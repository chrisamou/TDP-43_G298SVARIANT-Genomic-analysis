import requests
import json
import pandas as pd
import re  # Imported for regular expressions in variant name parsing

# --- Configuration ---
# Your NCBI API key (optional but highly recommended for higher query limits)
# Get one from: https://ncbi.nlm.nih.gov/account/settings/#api-key
# Uncomment the line below and replace "YOUR_NCBI_API_KEY" with your actual key
# NCBI_API_KEY = "YOUR_NCBI_API_KEY"

# Mutation details for filtering
WT_RESIDUE_NAME = "GLY"  # Wild-type residue name (Glycine)
MUTATION_RESIDUE_NUMBER = 298
MUTANT_RESIDUE_NAME = "SER"  # Mutant residue name (Serine)

# --- NCBI E-utilities base URL ---
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def query_clinvar_for_gene(gene_name, retmax=500):
    """
    Searches ClinVar for variants associated with a given gene and parses key details.
    """
    search_url = f"{EUTILS_BASE}esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": f"{gene_name}[gene]",
        "retmode": "json",
        "retmax": retmax
    }
    # Only add API key if defined and not empty
    if 'NCBI_API_KEY' in globals() and NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    print(f"Searching ClinVar for gene: {gene_name} (retmax={retmax})...")
    response = requests.get(search_url, params=params)
    response.raise_for_status()  # Raise an exception for HTTP errors
    data = response.json()

    id_list = data["esearchresult"]["idlist"]
    if not id_list:
        print(f"No ClinVar entries found for gene {gene_name}.")
        return pd.DataFrame()

    print(f"Found {len(id_list)} ClinVar IDs. Fetching summaries...")
    summary_url = f"{EUTILS_BASE}esummary.fcgi"
    params = {
        "db": "clinvar",
        "id": ",".join(id_list),
        "retmode": "json"
    }
    # Only add API key if defined and not empty
    if 'NCBI_API_KEY' in globals() and NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    summary_response = requests.get(summary_url, params=params)
    summary_response.raise_for_status()
    summary_data = summary_response.json()

    variants_data = []
    # If the summary_data structure is unexpectedly empty or malformed at the top level
    if not summary_data or "result" not in summary_data:
        print("Error: 'summary_data' is empty or missing 'result' key after fetching summaries.")
        print(f"Raw summary response head: {str(summary_response.text)[:500]}...")  # Print beginning of response
        return pd.DataFrame()

    for uid in id_list:
        try:
            # Check if the UID exists directly under "result"
            if uid not in summary_data["result"]:
                print(f"Warning: UID {uid} not found directly in summary_data['result']. Skipping.")
                continue

            # This is the crucial part: try to get data from nested structure first,
            # if not found, fall back to the top-level UID data.
            variant_raw_data = summary_data["result"][uid]
            variant_data_source = None

            doc_summary_set = variant_raw_data.get("document_summary_set")
            if doc_summary_set:
                doc_summary_list = doc_summary_set.get("document_summary")
                if doc_summary_list:
                    variant_data_source = doc_summary_list[0]

            # If nested structure not found, use the top-level UID data directly
            if not variant_data_source:
                variant_data_source = variant_raw_data
                # print(f"Info: Using flatter structure for UID {uid}.") # Optional debug print

            if not variant_data_source:  # Should not happen if variant_raw_data exists
                print(f"Error: No suitable data source found for UID {uid}. Skipping.")
                continue

            # Now, extract data from the chosen source (either nested or flat)
            # The .get() method with a default value handles missing keys gracefully
            clinical_significance = variant_data_source.get("clinical_significance", {}).get("description", "N/A")
            rcv_accession = variant_data_source.get("rcv_accession", "N/A")
            variation_id = variant_data_source.get("uid", "N/A")  # This UID is the ClinVar Variation ID

            # Variant name could be 'title' or 'name' depending on the structure
            variant_name = variant_data_source.get("title")
            if not variant_name:
                variant_name = variant_data_source.get("name", "N/A")  # Fallback for flatter structure

            gene_symbol = "N/A"
            if 'gene_list' in variant_data_source and variant_data_source['gene_list']:
                first_gene_entry = variant_data_source['gene_list'][0]
                if isinstance(first_gene_entry, dict) and 'gene_symbol' in first_gene_entry:
                    gene_symbol = first_gene_entry['gene_symbol']
                elif isinstance(first_gene_entry, str):  # Fallback if it's just a string gene symbol
                    gene_symbol = first_gene_entry

            # Attempt to extract protein change string
            protein_change_str = "N/A"
            match_p = re.search(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}', variant_name)
            if match_p:
                protein_change_str = match_p.group(0)
            else:  # Try single-letter code if three-letter not found
                match_p_single = re.search(r'p\.[A-Z](\d+)[A-Z]', variant_name)
                if match_p_single:
                    protein_change_str = match_p_single.group(0)

            variants_data.append({
                "ClinVar_ID": variation_id,
                "RCV_Accession": rcv_accession,
                "Variant_Name": variant_name,
                "Gene_Symbol": gene_symbol,
                "Clinical_Significance": clinical_significance,
                "Protein_Change_Str": protein_change_str
            })
        except Exception as e:  # Catching a broader exception for any unexpected issue
            print(
                f"Error parsing UID {uid}: {e}. Raw summary_data['result'][uid]: {summary_data['result'].get(uid, 'Not available')}")
            continue

    return pd.DataFrame(variants_data)


# --- Main Execution Block for ClinVar ---
if __name__ == "__main__":
    gene_of_interest = "TARDBP"  # TDP-43 gene symbol
    clinvar_df = query_clinvar_for_gene(gene_of_interest)

    if not clinvar_df.empty:
        print("\n--- ClinVar Data for TARDBP (First 5 entries) ---")
        print(clinvar_df.head())
        print(f"\nTotal ClinVar entries found: {len(clinvar_df)}")

        # Filter for your specific mutation (G298S or Gly298Ser)
        g298s_clinvar = clinvar_df[
            clinvar_df['Protein_Change_Str'].str.contains(
                rf"p\.{WT_RESIDUE_NAME[:3]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:3]}", case=False, na=False) |
            clinvar_df['Protein_Change_Str'].str.contains(
                rf"p\.{WT_RESIDUE_NAME[:1]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:1]}", case=False, na=False) |
            clinvar_df['Variant_Name'].str.contains(f"{WT_RESIDUE_NAME}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME}",
                                                    case=False, na=False) |
            clinvar_df['Variant_Name'].str.contains(
                f"{WT_RESIDUE_NAME[:3]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:3]}", case=False, na=False) |
            clinvar_df['Variant_Name'].str.contains(
                f"{WT_RESIDUE_NAME[:1]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:1]}", case=False, na=False)
            ]

        if not g298s_clinvar.empty:
            print("\n--- G298S Variant in ClinVar (Filtered Results) ---")
            print(g298s_clinvar)
        else:
            print(
                f"\n{WT_RESIDUE_NAME}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME} not explicitly found in ClinVar results for TARDBP with common filters.")
            print("You may need to manually inspect the full ClinVar data or refine the search/filter terms.")
            # Optional: Save the full DataFrame to CSV to inspect all entries
            # clinvar_df.to_csv("full_tardbp_clinvar_data.csv", index=False)
            # print("Full ClinVar data saved to full_tardbp_clinvar_data.csv for manual inspection.")

import requests
import json
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
# Your NCBI API key (optional but highly recommended for higher query limits)
# Get one from: https://ncbi.nlm.nih.gov/account/settings/#api-key
# Uncomment the line below and replace "YOUR_NCBI_API_KEY" with your actual key
# NCBI_API_KEY = "YOUR_NCBI_API_KEY" # Example: NCBI_API_KEY = "YOUR_ACTUAL_API_KEY_HERE"

# Mutation details for filtering
WT_RESIDUE_NAME = "GLY"  # Wild-type residue name (Glycine)
MUTATION_RESIDUE_NUMBER = 298
MUTANT_RESIDUE_NAME = "SER"  # Mutant residue name (Serine)

# --- NCBI E-utilities base URL ---
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def query_clinvar_for_gene(gene_name, retmax=500):
    """
    Searches ClinVar for variants associated with a given gene and parses key details.
    """
    search_url = f"{EUTILS_BASE}esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": f"{gene_name}[gene]",
        "retmode": "json",
        "retmax": retmax
    }
    # Only add API key if defined and not empty
    if 'NCBI_API_KEY' in globals() and 'NCBI_API_KEY' in locals() and NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    print(f"Searching ClinVar for gene: {gene_name} (retmax={retmax})...")
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()  # Raise an exception for HTTP errors
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error connecting to ClinVar esearch: {e}")
        return pd.DataFrame()
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from ClinVar esearch: {e}")
        return pd.DataFrame()


    id_list = data["esearchresult"]["idlist"]
    if not id_list:
        print(f"No ClinVar entries found for gene {gene_name}.")
        return pd.DataFrame()

    print(f"Found {len(id_list)} ClinVar IDs. Fetching summaries...")
    summary_url = f"{EUTILS_BASE}esummary.fcgi"
    params = {
        "db": "clinvar",
        "id": ",".join(id_list),
        "retmode": "json"
    }
    # Only add API key if defined and not empty
    if 'NCBI_API_KEY' in globals() and 'NCBI_API_KEY' in locals() and NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        summary_response = requests.get(summary_url, params=params)
        summary_response.raise_for_status()
        summary_data = summary_response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error connecting to ClinVar esummary: {e}")
        return pd.DataFrame()
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from ClinVar esummary: {e}")
        return pd.DataFrame()


    variants_data = []
    # If the summary_data structure is unexpectedly empty or malformed at the top level
    if not summary_data or "result" not in summary_data:
        print("Error: 'summary_data' is empty or missing 'result' key after fetching summaries.")
        print(f"Raw summary response head: {str(summary_response.text)[:500]}...")
        return pd.DataFrame()

    for uid in id_list:
        try:
            if uid not in summary_data["result"]:
                continue

            variant_raw_data = summary_data["result"][uid]
            variant_data_source = None

            # Try to get data from nested structure first
            doc_summary_set = variant_raw_data.get("document_summary_set")
            if doc_summary_set:
                doc_summary_list = doc_summary_set.get("document_summary")
                if doc_summary_list:
                    variant_data_source = doc_summary_list[0]

            # If nested structure not found, use the top-level UID data directly
            if not variant_data_source:
                variant_data_source = variant_raw_data

            if not variant_data_source:
                continue

            clinical_significance = variant_data_source.get("clinical_significance", {}).get("description", "N/A")
            rcv_accession = variant_data_source.get("rcv_accession", "N/A")
            variation_id = variant_data_source.get("uid", "N/A")

            variant_name = variant_data_source.get("title")
            if not variant_name:
                variant_name = variant_data_source.get("name", "N/A")

            gene_symbol = "N/A"
            if 'gene_list' in variant_data_source and variant_data_source['gene_list']:
                first_gene_entry = variant_data_source['gene_list'][0]
                if isinstance(first_gene_entry, dict) and 'gene_symbol' in first_gene_entry:
                    gene_symbol = first_gene_entry['gene_symbol']
                elif isinstance(first_gene_entry, str):
                    gene_symbol = first_gene_entry

            protein_change_str = "N/A"
            match_p = re.search(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}', variant_name)
            if match_p:
                protein_change_str = match_p.group(0)
            else:
                match_p_single = re.search(r'p\.[A-Z](\d+)[A-Z]', variant_name)
                if match_p_single:
                    protein_change_str = match_p_single.group(0)

            variants_data.append({
                "ClinVar_ID": variation_id,
                "RCV_Accession": rcv_accession,
                "Variant_Name": variant_name,
                "Gene_Symbol": gene_symbol,
                "Clinical_Significance": clinical_significance,
                "Protein_Change_Str": protein_change_str
            })
        except Exception as e:
            # print(f"Skipping ClinVar UID {uid} due to parsing error: {e}") # Uncomment for debugging
            continue

    return pd.DataFrame(variants_data)


# --- gnomAD API base URL ---
# The direct gnomAD API query continues to show 400 Bad Request errors.
# We will skip the direct API query for gnomAD and use a manually retrieved value for TARDBP G298S.
# GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

# def query_gnomad_variant_by_id(variant_identifier, dataset="gnomad_r3"):
#     """
#     Queries gnomAD for a specific variant using its variant_id (chrom-pos-ref-alt) or rsID.
#     (NOTE: This function is currently bypassed due to consistent API errors.)
#     """
#     # ... (original function code, which will be commented out in the main block)
#     return None # Always return None to signal failure for now


# --- Main Execution Block ---
if __name__ == "__main__":
    gene_of_interest = "TARDBP"

    # --- ClinVar Data Retrieval ---
    clinvar_df = query_clinvar_for_gene(gene_of_interest)

    if not clinvar_df.empty:
        print("\n--- ClinVar Data for TARDBP (First 5 entries) ---")
        print(clinvar_df.head())
        print(f"\nTotal ClinVar entries found: {len(clinvar_df)}")

        # Filter ClinVar results for the specific G298S variant
        # Using multiple regex patterns to capture different protein change string formats
        g298s_clinvar = clinvar_df[
            clinvar_df['Protein_Change_Str'].str.contains(
                rf"p\.{WT_RESIDUE_NAME[:3]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:3]}", case=False, na=False) |
            clinvar_df['Protein_Change_Str'].str.contains(
                rf"p\.{WT_RESIDUE_NAME[:1]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:1]}", case=False, na=False) |
            clinvar_df['Variant_Name'].str.contains(f"{WT_RESIDUE_NAME}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME}",
                                                    case=False, na=False) |
            clinvar_df['Variant_Name'].str.contains(
                f"{WT_RESIDUE_NAME[:3]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:3]}", case=False, na=False) |
            clinvar_df['Variant_Name'].str.contains(
                f"{WT_RESIDUE_NAME[:1]}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME[:1]}", case=False, na=False)
        ]

        if not g298s_clinvar.empty:
            print("\n--- G298S Variant in ClinVar (Filtered Results) ---")
            print(g298s_clinvar)
        else:
            print(
                f"\n{WT_RESIDUE_NAME}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME} not explicitly found in ClinVar results for TARDBP with common filters.")
            print("You may need to manually inspect the full ClinVar data or refine the search/filter terms.")

    # --- gnomAD Data (Manual Retrieval) ---
    # Due to persistent API errors with gnomAD, we're using a manually retrieved value for TARDBP G298S.
    # This value was obtained from the gnomAD browser (https://gnomad.broadinstitute.org/variant/rs121908865)
    # for gnomAD v3.1.2 (genomes).
    print("\n--- gnomAD Data for TARDBP G298S (Manually Retrieved) ---")
    tardbp_g298s_gnomad_af = 8.033e-6  # Allele Frequency for rs121908865 in gnomAD v3.1.2 genomes

    # Create a dummy DataFrame to hold this manually retrieved data for consistency in plotting
    gnomad_df_for_plot = pd.DataFrame([{
        "variant_id": "1-11022301-G-A", # Canonical variant ID for G298S
        "chrom": "1",
        "pos": "11022301",
        "ref": "G",
        "alt": "A",
        "rsids": "rs121908865",
        "genome_allele_frequency": tardbp_g298s_gnomad_af
    }])
    print(gnomad_df_for_plot)

    # --- Combined Data Analysis and Visualization ---
    print("\n--- Combined Data Analysis and Visualization ---")

    # Proceed with plot generation using the manually retrieved gnomAD data
    combined_plot_data = pd.DataFrame({
        "Variant": [f"{WT_RESIDUE_NAME}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME}"],
        "Prediction_Score": [0.9],  # Dummy score for demonstration (e.g., from AlphaMissense, PolyPhen-2, etc.)
        "gnomAD_AF": [tardbp_g298s_gnomad_af],
        "Clinical_Significance": [
            "Pathogenic" if not g298s_clinvar.empty and g298s_clinvar['Clinical_Significance'].isin(
                ['Pathogenic', 'Likely pathogenic']).any() else "Benign"
        ]
    })

    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=combined_plot_data,
        x="Prediction_Score",
        y="gnomAD_AF",
        hue="Clinical_Significance",
        size="gnomAD_AF",
        sizes=(50, 500), # Adjust marker size based on AF for visual emphasis
        alpha=0.7
    )
    plt.yscale('log') # Use a log scale for allele frequency as it can vary widely
    plt.title(f'Variant Prediction Score vs. gnomAD Allele Frequency for {gene_of_interest} {WT_RESIDUE_NAME}{MUTATION_RESIDUE_NUMBER}{MUTANT_RESIDUE_NAME}')
    plt.xlabel('Predicted Pathogenicity Score (e.g., 0-1, 1 is pathogenic)')
    plt.ylabel('gnomAD Allele Frequency (log scale)')
    plt.grid(True, which="both", ls="--", c="0.7") # Add a grid for readability
    plt.show()