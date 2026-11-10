import re
import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
#### Requires: pandas, numpy and Biopython


# ==================================================================================================================================================================================================
# Project: Proteomics Pipeline Function Collection
# Original Author: Patrick Ruether (Poriginal R code, updated by Claire Koenig, published in "Enamel proteins reveal biological sex and genetic variability in southern African Paranthropus" paper)
# Updated code and port to python: Ioanins Patramanis
# Script purpose: Function collection for FASTA and MaxQuant data handling
# ===================================================================================================================================================================================================


def fasta_to_datatable(path=".", paths=None):
    """
    Reads one or more FASTA protein sequence files, extracts standard UniProt/FASTA 
    header annotations (like Gene Name GN=, Organism OS=, and Organism ID OX=), 
    and compiles everything into a clean, well-structured pandas DataFrame.
    """
    
    
    ##### Create Fasta list 
    if paths is None:
        paths = []
    
    if len(paths) > 0:
        fasta_paths = paths
    else:
        fasta_paths = []
        if os.path.isfile(path): ### if input is a single file
            fasta_paths = [path]
        elif os.path.isdir(path):
            
            #### if input is a folder, search for .fa* files in the directory
            pattern = os.path.join(path, "*.fa*")
            fasta_paths = glob.glob(pattern)
            
            ### Also check uppercase extensions just in case
            if len(fasta_paths) == 0:
                fasta_paths = glob.glob(os.path.join(path, "*.FA*"))
        
        if len(fasta_paths) == 0:
            raise FileNotFoundError("No FASTA file in specified path")




    #### load fasta files using Biopython
    output_list = []
    #### For each fasta file in list, load it as an entry into a list of dictionaries
    for file_path in fasta_paths:
        records = []
        for record in SeqIO.parse(file_path, "fasta"):
            records.append({
                "name": record.id,
                "Annot": record.description,
                "Sequence": str(record.seq),
                "File": file_path
            })
        
        ### Convert list of dictionary of fasta entries to Pandas dataframe
        df_temp = pd.DataFrame(records)
        output_list.append(df_temp)
    
    if not output_list:
        raise ValueError("No sequences found in the provided FASTA files.")
    
    ### Join all fasta file dataframes into one dataframe
    output = pd.concat(output_list, ignore_index=True)


    #### Ensure required columns exist, if not add them as empy
    cols_to_init = ["Database", "Protein", "Gene_Species", "Gene", "Species", "SpeciesID"]
    for c_name in cols_to_init:
        if c_name not in output.columns:
            output[c_name] = None


    ##### Detect UniProt format, if it exists (if it contains '|' OR standard tags like OS= / OX=)
    output["is_uniptr"] = output["Annot"].str.contains(r"\|", na=False) | (output["Annot"].str.contains(r"OS=", na=False) & output["Annot"].str.contains(r"OX=", na=False))

    ##### --- MODE 1: UniProt Format --- #####
    uniprot_mask = output["is_uniptr"] == True
    if uniprot_mask.any():
        u_parts = output.loc[uniprot_mask, "Annot"].str.split(r"\|", n=2, expand=True)
        if u_parts.shape[1] >= 3: ### if 3 or more |parts|  Database is the first, protein the second and Gene species the third
            output.loc[uniprot_mask, "Database"] = u_parts[0]
            output.loc[uniprot_mask, "Protein"] = u_parts[1]
            output.loc[uniprot_mask, "Gene_Species"] = u_parts[2]
        
        ### Gene is grabbed from GN, Species from OS and SpeciesID from OX
        output.loc[uniprot_mask, "Gene"] = output.loc[uniprot_mask, "Annot"].str.extract(r"GN=([^\s]+)", expand=False)
        output.loc[uniprot_mask, "Species"] = output.loc[uniprot_mask, "Annot"].str.extract(r"OS=([^\s]+\s[^\s]+)", expand=False)
        output.loc[uniprot_mask, "SpeciesID"] = output.loc[uniprot_mask, "Annot"].str.extract(r"OX=([^\s]+)", expand=False)

        #### Fallback for missing Gene names in UniProt format
        missing_g = uniprot_mask & (output["Gene"].isna() | (output["Gene"] == ""))
        if missing_g.any():
            fallback_labels = output.loc[missing_g, "Annot"]
            ### This will grab whatever text appears after the last underscore in the header and assigns it as the Gene name
            output.loc[missing_g, "Gene"] = fallback_labels.str.replace(r"^.*_(.+)$", r"\1", regex = True) 



    ##### --- MODE 2: Non-UniProt Format --- #####
    custom_mask = output["is_uniptr"] == False
    if custom_mask.any():
        custom_labels = output.loc[custom_mask, "Annot"]
        prefix_part = custom_labels.str.replace(r"_(?!.*_).+$", "", regex = True)
        suffix_part = custom_labels.str.replace(r"^.*_(.+)$", r"\1", regex = True)
        
        ##### Everything before last '_' will be used for name, annotation, gene species, SpeciesID, Dataabase BUT everything after last '_' will be used for GENE name
        output.loc[custom_mask, "name"] = prefix_part
        output.loc[custom_mask, "Annot"] = prefix_part
        output.loc[custom_mask, "Gene_Species"] = prefix_part
        output.loc[custom_mask, "SpeciesID"] = prefix_part
        output.loc[custom_mask, "Database"] = "Custom_Database"
        output.loc[custom_mask, "Gene"] = suffix_part
    #### Remove column 'is uniprot'
    output.drop(columns=["is_uniptr"], inplace = True)
    
    #### Convert Gene to upper case for consistency 
    output["Gene"] = output["Gene"].str.upper()

    ##### If data has not 'Protein' labels, Handle empty Protein column entries with unique gene-based suffixes ####
    empty_prot = output["Protein"].isna() | (output["Protein"] == "")
    if empty_prot.any():
        for idx in output.index[empty_prot]:
            base_gene = output.loc[idx, "Gene"]
            if pd.isna(base_gene) or base_gene == "":
                base_gene = "UNKNOWN_GENE"
            counter = 1
            candidate = f"{base_gene}_{counter}"
            while candidate in output["Protein"].values:
                counter += 1
                candidate = f"{base_gene}_{counter}"
            output.loc[idx, "Protein"] = candidate
    
    ##########################
    #### Global uniqueness check for Protein column across all entries
    dup_prot = output["Protein"].duplicated(keep = False) ### marks all rows that have duplicate entry in Protein (keeping all occurances, not just first)
    if dup_prot.any():
        for idx in output.index[dup_prot]: ## go through duplicate entries
            base_prot = output.loc[idx, "Protein"]
            counter = 1
            candidate = f"{base_prot}_{counter}" ## new name (e.g. protein_1)
            while candidate in output["Protein"].values: ## if new name exists as well
                counter += 1
                candidate = f"{base_prot}_{counter}" ### increase counter and try again
            output.loc[idx, "Protein"] = candidate

    return output ## return dataframe







def mq_intensity_64(df):
    """
    Ensures that large quantitative columns from MaxQuant outputs (such as protein intensities, 
    iBAQ, and LFQ values) are stored as 64-bit integers (Int64) to prevent overflow errors, 
    and converts percentage columns to double precision (float).
    """
    df_copy = df.copy()
    intensity_cols = [col for col in df_copy.columns if re.search(r"(Intensity)|(iBAQ)|(LFQ)", col)]
    percent_cols = [col for col in df_copy.columns if "%" in col]

    for col in intensity_cols:
        if df_copy[col].dtype != "Int64":
            try:
                df_copy[col] = pd.to_numeric(df_copy[col], errors="coerce").astype("Int64")
            except Exception:
                pass

    for col in percent_cols:
        if df_copy[col].dtype != "float64":
            try:
                df_copy[col] = pd.to_numeric(df_copy[col], errors="coerce").astype("float64")
            except Exception:
                pass

    return df_copy









def most_common(x):
    """
    Finds the most frequently occurring element in a character/string sequence or list. 
    If there is a tie, it picks the shortest string; if lengths are also equal, 
    it selects the first one. Ignores empty strings ('').
    """
    if not isinstance(x, (list, pd.Series, np.ndarray)):
        raise TypeError("Input has to be a list, Series, or array of characters.")
    
    # Filter out empty strings and NA/NaN
    cleaned = [str(item) for item in x if item is not None and str(item) != "" and str(item).lower() != "nan"]
    
    if len(cleaned) == 0:
        return np.nan
    
    # Count frequencies
    counts = pd.Series(cleaned).value_counts()
    max_freq = counts.max()
    
    # Get all items with max frequency
    top_items = counts[counts == max_freq].index.tolist()
    
    if len(top_items) == 1:
        return top_items[0]
    else:
        # Tie-breaker 1: Shortest string length
        min_len = min(len(item) for item in top_items)
        shortest_items = [item for item in top_items if len(item) == min_len]
        # Tie-breaker 2: First one encountered
        return shortest_items[0]


def all_contain(x, pattern="CON__", sep=";"):
    """
    Checks multi-value delimited strings (e.g., protein IDs separated by semicolons) 
    to see if all individual entries within a string match a specific search pattern 
    (default: 'CON__').
    """
    if isinstance(x, str):
        x = [x]
    
    results = []
    for item in x:
        if pd.isna(item):
            results.append(False)
            continue
        if sep not in str(item):
            ### Warning equivalent
            pass
        parts = str(item).split(sep)
        matches = [bool(re.search(pattern, p)) for p in parts]
        results.append(all(matches) if len(matches) > 0 else False)
        
    return pd.Series(results) if len(results) > 1 else results[0]



def paste_list_cols(df, sep=";"):
    """
    Identifies list-type columns inside a pandas DataFrame and collapses their 
    elements into standard character strings by pasting them together with a separator.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input is not a pandas DataFrame")
    
    df_copy = df.copy()
    list_cols = [col for col in df_copy.columns if df_copy[col].apply(isinstance, args=(list,)).any()]
    
    for col in list_cols:
        df_copy[col] = df_copy[col].apply(lambda val: sep.join(str(v) for v in val) if isinstance(val, list) else val)
        
    return df_copy




def melt_mq_prot(df, appendix="_long"):
    """
    Reshapes a wide-format MaxQuant proteinGroups.txt table into a long-format table. 
    Melts experimental columns by experiment name, pairing them neatly with core protein metadata.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input has to be a pandas DataFrame (read from proteinGroups.txt)")

    has_ibaq = any(re.match(r"^iBAQ", col) for col in df.columns)
    has_lfq = any(re.match(r"^LFQ", col) for col in df.columns)
    has_msms = any(re.match(r"^MS/MS count\s.+", col) for col in df.columns)

    melt_groups = {
        "Peptides": [c for c in df.columns if re.match(r"^Peptides .+", c)],
        "Razor + unique peptides": [c for c in df.columns if re.match(r"^Razor \+ unique peptides .+", c)],
        "Unique peptides": [c for c in df.columns if re.match(r"^Unique peptides .+", c)],
        "Sequence coverage": [c for c in df.columns if re.match(r"^Sequence coverage [^\[]", c)],
        "Intensity": [c for c in df.columns if re.match(r"^Intensity .+", c)],
    }

    if has_msms:
        melt_groups["MS/MS count"] = [c for c in df.columns if re.match(r"^MS/MS count\s.+", c)]
    if has_lfq:
        melt_groups["LFQ"] = [c for c in df.columns if re.match(r"^LFQ .+", c)]
    if has_ibaq:
        melt_groups["iBAQ"] = [c for c in df.columns if re.match(r"^iBAQ .+", c)]

    # Check that all melt groups have the exact same size
    lengths = [len(v) for v in melt_groups.values()]
    if len(set(lengths)) != 1:
        raise ValueError("Something wrong with your column names: mismatched metric group sizes.")

    # Extract experiment names from the 'Peptides ' columns
    pep_cols = melt_groups["Peptides"]
    experiments = [re.search(r"Peptides\s(.+)", c).group(1) for c in pep_cols if re.search(r"Peptides\s(.+)", c)]

    # Base metadata columns to retain
    id_vars = [c for c in ["Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "Fasta headers", "Number of proteins"] if c in df.columns]

    # Build long format by melting each category and merging or aligning
    dfs_to_melt = []
    for metric, cols in melt_groups.items():
        temp_df = df[id_vars + cols].copy()
        # Rename measurement columns to just the experiment name
        col_rename = {c: re.sub(rf"^{re.escape(metric)}\s+", "", c) for c in cols}
        temp_df = temp_df.rename(columns=col_rename)
        
        melted = temp_df.melt(
            id_vars=id_vars,
            value_vars=list(col_rename.values()),
            var_name="Experiment",
            value_name=metric + appendix
        )
        dfs_to_melt.append(melted)

    # Merge all melted metric dataframes together on [id_vars + ['Experiment']]
    long_df = dfs_to_melt[0]
    for other_df in dfs_to_melt[1:]:
        long_df = pd.merge(long_df, other_df, on=id_vars + ["Experiment"], how="outer")

    return long_df