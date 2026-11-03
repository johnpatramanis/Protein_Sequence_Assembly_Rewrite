import os
from pathlib import Path
import pandas as pd
import numpy as np
import re
import glob
from Bio import SeqIO
#### Requires: pandas, numpy and Biopython

# ==========================================================
# Project: MaxQuant Analysis Pipeline
# Original Author: Patrick Ruether (updated by Claire Koenig)
# Converted to Python 3 by Ioannis Patramanis
# Code Updates by Ioannis Patramanis in the year of our lord 2026
# Script purpose: Environment setup, data loading, and FASTA processing
# ==========================================================

#### --- Environment Setup & Working Directory --- ####
#### Safely get the directory where this script is located (handles command-line & interactive execution)
try:
    script_dir = Path(__file__).resolve().parent
except NameError:
    # Fallback for interactive environments like Jupyter notebooks
    script_dir = Path.cwd()

##### Set working directory to the script's location
os.chdir(script_dir)

## Import custom functions from your functions module
from Functions import (fasta_to_datatable, mq_intensity_64, most_common, all_contain, paste_list_cols, melt_mq_prot)


## ==========================================================
## 0. Setting up paths and loading files from MaxQuant
## ==========================================================

maxquant_dir = "../Test_Input_Data/txt-ParanthropusSearch"
samples_to_experiment_path = "../Test_Input_Data/Sample_to_Experiment.tsv"

fasta_path = "../Test_Input_Data/PR220121_GreatApes_Enamelome_aligned.fasta"
exclusion_path = "../Test_Input_Data/CK_ExclusionListSample.csv"

output_dir = "../Output_and_Input_Paranthropus_Yiannis_Reanalysis_Python"

# Create output directory if it does not exist (equivalent to dir.create with recursive=TRUE)
os.makedirs(output_dir, exist_ok = True)


## ==========================================================
## 0) LOAD TABLES
## ==========================================================

### a) Load MaxQuant Evidence File
evidence = pd.read_csv(os.path.join(maxquant_dir, "evidence.txt"), sep = "\t", low_memory = False)

### b) Load MaxQuant Summary File and filter out rows where Experiment is empty
summary = pd.read_csv(os.path.join(maxquant_dir, "summary.txt"), sep = "\t", low_memory = False)
summary = summary[summary["Experiment"].notna() & (summary["Experiment"] != "")]

### c) Load MaxQuant MSMS File
MSMS = pd.read_csv(os.path.join(maxquant_dir, "msms.txt"), sep = "\t", low_memory = False)




### d) Load FASTA & Generate Site-Level Table

### Load FASTA using our previously defined custom function 
### (note: pass as single file or a folder with fasta)
FASTA = fasta_to_datatable(fasta_path)

### Generate gapless sequences (remove all hyphens '-' and 'X')
FASTA["Gapless"] = FASTA["Sequence"].str.replace("-", "", regex = False)
FASTA["Gapless"] = FASTA["Gapless"].str.replace("X", "", regex = False)

### Add a sequential identifier (1 to N) for each sequence row
FASTA["id"] = range(1, len(FASTA) + 1)


### f) Generate Site-Level FASTA Table
### Splits each protein sequence character-by-character into individual amino acid rows.

sites_records = []

for INDX, row in FASTA.iterrows(): ## Row is the whole row as dataframe
    seq = row["Sequence"]
    aa_vec = list(seq) ### every AA is an item in this list
    aligned_pos = list(range(1, len(aa_vec) + 1)) ### numbers with the position of each AA of aa_vec
    
    ## Calculate gapless positions: 
    ## Non-gap characters get a sequential index (1, 2, 3...), gaps get NaN
    gapless_pos = []
    curr_gapless = 0
    for aa in aa_vec:
        if (aa != "-") and (aa != "X"):
            curr_gapless += 1
            gapless_pos.append(curr_gapless)
        else:
            gapless_pos.append(np.nan)
            
    ### Bundle into records for each position
    for aa, apos, gpos in zip(aa_vec, aligned_pos, gapless_pos):
        sites_records.append({
            "Protein": row["Protein"],
            "Annot": row["Annot"],
            "Gene": row["Gene"],
            "Species": row["Species"],
            "id": row["id"],
            "AA": aa,
            "AlignedPos": apos,
            "GaplessPos": gpos
        })

#### Create the expanded site-level DataFrame
FASTAsites = pd.DataFrame(sites_records)

### Add identifier for each amino acid site: Gene_AlignedPosition
FASTAsites["GeneSite"] = FASTAsites["Gene"].astype(str) + "_" + FASTAsites["AlignedPos"].astype(str)










# ==========================================================
# 1. Experimental annotation from summary table
# ==========================================================

### Create new DataFrame from Summary where 'Enzyme mode' is non-empty
### Equivalent to: Summary[nzchar(`Enzyme mode`), .(Experiment, `Raw file`, Fraction)]
enzyme_mode_valid = summary["Enzyme mode"].notna() & (summary["Enzyme mode"].astype(str).str.strip() != "") ### List of rows that have 'enzyme mode'
ExpAnnot = summary.loc[enzyme_mode_valid, ["Experiment", "Raw file", "Fraction"]].copy() ### Get experiment, raw file and fraction for rows that have enzyme mode

### Sort by 'Raw file'
ExpAnnot = ExpAnnot.sort_values(by="Raw file").reset_index(drop=True)


### Load samples to experiments TSV file, if it exists , map experiments to samples (e.g. if multiple experiements analysed the same sample) ###
if os.path.exists(samples_to_experiment_path):
    ### Read samples to Experiment file
    samples_to_experiment = pd.read_csv(samples_to_experiment_path, sep="\t")
    
    
    ### Ensure standard column names (1st col = Sample, 2nd col = Experiment)
    if len(samples_to_experiment.columns) >= 2:
        samples_to_experiment.columns = ["Sample", "Experiment"] + list(samples_to_experiment.columns[2:])
    
    ### Merge mapping into ExpAnnot by Experiment (equivalent to all.x = TRUE / left join)
    ExpAnnot = pd.merge(
        ExpAnnot, 
        samples_to_experiment[["Experiment", "Sample"]], 
        on="Experiment", 
        how="left"
    )
    
    #### Fallback for any experiments missing from the mapping file (fill NA spots for 'sample' with Experiment name)
    ExpAnnot["Sample"] = ExpAnnot["Sample"].fillna(ExpAnnot["Experiment"])
else:
    #### Fallback if no mapping file is provided
    ExpAnnot["Sample"] = ExpAnnot["Experiment"]


### Convert Sample Column to Category
ExpAnnot["Sample"] = ExpAnnot["Sample"].astype("category")


## Add fractionation information (single shot vs. fractionated)
## New column 
## Count entries per Experiment (equivalent to .N by = Experiment)
ExpAnnot["n_entries"] = ExpAnnot.groupby("Experiment")["Experiment"].transform("count")

## New column
## Ensure Fraction can be evaluated numerically
ExpAnnot["Fraction_num"] = pd.to_numeric(ExpAnnot["Fraction"], errors = "coerce")

## Determine if fractionated based on multiple entries or Fraction > 1
## List of booleans
is_fractionated = (ExpAnnot["n_entries"] > 1) | (ExpAnnot["Fraction_num"].fillna(0) > 1)

## Where is_fractionated == True do frc, where Flase Single, then convert to category
ExpAnnot["Fractionated"] = np.where(is_fractionated, "frc", "single")
ExpAnnot["Fractionated"] = ExpAnnot["Fractionated"].astype("category")

### For fractionated samples, cast Fraction to integer
frac_mask = ExpAnnot["Fractionated"] == "frc"
if frac_mask.any():
    ExpAnnot.loc[frac_mask, "Fraction"] = ExpAnnot.loc[frac_mask, "Fraction_num"].astype("Int64")

### Clean up temporary helper columns
ExpAnnot = ExpAnnot.drop(columns=["n_entries", "Fraction_num"])

### Convert Experiment column to a category
ExpAnnot["Experiment"] = ExpAnnot["Experiment"].astype("category")





# ==========================================================
# 2. Annotate tables, by merging them with Exp Annot
# ==========================================================



### Merge ExpAnnot into MSMS on 'Raw file' (equivalent to R's default inner join merge)
MSMS = pd.merge(MSMS, ExpAnnot, on="Raw file", how="inner")

### Drop 'Experiment' and 'Fraction' columns from evidence (if they exist), 
### then merge with ExpAnnot on 'Raw file'
evidence_subset = evidence.drop(columns=["Experiment", "Fraction"], errors="ignore")
evidence = pd.merge(evidence_subset, ExpAnnot, on="Raw file", how="inner")



# ==========================================================
# 3. Annotate tables, by merging them with Exp Annot
# ==========================================================


###
### Check if PepSeq file already exists in the output directory
pepseq_files = [f for f in os.listdir(output_dir) if "PepSeq" in f] if os.path.exists(output_dir) else []

### If it does, simply load it
if len(pepseq_files) > 0:
    target_file = "PepSeqs.csv" if "PepSeqs.csv" in pepseq_files else pepseq_files[0]
    PepSeqs = pd.read_csv(os.path.join(output_dir, target_file))
### Else, we need to create it
else:
    
    ### Get all unique peptide sequences from MSMS table
    unique_peptides = MSMS["Sequence"].dropna().unique()
    
    ###
    pep_records = []
    id_counter = 1
    
    ### For each peptide sequence,
    ### Find all proteins in FASTA (looking at gapless sequences) containing this peptide
    for seq in unique_peptides:
        
        ### Get true or false, which sequences from FASTA contain this peptide 
        mask = FASTA["Gapless"].str.contains(re.escape(seq), na=False) 
        matched_rows = FASTA[mask] ### subset FASTA to those that contain it
        
        ########### loop through matched rows (if they exist)
        for _, f_row in matched_rows.iterrows():
            
            prot = f_row["Protein"] ## name of protein (e.g. protein_1) of this row
            gapless_seq = f_row["Gapless"] ## Get the gapless sequence
            
            ### Find all start positions of the peptide in the gapless sequence
            for match in re.finditer(re.escape(seq), gapless_seq):
                start_pos = match.start() + 1  # 1-based index        ###############################################################################################################                         
                
                #### Expand to site level: 1 row per amino acid for each peptide-protein match
                for offset in range(len(seq)):
                    pep_records.append({
                        "id": id_counter,
                        "Sequence": seq,
                        "Protein": prot,
                        "GaplessPos": start_pos + offset
                    })
                id_counter += 1

    PepSeqs_sites = pd.DataFrame(pep_records) ### Dataframe, one AA per row. AA are grouped together by ID, if they belong to the same peptide
    
    if not PepSeqs_sites.empty: ## if there is any matches
        
        ### Filter FASTAsites where AA != "-" and select key columns
        fasta_sites_filtered = FASTAsites[FASTAsites["AA"] != "-"][["Gene", "Protein", "AlignedPos", "GaplessPos"]]
        
        ### Merge with FASTAsites to add AlignedPos and Gene
        PepSeqs_merged = pd.merge(
            PepSeqs_sites, 
            fasta_sites_filtered, 
            on=["Protein", "GaplessPos"], 
            how="left"
        )
        
        ### Go back to peptide level: collapse GaplessPos and AlignedPos into comma-separated strings
        grouped1 = PepSeqs_merged.groupby(["Gene", "Protein", "Sequence", "id"]).agg(
            GaplessPos=('GaplessPos', lambda x: ','.join(map(str, x))),
            AlignedPos=('AlignedPos', lambda x: ','.join(map(str, x)))
        ).reset_index()
        
        ### Filter down to unique peptide sequence, Gene, and AlignedPos
        grouped2 = grouped1.groupby(["Sequence", "Gene", "AlignedPos"]).agg(
            Proteins=('Protein', lambda x: ';'.join(sorted(x.unique()))),
            GaplessPositions=('GaplessPos', lambda x: ';'.join(x))
        ).reset_index()
        
        ### Check for peptides that match multiple proteins per gene
        ### This is counting how often a peptide is found within the same protein group (e.g. matches ALB_1, ALB_2 and ALB_3 multimatch = 3)
        grouped2["Multimatch"] = grouped2.groupby(["Sequence", "Gene"])["Sequence"].transform("count")
        
        ### Reorder and get rid of unique alignment artifacts based on comma count in Proteins
        grouped2["prot_comma_count"] = grouped2["Proteins"].str.count(",")
        grouped2 = grouped2.sort_values(by=["Sequence", "Gene", "prot_comma_count"])
        
        # Collapse by Gene and Sequence
        PepSeqs = grouped2.groupby(["Gene", "Sequence"]).agg(
            Proteins=('Proteins', lambda x: ';'.join(x)),
            AlignedPos=('AlignedPos', lambda x: x.iloc[0]),
            GaplessPositions=('GaplessPositions', lambda x: x.iloc[0])
        ).reset_index()
        
        #### Calculate how often a peptide sequence is repeated in the table (Multiplicity)
        ### This is counting how often a peptide is found between different protein groups (e.g. matches COL1_A1 and COL1_A2)
        PepSeqs["Multiplicity"] = PepSeqs.groupby("Sequence")["Sequence"].transform("count")
        
        #### Apply razor principle (winner takes all) to peptides mapping to multiple genes
        multi_mask = PepSeqs["Multiplicity"] > 1
        if multi_mask.any():
            ### Flag genes that have any multi-mapping peptides
            PepSeqs["AnyMultiPep"] = PepSeqs.groupby("Gene")["Multiplicity"].transform(lambda x: (x > 1).any())
            
            ### Subset multi-peptides and calculate peptide count per gene across the dataset
            multi_subset = PepSeqs[PepSeqs["Multiplicity"] > 1].copy()
            pep_count_df = PepSeqs[PepSeqs["AnyMultiPep"] == True].groupby("Gene").size().reset_index(name="PepCount")
            
            multi_subset = pd.merge(multi_subset, pep_count_df, on="Gene", how="left")
            
            ### Find max peptide count per sequence
            multi_subset["MaxPepCount"] = multi_subset.groupby("Sequence")["PepCount"].transform("max")
            
            ### Keep razor gene with the most peptides in the entire dataset
            MultiPeps = multi_subset[multi_subset["PepCount"] == multi_subset["MaxPepCount"]].copy()
            
            ### Combine single multiplicity peptides with filtered multi-peptides
            single_peps = PepSeqs[PepSeqs["Multiplicity"] == 1].copy()
            
            cols_to_drop = ["MaxPepCount", "PepCount", "AnyMultiPep"]
            MultiPeps = MultiPeps.drop(columns=[c for c in cols_to_drop if c in MultiPeps.columns])
            single_peps = single_peps.drop(columns=[c for c in cols_to_drop if c in single_peps.columns])
            
            PepSeqs = pd.concat([single_peps, MultiPeps], ignore_index=True)
          
          
        #### Clean up temporary helper columns if present
        if "prot_comma_count" in PepSeqs.columns:
            PepSeqs = PepSeqs.drop(columns=["prot_comma_count"])
            
        # Export site table as CSV into the output folder
        os.makedirs(output_dir, exist_ok=True)
        PepSeqs.to_csv(os.path.join(output_dir, "PepSeqs.csv"), index=False)
        
        
        
        
        



# ==========================================================
## 4. Delta Score recalculation and cutoff optimization
# ==========================================================
