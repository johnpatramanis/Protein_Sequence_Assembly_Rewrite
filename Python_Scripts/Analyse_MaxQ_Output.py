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


## ================
## LOAD TABLES
## ================


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
## 3. Annotate tables, by merging them with Exp Annot
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
    pep_records = [] ### list of dictionaries, one entry per AA, per sequence match, per peptide (so peptides with no matches don't get here)
    id_counter = 1
    
    ### For each peptide sequence,
    ### Find all proteins in FASTA (looking at gapless sequences) containing this peptide
    for seq in unique_peptides:
        
        ### Get true or false, which sequences from FASTA contain this peptide 
        mask = FASTA["Gapless"].str.contains(re.escape(seq), na = False) 
        matched_rows = FASTA[mask] ### subset FASTA to those that contain it
        
        ########### loop through matched rows (if they exist)
        for _, f_row in matched_rows.iterrows():
            
            prot = f_row["Protein"] ## name of protein (e.g. protein_1) of this row
            gapless_seq = f_row["Gapless"] ## Get the gapless sequence
            
            ### Find all start positions of the peptide in the gapless sequence
            for match in re.finditer(re.escape(seq), gapless_seq):
                start_pos = match.start() + 1  ### +1 to convert python counting to real, protein length counting                     
                
                #### Expand to site level: 1 row per amino acid for each peptide-protein match
                for offset in range(len(seq)):
                    pep_records.append({
                        "id": id_counter,
                        "Sequence": seq,
                        "Protein": prot,
                        "GaplessPos": start_pos + offset
                    })
                id_counter += 1

    PepSeqs_sites = pd.DataFrame(pep_records) ### Convert Dataframe, one AA per row. AA are grouped together by ID, if they belong to the same peptide
    
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
        grouped1 = PepSeqs_merged.groupby(["Gene", "Protein", "Sequence", "id"]).agg( #### For AA that share Gene, Protein, Sequence and ID (so result should be one entry = 1 peptide match on one protein)
            GaplessPos=('GaplessPos', lambda x: ','.join(map(str, x))),
            AlignedPos=('AlignedPos', lambda x: ','.join(map(str, x)))
        ).reset_index()
        
        ### Filter down to unique peptide sequence, Gene, and AlignedPos, (so result should be one entry = 1 peptide and all of its matches together)
        grouped2 = grouped1.groupby(["Sequence", "Gene", "AlignedPos"]).agg(
            Proteins=('Protein', lambda x: ';'.join(sorted(x.unique()))),
            GaplessPositions=('GaplessPos', lambda x: ';'.join(x))
        ).reset_index()
        
        ### Check for peptides that match multiple proteins per gene
        ### This is counting how often a peptide is found within the same protein group (e.g. matches ALB_1, ALB_2 and ALB_3 multimatch = 3)
        grouped2["Multimatch"] = grouped2.groupby(["Sequence", "Gene"])["Sequence"].transform("count")
        
        ### Reorder and get rid of unique alignment artifacts based on comma count in Proteins
        grouped2["prot_comma_count"] = grouped2["Proteins"].str.count(",")
        grouped2 = grouped2.sort_values(by = ["Sequence", "Gene", "prot_comma_count"])
        
        ### Collapse by Gene and Sequence, (so result should be one entry = 1 peptide and all of its matches together)
        PepSeqs = grouped2.groupby(["Gene", "Sequence"]).agg(
            Proteins=('Proteins', lambda x: ';'.join(x)), ### Different proteins from the group that matched, merged with ;
            AlignedPos=('AlignedPos', lambda x: x.iloc[0]), #### The first position (1 for the whole group) in the alignment that matched
            GaplessPositions=('GaplessPositions', lambda x: x.iloc[0]) ### The first position for each gapless match (1 per protein in the protein group)
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


########## Control whether or not Delta Score will be recalculated
recalculate_delta_scores = True


#### Check if 'Reverse' or 'Decoy' column exists in MSMS (older version of MaxQuant seem to have 'Reverse', newer have 'Decoy'
if "Reverse" in MSMS.columns:
    decoy_col = "Reverse"
elif "Decoy" in MSMS.columns:
    decoy_col = "Decoy"
else:
    raise ValueError("Neither 'Reverse' nor 'Decoy' column found in the MSMS table.")




###### --- Optional: Recalculate delta scores if enabled ---
############################################################
if recalculate_delta_scores:
    
    
    ### Bckup MSMS Deltascores, just in case it is needed
    MSMS["DeltaBackup"] = MSMS["Delta score"].copy()


    #### Split top sequences, modified sequences, and scores (handled as lists by splitting on ';')
    MSMS["Sequences"] = MSMS["All sequences"].astype(str).str.split(";")
    MSMS["ModSequences"] = MSMS["All modified sequences"].astype(str).str.split(";")
    
    def parse_scores(val): ### Quick function to check if "All scores" column (delta scores) can be parsed
        try:
            return [float(x) for x in str(val).split(";") if x.strip() != ""]
        except Exception:
            return []
            
    MSMS["AllScores"] = MSMS["All scores"].apply(parse_scores)

    ##### 1) Identify rank 1 and 2 hits identical except 1 Leucine / Isoleucine change, add it as a new column in MSMS
    ##### For now we cannot distinguish I and L, so if this is the only difference, then there is no difference at all
    il_hits = []
    for seqs in MSMS["Sequences"]:
        if isinstance(seqs, list) and len(seqs) >= 3 and pd.notna(seqs[0]) and pd.notna(seqs[1]):
            s1 = seqs[0].replace("I", "L")
            s2 = seqs[1].replace("I", "L")
            il_hits.append(s1 == s2)
        else:
            il_hits.append(False)
    MSMS["IL_hit"] = il_hits

    ##### Update delta score for IL hits, use the score difference between 2nd and 3rd hit instead of 1st and 2nd
    new_delta = MSMS["Delta score"].copy()
    for idx in MSMS.index[MSMS["IL_hit"]]:
        scores = MSMS.loc[idx, "AllScores"]
        if isinstance(scores, list) and len(scores) >= 3 and pd.notna(scores[1]) and pd.notna(scores[2]):
            new_delta.loc[idx] = scores[1] - scores[2]
    MSMS["Delta score"] = new_delta




    #### 2) Identify rank 1 and 2 hits identical except Q/E or N/D changes
    ##### Deamidation are so common, so Q to E or N to D changes, should not count for difference between hits
    qe_hits = []
    for mod_seqs in MSMS["ModSequences"]:
        if isinstance(mod_seqs, list) and len(mod_seqs) >= 3 and pd.notna(mod_seqs[0]) and pd.notna(mod_seqs[1]):
            s1 = re.sub(r"Q\(de\)", "E", re.sub(r"N\(de\)", "D", mod_seqs[0]))
            s2 = re.sub(r"Q\(de\)", "E", re.sub(r"N\(de\)", "D", mod_seqs[1]))
            qe_hits.append(s1 == s2)
        else:
            qe_hits.append(False)
    MSMS["QE_hit"] = qe_hits
    
    ##### Again, update delta score for Q-E/N-D hits, use the score difference between 2nd and 3rd hit instead of 1st and 2nd
    new_delta = MSMS["Delta score"].copy()
    for idx in MSMS.index[MSMS["QE_hit"]]:
        scores = MSMS.loc[idx, "AllScores"]
        if isinstance(scores, list) and len(scores) >= 3 and pd.notna(scores[1]) and pd.notna(scores[2]):
            new_delta.loc[idx] = scores[1] - scores[2]
    MSMS["Delta score"] = new_delta






###===========================================================
### Delta Score Optimization Plots & FDR Cutoffs

## Standard FDR is 0.01, unless we want to change it
FDR = 0.01

#### Sort MSMS by Sample (== Subset)
MSMS["Subset"] = MSMS["Sample"]
MSMS = MSMS.sort_values(by = ["Subset", "Score"], ascending = [True, False]).reset_index(drop = True)



###### --- Optional: Recalculate delta scores if enabled ---
############################################################

if recalculate_delta_scores:
    
    #### Prepare iterative Delta cutoff list by 1 (0 to 30), then every 10 (40, 50, ..., 100)
    min_deltas = list(range(31)) + [i * 10 for i in range(4, 11)]
    delta_opt_records = []


    ### Cycle through cutoffs
    for min_d in min_deltas:
        tmp = MSMS[MSMS["Delta score"] >= min_d].copy() ### get peptides that are above this delta score
        tmp = tmp.sort_values(by=["Subset", "Score"], ascending=[True, False]).reset_index(drop = True) ### Sort them by sample and score (important, as we are going from best score to worse)
        
        ###### Calculate q-value per subset (sample), using the name of the decoy column.
        ###### The function goes from top to bottom of the table (which is sorted by Score)
        ###### At every entry, it counts how many false negatives (decoys/reverse) have been found by that point (negatives with higher score) and divides it by the number of total hits by that point (hits with higher score)
        ###### That number is the new Qvalue to be compared with the FDR
        tmp["qValue"] = tmp.groupby("Subset", group_keys = False, observed=False).apply(lambda group: (group[decoy_col] == "+").cumsum() / np.arange(1, len(group) + 1),include_groups = False)
  
        
        #### Select all that have qvalue <= FDR
        passed_ids = set(tmp.loc[tmp["qValue"] <= FDR, "Evidence ID"])
        evidence["PassedMSMS_FDR"] = evidence["id"].isin(passed_ids)
        evidence["SampleGroup"] = evidence["Sample"].astype(str)
        
        ##### Aggregate metrics for passed items
        passed_ev = evidence[evidence["PassedMSMS_FDR"] == True]
        if not passed_ev.empty:
            agg_df = passed_ev.groupby("SampleGroup").agg(
                EvidenceHits=("id", "count"),
                PSMs=("MS/MS count", "sum"),
                Peptides=("Sequence", "nunique")
            ).reset_index()
            agg_df["MinDelta"] = min_d
            delta_opt_records.append(agg_df)

    #### if not empty, prepare plot data
    if delta_opt_records:
        DeltaOptPlotting = pd.concat(delta_opt_records, ignore_index=True)
        
        # Melt dataframe for plotting
        plot_data = DeltaOptPlotting.melt(
            id_vars=["MinDelta", "SampleGroup"],
            value_vars=["Peptides", "EvidenceHits", "PSMs"],
            var_name="Identification",
            value_name=F"IDs-at-{FDR}-FDR"
        ).rename(columns={"SampleGroup": "Sample"})
        
        # Find optimal points (max IDs at 1% FDR for each Identification and Sample)
        optimal_points = plot_data.loc[
            plot_data.groupby(["Identification", "Sample"])[F"IDs-at-{FDR}-FDR"].idxmax()
        ].copy()

        # Calculate median optimal MinDelta across metrics for each sample
        Optimal_Delta_scores = optimal_points.groupby("Sample")["MinDelta"].median().reset_index(name="Optimal_Delta")
    else:
        Optimal_Delta_scores = pd.DataFrame(columns=["Sample", "Optimal_Delta"])
else:
    Optimal_Delta_scores = pd.DataFrame(columns=["Sample", "Optimal_Delta"])




### ==========-=----=-=-=-=-=-=-=============================
### Optional Plotting of FDR Cutoff effects

## Trigger to control plotting
Plotting_is_a_go = True

if recalculate_delta_scores and not plot_data.empty and Plotting_is_a_go:
    
    ### make sure these are installed before running
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    identifications = plot_data["Identification"].unique()
    samples = plot_data["Sample"].unique()
    
    num_rows = len(identifications)
    num_cols = len(samples)
    
    # Dynamically calculate width based on number of samples
    dynamic_width = max(10, num_cols * 2.5)
    
    fig, axes = plt.subplots(
        num_rows, num_cols, 
        figsize=(dynamic_width, num_rows * 3), 
        sharey="row", 
        sharex=True
    )
    
    # Ensure axes is always a 2D array even if rows/cols == 1
    if num_rows == 1 and num_cols == 1:
        axes = np.array([[axes]])
    elif num_rows == 1:
        axes = axes[np.newaxis, :]
    elif num_cols == 1:
        axes = axes[:, np.newaxis]

    # Color palette matching identifications
    palette = sns.color_palette("tab10", num_rows)
    color_dict = {ident: palette[i] for i, ident in enumerate(identifications)}

    for r_idx, ident in enumerate(identifications):
        for c_idx, sample in enumerate(samples):
            ax = axes[r_idx, c_idx]
            subset = plot_data[(plot_data["Identification"] == ident) & (plot_data["Sample"] == sample)]
            
            if not subset.empty:
                col_color = color_dict[ident]
                # Plot line
                ax.plot(subset["MinDelta"], subset[F"IDs-at-{FDR}-FDR"], color=col_color, linewidth=1.2)
                
                # Find optimal point for this facet
                opt_row = optimal_points[(optimal_points["Identification"] == ident) & (optimal_points["Sample"] == sample)]
                if not opt_row.empty:
                    opt_delta = opt_row["MinDelta"].values[0]
                    # Vertical dashed line at optimal MinDelta
                    ax.axvline(x=opt_delta, color=col_color, linestyle="--", alpha=0.6)
                    # Text annotation at the bottom
                    ax.text(opt_delta, 0, str(int(opt_delta)), rotation=90, va="bottom", ha="right", fontsize=9)
            
            # Label columns on top row, rows on left column
            if r_idx == 0:
                ax.set_title(str(sample), fontsize=10, fontweight="bold")
            if c_idx == 0:
                ax.set_ylabel(str(ident), fontsize=10, fontweight="bold")
            if r_idx == num_rows - 1:
                ax.set_xlabel("MinDelta", fontsize=9)
            
            ax.tick_params(axis="x", rotation=45, labelsize=8)
            ax.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    
    # Save plot as PDF (equivalent to ggsave)
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, "DeltaScore-Opt.pdf"), bbox_inches="tight", dpi=300)
    plt.close()


###==========================================================
### Apply PSM FDR Control

###### --- Optional: Recalculate delta scores if enabled ---
############################################################
if recalculate_delta_scores and not Optimal_Delta_scores.empty:
    MSMS_F = MSMS.copy()
    MSMS_F["Subset"] = MSMS_F["Sample"]
    
    # Merge optimal delta scores
    MSMS_F = pd.merge(MSMS_F, Optimal_Delta_scores, on="Sample", how="left")
    MSMS_F["Optimal_Delta"] = MSMS_F["Optimal_Delta"].fillna(15)
    
    # Filter by threshold
    MSMS_F = MSMS_F[MSMS_F["Delta score"] >= MSMS_F["Optimal_Delta"]]
    MSMS_F = MSMS_F.drop(columns=["Optimal_Delta"])
else:
    MSMS_F = MSMS.copy()




#### Sort and compute final q-value
MSMS_F = MSMS_F.sort_values(by=["Subset", "Score"], ascending=[True, False]).reset_index(drop=True)
MSMS_F["qValue"] = MSMS_F.groupby("Subset", group_keys = False, observed = False).apply(lambda group: (group[decoy_col] == "+").cumsum() / np.arange(1, len(group) + 1), include_groups = False) #### See previous identical function for explanation a bit above


#### Mini Function Select PSMs under q-value threshold
def filter_passed_fdr(group, FDR = 0.01):
    valid_idx = group["qValue"] <= FDR
    if not valid_idx.any():
        return pd.Series([False] * len(group), index=group.index)
    last_valid = np.where(valid_idx)[0][-1]
    
    return pd.Series(np.arange(len(group)) <= last_valid, index=group.index)



MSMS_F["PassedFDR"] = MSMS_F.groupby("Subset", group_keys = False, observed = False).apply(filter_passed_fdr, FDR = 0.01, include_groups = False)
MSMS_F = MSMS_F[MSMS_F["PassedFDR"] == True].copy()














# ==========================================================
### 5. Transfer the filtered sequences from the MSMS table to the evidence table and
### collapse to site level to combine information from overlapping peptides on 1 site
# ==========================================================












# ==========================================================
## 6. Generate the consensus sequence
# ==========================================================













# ==========================================================
## 7. Generate output Files
# ==========================================================
