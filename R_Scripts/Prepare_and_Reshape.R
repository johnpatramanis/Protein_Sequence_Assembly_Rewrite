##################################################
# Project: MaxQuant Analysis Pipeline
# Author(s): Patrick Ruether (Generalized)
# Script purpose: Environment setup, data loading, and FASTA processing
# Code Updates by Ioannis Patramanis in the year of our lord 2026
##################################################

### Load Libraries
library(data.table)
library(stringr)
library(stringi)
library(seqinr)
library(bit64)


###### Paths and Environment setup
### Make paths easily adjustable or pass them via arguments/config
### Safe working directory setting with command-line fallback
script_dir <- NULL
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    script_dir <- getwd()
  }
}

# Set working directory to the script's location
setwd(script_dir)

# Source Functions.R located in the same path as the script itself
source("Functions.R")




############## Prepare / Load Files ############################  

### 0. Setting up paths and loading files from MaxQuant

## Set paths as variables so the script can be easily adapted to any MaxQuant run folder and FASTA file
maxquant_dir <- "../Paranthropus_Yiannis_Reanalysis/txt"
fasta_path <- "../Test_Data_Input/PR220121_GreatApes_Enamelome_aligned.fasta"
exclusion_path <- "../Paranthropus_Yiannis_Reanalysis/IP_Exclusion_List.csv"
samples_to_experiment_path <- "../Paranthropus_Yiannis_Reanalysis/Sample_to_Experiment.tsv"
output_dir <- "../Paranthropus_Yiannis_Reanalysis"




if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
## LOAD TABLES

## 1) Load MaxQuant Evidence File
evidence <- fread(file.path(maxquant_dir, "evidence.txt"))

## 2) Load MaxQuant Summary File
Summary <- fread(file.path(maxquant_dir, "summary.txt"))
Summary <- Summary[Experiment != ""]

## 3) Load MaxQuant MSMS File
MSMS <- fread(file.path(maxquant_dir, "msms.txt"))


## 4) Load FASTA
#### Information from Fasta labels are loaded as columns
FASTA <- fasta_to_datatable(paths = fasta_path)

## 5) Generate gapless sequences, add a new column, which is just the same sequence, but without gaps
FASTA[, Gapless := str_remove_all(Sequence, "-")]

## Add identifier (number 1 to N), add new column that is just and ID per sequence
FASTA[,id := 1:.N]

## Create new table
## Generate site level fasta site by cutting the sequence after each a.a. and making it a new row in a new table
FASTAsites <- FASTA[, .(AA = unlist(str_split(Sequence, "")),
                        AlignedPos = 1:nchar(Sequence)),
                    by = .(Protein, Annot, Gene, Species, id)]

## Connect gapless and normal sequence position
## Add gapless position (add +1 for all "_")
FASTAsites[, GaplessPos := cumsum(AA != "-"), by = .(Annot)]

## 
## Add identifier for each amino acid: gene_aligned position
FASTAsites[, GeneSite := paste(Gene, AlignedPos, sep = "_")]






#################### Reshape and Output Files #############################



#####################################################
### 1. Experimental annotation from summary table 

#### Create new data table from Summary that has Experiment, Raw file, and Fraction columns
ExpAnnot <- Summary[nzchar(`Enzyme mode`), .(Experiment, `Raw file`, Fraction)]
setorder(ExpAnnot, `Raw file`)


## load samples to experiments tab separated file, if it exists
if (file.exists(samples_to_experiment_path)) {
  samples_to_experiment <- fread(samples_to_experiment_path)
  
  # Ensure standard column names (1st col = Sample, 2nd col = Experiment)
  setnames(samples_to_experiment, c("Sample", "Experiment"))
  
  # Merge mapping into ExpAnnot by Experiment
  ExpAnnot <- merge(ExpAnnot, samples_to_experiment, by = "Experiment", all.x = TRUE)
  
  # Fallback for any experiments missing from the mapping file
  ExpAnnot[is.na(Sample), Sample := Experiment]
} else {
  ## Fallback if no mapping file is provided
  ExpAnnot[, Sample := Experiment]
}
## Cleanup, create new column named "Sample"



## Convert Sample Column to Factor
ExpAnnot[, Sample := as.factor(Sample)]

## Add info whether or not they are fractionated or single shot
# Based on whether an Experiment has multiple entries or a non-empty/non-zero Fraction value
ExpAnnot[, n_entries := .N, by = Experiment]
ExpAnnot[, Fractionated := as.factor(n_entries > 1 | (is.numeric(Fraction) & Fraction > 1))]
levels(ExpAnnot$Fractionated) <- c("single", "frc")

### For fractionated samples, add number of the fraction using the Summary Fraction column
ExpAnnot[Fractionated == "frc", Fraction := as.integer(Fraction)]
ExpAnnot[, n_entries := NULL] ## Clean up temporary helper column
ExpAnnot[, Experiment := as.factor(Experiment)]







#########################################
##### 2. Annotate tables

#### Add new info to MSMS table and evidence table
MSMS <- merge(MSMS, ExpAnnot, by = "Raw file")
evidence <- merge(evidence[, !c("Experiment", "Fraction")], ExpAnnot, by = c("Raw file"))








#########################################
## 3. Annotate Peptide sequences
## Here the goal is to create a csv file with 7 columns: Gene(gene name), Peptide(1 short aa sequence), Proteins(A collection of protein sequences that correspond to that gene), AlignedPos(for each alignment of proteins sequences, the positions of that match, 1 sequence for the entire alignemnt), GaplessPositions(for each one protein sequence, which position in the gapless version it matches), Multiplicity(if it is matching multiple genes), AnyMultiPep()
#### if PepSeq File Exists, read it

pepseq_files <- list.files(output_dir, pattern = "PepSeq")

if(length(pepseq_files) > 0 && any(grepl("PepSeq", pepseq_files))) {
  target_file <- ifelse("PepSeqs.csv" %in% pepseq_files, "PepSeqs.csv", pepseq_files[1])
  PepSeqs <- fread(file.path(output_dir, target_file))
} else { #### Else we have to make it
  
  ### Get all unique peptide sequence from MSMS.txt -> list of unique peptide sequences
  PepSeqs <- unique(MSMS[, .(Sequence)])
  
  ### Find all proteins in the fasta (gapless sequences) that contain this peptide -> add second column with all the proteins containing that sequence
  PepSeqs[, Proteins := lapply(Sequence, function(x) FASTA[stri_detect(Gapless,  regex = x), ]$Protein)]
  
  ### Take only peptides found in the FASTA and Make 1 row for each PepSequence-Protein pair (= expand table to long format)
  PepSeqs <- PepSeqs[, .(Protein = unlist(Proteins)), by = .(Sequence)]
  
  ### Locate the peptide sequence in the respective FASTA gapless sequence -> add start and end position of the peptide from the gapless sequence
  PepSeqs[, GaplessPos := mapply(function(x, y) str_locate_all(FASTA[match(x, Protein)]$Gapless, y), x = Protein, y = Sequence)]
  
  ### Only keep the start position of each peptide match
  PepSeqs[, GaplessPos := lapply(GaplessPos, function(x) x[,1])]
  
  ### Repeat peptide with multiple matches in the same protein sequence (probably irrelevant for enamel)
  PepSeqs <- PepSeqs[, .(GaplessPos = unlist(GaplessPos)), by = .(Sequence, Protein)]
  
  ### Add an identifier before conversion to site level
  PepSeqs[, id := 1:.N]
  
  #### Transform to site level for mapping aligned positions -> expand table: 1 line per amino acid for each peptide-prot sequence
  PepSeqs <- PepSeqs[, .(GaplessPos = GaplessPos:(GaplessPos+nchar(Sequence)-1)), by = .(Sequence, Protein, id)]
  
  #### Add aligned positions and gene
  PepSeqs <- merge(PepSeqs, FASTAsites[AA != "-", .(Gene, Protein, AlignedPos, GaplessPos)], by = c("Protein", "GaplessPos"), all.x = T)
  
  #### Sanity check: Any peptides without aligned position?
  PepSeqs[is.na(AlignedPos)]
  
  ##### Go back to peptide level
  PepSeqs <- PepSeqs[, .(GaplessPos = paste(GaplessPos, collapse = ","),
                         AlignedPos = paste(AlignedPos, collapse = ",")),
                     by = .(Gene, Protein, Sequence, id)]
  
  ##### Filter down to unique peptide sequence and Gene
  PepSeqs <- PepSeqs[, .(Proteins = paste(Protein, collapse = ";"),
                         GaplessPositions = paste(GaplessPos, collapse = ";")),
                     .(Sequence, Gene, AlignedPos)]
  
  ##### Check for peptides that match multiple proteins per gene
  PepSeqs[, Multimatch := .N, by = .(Sequence, Gene)]
  
  ##### Reorder things, get rid of unique weird alignment artifact
  PepSeqs <- PepSeqs[order(Sequence, Gene, str_count(Proteins, ","))]
  PepSeqs <- PepSeqs[, .(Proteins = paste(Proteins, collapse = ";"),
                         AlignedPos = AlignedPos[1],
                         GaplessPositions = GaplessPositions[1]), by = .(Gene, Sequence)]
  
  
  ##### Calculate how often a peptide sequence is repeated in the table (generalized from == 2 to > 1)
  PepSeqs[, Multiplicity := .N, by = .(Sequence)]
  PepSeqs[Multiplicity > 1, .(Genes = paste(Gene, collapse = ";")), by = .(Sequence, Multiplicity)]
  
  ##### Apply razor principle (winner takes all) to peptides mapping to multiple genes
  PepSeqs[, AnyMultiPep := any(Multiplicity > 1), by = .(Gene)]
  
  ###### Subset only peptides matching multiple genes and count how many peptides per gene in the entire dataset
  MultiPeps <- merge(PepSeqs[Multiplicity > 1], PepSeqs[AnyMultiPep == T, .(PepCount = .N), by = .(Gene)],
                     by = "Gene")
  
  ###### Find and keep "razor" gene with most peptides in the entire dataset
  MultiPeps[, MaxPepCount := max(PepCount), by = .(Sequence)]
  MultiPeps <- MultiPeps[PepCount == MaxPepCount]
  MultiPeps[, .N, by = .(Sequence)][N > 1]
  MultiPeps[, .N, by = .(Gene)]
  
  #### Combine the peptide sequences with unique gene mapping
  PepSeqs <- rbind(PepSeqs[Multiplicity == 1], MultiPeps[, !c("MaxPepCount", "PepCount")])
  
  #### Export site table as csv into the output folder
  fwrite(PepSeqs, file.path(output_dir, "PepSeqs.csv"))
}








##################################################################################
##### 4. Delta Score recalculation and cutoff optimization

#### Recalculate delta scores
MSMS[, DeltaBackup := copy(`Delta score`)]

### Split top 3 sequences from 'all sequences' tab of MQ
MSMS[, Sequences := str_split(`All sequences`, ";")]

### Split top 3 modified sequences
MSMS[, ModSequences := str_split(`All modified sequences`, ";")]

#### Split top 3 scores
MSMS[, AllScores := str_split(`All scores`, ";")]
MSMS[, AllScores := lapply(AllScores, as.numeric)]

#### Identify rank 1 and 2 hits that are identical except 1 leucine/isoleucine change (with safety check for vector length)
MSMS[, IL_hit := unlist(lapply(Sequences, function(x) {
  if(length(x) >= 3 && !is.na(x[1]) && !is.na(x[2])) {
    str_replace_all(x[1], "I", "L") == str_replace_all(x[2], "I", "L")
  } else {
    FALSE
  }
}))]

MSMS[IL_hit == T, `Delta score` := mapply(function(scores, current_ds) {
  if(length(scores) >= 3 && !is.na(scores[2]) && !is.na(scores[3])) {
    scores[2] - scores[3]
  } else {
    current_ds
  }
}, AllScores[IL_hit], `Delta score`[IL_hit])]

##### Identify rank 1 and 2 hits that are identical except 1 Q/E or N/D (with safety check for vector length)
MSMS[, QE_hit := unlist(lapply(ModSequences, function(x) {
  if(length(x) >= 3 && !is.na(x[1]) && !is.na(x[2])) {
    str_replace_all(x[1], c("Q\\(de\\)" = "E", "N\\(de\\)" = "D")) == str_replace_all(x[2], c("Q\\(de\\)" = "E", "N\\(de\\)" = "D"))
  } else {
    FALSE
  }
}))]

MSMS[QE_hit == T, `Delta score` := mapply(function(scores, current_ds) {
  if(length(scores) >= 3 && !is.na(scores[2]) && !is.na(scores[3])) {
    scores[2] - scores[3]
  } else {
    current_ds
  }
}, AllScores[QE_hit], `Delta score`[QE_hit])]







##### Delta Score optimization plots - Sample-wise FDR at different delta score cutoffs - peptides, evidence hits, PSMs



FDR <- 0.01
MSMS[, Subset := Sample]
setorder(MSMS, Subset, -Score)



### Check if 'Reverse' or 'Decoy' column exist

decoy_col <- if ("Reverse" %in% names(MSMS)) {
  "Reverse"
} else if ("Decoy" %in% names(MSMS)) {
  "Decoy"
} else {
  stop("Neither 'Reverse' nor 'Decoy' column found in the MSMS table.")
}


### Prepare iterative Delta cutoff list
DeltaOpt <- data.table(MinDelta = c(0:30, 10 * (4:10)))


#### Go over each Delta cutoff
DeltaOpt[, IDs := lapply(MinDelta, function(x) {
  tmp <- MSMS[`Delta score` >= x]
  setorder(tmp, Subset, -Score)
  
  # Calculate q-value using whichever decoy column is available
  if (decoy_col == "Reverse") {
    tmp[, qValue := cumsum(Reverse == "+")/1:.N, by = .(Subset)]
  } else {
    tmp[, qValue := cumsum(Decoy == "+")/1:.N, by = .(Subset)]
  }
  
  #### Then 
  evidence[, .(EvidenceHits = .N,
               PSMs = sum(`MS/MS count`),
               Peptides = length(unique(Sequence)),
               Sequences = list(unique(Sequence))), 
           keyby = .(PassedMSMS_FDR = id %in% tmp[qValue <= FDR]$`Evidence ID`,
                     SampleGroup = as.character(Sample))][PassedMSMS_FDR == T]
})]


DeltaOpt[, Samples := lapply(IDs, function(x) x$SampleGroup)]
DeltaOpt[, Peptides := lapply(IDs, function(x) x$Peptides)]
DeltaOpt[, EvidenceHits := lapply(IDs, function(x) x$EvidenceHits)]
DeltaOpt[, PSMs := lapply(IDs, function(x) x$PSMs)]
DeltaOptPlotting <- DeltaOpt[, .(Sample = unlist(Samples),
                                 Peptides = unlist(Peptides),
                                 EvidenceHits = unlist(EvidenceHits),
                                 PSMs = unlist(PSMs)),
                             by = .(MinDelta)]

# Prepare melted data and calculate optimal MinDelta (max IDs) for each facet
plot_data <- melt(DeltaOptPlotting[!is.na(Sample), ], measure.vars = c("Peptides", "EvidenceHits", "PSMs"),
                  variable.name = "Identification", value.name = "IDs-at-1%-FDR")

optimal_points <- plot_data[, .SD[which.max(`IDs-at-1%-FDR`)], by = .(Identification, Sample)]

# Dynamically calculate width to ensure subplots are doubled in width per sample
num_samples <- length(unique(plot_data$Sample))
dynamic_width <- max(20, num_samples * 8)



######### Plotting ############

library(ggplot2)


ggplot(plot_data, aes(x = MinDelta, y = `IDs-at-1%-FDR`, color = Identification)) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = optimal_points, aes(xintercept = MinDelta, color = Identification), linetype = "dashed", alpha = 0.6) +
  geom_text(data = optimal_points, aes(x = MinDelta, y = 0, label = MinDelta), 
            angle = 90, hjust = -0.1, vjust = -0.5, size = 3.5, show.legend = FALSE) +
  facet_grid(Identification ~ Sample, scales = "free") +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.spacing = unit(1, "lines"))

ggsave(file.path(output_dir, "DeltaScore-Opt.pdf"),
       width = dynamic_width, height = 10, limitsize = FALSE)


################################



## Calculate the median optimal MinDelta across the three identification metrics for each sample
Optimal_Delta_scores <- optimal_points[, .(Optimal_Delta = median(MinDelta)), by = Sample]


##### Apply PSM FDR control with optimized delta score cutoff

FDR <- 0.01

# Copy MSMS and set Subset
MSMS_F <- copy(MSMS)
MSMS_F[, Subset := Sample]

# Merge with Optimal_Delta_scores to get sample-specific thresholds
MSMS_F <- merge(MSMS_F, Optimal_Delta_scores, by.x = "Sample", by.y = "Sample", all.x = TRUE)

# Fallback: if any sample lacks an optimal delta score, default to 15
MSMS_F[is.na(Optimal_Delta), Optimal_Delta := 15]

# Filter using the sample-specific Delta score cutoff
MSMS_F <- MSMS_F[`Delta score` >= Optimal_Delta]

# Clean up the temporary column
MSMS_F[, Optimal_Delta := NULL]




##### Order the MSMS table by 1. sample 2. decreasing score
setorder(MSMS_F, Subset, -Score)

#### Calculate q-Value
if (decoy_col == "Reverse") {
MSMS_F <- MSMS_F[, qValue := cumsum(Reverse == "+")/1:.N, by = .(Subset)]
} else{
  
MSMS_F <- MSMS_F[, qValue := cumsum(Decoy == "+")/1:.N, by = .(Subset)]
  
}
#### Select the last PSM that's under the q-value threshold
MSMS_F[, PassedFDR := 1:.N <= which.max(which(qValue <= FDR)), by = .(Subset)]
MSMS_F[, .N, by = .(Subset, PassedFDR)]
MSMS_F <- MSMS_F[PassedFDR == T]





####################################################################################
### 5. Transfer the filtered sequences from the MSMS table to the evidence table and
### collapse to site level to combine information from overlapping peptides on 1 site

#### Filter the evidence table based on re-calculated delta scores and FDR
evidenceF <- evidence[id %in% MSMS_F$`Evidence ID`]

##### FDR from the 1% applied on the MSMS table. Still need to remove the 1% reverse hits. Remove matches to Contaminants
if (decoy_col == "Reverse") {
evidenceF <- evidenceF[Reverse != "+"]
} else{
evidenceF <- evidenceF[Decoy != "+"]
}
evidenceF <- evidenceF[!all.contain(Proteins, pattern = "CON"),]

#### Annotate evidence with pep locations and filter based on exclusion list
evidenceF <- merge(evidenceF[, !c("Proteins")], 
                   unique(PepSeqs[, .(Sequence, Gene, Proteins, AlignedPos, GaplessPositions)]), 
                   "Sequence", all.x = T, allow.cartesian = T)

##### Facultative - if exclusion list is available (made conditional for general datasets)

if (file.exists(exclusion_path)) {
  ExclusionList <- fread(exclusion_path)
  ExclusionList[, Excl := paste0(Sequence, "_", Sample)]
  evidenceF[, Excl := paste0(Sequence, "_", Sample)]
  evidenceF <- evidenceF[!Excl %in% ExclusionList$Excl]
  evidenceF[, Excl := NULL] # Clean up temporary column
}

#### Export filtered evidence into the output folder
fwrite(evidenceF, file.path(output_dir, "EvidenceFiltered.csv"))





##### IDs - generate table with PSMs, mod peptides and peptides count per experiment (ie. sample + frac or SS)
IDbyExperiment <- evidenceF[, .(PSMs = sum(`MS/MS count`),
                                ModPeps = length(unique(`Modified sequence`)),
                                Peptides = length(unique(Sequence))),
                            keyby = .(Sample, Fractionated, Experiment)]


#### IDs to long - generate long table with PSMs, mod peptides and peptides count per experiment (ie. sample + frac or SS)
IDbyExperimentLong <- melt(IDbyExperiment, 
                           measure.vars = c("PSMs", "ModPeps", "Peptides"), 
                           variable.name = "Identifications", 
                           value.name = "IDs")

#### IDs by sample - generate table with PSMs, mod peptides and peptides count per sample
IDbySample <- evidenceF[, .(PSMs = sum(`MS/MS count`),
                            ModPeps = length(unique(`Modified sequence`)),
                            Peptides = length(unique(Sequence))),
                        keyby = .(Sample)]
#### IDs to long - generate long table with PSMs, mod peptides and peptides count per sample
IDbySampleLong <- melt(IDbySample, 
                       measure.vars = c("PSMs", "ModPeps", "Peptides"), 
                       variable.name = "Identifications", 
                       value.name = "IDs")

#### Export tables into the output directory
fwrite(IDbyExperiment, file.path(output_dir, "IDbyExp.csv"))
fwrite(IDbyExperimentLong, file.path(output_dir, "IDbyExpLong.csv"))
fwrite(IDbySample, file.path(output_dir, "IDbySample.csv"))
fwrite(IDbySampleLong, file.path(output_dir, "IDbySampleLong.csv"))




##### Collapse unique modified peptides = Collapse charge state
evidenceF[is.na(Intensity), Intensity := 0]
modpepsUnique <- evidenceF[!is.na(Gene), 
                           .(PSMcount = sum(`MS/MS count`),
                             Charges = paste(unique(Charge), collapse = "/"),
                             MedScore = median(Score, na.rm = T),
                             MaxScore = max(Score, na.rm = T),
                             AvgMass = mean(`m/z`),
                             AvgMassErrorPPM = mean(`Mass error [ppm]`),
                             MaxLogInt = max(log10(Intensity), na.rm = T)),
                           by = .(Sequence, Proteins, AlignedPos, `Modified sequence`, Gene, 
                                  Sample, Fractionated, Experiment)]
modpepsUnique[MaxLogInt < 0, MaxLogInt := 0]


###### Expand to site level - add identifier
modpepsUnique[, id := 1:.N]

# For (i) all peptides mapped to the alignment, split the peptide sequence into AAs and split the aligned positions
sites <- modpepsUnique[!is.na(AlignedPos), 
                       .(AA = unlist(str_split(Sequence, "")),
                         Site = as.integer(unlist(str_split(AlignedPos, ",")))),
                       by = .(Proteins, Gene, MedScore, MaxScore, AvgMassErrorPPM, MaxLogInt, 
                              PSMcount, id, Sequence, `Modified sequence`, Sample, Fractionated, Experiment)]

##### Annotate individual sites
sites[, GeneSite := paste(Gene, Site, sep = "_")]

##### Calculate sequence coverage per raw file
AACoverageRawFile <- sites[, .(AbsCoverage = length(unique(GeneSite))), by = .(Experiment)]
fwrite(AACoverageRawFile, file.path(output_dir, "AA_Coverage_by_Experiment.csv"))





##### Collapse sites to have unique unmodified amino acids
sites <- sites[, .(MedScore = median(MedScore, na.rm = T),
                   MaxScore = max(MaxScore, na.rm = T),
                   AvgMassErrorPPM = mean(AvgMassErrorPPM, na.rm = T),
                   MaxLogInt = max(MaxLogInt, na.rm = T),
                   PSMcount = sum(PSMcount, na.rm = T),
                   Sequences = paste(unique(Sequence), collapse = ";")),
               by = .(Gene, Site, GeneSite, AA, Sample)]

##### Absolute sequence coverage per sample
AACoverageSample <- sites[, .(AbsCoverage = length(unique(GeneSite))), by = .(Sample)]
fwrite(AACoverageSample, file.path(output_dir, "AA_Coverage_by_Sample.csv"))








##################################################
### 6. Generate the consensus sequence

#### Calculate relative variant coverage
sites[, SiteCoverage := sum(PSMcount), by = .(GeneSite, Sample)]
sites[, SiteInt := sum(MaxLogInt), by = .(GeneSite, Sample)]

fwrite(sites, file.path(output_dir, "Site table.csv"))

#### Consensus based on most common variant
sites[, RelCoverage := PSMcount/SiteCoverage, by = .(GeneSite, Sample)]
sites[, RelInt := MaxLogInt/SiteInt, by = .(GeneSite, Sample)]
sites[, MaxCoverage := max(RelCoverage), by = .(GeneSite, Sample)]
sites[, Multiplicity := length(unique(AA[RelCoverage > 0.1 & RelInt > 0.1 & PSMcount > 1])), by = .(GeneSite, Sample)]
VarSites <- sites[Multiplicity > 1, .(Sample, GeneSite, RelCoverage, RelInt, 
                                      PSMcount, SiteCoverage, AA)][order(Sample, GeneSite, RelCoverage, RelInt, 
                                                                         PSMcount, SiteCoverage, AA)]
#### Export table
fwrite(VarSites, file.path(output_dir, "VarSites.csv"))

#### Generate consensus table with the most common variant
consensus <- sites[RelCoverage == MaxCoverage]
consensus[, .N, by = .(GeneSite, Sample)][N > 1]
consensus[, MaxSiteInt := max(MaxLogInt), by = .(Sample, GeneSite)]
consensus <- consensus[MaxLogInt == MaxSiteInt]
consensus[, .N, by = .(GeneSite, Sample)][N > 1]

### All covered position per sample
AllPos <- copy(FASTA)
AllPos <- AllPos[Gene %in% sites$Gene, .(Sequence = Sequence[1]), by = .(Gene)]
AllPos <- AllPos[, .(Site = 1:nchar(Sequence)), by = .(Gene)]
AllPos <- AllPos[, .(Sample = unique(sites$Sample)), by = .(Gene, Site)]

##### Combine the identified sites/AAs with all possible alignment positions
consensus <- merge(consensus, AllPos, by = c("Sample", "Gene", "Site"), all = T)

##### Fill up gaps with dashes
consensus[is.na(AA), AA := "-"]
consensus[, .N, by = .(Site, Gene, Sample)][N > 1]


##### Sort for sequence assembly
setorder(consensus, Sample, Gene, Site)

##### Paste together the sites into a consensus sequence
consensus <- consensus[, .(Consensus = paste(AA, collapse = "")), by = .(Sample, Gene)]
consensus[, AbsCoverage := str_count(Consensus, "[^-]")]
consensus[, RelCoverage := AbsCoverage/nchar(Consensus)]
consensus[, Annot := paste(Sample, Gene, AbsCoverage, sep = "_")]
consensus <- merge(consensus, unique(ExpAnnot[, .(Sample)]), by = "Sample")




##################################################
####### 7. Generate output

#### Add reference sequences - 000 --> add reference sequences with sample name "000"
consensus <- rbind(consensus[, .(Consensus, Annot, Gene, Sample = as.character(Sample))], 
                   FASTA[Gene %in% consensus$Gene, .(Consensus = Sequence, Annot, Gene, Sample = "000")])
sort(c("000", as.character(unique(sites$Sample))))
setorder(consensus, Gene, Sample)

#### Export as FASTA Files into the output directory
#### Option 1: Consensus sequences for all the samples aligned with the reference sequences
write.fasta(as.list(consensus[]$Consensus),
            consensus[]$Annot,
            file.path(output_dir, paste0("PR", format(Sys.Date(), "%y%m%d"), "_", "Consensus_Ref_MQ.fasta")), as.string = T)

#### Option 2: Consensus sequences of all samples only (without reference "000")
write.fasta(as.list(consensus[Sample != "000"]$Consensus),
            consensus[Sample != "000"]$Annot,
            file.path(output_dir, paste0("PR", format(Sys.Date(), "%y%m%d"), "_", "Consensus_MQ.fasta")), as.string = T)

