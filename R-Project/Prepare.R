##################################################
# Project:
# Author(s): Patrick Ruether
# Date: 03.03.2020
# Version:
# Script purpose: Function collection
##
##
##################################################

##########################
#  Date: 190404
#  User: PR
#  Purpose: Load files
##########################

library(data.table)
library(stringr)
library(bit64)
library(seqinr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
source("R-Project/Functions.R")

##LOAD TABLES

## 1. Load Evidence
evidence <- fread("./Test_Data_Input/txt-ParanthropusSearch/evidence.txt")

## 2. Load Summary
Summary <- fread("./Test_Data_Input/txt-ParanthropusSearch/summary.txt")
Summary <- Summary[Experiment != ""]

## 3. Load MSMS 
MSMS <- fread("./Test_Data_Input/txt-ParanthropusSearch/msms.txt")


## LOAD FASTA
#### Infor from Fasta labels are loaded as columns
FASTA <- fasta_to_datatable(paths = "./Test_Data_Input/PR220121_GreatApes_Enamelome_aligned.fasta")

## Generate gapless sequences, add a new column, which is just the same sequence, but without gaps
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

