#Author: Rio Kashinmoto
#Last edit: April 2nd, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description: DADA2 pipeline from filtered reads to initial Phyloseq object.

# ---- 1. Load Required Packages ----
library(dada2)
library(phyloseq)
library(Biostrings)
library(tidyverse)

# ---- 2. Define Paths & Sample Names ----
# Set path to the cutadapt filtered fastq files
filteredReadsPath <- "path"
taxaDBPath <- "path"
