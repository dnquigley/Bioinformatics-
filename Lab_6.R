library(Biostrings)
library(seqinr)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")

library(msa)

setwd("/Users/devinquigley/Desktop/Bioinformatics/Bioinformatics-")

KU18 <- readDNAStringSet("VAME_Lab6/VAME_KU18.fasta")
KU18

KU15 <- readDNAStringSet("VAME_Lab6/VAME_KU15.fasta")
KU15

KU14 <- readDNAStringSet("VAME_Lab6/VAME_KU14.fasta")
KU14

KU13 <- readDNAStringSet("VAME_Lab6/VAME_KU13.fasta")
KU13

KU16 <- readDNAStringSet("VAME_Lab6/VAME_KU16.fasta")
KU16











