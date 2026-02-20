library(Biostrings)
library(seqinr)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
BiocManager::install("muscle")

library(msa)

setwd("/Users/devinquigley/Desktop/Bioinformatics/Bioinformatics-")

#Read in DNA fasta files
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

#Combine fasta files
sequences <-c(KU18, KU15, KU14, KU13, KU16)
sequences

#Rename sequences in combined file
names(sequences) <- c("KU18", "KU15", "KU14", "KU13", "KU16")
sequences

#Options to align sequences 
alignment <- msa(sequences, method = "Muscle")
alignment

alignment2 <- msa(sequences, method = "ClustalW")
alignment2

print(alignment, show="complete")

#Finding gaps in the alignment
install.packages("ape")
library(ape)



#Redo alignment and gap check
alignment_dnabin <- msaConvert(alignment, type = "ape::DNAbin")

gap.stats <- gap.inspect(alignment_dnabin)

gap.stats

aligned_stringset <- as(alignment, "DNAStringSet")
alignment_matrix <- as.matrix(aligned_stringset)

sum(alignment_matrix == "-")

gap_columns <- which(colSums(alignment_matrix == "-") > 0)
gap_columns



#extract variable sites to check for substitutions

aligned_stringset <- as(alignment, "DNAStringSet")
alignment_matrix <- as.matrix(aligned_stringset)

is_variable <- function(column) {
  length(unique(column)) > 1
}

variable_sites <- which(apply(alignment_matrix, 2, is_variable))

variable_sites
length(variable_sites)

#Length of alignment

# If you have a DNAStringSet
aligned_stringset <- as(alignment, "DNAStringSet")
width(aligned_stringset)  # Returns the length of each sequence (all should be equal in an alignment)



#Finding GC content

library(Biostrings)

#Convert to string set to do count
aligned_stringset <- as(alignment, "DNAStringSet")

# Count bases for each sequence
counts <- alphabetFrequency(aligned_stringset)
counts

#Calculate GC content for all sequences
base_counts <- counts[, c("A","C","G","T")]
GC_content_all <- rowSums(base_counts[, c("G","C")]) / rowSums(base_counts) * 100
GC_content_all

#Convert to seqinr
seqinr_alignment <- msaConvert(alignment)

#Distance matrix
dist_matrix <- dist.alignment(seqinr_alignment, "identity")
dist_matrix

as.matrix(dist_matrix)



#Translate to amino acid
library(Biostrings)

# Extract the first sequence from the alignment
library(seqinr)

KU18_seq <- KU18[[1]]

# Convert DNAString to a plain character string
KU18_string <- as.character(KU18_seq)  

# Now convert the string to characters for seqinr
KU18_chars <- s2c(KU18_string)  

# Translate to amino acids
KU18_protein <- translate(KU18_chars)

# View as a string
paste(KU18_protein, collapse = "")



#Write alignment to a file
install.packages("phangorn")
library(phangorn)

# Convert MSA to phyDat object
Alignment_phyDat <- msaConvert(alignment, type = "phangorn::phyDat")

# Write to file in FASTA format
write.phyDat(Alignment_phyDat, file = "alignment.fasta", format = "fasta")





