#Devin Quigley

#library packages
library(Biostrings)
library(msa)
library(seqinr)
library(pwalign)
library(ape)

#set the working directory for my file
setwd("/Users/devinquigley/Desktop/Bioinformatics/Bioinformatics-/Midterm1")

#1. Read in DNA fasta file
seqs<-readDNAStringSet("sequences.fasta")
seqs

#1a. align sequences using muscle
alignment <- msa(seqs, method = "Muscle")
alignment



#2 measure how good my alignment is

#2a. Convert alignment to matrix
aln_matrix <- as.matrix(alignment)

#2a. Calculate conservation per column
conservation <- apply(aln_matrix, 2, function(col) {
  max(table(col)) / length(col)
})

#2a. Average conservation across alignment
mean(conservation)

#2 ANSWER: The value for my average conservation across the alignment was 0.9992212
#Continued: This suggests that the alignment is high quality. Values closer to 1 suggest strong
#similarity, whereas values closer to zero, (<0.85) suggest a lot of variation. 
#Only a small number of vartion was observed, which is expected when sequencing the same gene within a population. 
#http://thegrantlab.org/bio3d/reference/conserv.html#:~:text=Score%20Residue%20Conservation%20At%20Each,a%20matrix%20to%20score%20conservation.


#3 calculate the consensus sequence for the alignment
consString <- msaConsensusSequence(alignment, type="Biostrings")
print(consString)

#3a. consensus sequence
#"AACTCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATAAAAGTCAGG
#GCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCT
#CAAACAGACACCATGGTGCATCTGACTCCTGAGGAGAACTCTGCCGTTACTGCCCTG
#TGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCA
#AGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAA
#GACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCAC
#CCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGA
#TCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAA
#AGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTT
#TGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAG
#GGTGAGTCTATGGGACGCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAATTC
#ATGTCATAGGAAGGGG"



#4. finding GC content

#4a. Convert to string set to do count
aligned_stringset <- as(alignment, "DNAStringSet")

#4a. Count bases for each sequence
counts <- alphabetFrequency(aligned_stringset)
counts

#4a. Calculate GC content for all sequences
base_counts <- counts[, c("A","C","G","T")]

#4b. Calculate overall GC content across all sequences
total_GC <- sum(base_counts[, "G"]) + sum(base_counts[, "C"])
total_bases <- sum(rowSums(base_counts))
GC_content_overall <- total_GC / total_bases * 100

GC_content_overall

#4. ANSWER GC content value: 51.56944%. This is typical for vertebrates, 
#and does not suggest anything unusual



#5 Difference check between samples 

#5a. Gap check

#5a. Gap check without gap.inspect()
aligned_stringset <- as(alignment, "DNAStringSet")
alignment_matrix <- as.matrix(aligned_stringset)

#5a. Total number of gaps
total_gaps <- sum(alignment_matrix == "-")
total_gaps

#Total gaps=1

#5b. Which columns have gaps
gap_columns <- which(colSums(alignment_matrix == "-") > 0)
gap_columns

#Gap columns=3


#5c. extract variable sites to check for substitutions

is_variable <- function(column) {
  length(unique(column)) > 1
}

variable_sites <- which(apply(alignment_matrix, 2, is_variable))

variable_sites
length(variable_sites)

#5d. Table showing the nucleotides at variable positions
alignment_matrix[, variable_sites]

#5 ANSWER: Out of the 20 sequences, 9 positions were variable (positions 3, 39, 45, 47, 134, 145, 152, 586, and 623). 
#Position 3 contained a gap in one sequence, and the other variable positions represent single nucleotide substitutions. 



#6. Export consensus sequence to search in database

#6a. Wrap the consensus in a DNAStringSet
consensus_set <- DNAStringSet(consString)

#6b. Export to FASTA
writeXStringSet(consensus_set, filepath = "consensus_midterm1.fasta")

# Step 6: Identify gene via BLAST
# Consensus sequence BLASTn top hit:
# Gene: Homo sapiens hbb gene (beta globin), partial cds
# Accession number: LC121775.1
# Alignment: 642/642 identical (100%), 0 gaps, E-value = 0



#7. Which individual is most different-translate to amino acid

#7a. Which individual most different

# Use aligned sequences, which are all the same length
aligned_stringset <- as(alignment, "DNAStringSet")

# 7a. Calculate pairwise Hamming distance
dist_matrix <- stringDist(aligned_stringset, method = "hamming")

# 7a. Convert to a full symmetric matrix
dist_matrix_full <- as.matrix(dist_matrix)

# 7a. Calculate the average distance for each sequence
avg_dist <- rowMeans(dist_matrix_full)

# 7a. Identify the sequence with the highest average distance
most_different_index <- which.max(avg_dist)

# 7a. Get the name of the most different individual
names(aligned_stringset)[most_different_index]

#7a. ANSWER: The most different different individual is "Homo_sapiens_6"


#7b. Translate to amino acid

#7b. Extract the most different sequence from the alignment
most_diff_seq <- aligned_stringset[most_different_index]

#7b. Convert DNAString to a plain character string
most_diff_string <- as.character(most_diff_seq)

#7b. Convert the string to characters for seqinr
most_diff_chars <- s2c(most_diff_string)

#7b. Translate to amino acids
most_diff_protein <- translate(most_diff_chars)

#7b. View as a string
paste(most_diff_protein, collapse = "")

#7b. ANSWER: "XSTPRSREGRSQGWA*KSGQSHLLLTFASDTTVFTSNLKQTPWCT*LLWRSLPLLPC
#GAR*TWMKLVVRPWAGWYQGYKTGLRRPIETGHVETEKTLGFLIGTDSLCLLVYFPTLRLLVVYPWTQRF
#FESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRVSLWDP*
#CFLSPSFLWLSSCHRKG"

#7c. export to fasta file

#7c. Wrap the protein sequence in an AAStringSet
protein_set <- AAStringSet(paste(most_diff_protein, collapse = ""))

#7c. Give the sequence a name
names(protein_set) <- names(aligned_stringset)[most_different_index]

#7c. Export to FASTA
writeXStringSet(protein_set, filepath = "HS6_protein.fasta")



#8. Protein match

# Step 8: Protein BLAST
# Best match by BLAST: XP_025213810.1 (hemoglobin subunit beta isoform X1, Theropithecus gelada)
# Biologically relevant match: ACE80932.1 (mutant beta-globin, Homo sapiens)



#9. Disease associated with the mutated hbb gene
#9a. Sickle cell diseaase is associated with the mutated hbb gene which encodes the beta‑globin subunit of hemoglobin
#This is according to a search on the NIH website:
#(https://www.nih.gov/news-events/nih-research-matters/fixing-sickle-cell-disease-gene)

#9a. The most well-known disease-causing mutation in HBB (the sickle cell mutation) is at codon 6 (Glu → Val), which is not covered in my partial fragment.
#Reference for above claim: (https://www.nhlbi.nih.gov/research/sickle-cell-disease)

#9a. Based on the partial sequence data we have, we cannot determine whether Homo_sapiens_6
# carries the specific mutation(s) that cause disease. Therefore, there is not enough evidence
# to conclude that this individual actually has sickle cell disease.






