#install.packages("BiocManager")
#BiocManager::install("biomaRt")

library(Biostrings)
library(biomaRt)
library(ggplot2)

#setwd
setwd("/Users/devinquigley/Desktop/Bioinformatics/Bioinformatics-/Midterm2")

#8 AI/Previous code
#Load in file and read sequences
fasta_file <- "metazoa_alignment.gene.fasta"
seqs <- readDNAStringSet(fasta_file)

#Look at the names of sequences
names(seqs)

#Extract the Homo sapiens from sequence list
hs_seq <- seqs[grep("Homo_sapiens", names(seqs))]

#Remove gaps and ambiguus bases for translation due to error
hs_seq_clean <- DNAString(gsub("[-N]", "", as.character(hs_seq)))

#Translation of cleaned sequence for homo_sapiens
hs_protein <- translate(hs_seq_clean)

#Warning message:In .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet],
#:last 2 bases were ignored

#Wrap in AAStringSet due to error 
hs_protein_set <- AAStringSet(hs_protein)

#Naming the protein for exportation
names(hs_protein_set) <- "Homo_sapiens_protein"

#Write the protein to FASTA file
writeXStringSet(hs_protein_set, "Homo_sapiens_protein.fasta")

#10 Stack Overflow/AI
# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast")

#Retrive GO terms using Accession number
go_terms <- getBM(
  mart = ensembl,
  attributes = c("uniprotswissprot", "go_id", "namespace_1003", "name_1006"),
  filters = "uniprotswissprot",
  values = "P54098"
)

#View results 
head(go_terms)

# Get one term from each ontology
mf <- go_terms[go_terms$namespace_1003 == "molecular_function", ][1, ]
bp <- go_terms[go_terms$namespace_1003 == "biological_process", ][1, ]
cc <- go_terms[go_terms$namespace_1003 == "cellular_component", ][1, ]

#Stacks three GO terms into one table
selected_go <- rbind(mf, bp, cc)

#View results
selected_go

#Graph selected GO results: Barplot
barplot(rep(1, 3),
        names.arg = selected_go$namespace_1003,
        main = "GO terms for P54098",
        ylab = "Presence")

text(x = c(0.7, 1.9, 3.1), y = 0.5,
     labels = selected_go$name_1006,
     srt = 45, adj = 1, xpd = TRUE)






#test

library(ggplot2)

go_data <- data.frame(
  ONTOLOGY = c("MF", "BP", "CC"),
  Term = c("DNA binding",
           "DNA replication",
           "gamma DNA polymerase complex"),
  Count = c(1, 1, 1)
)

ggplot(go_data, aes(x = Term, y = Count, fill = ONTOLOGY)) +
  geom_col() +
  geom_text(aes(label = Count), hjust = -0.2) +
  coord_flip() +
  labs(title = "GO Terms for P54098", x = NULL, y = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom")











