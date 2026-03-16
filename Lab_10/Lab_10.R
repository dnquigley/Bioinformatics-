#setting the working directory 
setwd("/Users/devinquigley/Desktop/Bioinformatics/Bioinformatics-")

#install packages
#install.packages("UniprotR")
#install.packages("protti")
#BiocManager::install("GenomicAlignments")
#install.packages("r3dmol")
#install.packages("easyr")
#install.packages("shiny")


#Library packages
library(UniprotR)
library(protti)
library(GenomicAlignments)
library(r3dmol)
library(Biostrings)
library(seqinr)
library(easyr)
library(dplyr)
library(shiny)

#Read in DNA Fasta file
KU18 <- readDNAStringSet("Lab_10/VAME_KU18.fasta")
KU18

#Check class of read
class(KU18)

#Convert to amino acid sequence 
KU18_protein <- Biostrings::translate(KU18)
KU18_protein

#Export amino acid sequence to Fasta file-AI 
write.fasta(sequences = KU18_protein,
            names = "KU18",
            file.out = "KU18_protein.fasta")

#Read in text file with accession numbers-AI 
accessions <- read.csv("Lab_10_Accession.txt", header = FALSE)
accessions

#Check class 
class(accessions)

#Convert to character string-AI 
acc_string <-accessions[,"V1"]

class(acc_string)

#Alternate accession numbers
#Used provided accessions from number 2
accessions <- c("P0A799", "P08839")

#Get gene ontology terms
acc_go <- GetProteinGOInfo(accessions)
View(acc_go)

#Plot the gene ontology results for result 
PlotGoInfo(acc_go)

#Handy visualization for publications-Part AI part Internet
PlotGOAll(GOObj = acc_go, Top = 10, directorypath = getwd(), width = 10, height = 8)
getwd()


#Info on disease and pathologies associated with gene
#Pathology associates
pathology <- GetPathology_Biotech(accessions)
View(pathology)

# Disease associations-AI
Get.diseases(pathology)

#Access structural info
stuctural_info <-fetch_uniprot(accessions) 

#Pull available structural info from Protein database
#Used provided accessions from number 11
acc <- c("1ZMR", "2HWG")

structural_protein <-fetch_pdb(acc)

#Any available 3D structures 
fetch_alphafold_prediction("P0A799") 

fetch_alphafold_prediction("P08839") 







