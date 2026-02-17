# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("Biostrings") # installing biostrings
# install.packages("seqinr")
# install.packages("phangorn")
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("msa")
library("Biostrings") 
library("seqinr")
library("phangorn")
library("msa")
library("stringr")

# seting wd to sub-folder 
setwd("Labsix")

# doing msa on combined file
read.fasta("BWseqs.fasta")
mySequences <- readDNAStringSet("BWseqs.fasta")
names(mySequences) <- c("mitochondrion", "ADRB2", "OPN1SW", "DMP1", "EWF1415")
msaalign <- msaMuscle(mySequences)

# msa analysis: length, gaps, GC content
print(msaalign, show="complete")
convert <- as(msaalign, "DNAStringSet")
letterFrequency(convert, "-")
letterFrequency(convert, letters="GC", as.prob = TRUE)

# converting to seqinr format & computing distance matrix
msaseqinr <- msaConvert(msaalign, type="seqinr::alignment")
as.matrix(dist.alignment(msaseqinr, "identity"))

# translating 1 file into amino acid seq
read.fasta("BWseq2.fasta")
dna_seq <- readDNAStringSet("BWseq2.fasta")
aa_seq <- Biostrings::translate(dna_seq)
aa_seq

# writing alignment to a file
Alignment_phyDat <- msaConvert(msaalign, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")
msa