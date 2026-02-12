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
library(msa)
library("stringr")

# created subfolder 
setwd("Labsix")

seq1 <- read.fasta("BWseq1.fasta")

# msa
read.fasta("BWseqs.fasta")
mySequences <- readDNAStringSet("BWseqs.fasta")
names(mySequences) <- c("mitochondrion", "ADRB2", "OPN1SW", "DMP1", "EWF1415")

# msa analysis: length, gaps, GC content
msaalign <- msaMuscle(mySequences)
print(msaalign, show="complete")
convert <- as(msaalign, "DNAStringSet")
letterFrequency(convert, "-")
letterFrequency(convert, "GC")
letterFrequency(convert, letters="GC", as.prob = TRUE)

# converting to seqinr format & compute distance matrix
msaseqinr <- msaConvert(msaalign, type="seqinr::alignment")
as.matrix(dist.alignment(msaseqinr, "identity"))

# translating 1 file into amino acid sequence
read.fasta("BWseq2.fasta")
dna_seq <- readDNAStringSet("BWseq2.fasta")
dna_seq <- dna_seq[width(dna_seq) > 0]
aa_seq <- Biostrings::translate(dna_seq)
aa_seq

# write alignment to a file
Alignment_phyDat <- msaConvert(msaalign, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")

