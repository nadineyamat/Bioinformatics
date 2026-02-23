library("msa")
library("seqinr")
library("BiocManager")
library("Biostrings")

# Seting wd to sub-folder 
setwd("Midterm")

# Loading file in R and performing MSA
read.fasta("sequences.fasta")
mySeq <- readDNAStringSet("sequences.fasta")
myAlign <- msaMuscle(mySeq)

# Determining alignment is good or bad
   # Looking at MSA visually
print(myAlign, show="complete")
   # Looking for a highly conserved consensus sequence
myAlignMatrix <- msaConvert(myAlign, type="seqinr::alignment")
genmatrix <-as.matrix(myAlignMatrix)
consmatrix <- consensusMatrix(myAlign, baseOnly=TRUE)
cons_count <- apply(consmatrix, 2, max)
cons_proportion <- cons_count / nrow(genmatrix)
print(cons_proportion)

# Determined MSA is good alignment; looking at the msa visually shows nearly 
  # identical letters aligned vertically across sequences. Quantitative analysis 
  # of the consensus sequence supports this as well, with most positions showing 
  # full conservation, 1, or differences of about 0.10-0.15. Overall, I believe the 
  # msa is good.

# Analzying multiple sequence alignment
myAlignFreq <- as(msaalign, "DNAStringSet")
letterFrequency(myAlignFreq, "-")
msaConsensusSequence(myAlign)
letterFrequency(myAlignFreq, letters="GC", as.prob = TRUE)


# Calculating distant between sequences via matrix
as.matrix(dist.alignment(myAlignMatrix, "identity"))

# Notable samples: Samples 4 and 10 have moderate differences, but sample 6 
  # appears the most different, due to the large number in the matrix.


# The rest of the samples are identical to each other, evident of the "0" 


# Gene of interest appears to be the hemoglobin subunit beta gene.
# Best matchs & correlating accession numbers: 
# homo 6 is AY356351.1, homo 4 is LC121775.1, homo 10 is LC121775.1
# Achieved by pasting each sample from the fasta file as a query in blastn.


# Best protein match accession number: ACE80932.1
# Achieved via pasting sample as a query into blastx.


# Disease of interest appears to be beta-thalessemia.
# Achieved by searching through GenBank (nulceotide & protein), and Blastn 
  # and Blastx.