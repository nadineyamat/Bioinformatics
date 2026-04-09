library("protti")
library("UniprotR")
library("r3dmol")
library("Biostrings")

# Loading accession into R
uniprot_accessions <- readLines("AccessionUniProt")

# Assigning GO info to variable
geneonto <- GetProteinGOInfo(uniprot_accessions)

# Plotting assigned variable
plotgeneonto <- PlotGoInfo(geneonto)

