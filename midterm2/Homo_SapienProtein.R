library("Biostrings")

# Read fasta into R
metazoa_alignment <- readDNAStringSet("metazoa_alignment.gene.fasta")

# Pluck out Homo_Sapien gene
names(metazoa_alignment)
homo_sapien_gene <- ("Homo_sapiens")

# Extract sequence as a DNAStringSet object
target_dna <- metazoa_alignment[["Homo_sapiens"]]

# Replace non-base characters (gaps, periods) with 'N' to prevent translation errors
target_dna_n <- DNAString(gsub("[^ATGC]", "N", as.character(target_dna)))

# Translate to protein sequence
protein_seq <- translate(target_dna_n, if.fuzzy.codon="solve")

# Visualize protein seq
as.character(protein_seq)

# Write to fasta file
writeXStringSet(AAStringSet(protein_seq), filepath = "human_protein.fasta")
