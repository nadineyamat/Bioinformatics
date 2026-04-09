# install.packages ("ape")
library(ape)

# set working directory
setwd("midterm2")

# assigning tree to variable
my_tree <- read.tree("metazoa_alignment.5k.fasta.raxml.bestTree")

# Visualizing unrooted tree
plot(my_tree)

# rooting tree with outgroup and visualizing
rooted_tree <- root(my_tree, outgroup = "Grantia_compressa")
plot(rooted_tree)

# editing tree w/ less cramping
plot(rooted_tree, cex = 0.5)
