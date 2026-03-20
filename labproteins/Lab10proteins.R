# Set working directory
setwd("labproteins")
# install.packages("UniprotR")
# install.packages("protti")
# install.packages("r3dmol")
library("protti")
library("UniprotR")
library("r3dmol")
library("Biostrings")

# Reading text file into R & verifying
uniprot_vector <- readLines("Proteins_Accession")
head(uniprot_vector)

# Obtaining Gene Ontology terms
geneonto <- GetProteinGOInfo(uniprot_vector)

# Plotting results
plotgeneonto <- PlotGoInfo(geneonto)

# Commmand failed; Checking if Gene Ontology results include all 3 key terms
print(geneonto)

# N/A for molecular info; plotting Subcellular & Biological info
Plot.GOSubCellular(geneonto)
PlotGOBiological(geneonto)

# Saving plots to repo
PlotGOBiological(GOObj = geneonto, Top = 10, directorypath = getwd())
Plot.GOSubCellular(GOObj = geneonto, Top = 10, directorypath = getwd())

# Checking for associated diseases 
knowndiseases <- GetPathology_Biotech(uniprot_vector)
print(knowndiseases)
diseases <- Get.diseases(knowndiseases)
print(diseases)

# Accessing structural information 
untiprotdata <- fetch_uniprot(uniprot_vector)
print(untiprotdata, width = Inf)

# Defining PDB IDs
proteinid <- c("1ZMR", "2HWG")

# Fetching information for the IDs
pdb_data <- fetch_pdb(pdb_ids = proteinid)
print(pdb_data)

# Checking for available 3d structures
fetch_alphafold_prediction(uniprot_vector)
