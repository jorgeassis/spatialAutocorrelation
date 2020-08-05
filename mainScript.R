# -----------------------------------------------
# Read libraries

library(raster)

# -----------------------------------------------
# Read occurrence records from Assis et al. (2020)
# Assis, J., Fragkopoulou, E., Frade, D. et al. A fine-tuned global distribution dataset of marine forests. Sci Data 7, 119 (2020). https://doi.org/10.1038/s41597-020-0459-x

source("https://raw.githubusercontent.com/jorgeassis/marineforestsDB/master/sourceMe.R")
dataset <- extractDataset("brownAlgae",pruned=TRUE)
occurrenceRecords <- subsetDataset(dataset,taxa="Laminaria digitata",status="accepted")

# -----------------------------------------------
# Read raster from Assis et al. (2018)
# Assis, J, Tyberghein, L, Bosch, S, Verbruggen, H, Serrão, EA, De Clerck, O. Bio‐ORACLE v2.0: Extending marine data layers for bioclimatic modelling. Global Ecol Biogeogr. 2018; 27: 277– 284. https://doi.org/10.1111/geb.12693



