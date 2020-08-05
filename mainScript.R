# -----------------------------------------------
# Read libraries

library(raster)

# -----------------------------------------------
# Read functions

source("https://raw.githubusercontent.com/jorgeassis/spatialAutocorrelaiton/master/functions.R")

# -----------------------------------------------
# Read occurrence records from Assis et al. (2020)
# Assis, J., Fragkopoulou, E., Frade, D. et al. A fine-tuned global distribution dataset of marine forests. Sci Data 7, 119 (2020). https://doi.org/10.1038/s41597-020-0459-x

dataset <- extractDataset("brownAlgae",pruned=TRUE)
occurrenceRecords <- subsetDataset(dataset,taxa="Laminaria digitata",status="accepted")

# -----------------------------------------------
# Read raster from Assis et al. (2018)
# Assis, J, Tyberghein, L, Bosch, S, Verbruggen, H, Serrão, EA, De Clerck, O. Bio‐ORACLE v2.0: Extending marine data layers for bioclimatic modelling. Global Ecol Biogeogr. 2018; 27: 277– 284. https://doi.org/10.1111/geb.12693







distanceUncorr <- data.frame(Predictor=names(rasterLayers),Distance=NA)
for( i in 1:length(names(rasterLayers))) {
  distanceUncorr[i,2] <- spatialAutocorrelation(occurrenceRecords=occurrenceRecords,subset(rasterLayers,i),autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif)
}

distanceUncorrPlot <- ggplot(distanceUncorr[sort(distanceUncorr[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
  geom_bar( aes(x= reorder(Predictor, Distance) , y=Distance), stat="identity", fill="black", alpha=0.5) +
  coord_flip() + theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
  ) + labs(x = "Predictor") + 
  labs(y = "Spatial correlation (km)") + geom_hline(aes(yintercept=round(mean(distanceUncorr[,2]),digits=2)   ),color="Black", linetype="dashed", size=0.3) +
  annotate("text", y = round(mean(distanceUncorr[,2]),digits=2) + 2 , x = 1 , label = paste0( round(mean(distanceUncorr[,2]),digits=2)," Km") , hjust = 0)

pdf(file = paste0(resultsFolder.i,"/SpatialAutocorrelation.pdf"), width=12, height=8 )
print(distanceUncorrPlot)
dev.off()

meanCorrDistance <- mean(distanceUncorr[,2])
maxCorrDistance <- max(distanceUncorr[,2])
occurrenceRecords <- spatialThinning(occurrenceRecords,rasterLayers, ifelse(meanCorrDistance > 25 , 25 , meanCorrDistance) )

resultsSummary[1,"meanCorrDistance"] <- meanCorrDistance
resultsSummary[1,"maxCorrDistance"] <- maxCorrDistance
resultsSummary[1,"recordsFinal"] <- nrow(occurrenceRecords)


