############################################################################
# Script Purpose: Analysis of Flow Cytometry Annexin-V Apoptosis Data      #
# Date: ver 3, 27-04-20                                                    #
# Author: Peter D. Jones                                                   #
############################################################################

#load required libraries
library(flowCore)
library(flowViz)
library(flowStats)
#uses lattice packages

#load data files specified in annotation.txt - only loads "Lin" data, ignores Log and Area data
flowData <- read.flowSet(phenoData="annotation.txt", transformation=FALSE, column.pattern="Lin", alter.names=TRUE)

#replace sample names with those specified in annotation.txt
sampleNames(flowData) <- as.character(pData(flowData)[, "SampleID"])

#apply asinh transformation to fluorescent channels(log, logicle would also work)
tData <- transform(flowData, transformList(colnames(flowData)[3:8], asinh))

#create plot of FS vs SS
unGatedPlot <- xyplot(`FS.Lin`~`SS.Lin`, data = tData)

#Gate on FS and SS to remove debris, clumpy cells etc (crude rectangular gate)
rg <- rectangleGate("FS.Lin"=c(5000, 35000), "SS.Lin"=c(5000, 50000), filterId="rectangle")
gData <-Subset(tData, rg)

#create plot of gated FS and SS data
gatedPlot <- xyplot(`FS.Lin`~`SS.Lin`, data = gData)

#normalise data to correct for variability in staining
normGData <- warpSet(gData, "FL.1.Lin")

#plot ovelapping density plot of Annexin-V-FITC on un-normalised data
overlapPlot <- densityplot(~ `FL.1.Lin`, gData, main = "density plot of FITC")

#normalised overlap plot
normOverlapPlot <- densityplot(~ `FL.1.Lin`, normGData, main = "density plot of FITC")

#create gate on FITC fluorescence to show %apoptotic, plot graph showing %apo for each sample
apoGated <- rectangleGate("FL.1.Lin"=c(7.5, Inf), filterId="Apoptotic")
#create a gating set, apply the gate, create gated plot
gs<-GatingSet(normGData)
add(gs, apoGated, parent="root")
gs
getNodes(gs)
recompute(gs)
percentPlot <- plotGate(gs, "Apoptotic", type="densityplot")

#extract numerical results data 
result = filter(normGData, apoGated)

#set up a dataframe to store the results
finalResults <- data.frame(SampleID=character(), PositiveCells=integer(), TotalCells=integer(), Percentage=double())

#extract sample names
sampleNamesOnly <- as.character(pData(flowData)[, "SampleID"])

#loop through each sample, extracting data and storing it in the finalResults data frame, then save to csv
for (i in 1:length(sampleNamesOnly)) {
resultsTemp = data.frame(SampleID=sampleNamesOnly[i], PositiveCells=summary(result[[i]])$true, TotalCells=summary(result[[i]])$n, Percentage=summary(result[[i]])$p)
finalResults <- rbind(finalResults, resultsTemp)
}

write.csv(finalResults, "apoptosisResults.csv")

#create pdf, print graphs to it
pdf("graphs2.pdf")
unGatedPlot
gatedPlot
overlapPlot
normOverlapPlot
percentPlot
dev.off()