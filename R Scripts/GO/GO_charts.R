###########################################################
# GO_charts.R
#
# Utilizes the RRVGO package to create treemap charts of
# our RNA-Seq Gene Ontology data (output from GOrilla).
#
# Code by Anders Ohman
# Spring 2022
###########################################################

########################
# User Settings
# rrvgo::shiny_rrvgo()

# Which GOrilla result set to use? host or colony
whichSet = "colony"

pvalCutoff = 0.00001 #10E-6=0.000001, 10E-5=0.00001, 10E-3=0.001

GOthreshold = 0.7

# Which plot to display? all, scatter, heatmap, treemap, wordcloud
whichPlot = "treemap"

colony <- "Colony GO Output - pval 10^-3.txt"
host <- "Host GO Output - pval 10^-3.txt"
organism <- "org.Rn.eg.db" # or "org.Hs.eg.db"

########################
# Load Libraries
library(tidyverse)
library(rrvgo)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/GO", "")
gowd <- paste(githubwd,"Data Files/go/",sep="")

########################
# Read input files

setwd(gowd)
colonyGO <- read.csv2(colony, sep="\t")
hostGO <- read.csv2(host, sep="\t")

if(whichSet=="colony"){
  go_analysis <- colonyGO
}
if(whichSet=="host"){
  go_analysis <- hostGO
}

go_analysis <- transform(go_analysis, pValue = as.numeric(pValue))

go_analysis <- subset(go_analysis, pValue<pvalCutoff) 

#goTerm <- as.vector(colonyGO['GO.term'])
#pVal <- as.numeric(colonyGO['pValue'])

simMatrix <- calculateSimMatrix(go_analysis$GO.term,
                                orgdb=toString(organism),
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_analysis$pValue), go_analysis$GO.term)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=GOthreshold,
                                orgdb=toString(organism))

if(whichPlot=="all"||whichPlot=="heatmap"){
  heatmapPlot(simMatrix,
              reducedTerms,
              annotateParent=TRUE,
              annotationLabel="parentTerm",
              fontsize=6)
}
if(whichPlot=="all"||whichPlot=="scatter"){
  scatterPlot(simMatrix, reducedTerms)
}
if(whichPlot=="all"||whichPlot=="treemap"){
  treemapPlot(reducedTerms,
              inflate.labels = TRUE,      # set this to TRUE for space-filling group labels - good for posters
              #fontsize.labels = 100, # if inflate is not on
              aspRatio = 1.0, #1.0 for host, 1.9 for colony
              overlap.labels = 1, # 0 shows no overlap, 1 shows many
              lowerbound.cex.labels = 1,   # try to draw as many labels as possible (still, some small squares may not get a label)
              )
}
if(whichPlot=="all"||whichPlot=="wordcloud"){
  wordcloudPlot(reducedTerms, min.freq=1, colors="black")
}