###########################################################
# PCA_hPTM.R
#
# Generate Principal Component Analysis Plots from our
# hPTM-MS Data.
#
# Code by Anders Ohman and John Santiago
# Sanders Lab, Brown University, Pathobiology
# Spring 2020 - 2022
###########################################################

########################
# User Settings

# Graph All or Metadata samples?
samplesToGraph = "Metadata"

metadataColor <- "adultfetal"
metadataColorFactors <- c("Adult","Fetal")
colors <- c("green4","red")

metadataShape <- "Group"
metadataShapeFactors <- c("A-20","F-20","FEHS","FEHD","ADH","FEHL")
shapes <- c(17,17,19,15,19,18) # Identify shape; NOTE: see http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r

metadataLine <- "replicate"

inputMetadata <- NULL

# 2D or 3D plot?
typeOfPCA = "3D"

  # IF 2D, what style of plot? Shape or Label points?
  styleOfPCA = "Shape"

  # IF 2D, What Principal Components? (Default is X=1, Y=2)
  PCX = 1
  PCY = 2

  # IF 2D, draw lines between points by replicate?
  lines = FALSE

  # IF 3D, what angle to show? (Default is leftToRight=60, bottomToTop=20)
  ## Also try 40/20, 60/65
  leftToRight = 60
  bottomToTop = 40

  # What size should the points be? (Default is 2.5 for 3D)
  ## NOTE: Only working in 3D at the moment, I will try to fix for ggplot2.
  pointSize = 2.5

#Include LCM liver?
includeLiver <- FALSE
  
########################
# Load Libraries
library(tidyverse)
library(edgeR)
library(plot3D)
library(ggfortify)
library(directlabels)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/PCA", "")
protwd <- paste(githubwd,"Data Files/proteomics/",sep="")
setwd(paste(protwd,"hPTM/",sep=""))

########################
# Read input files
cellData <- read.csv("Cells - 020716.csv")
liverData <- read.csv("Liver - 080919.csv")
colonyData <- read.csv("Colonies - 012320.csv")

merged <- merge(cellData, colonyData, by="rnm", all=TRUE)

if(includeLiver){
merged <- merge(merged, liverData, by="rnm", all=TRUE)
}

row.names(merged) <- merged[,1]
merged <- merged[,-1]
merged <- merged[complete.cases(merged),]

groups <- read.csv("hptm_metadata.csv",row.names=1)

if(includeLiver==FALSE){
groups <- groups[-c(1:6),]
}

pcadata <- merged

col.order <- row.names(groups)
pcadata <- pcadata[,col.order]

########################
# Set up individual groups
if (samplesToGraph!="All"){
# Relabel for graph display
groups$adultfetal <- factor(groups$adultfetal, labels=metadataColorFactors)
#groups$viability <- factor(groups$viability, levels=c("V","N"), labels=c("Viable","Nonviable"))
#groups$timepoint <- factor(groups$timepoint, levels=c("pre","three","six"), labels=c("0h","3h","6h"))
groups$Type <- factor(groups$Type, levels=c("Cells","HostCol","Liver"), labels=c("Cells","HostCol","Liver"))
groups$Group <- factor(groups$Group, levels=c("A-20","F-20","FEHS","FEHD","ADH","FEHL"), labels=c("A-20","F-20","FEHS","FEHD","ADH","FEHL"))
groups$replicate <- factor(groups$replicate)

# Run the actual PCA calculation
pca <- prcomp(t(pcadata), scale.=TRUE)

# Create merged metadata/data file for use with 3D Plots
if(typeOfPCA=="3D"){
  merged3D <- data.frame(pca$x, groups)
  
  shapeIndex <- metadataShapeFactors
  shapeValues <- shapes
  merged3D$shape <- shapeValues[match(merged3D[[metadataShape]], shapeIndex)]
  
  colorIndex <- metadataColorFactors
  colorValues <- colors
  merged3D$color <- colorValues[match(merged3D[[metadataColor]], colorIndex)]
  
  #merged3D$unique <- c(1,2,3,1,2,3,1,2,3,4,5,6,4,5,6,4,5,6) #See line 179 - trying to identify group
}
}

########################
# ALL
if (samplesToGraph=="All"){
  pca <- prcomp(t(pcadata), scale.=TRUE)
  if (typeOfPCA=="3D") {
    plot <- scatter3D(merged3D[,1],merged3D[,2],merged3D[,3],
              colvar = NULL, col = merged3D$color, pch = merged3D$shape, theta = leftToRight, phi = bottomToTop,
              xlab = "PC1", ylab = "PC2", zlab = "PC3", cex = pointSize, col.grid = "black",
              cex.axis = 1.0, cex.lab = 1.5, main=paste(samplesToGraph," Samples (3D)",sep=""), cex.main=2.0,
              bty = "b2", ticktype = "detailed", d = 2)
  }
  if (typeOfPCA=="2D"){
    if (styleOfPCA=="Shape"){
      #Create plot with this criteria:
      #Colored based on replicate group (i.e. LV1)
      #Group labeled at center point of the 3 replicates
      #Shaped based on fatty (triagle is lean, round is fatty)
      plot <- autoplot(pca, x=PCX, y=PCY, data=groups, colour = toString(metadataColor), shape = toString(metadataShape), frame = TRUE, frame.type = 'norm') +
        labs(shape=toString(metadataShape)) +
        ggtitle(paste("PCA - ", samplesToGraph, " Samples", sep=""))
      # Use Directlabels to put the group names on the plot
      plot <- direct.label.ggplot(plot, method=smart.grid)
    }
    # Label each point for individual sample tracking
    if (styleOfPCA=="Label"){
      if (!is.null(inputMetadata)) {
        plot <- autoplot(pca, x=PCX, y=PCY, data = groups, colour = groups[[metadataColor]], shape = FALSE, label.size = 7)+
          labs(color = "metadataColor")+
          ggtitle(paste("PCA - ", samplesToGraph, " Samples (Labeled)", sep=""))
      }
      if (is.null(inputMetadata)) {
        plot <- autoplot(pca, x=PCX, y=PCY, shape = FALSE, label.size = 7)+
          ggtitle(paste("PCA - ", samplesToGraph, " Samples (Labeled)", sep=""))
      }
    }
  }
}

########################
# Metadata PCAs
if (samplesToGraph=="Metadata"){
  if (typeOfPCA=="3D") {
    plot <- scatter3D(merged3D[,1],merged3D[,2],merged3D[,3],

              # NOTE: This is a work in progress trying to add lines between points
              # in the 3D plot. Currently it is not working, and plot3D may
              # not be flexible enough to do so automatically by a group variable.
              # John has done this successfully with manual plotting in a different
              # script, but this requires a lot of adjustment for every plot.
              # Documentation: http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization#scatter-plots
              #
              # Type B means Both lines and points, but this just draws lines between everything...
              #
              #type="b", lwd = 4,
              #
              #scatter3D(x, y, z, phi = 0, bty = "g", type = "b",
              #          ticktype = "detailed", pch = 20,
              #          cex = c(0.5, 1, 1.5))
              # [order(merged3D$unique)] #See line 126 - trying to identify group

              colvar = NULL, col = merged3D$color, pch = merged3D$shape, theta = leftToRight, phi = bottomToTop,
              xlab = "PC1", ylab = "PC2", zlab = "PC3", cex = pointSize, col.grid = "black",
              cex.axis = 1.0, cex.lab = 1.5, main=paste(samplesToGraph," Samples (3D)",sep=""), cex.main=2.0,
              bty = "b2", ticktype = "detailed", d = 2)
  }
  if (typeOfPCA=="2D"){
    #Shapes (circle, triangle, square) for each point
    if (styleOfPCA=="Shape"){
      plot <- autoplot(pca, x=PCX, y=PCY, data=groups)+ # What data to use for the chart
        geom_point(aes(colour=groups[[metadataColor]], shape=groups[[metadataShape]], size=pointSize))+ # What style of graph to make
        #stat_ellipse(mapping=aes(color=groups[[metadataColor]]),type="t")+
        geom_dl(aes(label=NA,colour=groups[[metadataColor]]),method="smart.grid")+ # How to color the chart
        labs(color = metadataColor)+ # Label for viability legend
        labs(shape = metadataShape)+ # Label for timepoint legend
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples", sep=""))+ # Title of chart
        scale_color_manual(values=colors)+ # Manually set the colors of the nodes, sequentially
        theme_bw()+ # Change background theme
        scale_size(guide = 'none') # Hide the size from the legend
      # Draw lines between points in path order
      if (lines){
        plot <- plot + geom_path(aes(color=toString(metadataColor),
                                     linetype=toString(metadataLine),
                                     size=1.5)) + # Draw lines between replicate over timecourse
          
          #FIXME#######################################################################################################################
          scale_linetype_manual(values=c("solid", "longdash", "dotted","solid", "longdash"))+
          ########################################################################################################################
          
          labs(linetype = metadataLine) # Label for replicate line legend
      }
    }
    
    # Label each point for individual sample tracking
    if (styleOfPCA=="Label"){
      plot <- autoplot(pca, x=PCX, y=PCY, data = groups, colour = toString(metadataColor), shape = FALSE, label.size = 7)+
        labs(color = metadataColor)+
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples (Labeled)", sep=""))
    }
  }
}

#Print out the selected plot
plot


########################
# Calculate Euclidean Distances based on group

df.pca.x <- as.data.frame(pca$x)
df.pca.x$groups <- paste(groups$adultfetal,groups$Type,sep=".")
ptm.pca.centroids <- aggregate(df.pca.x[,1:13], list(Type = df.pca.x$groups), mean)

ugroups <- unique(df.pca.x$groups)
distances <- matrix(nrow=4, ncol=4)
colnames(distances) <- ugroups
rownames(distances) <- ugroups

#distances["Adult.Cells","Adult.HostCol"] <- dist(rbind(ptm.pca.centroids[ptm.pca.centroids$Type == "Adult.Cells",2:4],ptm.pca.centroids[ptm.pca.centroids$Type == "Adult.HostCol",2:4]), method = "euclidean")

for(i in 1:4){
  varOne <- colnames(distances)[i]
  for(j in 1:4){
    varTwo <- colnames(distances)[j]
    distances[varOne,varTwo] <- dist(rbind(ptm.pca.centroids[ptm.pca.centroids$Type == varOne,2:4],ptm.pca.centroids[ptm.pca.centroids$Type == varTwo,2:4]), method = "euclidean")
  }
}

distances


########################
# Calculate Euclidean Distances based on specific subgroup

df.pca.x <- as.data.frame(pca$x)
df.pca.x$groups <- groups$Group
ptm.pca.centroids <- aggregate(df.pca.x[,1:13], list(Type = df.pca.x$groups), mean)

ugroups <- unique(df.pca.x$groups)
distances <- matrix(nrow=6, ncol=6)
colnames(distances) <- ugroups
rownames(distances) <- ugroups

#distances["Adult.Cells","Adult.HostCol"] <- dist(rbind(ptm.pca.centroids[ptm.pca.centroids$Type == "Adult.Cells",2:4],ptm.pca.centroids[ptm.pca.centroids$Type == "Adult.HostCol",2:4]), method = "euclidean")

for(i in 1:6){
  varOne <- colnames(distances)[i]
  for(j in 1:6){
    varTwo <- colnames(distances)[j]
    distances[varOne,varTwo] <- dist(rbind(ptm.pca.centroids[ptm.pca.centroids$Type == varOne,2:4],ptm.pca.centroids[ptm.pca.centroids$Type == varTwo,2:4]), method = "euclidean")
  }
}

distances
