###########################################################
# PCA_RNA.R
#
# Generate Principal Component Analysis Plots from our
# RNA-Seq Data.
#
# Code by Anders Ohman and John Santiago
# Sanders Lab, Brown University, Pathobiology
# Spring 2020 - 2022
###########################################################

########################
# User Settings

#Graph Cells, LCM, or Both?
inputToGraph = "Both"

# Batch correct the data using SVA and ComBat-seq?
if (inputToGraph == "Both"){ comBatSeq = TRUE } else { comBatSeq = FALSE }

# Graph All or Metadata samples?
samplesToGraph = "Metadata"

# What metadata to use to define COLOR?
metadataColor <- "adultfetal"
metadataColorFactors <- c("Adult","Fetal")
colors <- c("blue","red")

# What metadata to use to define SHAPE?
if (inputToGraph == "Both" || inputToGraph == "Cells"){ metadataShape <- "group" }
if (inputToGraph == "LCM"){ metadataShape <- "adultfetal" }

# How to define SHAPE?
if (inputToGraph == "Both"){
  metadataShapeFactors <- c("Adult","Dual","LAP","Single","Host","Colony")
  shapes <- c(16,15,18,19,17,17) # Identify shape; NOTE: see http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  order <- c("A","B","C","D","E","F")
}

if (inputToGraph == "Cells"){
  metadataShapeFactors <- c("Adult","Dual","LAP","Single")
  shapes <- c(16,15,18,19) #Identify shape; NOTE: see http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  order <- c("A","B","C","D")
}

if (inputToGraph == "LCM"){
  metadataShapeFactors <- c("Adult","Fetal")
  shapes <- c(17,17) #Identify shape; NOTE: see http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  order <- c("A","B")
}


metadataLine <- "replicate"

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
  if (inputToGraph == "Both" || inputToGraph == "Cells"){
    leftToRight = 10
    bottomToTop = 65
  }
  if (inputToGraph == "LCM"){
    leftToRight = -30
    bottomToTop = 40
  }

  # What size should the points be? (Default is 2.5 for 3D)
  ## NOTE: Only working in 3D at the moment, I will try to fix for ggplot2.
  pointSize = 2.5

########################
# Load Libraries
library(tidyverse)
library(edgeR)
library(plot3D)
library(ggfortify)
library(directlabels)
if(comBatSeq){ library(sva) }
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(annotate)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/PCA", "")
rnawd <- paste(githubwd,"Data Files/rna-seq/",sep="")
setwd(paste(rnawd,"AF_Genewiz",sep=""))

########################
# Read input files
countdata_cells <- read.csv("AF_Cells_raw_counts.csv",row.names=1)
countdata_cells <- dplyr::select(countdata_cells, -Gene.name) #Remove gene name column from end
countdata_cells <- countdata_cells[rowSums(countdata_cells[])>0,] #Remove all rows that are only 0

countdata_lcm <- read.csv("AF_LCM_all_counts_no246.csv",row.names=1)

# # Merge cell and LCM tables
# if (inputToGraph=="Both"){
#   countdata_merged <- merge(countdata_cells, countdata_lcm, by="row.names", all=TRUE)
#   countdata <- countdata_merged
#   countdata[is.na(countdata)] <- 0
#   write.csv(countdata, "AF_Merged_raw_counts.csv")
# }

countdata_merged <- read.csv("AF_Merged_raw_counts.csv",row.names="Row.names")
#countdata_merged <- select(countdata_merged, -X)

if (inputToGraph=="Both"){
  countdata <- countdata_merged
  groups <- read.csv("AF_Merged_metadata.csv",row.names=1)
}

if (inputToGraph=="Cells"){
  countdata <- countdata_cells
  groups <- read.csv("AF_Cells_metadata.csv",row.names=1)
}

if (inputToGraph=="LCM"){
  countdata <- countdata_lcm
  groups <- read.csv("AF_LCM_metadata_no246.csv",row.names=1)
}

if (comBatSeq){
  # SVA Batch Correct
  count_matrix <- as.matrix(countdata)
  batch <- groups$batch
  adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)
  countdata <- adjusted
}

########################
# Normalize and remove low-count reads

#Saves normalized CPM file and a count table both with low count reads removed

countdata <- countdata[,row.names(groups)]
x <- DGEList(countdata)
#class(x)
geneid <- rownames(x)
x$samples <- merge(x$samples,groups,by="row.names")
row.names(x$samples) <- x$samples[,1]
x$samples <- x$samples[,-1]

x <- DGEList(countdata)
geneid <- rownames(x)

# Trim counts
keep.exprs <- filterByExpr(countdata)
filt1 <- x[keep.exprs,,keep.lib.sizes=FALSE]
# dim(filt1)
y <- filt1

#Normalize via TMM method
z <- calcNormFactors(y, method = "TMM")

##using z gives the normalized cpm, using y gives the not normalized cpm
cpmsc <- cpm(z)
#cpmsc <- cpm(y)
pcadata <- cpmsc

# col.order <- row.names(groups)
# pcadata <- pcadata[,col.order]

# Write annontated, normalized data
# speciesDB <- "org.Rn.eg.db"
# conFrom <- "ENSEMBL"
# conTo <- "SYMBOL"
# pcadata <- as.data.frame(pcadata)
# pcadata$ENSEMBL <- row.names(pcadata)
# pcadata$SYMBOL <- mapIds(get(speciesDB),
#                     keys=row.names(pcadata),
#                     column=conTo,
#                     keytype=conFrom,
#                     multiVals="first")
# write.csv(pcadata, "AF_All_filterByExpr-normalized_comBatSeq-corrected.csv")
# #DO NOT USE heatmaply(pcadata)

########################
# Set up individual groups
if (samplesToGraph!="All"){
# Relabel for graph display
groups$adultfetal <- factor(groups$adultfetal, labels=metadataColorFactors)
#groups$viability <- factor(groups$viability, levels=c("V","N"), labels=c("Viable","Nonviable"))
#groups$timepoint <- factor(groups$timepoint, levels=c("pre","three","six"), labels=c("0h","3h","6h"))
groups$replicate <- factor(groups$replicate)

# Run the actual PCA calculation
pca <- prcomp(t(pcadata), scale.=TRUE)

# Create merged metadata/data file for use with 3D Plots
if (typeOfPCA == "3D") {
  
  #shapes <- c(16,15,18,19,17,17) # Identify shape; NOTE: see http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  
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
      
      
      groups$order <- factor(groups[[metadataShape]], 
                             levels=metadataShapeFactors,
                             labels=order)
      
      plot <- autoplot(pca, x=PCX, y=PCY, data=groups)+ # What data to use for the chart
        geom_point(aes(colour=groups[[metadataColor]], shape=groups$order, size=pointSize))+ # What style of graph to make
        #stat_ellipse(mapping=aes(color=groups[[metadataColor]]),type="t")+
        geom_dl(aes(label=NA,colour=groups[[metadataColor]]),method="smart.grid")+ # How to color the chart
        labs(color = metadataColor)+ # Label for viability legend
        labs(shape = metadataShape)+ # Label for timepoint legend
        #ggtitle(paste("PCA - All ", samplesToGraph, " Samples", sep=""))+ # Title of chart
        scale_color_manual(values=colors)+ # Manually set the colors of the nodes, sequentially
        theme_bw()+ # Change background theme
        scale_shape_manual(values=shapes)+
        # To remove the legend:
        theme(legend.position="none") +
        # To theme the tick labels (0.0, etc):
        theme(axis.text.x = element_text(face="plain", color="black", 
                                         size=15, angle=0),
              axis.text.y = element_text(face="plain", color="black", 
                                         size=15, angle=0)) +
        # To remove the axis labels (PC#):
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
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
      plot <- autoplot(pca, x=PCX, y=PCY, data = groups, colour = merged3D$color, shape = FALSE, label.size = 7)+
        labs(color = metadataColor)+
        ggtitle(paste("PCA - All ", samplesToGraph, " Samples (Labeled)", sep=""))
    }
  }
}

#Print out the selected plot
plot


################################
# John's Code: top % genes

#pca <- prcomp(t(cpmdata), scale.=TRUE) # What is cpmdata?
pca <- prcomp(t(pcadata), scale.=TRUE)
gr <- factor(row.names(groups))
contribution=pca$rotation

##these add up to the total variation attributed to PC2 (or whatever PC you pick)
PCcontribute=contribution[,"PC2"]

##These are the percent of the total PC2 variation that each gene contributes.
percents=(PCcontribute*PCcontribute)

##Top of the list is the highest contributor to PC2 variation
percents=percents[order(percents,decreasing=T)]
