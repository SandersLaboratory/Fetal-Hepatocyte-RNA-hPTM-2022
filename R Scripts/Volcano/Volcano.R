###########################################################
# Volcano.R
#
# Generate Volcano Plots from our RNA-Seq Data.
#
# Code by Anders Ohman and John Santiago
# Sanders Lab, Brown University, Pathobiology
# Spring 2020 - 2022
###########################################################

########################
# User Settings

subfolder <- "AF_Genewiz"

# What analysis to run?
# 'Dual', 'Single', or 'LCM'
comparison <- "LCM"

if(comparison == "LCM"){
  inputCountdata <- "AF_LCM_all_counts_no246.csv"
  inputMetadata <- "AF_LCM_metadata_no246.csv"
  comparisonGroup <- "adultfetal"
  inputContrast <- "Fetal-Adult"
  xLimits <- c(-6,6)
  yLimits = c(0,6)
}

if(comparison == "Dual" || comparison == "Single"){
  inputCountdata <- "AF_Cells_raw_counts.csv"
  inputMetadata <- "AF_Cells_metadata.csv"
  comparisonGroup <- "group"
  if(comparison == "Dual"){inputContrast <- "Dual-Adult"}
  if(comparison == "Single"){inputContrast <- "Single-Adult"}
  xLimits <- c(-11,11)
  yLimits = c(0,150)
}

# Volcano style: Simple, Enhanced, or Minimal?
# Minimal is designed for figures, and has axes +/- 15 with no labels
style = "Minimal"

# What FDR cutoff? (Default is 0.05)
userFDR = 0.05

  # For Enhanced only, what Log2FC cutoff? (Default is 2.0)
  userLog2fc = 1.0
  
  # For Enhanced only, Human ('org.Hs.eg') or Rat ('org.Rn.eg') genome?
  genome <- 'org.Rn.eg'

  # For Minimal only, use Red for both +/-, or RedBlue?
  colorStyle="RedBlue"
  
  # Axis Limits
  #xLimits <- c(-20,20)
  #yLimits = c(0,20)

# What axis label size? (Default is 18)
axisLab = 25

########################
# Load Libraries
library(tidyverse)
library(edgeR)
library(annotate)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(org.Rn.eg.db)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/Volcano", "")
rnawd <- paste(githubwd,"Data Files/rna-seq/",sep="")
userwd <- paste(rnawd,subfolder,"/",sep="")
setwd(userwd)

########################
# Load Gene Sets
countdata <- read.csv(inputCountdata,row.names=1)
groups <- read.csv(inputMetadata,row.names=1)

#Saves normalized CPM file and a count table both with low count reads removed
countdata=countdata[,row.names(groups)]
x<-DGEList(counts=countdata)

# Design matrix
design<-model.matrix(~0+groups[[comparisonGroup]])
colnames(design) <- c(levels(as.factor(groups[[comparisonGroup]])))

# Trim counts
keep.exprs <- filterByExpr(countdata, design = design)
filt1 <- x[keep.exprs,,keep.lib.sizes=FALSE]
# dim(filt1)
y <- filt1

z <- calcNormFactors(y, method = "TMM")
##using z gives the normalized cpm, using y gives the not normalized cpm
cpmsc=cpm(z)

#write.csv(cpmsc,file="AF_Cells_trimmed_normalized_counts.csv")

##To run an edger comparison
z = estimateGLMCommonDisp(z,design, verbose=F)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)

# Plot the Biological Coefficient of Variation
plotBCV(z)

########################
# Generate Contrasts based on input
compare <- makeContrasts(inputContrast,levels=design)

lrt <- glmLRT(fit,contrast=as.vector(compare))
G_X_E<-topTags(lrt, n=nrow(cpmsc),adjust.method="BH", sort.by="PValue")
#hist(G_X_E$table$PValue, breaks=100,main=(G_X_E$comparison))
deanalysis=G_X_E$table

# Write DE analysis to file
#write.csv(deanalysis,file=paste("DEAnalysis",inputContrast,".csv",sep="_"))

# Keep only q-values of <=0.05
sigde=deanalysis[deanalysis[,"FDR"]<=.05,]

# Enhanced Volcano Plot
if (style=="Enhanced"){
  ENSEMBLList <- as.list(get(paste(genome,"ENSEMBL2EG",sep="")))

  #Look up Entrez IDs from row name ENSEMBL
  deanalysis$geneid <- ENSEMBLList[row.names(deanalysis)]

  #Look up gene symbol from Entrez ID
  
  # IF there are multiple entrez IDs for an ENSEMBL:
  # Set all mutli-entrez elements to the first entry in the list.
  deanalysis$geneid <- map(deanalysis$geneid, 1)
  
  deanalysis$genesymbol <- lookUp(as.character(deanalysis$geneid), genome, 'SYMBOL')
  
  # Set all NAs (no entrez mapped to ensembl) to " " (blank) before continuing.
  #is.na(deanalysis$genesymbol)
  deanalysis$genesymbol[is.na(deanalysis$genesymbol)] <- " "
  
  plot <- EnhancedVolcano(deanalysis,
                                 lab = deanalysis$genesymbol,
                                 x = 'logFC',
                                 y = 'FDR',
                                 xlim = xLimits,
                                 title = "Enhanced Plot",
                                 subtitle = paste(comparisonGroup,inputContrast,sep=" "),
                                 caption = paste('FC cutoff = ',sprintf("%.1f", round(userLog2fc,1)),', FDR cutoff = ',userFDR),
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 ylab = bquote(~-Log[10]~adjusted~italic(P)/FDR),
                                 axisLabSize = axisLab, #default = 18
                                 pCutoff = userFDR,
                                 FCcutoff = userLog2fc,
                                 labSize = 3,
                                 colAlpha = 1,
                                 #legend=c('NS','Log2 FC','Adjusted p-value (FDR)',
                                 #         'Adjusted p-value & Log2 FC'),
                                 legendPosition = 'top',
                                 legendLabSize = 10,
                                 legendIconSize = 3.0)
}

#Simple volcano plot (red and black)
if (style=="Simple"){
  plot <- EnhancedVolcano(deanalysis,
                                 lab = NA,
                                 x = 'logFC',
                                 y = 'FDR',
                                 xlim = xLimits,
                                 ylim = yLimits,
                                 title = "Simple Plot",
                                 subtitle = paste(comparisonGroup,inputContrast,sep=" "),
                                 caption = paste('FDR cutoff = ',userFDR),
                                 cutoffLineType = 'blank',
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 ylab = bquote(~-Log[10]~adjusted~italic(P)/FDR),
                                 axisLabSize = axisLab, #default = 18
                                 pCutoff = userFDR,
                                 FCcutoff = 0.0,
                                 labSize = 1.5,
                                 legendLabels =c('NS','NS','NS','Adjusted p-value (FDR) & Log2 FC'),
                                 legendPosition = 'top',
                                 col=c('black','black','black','red'),
                                 colAlpha = 1)
}

#Minimal volcano plot (red and black)
if (style=="Minimal"){

  # Using red for upregulation and blue for downregulation
  if (colorStyle=="RedBlue"){
    keyvals.colour <- ifelse(
      deanalysis$FDR > userFDR, 'black',
      ifelse(deanalysis$logFC < 0, 'royalblue',
      ifelse(deanalysis$logFC > 0, 'red',
             'black')))
    keyvals.colour[is.na(keyvals.colour)] <- 'black'
    names(keyvals.colour)[keyvals.colour == 'red'] <- 'up'
    names(keyvals.colour)[keyvals.colour == 'black'] <- 'ns'
    names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'down'
  }

  # Using red color for both up and downregulation
  if (colorStyle=="Red"){
    keyvals.colour <- ifelse(deanalysis$FDR > userFDR, 'black', 'red')
    keyvals.colour[is.na(keyvals.colour)] <- 'black'
    names(keyvals.colour)[keyvals.colour == 'red'] <- 'sig'
    names(keyvals.colour)[keyvals.colour == 'black'] <- 'ns'
  }

  plot <- EnhancedVolcano(deanalysis,
                                 lab = NA,
                                 x = 'logFC',
                                 y = 'FDR',
                                 xlim = xLimits,
                                 ylim = yLimits,
                                 title = NULL, #NULLs remove the section
                                 subtitle = NULL,
                                 caption = NULL,
                                 cutoffLineType = 'blank',
                                 xlab = NULL,
                                 ylab = NULL,
                                 axisLabSize = axisLab, #default = 18
                                 pCutoff = userFDR,
                                 FCcutoff = 0.0,
                                 labSize = 1.5,
                                 #legendLabels =c('NS','NS','NS','Adjusted p-value (FDR) & Log2 FC'),
                                 legendPosition = 'none',
                                 #colCustom = c('black','black','black','red'), # Manual colors
                                 colCustom = keyvals.colour,
                                 colAlpha = 1,
                                 border = 'full', #vs 'full'
                                 borderWidth = 2,
                                 borderColour = 'black')
}

# Print the selected plot to the viewer
plot

# Write DE analysis to file
deanalysis$ENSEMBL <- row.names(deanalysis)

# Number of Sig Fetal:
sum(deanalysis$logFC>0 & deanalysis$FDR<0.05)
sigUp <- deanalysis[(deanalysis$logFC>0 & deanalysis$FDR<0.05),]
# Number of Sig Adult:
sum(deanalysis$logFC<0 & deanalysis$FDR<0.05)
sigDown <- deanalysis[(deanalysis$logFC<0 & deanalysis$FDR<0.05),]

deanalysis2 <- apply(deanalysis,2,as.character)
write.csv(deanalysis2,file=paste("DEAnalysis_",comparisonGroup,"_",inputContrast,".csv",sep=""))
#ggsave(paste(style,"_",comparisonGroup,"_",inputContrast,".png",sep=""),width = 8, height = 8)

# Need to convert list element to character to save correctly, but have to do after the > and < statements
sigUp <- apply(sigUp,2,as.character)
sigDown <- apply(sigDown,2,as.character)

write.csv(sigUp,file=paste("SigUp_","DEAnalysis_",comparisonGroup,"_",inputContrast,".csv",sep=""))
write.csv(sigDown,file=paste("SigDown_","DEAnalysis_",comparisonGroup,"_",inputContrast,".csv",sep=""))

# ###################
# # Prune the master count data list based on what is significant here
# 
# #Describe how you're filtering
# Normalized <- TRUE
# 
# #Read in both files
# if (!Normalized){
#   countdata <- read.csv(inputCountdata, row.names=1) #Raw
#   datatype <- "Raw"
# }
# if (Normalized){
#   countdata <- read.csv(inputCountdata, row.names=1) #Raw
#   countdata=countdata[,row.names(groups)]
#   x<-DGEList(countdata)
#   class(x)
#   geneid <- rownames(x)
#   atleast=(10000000/min(x$samples$lib.size))
#   keep=rowSums(cpm(x)>=atleast)>=3
#   y=x[keep, ,keep.lib.sizes=FALSE]
#   z <- calcNormFactors(y, method = "TMM")
#   countdata=cpm(z)
#   datatype <- "Normalized"
# }
# 
# filtereddata <- read.csv("DEAnalysis_Fetal-Adult.csv")
# exclusion <- "IncludingA246"
# 
# Exclude246 <- FALSE
# if (Exclude246){
#   filtereddata <- read.csv("DEAnalysis_Fetal-Adult_withoutA246.csv")
#   exclusion <- "ExcludingA246"
# }
# 
# #Filter data by FDR and then by that ENSEMBL ID list
# filtereddata <- filtereddata[!(filtereddata$FDR >0.05),]
# sigRows <- filtereddata$ENSEMBL
# filteredcountdata <- countdata[ row.names(countdata) %in% sigRows, ]
# 
# write.csv(filteredcountdata,file=paste("Filtered_CountData",datatype,exclusion,inputContrast,".csv",sep="_"))