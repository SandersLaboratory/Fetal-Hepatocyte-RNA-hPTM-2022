###########################################################
# Volcano_Overlay.R
#
# Generate Volcano Plots from RNA-Seq Data.
#
# Code by Anders Ohman and John Santiago
# Sanders Lab, Brown University, Pathobiology
# Spring 2020 - 2022
###########################################################

########################
# User Settings

subfolder <- "AF_Genewiz"

## Which data will make the volcano?

#inputCountdata <- "AF_Merged_raw_counts.csv"
#inputMetadata <- "AF_Merged_metadata.csv"
inputCountdata <- "AF_Cells_raw_counts.csv"
inputMetadata <- "AF_Cells_metadata.csv"
#inputCountdata <- "AF_LCM_all_counts_no246.csv"
#inputMetadata <- "AF_LCM_metadata_no246.csv"

## String of genes to overlay on volcano
## Set to "LIST" to use InputList, or "GO" for GO Term

calloutList <- "LIST"
inputListUP <- "SigUp_DEAnalysis_adultfetal_Fetal-Adult.csv"
inputListDOWN <- "SigDown_DEAnalysis_adultfetal_Fetal-Adult.csv"
#inputListUP <- "SigUp_DEAnalysis_group_Dual-Adult.csv"
#inputListDOWN <- "SigDown_DEAnalysis_group_Dual-Adult.csv"
#inputListUP <- "SigUp_DEAnalysis_group_Single-Adult.csv"
#inputListDOWN <- "SigDown_DEAnalysis_group_Single-Adult.csv"
#inputList <- "SH_Genes.csv"
#inputList <- "Dual_Exclusive_Upregulated_Genes.csv"

listType <- "SigUpandDown_LCM"

#nameOfPathway <- "genesymbol"
nameOfPathway <- "ENSEMBL"
#nameOfPathway <- "SYMBOL"

topNum <- 50

# Draw connecting lines?
drawBool <- FALSE

removeThese <- c("", " ")

comparisonGroup <- "group"
inputContrast <- "LAP-Adult"
#inputContrast <- "Dual-Adult"
#inputContrast <- "Single-Adult"
#comparisonGroup <- "adultfetal"
#inputContrast <- "Fetal-Adult"

# Volcano style: Overlay
style = "Overlay"

# What FDR cutoff? (Default is 0.05)
userFDR = 0.05

  # What Log2FC cutoff? (Default is 2.0)
  userLog2fc = 0.0
  
  # Human ('org.Hs.eg') or Rat ('org.Rn.eg') genome?
  genome <- 'org.Rn.eg.db'

  # Use RedBlue or Custom coloring?
  colorStyle="RedGreen"
  
  # Axis Limits
  xLimits <- c(-11,11)
  yLimits = c(0,150)
  #xLimits <- c(-6,6)
  #yLimits = c(0,6)

# What axis label size? (Default is 18)
axisLab = 25

########################
# Load Libraries
library(tidyverse)
library(edgeR)
library(annotate)
library(EnhancedVolcano)
library(ggfortify)
library(org.Hs.eg.db)
library(org.Rn.eg.db)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/Volcano", "")
rnawd <- paste(githubwd,"Data Files/rna-seq/",sep="")
userwd <- paste(rnawd,subfolder,"/",sep="")

########################
# Load and Process Input List
if (calloutList == "LIST") {
  setwd(paste(rnawd,"DEGs",sep=""))
  pathwayListUP <- read.csv(paste(userwd,inputListUP,sep=""))
  pathwayListDOWN <- read.csv(paste(userwd,inputListDOWN,sep=""))
  pathwayListMERGE <- merge(pathwayListUP, pathwayListDOWN, all=TRUE)
  pathwayList <- pathwayListMERGE[,nameOfPathway]
  #pathwayList <- pathwayList[pathwayList != ""]
  pathwayList <- pathwayList[!pathwayList %in% removeThese]
  calloutList <- pathwayList
}
setwd(userwd)

########################
# Load Gene Sets
countdata <- read.csv(inputCountdata,row.names=1)
#countdata <- dplyr::select(countdata, -X)
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

##To run an edger comparison
z = estimateGLMCommonDisp(z,design, verbose=F)
z = estimateGLMTrendedDisp(z,design)
z = estimateGLMTagwiseDisp(z,design)
fit <- glmFit(z, design)

########################
# Generate Contrasts based on input
compare <- makeContrasts(inputContrast,levels=design)

lrt <- glmLRT(fit,contrast=as.vector(compare))
G_X_E<-topTags(lrt, n=nrow(cpmsc),adjust.method="BH", sort.by="PValue")
#hist(G_X_E$table$PValue, breaks=100,main=(G_X_E$comparison))
deanalysis=G_X_E$table

# Annotate with Gene Symbol
deanalysis$ENSEMBL <- row.names(deanalysis)
deanalysis$SYMBOL <- mapIds(get(genome),
                                keys=deanalysis$ENSEMBL,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")

# Set all NAs (no entrez mapped to ensembl) to " " (blank) before continuing.
#is.na(deanalysis$genesymbol)
deanalysis$SYMBOL[is.na(deanalysis$SYMBOL)] <- " "

# Keep only q-values of <=0.05
sigde=deanalysis[deanalysis[,"FDR"]<=.05,]

# # If looking for top #
# if (calloutList == "TOP") {
#   sigdeCOPY <- sigde
#   #sigdeCOPY <- sigdeCOPY %>% arrange(desc(sigdeCOPY$logFC)) # For logFC
#   symbolList <- sigdeCOPY$SYMBOL[1:topNum]
#   ensemblList <- sigdeCOPY$ENSEMBL[1:topNum]
#   symbolList <- symbolList[!symbolList %in% removeThese]
#   ensemblList <- ensemblList[!ensemblList %in% removeThese]
#   
#   pathwayList <- ensemblList
#   #pathwayList <- pathwayList[!pathwayList %in% removeThese]
#   calloutList <- pathwayList
# }

# Remove from name list if FDR >0.05
if (nameOfPathway=="ENSEMBL") {
  ensemblList <- calloutList[calloutList %in% sigde$ENSEMBL]
  symbolList <- mapIds(get(genome),
                        keys=ensemblList,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
  calloutList <- ensemblList
} else if (nameOfPathway=="SYMBOL"){
  symbolList <- calloutList[calloutList %in% sigde$SYMBOL]
  calloutList <- symbolList
}

# Remove from name list if logFC is between -1 and +1
sigde <- sigde[sigde$logFC > userLog2fc | sigde$logFC < -userLog2fc, ]         # Multiple conditions

if (nameOfPathway=="ENSEMBL") {
  calloutList <- calloutList[calloutList %in% sigde$ENSEMBL]
} else if (nameOfPathway=="SYMBOL"){
  calloutList <- calloutList[calloutList %in% sigde$SYMBOL]
}

#Minimal volcano plot (red and black)
if (style=="Overlay"){

  # Using red for upregulation and blue for downregulation
  if (colorStyle=="RedBlue"){
    
    #Copy in Colony logFC into table
    pathwayListMERGE_ClogFC <- pathwayListMERGE[,c("ENSEMBL","logFC")]
    colnames(pathwayListMERGE_ClogFC) <- c("ENSEMBL","ClogFC")
    deanalysis <- left_join(deanalysis, pathwayListMERGE_ClogFC, by="ENSEMBL")
    
    deanalysis$COLOR <- ifelse(
      deanalysis$FDR > userFDR, 'black',
      #ifelse(deanalysis$logFC < 0 & deanalysis$SYMBOL %in% calloutList, 'royalblue',
      #ifelse(deanalysis$logFC > 0 & deanalysis$SYMBOL %in% calloutList, 'red',
      ifelse(deanalysis$ClogFC < 0 & deanalysis$ENSEMBL %in% calloutList, 'royalblue',
      ifelse(deanalysis$ClogFC > 0 & deanalysis$ENSEMBL %in% calloutList, 'red',      
             'gray')))
    deanalysis$COLOR[is.na(deanalysis$COLOR)] <- 'gray'
    
    #Ugh. Now reorder by that so the colored dots go on top...
    
    deanalysis <- deanalysis %>% arrange(deanalysis$COLOR)
    keyvals.colour2 <- deanalysis$COLOR
    names(keyvals.colour2)[keyvals.colour2 == 'gray'] <- 'null'
    names(keyvals.colour2)[keyvals.colour2 == 'red'] <- 'up'
    names(keyvals.colour2)[keyvals.colour2 == 'black'] <- 'ns'
    names(keyvals.colour2)[keyvals.colour2 == 'royalblue'] <- 'down'
  }
  
  if (colorStyle=="RedGreen"){
    
    #Copy in Colony logFC into table
    pathwayListMERGE_ClogFC <- pathwayListMERGE[,c("ENSEMBL","logFC")]
    colnames(pathwayListMERGE_ClogFC) <- c("ENSEMBL","ClogFC")
    deanalysis <- left_join(deanalysis, pathwayListMERGE_ClogFC, by="ENSEMBL")
    
    deanalysis$COLOR <- ifelse(
      deanalysis$FDR > userFDR, 'black',
      #ifelse(deanalysis$logFC < 0 & deanalysis$SYMBOL %in% calloutList, 'royalblue',
      #ifelse(deanalysis$logFC > 0 & deanalysis$SYMBOL %in% calloutList, 'red',
      ifelse(deanalysis$ClogFC < 0 & deanalysis$ENSEMBL %in% calloutList, 'green4',
             ifelse(deanalysis$ClogFC > 0 & deanalysis$ENSEMBL %in% calloutList, 'red',      
                    'gray')))
    deanalysis$COLOR[is.na(deanalysis$COLOR)] <- 'gray'
    
    #Ugh. Now reorder by that so the colored dots go on top...
    
    deanalysis <- deanalysis %>% arrange(deanalysis$COLOR)
    keyvals.colour2 <- deanalysis$COLOR
    names(keyvals.colour2)[keyvals.colour2 == 'gray'] <- 'null'
    names(keyvals.colour2)[keyvals.colour2 == 'red'] <- 'up'
    names(keyvals.colour2)[keyvals.colour2 == 'black'] <- 'ns'
    names(keyvals.colour2)[keyvals.colour2 == 'green4'] <- 'down'
  }
  
  sigs <- deanalysis[deanalysis$FDR < 0.05,]
  
  plot <- EnhancedVolcano(deanalysis,
                                #lab = deanalysis$SYMBOL,
                                lab = NA,
                                #selectLab = symbolList,
                                x = 'logFC',
                                y = 'FDR',
                                xlim = xLimits,
                                ylim = yLimits,
                                #title = NULL, #NULLs remove the section
                                title = paste(inputContrast,": ",length(calloutList)," ",listType," overlap items in ",nrow(sigs)," sig total.",sep=""),
                                subtitle = NULL,
                                caption = NULL,
                                cutoffLineType = 'dashed',
                                #xlab = NULL,
                                #ylab = NULL,
                                xlab = bquote(~Log[2]~ 'fold change'),
                                ylab = bquote(~-Log[10]~adjusted~italic(P)/FDR),
                                axisLabSize = axisLab, #default = 18
                                pCutoff = userFDR,
                                FCcutoff = userLog2fc,
                                labSize = 4,
                                #labhjust = 0.5,
                                drawConnectors = drawBool,
                                arrowheads = FALSE,
                                widthConnectors = 0.5,
                                #boxedLabels = TRUE,
                                #legendLabels =c('NS','NS','NS','Adjusted p-value (FDR) & Log2 FC'),
                                legendPosition = 'none',
                                #col = c('gray','gray','gray','gray'),
                                #colCustom = c('gray','gray'), # Manual colors
                                colCustom = keyvals.colour2,
                                colAlpha = 1,
                                border = 'full', #vs 'full'
                                borderWidth = 2,
                                borderColour = 'black')
}


# Print the selected plot to the viewer
plot

# Can't create a file with a colon, replace in GO term with underscore for now.
nameOfOutputFile <- paste(inputContrast,"_",listType,sep="")

# Save output of plot
#setwd(paste(userwd,"/Volcanoes",sep=""))
setwd(paste(userwd))
#ggsave(paste(nameOfOutputFile,".pdf",sep=""),width = 8, height = 8)
ggsave(paste(nameOfOutputFile,".png",sep=""),width = 10, height = 10)

# Write DE analysis to file
#write.csv(deanalysis,file=paste("DEAnalysis",inputContrast,".csv",sep="_"))

# Write post-filter callout list to CSV
calloutList <- as.data.frame(calloutList,symbolList)
colnames(calloutList) <- nameOfPathway
write.csv(calloutList, paste(nameOfOutputFile,"_filtered_genes.csv",sep=""))

#join deanalysis(logFC and FDR) by ENSEMBL
df = merge(x=calloutList,y=deanalysis[ , c("SYMBOL","ENSEMBL","logFC","FDR")],by="ENSEMBL",all.x=TRUE)
names(df)[names(df) == "logFC"] <- paste(inputContrast,"logFC",sep="_")
names(df)[names(df) == "FDR"] <- paste(inputContrast,"FDR",sep="_")
#join inputListUP (logFC and FDR) by ENSEMBL
df = merge(x=df,y=pathwayListUP[ , c("ENSEMBL","logFC","FDR")],by="ENSEMBL",all.x=TRUE)
names(df)[names(df) == "logFC"] <- paste("Colony_Host","logFC",sep="_")
names(df)[names(df) == "FDR"] <- paste("Colony_Host","FDR",sep="_")
#join inputlistDOWN (logFC and FDR) by ENSEMBL
df = merge(x=df,y=pathwayListDOWN[ , c("ENSEMBL","logFC","FDR")],by="ENSEMBL",all.x=TRUE)
names(df)[names(df) == "logFC"] <- paste("Colony_Host","logFC",sep="_")
names(df)[names(df) == "FDR"] <- paste("Colony_Host","FDR",sep="_")

write.csv(df, paste(nameOfOutputFile,"_overlap_gene_expression.csv",sep=""))

# blueDown <- sum(deanalysis$logFC<0 & deanalysis$COLOR=='royalblue')
# blueUp <- sum(deanalysis$logFC>0 & deanalysis$COLOR=='royalblue')
# redDown <- sum(deanalysis$logFC<0 & deanalysis$COLOR=='red')
# redUp <- sum(deanalysis$logFC>0 & deanalysis$COLOR=='red')
