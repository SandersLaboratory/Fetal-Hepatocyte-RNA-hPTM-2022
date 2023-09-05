########################
# Load Libraries
library(VennDiagram)
library(tidyverse)
library(heatmaply)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(annotate)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/Venn", "")
vennwd = paste(githubwd,"Data Files/venn","/",sep="")
rnawd <- paste(githubwd,"Data Files/rna-seq/AF_Genewiz/",sep="")

# Import RNA-Seq DE Dataset
setwd(rnawd)
colonyIn <- read.csv("DEAnalysis_adultfetal_Fetal-Adult.csv", row.names="ENSEMBL")
dualIn <- read.csv("DEAnalysis_group_Dual-Adult.csv", row.names="ENSEMBL")
singleIn <- read.csv("DEAnalysis_group_Single-Adult.csv", row.names="ENSEMBL")


# Annotate with Gene Symbol
colonyIn$ENSEMBL <- row.names(colonyIn)
colonyIn$SYMBOL <- mapIds(org.Rn.eg.db,
                            keys=colonyIn$ENSEMBL,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
colonyIn$UNIPROT <- mapIds(org.Rn.eg.db,
                          keys=colonyIn$ENSEMBL,
                          column="UNIPROT",
                          keytype="ENSEMBL",
                          multiVals="first")

dualIn$ENSEMBL <- row.names(dualIn)
dualIn$SYMBOL <- mapIds(org.Rn.eg.db,
                          keys=dualIn$ENSEMBL,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
dualIn$UNIPROT <- mapIds(org.Rn.eg.db,
                           keys=dualIn$ENSEMBL,
                           column="UNIPROT",
                           keytype="ENSEMBL",
                           multiVals="first")

singleIn$ENSEMBL <- row.names(singleIn)
singleIn$SYMBOL <- mapIds(org.Rn.eg.db,
                        keys=singleIn$ENSEMBL,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
singleIn$UNIPROT <- mapIds(org.Rn.eg.db,
                         keys=singleIn$ENSEMBL,
                         column="UNIPROT",
                         keytype="ENSEMBL",
                         multiVals="first")

# Reduce all to only Signif q<0.05
dual <- dualIn[!(dualIn$FDR > 0.05),]
single <- singleIn[!(singleIn$FDR > 0.05),]
colony <- colonyIn[!(colonyIn$FDR > 0.05),]

# Count and Purge NA names
sum(is.na(dual$SYMBOL))
sum(is.na(single$SYMBOL))
sum(is.na(colony$SYMBOL))
#sum(colony$genesymbol==" ")

dual <- dual[!is.na(dual$SYMBOL),]
single <- single[!is.na(single$SYMBOL),]
colony <- colony[!is.na(colony$SYMBOL),]
#colony <- colony[!colony$genesymbol==" ",]

# Export only gene names from each list
dual <- dual$SYMBOL
single <- single$SYMBOL
colony <- colony$SYMBOL
#colony <- colony$genesymbol

# sums
length(dual)
length(single)
length(colony)

# Venn!
set1="dual"
set2="single"
set3="colony"
s1=dual
s2=single
s3=colony
grid.newpage()
draw.triple.venn(area1=length(s1),area2=length(s2),area3=length(s3),n12=length(intersect(s1,s2)),n23 = length(intersect(s2,s3)),n13 = length(intersect(s1,s3)),n123 = length(intersect(intersect(s1,s2),s3)),category=c(set1,set2,set3),fill = c("green", "blue","red"),cex=1.5,cat.cex = 1)

all3genes <- intersect(intersect(s1,s2),s3)
dualandcolonygenes <- intersect(s1,s3)
singleandcolony <- intersect(s2,s3)

setwd(vennwd)
write.csv(all3genes, "Genes_In_Dual&Single&Colony.csv")
write.csv(dualandcolonygenes, "Genes_In_Dual&Colony.csv")
write.csv(singleandcolony, "Genes_In_Single&Colony.csv")

dualcolonyonly <- dualandcolonygenes [! dualandcolonygenes %in% singleandcolony]
write.csv(dualcolonyonly, "Genes_Only_In_Dual&Colony_NotSingle.csv")

set1="dual"
set2="colony"
s1=dual
s2=colony
grid.newpage()
draw.pairwise.venn(area1=length(s1),area2=length(s2),cross.area=length(intersect(s1,s2)),fill = c("blue","red"),cex=1.5,cat.cex = 1,category = c(set1,set2))
