###########################################################
# Heatmap_hPTM.R
#
# Generates Heatmaps, Dotplots, and Bar Graphs from our
# hPTM data.
#
# Code by Anders Ohman
# Sanders Lab, Brown University, Pathobiology
# Spring 2020 - 2022
###########################################################

########################
# User Settings
########################

filename <- "hPTM Merged.csv"
names <- c("Colonies","Colonies.q","Liver","Liver.q","LAP+OC2+","LAP+OC2+.q","LAP+OC2-","LAP+OC2-.q","LAP+","LAP+.q")
significance <- 0.05

# Drop Liver column?
dropLiver = FALSE

########################
# Load Libraries
########################
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggpattern)
library(plotly)
library(heatmaply)

########################
# Automatically set working directory
currentwd <- dirname(rstudioapi::getSourceEditorContext()$path)
githubwd <- str_replace(currentwd, "R Scripts/Heatmap", "")
protwd = paste(githubwd,"Data Files/proteomics/hPTM",sep="")
heatmapwd = paste(githubwd,"Data Files/proteomics/Heatmaps",sep="")
setwd(protwd)

df <- read.csv(filename, row.names=1)
colnames(df) <- names

########################
# Heatmap Dot

df1 <- df
df1 <- df1[complete.cases(df1), ]
df1 <- df1[!(df1$Colonies.q > significance | df1$"LAP+OC2+.q" > significance | df1$"LAP+OC2-.q" > significance),]
dfq <- df1 %>% dplyr::select(contains("q"))

# Change Liver NSes to NA
#df1[(df1$Liver.q > significance),]$Liver <- NA

dfc <- df1 %>% dplyr::select(-contains("q"))

# Drop Liver column
if(dropLiver==TRUE){
  dfc <- dfc %>% dplyr::select(-contains("Liver"))
  dfq <- dfq %>% dplyr::select(-contains("Liver"))
}

# Cell Note the non-sigs
cn <- as.data.frame(dfq)
#cn <- round(cn, digits = 3)
cn[cn >= 0.05] <- "ns"
dat <- as.data.frame(sapply(cn, as.numeric)) #<- sapply is here
dat[dat <= 0.05] <- ""
dat[is.na(dat)] <- "ns"

setwd(heatmapwd)

heatmaply_cor(dfc, 
              dendrogram='row', 
              #dendrogram='none',
              #seriate = "mean",
              limits = c(-32, 32),
              node_type = "scatter",
              point_size_mat = -log10(dfq), 
              point_size_name = "-log10(q)",
              label_names = c("hPTM", "Group", "%"),
              cellnote = dat,
              #scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colours = c("blue", "white", "red"),
              #                                                        oob=scales::squish,
              #                                                        na.value = "white",
              #                                                        values = scales::rescale(c(-1, -0.05, 0, 0.05, 1))),
              file = paste("hPTM_heatmaply_cor.pdf", sep=""), width=600,height=500)

#write.csv(dfc, "hPTM_Significant_13_Marks.csv")

setwd(protwd)

########################
# Heatmap Full Dot

df1 <- df
df1 <- df1[complete.cases(df1), ]
#df1 <- df1[!(df1$Colonies.q > significance | df1$"LAP+OC2+.q" > significance | df1$"LAP+OC2-.q" > significance),]
dfc <- df1 %>% select(-contains("q"))
dfq <- df1 %>% select(contains("q"))

cn <- as.data.frame(dfq)
#cn <- round(cn, digits = 3)
cn[cn >= 0.1 | is.na(cn)] <- ""
dat <- as.data.frame(sapply(cn, as.numeric)) #<- sapply is here
dat[dat <= 0.05] <- "*"
dat[dat <= 0.1 & dat > 0.05] <- "^"
dat[is.na(dat)] <- ""

setwd(heatmapwd)

heatmaply_cor(dfc, dendrogram='row', limits = c(-35, 35),
  node_type = "scatter",
  point_size_mat = -log10(dfq), 
  point_size_name = "-log10(q)",
  label_names = c("hPTM", "Group", "% Change"),
  cellnote = dat,
  #scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colours = c("blue", "white", "red"),
  #                                                        oob=scales::squish,
  #                                                        na.value = "white",
  #                                                        values = scales::rescale(c(-1, -0.05, 0, 0.05, 1))),
  file = paste("hPTM_heatmaply_cor_full.png", sep=""), width=500,height=1200,
  
)

setwd(protwd)


########################
# Heatmap Regular

####workaround for Error in hclustfun(dist) : NA/NaN/Inf in foreign function call (arg 10)
df2 <- df
#df2[is.na(df2)] <- 0
df2 <- df2[complete.cases(df2), ]
df2 <- df2[!(df2$Colonies.q > significance | df2$"LAP+OC2+.q" > significance | df2$"LAP+OC2-.q" > significance),] 
#df2 <- df2[!(df2$Colonies.q > significance),] 
df2 <- df2 %>% select(-contains("q"))

setwd(heatmapwd)

heatmaply(df2, dendrogram='row', main="",
          showticklabels = c(TRUE,TRUE),
          show_dendrogram = c(TRUE,TRUE),
          #hide_colorbar = TRUE,
          #symm = TRUE,
          #point_size_mat = df,
          key.title="% Fetal Diff",
          label_names = c("hPTM", "Group", "logFC"), xlab = NULL, ylab = NULL, fontsize_row = NULL, fontsize_col = 15,
          scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(colours = c("blue", "white", "red"),limits=c(-30,30), oob=scales::squish,
                                                                  na.value = "white",
                                                                  values = scales::rescale(c(-1, -0.05, 0, 0.05, 1))),
          file = paste("hPTM_heatmaply.png", sep=""), width=500,height=1200,
)

graphics.off()
setwd(protwd)

########################
# Bar graphs

df <- read.csv("hPTM Significant.csv")
setDT(df)
dfMelt <- melt(df, id.vars = "Sample")

dfMelt$Sample <- factor(dfMelt$Sample,
                        levels=c("ADH","FEHS","FEHD",
                                 "A-19","F-19",
                                 "A-20","F-20"),
                        labels=c("Adult Hepatocytes","Fetal Hepatocytes (S)","Fetal Hepatocytes (D)",
                                 "Adult Liver","Fetal Liver",
                                 "Adult Host","Fetal Colonies"))

agg <- dfMelt[, .(mean = mean(value), 
                  se = sd(value)/.N), 
              by = .(Sample, variable)]

agg <- arrange(agg, variable, Sample)

########################
# Plot by mark

agg$variable <- factor(agg$variable,
                       levels=c("H3.1_K27UN",
                                "H3.3_K27UN",
                                "H3.3_K36ME2",
                                "H3.1_K27ME1",
                                "H3_K9UN",
                                "H3.1_K36UN",
                                "H3_K18ME1",
                                "H3.3_K27ME3",
                                "H3.3_K27ME2",
                                "H3.1_K27ME2",
                                "H3.3_K36ME1",
                                "H3.1_K27ME3",
                                "H3.1_K36ME1"))

plot <- ggplot(agg, 
       aes(x = variable, y = mean, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Histone Modifications") +
  # Add error bars (here +/- 1.96 SE)
  geom_errorbar(aes(ymax = mean + 1.96*se, 
                    ymin = mean - 1.96*se),
                position = "dodge") +
  xlab(NULL) + 
  ylab("% Relative Abundance") +
  labs(fill = "") +
  theme_light() +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(face="bold"))

plot

# Plot by mark (B&W)

agg$color[agg$Sample == "Adult Hepatocytes"] <- "black"
agg$color[agg$Sample == "Fetal Hepatocytes (S)"] <- "white"
agg$color[agg$Sample == "Fetal Hepatocytes (D)"] <- "white"
agg$color[agg$Sample == "Adult Liver"] <- "royalblue"
agg$color[agg$Sample == "Fetal Liver"] <- "skyblue"
agg$color[agg$Sample == "Adult Host"] <- "darkred"
agg$color[agg$Sample == "Fetal Colonies"] <- "darksalmon"

ggplot(agg, 
       aes(x = variable, y = mean, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  #ggtitle("Histone Modifications") +
  # Add error bars (here +/- 1.96 SE)
  geom_errorbar(aes(ymax = mean + 1.96*se, 
                    ymin = mean - 1.96*se),
                position = "dodge") +
  xlab(NULL) + 
  ylab("% Relative Abundance") +
  labs(fill = "") +
  #theme_light() +
  scale_fill_manual(values = agg$color)+
  theme(axis.text.x = element_text(face="bold"))

# Plot by mark (pattern)

agg$color[agg$Sample == "Adult Hepatocytes"] <- "black"
agg$color[agg$Sample == "Fetal Hepatocytes (S)"] <- "red3"
agg$color[agg$Sample == "Fetal Hepatocytes (D)"] <- "red3"
agg$color[agg$Sample == "Adult Liver"] <- "gray0"
agg$color[agg$Sample == "Fetal Liver"] <- "red3"
agg$color[agg$Sample == "Adult Host"] <- "gray1"
agg$color[agg$Sample == "Fetal Colonies"] <- "red3"

agg$pattern[agg$Sample == "Adult Hepatocytes"] <- "none"
agg$color[agg$Sample == "Fetal Hepatocytes (S)"] <- "none"
agg$color[agg$Sample == "Fetal Hepatocytes (D)"] <- "none"
agg$pattern[agg$Sample == "Adult Liver"] <- "stripe"
agg$pattern[agg$Sample == "Fetal Liver"] <- "stripe"
agg$pattern[agg$Sample == "Adult Host"] <- "crosshatch"
agg$pattern[agg$Sample == "Fetal Colonies"] <- "crosshatch"

agg$patternCol[agg$Sample == "Adult Hepatocytes"] <- "black"
agg$color[agg$Sample == "Fetal Hepatocytes (S)"] <- "black"
agg$color[agg$Sample == "Fetal Hepatocytes (D)"] <- "black"
agg$patternCol[agg$Sample == "Adult Liver"] <- "white"
agg$patternCol[agg$Sample == "Fetal Liver"] <- "black"
agg$patternCol[agg$Sample == "Adult Host"] <- "white"
agg$patternCol[agg$Sample == "Fetal Colonies"] <- "black"

# plot <- ggplot(agg, 
#        aes(x = variable, y = mean, fill = Sample)) + 
#   geom_bar(stat = "identity", position = "dodge", color = "black") +
#   ggtitle("Histone H3.1 Modifications") +
#   # Add error bars (here +/- 1.96 SE)
#   geom_errorbar(aes(ymax = mean + 1.96*se, 
#                     ymin = mean - 1.96*se),
#                 position = "dodge") +
#   xlab(NULL) + 
#   ylab("% Relative Abundance") +
#   labs(fill = "") +
#   #theme_light() +
#   scale_fill_manual(values = agg$color)+
#   theme(axis.text.x = element_text(face="bold"))
# 
# plot +
#   geom_bar_pattern(
#     aes(variable, mean, pattern_fill = variable), 
#     pattern = 'stripe',
#     fill    = 'white',
#     colour  = 'black'
#   )

ggplot(data = agg, aes(x = variable, y = mean, fill = Sample, pattern = Sample)) +
  geom_bar_pattern(stat='identity',
                   position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = agg$patternCol,
                   pattern_angle = 45,
                   pattern_density = 0.35,
                   pattern_spacing = 0.008,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = agg$color) +
  geom_errorbar(aes(ymax = mean + 1.96*se, 
                    ymin = mean - 1.96*se),
                position = "dodge") +
  scale_pattern_manual(values = agg$pattern)+
  xlab(NULL) + 
  ylab("% Relative Abundance")

# df2 <- df
# df2 %>% add_row(hello = "hola", goodbye = "ciao")
# 
# agg[3]
# df2$flMinusal <- agg[2,mean]-agg[1,mean] #F-A Liver
#agg[4,mean]-agg[3,mean] #F-A Colonies
# 
# # Plot by sample type
# ggplot(agg, 
#        aes(x = Sample, y = mean, fill = variable)) + 
#   geom_bar(stat = "identity", position = "dodge") +
#   ggtitle("Histone H3.1 Modifications") +
#   # Add error bars (here +/- 1.96 SE)
#   geom_errorbar(aes(ymax = mean + 1.96*se, 
#                     ymin = mean - 1.96*se),
#                 position = "dodge") +
#   xlab(NULL) + 
#   ylab("% Relative Abundance") +
#   labs(fill = "") +
#   theme_light() +
#   scale_fill_brewer(palette="Set1") +
#   theme(axis.text.x = element_text(face="bold"))
# 
# # Plot by
# h3 <- read.csv("H3 hPTMs.csv", row.names = 1)
# colnames(h3) <- c("Colonies", "Liver", "LAP+OC2+","LAP+OC2-")
# h3 <- t(h3)
# g <- ggplot(h3, aes(row.names(h3))) + geom_bar()

########################
# PLOT ONLY ONE MARK

df <- read.csv("hPTM Significant.csv")
setDT(df)
dfMelt <- melt(df, id.vars = "Sample")

dfMelt$Sample <- factor(dfMelt$Sample,
                        levels=c("ADH","FEHS","FEHD",
                                 "A-19","F-19",
                                 "A-20","F-20"),
                        labels=c("Adult Hepatocytes","Fetal Hepatocytes (S)","Fetal Hepatocytes (D)",
                                 "Adult Liver","Fetal Liver",
                                 "Adult Host","Fetal Colonies"))

agg <- dfMelt[, .(mean = mean(value), 
                  se = sd(value)/.N), 
              by = .(Sample, variable)]

agg <- arrange(agg, variable, Sample)

########################
# Plot by mark

agg <- as.data.frame(agg)
agg[agg[, "variable"] == "H3.1_K27UN",]
agg[agg[,"variable"] == "H3.1_K27UN",,drop=FALSE]


agg$variable <- factor(agg$variable,
                       levels=c("H3.1_K27UN",
                                "H3.3_K27UN",
                                "H3.3_K36ME2",
                                "H3.1_K27ME1",
                                "H3_K9UN",
                                "H3.1_K36UN",
                                "H3_K18ME1",
                                "H3.3_K27ME3",
                                "H3.3_K27ME2",
                                "H3.1_K27ME2",
                                "H3.3_K36ME1",
                                "H3.1_K27ME3",
                                "H3.1_K36ME1"))

plot <- ggplot(agg, 
               aes(x = variable, y = mean, fill = Sample)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Histone Modifications") +
  # Add error bars (here +/- 1.96 SE)
  geom_errorbar(aes(ymax = mean + 1.96*se, 
                    ymin = mean - 1.96*se),
                position = "dodge") +
  xlab(NULL) + 
  ylab("% Relative Abundance") +
  labs(fill = "") +
  theme_light() +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(face="bold"))

plot

