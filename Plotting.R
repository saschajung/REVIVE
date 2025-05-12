library(ggplot2)
library(viridis)
library(ggalluvial)
library(treemapify)
library(ggbreak)
library(h2o)

setwd("~/Documents/RejuvenationPerturbationPaper/PlottingData/")

my_theme <- function(){
  theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = rel(1.25), color = "black", family="Arial"),
          legend.title = element_text(size=rel(1.25), color = "black", family="Arial"),
          legend.text = element_text(size=rel(1.25), color = "black", family="Arial"),
          axis.title = element_text(size=rel(2), color = "black", family="Arial"),
          axis.title.y = element_text(margin = margin(r = 10)),
          axis.text.x = element_text(size=rel(1.5), angle = 0, hjust = 0.5, vjust=0.8, color = "black", family="Arial"),
          axis.text.y = element_text(size=rel(2), hjust = 1, vjust=0.5, color = "black", family="Arial"),
          strip.text.y = element_text(size=rel(2), color = "black", family="Arial"),
          strip.text.x = element_text(size=rel(2), color = "black", family="Arial"),
          strip.background = element_rect(fill="white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
}

load("~/Documents/RejuvenationPerturbationPaper/PlottingData/mbotc.rda",verbose = T)
symbToEns <- readRDS("~/Documents/RejuvenationPerturbationPaper/PlottingData/mbotc_symbToEns.rds")

mbotc <- merge(mbotc,symbToEns,by.x = "Symbol",by.y = "hgnc_symbol")
mbotc <- mbotc[mbotc$ProcessLevel == 1,]

h2o.init()
mod <- h2o.loadModel("./Documents/GLM_model_R_1708656074038_1")

predictors <- mod@model[["coefficients_table"]]
predictors <- predictors[predictors$coefficients != 0,]

#Check although it's clear - Do not check intercept
all(predictors$names[-1] %in% mbotc$ensembl_gene_id)

mbotc <- mbotc[mbotc$ensembl_gene_id %in% predictors$names,]
mbotc <- mbotc[,c("ensembl_gene_id","Symbol","ProcessName")]

suptab <- aggregate(mbotc$ProcessName,by = list(EnsemblID = mbotc$ensembl_gene_id,Symbol = mbotc$Symbol), FUN = paste0,collapse = "; ")

###Supplementary Table S1
write.table(suptab, file = "SupplementaryTableS1_Input.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#Supplementary Fig. 1a
genesPerProc <- as.data.frame(table(mbotc$ProcessName))
genesPerProc <- genesPerProc[order(genesPerProc$Freq,decreasing = T),]
genesPerProc$Var1 <- factor(genesPerProc$Var1,levels = genesPerProc$Var1)

p <- ggplot(genesPerProc, aes(x = Var1, y = Freq, fill= Var1)) +
  geom_bar(stat="identity",show.legend = F) +
  scale_fill_viridis(option="plasma", name="Process", discrete=T) +
  my_theme() +
  ylab("Number of Genes") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("SupplementaryFig1a.png",p,dpi = 300, width = 25, height = 10)

mbotc_wScores <- merge(mbotc,mod@model$coefficients_table,by.x = "ensembl_gene_id",by.y = "names")

aggScores <- aggregate(abs(mbotc_wScores$standardized_coefficients),by = list(Process = mbotc_wScores$ProcessName),FUN = sum)
aggScores$Process <- factor(aggScores$Process, levels = genesPerProc$Var1)

p <- ggplot(aggScores, aes(x = Process, y = x, fill= Process)) +
  geom_bar(stat="identity",show.legend = F) +
  scale_fill_viridis(option="plasma", name="Process", discrete=T) +
  my_theme() +
  ylab("Process Importance") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("SupplementaryFig1b.png",p,dpi = 300, width = 25, height = 10)

cor(genesPerProc$Freq,aggScores$x[match(genesPerProc$Var1,aggScores$Process)], method = "spearman") #0.9403069

###### Fig. 2a ######
inp <- data.frame(Platform = c("LINCS","Other"),
                 Value = c(228681/(1734+228681)*100,1734/(1734+228681)*100))
inp$Platform <- factor(inp$Platform,levels = (inp$Platform))

p <- ggplot(inp, aes(x = '', y = Value, fill= Platform)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_viridis(option="plasma", begin = 0.15, end = 0.85,name="Platform", discrete=T) +
  theme_void() +
  my_theme() +
  ylab("") +
  xlab("")

ggsave("Fig2a.png",p,dpi = 300, width = 7, height = 7)

###### Fig. 2b ######
mappingofcts <- read.table("../MappingOfCTs.txt",header = F, sep = "\t")

cts <- unique(sapply(rownames(res_df),function(x){strsplit(x,";")[[1]][1]}))

mappingofcts <- mappingofcts[mappingofcts$V2 %in% cts,]
mappingofcts$V2 <- factor(mappingofcts$V2,levels = mappingofcts$V2[order(mappingofcts$V3,decreasing = T)])

p <- ggplot(mappingofcts, aes(x = V2, y = V3, fill= V2)) +
  geom_bar(stat="identity",show.legend = F) +
  scale_fill_viridis(option="plasma", name="Process", discrete=T) +
  my_theme() +
  ylab("Number of Cell/Tissue Types") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("Fig2b.png",p,dpi = 300, width = 10, height = 7)

###### Fig. 2c #######
inp <- as.data.frame(table((sapply(rownames(res_df),function(x){strsplit(x,";")[[1]][1]}))))
inp$Label <- paste0(as.character(inp$Var1),"\n(",inp$Freq,")")
inp$Label[inp$Freq < 4000] <- ""

p <- ggplot(inp, aes(area = Freq, fill = Var1, label = Label)) +
  geom_treemap() +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    grow = T) +
  scale_fill_viridis(option="plasma", name="Cell Type", discrete=T) +
  my_theme() +
  theme(legend.position = "none")

ggsave("Fig2c.png",p,dpi = 300, width = 9, height = 7)


###### Fig. 2d #######
inp <- as.data.frame(table(res_df$PertType))
inp$Var1 <- as.character(inp$Var1)
inp$Var1[inp$Var1 == "ctl_untrt"] <- "Untreated"
inp$Var1[inp$Var1 == "ctl_vector"] <- "Mock\nVector"
inp$Var1[inp$Var1 == "ctl_vehicle"] <- "DMSO"
inp$Var1[inp$Var1 == "trt_cp"] <- "Compound"
inp$Var1[inp$Var1 == "trt_lig"] <- "Ligand"
inp$Var1[inp$Var1 == "trt_mut"] <- "Mutation"
inp$Var1[inp$Var1 == "trt_oe"] <- "Gene Overexpression"
inp$Var1[inp$Var1 == "trt_sh"] <- "shRNA"

inp$Var1 <- factor(as.character(inp$Var1),levels = as.character(inp$Var1[order(inp$Freq)]))

p <- ggplot(inp, aes(x = Var1, y = log10(Freq), fill = Var1)) +
  geom_bar(stat="identity") + 
  coord_polar(start = 0) +
  annotate("text", x = rep("Mutation",5), y = c(1, 2, 3, 4,5), label = c("1", "2", "3", "4","5") , color="black", size=5 , angle=0, fontface="bold", hjust=1) +
  scale_fill_viridis(option="plasma", name="Perturbation Type", discrete=T) +
  theme(legend.position = "none") +
  my_theme() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.minor.y = element_line(color = "black"),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill="white"),
    panel.border = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  )

ggsave("Fig2d.png",p,dpi = 300, width = 8, height = 7)

###### Fig. 3? #######
library(stringr)
res_df <- readRDS("~/Documents/RejuvenationPerturbationPaper/PredictedAgeDifferences.rds")
res_df <- res_df[which(!is.na(res_df$AgeDiff)),]
res_df <- res_df[which(grepl(";[0-9]+\ h$",rownames(res_df))),]
res_df$Timepoint <- str_extract(rownames(res_df),pattern = "[0-9]+\ h$")

perts_notimepoint <- gsub(";[0-9]+\ h$",replacement = "",x = rownames(res_df))
perts_notimepoint_stats <- as.data.frame(table(perts_notimepoint))
perts_notimepoint_stats <- perts_notimepoint_stats[perts_notimepoint_stats$Freq > 1,]

vals_list <- lapply(as.character(perts_notimepoint_stats$perts_notimepoint),function(x){
  if(min(res_df$zScore[which(grepl(x,rownames(res_df),fixed = T))]) > -2){
    return(NULL)
  }
  df <- data.frame(Perturbation = x,
             AgeDiff = res_df$AgeDiff[which(grepl(x,rownames(res_df),fixed = T))],
             zScore = res_df$zScore[which(grepl(x,rownames(res_df),fixed = T))],
             Timepoint = res_df$Timepoint[which(grepl(x,rownames(res_df),fixed = T))],
             stringsAsFactors = F)
  df$Timepoint <- as.numeric(gsub("\ h","",df$Timepoint))
  df$Timepoint <- rank(df$Timepoint)
  return(df)
})

vals_df <- do.call("rbind.data.frame",vals_list)
vals_df$Timepoint <- as.factor(vals_df$Timepoint)

p <- ggplot(vals_df, aes(y = Perturbation, x = Timepoint, fill= AgeDiff)) + 
  geom_tile() +
  scale_fill_viridis(option="plasma", name="Age Difference versus Control",discrete = F) +
  my_theme() +
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Timepoint") +
  ylab("")

ggsave("Fig3e.png",p,dpi = 300, width = 15, height = 5)

###### Fig. 3a #######
inp <- data.frame(Condition = rownames(res_df),
                  ZScore = res_df$zScore,
                  stringsAsFactors = F)
inp <- inp[!is.na(inp$ZScore),]
inp$Condition <- factor(inp$Condition,levels = inp$Condition[order(inp$ZScore)])
inp$x <- rank(inp$ZScore)

p <- ggplot(inp, aes(x=x, y=ZScore, color = ZScore)) +
  geom_point() +
  scale_color_viridis(option="plasma", name="Z-score",discrete = F) +
  #scale_y_cut(breaks=c(-3.5,7),scales = c(1,8,1),expand = expansion(mult = c(0.05,0.05))) +
  scale_y_cut(breaks=c(-4,7),scales = 0.25,space = 0.2) +
  my_theme() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  coord_cartesian(clip = "off")

ggsave("Fig3a.png",p,dpi = 300, width = 10, height = 7)

##### Fig 3b #####
inp <- res_df
inp_younger <- inp[which(inp$zScore <= -2),]
inp_older <- inp[which(inp$zScore >= 2),]

inp_younger <- as.data.frame(table(inp_younger$PertType))
inp_younger$Var1 <- as.character(inp_younger$Var1)
inp_younger$Var1[inp_younger$Var1 == "trt_cp"] <- "Compound"
inp_younger$Var1[inp_younger$Var1 == "trt_sh"] <- "shRNA"
inp_younger$Condition <- "Rejuvenating"

inp_older <- as.data.frame(table(inp_older$PertType))
inp_older$Var1 <- as.character(inp_older$Var1)
inp_older$Var1[inp_older$Var1 == "trt_cp"] <- "Compound"
inp_older$Var1[inp_older$Var1 == "trt_sh"] <- "shRNA"
inp_older$Var1[inp_older$Var1 == "trt_oe"] <- "Gene Overexpression"
inp_older$Condition <- "Pro-Aging"

inp_unchanged <- inp[which(inp$zScore > -2 & inp$zScore < 2),]
inp_unchanged <- as.data.frame(table(inp_unchanged$PertType))
inp_unchanged$Var1 <- as.character(inp_unchanged$Var1)
inp_unchanged$Var1[inp_unchanged$Var1 == "trt_cp"] <- "Compound"
inp_unchanged$Var1[inp_unchanged$Var1 == "trt_sh"] <- "shRNA"
inp_unchanged$Var1[inp_unchanged$Var1 == "trt_oe"] <- "Gene Overexpression"
inp_unchanged$Var1[inp_unchanged$Var1 == "trt_lig"] <- "Ligand"
inp_unchanged$Condition <- "No Effect"

inp_p <- rbind(inp_older,inp_younger,inp_unchanged)

p <- ggplot(inp_p, aes(fill=Var1, y=Freq, x=Condition)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(option="plasma", name="Perturbation Type", begin = 0.15, end = 0.85,discrete = T) +
  scale_y_cut(breaks=c(1000,12500),scales = 0.25) +
  my_theme() +
  ylab("Frequency")
  
ggsave("Fig3b.png",p,dpi = 300, width = 8, height = 7)


##### Fig 3c #####
inp <- read.table("~/Documents/RejuvenationPerturbationPaper/Fig3c_input.txt",sep = "\t",header =T, stringsAsFactors = F)
inp <- as.data.frame(table(inp$DrugAge.Prediction,inp$REVIVE.Prediction))
inp$Label <- paste0(round(inp$Freq/sum(inp$Freq)*100,2),"%")

mat <- matrix(inp$Freq,nrow = 2,ncol = 2)
chisq.test(t(mat)) #Significance test for association

p <- ggplot(inp, aes(x = Var1, y = Var2, fill= Freq)) + 
  geom_tile() +
  geom_text(aes(label=Label),size=6,color = "white") +
  scale_fill_viridis(option="plasma", name="",discrete = F,begin = 0.15, end = 0.85) +
  my_theme() +
  xlab("Significance (DrugAge)") +
  ylab("Significance (REVIVE)")

ggsave("Fig3c.png",p,dpi = 300, width = 7, height = 7)

##### Fig 3d #####
res_df <- readRDS("~/Documents/RejuvenationPerturbationPaper/PredictedAgeDifferences.rds")

part1 <- readRDS("../Hallmarks/1_10000.rds")
part2 <- readRDS("../Hallmarks/10001_20000.rds")
part3 <- readRDS("../Hallmarks/20001_30000.rds")
part4 <- readRDS("../Hallmarks/30001_70021.rds")

load("../Hallmarks/ProcessNames.rds",verbose =T)

hallmark_list <- c(part1,part2,part3,part4)

hallmark_list_younger <- hallmark_list[rownames(res_df)[which(res_df$zScore <= -2)]]
hallmark_younger_mat <- do.call("rbind",hallmark_list_younger)
colnames(hallmark_younger_mat) <- allprocs

expsToConsider <- c()
for(proc in 1:ncol(hallmark_younger_mat)){
  exps <- order(hallmark_younger_mat[,proc],decreasing = F)[1:5]
  idx <- which(hallmark_younger_mat[exps,proc] < 0)
  expsToConsider <- c(expsToConsider,rownames(hallmark_younger_mat)[exps[idx]])
}

expsToConsider <- unique(expsToConsider)

hallmark_younger_mat <- hallmark_younger_mat[expsToConsider,]

library(pheatmap)
p_heat <- pheatmap(hallmark_younger_mat,
                   scale = "none",
                   labels_row = "")

#write.table(res_df[which(res_df$zScore <= -2),],file = "../SupplementaryTableS6_Raw.txt", row.names = T, col.names = T, sep = "\t", quote = F)

inp <- melt(hallmark_younger_mat)
inp$Var1 <- factor(inp$Var1,levels = p_heat$tree_row$labels[p_heat$tree_row$order])
inp$Var2 <- factor(inp$Var2,levels = rev(p_heat$tree_col$labels[p_heat$tree_col$order]))

library(RColorBrewer)
color_x <- RColorBrewer::brewer.pal(12,"Paired")
names(color_x) <- rev(p_heat$tree_col$labels[p_heat$tree_col$order])

color_y <- rep("black",length(p_heat$tree_row$labels[p_heat$tree_row$order]))
names(color_y) <- p_heat$tree_row$labels[p_heat$tree_row$order]

label_y <- rep("",length(p_heat$tree_row$labels[p_heat$tree_row$order]))
names(label_y) <- p_heat$tree_row$labels[p_heat$tree_row$order]

for(intervention in 1:nrow(hallmark_younger_mat)){
  procs <- names(which.min(hallmark_younger_mat[intervention,]))
  
  color_y[rownames(hallmark_younger_mat)[intervention]] <- unname(color_x[procs])
  label_y[rownames(hallmark_younger_mat)[intervention]] <- rownames(hallmark_younger_mat)[intervention]
  
}

p <- ggplot(inp, aes(y = Var1, x = Var2, fill= value)) + 
  geom_tile() +
  scale_fill_viridis(option="plasma", name="Process Enrichment",discrete = F) +
  scale_y_discrete(labels = label_y) + 
  my_theme() +
  theme(axis.text.y = element_text(colour = color_y),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,colour = color_x)) +
  xlab("") +
  ylab("")

ggsave("Fig3d.png",p,dpi = 300, width = 15, height = 15)


######### Fig 4a #########
files <- list.files("~/Documents/RejuvenationPerturbationPaper/DrugRepurposing/AgingDEGs", pattern = ".rds", full.names = T)

inp_list <- list()
for(f in files){
  n <- gsub("/Users/sjung/Documents/RejuvenationPerturbationPaper/DrugRepurposing/AgingDEGs/","",f,fixed = T)
  n <- gsub(".rds","",n)
  spl <- strsplit(n,"_")[[1]]
  
  df <- readRDS(f)
  num_up <- length(which(df$avg_log2FC > 0))
  num_down <- length(which(df$avg_log2FC < 0))
  
  inp_list[[n]] <- data.frame(Tissue = c(spl[1],spl[1]),
                              CT = c(spl[2],spl[2]),
                              NumDiff = c(num_up,num_down),
                              Direction = c("UP","DOWN"),
                              stringsAsFactors = F)
}

names(inp_list) <- NULL

inp_df <- do.call("rbind",inp_list)

p <- ggplot(inp_df, aes(y = CT, x = NumDiff, fill= Direction)) +
geom_col() +
scale_fill_viridis(option="plasma", name="",discrete = T,begin = 0.15, end = 0.85) +
facet_grid(.~Tissue) +
my_theme() +
theme(panel.grid.major.y = element_line(color = "black", linewidth = 0.5)) +
#theme(legend.position = "none") +
xlab("") +
ylab("")

ggsave("Fig4a.png",p,dpi = 300, width = 15, height = 10)


#### Fig 4b ####
repurposingData <- readRDS("../DrugRepurposing/GSEA_Results_Combined.rds")

scoringFun <- function(x){
  idx_DOWN <- which(grepl("_DOWN",x$pathway))
  idx_UP <- which(grepl("_UP",x$pathway))
  
  return(x$NES[idx_DOWN]-x$NES[idx_UP])
}

conditions  <- repurposingData[[1]]$pathway
conditions <- gsub("_DOWN","",conditions)
conditions <- gsub("_UP","",conditions)
conditions <- unique(conditions)

allScores_list <- lapply(conditions,function(cond){
  df <- lapply(repurposingData,function(x){
    x <- as.data.frame(x)
    return(x[which(x[,"pathway"] %in% c(paste0(cond,"_DOWN"),paste0(cond,"_UP"))),c(1,2,3,6)])
  })
  scores <- sapply(df,scoringFun)
  names(scores) <- gsub("^[0-9]+_[0-9]+\\.","",names(scores))
  scores_rejuv <- scores[rownames(res_df)[which(res_df$zScore <= -2)]]
  return(scores_rejuv)
  
})

names(allScores_list) <- conditions

allScores_df <- do.call("cbind",allScores_list)

write.table(allScores_df,file = "../SupplementaryTableS7_Raw.txt",col.names = T, row.names = T, sep = "\t", quote =F)

allScores_mat <- as.matrix(allScores_df)

library(pheatmap)
p_heat <- pheatmap(allScores_mat,
                   scale = "row",
                   labels_col = "")

inp <- melt(allScores_mat)
inp$Var1 <- factor(inp$Var1,levels = p_heat$tree_row$labels[p_heat$tree_row$order])
inp$Var2 <- factor(inp$Var2,levels = rev(p_heat$tree_col$labels[p_heat$tree_col$order]))

p <- ggplot(inp, aes(y = Var1, x = Var2, fill= value)) + 
  geom_tile() +
  scale_fill_viridis(option="plasma", name="Perturbation Enrichment",discrete = F) +
  my_theme() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("") +
  ylab("")



### Fibroblast
fibroblast_input <- lapply(repurposingData,function(x){
  x <- as.data.frame(x)
  return(x[which(x[,"pathway"] %in% c("Skin_Fibroblast cell_DOWN","Skin_Fibroblast cell_UP")),c(1,2,3,6)])
})

fibroblast_scores <- sapply(fibroblast_input,scoringFun)
names(fibroblast_scores) <- gsub("^[0-9]+_[0-9]+\\.","",names(fibroblast_scores))
View(as.data.frame(sort(fibroblast_scores,decreasing = T)))

scores_rejuv <- fibroblast_scores[rownames(res_df)[which(res_df$zScore <= -2)]]
head(sort(scores_rejuv,decreasing = T))
View(as.data.frame(sort(scores_rejuv,decreasing = T)))

### Excitatory Neuron
exNeuron_input <- lapply(repurposingData,function(x){
  x <- as.data.frame(x)
  return(x[which(x[,"pathway"] %in% c("Brain_Excitatory neurons cell_DOWN","Brain_Excitatory neurons cell_UP")),c(1,2,3,6)])
})

exNeuron_scores <- sapply(exNeuron_input,scoringFun)
names(exNeuron_scores) <- gsub("^[0-9]+_[0-9]+\\.","",names(exNeuron_scores))
View(as.data.frame(sort(exNeuron_scores,decreasing = T)))

scores_rejuv <- exNeuron_scores[rownames(res_df)[which(res_df$zScore <= -2)]]
head(sort(scores_rejuv,decreasing = T))
View(as.data.frame(sort(scores_rejuv,decreasing = T)))



