library(h2o)
model <- h2o.loadModel("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/GLM_model_R_1708656074038_1")
h2o.rm(trainData_h2o)

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs1Data.rds")
numParts <- ceiling(ncol(lincs)/1000)
preds_list <- list()
for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  lincs_small <- t(lincs[,from:to])
  lincs_h2o <- as.h2o(lincs_small)
  pred_h2o <- h2o.predict(model,lincs_h2o)
  preds_list[[i]] <- as.data.frame(pred_h2o)
  h2o.rm(lincs_h2o)
}

for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  rownames(preds_list[[i]]) <- colnames(lincs)[from:to]
}

preds_df <- do.call("rbind",preds_list)

saveRDS(preds_df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Preds_Lincs1.rds")

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs2Data.rds")
numParts <- ceiling(ncol(lincs)/1000)
preds_list <- list()
for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  lincs_small <- t(lincs[,from:to])
  lincs_h2o <- as.h2o(lincs_small)
  pred_h2o <- h2o.predict(model,lincs_h2o)
  preds_list[[i]] <- as.data.frame(pred_h2o)
  h2o.rm(lincs_h2o)
}

for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  rownames(preds_list[[i]]) <- colnames(lincs)[from:to]
}

preds_df <- do.call("rbind",preds_list)

saveRDS(preds_df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Preds_Lincs2.rds")

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs3Data.rds")
numParts <- ceiling(ncol(lincs)/1000)
preds_list <- list()
for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  lincs_small <- t(lincs[,from:to])
  lincs_h2o <- as.h2o(lincs_small)
  pred_h2o <- h2o.predict(model,lincs_h2o)
  preds_list[[i]] <- as.data.frame(pred_h2o)
  h2o.rm(lincs_h2o)
}

for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  rownames(preds_list[[i]]) <- colnames(lincs)[from:to]
}

preds_df <- do.call("rbind",preds_list)

saveRDS(preds_df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Preds_Lincs3.rds")

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs4Data.rds")
numParts <- ceiling(ncol(lincs)/1000)
preds_list <- list()
for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  lincs_small <- t(lincs[,from:to])
  lincs_h2o <- as.h2o(lincs_small)
  pred_h2o <- h2o.predict(model,lincs_h2o)
  preds_list[[i]] <- as.data.frame(pred_h2o)
  h2o.rm(lincs_h2o)
}

for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(lincs))
  rownames(preds_list[[i]]) <- colnames(lincs)[from:to]
}

preds_df <- do.call("rbind",preds_list)

saveRDS(preds_df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Preds_Lincs4.rds")

rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/RNAData.rds")
numParts <- ceiling(ncol(rna)/1000)
preds_list <- list()
for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(rna))
  lincs_small <- t(rna[,from:to])
  lincs_h2o <- as.h2o(lincs_small)
  pred_h2o <- h2o.predict(model,lincs_h2o)
  preds_list[[i]] <- as.data.frame(pred_h2o)
  h2o.rm(lincs_h2o)
}

for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(rna))
  rownames(preds_list[[i]]) <- colnames(rna)[from:to]
}

preds_df <- do.call("rbind",preds_list)

saveRDS(preds_df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Preds_RNA.rds")

micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/MicData.rds")
numParts <- ceiling(ncol(micro)/1000)
preds_list <- list()
for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(micro))
  lincs_small <- t(micro[,from:to])
  lincs_h2o <- as.h2o(lincs_small)
  pred_h2o <- h2o.predict(model,lincs_h2o)
  preds_list[[i]] <- as.data.frame(pred_h2o)
  h2o.rm(lincs_h2o)
}

for(i in 1:numParts){
  print(i)
  from <- (i-1)*1000 + 1
  to <- min(i*1000,ncol(micro))
  rownames(preds_list[[i]]) <- colnames(micro)[from:to]
}

preds_df <- do.call("rbind",preds_list)

saveRDS(preds_df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Preds_Mic.rds")


unique(sapply(unname(group_all),function(x){
  strsplit(x,";")[[1]][2]
}))

instances <- read.table("~/Documents/LINCS/instinfo_beta.txt",sep = "\t", header = T, stringsAsFactors = F,comment.char = "", quote = "")
dim(instances)
cellinfo <- read.table("~/Documents/LINCS/cellinfo_beta.txt",sep = "\t", header = T, stringsAsFactors = F,comment.char = "",quote = "")
cellinfo <- cellinfo[cellinfo$cell_type == "normal",]

instances <- instances[instances$qc_pass == 1,]
instances <- instances[instances$cell_iname %in% cellinfo$cell_iname,]

instances$cell_iname[instances$cell_iname %in% c("HEK293T","HEK293","HEKTE","AALE")] <- "Kidney Cell"
instances$cell_iname[instances$cell_iname %in% c("HA1E","HPTEC","HME","LHSAR","NL20")] <- "Epithelial Cell"
instances$cell_iname[instances$cell_iname %in% c("NPC")] <- "Neural Progenitor Cell"
instances$cell_iname[instances$cell_iname %in% c("HUVEC")] <- "Endothelial Cell"
instances$cell_iname[instances$cell_iname %in% c("NAMEC8")] <- "Mesenchymal Cell"
instances$cell_iname[instances$cell_iname %in% c("PHH")] <- "Hepatocyte"
instances$cell_iname[instances$cell_iname %in% c("CD34")] <- "HSC"
instances$cell_iname[instances$cell_iname %in% c("WA09","HUES3")] <- "Embryonic Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("ASC")] <- "Adipose-derived Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("SKL")] <- "Bone Marrow Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("NEU")] <- "Neural Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("WI38","1HAE","HFL1","IMR90")] <- "Fibroblast"
instances$cell_iname[instances$cell_iname %in% c("CD34")] <- "HSC"

instances$pert_itime[which(instances$cmap_name == "UnTrt")] <- "0 h"

instances$cmap_name[which(instances$cmap_name %in% c("GFP","LUCIFERASE","EGFP","BFP","HCRED","LACZ","EMPTY_VECTOR","PGW","RFP"))] <- "CONTROL_VECTOR"

preds <- rbind(readRDS("~/Documents/Preds_Lincs1.rds"),
               readRDS("~/Documents/Preds_Lincs2.rds"),
               readRDS("~/Documents/Preds_Lincs3.rds"),
               readRDS("~/Documents/Preds_Lincs4.rds"),
               readRDS("~/Documents/Preds_RNA.rds"),
               readRDS("~/Documents/Preds_Mic.rds"))
load("~/Documents/SampleGroups.RData",verbose = T)

group_all <- gsub(";GFP;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";LUCIFERASE;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";EGFP;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";BFP;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";HCRED;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";LACZ;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";EMPTY_VECTOR;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";PGW;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";RFP;",";CONTROL_VECTOR;",group_all)

all(names(group_all) %in% rownames(preds))
all(rownames(preds) %in% names(group_all))

save(group_all,preds,file = "~/Downloads/AgePredData.RData")

groups <- group_all[rownames(preds)]

library(foreach)
library(doParallel)

cores=16
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#numSamps_control <- c()
#nummSamps_exp <- c()
#diffAge <- c()
ugroups <- unique(unname(groups))
clusterExport(cl,varlist = c("groups","preds"))
res <- parLapply(cl =cl,
                 X = ugroups,
                 fun = function(g){
                   g_split <- strsplit(g,";")[[1]]
                   if(g_split[3] == "trt_cp"){
                     g_control <- paste0(g_split[1],";DMSO;ctl_vehicle;\"\";",g_split[5])
                   }else if(g_split[3] == "trt_oe"){
                     if(grepl("^CMAP\\-ERGK[0-9]{3}$",g_split[2])){
                       g_control <- paste0(g_split[1],";CMAP\\-ERGK[0-9]{3};ctl_vector;\"\";",g_split[5])
                     }else if(grepl("^CMAP\\-ERGK[0-9]{3}\\.DOX$",g_split[2])){
                       g_control <- paste0(g_split[1],";CMAP\\-ERGK[0-9]{3}\\.DOX;ctl_vector;\"\";",g_split[5])
                     }else{
                       g_control <- paste0(g_split[1],";CONTROL_VECTOR;ctl_vector;",g_split[4],";",g_split[5])
                     }
                   }else if(g_split[3] == "trt_sh"){
                     g_control <- paste0(g_split[1],";UnTrt;ctl_untrt;\"\";","0 h")
                   }else if(g_split[3] == "trt_lig"){
                     g_control <- paste0(g_split[1],";UnTrt;ctl_untrt;\"\";","0 h")
                   }else if(g_split[3] == "trt_xpr"){
                     g_control <- paste0(g_split[1],";UnTrt;ctl_untrt;\"\";","0 h")
                   }else if(g_split[3] %in% c("ctl_vector","ctl_vehicle","ctl_untrt")){
                     return(c(-500,-500,-500))
                   }else{
                     #print(paste0("ERRROR:::",g))
                     return(c(-1000,-1000,-1000))
                     #break
                   }
                   
                   exp_samps <- names(which(groups == g))
                   con_samps <- names(groups)[which(grepl(g_control,groups,fixed = T))]
                   
                   return(c(length(con_samps),length(exp_samps),mean(preds[exp_samps,"predict"]) - mean(preds[con_samps,"predict"])))
                   
                 })

names(res) <- ugroups

perturbagen <- unname(sapply(ugroups,function(x){strsplit(x,";")[[1]][2]}))
perttype <- unname(sapply(ugroups,function(x){strsplit(x,";")[[1]][3]}))

res_df <- do.call("rbind",res)
colnames(res_df) <- c("NumConSamps","NumExpSamps","AgeDiff")
res_df[res_df == -500] <- NA
res_df$AgeDiff[is.nan(res_df$AgeDiff)] <- NA
res_df[res_df == -1000] <- NA
res_df$AgeDiff[res_df$Perturbagen %in% c("shRNA","siRNA")] <- NA
res_df <- as.data.frame(res_df)
res_df$Perturbagen <- perturbagen
res_df$PertType <- perttype

res_df$zScore <- (res_df$AgeDiff-mean(res_df$AgeDiff,na.rm = T))/sd(res_df$AgeDiff,na.rm = T)

saveRDS(res_df,file = "~/Documents/RejuvenationPerturbationPaper/PredictedAgeDifferences.rds")

