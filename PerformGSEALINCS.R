library(sva)
library(foreach)
library(parallel)
library(doParallel)
library(tools)
library(fgsea)
library(GSA)

logfcs_df <- readRDS("./data/logfcs_normalized_45svs.rds")

#Preprocess hallmark genes from the aging atlas
test <- GSA::GSA.read.gmt("./data/AgingAtlas_HallmarkGenes_Complete.gmt")
out_list <- test$genesets
names(out_list) <- test$geneset.names
out_list <- lapply(out_list,function(x){
  setdiff(unique(x),"")
})
saveRDS(out_list,file = "./data/AgingAtlas_HallmarkGenes_Complete.rds")

genesets_hallmarks <- readRDS("./data/AgingAtlas_HallmarkGenes_Complete.rds")
#Load prediction table
predsummary <- read.table("./data/PredictedAgeDifferences_sign_revision.txt",row.names = 1, header = T, sep = "\t", stringsAsFactors = F, quote = "")
parts_list <- list()
for(n in rownames(predsummary)){
  print(which(rownames(predsummary) == n))
  stat <- logfcs_df[n,]
  gseares <- fgsea::fgsea(pathways = genesets_hallmarks,
                          stats = stat,
                          nPermSimple = 1000,
                          nproc = 1)
  
  parts_list[[n]] <- gseares
  
}

saveRDS(parts_list,file = "./data/GSEA_logfc_PredictedAgeDifferences_sign_revision_Hallmarks.rds")

##### GSEA of aging DEGs #####
#Preprocess aging DEGs from Aging Atlas
files_agingdegs <- list.files("./AgingDEGs",pattern = ".rds", full.names = F)
agingDEG_list <- list()
for(f in files_agingdegs){
  n <- gsub(".rds","",n,fixed = T)
  dat <- readRDS(f)
  dat <- dat[dat$gene %in% colnames(logfcs_df),]
  
  agingDEG_list[[paste0(n,"_UP")]] <- dat$gene[dat$avg_log2FC > 0]
  agingDEG_list[[paste0(n,"_DOWN")]] <- dat$gene[dat$avg_log2FC < 0]
  
}

#Get significant conditions from the predictions
predsummary <- read.table("./PredictedAgeDifferences_sign_revision.txt",row.names = 1, header = T, sep = "\t", stringsAsFactors = F, quote = "")
predsummary <- predsummary[which(predsummary$AgeDiff < 0 & predsummary$t.qval < 0.01),]

logfcs_df <- logfcs_df[rownames(predsummary),]
cores=30
cl <- makeCluster(cores)
registerDoParallel(cl)
clusterExport(cl = cl,c("predsummary","agingDEG_list","logfcs_df"))

parts_list <- foreach(i = 1:nrow(predsummary),.combine = "c",.inorder = T, .packages = "fgsea") %dopar% {
  stat <- logfcs_df[rownames(predsummary)[i],]
  gseares <- fgsea::fgsea(pathways = agingDEG_list,
                          stats = stat,
                          nPermSimple = 1000,
                          nproc = 1)
  list(gseares)
}

names(parts_list) <- rownames(predsummary)
saveRDS(parts_list,file = "./data/GSEA_logfc_PredictedAgeDifferences_sign_rejuvenating_revision_AgingDEGs.rds")

stopCluster(cl)

