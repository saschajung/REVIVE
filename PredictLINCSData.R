library(cmapR)
library(sva)
library(parallel)
library(doParallel)
library(tools)

commonGenes <- readRDS("./data/commonGenes.rds")

setwd("./data")
samps <- read.table("TableS7_raw.txt",sep = "\t",header =F,stringsAsFactors = F) #Load final samples to save time
samps <- samps[,1] 
fs <- c(list.files(pattern = "cp_predicted_RNAseq_profiles_[0-9]+_[0-9+e]+.rds"),
        list.files(pattern = "ctl_predicted_RNAseq_profiles[0-9]+_[0-9]+.rds"),
        list.files(pattern = "lig_predicted_RNAseq_profiles[0-9]+_[0-9]+.rds"),
        list.files(pattern = "oe_predicted_RNAseq_profiles[0-9]+_[0-9]+.rds"),
        list.files(pattern = "shRNA_predicted_RNAseq_profiles[0-9]+_[0-9]+.rds"))

load("./data/WS_45svs.RData",verbose = T)
cl <- makeCluster(30) #Change if more resources are available
registerDoParallel(cl)
clusterExport(cl = cl,varlist = c("origData","mod","out"))
for(f in fs){
  print(f)
  if(file.exists(paste0(tools::file_path_sans_ext(f),"_preds_45svs.rds"))){
    next
  }
  testData <- readRDS(f)
  
  testData_symb <- testData@mat
  testData_normalized <- testData_symb[commonGenes,]
  rm(testData)
  
  testData_normalized <- foreach(i = 1:ncol(testData_normalized),.combine = "cbind",.inorder = T, .packages = "sva",.export = c("origData","mod","out")) %dopar% {
    newout <- fsva(dbdat = origData,
                   mod = mod[,1:2],
                   sv = out,
                   newdat = testData_normalized[,i,drop = F],
                   method = "exact")
    
    newout$new
  }
  
  rm(testData_symb)
  
  testData_h2o <- as.h2o(t(testData_normalized))
  
  preds <- h2o.predict(h2o_model,testData_h2o)
  
  preds_df <- as.data.frame(preds)
  rownames(preds_df) <- colnames(testData_normalized)
  preds_df$predict <- horvath_inverse_transform(preds_df$predict)
  
  saveRDS(preds_df,file = paste0(tools::file_path_sans_ext(f),"_preds_45svs.rds"))
  h2o.rm(testData_h2o)
  rm(testData_normalized)
}

stopCluster(cl)
