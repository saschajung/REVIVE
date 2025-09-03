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
counter <- 1
for(f in fs){
  print(counter)
  if(file.exists(paste0(tools::file_path_sans_ext(f),"_normalized_45svs.rds"))){
    counter <- counter + 1
    next
  }
  testData <- readRDS(f)
  
  testData_symb <- testData@mat
  testData_normalized <- testData_symb[commonGenes,]
  #Do not perform log normalization of LINCS pseudoRNA-seq data, since it is already log normalized
  
  idx <- intersect(samps,colnames(testData_normalized))
  print(length(idx))
  rm(testData)
  if(length(idx) == 0){
    next
  }
  testData_normalized <- testData_normalized[,idx]
  
  testData_normalized <- foreach(i = 1:ncol(testData_normalized),.combine = "cbind",.inorder = T, .packages = "sva") %dopar% {
    newout <- fsva(dbdat = origData,
                   mod = mod[,1:2],
                   sv = out,
                   newdat = testData_normalized[,i,drop = F],
                   method = "exact")
    
    newout$new
  }
  
  saveRDS(testData_normalized,file = paste0(tools::file_path_sans_ext(f),"_normalized_45svs.rds"))
  counter <- counter + 1
}

stopCluster(cl)
