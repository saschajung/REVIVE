library(DescTools)
library(sva)
library(psych)
library(h2o)
library(biomaRt)
library(parallel)
library(doParallel)

#Get dataframe for mapping ensembl Ids to gene symbols 
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 mart = hsmart)

#Function to test new RNA-seq data (expects ensembl Ids)
testNewData <- function(newDataPath,origData,h2o_model,groundTruthAgesPath, cores = 8, ensToSymbMapping){
  
  #Read in the test data
  testData <- readRDS(newDataPath)
  
  # Step 1: Filter mapping table
  ensToSymbMapping <- ensToSymbMapping[
    ensToSymbMapping$ensembl_gene_id %in% rownames(testData) & 
      ensToSymbMapping$hgnc_symbol != "", 
  ]
  
  # Step 2: Match order of Ensembl IDs in testData and mapping
  idx <- match(rownames(testData), ensToSymbMapping$ensembl_gene_id)
  symbol_mapping <- ensToSymbMapping$hgnc_symbol[idx]
  keep <- !is.na(symbol_mapping)
  
  # Step 3: Collapse rows by gene symbol
  testData_symb <- rowsum(testData[keep, ], group = symbol_mapping[keep])
  
  # Optional: Sorting rows by symbol name
  testData_symb <- testData_symb[rownames(origData), ]
  
  #Remove samples with less than 1 million counts. Less depths may work well, but keeps it similar to training RNA-seq data
  toKeep <- which(colSums(testData_symb)>=10e+6)
  testData_symb <- testData_symb[,toKeep]
  
  #Perform log-normalization
  testData_normalized <- testData_symb[commonGenes,]
  testData_normalized <- log(testData_normalized+1)
  
  #Create cluster with specified number of cores
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  #Perform frozen surrogate variable analysis for each sample 
  testData_normalized <- foreach(i = 1:ncol(testData_normalized),.combine = "cbind",.inorder = T, .packages = "sva",.export = c("origData","mod","out")) %dopar% {
    newout <- fsva(dbdat = origData,
                   mod = mod[,1:2],
                   sv = out,
                   newdat = testData_normalized[,i,drop = F],
                   method = "exact")
    
    newout$new
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  rm(testData_symb)
  
  #Create predictions
  testData_h2o <- as.h2o(t(testData_normalized))
  preds <- h2o.predict(h2o_model,testData_h2o)
  
  preds_df <- as.data.frame(preds)
  rownames(preds_df) <- colnames(testData_normalized)
  
  #Read ground truth chronological ages
  groundTruth_df <- read.table(groundTruthAgesPath,header = T, sep = "\t", stringsAsFactors = F)
  rownames(groundTruth_df) <- groundTruth_df$GSM
  
  #Perform correlation analysis
  #Only consider samples above the age of 20
  comb <- cbind(groundTruth_df[colnames(testData_normalized),"age"],horvath_inverse_transform(preds_df[colnames(testData_normalized),"predict"]))
  rownames(comb) <- colnames(testData_normalized)
  comb <- comb[which(comb[,1] >= 20),]
  
  print(cor(comb))
  #Return predictions and normalized data
  return(list("Predictions"=comb,
              "Data"=testData_normalized))
  
}

load("./data/WS_45svs.RData") #

epi <- testNewData(newDataPath = "./data/celltype_data/Epithelial_raw_counts.rds",
                   origData = origData,
                   h2o_model = h2o_model,
                   groundTruthAgesPath = "./data/Epithelial_ages.tsv",
                   cores = 16,
                   ensToSymbMapping = mapping
)

hep <- testNewData(newDataPath = "./data/HepatocyteSamples_raw_counts.rds",
                   origData = origData,
                   h2o_model = h2o_model,
                   groundTruthAgesPath = "./data/Hepatocytes_ages.tsv",
                   cores = 16,
                   ensToSymbMapping = mapping
)

neu <- testNewData(newDataPath = "./data/Neurons_raw_counts.rds",
                   origData = origData,
                   h2o_model = h2o_model,
                   groundTruthAgesPath = "./data/Neurons_ages.tsv",
                   cores = 16,
                   ensToSymbMapping = mapping
)

hsc <- testNewData(newDataPath = "./data/HSC_raw_counts.rds",
                   origData = origData,
                   h2o_model = h2o_model,
                   groundTruthAgesPath = "./data/HSC_ages.tsv",
                   cores = 16,
                   ensToSymbMapping = mapping
)
