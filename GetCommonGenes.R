library(cmapR)
library(biomaRt)
setwd("./data/") #Assuming that all datasets are in the 'data' subfolder 

#Select common genes between LINCS, training data and test data
lincs <- readRDS("cp_predicted_RNAseq_profiles_1_1000.rds") #Read one of pseudoRNA-seq subsets
load("ArchS4_raw_manCur.RData",verbose = T) #Load training data
testData <- readRDS("HSC_raw_counts.rds") #Load one of the test datasets

#The test data has ensembl gene IDs. Transform to HGNC symbols
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mapping1 <- getBM(
  attributes = c('ensembl_gene_id',"hgnc_symbol"), 
  filters = 'ensembl_gene_id',
  values = rownames(testData),
  mart = hsmart
)
mapping1 <- mapping1[which(mapping1$hgnc_symbol != ""),]
testData_symb <- matrix(0,nrow = length(unique(mapping1$hgnc_symbol)),ncol = ncol(testData))
rownames(testData_symb) <- sort(unique(mapping1$hgnc_symbol))
colnames(testData_symb) <- colnames(testData)
for(i in 1:nrow(testData_symb)){
  ens <- mapping1$ensembl_gene_id[which(mapping1$hgnc_symbol == rownames(testData_symb)[i])]
  testData_symb[i,] <- unname(colSums(testData[ens,,drop = F]))
}

#Find common genes
commonGenes <- intersect(lincs@rid,rownames(expr_train_manCur)) #'expr_train_manCur' is the training data object
commonGenes <- intersect(rownames(testData_symb),commonGenes)

#Remove all the data except the common genes
rm(lincs)
rm(expr_train_manCur)
rm(props_train_manCur)
rm(testData_symb)
rm(testData)
rm(mapping1)
saveRDS(commonGenes,file = "./commonGenes.rds")
