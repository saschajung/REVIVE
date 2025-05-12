library(biomaRt)
library(DGEobj.utils)
load("/vols/GBIAntonio_bigdata/ArchS4_raw_manCur.RData",verbose = T)

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org", biomart = "ensembl")
mapping1 <- getBM(
  attributes = c('ensembl_gene_id',"hgnc_symbol"), 
  filters = 'hgnc_symbol',
  values = rownames(expr_train_manCur),
  mart = hsmart
)

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mapping2 <- getBM(
  attributes = c('ensembl_gene_id',"hgnc_symbol"), 
  filters = 'hgnc_symbol',
  values = rownames(expr_train_manCur),
  mart = hsmart
)

mapping <- rbind(mapping1,mapping2)
mapping <- mapping[!duplicated(mapping),]

expr_train_manCur <- expr_train_manCur[which(rownames(expr_train_manCur) %in% mapping$hgnc_symbol),]
mapping <- mapping[mapping$hgnc_symbol %in% rownames(expr_train_manCur),]
mapped_list <- lapply(rownames(expr_train_manCur),function(x){
  idx <- which(mapping$hgnc_symbol == x)
  df <- expr_train_manCur[rep(x,length(idx)),,drop = F]
  rownames(df) <- mapping$ensembl_gene_id[idx]
  return(df)
})

mat <- do.call("rbind",mapped_list)
mat <- mat[genes,]
any(is.na(mat))
dim(mat)

saveRDS(mat, file = "~/Documents/TrainingSamples_Ensembl_mat.rds")

mat <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/TrainingSamples_Ensembl_mat.rds")

ensToLength <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/EnsemblToLength.rds")
ensToLength <- ensToLength[ensToLength$Ensembl %in% rownames(mat),]
mat <- convertCounts(countsMatrix = mat,
                     unit = "TPM",
                     geneLength = ensToLength$Length,
                     normalize = "none")

saveRDS(mat, file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/TrainingSamples_Ensembl_mat.rds")
