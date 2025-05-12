#Create adjusted rank matrices
setwd("/vols/GBIAntonio_bigdata/LINCS")

bin_equal = function(x, nbin = 5) {
  breaks = quantile(x, probs = seq(0, 1, length.out = nbin + 1), na.rm = TRUE)
  return(findInterval(x, breaks, all.inside = TRUE))
}

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/lincs/AllSamples_Ensembl_mat.rds")
lincs <- lincs[,1:114340]
gc()
rankmat_lincs <- apply(lincs,2,function(x){
  r <- base::rank(x,ties.method = "min")
  return(ceiling((r*100)/max(r)))
})
#rankmat_lincs <- rankmat_lincs[1:nrow(lincs),]
rownames(rankmat_lincs) <- rownames(lincs)
rankmat_lincs[1:5,1:5]
rm(lincs)

saveRDS(rankmat_lincs,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part1.rds")

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/lincs/AllSamples_Ensembl_mat.rds")
lincs <- lincs[,114341:ncol(lincs)]
gc()
rankmat_lincs <- apply(lincs,2,function(x){
  r <- base::rank(x,ties.method = "min")
  return(ceiling((r*100)/max(r)))
})
#rankmat_lincs <- rankmat_lincs[1:nrow(lincs),]
rownames(rankmat_lincs) <- rownames(lincs)
rankmat_lincs[1:5,1:5]
rm(lincs)

saveRDS(rankmat_lincs,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part2.rds")


rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/AllSamples_RNA_Ensembl_mat.rds")
rna <- log2(rna+1)
rankmat_rna <- apply(rna,2,function(x){
  r <- base::rank(x,ties.method = "min")
  return(ceiling((r*100)/max(r)))
})
#rankmat_rna <- rankmat_rna[1:nrow(rna),]
rownames(rankmat_rna) <- rownames(rna)
rankmat_rna[1:5,1:5]
rm(rna)

saveRDS(rankmat_rna,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/RNA_Ranks_mat.rds")


micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/AllSamples_Ensembl_mat.rds")
rankmat_micro <- apply(micro,2,function(x){
  r <- base::rank(x,ties.method = "min")
  return(ceiling((r*100)/max(r)))
})
#rankmat_micro <- rankmat_micro[1:nrow(micro),]
rownames(rankmat_micro) <- rownames(micro)
rankmat_micro[1:5,1:5]
rm(micro)
any(is.na(micro))
saveRDS(rankmat_micro,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Micro_Ranks_mat.rds")


train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/TrainingSamples_Ensembl_mat.rds")
train <- log2(train+1)
rankmat_train <- apply(train,2,function(x){
  r <- base::rank(x,ties.method = "min")
  return(ceiling((r*100)/max(r)))
})
#rankmat_train <- rankmat_train[1:nrow(train),]
rownames(rankmat_train) <- rownames(train)
rankmat_train[1:5,1:5]
rm(train)

saveRDS(rankmat_train,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Training_Ranks_mat.rds")


lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/lincs/AllSamples_Ensembl_mat.rds")
rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/AllSamples_RNA_Ensembl_mat.rds")
micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/AllSamples_Ensembl_mat.rds")
train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/TrainingSamples_Ensembl_mat.rds")
comb <- cbind(lincs,rna,micro,train)
rm(lincs)
rm(rna)
rm(micro)
rm(train)
gc()

expvec <- as.vector(comb)
rm(comb)
saveRDS(expvec,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_ExpVec.rds")
rm(expvec)

lincs1 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part1.rds")
lincs2 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part2.rds")
rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/RNA_Ranks_mat.rds")
micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Micro_Ranks_mat.rds")
train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Training_Ranks_mat.rds")
comb <- cbind(lincs1,lincs2,rna,micro,train)
rm(lincs1)
rm(lincs2)
rm(rna)
rm(micro)
rm(train)
gc()
rankvec <- as.vector(comb)
rm(comb)
gc()
rankvec2 <- rankvec^2

saveRDS(rankvec,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_RankVec.rds")
saveRDS(rankvec2,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_RankVec2.rds")

rm(list = ls())

set.seed(42)
e <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_ExpVec.rds")
idx <- floor(runif(28157893)*length(e))
idx <- sapply(idx,function(x){max(1,x)})
e <- e[idx]
gc()
df <- data.frame("e" = e,
                 stringsAsFactors = F)
rm(e)
gc()
r2 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_RankVec2.rds")
r2 <- r2[idx]
gc()
df$r2 <- r2
rm(r2)
gc()
r <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_RankVec.rds")
r <- r[idx]
gc()
df$r <- r
rm(r)
gc()

mod <- lm(e ~ r2 + r, data = df)
summary(mod)

lincs1 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part1.rds")
lincs2 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part2.rds")
rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/RNA_Ranks_mat.rds")
micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Micro_Ranks_mat.rds")
train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Training_Ranks_mat.rds")
comb <- cbind(lincs1,lincs2,rna,micro,train)
rm(lincs1)
rm(lincs2)
rm(rna)
rm(micro)
rm(train)
gc()
weightmat = 2*unname(mod$coefficients[2])*comb + unname(mod$coefficients[3])
saveRDS(weightmat,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Weights_mat.rds")

rankmat_adj <- comb * weightmat
rankmat_adj[1:5,1:5]

saveRDS(rankmat_adj,"/vols/GBIAntonio_bigdata/rejuvenationpredata/RanksAdj_mat.rds")


#Split into individual rankmat_adj per type
lincs1 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part1.rds")
lincs2 <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_Ranks_mat_Part2.rds")
all(c(colnames(lincs1),colnames(lincs2)) %in% colnames(rankmat_adj))
rankmat_adj_l <- rankmat_adj[,c(colnames(lincs1),colnames(lincs2))]
saveRDS(rankmat_adj_l,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_RanksAdj_mat.rds")
rm(lincs1)
rm(lincs2)
rm(rankmat_adj_l)
gc()

rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/RNA_Ranks_mat.rds")
all(colnames(rna) %in% colnames(rankmat_adj))
rankmat_adj_l <- rankmat_adj[,colnames(rna)]
saveRDS(rankmat_adj_l,"/vols/GBIAntonio_bigdata/rejuvenationpredata/RNA_RanksAdj_mat.rds")
rm(rna)
rm(rankmat_adj_l)
gc()

micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Micro_Ranks_mat.rds")
all(colnames(micro) %in% colnames(rankmat_adj))
rankmat_adj_l <- rankmat_adj[,colnames(micro)]
saveRDS(rankmat_adj_l,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Micro_RanksAdj_mat.rds")
rm(micro)
rm(rankmat_adj_l)
gc()

train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Training_Ranks_mat.rds")
all(colnames(train) %in% colnames(rankmat_adj))
rankmat_adj_l <- rankmat_adj[,colnames(train)]
saveRDS(rankmat_adj_l,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Training_RanksAdj_mat.rds")
rm(train)
rm(rankmat_adj_l)
gc()

#Create M matrix
load("/vols/GBIAntonio_bigdata/rejuvenationpredata/SampleGroups.RData")

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/lincs/AllSamples_Ensembl_mat.rds")
rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/AllSamples_RNA_Ensembl_mat.rds")
micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/AllSamples_Ensembl_mat.rds")
micro <- micro[,which(!colnames(micro) %in% c(colnames(rna),colnames(lincs)))]
train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/TrainingSamples_Ensembl_mat.rds")
train <- train[,which(!colnames(train) %in% c(colnames(rna),colnames(lincs),colnames(micro)))]
comb <- cbind(lincs,rna,micro,train)
rm(lincs)
rm(rna)
rm(micro)
rm(train)

group_all <- group_all[intersect(names(group_all),colnames(comb))]

groupsUnique <- unique(unname(group_all))
for(i in 1:length(groupsUnique)){
  group <- groupsUnique[i]
  print(group)
  samps <- names(which(group_all == group))
  mat <- comb[,samps,drop = F]
  m <- rowMeans(mat)
  saveRDS(m,file = paste0("/vols/GBIAntonio_bigdata/rejuvenationpredata/M_Matrices/",i,".rds"))
  gc()
}

saveRDS(groupsUnique,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/M_Matrices/UniqueGroups.rds")
rm(comb)

lincs <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Lincs_RanksAdj_mat.rds")
rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/RNA_RanksAdj_mat.rds")
micro <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Micro_RanksAdj_mat.rds")
micro <- micro[,which(!colnames(micro) %in% c(colnames(rna),colnames(lincs)))]
train <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Training_RanksAdj_mat.rds")
train <- train[,which(!colnames(train) %in% c(colnames(rna),colnames(lincs),colnames(micro)))]
comb <- cbind(lincs,rna,micro,train)
rm(lincs)
rm(rna)
rm(micro)
rm(train)

groupsUnique <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/M_Matrices/UniqueGroups.rds")
for(i in 1:length(groupsUnique)){
  group <- groupsUnique[i]
  print(i)
  m <- readRDS(paste0("/vols/GBIAntonio_bigdata/rejuvenationpredata/M_Matrices/",i,".rds"))
  samps <- names(which(group_all == group))
  mat <- comb[,samps,drop = F]
  mat <- sweep(mat,1,m)
  saveRDS(mat,file = paste0("/vols/GBIAntonio_bigdata/rejuvenationpredata/MeanAdjRanks/",i,".rds"))
}

#Combine adjusted matrix
cols <- colnames(comb)
rm(comb)

groupsUnique <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/M_Matrices/UniqueGroups.rds")
mat_list <- list()
for(i in 1:length(groupsUnique)){
  mat_list[[i]] <- readRDS(paste0("/vols/GBIAntonio_bigdata/rejuvenationpredata/MeanAdjRanks/",i,".rds"))
  print(dim(mat_list[[i]]))
}

adj_mat <- do.call("cbind",mat_list)
adj_mat[1:5,1:5]
rm(mat_list)
gc()
saveRDS(adj_mat,"/vols/GBIAntonio_bigdata/rejuvenationpredata/MeanAdjRanks_mat.rds")

#install.packages("bigstatsr")
library(bigstatsr)

adj_mat <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/MeanAdjRanks_mat.rds")
cnames <- colnames(adj_mat)
saveRDS(cnames, file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Colnames_MeanAdjRanks.rds")

comb_FBM <- as_FBM(x = adj_mat,
                   type = "double",
                   backingfile = "/vols/GBIAntonio_bigdata/rejuvenationpredata/MeanAdjRanks_FBM_new")
rm(adj_mat)
gc()

res <- big_randomSVD(X = comb_FBM,
                     k = 50,
                     ncores = 1)

saveRDS(res,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Results_RSVD_K50.rds")

#Determine number of singular values to consider
res <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Results_RSVD_K50.rds")
cnames <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Colnames_MeanAdjRanks.rds")
load("~/Documents/SampleGroups.RData",verbose = T)

type_lincs <- rep("LINCS",length(group_lincs))
names(type_lincs) <- names(group_lincs)
type_mic <- rep("MIC",length(group_mic))
names(type_mic) <- names(group_mic)
type_rna <- rep("RNA",length(group_rna))
names(type_rna) <- names(group_rna)
type_train <- rep("TRAIN",length(group_train))
names(type_train) <- names(group_train)

type_all <- c(type_lincs,type_mic,type_rna,type_train)

all(cnames %in% names(group_all))

type_all <- type_all[cnames]
group_all <- group_all[cnames]
any(is.na(group_all))

res <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Results_RSVD_K50.rds")
idx <- c(5,6,7,8,20,22,23,25,27,28,30,31,32,33,34,35,36,37,38,39,41,42,43,44,45,46,47,48,49,50)
res$d <- res$d[idx]
res$u <- res$u[,idx]
res$v <- res$v[,idx]

svdmat <- res$u %*% diag(res$d) %*% t(res$v)

rankmat <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/RanksAdj_mat.rds")
cnames <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Colnames_MeanAdjRanks.rds")
rankmat <- rankmat[,cnames]

#save(cnames,"/vols/GBIAntonio_bigdata/rejuvenationpredata/AllSamples_Finalmat.rds")

outmat <- rankmat - svdmat

saveRDS(outmat, "/vols/GBIAntonio_bigdata/rejuvenationpredata/AdjustedExpMat_K30.rds")
