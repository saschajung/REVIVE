files <- list.files(path = "/vols/GBIAntonio_bigdata/LINCS",pattern = "*[0-9]+_[0-9]+.rds",full.names = T)
mapping <- readRDS("/vols/GBIAntonio_bigdata/LINCS/Rownames_SymbolToEnsembl.rds")
genes <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Genes_Lincs_RNA_Microarray.rds")
for(f in files){
  print(f)
  df <- readRDS(f)
  mapping$hgnc_symbol <- toupper(mapping$hgnc_symbol)
  idx <- match(mapping$hgnc_symbol,df@rid)
  mat <- df@mat[idx,]
  rownames(mat) <- mapping$ensembl_gene_id
  mat <- mat[order(rownames(mat),decreasing = F),]
  mat <- mat[rownames(mat) %in% genes,]
  f_new <- gsub("\\.rds","_ensembl_mat_commonGenes.rds",f)
  saveRDS(mat,file = f_new)
}

### Get largest possible intersection of genes
library(stringr)
genes_lincs_rna <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Genes_Lincs_RNA.rds")
files <- c(list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/Javier/rankordered_javier",pattern = "*.csv",full.names = T),
           list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/Sascha/rankordered_sascha",pattern = "*.csv",full.names = T),
           list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/chempert/rankordered_chempert",pattern = "*.csv",full.names = T))

genes_micro_list <- list()
for(f in files){
  df <- read.csv(f,row.names = 1)
  n <- str_extract(f,"GSE[0-9]+\\.[0-9]+")
  genes_micro_list[[n]] <- intersect(genes_lincs_rna,rownames(df))
}

saveRDS(genes_micro_list,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Genes_Lincs_RNA_Microarray_PerGSE.rds")

files <- c(list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/Javier/rankordered_javier",pattern = "*.csv",full.names = T),
           list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/Sascha/rankordered_sascha",pattern = "*.csv",full.names = T),
           list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/chempert/rankordered_chempert",pattern = "*.csv",full.names = T))

numSamples_micro_list <- list()
for(f in files){
  df <- read.csv(f,row.names = 1)
  n <- str_extract(f,"GSE[0-9]+\\.[0-9]+")
  numSamples_micro_list[[n]] <- ncol(df)
}

saveRDS(numSamples_micro_list,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Samples_Microarray_PerGSE.rds")


genes_micro_list <- readRDS("~/Documents/Genes_Lincs_RNA_Microarray_PerGSE.rds")
genes_micro_list <- genes_micro_list[which(sapply(genes_micro_list,length) >= 10000)]

visited <- rep(F,length(genes_micro_list))
names(visited) <- names(genes_micro_list)
equivClasses_list <- list()
for(i in 1:length(genes_micro_list)){
  print(i)
  if(visited[i]){
    next
  }
  samples <- names(which(sapply(genes_micro_list[(i+1):length(genes_micro_list)],function(x){setequal(x,genes_micro_list[[i]])})))
  samples <- c(samples,names(genes_micro_list)[i])
  equivClasses_list <- c(equivClasses_list,list(samples))
  visited[samples] <- T
}

equivClasses_genes_list <- list()
for(i in 1:length(equivClasses_list)){
  n <- paste0(equivClasses_list[[i]],collapse = "_")
  equivClasses_genes_list[[n]] <- genes_micro_list[[equivClasses_list[[i]][1]]]
}

numExp_PerEquivClass <- sapply(equivClasses_list,length)
numGenes_PerEquivClass <- sapply(equivClasses_genes_list,length)

library(rmoo)

max(sapply(equivClasses_genes_list,length))
sum(numExp_PerEquivClass)

fitfun <- function(x){
  #x <- strsplit(x,"")[[1]]
  idx <- which(x == 1)
  obj1 <- length(Reduce(intersect,equivClasses_genes_list[idx]))/16913
  obj2 <- sum(numExp_PerEquivClass[idx])/919
  return(c(-obj1,-obj2,abs(obj1-obj2)))
}

res <- rmoo(type = "binary",
            fitness = fitfun,
            nBits = 39,
            popSize = 1000,
            maxiter = 1000,
            strategy = "NSGA-II",
            nObj = 3,
            monitor = F,
            summary = F,
            seed = 1)

fit <- getFitness(res)
fit[,1] <- -fit[,1]*16913
fit[,2] <- -fit[,2]*919
fit$Idx <- 1:nrow(fit)
fit <- fit[fit$Fit_1 >= 10000 & fit$Fit_2 >= 400,]

library(plotly)
plot_ly(x = fit$Fit_1,y = fit$Fit_2)

#Take index 8
pop <- res@population[8,]

selectedEquivClasses <- equivClasses_genes_list[which(pop == 1)]

selectedSamples <- Reduce(union,strsplit(names(selectedEquivClasses),"_"))
selectedGenes <- Reduce(intersect,selectedEquivClasses)

saveRDS(res,file = "~/Documents/Microarray_SampleOptimization.rds")
saveRDS(selectedSamples, file = "~/Documents/SelectedMicroarraySamples.rds")
saveRDS(selectedGenes, file = "~/Documents/Genes_Lincs_RNA_Microarray.rds")
