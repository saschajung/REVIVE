library(codelink)
library(crossmeta)
library(GEOquery)
setwd("~/Documents/Microarray/")

samples <- readRDS("~/Documents/SelectedMicroarraySamples.rds")
samples <- unique(gsub("\\.[0-9]+","",samples))

saveRDS(samples,"~/Documents/GSEs_SelectedMicroarraySamples.rds")


gpllist <- list()
for(s in samples){
  print(s)
  gpl <- GEOquery::getGEO(GEO = s,
                          destdir = "./tmp",
                          GSEMatrix = F,
                          AnnotGPL = F,
                          getGPL = T)
  gpllist[[s]] <- names(gpl@gpls)
}

gpls <- unique(do.call("c",gpllist))
saveRDS(gpls,"~/Documents/GPLs_SelectedMicroarraySamples.rds")


library(Biobase)
for(s in samples){
  print(s)
  crossmeta::get_raw(s)
  raw <- crossmeta::load_raw(s,gpl_dir = "./GPL")
  mat <- exprs(raw[[1]])
  saveRDS(mat,paste0("./",s,".rds"))
}


###
library(stringr)
samples <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/SelectedMicroarraySamples.rds")
genes <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Genes_Lincs_RNA_Microarray.rds")

files <- c(list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/Javier/rankordered_javier",pattern = "*.csv",full.names = T),
           list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/Sascha/rankordered_sascha",pattern = "*.csv",full.names = T),
           list.files(path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/chempert/rankordered_chempert",pattern = "*.csv",full.names = T))

inp_list <- list()
for(s in samples){
  print(s)
  f <- files[which(grepl(paste0(s,".csv"),files))]
  if(length(f) > 1){
    print("MORE THANT ONE FILE FOUND")
    fnames <- unname(sapply(f,function(x){
      return(str_extract(x,"GSE[0-9]+\\.[0-9]+"))
    }))
    if(!all(fnames == fnames[1])){
      print("PROBLEM!")
      break
    }else{
      print("OK. JUST DUPLICATED FILE!")
      f <- f[1]
    }
    
  }
  inp_list[[f]] <- read.csv(f,row.names = 1)
  if(any(is.na(inp_list[[f]]))){
    print(paste0("NA found... Skipping sample ",f))
    inp_list[[f]] <- NULL
    next
  }
  if(!all(genes %in% rownames(inp_list[[f]]))){
    print("NOT ALL GENES FOUND!!! Skipping...")
    inp_list[[f]] <- NULL
    next
  }
  inp_list[[f]] <- inp_list[[f]][genes,]
}

mat <- do.call("cbind",inp_list)
cnames <- str_extract(colnames(mat),"GSM[0-9]+")
mat <- mat[,which(!duplicated(cnames))]
colnames(mat) <- str_extract(colnames(mat),"GSM[0-9]+")
mat <- as.matrix(mat)

saveRDS(mat,"/vols/GBIAntonio_bigdata/rejuvenationpredata/microarray/AllSamples_Ensembl_mat.rds")
