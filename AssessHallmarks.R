library(fgsea)
library(stringr)
genesets_hallmarks <- readRDS("~/Documents/RejuvenationPerturbationPaper/DrugRepurposing/AgingAtlas_HallmarkGenes_Complete_Ensembl.rds")
allgenes <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Repurposing/AllGenes_Ensembl.rds")

files_agingdegs <- list.files("/vols/GBIAntonio_bigdata/rejuvenationpredata/Repurposing/AgingDEGs",pattern = ".rds", full.names = T)
agingDEG_list <- list()
for(f in files_agingdegs){
  n <- gsub("/vols/GBIAntonio_bigdata/rejuvenationpredata/Repurposing/AgingDEGs/","",f,fixed = T)
  n <- gsub(".rds","",n,fixed = T)
  dat <- readRDS(f)
  dat <- dat[dat$ensembl %in% allgenes,]
  
  agingDEG_list[[paste0(n,"_UP")]] <- dat$ensembl[dat$avg_log2FC > 0]
  agingDEG_list[[paste0(n,"_DOWN")]] <- dat$ensembl[dat$avg_log2FC < 0]
  
}

files <- list.files("/vols/GBIAntonio_bigdata/rejuvenationpredata/Log2FCs",pattern = "[0-9]+_[0-9]+.rds", full.names = T)

library(foreach)
library(doParallel)

cores=12
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)

parts_list <- list()
for(f in files){
  print(Sys.time())
  load(f,verbose = T)
  res <- res[which(!sapply(res,is.null))]
  n <- gsub("/vols/GBIAntonio_bigdata/rejuvenationpredata/Log2FCs/","",f,fixed = T)
  n <- gsub(".rds","",n,fixed = T)
  
  enr_list <- list()
  clusterExport(cl,varlist = c("res","agingDEG_list","allgenes"))
  enr_list <- parLapply(cl = cl,
                        X = 1:length(res),
                        fun = function(i){
                          stat <- res[[i]]
                          names(stat) <- allgenes
                          gseares <- fgsea::fgsea(pathways = agingDEG_list,
                                                  stats = stat,
                                                  nproc = 1)
                          return(gseares)
                        })
  
  names(enr_list) <- names(res)
  
  #for(i in 1:length(res)){
  #  stat <- res[[i]]
  #  names(stat) <- allgenes
  #  gseares <- fgsea(pathways = agingDEG_list,
  #                   stats = stat,
  #                   nproc = 1)
  #  enr_list[[names(res)[i]]] <- gseares[,c(1,2,3,6)]
  #}
  
  parts_list[[n]] <- enr_list
  
}
