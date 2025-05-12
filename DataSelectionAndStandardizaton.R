library(reshape2)
setwd("~/Documents/LINCS/")

instances <- read.table("~/Documents/LINCS/instinfo_beta.txt",sep = "\t", header = T, stringsAsFactors = F,comment.char = "", quote = "")
dim(instances)
cellinfo <- read.table("~/Documents/LINCS/cellinfo_beta.txt",sep = "\t", header = T, stringsAsFactors = F,comment.char = "",quote = "")
cellinfo <- cellinfo[cellinfo$cell_type == "normal",]

instances <- instances[instances$qc_pass == 1,]
instances <- instances[instances$cell_iname %in% cellinfo$cell_iname,]
table(instances$pert_time)
table(instances$cell_iname,instances$pert_type)

trt_cps <- instances$sample_id[instances$pert_type == "trt_cp" & instances$cell_iname %in% c("1HAE","ASC","CD34","HA1E","HEK293","HEK293T","HFL1","HPTEC","HUES3","HUVEC","IMR90","NEU","NL20","NPC","PHH","SKL","WA09","WI38")]
trt_vehicles <- instances$sample_id[instances$pert_type == "ctl_vehicle" & instances$cell_iname %in% c("1HAE","ASC","CD34","HA1E","HEK293","HEK293T","HFL1","HPTEC","HUES3","HUVEC","IMR90","NEU","NL20","NPC","PHH","SKL","WA09","WI38")]

trt_oe <- instances$sample_id[instances$pert_type == "trt_oe" & instances$cell_iname %in% c("AALE","HA1E","HEK293","HEK293T","LHSAR")]
trt_vectors <- instances$sample_id[instances$pert_type == "ctl_vector" & instances$cell_iname %in% c("AALE","HA1E","HEK293","HEK293T","LHSAR")]

trt_lig <- instances$sample_id[instances$pert_type == "trt_lig" & instances$cell_iname %in% c("HA1E","HEK293T")]
trt_sh <- instances$sample_id[instances$pert_type == "trt_sh" & instances$cell_iname %in% c("ASC","HA1E","HEK293T","HEKTE","NPC","SKL")]
trt_untrt <- instances$sample_id[instances$pert_type == "ctl_untrt" & instances$cell_iname %in% c("ASC","HA1E","HEK293T","HEKTE","NPC","SKL")]

saveRDS(trt_cps,"selectedSamples_LINCS_cp.rds")
saveRDS(trt_vehicles,"selectedSamples_LINCS_vehicle.rds")
saveRDS(trt_oe,"selectedSamples_LINCS_oe.rds")
saveRDS(trt_vectors,"selectedSamples_LINCS_vector.rds")
saveRDS(trt_lig,"selectedSamples_LINCS_lig.rds")
saveRDS(trt_sh,"selectedSamples_LINCS_sh.rds")
saveRDS(trt_untrt,"selectedSamples_LINCS_untrt.rds")


saveRDS(instances,"selectedSamples_LINCS.rds")

cps <- instances$sample_id[instances$pert_type == "trt_cp"]
lig <- instances$sample_id[instances$pert_type == "trt_lig"]
oe <- instances$sample_id[instances$pert_type == "trt_oe"]
ctl <- instances$sample_id[instances$pert_type == "trt_ctl"]

cps <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_cp.rds")
numiters <- ceiling(length(cps)/1000)
for(i in 1:numiters){
  print(i)
  s <- (i-1)*1000 + 1
  e <- min(i*1000,length(cps))
  my_ds_1000_columns <- parse_gctx("/vols/GBIAntonio_bigdata/LINCS/cp_predicted_RNAseq_profiles.gctx", cid=cps[s:e])
  saveRDS(my_ds_1000_columns,paste0("/vols/GBIAntonio_bigdata/LINCS/cp_predicted_RNAseq_profiles_",s,"_",e,".rds"))
}

lig <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_lig.rds")
numiters <- ceiling(length(lig)/1000)
for(i in 1:numiters){
  print(i)
  s <- (i-1)*1000 + 1
  e <- min(i*1000,length(lig))
  my_ds_1000_columns <- parse_gctx("/vols/GBIAntonio_bigdata/LINCS/lig_predicted_RNAseq_profiles.gctx", cid=lig[s:e])
  saveRDS(my_ds_1000_columns,paste0("/vols/GBIAntonio_bigdata/LINCS/lig_predicted_RNAseq_profiles",s,"_",e,".rds"))
}

oe <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_oe.rds")
numiters <- ceiling(length(oe)/1000)
for(i in 1:numiters){
  print(i)
  s <- (i-1)*1000 + 1
  e <- min(i*1000,length(oe))
  my_ds_1000_columns <- parse_gctx("/vols/GBIAntonio_bigdata/LINCS/oe_predicted_RNAseq_profiles.gctx", cid=oe[s:e])
  saveRDS(my_ds_1000_columns,paste0("/vols/GBIAntonio_bigdata/LINCS/oe_predicted_RNAseq_profiles",s,"_",e,".rds"))
}

sh <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_sh.rds")
numiters <- ceiling(length(sh)/1000)
for(i in 1:numiters){
  print(i)
  s <- (i-1)*1000 + 1
  e <- min(i*1000,length(sh))
  my_ds_1000_columns <- parse_gctx("/vols/GBIAntonio_bigdata/LINCS/shRNA_predicted_RNAseq_profiles.gctx", cid=sh[s:e])
  saveRDS(my_ds_1000_columns,paste0("/vols/GBIAntonio_bigdata/LINCS/shRNA_predicted_RNAseq_profiles",s,"_",e,".rds"))
}

untr <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_untrt.rds")
vehicle <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_vehicle.rds")
vector <- readRDS("/vols/GBIAntonio_bigdata/LINCS/selectedSamples_LINCS_vector.rds")
ctl <- c(untr,vehicle,vector)
numiters <- ceiling(length(ctl)/1000)
for(i in 1:numiters){
  tryCatch({
    print(i)
    s <- (i-1)*1000 + 1
    e <- min(i*1000,length(ctl))
    my_ds_1000_columns <- parse_gctx("/vols/GBIAntonio_bigdata/LINCS/ctl_predicted_RNAseq_profiles.gctx", cid=ctl[s:e])
    saveRDS(my_ds_1000_columns,paste0("/vols/GBIAntonio_bigdata/LINCS/ctl_predicted_RNAseq_profiles",s,"_",e,".rds"))
  },warning = function(x){print(x)},
  error = function(x){print(x)
    print(paste0("Iteration ",i," failed."))})
  
}

#Map LINCS data to ensembl ids
library(biomaRt)
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
gs <- readRDS("~/Documents/Rownames_Symbol.rds")

mapping <- getBM(
  attributes = c("ensembl_gene_id","hgnc_symbol"), 
  filters = "hgnc_symbol",
  values = gs,
  mart = hsmart
)

saveRDS(mapping,"~/Documents/Rownames_SymbolToEnsembl.rds")
