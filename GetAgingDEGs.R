library(biomaRt)

ageDEGs_df <- read.table("~/Documents/RejuvenationPerturbationPaper/DrugRepurposing/AgingDEGs/Aging-related DEGs.txt", header = T, sep = ",", stringsAsFactors = F)
ageDEGs_df$cell_type[which(ageDEGs_df$cell_type == "\xa6æ\xc4 T cell")] <- "Gamma delta T cell"
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

for(ti in unique(ageDEGs_df$Tissue)){
  for(ct in unique(ageDEGs_df$cell_type[which(ageDEGs_df$Tissue == ti)])){
    print(paste0(ti,"_",ct))
    df <- ageDEGs_df[which(ageDEGs_df$Tissue == ti & ageDEGs_df$cell_type == ct & ageDEGs_df$group == "old vs youth"),]
    df <- df[,c("gene","avg_log2FC","p_val","p_val_adj")]
    
    if(nrow(df) == 0){
      next
    }
    
    mapping <- getBM(
      attributes = c('ensembl_gene_id',"hgnc_symbol"), 
      filters = 'hgnc_symbol',
      values = unique(df$gene),
      mart = hsmart
    )
    
    df <- merge(df,mapping,by.x = "gene",by.y = "hgnc_symbol")
    df <- df[,c(1,5,2,3,4)]
    
    colnames(df)[2] <- "ensembl"
    
    saveRDS(df,file = paste0("~/Documents/RejuvenationPerturbationPaper/DrugRepurposing/AgingDEGs/",ti,"_",ct,".rds"))
  }
}
