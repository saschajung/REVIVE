mapping <- read.csv("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/tx2gene.human.GRCh38.rna.cdnaplusnrna.csv")
exsamp <- read.table("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/processed_final/sascha_illumina/GSM3191786/kallisto_output/SRR7346994/abundance.tsv",header = T,stringsAsFactors=F)

m <- merge(mapping,exsamp, by.x = "TXNAME", by.y = "target_id")
m <- m[,c("GENEID","eff_length")]
m <- aggregate(m$eff_length,by = list(Ensembl = m$GENEID),FUN = median)
colnames(m) <- c("Ensembl","Length")
m$Ensembl <- gsub("\\.[0-9]+","",m$Ensembl)

saveRDS(m,"/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/EnsemblToLength.rds")

df_j <- read.csv("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/javier_111_counts.csv",row.names = 1)
df_s <- read.csv("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/chempert_111_counts.csv",row.names = 1)
df_c <- read.csv("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/111_sascha_complete_counts_df.csv",row.names = 1)

rownames(df_j) <- gsub("\\.[0-9]+","",rownames(df_j))
rownames(df_s) <- gsub("\\.[0-9]+","",rownames(df_s))
rownames(df_c) <- gsub("\\.[0-9]+","",rownames(df_c))

df_j <- as.matrix(df_j)
df_s <- as.matrix(df_s)
df_c <- as.matrix(df_c)

m <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/EnsemblToLength.rds")

library(DGEobj.utils)
df_j <- convertCounts(countsMatrix = df_j,
                      unit = "TPM",
                      geneLength = m$Length,
                      normalize = "none")

df_s <- convertCounts(countsMatrix = df_s,
                      unit = "TPM",
                      geneLength = m$Length,
                      normalize = "none")

df_c <- convertCounts(countsMatrix = df_c,
                      unit = "TPM",
                      geneLength = m$Length,
                      normalize = "none")

genes <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Genes_Lincs_RNA_Microarray.rds")

all(genes %in% rownames(df_j))
all(genes %in% rownames(df_s))
all(genes %in% rownames(df_c))

df_j <- df_j[genes,]
df_s <- df_s[genes,]
df_c <- df_c[genes,]

df <- cbind(df_j,df_s,df_c)

saveRDS(df,file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/bulk/AllSamples_RNA_Ensembl_mat.rds")
