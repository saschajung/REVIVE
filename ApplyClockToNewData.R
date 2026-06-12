# Apply the trained REVIVE transcriptional aging clock to new RNA-seq data.
#
# Requires BOTH artifacts produced by TrainModel.R (neither is sufficient alone):
#   data/FinalModel     - trained h2o GLM (the only object with fitted gene
#                         coefficients; this is the clock)
#   data/WS_45svs.RData - frozen-SVA workspace; loads origData, mod, out,
#                         horvath_inverse_transform, commonGenes
#
# Input: raw counts, Ensembl gene IDs in rows, samples in columns.
# Output: a named numeric vector of predicted ages (years), one per sample.

library(sva)
library(h2o)
library(biomaRt)
library(parallel)
library(doParallel)

# ---- configure ----
newDataPath <- "your_new_raw_counts.rds"   # raw counts, Ensembl IDs x samples
modelPath   <- "data/FinalModel"
svaWSPath   <- "data/WS_45svs.RData"
cores       <- 8

# ---- load model + frozen-SVA workspace ----
h2o.init(max_mem_size = "16G", nthreads = cores)
model <- h2o.loadModel(modelPath)
load(svaWSPath)   # origData, mod, out, horvath_inverse_transform, commonGenes

# ---- read new data and map Ensembl IDs -> HGNC symbols ----
newCounts <- readRDS(newDataPath)

mart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
map  <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = mart)
map  <- map[map$ensembl_gene_id %in% rownames(newCounts) & map$hgnc_symbol != "", ]
sym  <- map$hgnc_symbol[match(rownames(newCounts), map$ensembl_gene_id)]
keep <- !is.na(sym)

expr <- rowsum(newCounts[keep, ], group = sym[keep])   # collapse to one row per symbol
expr <- expr[rownames(origData), ]                     # align to the training gene set/order

# ---- depth filter + log-normalise ----
# RAW RNA-seq only. LINCS pseudo-RNA-seq is already log-scale: skip the log() step.
expr <- expr[, colSums(expr) >= 10e6]                  # keep samples with >= 10M counts
expr <- log(expr[commonGenes, ] + 1)

# ---- frozen SVA: project each new sample onto the training SVA solution ----
# One sample at a time, method = "exact" (by design).
# NB: sv = out  (the WHOLE sva object, NOT out$sv).  mod[,1:2] = Intercept + Age.
cl <- makeCluster(cores)
registerDoParallel(cl)
adj <- foreach(i = seq_len(ncol(expr)), .combine = "cbind", .inorder = TRUE,
               .packages = "sva", .export = c("origData", "mod", "out")) %dopar% {
  fsva(dbdat = origData, mod = mod[, 1:2], sv = out,
       newdat = expr[, i, drop = FALSE], method = "exact")$new
}
stopCluster(cl)

# ---- predict + back-transform to years ----
pred <- as.data.frame(h2o.predict(model, as.h2o(t(adj))))$predict
predicted_age <- horvath_inverse_transform(pred)
names(predicted_age) <- colnames(adj)

print(predicted_age)
