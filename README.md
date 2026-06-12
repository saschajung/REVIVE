# REVIVE

This repository contains all scripts used to generate the results in the manuscript: "REVIVE: A Computational Platform for Systematically Identifying Rejuvenating Chemical and Genetic Perturbations"

## Description

- PerformGSEALINCS.R: Performs fast GSEA of the age-related log foldchanges with respect to the hallmark gene sets and aging DEGs.
- ChronologicalAgeValidation.R: Run the transcriptional clock on the chronological age validation datasets and compute pearson correlations.
- ComputeLogFCsLINCS.R: Compute log-fold changes for each conditions (treated vs control) of the sva normalized LINCS pseudoRNA-seq samples.
- GetCommonGenes.R: Obtain a list of genes common to all datasets.
- NormalizeLINCSData.R: Create sva corrected and normalized LINCS samples.
- PredictLINCSData.R: Apply the transcriptional clock to all LINCS pseudoRNA-seq samples after fsva correction.
- TrainModel.R: Train the transcriptional aging clock on the normalized and sva corrected training data.

## Applying the trained clock to new RNA-seq data

To score a new RNA-seq sample you need **both** artifacts produced by `TrainModel.R`
(neither is sufficient alone):

- `data/FinalModel` — the trained h2o GLM. This is the only object that holds the
  fitted gene coefficients (the clock itself).
- `data/WS_45svs.RData` — the frozen-SVA workspace. Loading it brings in `origData`,
  `mod`, `out`, `horvath_inverse_transform`, and `commonGenes`.

The procedure is: collapse the new data to the training gene space, log-normalise,
run frozen SVA per sample (`fsva`, `method = "exact"`), `h2o.predict`, then invert
the Horvath age transform. The full runnable script is
[`ApplyClockToNewData.R`](ApplyClockToNewData.R); the core is:

```r
library(sva)
library(h2o)
library(biomaRt)
library(parallel)
library(doParallel)

h2o.init(max_mem_size = "16G", nthreads = 8) #change according to your hardware
model <- h2o.loadModel("data/FinalModel")
load("data/WS_45svs.RData")     #origData, mod, out, horvath_inverse_transform, commonGenes

#new data: raw counts, Ensembl gene IDs in rows, samples in columns
newCounts <- readRDS("your_original_raw_counts.rds")

#map Ensembl IDs -> HGNC symbols (the clock uses HGNC symbol annotation)
mart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
map  <- getBM(c("ensembl_gene_id", "hgnc_symbol"), mart = mart)
map  <- map[map$ensembl_gene_id %in% rownames(newCounts) & map$hgnc_symbol != "", ]
sym  <- map$hgnc_symbol[match(rownames(newCounts), map$ensembl_gene_id)]
keep <- !is.na(sym)
expr <- rowsum(newCounts[keep, ], group = sym[keep])   #one row per gene symbol
expr <- expr[rownames(origData), ]                     #align to training gene set

#depth filter + log-normalise (RAW RNA-seq only; LINCS data is already log-scale)
expr <- expr[, colSums(expr) >= 10e6]                  #keep samples with >= 10M counts
expr <- log(expr[commonGenes, ] + 1)

#frozen SVA, one sample at a time (method = "exact").
#NB: sv = out (the whole sva object, NOT out$sv); mod[,1:2] = Intercept + Age.
cl <- makeCluster(8); registerDoParallel(cl)
adj <- foreach(i = seq_len(ncol(expr)), .combine = "cbind", .inorder = TRUE,
               .packages = "sva", .export = c("origData", "mod", "out")) %dopar% {
  fsva(dbdat = origData, mod = mod[, 1:2], sv = out,
       newdat = expr[, i, drop = FALSE], method = "exact")$new
}
stopCluster(cl)

#predict, back-transform to years
pred <- as.data.frame(h2o.predict(model, as.h2o(t(adj))))$predict
predicted_age <- horvath_inverse_transform(pred)
names(predicted_age) <- colnames(adj)
predicted_age
```
