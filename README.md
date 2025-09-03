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
