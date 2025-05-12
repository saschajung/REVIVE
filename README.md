# REVIVE

This repository contains all scripts used to generate the results in the manuscript: "REVIVE: A Computational Platform for Systematically Identifying Rejuvenating Chemical and Genetic Perturbations"

## Description
- AssessHallmarks.R: Performs fast GSEA of the age-related log2 foldchanges with respect to the hallmark gene sets.
- CountsToTPM.R: Converts all RNAseq of the intervention data to TPMs.
- DataSelectionAndStandardization.R: Select LINCS samples and transform to ENSEMBL ids.
- GetAgingDEGs.R: Obtain a list of age-related DEGs from AgeAnno data.
- MakeClockPredictions.R: Apply the transcriptional clock to all intervention data.
- NormalizeAllData.R: Normalize all datasets (i.e. intervention data and training samples) together using a custom implementation of the Rank-in method.
- Plotting: Generate the plots in the manuscript.
- ProcessMicroarray.R: Process and obtain the microarray data.
- ProcessOriginalTrainingData.R: Process the original training data from the MultiTIMER paper and convert to TPMs.
- SelectMaximumSamples.R: Find a pareto-optimum that maximizes the number of samples and the number of common genes.
- StandardizeMetadata.R: Map cell types and interventions to a standard vocabulary to ease comparison.
- TrainClock.R: Train the transcriptional aging clock on the normalized training data.
