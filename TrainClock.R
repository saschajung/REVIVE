mat <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/AdjustedExpMat_K30.rds")
load("/vols/GBIAntonio_bigdata/rejuvenationpredata/SampleGroups.RData",verbose = T)
group_lincs <- group_lincs[which(names(group_lincs) %in% colnames(mat))]
group_rna <- group_rna[which(names(group_rna) %in% colnames(mat))]
group_mic <- group_mic[which(names(group_mic) %in% colnames(mat))]
group_train <- group_train[which(names(group_train) %in% colnames(mat))]

train <- mat[,names(group_train)]
saveRDS(train,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/TrainingData.rds")
mic <- mat[,names(group_mic)]
saveRDS(mic,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/MicData.rds")
rna <- mat[,names(group_rna)]
saveRDS(rna,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/RNAData.rds")
lincs <- mat[,names(group_lincs)]

lincs1 <- lincs[,1:57171]
lincs2 <- lincs[,57172:114343]
lincs3 <- lincs[,114344:171515]
lincs4 <- lincs[,171516:228681]
saveRDS(lincs1,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs1Data.rds")
saveRDS(lincs2,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs2Data.rds")
saveRDS(lincs3,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs3Data.rds")
saveRDS(lincs4,"/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/Lincs4Data.rds")

load("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/mbotc.rda")
load("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/ArchS4_raw_manCur.RData",verbose = T)
trainExpr <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/TrainingData.rds")
trainAge <- props_train_manCur$Value[match(colnames(trainExpr),props_train_manCur$GSM)]
all_procs <- setdiff(unique(mbotc$ProcessName),"Background genes")
all_genes <- unique(mbotc$Symbol)

mapping <- readRDS("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/mbotc_symbToEns.rds")

genes_to_procs <- lapply(all_procs,function(x){unique(mbotc$Symbol[mbotc$ProcessName == x])})
names(genes_to_procs) <- all_procs
genes_to_procs <- lapply(genes_to_procs,function(x){unique(mapping$ensembl_gene_id[which(mapping$hgnc_symbol %in% x)])})
genes_to_procs <- lapply(genes_to_procs,function(x){intersect(rownames(trainExpr),x)})
genes_to_procs <- genes_to_procs[unname(which(sapply(genes_to_procs,length) >= 3))]
level_to_procs <- sapply(names(genes_to_procs),function(x){unique(mbotc$ProcessLevel[mbotc$ProcessName == x])})
levels_to_consider <- c(1)
level_to_procs <- level_to_procs[level_to_procs %in% levels_to_consider]

processes_to_consider <- names(level_to_procs)

genesInLevel <- intersect(unname(do.call(c,genes_to_procs[processes_to_consider])),rownames(trainExpr))

trainSamples <- colnames(trainExpr)

trainData <- trainExpr[genesInLevel,]
trainData <- trainData[,trainSamples]

trainData <- t(trainData)
trainData <- cbind(trainAge,trainData)
colnames(trainData)[1] <- "Age"

library(h2o)

conn <- h2o.init(max_mem_size="20G")
trainData_h2o <- as.h2o(trainData)

lambda = NULL
lambda_search = TRUE

model <- h2o.glm(y = "Age",
                 training_frame = trainData_h2o,
                 nfolds = 10,
                 fold_assignment = "Random",
                 family = "AUTO",
                 link = "family_default",
                 lambda_search = lambda_search,
                 lambda = lambda,
                 standardize = T,
                 alpha = 0.5,
                 seed = 100,
                 max_active_predictors = ncol(trainData_h2o),
                 solver = "IRLSM"
)

h2o.saveModel(model,path = "/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/")
