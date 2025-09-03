library(DescTools)
library(sva)
library(psych)
library(h2o)

setwd("./")

load("./data/commonGenes.rds")
#Some age transformation functions
horvath_transform <- function(age) {
  ifelse(age <= 20,
         log(age + 1) - log(21),
         (age - 20) / 21)
}

horvath_inverse_transform <- function(transformed_age) {
  ifelse(transformed_age <= 0,
         exp(transformed_age + log(21)) - 1,
         21 * transformed_age + 20)
}

###Curate the training data
ctmap <- read.table("~/TrainSamplesToCellTypes.txt",header = F, sep = "\t",stringsAsFactors = F)
load("/vols/GBIAntonio_bigdata/rejuvenationpredata/Clock/ArchS4_raw_manCur.RData",verbose = T)
trainExpr <- expr_train_manCur
trainExpr <- trainExpr[commonGenes,]
normgenes <- intersect(commonGenes,agestable)

#Remove cancer samples
toKeep <- ctmap$V1[which(ctmap$V2 != "Cancer")]
ctmap <- ctmap[which(ctmap$V2 != "Cancer"),]
props_train_manCur <- props_train_manCur[which(props_train_manCur$GSM %in% toKeep),]
trainExpr <- trainExpr[,which(colnames(trainExpr) %in% toKeep)]
trainAge <- props_train_manCur$Value[match(colnames(trainExpr),props_train_manCur$GSM)]
names(trainAge) <- colnames(trainExpr)

#Remove samples below 20 years
toKeep <- which(trainAge >= 20)
trainExpr <- trainExpr[,toKeep]
trainAge <- trainAge[toKeep]
trainAge <- round(trainAge)

#Remove samples with less than 10 million counts
toKeep <- which(colSums(trainExpr)>=10e+6)
trainExpr <- trainExpr[,toKeep]
trainAge <- trainAge[toKeep]

#Compute weights for cell types
freq_df <- as.data.frame(table(ctmap[match(colnames(trainExpr),ctmap$V1),"V2"]))
freq_df$Freq <- ceiling(freq_df$Freq / 10) * 10
nom <- 15000 #In principle, the least common multiple would be the best, but that becomes too large. Set 15000 instead for approximative weights.
freq_df$Weight <- round(nom/freq_df$Freq)

weights_ct <- freq_df$Weight[match(ctmap[match(colnames(trainExpr),ctmap$V1),"V2"],freq_df$Var1)]
names(weights_ct) <- colnames(trainExpr)

#Compute weights for ages
freq_age_df <- as.data.frame(table(trainAge))
freq_age_df$Freq <- ceiling(freq_age_df$Freq / 10) * 10
nom <- 15000
freq_age_df$Weight <- round(nom/freq_age_df$Freq)

weights_age <- freq_age_df$Weight[match(trainAge,freq_age_df$trainAge)]
names(weights_age) <- names(trainAge)

all(names(weights_age) == names(weights_ct)) #Needs to be true

weights <- weights_ct + weights_age #Take the sum as the final training weights

#Log-normalize the training data
trainExpr <- log(trainExpr+1)

#Create metadata frame
age_df <- props_train_manCur[match(colnames(trainExpr),props_train_manCur$GSM),]
age_df$TransValue <- horvath_transform(age_df$Value) #Perform the age-transformation of Horvath. It's basically linear for the training samples since ages >= 20 

trainSamples <- colnames(trainExpr)
trainAge <- age_df$TransValue[match(trainSamples,age_df$GSM)]

trainData <- trainExpr[,trainSamples]

#Perform surrogate variable analysis
meta <- data.frame(Group = as.factor(rep("1",ncol(trainData))),
                   Age = trainAge,
                   row.names = colnames(trainData))
mod <- model.matrix(~Age,meta) #Keep age as the covariate of interest
nsv_be <- num.sv(trainData, mod, method = "be")
out <- sva(dat = trainData,
           mod = mod,
           n.sv = nsv_be)

sapply(1:ncol(out$sv),function(i){cor(trainAge,out$sv[,i])}) #Look at the correlation of each surrogate variable with age. It should be low if sva worked.

#Adjust the training data. Same adjustment that is used in frozen surrogate variable analysis (fsva function)
nmod = dim(mod)[2]
mod1 = cbind(mod,out$sv)
gammahat = (trainData %*% mod1 %*% solve(t(mod1) %*% mod1))[,(nmod+1):(nmod + out$n.sv)]
trainData = trainData - gammahat %*% t(out$sv)

#Remove correlated genes to avoid multicollinearity. Threshold is 0.7 pearson correlation
#Compute correlation of each gene with age
cors <- sapply(rownames(trainData),function(x){cor(trainData[x,],trainAge)})
cors <- sort(abs(cors),decreasing = T)
head(cors)

#Efficiently compute the weighted correlation matrix with the weights computed before
corMat <- cor.wt(t(trainData),vars=NULL, w=unname(weights[colnames(trainData)]),sds=NULL, cor=TRUE)
corMat <- corMat$r
corMat[is.na(corMat)] <- 0 #Should not happen

#Perform a greedy selection of genes by looping through all genes based on their absolute correlation with age
selectedGenes <- c()
for(g in names(cors)){
  print(g)
  subgenes <- colnames(corMat)[which(corMat[g,] >= 0.7)] #Get all correlated genes
  if(length(selectedGenes) > 0){
    subgenes <- subgenes[which(rowSums(corMat[subgenes,selectedGenes,drop=F] >= 0.7) == 0)] #Select only genes that are not correlated with any of the already selected genes
  }
  if(length(subgenes) == 0){
    print("............skipped")
    next
  }
  agecors <- sapply(subgenes,function(x){
    tmp <- cor.wt(cbind(trainData[x,],trainAge),w=unname(weights[colnames(trainData)]))$r
    abs(tmp[1,2])
  })
  selectedGenes <- c(selectedGenes,names(agecors)[which(agecors == max(agecors))]) #Take the gene that is most informative about age
}

length(selectedGenes) #Just see how many genes have been finally selected

#Compose the final training data frame including age and sample weights
trainData <- trainData[selectedGenes,]

trainData <- t(trainData)
trainData <- cbind(trainAge,trainData)
colnames(trainData)[1] <- "Age"

w <- weights[rownames(trainData)]
trainData <- cbind(unname(w),trainData)
colnames(trainData)[1] <- "Weight"

#Only retain those genes that are differentially expressed between young and old (60-69 or 70-79) in at least one tissue
attr_mat <- read.table("./data/DEGenes_GTEx.txt",sep = "\t",header = T)
rownames(attr_mat) <- attr_mat$Gene
attr_mat$Gene <- NULL

attr_mat <- attr_mat[,which(grepl("20\\.29\\.vs\\.[67]0\\.[67]9",colnames(attr_mat)))]
attr_mat <- attr_mat[names(which(!apply(attr_mat,1,function(x){all(x == 0)}))),]

attr_mat <- attr_mat[which(apply(attr_mat,1,function(x){
  pos <- sum(x == 1)
  neg <- sum(x == -1)
  return(pos == 0 | neg == 0)
})),]

x_genes <- intersect(selectedGenes,rownames(attr_mat))

#Train the transcriptional aging clock
conn <- h2o.init(max_mem_size="16G",nthreads = 8)
trainData_h2o <- as.h2o(trainData)

lambda = NULL
lambda_search = TRUE

#Run first model inference using an elastic net
model <- h2o.glm(y = "Age",
                 x = x_genes,
                 training_frame = trainData_h2o,
                 nfolds = 10,
                 fold_assignment = "Modulo",
                 early_stopping = F,
                 family = "AUTO",
                 link = "family_default",
                 lambda_search = lambda_search,
                 lambda = lambda,
                 standardize = T,
                 weights_column = "Weight",
                 alpha = 0.25,
                 seed = 100,
                 solver = "IRLSM"
)

#Refine model with strict lasso regression
activePredictors <- setdiff(names(which(model@model$coefficients != 0)),"Intercept")
model1 <- h2o.glm(y = "Age",
                  x = activePredictors,
                  training_frame = trainData_h2o,
                  nfolds = 10,
                  fold_assignment = "Modulo",
                  early_stopping = F,
                  family = "AUTO",
                  link = "family_default",
                  lambda_search = lambda_search,
                  #nlambdas = 200,
                  lambda = lambda,
                  standardize = T,
                  weights_column = "Weight",
                  alpha = 1,
                  seed = 100,
                  max_active_predictors = nrow(trainData_h2o),
                  solver = "IRLSM"
)

#Save model and data needed for applying the model to new data
h2o.saveModel(model1, path = "./data/", force = T,filename = "FinalModel")
save(origData,mod,out,horvath_inverse_transform,commonGenes,file = "./data/WS_45svs.RData")






