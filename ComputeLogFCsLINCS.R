library(cmapR)
library(sva)
library(parallel)
library(doParallel)
library(tools)

fs <- list.files(path = "./data/",pattern = "*normalized_45svs.rds", full.names = T)
data_list <- list()
k <- 1
for(f in fs){
  data_list[[k]] <- readRDS(f)
  k <- k + 1
}
data <- do.call("cbind",data_list)
rm(data_list)

load("./data/SampleGroups.RData",verbose = T)

group_all <- gsub(";GFP;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";LUCIFERASE;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";EGFP;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";BFP;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";HCRED;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";LACZ;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";EMPTY_VECTOR;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";PGW;",";CONTROL_VECTOR;",group_all)
group_all <- gsub(";RFP;",";CONTROL_VECTOR;",group_all)

all(colnames(data) %in% names(group_all))
group_all <- group_all[colnames(data)]

out_list <- list()
k <- 1
#Load the prediction table of the LINCS data
predsummary <- read.table("./PredictedAgeDifferences_sign_revision.txt",row.names = 1, header = T, sep = "\t", stringsAsFactors = F, quote = "") 
for(n in rownames(predsummary)){
  print(k)
  g_split <- strsplit(n,";")[[1]]
  if(g_split[3] == "trt_cp"){
    g_control <- paste0(g_split[1],";DMSO;ctl_vehicle;\"\";",g_split[5])
  }else if(g_split[3] == "trt_oe"){
    if(grepl("^CMAP\\-ERGK[0-9]{3}$",g_split[2])){
      g_control <- paste0(g_split[1],";CMAP\\-ERGK[0-9]{3};ctl_vector;\"\";",g_split[5])
    }else if(grepl("^CMAP\\-ERGK[0-9]{3}\\.DOX$",g_split[2])){
      g_control <- paste0(g_split[1],";CMAP\\-ERGK[0-9]{3}\\.DOX;ctl_vector;\"\";",g_split[5])
    }else{
      g_control <- paste0(g_split[1],";CONTROL_VECTOR;ctl_vector;",g_split[4],";",g_split[5])
    }
  }else if(g_split[3] == "trt_sh"){
    g_control <- paste0(g_split[1],";UnTrt;ctl_untrt;\"\";","0 h")
  }else if(g_split[3] == "trt_lig"){
    g_control <- paste0(g_split[1],";UnTrt;ctl_untrt;\"\";","0 h")
  }else if(g_split[3] == "trt_xpr"){
    g_control <- paste0(g_split[1],";UnTrt;ctl_untrt;\"\";","0 h")
  }else if(g_split[3] %in% c("ctl_vector","ctl_vehicle","ctl_untrt")){
    next
  }else{
    next #Should not happen!
  }
  
  exp_samps <- names(which(group_all == n))
  con_samps <- names(group_all)[which(grepl(g_control,group_all,fixed = T))]
  
  out_list[[n]] <- rowMeans(data[,exp_samps]) - rowMeans(data[,con_samps]) #Since the data is log-normalized, just calculate the difference
  k <- k + 1
  
}

rm(data)

logfcs_df <- do.call("rbind",out_list)

saveRDS(logfcs_df,file = "./data/logfcs_normalized_45svs.rds")
