## LINCS GROUPS ##
instances <- read.table("~/Documents/LINCS/instinfo_beta.txt",sep = "\t", header = T, stringsAsFactors = F,comment.char = "", quote = "")
dim(instances)
cellinfo <- read.table("~/Documents/LINCS/cellinfo_beta.txt",sep = "\t", header = T, stringsAsFactors = F,comment.char = "",quote = "")
cellinfo <- cellinfo[cellinfo$cell_type == "normal",]

instances <- instances[instances$qc_pass == 1,]
instances <- instances[instances$cell_iname %in% cellinfo$cell_iname,]

instances$cell_iname[instances$cell_iname %in% c("HEK293T","HEK293","HEKTE","AALE")] <- "Kidney Cell"
instances$cell_iname[instances$cell_iname %in% c("HA1E","HPTEC","HME","LHSAR","NL20")] <- "Epithelial Cell"
instances$cell_iname[instances$cell_iname %in% c("NPC")] <- "Neural Progenitor Cell"
instances$cell_iname[instances$cell_iname %in% c("HUVEC")] <- "Endothelial Cell"
instances$cell_iname[instances$cell_iname %in% c("NAMEC8")] <- "Mesenchymal Cell"
instances$cell_iname[instances$cell_iname %in% c("PHH")] <- "Hepatocyte"
instances$cell_iname[instances$cell_iname %in% c("CD34")] <- "HSC"
instances$cell_iname[instances$cell_iname %in% c("WA09","HUES3")] <- "Embryonic Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("ASC")] <- "Adipose-derived Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("SKL")] <- "Bone Marrow Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("NEU")] <- "Neural Stem Cell"
instances$cell_iname[instances$cell_iname %in% c("WI38","1HAE","HFL1","IMR90")] <- "Fibroblast"
instances$cell_iname[instances$cell_iname %in% c("CD34")] <- "HSC"

instances$pert_itime[which(instances$cmap_name == "UnTrt")] <- "0 h"

group_lincs <- paste0(instances$cell_iname,";",instances$cmap_name,";",instances$pert_type,";",instances$pert_idose,";",instances$pert_itime)
names(group_lincs) <- instances$sample_id
head(group_lincs)

tab_lincs <- instances[,c("sample_id","cell_iname","cmap_name","pert_type")]
colnames(tab_lincs) <- c("SampleID", "CellType", "Perturbagen","PerturbationType")

tab <- table(unname(group_lincs))
head(sort(tab,decreasing = T))

rna1 <- read.csv("~/Documents/javier_metadata.csv")
rna1$Celltype[rna1$Celltype %in% c(" ARPE-19","ARPE19")] <- "Epithelial Cell"
rna1$Celltype[rna1$Celltype %in% c("") & rna1$Intervention == "CPEB1 "] <- "Fibroblast"
rna1$Celltype[rna1$Celltype %in% c("") & rna1$Intervention == "Myt1L"] <- "Fibroblast"
rna1$Celltype[rna1$Celltype %in% c(" human coronary artery smooth muscle cells", "HASMC","human aortic smooth muscle cells","smooth muscle cells","vascular smooth muscle cell","vascular smooth muscle cells","Vascular Smooth Muscle cells","VSMCs")] <- "Smooth Muscle Cell"
rna1$Celltype[rna1$Celltype %in% c("adipocytes","ASC-adipocytes")] <- "Adipocyte"
rna1$Celltype[rna1$Celltype %in% c("ADSCs","human adipose derived stem cells ")] <- "Adipose-derived Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("alveolar epithelial cells ")] <- "Epithelial Cell"
rna1$Celltype[rna1$Celltype %in% c("ADSCs")] <- "Adipose-derived Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("assumes ev control","myoblasts")] <- "Myoblast"
rna1$Celltype[rna1$Celltype %in% c("astrocyte")] <- "Astrocyte"
rna1$Celltype[rna1$Celltype %in% c("BEAS-2B cells","BEAS2B")] <- "BEAS-2B"
rna1$Celltype[rna1$Celltype %in% c("BMSCs")] <- "Mesenchymal Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("cardiomyocytes (AC16)","human-derived primary cardiomyocytes")] <- "Cardiomyocyte"
rna1$Celltype[rna1$Celltype %in% c("BM-HSC","LT_HSC","Lind-CD34+ cord-blood cells","human hematopoietic stem cells ","human haematopoietic stem cells","human cord blood HSPCs","HUDEP-2 cells and primary human CD34+","HSC","HSC ","HSCs","erythroid differentiation  CD34+ cells","EB-HSPCs ","CD34+ ","CD34+ cells","CD34+ hematopoietic progenitor cells","CD34+","CD34+ hematopoietic stem/progenitor cells","CD34+ HPCs","cord blood-derived CD34+ hematopoietic stem","hematopoietic progenitor cells ","human hematopoietic stem and progenitor cells","human hematopoietic stem cells","HSPCs","human cord blood HSPCs")] <- "HSC"
rna1$Celltype[rna1$Celltype %in% c("CD56plus")] <- "CD56+ Cells"
rna1$Celltype[rna1$Celltype %in% c("Cultured cortical neurons,","GABAergic neurons ","human neurons")] <- "Neuron"
rna1$Celltype[rna1$Celltype %in% c("MRC5","FS087 ","human lung fibroblasts","Rheumatoid arthritis Synovial Fibroblasts","quiescent human fibroblasts","primary fibroblasts from neonatal foreskin ","neonatal dermal fibroblasts cell lines","N78 ","Lung fibroblasts","Kidney fibroblasts","IMR90 fibroblasts","IMR90","human primary fibroblast","Human gingival fibroblast ","human foreskin fibroblast","human fibroblasts","Human fibroblast","Human Dermal Fibroblasts ","human dermal fibroblasts","HFF cells","HDFs","HCF","DCM","dermal fibroblasts","FB","fibroblasts","Fibroblasts","FS087")] <- "Fibroblast"
rna1$Celltype[rna1$Celltype %in% c("dental follicle cells","dental follicle cells ")] <- "Dental Follicle Cell"
rna1$Celltype[rna1$Celltype %in% c("dermal papilla cells")] <- "Dermal Papilla Cell"
rna1$Celltype[rna1$Celltype %in% c("developing mammalian chondrocyte")] <- "Chondrocyte"
rna1$Celltype[rna1$Celltype %in% c("Diff","NSC ","Undiff")] <- "Neural Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("RWPE1","pulmonary epithelial ","Retinal Epithelial Cell","retinal epithelial cells","proliferating human airway epithelial cells","normal human bronchoepithelial ","MC10A","MCF-10","MCF-10A","MCF10a","MCF10A","MCF10A cells","human mammary epithelial cells","human corneal epithelial cell","human airway epithelial cells","HMLE","endometrial epithelium","esophageal epithelial cells","gastric epithelial cells","H9F")] <- "Epithelial Cell"
rna1$Celltype[rna1$Celltype %in% c("Blood Endothelial Cell","huvec","HUVEC","HUVEC ","HUVEC cells ","HUVECs","hPAECs","HAOEC","HAEC","Endothelial","endothelial cells","endothelial cells ","blood endothelial cells","Human Aortic Endothelial Cells","human coronary artery endothelial cells","human dermal lymphatic endothelial cells ","human endothelial cells","Human pulmonary artery endothelial cells","human retinal vascular endothelial cells","human umbilical vein endothelial","Human Umbilical Vein Endothelial Cells","lymphatic endothelial cells","microvascular endothelial","rimary human lung microvascular endothelial cells","vascular endothelial cells")] <- "Endothelial Cell"
rna1$Celltype[rna1$Celltype %in% c("fibroblast-like synoviocytes")] <- "Synoviocytes"
rna1$Celltype[rna1$Celltype %in% c("hiHep")] <- "Hepatocyte"
rna1$Celltype[rna1$Celltype %in% c("MSC","MSCs","hMSC","hMSCs","hMSCs Transcriptomes","hMSC ","human MSCs ","in vivo MSCs","mesenchymal stem cells","Mesenchymal Stem Cells ","BM-MSC")] <- "Mesenchymal Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("hskm","Human Skeletal Muscle Cells")] <- "Skeletal Muscle Cell"
rna1$Celltype[rna1$Celltype %in% c("HUDEP2 cell")] <- "Erythroid Progenitor Cell"
rna1$Celltype[rna1$Celltype %in% c("human brain ")] <- "Brain"
rna1$Celltype[rna1$Celltype %in% c("human CD14+ monocytes")] <- "CD14+ Monocytes"
rna1$Celltype[rna1$Celltype %in% c("Human cord blood cells")] <- "Cord Blood"
rna1$Celltype[rna1$Celltype %in% c("human dental pulp stem cells")] <- "Dental Pulp Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("primary human keratinocytes","human epidermal keratinocytes","neonatal foreskin keratinocytes","neonatal human foreskin keratinocytes ","Human epidermal keratinocytes ","human keratinocytes","KC","keratinocyte (psoriasis)","keratinocytes")] <- "Keratinocyte"
rna1$Celltype[rna1$Celltype %in% c("human islets","Human islets","human islets ")] <- "Islet Cell"
rna1$Celltype[rna1$Celltype %in% c("human Natural Killer ","Natural Killer","NK cells")] <- "Natural Killer Cell"
rna1$Celltype[rna1$Celltype %in% c("human pancreatic β cells ")] <- "Beta Cell"
rna1$Celltype[rna1$Celltype %in% c("human podocyte ","podocyte","Podocyte ","podocytes")] <- "Podocyte"
rna1$Celltype[rna1$Celltype %in% c("human primordial germ cell")] <- "Germ Cell"
rna1$Celltype[rna1$Celltype %in% c("human Sertoli ")] <- "Sertoli Cell"
rna1$Celltype[rna1$Celltype %in% c("human umbilical vein cells")] <- "Umbilical Vein Cell"
rna1$Celltype[rna1$Celltype %in% c("LX2","TWNT4")] <- "Stellate Cell"
rna1$Celltype[rna1$Celltype %in% c("macrophag","macrophages ","human macrophages","Macrophages")] <- "Macrophage"
rna1$Celltype[rna1$Celltype %in% c("mature neutrophil like cells")] <- "Neutrophil"
rna1$Celltype[rna1$Celltype %in% c("mesenchymal stromal cell ")] <- "Mesenchymal Stromal Cell"
rna1$Celltype[rna1$Celltype %in% c("monocyte-derived dendritic cells","monocytes-derived dendritic cells ")] <- "Dendritic Cell"
rna1$Celltype[rna1$Celltype %in% c("Normal esophageal squamous cells,")] <- "Squamous Cell"
rna1$Celltype[rna1$Celltype %in% c("periodontal ligament cells","Human Periodontal Ligament Cells")] <- "Periodontal Ligament Cell"
rna1$Celltype[rna1$Celltype %in% c("Pre-osteoblasts")] <- "Pre-Osteoblasts"
rna1$Celltype[rna1$Celltype %in% c("Primary erythroblasts ")] <- "Erythroblast"
rna1$Celltype[rna1$Celltype %in% c("primary human platelets")] <- "Platelet Cell"
rna1$Celltype[rna1$Celltype %in% c("regenerated human epidermis")] <- "Epidermis"
rna1$Celltype[rna1$Celltype %in% c("s megakaryocytE")] <- "Megakaryocyte"
rna1$Celltype[rna1$Celltype %in% c("s syncytiotrophoblast ")] <- "Syncytiotrophoblast"
rna1$Celltype[rna1$Celltype %in% c("SHED")] <- "Dental Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("SRA01/04 cells")] <- "SRA01/04 Cell"
rna1$Celltype[rna1$Celltype %in% c("tagged, ADPKD patient","WT 9-12 ADPKD patient. ")] <- "Cystic Cell"
rna1$Celltype[rna1$Celltype %in% c("THP1 human","THP1-MD2-CD14 cells ")] <- "Monocyte"
rna1$Celltype[rna1$Celltype %in% c("transformed MSCs")] <- "Transformed Mesenchymal Stem Cell"
rna1$Celltype[rna1$Celltype %in% c("transformed normal breast cells")] <- "Transformed Breast Cell"
rna1$Celltype[rna1$Celltype %in% c("uman c-kit+ cardiac interstitial cells")] <- "Interstitial Cell"
rna1$Celltype[rna1$Celltype %in% c("human skin cells")] <- "Skin"
rna1$Celltype[rna1$Celltype %in% c("T-cell ")] <- "T-cell"
rna1$Celltype[rna1$Celltype %in% c("medial ganglionic eminence (MGE) progenitors")] <- "Medial Ganglion Progenitor Cell"
rna1$Celltype[startsWith(rna1$Celltype,"RA fibroblast")] <- "Skin"

tab_rna1 <- do.call("rbind",apply(rna1,1,function(x){
  sampsExp <- strsplit(x[3],";")[[1]]
  sampsControl <- strsplit(x[4],";")[[1]]
  df <- data.frame(SampleID = c(sampsControl,sampsExp),
                   CellType = x[5],
                   Perturbagen = c(rep("UnTrt",length(sampsControl)),rep(x[2],length(sampsExp))),
                   PerturbationType = c(rep("ctl_untrt",length(sampsControl)),rep("trt_oe",length(sampsExp))),
                   stringsAsFactors = F)
  return(df)
}))

group_rna1 <- unlist(unname(apply(rna1,1,function(x){
  sampsExp <- strsplit(x[3],";")[[1]]
  sampsControl <- strsplit(x[4],";")[[1]]
  map1 <- rep(paste0(x[5],";UnTrt",";ctl_untrt;\"\";0 h"),length(sampsControl))
  names(map1) <- sampsControl
  map2 <- rep(paste0(x[5],";",x[2],";trt_oe;\"\";24 h"),length(sampsExp))
  names(map2) <- sampsExp
  return(c(map1,map2))
})))

group_rna1 <- group_rna1[which(!duplicated(names(group_rna1)))]

rna2 <- read.csv("~/Documents/Sascha_metadata.csv")
rna2$Biosample.Name[rna2$Biosample.Name %in% c("1BR3")] <- "Skin"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Adipose stem cells")] <- "Adipose Stem Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("AoSMC","HASMCs","HCA-SMCs","HCASMC","HUASMCs","HUtSMC","PASMC","Smooth muscle cells","UtSMC","Vascular smooth muscle cells","hVSMC")] <- "Smooth Muscle Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("SAECs","RWPE-1","Primary human bronchial epithelial cells","ARPE-19","Corneal epithelial cells","endometrial epithelial cells","ES-derived RPE cells","HBECs","HCE","HCEC","HMLE","HPrEC","hSAEC","hVEC","MCF10A","MCF10A-HER2-FKBP-HA")] <- "Epithelial Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("ARPE-19 + Hypoxia")] <- "Epithelial Cell + Hypoxia"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("CAR T Cell")] <- "CAR cells"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("CCD-8Lu","Coriell Line ND29510","Coriell Line ND38530","fibroblasts fron neonatal foreskin","HAF","HDF","primary dermal fibroblasts","heterozygous fibroblast","HFF","HFF-1","HLF","homozygous fibroblast","NHLF","IMR-90","IMR-90 (Proliferating)","IMR90","NHDF")] <- "Fibroblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("CD34-derived erythroid progenitors")] <- "Erythroid Progenitor Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("CD34+","CD34+ cells","CD34+ progenitor cells","human primary CD34+","Primary CD34+ hematopoietic stem/progenitor cells (HSPCs)","Lin-CD34+ HSPCs")] <- "HSC"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("CD34+ cells + UNC0638")] <- "HSC + UNC0638"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("CD34+ HSPC-derived proerythroblasts")] <- "Proerythroblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Decidual stromal cells","endometrial stromal cells","Endometrial stromal cells","Primary eutopic stromal cells")] <- "Mesenchymal Stromal Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Dental Pulp Cells")] <- "Dental Pulp Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Valvular endothelial cell","primary human dermal lymphatic endothelial cells","ECFCs","ECs","Endothelial cells","HAEC","HCAECs","HCEnC","HCPECs","HDLEC","HMVEC","HPMVEC","human primary microvascular endothelial cells","Lymphatic endothelial cells","Human Umbilical Vein Endothelial Cell","Huvec","HUVEC","HUVECs")] <- "Endothelial Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Endometrial stromal cells (24h)")] <- "Mesenchymal Stromal Cell (24h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Endometrial stromal cells (48h)")] <- "Mesenchymal Stromal Cell (48h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Endometrial stromal cells (72h)")] <- "Mesenchymal Stromal Cell (72h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("erythroblasts")] <- "Erythroblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Erythroid cells","primary erythroid cells","primary human erythroid cells")] <- "Erythroid"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("FLS","Primary fibroblast-like synoviocytes")] <- "Synoviocytes"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("FLS + Halofuginone")] <- "Synoviocytes + Halofuginone"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("FLS + TNF")] <- "Synoviocytes + TNF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("FLS + TNF + Halofuginone")] <- "Synoviocytes + TNF + Halofuginone"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Fully differentiated adipocytes","White Adipocytes (hWAPC)")] <- "Adipocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary human macrophages","GM-CSF Macrophage","human macrophages","M-CSF-derived macrophages","M0 Macrophages","M1 Macrophages","M2 Macrophages","macrophage","Macrophages","MDMs","Monocyte-derived M-CSF macrophages","PBMC derived macrophage")] <- "Macrophage"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("GM12878")] <- "Lymphoblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("GM12878 + DNMT inhibitors")] <- "Lymphoblast + DNMT inhibitors"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("primary human neonatal keratinocytes","Primary human keratinocytes","HaCaT keratinocytes","human primary keratinocytes","Neonatal foreskin keratinocytes","NHEK","NHK")] <- "Keratinocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hADSC","hASC","human adipose tissue derived stem cells")] <- "Adipose-derived Stem Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hASC + rosiglitazone")] <- "Adipose-derived Stem Cell + rosiglitazone"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HCM")] <- "Myocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hDFCs")] <- "Dental Follicle Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hDPSCs")] <- "Dental Pulp Stem Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HEM","Melanocytes")] <- "Melanocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hepatocyte progenitors")] <- "Hepatocyte Progenitor"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hepatocytes","human hepatocytes")] <- "Hepatocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HFF + doxorubicin")] <- "Fibroblast + doxorubicin"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hMSC","HMSC","Mesenchymal stem cells","MSCs","MSC")] <- "Mesenchymal Stem Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hMSC + IL1B treatment")] <- "Mesenchymal Stem Cell + IL1B"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hMSC + TNF treatment")] <- "Mesenchymal Stem Cell + TNF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HSPC","HSPCs")] <- "HSC"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("human keratinocyte precursor cells")] <- "Keratinocyte Progenitor Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("human motor neurons","Motor Neuron")] <- "Neuron"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + 0h VEGF")] <- "Endothelial Cell + 0h VEGF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + 1h VEGF")] <- "Endothelial Cell + 1h VEGF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + 4h VEGF")] <- "Endothelial Cell + 4h VEGF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + BMP9")] <- "Endothelial Cell + BMP9"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + Fluid Shear Stress")] <- "Endothelial Cell + Fluid Shear Stress"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Huvec + Hypoxia")] <- "Endothelial Cell + Hypoxia"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + IL4")] <- "Endothelial Cell + IL4"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("HUVEC + TNF")] <- "Endothelial Cell + TNF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hVEC + TNF treatment")] <- "Epithelial Cell + TNF"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("iEEC16")] <- "Endometriosis"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("IMR-90 (Quiescent)")] <- "Quiescent Fibroblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("IMR-90 (Senescent)")] <- "Senescent Fibroblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("IMR-90 + Etoposide")] <- "Fibroblast + Etoposide"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("IRM90 + 4OHT")] <- "Fibroblast + 40HT"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Keratinocytes (1h)")] <- "Keratinocyte (1h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Keratinocytes (4h)")] <- "Keratinocyte (4h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Keratinocytes (24h)")] <- "Keratinocyte (24h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Keratinocytes (48)")] <- "Keratinocyte (48h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Lineage-depleted chord blood cells","UCBC")] <- "Cord Blood"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("low passage primary melanoma cultures")] <- "Melanoma"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Lymphatic endothelial cells + Shear Stress")] <- "Endothelial Cell + Shear Stress"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Macrophages + MALP-2 (2h)")] <- "Macrophage + MALP-2 (2h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Macrophages + MALP-2 (4h)")] <- "Macrophage + MALP-2 (4h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Macrophages + MALP-2 (8h)")] <- "Macrophage + MALP-2 (8h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Macrophages + MALP-2 4h)")] <- "Macrophage + MALP-2 (4h)"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("MCF10A + collagen")] <- "Epithelial Cell + collagen"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("monocyte")] <- "Monocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Monocyte-derived M-CSF macrophages + LPS")] <- "Macrophage + LPS"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("monocyte-derived Mo-GMCSF[IL4 (0-72h)] cell")] <- "Macrophage + IL4"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Neonatal foreskin keratinocytes + Adriamycin")] <- "Keratinocyte + Adriamycin"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Neonatal foreskin keratinocytes + Cisplatin")] <- "Keratinocyte + Cisplatin"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Neural precursor cells","NPC","NPCs","U5-NPC (G1 Phase)")] <- "Neural Progenitor Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("neuroepithelial stem cells")] <- "Neuroepithelial Stem Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("neutrophil")] <- "Neutrophil"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("NHA")] <- "Astrocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("NHLF + IL-13")] <- "Fibroblast + IL13"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("NHLF + TGFB")] <- "Fibroblast + TGFB"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("NK")] <- "Natural Killer Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("NK cell progenitors")] <- "Natural Killer Cell Progenitor"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("NSC","NSCs")] <- "Neural Stem Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("oenocytes")] <- "Oenocytes"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("pancreatic islets")] <- "Islet Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("podocytes")] <- "Podocyte"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary CD14-myeloblasts")] <- "Myeloblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary CD4+ T cells (Activated)","Primary CD4+ T cells (Resting)")] <- "CD4+ T cells"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary human macrophages + Hypoxia + IL10")] <- "Macrophage + Hypoxia + IL10"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary human macrophages + Hypoxia")] <- "Macrophage + Hypoxia"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary human macrophages + IL10")] <- "Macrophage + IL10"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary Skeletal Muscle Cells","SkMc")] <- "Skeletal Muscle Cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Primary Skeletal Muscle Cells + Ethanol")] <- "Skeletal Muscle Cell + Ethanol"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("T cell precursors CD1a-")] <- "T Cell Progenitor"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("T cell precursors CD1a+")] <- "T Cell Progenitor"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("T lymphocyte")] <- "T-cell"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Trophoblast cells")] <- "Trophoblast"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("Valvular endothelial cell + OxLdl")] <- "Endothelial Cell + OxLdl"
rna2$Biosample.Name[rna2$Biosample.Name %in% c("hVSMC + TNF treatment")] <- "Smooth Muscle Cell + TNF"

methodmap <- c("ko" = "trt_sh",
               "siRNA" = "trt_sh",
               "shRNA" = "trt_sh",
               "inhibit" = "trt_cp",
               "Overexpression" = "trt_oe",
               "CrISPRi" = "trt_sh",
               "knockin" = "trt_oe",
               "mut" = "trt_mut",
               "CRISPRi" = "trt_sh",
               "mut, knock" = "trt_mut",
               "mut, knock, knock" = "trt_mut",
               "CRISPR" = "trt_sh",
               "shRNA, Overexpression" = "trt_sh",
               "knockdown" = "trt_sh",
               "drug" = "trt_cp",
               "mutant" = "trt_mut",
               "knockout" = "trt_sh",
               "overexpression" = "trt_oe",
               "activemutant" = "trt_mut",
               "druginhibition" = "trt_cp",
               "drugactivation" = "trt_cp")
rna2$PertType <- methodmap[rna2$Method]

tab_rna2 <- do.call("rbind",apply(rna2,1,function(x){
  sampsExp <- strsplit(x[6],";")[[1]]
  sampsControl <- strsplit(x[7],";")[[1]]
  df <- data.frame(SampleID = c(sampsControl,sampsExp),
                   CellType = x[2],
                   Perturbagen = c(rep("UnTrt",length(sampsControl)),rep(x[4],length(sampsExp))),
                   PerturbationType = c(rep("ctl_untrt",length(sampsControl)),rep(x[8],length(sampsExp))),
                   stringsAsFactors = F)
  return(df)
}))

group_rna2 <- unlist(unname(apply(rna2,1,function(x){
  sampsExp <- strsplit(x[6],";")[[1]]
  sampsControl <- strsplit(x[7],";")[[1]]
  map1 <- rep(paste0(x[2],";UnTrt",";ctl_untrt;\"\";0 h"),length(sampsControl))
  names(map1) <- sampsControl
  map2 <- rep(paste0(x[2],";",x[3],";",x[8],";\"\";24 h"),length(sampsExp))
  names(map2) <- sampsExp
  return(c(map1,map2))
})))

group_rna2 <- group_rna2[which(!duplicated(names(group_rna2)))]

group_rna <- c(group_rna1,group_rna2[which(!names(group_rna2) %in% names(group_rna1))])

mic1 <- read.csv("~/Documents/Cmicroarrays.csv", row.names = 1)
mic1 <- mic1[mic1$Species == "Human" & startsWith(mic1$Data_Accession,"GSE") & mic1$Control_GSM != "",]
mic1$Cell_Source[mic1$Cell_Source %in% c("2fTGH (uninfected)")] <- "Fibrosarcoma"
mic1$Cell_Source[mic1$Cell_Source %in% c("Acinar cell","Acinar like cell")] <- "Acinar Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Adipocytes")] <- "Adipocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("Adult pDC","Conventional dendritic cells(neonate)","Dendritic cells","Monocyte-derived denditic cells","Plasmacytoid dendritic cells")] <- "Dendritic Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Airway smooth muscle(expressing GFP)","Umbilical vein smooth muscle cells","Umbilical vein smooth muscle cell","Primary bladder smooth muscle cells","Pulmonary arterial smooth muscle cells","Smooth muscle cell")] <- "Smooth Muscle Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Alpha cell")] <- "Alpha Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Amnion mesenchymal cells")] <- "Mesenchymal Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("HMEC","Nasal epithelial cells","Mammary epithelial cells","HMLE-shGFP breast cells","HMECs","Epithelial cell","Ect1 ectocervical epithelial cells","ARPE-19","Bronchial epithelia","Bronchial Epithelial Cells","Bronchial epithelium cells")] <- "Epithelial Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("B cells")] <- "B Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("BE(2)-C neuronal cells(immature)")] <- "Neuroblast"
mic1$Cell_Source[mic1$Cell_Source %in% c("BEAS-2B cells")] <- "BEAS-2B"
mic1$Cell_Source[mic1$Cell_Source %in% c("Beta cell")] <- "Beta Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("BJ fibroblasts","Synovial fibroblasts","Primary fibroblasts","Normal ovarian fibroblasts cells","Normal colon mucosa-derived fibroblasts","MRC5","Dermal fibroblasts","Fibroblasts","Foreskin fibroblasts","Gingival fibroblast","Gingival fibroblasts","IMR-90 cellslung mesenchymal cells","Lung fibroblasts","Lung Fibroblasts")] <- "Fibroblast"
mic1$Cell_Source[mic1$Cell_Source %in% c("Blood monocyte","Blood monocytes","CD14+ monocytes","Monocytes","Peripheral blood monocytes","Primary monocytes","Primary Monocytes")] <- "Monocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("Brain vascular pericytes")] <- "Pericyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("Cardiomyocytes")] <- "Cardiomyocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("CCR7-CD45RA+ CD8+ T cells")] <- "CD8+ T cells"
mic1$Cell_Source[mic1$Cell_Source %in% c("CD34+ cells","Hematopoietic stem cell","Hematopoietic stem/progenitor cells")] <- "HSC"
mic1$Cell_Source[mic1$Cell_Source %in% c("CD4+ cells","CD4+ T cells")] <- "CD4+ T cells"
mic1$Cell_Source[mic1$Cell_Source %in% c("Decidual cells")] <- "Decidual Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Delta cell")] <- "Delta Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Endothelial cell","Umbilical vein endothelial cells","Retinal endothelial cells","HUVEC","HUVECs","HUVECs (nomoxia)","Immortalized HUVEC cell line EA.hy926","HMVECs","Dermal lymphatic endothelia","Dermal lymphatic endothelial cells","Lung microvascular endothelial cells")] <- "Endothelial Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Ductal cell")] <- "Ductal Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Ecto-cervical cells")] <- "Cervical Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Embryonic kidney cell")] <- "Kidney Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Epidermal keratinocytes","Vaginal Keratinocytes","Skinethic epidermal substitute Keratinocytes","Foreskin epidermal keratinocytes","Foreskin keratinocyte","Immortalized primary keratinocytes (HaCaT)","Keratinocytes","Keratinocytes (HaCaT)","Keratinocytes(differentiated)")] <- "Keratinocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("Epidermal melanocytes")] <- "Melanocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("ESCs")] <- "Embryonic Stem Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Gamma cell")] <- "Gamma Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Hepatocytes","Primary hepatocytes")] <- "Hepatocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("Huh6","Huh7","Huh7.5.1")] <- "Hepatoma"
mic1$Cell_Source[mic1$Cell_Source %in% c("Lamina cribrosa cells")] <- "Lamina Cribrosa Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Leukocytes")] <- "Leukocyte"
mic1$Cell_Source[mic1$Cell_Source %in% c("Lymphoblastoid")] <- "Lymphoblast"
mic1$Cell_Source[mic1$Cell_Source %in% c("Mesenchymal stem cells")] <- "Mesenchymal Stem cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Macrophages","Monocyte-derived GM-CSF macrophages","Monocyte-derived macrophages")] <- "Macrophage"
mic1$Cell_Source[mic1$Cell_Source %in% c("Mononuclear cells")] <- "Mononuclear Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Neural cells")] <- "Neural Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Neutrophils")] <- "Neutrophil"
mic1$Cell_Source[mic1$Cell_Source %in% c("NK cells")] <- "Natural Killer Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Normal primary thyroid cells")] <- "Thyroid Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("NPCs")] <- "Neural Progenitor Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Osteoblasts")] <- "Osteoblast"
mic1$Cell_Source[mic1$Cell_Source %in% c("PBMC","PBMCs")] <- "Blood"
mic1$Cell_Source[mic1$Cell_Source %in% c("Placental cytotrophoblasts")] <- "Trophoblast"
mic1$Cell_Source[mic1$Cell_Source %in% c("Primary microglial cells")] <- "Microglia"
mic1$Cell_Source[mic1$Cell_Source %in% c("Regulatory T cells")] <- "Treg Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Renal stem cells")] <- "Renal Stem Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Skin-resident T cells")] <- "T-cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Stromal cells")] <- "Mesenchymal Stromal Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Umbilical vein cells")] <- "Umbilical Vein Cell"
mic1$Cell_Source[mic1$Cell_Source %in% c("Whole blood")] <- "Blood"

mic1$Duration <- gsub("([0-9])(h)","\\1 \\2",mic1$Duration)
mic1$Concentration[mic1$Concentration == "/"] <- ""
mic1$Concentration <- gsub("([0-9])(mg/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(ng/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(ug/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(mug/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(IU/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(ugram/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(U/ml)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(uM)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(mM)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(nM)","\\1 \\2",mic1$Concentration)
mic1$Concentration <- gsub("([0-9])(ng)","\\1 \\2",mic1$Concentration)

tab_mic1 <- do.call("rbind",apply(mic1,1,function(x){
  sampsExp <- strsplit(x[12],";")[[1]]
  sampsControl <- strsplit(x[13],";")[[1]]
  df <- data.frame(SampleID = c(sampsControl,sampsExp),
                   CellType = x[5],
                   Perturbagen = c(rep("UnTrt",length(sampsControl)),rep(x[2],length(sampsExp))),
                   PerturbationType = c(rep("ctl_untrt",length(sampsControl)),rep("trt_cp",length(sampsExp))),
                   stringsAsFactors = F)
  return(df)
}))

group_mic1 <- unlist(unname(apply(mic1,1,function(x){
  sampsExp <- strsplit(x[12],";")[[1]]
  sampsControl <- strsplit(x[13],";")[[1]]
  map1 <- rep(paste0(x[5],";UnTrt",";ctl_untrt;\"\";0 h"),length(sampsControl))
  names(map1) <- sampsControl
  map2 <- rep(paste0(x[5],";",x[2],";","trt_cp",";\"\";",x[7]),length(sampsExp))
  names(map2) <- sampsExp
  return(c(map1,map2))
})))

group_mic1 <- group_mic1[which(!duplicated(names(group_mic1)))]

mic2 <- read.csv("~/Documents/Smicroarrays.csv",row.names = 1)
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Adipose stem cells","hADSC")] <- "Adipose Stem Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("CD34-derived erythroid progenitors")] <- "Erythroid Progenitor Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("CD34+","CD34+ cells","CD34+ progenitor cells","HSPCs","Primary CD34+ hematopoietic stem/progenitor cells (HSPCs)")] <- "HSC"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Corneal epithelial cells","HCE","MCF10A","MCF10A-HER2-FKBP-HA","Primary human bronchial epithelial cells","SAECs")] <- "Epithelial Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Decidual stromal cells","endometrial stromal cells","Endometrial stromal cells","Primary eutopic stromal cells")] <- "Mesenchymal Stromal Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Dental Pulp Cells")] <- "Dental Pulp Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("ECFCs","Endothelial cells","HDLEC","HMVEC","HUVEC","Lymphatic endothelial cells","primary human dermal lymphatic endothelial cells","Primary lymphatic endothelial cells")] <- "Endothelial Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Endometrial stromal cells (24h)")] <- "Mesenchymal Stromal Cell (24h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Endometrial stromal cells (48h)")] <- "Mesenchymal Stromal Cell (48h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Endometrial stromal cells (72h)")] <- "Mesenchymal Stromal Cell (72h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("fibroblasts fron neonatal foreskin","heterozygous fibroblast","homozygous fibroblast","IMR-90","NHLF")] <- "Fibroblast"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Fully differentiated adipocytes")] <- "Adipocyte"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("GM-CSF Macrophage","Macrophages","Monocyte-derived M-CSF macrophages","Primary human macrophages")] <- "Macrophage"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HaCaT keratinocytes","Neonatal foreskin keratinocytes","NHEK","Primary human keratinocytes","Primary neonatal keratinocytes")] <- "Keratinocyte"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HUASMCs","Smooth muscle cells","Vascular smooth muscle cells")] <- "Smooth Muscle Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("human hepatocytes")] <- "Hepatocyte"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HUVEC + 0h VEGF")] <- "Endothelial Cell + VEGF (0h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HUVEC + 1h VEGF")] <- "Endothelial Cell + VEGF (1h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HUVEC + 4h VEGF")] <- "Endothelial Cell + VEGF (4h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HUVEC + BMP9")] <- "Endothelial Cell + BMP9"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("HUVEC + IL4")] <- "Endothelial Cell + IL4"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("iEEC16")] <- "Endometriosis"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Keratinocytes (24h)")] <- "Keratinocyte (24h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Keratinocytes (1h)")] <- "Keratinocyte (1h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Keratinocytes (4h)")] <- "Keratinocyte (4h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Keratinocytes (48h)")] <- "Keratinocyte (48h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Lineage-depleted chord blood cells")] <- "Cord Blood"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Lymphatic endothelial cells + Shear Stress")] <- "Endothelial Cell + Shear Stress"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Macrophages + MALP-2 (2h)")] <- "Macrophage + MALP-2 (2h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Macrophages + MALP-2 (4h)")] <- "Macrophage + MALP-2 (4h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Macrophages + MALP-2 (8h)")] <- "Macrophage + MALP-2 (8h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Macrophages + MALP-2 4h)")] <- "Macrophage + MALP-2 (4h)"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("MCF10A + collagen")] <- "Endothelial Cell + collagen"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Melanocytes")] <- "Melanocyte"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Mesenchymal stem cells","MSCs")] <- "Mesenchymal Stem Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Neonatal foreskin keratinocytes + Adriamycin")] <- "Keratinocyte + Adriamycin"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Neonatal foreskin keratinocytes + Cisplatin")] <- "Keratinocyte + Cisplatin"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Neural precursor cells")] <- "Neural Progenitor Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("NHLF + IL-13")] <- "Fibroblast + IL13"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("NHLF + TGFB")] <- "Fibroblast + TGFB"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("NSC")] <- "Neural Stem Cell"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Primary CD14-myeloblasts")] <- "Myeloblast"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Primary fibroblast-like synoviocytes")] <- "Synoviocyte"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Primary human macrophages + Hypoxia + IL10")] <- "Macrophage + Hypoxia + IL10"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Primary human macrophages + Hypoxia")] <- "Macrophage + Hypoxia"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Primary human macrophages + IL10")] <- "Macrophage + IL10"
mic2$Biosample.Name[mic2$Biosample.Name %in% c("Trophoblast cells")] <- "Trophoblast"

methodmap <- c("ko" = "trt_sh",
               "siRNA" = "trt_sh",
               "shRNA" = "trt_sh",
               "inhibit" = "trt_cp",
               "Overexpression" = "trt_oe",
               "CrISPRi" = "trt_sh",
               "knockin" = "trt_oe",
               "mut" = "trt_mut",
               "CRISPRi" = "trt_sh",
               "mut, knock" = "trt_mut",
               "mut, knock, knock" = "trt_mut",
               "CRISPR" = "trt_sh",
               "shRNA, Overexpression" = "trt_sh",
               "knockdown" = "trt_sh",
               "drug" = "trt_cp",
               "mutant" = "trt_mut",
               "knockout" = "trt_sh",
               "overexpression" = "trt_oe",
               "activemutant" = "trt_mut",
               "druginhibition" = "trt_cp",
               "drugactivation" = "trt_cp")
mic2$PertType <- methodmap[mic2$Method]

tab_mic2 <- do.call("rbind",apply(mic2,1,function(x){
  sampsExp <- strsplit(x[6],";")[[1]]
  sampsControl <- strsplit(x[7],";")[[1]]
  df <- data.frame(SampleID = c(sampsControl,sampsExp),
                   CellType = x[2],
                   Perturbagen = c(rep("UnTrt",length(sampsControl)),rep(x[3],length(sampsExp))),
                   PerturbationType = c(rep("ctl_untrt",length(sampsControl)),rep(x[8],length(sampsExp))),
                   stringsAsFactors = F)
  return(df)
}))

group_mic2 <- unlist(unname(apply(mic2,1,function(x){
  sampsExp <- strsplit(x[6],";")[[1]]
  sampsControl <- strsplit(x[7],";")[[1]]
  map1 <- rep(paste0(x[2],";UnTrt",";ctl_untrt;\"\";0 h"),length(sampsControl))
  names(map1) <- sampsControl
  map2 <- rep(paste0(x[2],";",x[3],";",x[8],";\"\";24 h"),length(sampsExp))
  names(map2) <- sampsExp
  return(c(map1,map2))
})))

group_mic2 <- group_mic2[which(!duplicated(names(group_mic2)))]

mic3 <- read.csv("~/Documents/Jmicroarrays.csv", row.names = 1)
mic3$Celltype[mic3$Celltype %in% c(" ARPE-19","ARPE19")] <- "Epithelial Cell"
mic3$Celltype[mic3$Celltype %in% c("") & mic3$Intervention == "CPEB1 "] <- "Fibroblast"
mic3$Celltype[mic3$Celltype %in% c("") & mic3$Intervention == "Myt1L"] <- "Fibroblast"
mic3$Celltype[mic3$Celltype %in% c(" human coronary artery smooth muscle cells", "HASMC","human aortic smooth muscle cells","smooth muscle cells","vascular smooth muscle cell","vascular smooth muscle cells","Vascular Smooth Muscle cells","VSMCs")] <- "Smooth Muscle Cell"
mic3$Celltype[mic3$Celltype %in% c("adipocytes","ASC-adipocytes")] <- "Adipocyte"
mic3$Celltype[mic3$Celltype %in% c("ADSCs","human adipose derived stem cells ")] <- "Adipose-derived Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("alveolar epithelial cells ")] <- "Epithelial Cell"
mic3$Celltype[mic3$Celltype %in% c("ADSCs")] <- "Adipose-derived Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("assumes ev control","myoblasts")] <- "Myoblast"
mic3$Celltype[mic3$Celltype %in% c("astrocyte")] <- "Astrocyte"
mic3$Celltype[mic3$Celltype %in% c("BEAS-2B cells","BEAS2B")] <- "BEAS-2B"
mic3$Celltype[mic3$Celltype %in% c("BMSCs")] <- "Mesenchymal Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("cardiomyocytes (AC16)","human-derived primary cardiomyocytes")] <- "Cardiomyocyte"
mic3$Celltype[mic3$Celltype %in% c("BM-HSC","LT_HSC","Lind-CD34+ cord-blood cells","human hematopoietic stem cells ","human haematopoietic stem cells","human cord blood HSPCs","HUDEP-2 cells and primary human CD34+","HSC","HSC ","HSCs","erythroid differentiation  CD34+ cells","EB-HSPCs ","CD34+ ","CD34+ cells","CD34+ hematopoietic progenitor cells","CD34+","CD34+ hematopoietic stem/progenitor cells","CD34+ HPCs","cord blood-derived CD34+ hematopoietic stem","hematopoietic progenitor cells ","human hematopoietic stem and progenitor cells","human hematopoietic stem cells","HSPCs","human cord blood HSPCs")] <- "HSC"
mic3$Celltype[mic3$Celltype %in% c("CD56plus")] <- "CD56+ Cells"
mic3$Celltype[mic3$Celltype %in% c("Cultured cortical neurons,","GABAergic neurons ","human neurons")] <- "Neuron"
mic3$Celltype[mic3$Celltype %in% c("MRC5","FS087 ","human lung fibroblasts","Rheumatoid arthritis Synovial Fibroblasts","quiescent human fibroblasts","primary fibroblasts from neonatal foreskin ","neonatal dermal fibroblasts cell lines","N78 ","Lung fibroblasts","Kidney fibroblasts","IMR90 fibroblasts","IMR90","human primary fibroblast","Human gingival fibroblast ","human foreskin fibroblast","human fibroblasts","Human fibroblast","Human Dermal Fibroblasts ","human dermal fibroblasts","HFF cells","HDFs","HCF","DCM","dermal fibroblasts","FB","fibroblasts","Fibroblasts","FS087")] <- "Fibroblast"
mic3$Celltype[mic3$Celltype %in% c("dental follicle cells","dental follicle cells ")] <- "Dental Follicle Cell"
mic3$Celltype[mic3$Celltype %in% c("dermal papilla cells")] <- "Dermal Papilla Cell"
mic3$Celltype[mic3$Celltype %in% c("developing mammalian chondrocyte")] <- "Chondrocyte"
mic3$Celltype[mic3$Celltype %in% c("Diff","NSC ","Undiff")] <- "Neural Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("RWPE1","pulmonary epithelial ","Retinal Epithelial Cell","retinal epithelial cells","proliferating human airway epithelial cells","normal human bronchoepithelial ","MC10A","MCF-10","MCF-10A","MCF10a","MCF10A","MCF10A cells","human mammary epithelial cells","human corneal epithelial cell","human airway epithelial cells","HMLE","endometrial epithelium","esophageal epithelial cells","gastric epithelial cells","H9F")] <- "Epithelial Cell"
mic3$Celltype[mic3$Celltype %in% c("Blood Endothelial Cell","huvec","HUVEC","HUVEC ","HUVEC cells ","HUVECs","hPAECs","HAOEC","HAEC","Endothelial","endothelial cells","endothelial cells ","blood endothelial cells","Human Aortic Endothelial Cells","human coronary artery endothelial cells","human dermal lymphatic endothelial cells ","human endothelial cells","Human pulmonary artery endothelial cells","human retinal vascular endothelial cells","human umbilical vein endothelial","Human Umbilical Vein Endothelial Cells","lymphatic endothelial cells","microvascular endothelial","rimary human lung microvascular endothelial cells","vascular endothelial cells")] <- "Endothelial Cell"
mic3$Celltype[mic3$Celltype %in% c("fibroblast-like synoviocytes")] <- "Synoviocytes"
mic3$Celltype[mic3$Celltype %in% c("hiHep")] <- "Hepatocyte"
mic3$Celltype[mic3$Celltype %in% c("MSC","MSCs","hMSC","hMSCs","hMSCs Transcriptomes","hMSC ","human MSCs ","in vivo MSCs","mesenchymal stem cells","Mesenchymal Stem Cells ","BM-MSC")] <- "Mesenchymal Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("hskm","Human Skeletal Muscle Cells")] <- "Skeletal Muscle Cell"
mic3$Celltype[mic3$Celltype %in% c("HUDEP2 cell")] <- "Erythroid Progenitor Cell"
mic3$Celltype[mic3$Celltype %in% c("human brain ")] <- "Brain"
mic3$Celltype[mic3$Celltype %in% c("human CD14+ monocytes")] <- "CD14+ Monocytes"
mic3$Celltype[mic3$Celltype %in% c("Human cord blood cells")] <- "Cord Blood"
mic3$Celltype[mic3$Celltype %in% c("human dental pulp stem cells")] <- "Dental Pulp Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("primary human keratinocytes","human epidermal keratinocytes","neonatal foreskin keratinocytes","neonatal human foreskin keratinocytes ","Human epidermal keratinocytes ","human keratinocytes","KC","keratinocyte (psoriasis)","keratinocytes")] <- "Keratinocyte"
mic3$Celltype[mic3$Celltype %in% c("human islets","Human islets","human islets ")] <- "Islet Cell"
mic3$Celltype[mic3$Celltype %in% c("human Natural Killer ","Natural Killer","NK cells")] <- "Natural Killer Cell"
mic3$Celltype[mic3$Celltype %in% c("human pancreatic β cells ")] <- "Beta Cell"
mic3$Celltype[mic3$Celltype %in% c("human podocyte ","podocyte","Podocyte ","podocytes")] <- "Podocyte"
mic3$Celltype[mic3$Celltype %in% c("human primordial germ cell")] <- "Germ Cell"
mic3$Celltype[mic3$Celltype %in% c("human Sertoli ")] <- "Sertoli Cell"
mic3$Celltype[mic3$Celltype %in% c("human umbilical vein cells")] <- "Umbilical Vein Cell"
mic3$Celltype[mic3$Celltype %in% c("LX2","TWNT4")] <- "Stellate Cell"
mic3$Celltype[mic3$Celltype %in% c("macrophag","macrophages ","human macrophages","Macrophages")] <- "Macrophage"
mic3$Celltype[mic3$Celltype %in% c("mature neutrophil like cells")] <- "Neutrophil"
mic3$Celltype[mic3$Celltype %in% c("mesenchymal stromal cell ")] <- "Mesenchymal Stromal Cell"
mic3$Celltype[mic3$Celltype %in% c("monocyte-derived dendritic cells","monocytes-derived dendritic cells ")] <- "Dendritic Cell"
mic3$Celltype[mic3$Celltype %in% c("Normal esophageal squamous cells,")] <- "Squamous Cell"
mic3$Celltype[mic3$Celltype %in% c("periodontal ligament cells","Human Periodontal Ligament Cells")] <- "Periodontal Ligament Cell"
mic3$Celltype[mic3$Celltype %in% c("Pre-osteoblasts")] <- "Pre-Osteoblasts"
mic3$Celltype[mic3$Celltype %in% c("Primary erythroblasts ")] <- "Erythroblast"
mic3$Celltype[mic3$Celltype %in% c("primary human platelets")] <- "Platelet Cell"
mic3$Celltype[mic3$Celltype %in% c("regenerated human epidermis")] <- "Epidermis"
mic3$Celltype[mic3$Celltype %in% c("s megakaryocytE")] <- "Megakaryocyte"
mic3$Celltype[mic3$Celltype %in% c("s syncytiotrophoblast ")] <- "Syncytiotrophoblast"
mic3$Celltype[mic3$Celltype %in% c("SHED")] <- "Dental Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("SRA01/04 cells")] <- "SRA01/04 Cell"
mic3$Celltype[mic3$Celltype %in% c("tagged, ADPKD patient","WT 9-12 ADPKD patient. ")] <- "Cystic Cell"
mic3$Celltype[mic3$Celltype %in% c("THP1 human","THP1-MD2-CD14 cells ")] <- "Monocyte"
mic3$Celltype[mic3$Celltype %in% c("transformed MSCs")] <- "Transformed Mesenchymal Stem Cell"
mic3$Celltype[mic3$Celltype %in% c("transformed normal breast cells")] <- "Transformed Breast Cell"
mic3$Celltype[mic3$Celltype %in% c("uman c-kit+ cardiac interstitial cells")] <- "Interstitial Cell"
mic3$Celltype[mic3$Celltype %in% c("human skin cells")] <- "Skin"
mic3$Celltype[mic3$Celltype %in% c("T-cell ")] <- "T-cell"
mic3$Celltype[mic3$Celltype %in% c("medial ganglionic eminence (MGE) progenitors")] <- "Medial Ganglion Progenitor Cell"
mic3$Celltype[startsWith(mic3$Celltype,"RA fibroblast")] <- "Skin"

tab_mic3 <- do.call("rbind",apply(mic3,1,function(x){
  sampsExp <- strsplit(x[3],";")[[1]]
  sampsControl <- strsplit(x[4],";")[[1]]
  df <- data.frame(SampleID = c(sampsControl,sampsExp),
                   CellType = x[5],
                   Perturbagen = c(rep("UnTrt",length(sampsControl)),rep(x[2],length(sampsExp))),
                   PerturbationType = c(rep("ctl_untrt",length(sampsControl)),rep("trt_oe",length(sampsExp))),
                   stringsAsFactors = F)
  return(df)
}))

group_mic3 <- unlist(unname(apply(mic3,1,function(x){
  sampsExp <- strsplit(x[3],";")[[1]]
  sampsControl <- strsplit(x[4],";")[[1]]
  map1 <- rep(paste0(x[5],";UnTrt",";ctl_untrt;\"\";0 h"),length(sampsControl))
  names(map1) <- sampsControl
  map2 <- rep(paste0(x[5],";",x[2],";","trt_oe",";\"\";24 h"),length(sampsExp))
  names(map2) <- sampsExp
  return(c(map1,map2))
})))

group_mic3 <- group_mic3[which(!duplicated(names(group_mic3)))]

group_mic <- c(group_mic1,group_mic2[which(!names(group_mic2) %in% names(group_mic1))])
group_mic <- c(group_mic,group_mic3[which(!names(group_mic3) %in% names(group_mic))])
group_mic <- group_mic[which(!names(group_mic) %in% names(group_rna))]

load("~/Documents/ArchS4_raw_manCur.RData",verbose = T)
meta_list <- list()
for(gsm in props_train_manCur$GSM){
  print(gsm)
  test <- getGEO(gsm,GSEMatrix=F,getGPL = F)
  meta_list[[gsm]] <- test@header$characteristics_ch1
}

for(i in 1:length(meta_list)){
  write(c(names(meta_list)[i],meta_list[[i]]),file = "~/Documents/test.txt",ncolumns = 22,append = T,sep = "\t")
  #write("\n",file = "~/Documents/test.txt",append = T)
}

###Manual Curation in Excel ###
#Read results and create map
df <- read.table("~/Documents/TrainSamplesToCellTypes.txt", sep = "\t", header = F)

tab_train <-  data.frame(SampleID = df$V1,
                         CellType = df$V2,
                         Perturbagen = "UnTrt",
                         PerturbationType = "ctl_untrt",
                         stringsAsFactors = F)

group_train <- paste0(df$V2,";UnTrt",";ctl_untrt;\"\";0 h")
names(group_train) <- df$V1

group_train <- group_train[which(!names(group_train) %in% names(c(group_lincs,group_rna,group_mic)))]

group_all <- c(group_lincs,group_rna,group_mic,group_train)

save(group_lincs,group_rna,group_mic,group_train,group_all, file = "/vols/GBIAntonio_bigdata/rejuvenationpredata/SampleGroups.RData")

tab_combined <- rbind.data.frame(tab_lincs,tab_rna1,tab_rna2,tab_mic1,tab_mic2,tab_mic3,tab_train)
tab_combined <- tab_combined[!duplicated(tab_combined),]
tab_combined <- tab_combined[!duplicated(tab_combined$SampleID),]
tab_combined <- tab_combined[tab_combined$SampleID %in% rownames(preds),]

save(tab_lincs,tab_rna1,tab_rna2,tab_mic1,tab_mic2,tab_mic3,tab_train,tab_combined, file = "~/Documents/RejuvenationPerturbationPaper/SupplementaryTableS4_raw.RData")
