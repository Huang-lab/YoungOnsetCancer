##### ImmuneFeat_YoungOnset_Assoc.R #####
# William Lee @ September 2019, updated October 2021

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/immune_features"
setwd(bdir)

# read in statistical functions for glm wrapper
source("../stat_functions.R")
system("mkdir out")

library(readxl)
library(tidyverse)

# reads excel spreadsheet into R dataframe
immune_feat_xlfile = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/immune_features/1-s2.0-S1074761318301213-mmc2.xlsx"
immune_feat_df = read_xlsx(immune_feat_xlfile, sheet = NULL, range = NULL, col_names = TRUE,
                          col_types = NULL, na = "NA", trim_ws = TRUE, skip = 0,
                          progress = readxl_progress(), .name_repair = "unique")
unique_immune_samples = unique(substr(immune_feat_df$"TCGA Participant Barcode",1,12))

colnames(immune_feat_df)[1] <- 'bcr_patient_barcode'
colnames(immune_feat_df)[2] <- 'TCGA_Study'
colnames(immune_feat_df)[5:length(colnames(immune_feat_df))] <- gsub('[[:space:]]', '_', colnames(immune_feat_df)[5:length(colnames(immune_feat_df))])
colnames(immune_feat_df)[5:length(colnames(immune_feat_df))] <- gsub('[[:punct:]]', '_', colnames(immune_feat_df)[5:length(colnames(immune_feat_df))])

immune_features = colnames(immune_feat_df)[5:length(colnames(immune_feat_df))]

##### clinical files #####
# clinical file #
clin_complete_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv"
clin_complete = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_complete_f, stringsAsFactors=FALSE)

# Principal Component (PC) file #
PCs_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/2017-04-24_GDAN_AIM_PCA_ethnicity_assigned_WashU.tsv"
PCs = read.table(header=T, quote = "", sep="\t", fill =T, file = PCs_f, stringsAsFactors=FALSE)

clin_brief = clin_complete[,c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender")]
colnames(PCs)[1] = "bcr_patient_barcode"
PCs = PCs[!duplicated(PCs$bcr_patient_barcode),]
clin_merge = merge(clin_brief,PCs, by= "bcr_patient_barcode",all.x=T,all.y=F)

# subtype file #
subtype_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/Subtype\ Assignments.xlsx"
subtype = data.frame(readxl::read_xlsx(subtype_f))
colnames(subtype)[1] = "bcr_patient_barcode"
clin_merge_subtype = merge(clin_merge,subtype, by= "bcr_patient_barcode")

# create binary age variable; typically <=50 considered as young onset cancer
clin_merge_subtype$age_binary = clin_merge_subtype$age_at_initial_pathologic_diagnosis <= 50

# only include the ones with available molecular data
clin_merge_subtype_avail = clin_merge_subtype[clin_merge_subtype$bcr_patient_barcode %in% unique_immune_samples,]
clin_merge_subtype_avail_complete = clin_merge_subtype_avail[complete.cases(clin_merge_subtype_avail[,2:4]) & clin_merge_subtype_avail$gender != "",]

# quick check of merged data
cat("TCGA cases count available for this analysis\n")
nrow(clin_merge_subtype_avail_complete)
cat("Cancer type distribution of TCGA cases with young onset cancer\n")
table(clin_merge_subtype_avail_complete$age_at_initial_pathologic_diagnosis<=50,clin_merge_subtype_avail_complete$acronym)

# cleaning of merged df
clin_merge_subtype_avail_complete$cancer = NULL
clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

# initiate empty list and index = 1 to store results
immune_results_list = as.list(NULL)
immune_results_list_binary = as.list(NULL)
i=1
# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)) {
  
  clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
  clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
  
  if ((sum(clin_merge_subtype_avail_complete_c$age_binary==TRUE) >= 40) &
      (sum(clin_merge_subtype_avail_complete_c$age_binary==FALSE) >= 40)) {
  
  # conduct the test by immune feature
  for (immune_feat in immune_features) {
    
   if (immune_feat == "Th1_Cells" | immune_feat == "Th2_Cells" | immune_feat == "Th17_Cells") {
     
   # if (immune_feat == "Lymphocyte_Infiltration_Signature_Score" | immune_feat == "Wound_Healing" | immune_feat == "Proliferation" |
   #     immune_feat == "Macrophage_Regulation" | immune_feat == "IFN_gamma_Response" | immune_feat == "TGF_beta_Response") {
     
   # if (immune_feat == "Macrophages" | immune_feat == "Macrophages_M1" | immune_feat == "Macrophages_M2" |
   #     immune_feat == "Dendritic_Cells" | immune_feat == "B_Cells_Naive" | immune_feat == "T_Cells_CD4_Naive" |
   #     immune_feat == "T_Cells_CD8" | immune_feat == "NK_Cells_Activated" | immune_feat == "NK_Cells_Resting" |
   #     immune_feat == "T_Cells_Regulatory_Tregs") {
    
   # if (immune_feat == "Indel_Neoantigens" | immune_feat == "SNV_Neoantigens") {
   # if (cancer == "SKCM" & immune_feat == "SNV_Neoantigens") {next}
   # if (cancer == "SKCM" & immune_feat == "Indel_Neoantigens") {next}
   # if (cancer == "OV" & immune_feat == "Indel_Neoantigens") {next}
    
   #if (immune_feat == "Silent_Mutation_Rate" | immune_feat == "Nonsilent_Mutation_Rate") {
       
      immune_feat_df_current = select(immune_feat_df, bcr_patient_barcode, immune_feat)
      immune_feat_df_current = na.omit(immune_feat_df_current)

      clin_merge_subtype_avail_complete_c_merge_if = merge(clin_merge_subtype_avail_complete_c,immune_feat_df_current, by= "bcr_patient_barcode", all.x = T, all.y = F)
      clin_merge_subtype_avail_complete_c_merge_if = clin_merge_subtype_avail_complete_c_merge_if[complete.cases(clin_merge_subtype_avail_complete_c_merge_if[29]),]

      if (dim(clin_merge_subtype_avail_complete_c_merge_if)[1] == 0) {next}

      if (all(clin_merge_subtype_avail_complete_c_merge_if$age_binary) | 
          any(clin_merge_subtype_avail_complete_c_merge_if$age_binary) == FALSE) {next}
      
      else {
        
      # model onset age as a linear variable
      model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_if, yi = immune_feat, xi = "age_at_initial_pathologic_diagnosis", ytype = "Continuous", covi = c("SUBTYPE","gender","PC1","PC2"))
      cancer_stat = data.frame(cbind(cancer, immune_feat, model_results))
      immune_results_list[[i]] = cancer_stat
        
      # model onset age as a binary variable
      model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_if, yi = immune_feat, xi = "age_binary", ytype = "Continuous", covi = c("SUBTYPE","gender","PC1","PC2"))
      cancer_stat = data.frame(cbind(cancer, immune_feat, model_results))
      immune_results_list_binary[[i]] = cancer_stat
      
      }
    }
    
    # increment index
    i = i + 1
  }
 }
}

# compile linear result
tt = do.call(rbind,immune_results_list)
colnames(tt) = c("cancer","immune_feature","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F", "Pr_F","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"Pr_F"], method="fdr") 
tt=tt[order(tt$Pr_F, decreasing=FALSE),]
tn = "out/Th_Cells_onsetAge_assoc.txt"
#tn = "out/IGS_onsetAge_assoc.txt"
#tn = "out/IIC_onsetAge_assoc.txt"
#tn = "out/Neoantigens_onsetAge_assoc.txt"
#tn = "out/Mut_Rate_onsetAge_assoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

# compile binary result
tt = do.call(rbind,immune_results_list_binary)
colnames(tt) = c("cancer","immune_feature","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F", "Pr_F","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"Pr_F"], method="fdr") 
tt=tt[order(tt$Pr_F, decreasing=FALSE),]
tn = "out/Th_Cells_onsetAgeBinary_assoc.txt"
#tn = "out/IGS_onsetAgeBinary_assoc.txt"
#tn = "out/IIC_onsetAgeBinary_assoc.txt"
#tn = "out/Neoantigens_onsetAgeBinary_assoc.txt"
#tn = "out/Mut_Rate_onsetAgeBinary_assoc.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)