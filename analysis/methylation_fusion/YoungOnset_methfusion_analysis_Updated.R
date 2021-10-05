##### YoungOnset_methfusion_analysis_Updated.R #####
# Updated by Will Lee, October 2021

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion"
setwd(bdir)

# read in statistical functions for glm wrapper
source("~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/stat_functions_MF.R")
system("mkdir out")

library(readxl)
library(tidyverse)

onc_sig_file = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/oncosigfilefixed.txt"
oncosig_df = read.table(sep="\t",header=T, file=onc_sig_file , stringsAsFactors=FALSE)

oncosig_df$SAMPLE_BARCODE <- substr(oncosig_df$SAMPLE_BARCODE, 1,12) # make sample barcode into patient barcode

names(oncosig_df)[1] <- "bcr_patient_barcode"
oncosig_df$bcr_patient_barcode= as.character(oncosig_df$bcr_patient_barcode)
unique_oncosig_samples <- unique(oncosig_df$bcr_patient_barcode) # 9125

meth_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("EPISIL"))) # 17
names(meth_df)[1]='bcr_patient_barcode'
fusion_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("FUSION"))) # 125
names(fusion_df)[1]='bcr_patient_barcode'
cnv_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("AMP")), select(oncosig_df, starts_with("DEL"))) # 103
names(cnv_df)[1]='bcr_patient_barcode'
mut_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("MUT"))) # 169
names(mut_df)[1]='bcr_patient_barcode'

## clinical files
clin_complete_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv" #for analysis locally
clin_complete = read.table(header=T, quote = "", sep="\t", fill=T, file = clin_complete_f, stringsAsFactors=FALSE)

## PCA files
PCs_f = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/2017-04-24_GDAN_AIM_PCA_ethnicity_assigned_WashU.tsv" #for analysis locally
PCs = read.table(header=T, quote = "", sep="\t", fill =T, file = PCs_f, stringsAsFactors=FALSE)

clin_brief = clin_complete[,c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender")]
colnames(PCs)[1] = "bcr_patient_barcode"
PCs = PCs[!duplicated(PCs$bcr_patient_barcode),]
clin_merge = merge(clin_brief,PCs, by= "bcr_patient_barcode",all.x=T,all.y=F) # 10956 x 28

subtype_f = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/Subtype\ Assignments.xlsx" #for analysis locally
subtype = data.frame(readxl::read_xlsx(subtype_f))
colnames(subtype)[1] = "bcr_patient_barcode"
clin_merge_subtype = merge(clin_merge,subtype, by= "bcr_patient_barcode") # 8943 x 31

clin_merge_subtype$age_binary = clin_merge_subtype$age_at_initial_pathologic_diagnosis <= 50

# only include the ones with available data in oncosig file
clin_merge_subtype_avail = clin_merge_subtype[clin_merge_subtype$bcr_patient_barcode %in% unique_oncosig_samples,] # 8943 x 32 
clin_merge_subtype_avail_complete = clin_merge_subtype_avail[complete.cases(clin_merge_subtype_avail[,2:4]) & clin_merge_subtype_avail$gender != "",] # 8943 x 32

# run line depending on which type you are looking at (meth, fusion, mut, cnv) 
alt_events = names(meth_df)[2:length(names(meth_df))] # for meth
# alt_events = names(fusion_df)[2:length(names(fusion_df))] # for fusion
# alt_events = names(cnv_df)[2:length(names(cnv_df))] # for cnv
# alt_events = names(mut_df)[2:length(names(mut_df))] # for mut

# cleaning of merged df
clin_merge_subtype_avail_complete$cancer = NULL
clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

# initiate empty list and index = 1 to store results
results_list_oncosig = as.list(NULL)
results_list_oncosig_binary = as.list(NULL)
i=1
# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
    
  clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
  clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
  
  # only proceed if cancer has >= 40 young adult & >= 40 later-onset cases post-cleaning
  if ((sum(clin_merge_subtype_avail_complete_c$age_binary==TRUE) >= 40) &
      (sum(clin_merge_subtype_avail_complete_c$age_binary==FALSE) >= 40)) {

  # conduct the test by each alteration event
  for (mf in alt_events){
    
    mf_df_current = select(meth_df, bcr_patient_barcode, mf) # change line to specify which type you are looking at (meth, fusion, mut, cnv) 
    mf_df_current = mf_df_current[mf_df_current[,mf] != 'NA',]
    if (!is.logical(mf_df_current[,mf])) {
      mf_df_current[,mf] = as.logical(mf_df_current[,mf])
    }
    
    clin_merge_subtype_avail_complete_c_merge_mf = merge(clin_merge_subtype_avail_complete_c, mf_df_current, by= "bcr_patient_barcode", all.x=T, all.y = F)
    clin_merge_subtype_avail_complete_c_merge_mf = clin_merge_subtype_avail_complete_c_merge_mf[clin_merge_subtype_avail_complete_c_merge_mf[,mf] != 'NA',]
    
    if (sum(clin_merge_subtype_avail_complete_c_merge_mf[,mf]==TRUE, na.rm=TRUE) <5) {next}
    
    model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_mf, yi = mf, xi = "age_at_initial_pathologic_diagnosis", ytype = "Binary", covi = c("SUBTYPE","gender","PC1","PC2"))
    cancer_stat = data.frame(cbind(cancer, mf, model_results))
    results_list_oncosig[[i]] = cancer_stat
    
    ## model onset age as a binary variable
    model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_mf, yi = mf, xi = "age_binary", ytype = "Binary", covi = c("SUBTYPE","gender","PC1","PC2"))
    cancer_stat = data.frame(cbind(cancer, mf, model_results))
    results_list_oncosig_binary[[i]] = cancer_stat
    
    # increment index
    i = i + 1
  }  
 }
}

##### Compile and store results #####
# make sure to change col names (to meth_event, cnv_event, etc.) to reflect which analysis is being run

# compile linear result
tt = do.call(rbind,results_list_oncosig)
colnames(tt) = c("cancer","cnv_event","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                  "P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="fdr") 
tt=tt[order(tt$P_Chi, decreasing=FALSE),]
tn = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/out/meth_events_onsetAge_assoc_Updated.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)

# compile binary result
tt = do.call(rbind,results_list_oncosig_binary)
colnames(tt) = c("cancer","cnv_event","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                  "P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="fdr") 
tt=tt[order(tt$P_Chi, decreasing=FALSE),]
tn = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/out/meth_events_onsetAgeBinary_assoc_Updated.txt"
write.table(tt, quote=F, sep="\t", file = tn, row.names = F)