##### YoungAdult_Amplification_Clin_Assoc.R #####
# William Lee @ January 2020, updated October 2021

library(readxl)
library(tidyverse)
library(stringr)

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion"
setwd(bdir)

onc_sig_file = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/methylation_fusion/oncosigfilefixed.txt"
oncosig_df = read.table(sep="\t",header=T, file=onc_sig_file , stringsAsFactors=FALSE)

oncosig_df$SAMPLE_BARCODE <- substr(oncosig_df$SAMPLE_BARCODE, 1,12) # make sample barcode into patient barcode

names(oncosig_df)[1] <- "bcr_patient_barcode"
oncosig_df$bcr_patient_barcode= as.character(oncosig_df$bcr_patient_barcode)
unique_oncosig_samples <- unique(oncosig_df$bcr_patient_barcode) # 9125

amp_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("AMP"))) 
names(amp_df)[1]='bcr_patient_barcode'

## clinical files
clin_complete_f = '~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv' #for analysis locally
clin_complete = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_complete_f, stringsAsFactors=FALSE)

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

amp_events = names(amp_df)[2:length(names(amp_df))]

clin_merge_subtype_avail_complete$cancer = NULL; clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL; clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL

clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="BRCA"] = "1612" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="CESC"] = "4362" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="COAD"] = "9256" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="HNSC"] = "5520" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="KIRC"] = "4467" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="KIRP"] = "4465" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="LGG"]  = "60108" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="LIHC"] = "684" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="OV"]   = "2394" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="PCPG"] = "50773" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="SARC"] = "1115" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="SKCM"] = "8923" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="THCA"] = "1781" ###
clin_merge_subtype_avail_complete$doid[clin_merge_subtype_avail_complete$acronym=="UCEC"] = "1380"###

### analytic pipeline ###

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/clinical_associations" ###
setwd(bdir)

# civic_raw_data updated 06/04/2020
civic_raw_data = "01-May-2020-ClinicalEvidenceSummaries.tsv" # 3431 obs. of 41 variables
civic_df = read.table(header=T, quote = "", sep="\t", fill=T, file = civic_raw_data, stringsAsFactors=FALSE)

### civic data preparation 
civic_df = civic_df[str_detect(civic_df$variant, "AMPLIFICATION"),]

civic_df$ampEvnt = paste(civic_df$gene, "AMP" ,sep=":")
civic_df$doid[is.na(civic_df$doid)] = "XXXXX"
### end of civic data preparation 

# cgi_raw_data updated 06/04/2020
cgi_raw_data = "cgi_biomarkers_per_variant.tsv"
cgi_df = read.table(header=T, quote = "", sep="\t", fill=T, file = cgi_raw_data, stringsAsFactors=FALSE)

### cgi data preparation
cgi_df = cgi_df[str_detect(cgi_df$Biomarker, "amplification"),]

cgi_df$Evidence.level[cgi_df$Evidence.level=="NCCN guidelines"] = "A"
cgi_df$Evidence.level[cgi_df$Evidence.level=="FDA guidelines"] = "A" 
cgi_df$Evidence.level[cgi_df$Evidence.level=="European LeukemiaNet guidelines"] = "A"
cgi_df$Evidence.level[cgi_df$Evidence.level=="NCCN/CAP guidelines"] = "A"
cgi_df$Evidence.level[cgi_df$Evidence.level=="CPIC guidelines"] = "A"

cgi_df$Evidence.level[cgi_df$Evidence.level=="Clinical trials"] = "B"
cgi_df$Evidence.level[cgi_df$Evidence.level=="Late trials"] = "B"
cgi_df$Evidence.level[cgi_df$Evidence.level=="Early trials"] = "B"
cgi_df$Evidence.level[cgi_df$Evidence.level=="Early Trials,Case Report"] = "B"

cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Glioma")] = "60108"
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "glioma")] = "60108" 
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Sarcoma")] = "1115"
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "sarcoma")] = "1115"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Breast adenocarcinoma")] = "1612"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Cervix")] = "4362"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Colorectal adenocarcinoma")] = "9256"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Head an neck squamous")] = "5520"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Hepatic carcinoma")] = "684"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Ovary")] = "2394"
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Cutaneous melanoma")] = "8923" 
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Thyroid")] = "1781" 
#
cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, "Endometrium")] = "1380"

cgi_df$doid[str_detect(cgi_df$Primary.Tumor.type, ";")] = "XXXXX"
cgi_df$doid[is.na(cgi_df$doid)] = "XXXXX"

cgi_df$Drug = str_replace_all(cgi_df$Drug, fixed("["), "")
cgi_df$Drug = str_replace_all(cgi_df$Drug, fixed("]"), "")

cgi_df$ampEvnt = paste(cgi_df$Gene, "AMP", sep=":")
### end of cgi data preparation 

# onco_raw_data updated 06/04/2020
onco_raw_data = "oncokb_biomarker_drug_associations.tsv"
onco_df = read.table(header=T, quote = "", sep="\t", fill=T, file = onco_raw_data, stringsAsFactors=FALSE)

### onco kb data preparation
onco_df = onco_df[str_detect(onco_df$Alterations, "Amplification"),]

onco_df$Evidence.Level[onco_df$Level=="1"] = "A"
onco_df$Evidence.Level[onco_df$Level=="2"] = "A" 
onco_df$Evidence.Level[onco_df$Level=="3"] = "B"

onco_df$Evidence.Level[onco_df$Level=="R1"] = "A" 
onco_df$Evidence.Level[onco_df$Level=="R2"] = "B"

onco_df$doid[str_detect(onco_df$Tumor.Type, "Glioma")] = "60108" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Breast Cancer")] = "1612" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Colorectal Cancer")] = "9256" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Head and Neck Squamous Cell Carcinoma")] = "5520" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Melanoma")] = "8923" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Liposarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Ewing Sarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Ewing Sarcoma of Soft Tissue")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Dedifferentiated Liposarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Well-Differentiated Liposarcoma")] = "1115" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Dermatofibrosarcoma Protuberans")] = "1115" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Low-Grade Serous Ovarian Cancer")] = "2394" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Ovarian Cancer")] = "2394" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Anaplastic Thyroid Cancer")] = "1781" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Thyroid Cancer")] = "1781" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Medullary Thyroid Cancer")] = "1781" ###
#
onco_df$doid[str_detect(onco_df$Tumor.Type, "Uterine Serous Carcinoma/Uterine Papillary Serous Carcinoma")] = "1380" ###
onco_df$doid[str_detect(onco_df$Tumor.Type, "Endometrial Cancer")] = "1380" ###

onco_df$doid[is.na(onco_df$doid)] = "XXXXX"

onco_df$ampEvnt = paste(onco_df$Gene, "AMP", sep=":")

onco_df$clin.sig[onco_df$Level=="1" | onco_df$Level=="2" | onco_df$Level=="3"] = "Responsive"
onco_df$clin.sig[onco_df$Level=="R1" | onco_df$Level=="R2"] = "Resistant"
### end of onco kb data preparation

##### the following code generates combination df fed into analytic pipeline ##### 
doid1 = civic_df$doid
doid2 = cgi_df$doid
doid3 = onco_df$doid
doid = c(doid1, doid2, doid3)

drug1 = civic_df$drugs
drug2 = cgi_df$Drug
drug3 = onco_df$Drugs
drugs = c(drug1, drug2, drug3)

ampEvnt1 = civic_df$ampEvnt
ampEvnt2 = cgi_df$ampEvnt
ampEvnt3 = onco_df$ampEvnt
ampEvnt = c(ampEvnt1, ampEvnt2, ampEvnt3)

elvl1 = civic_df$evidence_level
elvl2 = cgi_df$Evidence.level
elvl3 = onco_df$Evidence.Level
evidence_level = c(elvl1, elvl2, elvl3)

clin1 = civic_df$clinical_significance
clin2 = cgi_df$Association
clin3 = onco_df$clin.sig
clinical_significance = c(clin1, clin2, clin3)

poten_clin_act_3db = as.data.frame(cbind(doid, drugs, ampEvnt, evidence_level, clinical_significance))

######################################################################################################

unique_ampevnt_3db = unique(poten_clin_act_3db$ampEvnt)

clin_amplification_therapeutics <- data.frame(matrix(ncol = 32, nrow = 0))
colnames(clin_amplification_therapeutics) <- colnames(clin_merge_subtype_avail_complete)
colnames(clin_amplification_therapeutics)[30:32] <- c("amp_binary","ampEvnt","ab_drug")

for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    if ((sum(clin_merge_subtype_avail_complete_c$age_binary==TRUE) >= 40) &
       (sum(clin_merge_subtype_avail_complete_c$age_binary==FALSE) >= 40)) {
    
    # conduct the test by each alteration event
    for (amp in amp_events){
      
      amp_df_current = select(amp_df, bcr_patient_barcode, amp) ###
      amp_df_current = amp_df_current[amp_df_current[,amp] != 'NA',] ###
      
      clin_merge_subtype_avail_complete_c_merge_amp = merge(clin_merge_subtype_avail_complete_c, amp_df_current, by= "bcr_patient_barcode", all.x=T, all.y=F)
      clin_merge_subtype_avail_complete_c_merge_amp = clin_merge_subtype_avail_complete_c_merge_amp[clin_merge_subtype_avail_complete_c_merge_amp[,amp] != 'NA',]
      colnames(clin_merge_subtype_avail_complete_c_merge_amp)[30] <- "amp_binary"
      
      clin_merge_subtype_avail_complete_c_merge_amp$ampEvnt[clin_merge_subtype_avail_complete_c_merge_amp$amp_binary==TRUE] = paste((unlist(strsplit(amp, "[.]")))[2], "AMP", sep=":")
      clin_merge_subtype_avail_complete_c_merge_amp = clin_merge_subtype_avail_complete_c_merge_amp[complete.cases(clin_merge_subtype_avail_complete_c_merge_amp),]
      
      # skips to next loop iteration if clin_merge_subtype_avail_complete_c_merge_amp has no obs.
      if (dim(clin_merge_subtype_avail_complete_c_merge_amp)[1] == 0) {next}
      
      current_ampevnt = unique(clin_merge_subtype_avail_complete_c_merge_amp$ampEvnt)
      
      if (current_ampevnt %in% unique_ampevnt_3db) {
          
        # poten_clin_act_evnt contains only rows that have amplification event of interest 
        poten_clin_act_evnt = poten_clin_act_3db[poten_clin_act_3db$ampEvnt==current_ampevnt,]
        
        # replaces blank cells in poten_clin_act_evnt with NA (so complete.cases() will work)
        poten_clin_act_evnt[ poten_clin_act_evnt == "" ] <- NA
        
        # confirm that all rows in poten_clin_act_evnt have evidence levels and drugs
        poten_clin_act_evnt = poten_clin_act_evnt[complete.cases(poten_clin_act_evnt$evidence_level),]
        poten_clin_act_evnt = poten_clin_act_evnt[complete.cases(poten_clin_act_evnt$drugs),]
        
        # skips if poten_clin_act_evnt has no observations
        if (dim(poten_clin_act_evnt)[1] == 0) {
          clin_merge_subtype_avail_complete_c_merge_amp$ab_drug = ""
          clin_amplification_therapeutics <- rbind(clin_amplification_therapeutics, clin_merge_subtype_avail_complete_c_merge_amp)
          next
        }
        
        # instantiates empty list (elements to be appended)
        current_list = list()
          
        # initiates counter at 1
        i = 1
        
        for (elvl in poten_clin_act_evnt$evidence_level) {
          
          if (elvl == "A" & (poten_clin_act_evnt[i, "clinical_significance"] == "Sensitivity/Response" |
                             poten_clin_act_evnt[i, "clinical_significance"] == "Responsive")) {
            
            if (poten_clin_act_evnt[i, "doid"] == unique(clin_merge_subtype_avail_complete_c_merge_amp$doid)) {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "A on-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            else {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "A off-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            i = i + 1
            
          }
          
          else if (elvl == "B" & (poten_clin_act_evnt[i, "clinical_significance"] == "Sensitivity/Response" |
                                  poten_clin_act_evnt[i, "clinical_significance"] == "Responsive")) {
            
            if (poten_clin_act_evnt[i, "doid"] == unique(clin_merge_subtype_avail_complete_c_merge_amp$doid)) {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "B on-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            else {
              current_drug = paste(poten_clin_act_evnt[i, "drugs"], "B off-label", sep=": ")
              current_list = append(current_list, current_drug)
            }
            
            i = i + 1
            
          }
          
          else {i = i + 1}
          
        }
        
        current_list = unique(current_list)
        therapeutic_drugs = paste(current_list, collapse = " & ")
        clin_merge_subtype_avail_complete_c_merge_amp$ab_drug = therapeutic_drugs
        
        clin_amplification_therapeutics <- rbind(clin_amplification_therapeutics, clin_merge_subtype_avail_complete_c_merge_amp)
        
      }
        
      else {
        
        clin_merge_subtype_avail_complete_c_merge_amp$ab_drug = ""
        clin_amplification_therapeutics <- rbind(clin_amplification_therapeutics, clin_merge_subtype_avail_complete_c_merge_amp)
        
      }
    }
  }
}

# get df for actionable ampEvnt analysis later (need duplicates) #
clin_amplification_actionable = clin_amplification_therapeutics
clin_amplification_actionable[clin_amplification_actionable == ""] <- NA
clin_amplification_actionable = clin_amplification_actionable[!is.na(clin_amplification_actionable$ab_drug),]
clin_amplification_actionable = clin_amplification_actionable[
clin_amplification_actionable$acronym=="BRCA" | clin_amplification_actionable$acronym=="CESC" |
clin_amplification_actionable$acronym=="HNSC" | clin_amplification_actionable$acronym=="LGG" |
clin_amplification_actionable$acronym=="OV" | clin_amplification_actionable$acronym=="SKCM",]

duplicated_patient_samples = clin_amplification_therapeutics[duplicated(clin_amplification_therapeutics$bcr_patient_barcode),1]

for (duplicate in duplicated_patient_samples) {
  
  current_duplicates_df = clin_amplification_therapeutics[clin_amplification_therapeutics$bcr_patient_barcode==duplicate,]
  clin_amplification_therapeutics = clin_amplification_therapeutics[clin_amplification_therapeutics$bcr_patient_barcode!=duplicate,]
  
  drugs_as_string = paste(current_duplicates_df$ab_drug, collapse = " ")
  
  if (str_detect(drugs_as_string, "A on-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "A on-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "A off-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "A off-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "B on-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "B on-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else if (str_detect(drugs_as_string, "B off-label")) {
    current_duplicates_df = current_duplicates_df[str_detect(current_duplicates_df$ab_drug, "B off-label"),]
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
  else {
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics, current_duplicates_df[1,])
  }
  
}

clin_merge_amp_temp = data.frame(matrix(ncol = 29, nrow = 0))
colnames(clin_merge_amp_temp) <- colnames(clin_merge_subtype_avail_complete)

for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  temp = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
  temp = temp[complete.cases(temp),]
  
  if (nrow(temp) >= 40) {
    
    clin_merge_amp_temp = rbind(clin_merge_amp_temp,temp)
    
  }
}

clin_merge_amp_temp["amp_binary"] = NA; clin_merge_amp_temp["ampEvnt"] = NA; clin_merge_amp_temp["ab_drug"] = NA

for (barcode in clin_merge_amp_temp$bcr_patient_barcode) {
  if (!(barcode %in% clin_amplification_therapeutics$bcr_patient_barcode)) {
    clin_amplification_therapeutics = rbind(clin_amplification_therapeutics,
                                            clin_merge_amp_temp[clin_merge_amp_temp$bcr_patient_barcode==barcode,])
  }
}

### the following code produces the ggplot-ready dataframe ###

# replaces blank cells in clin_amplification_therapeutics with NA   
clin_amplification_therapeutics[clin_amplification_therapeutics == ""] <- NA

clin_amplification_therapeutics_ggplot <- data.frame(matrix(ncol = 32, nrow = 0))
colnames(clin_amplification_therapeutics_ggplot) <- colnames(clin_amplification_therapeutics)

for (cancer in unique(clin_amplification_therapeutics$acronym)){
  
  i = 1
  
  clin_amplification_therapeutics_c = clin_amplification_therapeutics[clin_amplification_therapeutics$acronym==cancer,]
  
  total_reg_cases = sum(clin_amplification_therapeutics_c$age_binary==FALSE & !is.na(clin_amplification_therapeutics_c$ab_drug))
  total_yng_cases = sum(clin_amplification_therapeutics_c$age_binary==TRUE & !is.na(clin_amplification_therapeutics_c$ab_drug))
  
  if (total_reg_cases >= 10 & total_yng_cases >= 10) {
    
    for (element in clin_amplification_therapeutics_c$ab_drug) {
      
      if (is.na(element)) {i = i + 1}
      
      else if (str_detect(element, "A on-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "A on-label"; i = i + 1}
      
      else if (str_detect(element, "A off-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "A off-label"; i = i + 1}
      
      else if (str_detect(element, "B on-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "B on-label"; i = i + 1}
      
      else if (str_detect(element, "B off-label")) {clin_amplification_therapeutics_c[i, "ab_drug"] =  "B off-label"; i = i + 1}
    }
    
    clin_amplification_therapeutics_ggplot <- rbind(clin_amplification_therapeutics_ggplot, clin_amplification_therapeutics_c)
    
  }
}

clin_amplification_therapeutics_ggplot$ab_drug[is.na(clin_amplification_therapeutics_ggplot$ab_drug)] = "None"

clin_amplification_therapeutics_ggplot_rmNone = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$ab_drug!="None",]
clin_amplification_therapeutics_ggplot_None = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$ab_drug=="None",]

clin_amplification_therapeutics_ggplot$ab_drug = factor(clin_amplification_therapeutics_ggplot$ab_drug, levels=c("None", "B off-label", "B on-label", "A off-label", "A on-label"))

clin_amplification_therapeutics_ggplot$plot_age[clin_amplification_therapeutics_ggplot$age_binary==TRUE] = "  50"
clin_amplification_therapeutics_ggplot$plot_age[clin_amplification_therapeutics_ggplot$age_binary==FALSE] = "> 50"

################################ RESPONSIVE FIGURES ################################ 

###GGPLOT DRUG PERCENT###
p = ggplot(data=clin_amplification_therapeutics_ggplot, aes(x = plot_age, alpha = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "fill", fill = "dodgerblue4")
p = p + labs(y = "cases with treatment option(s) (%)", x = "age at pathologic diagnosis (yrs.)",
             alpha = "highest predicted\nactionability level") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + scale_alpha_manual(values=c("A on-label" = 1, "A off-label" = 0.7, "B on-label" = 0.35, "B off-label" = 0.1, "None" = 0),
                           breaks=c("A on-label", "A off-label", "B on-label", "B off-label"))
p = p + theme(legend.position = "top")
p
fn = "out/YoungAdult_Amplification_ClinAct_BarplotBinary_Percent.pdf"
ggsave(fn, w=7, h=4, useDingbats=FALSE) ###

###GGPLOT DRUG COUNT###
p = ggplot(data=clin_amplification_therapeutics_ggplot, aes(x = plot_age, fill = ab_drug))
p = p + facet_wrap( .~acronym, drop=T,scale="fixed", nrow=1)
p = p + geom_bar(position = "stack") 
p = p + labs(fill = "highest predicted\nactionability level", y = "individual cases",
             x = "age at pathologic diagnosis (yrs.)") + theme_bw()
p = p + scale_fill_manual(values=c("None" = "lightslategray", "A on-label" = "dodgerblue1", "A off-label" = "seagreen3",
                                   "B on-label" = "lightgoldenrod", "B off-label" = "lightcoral"),
                          breaks=c("A on-label", "A off-label", "B on-label", "B off-label", "None"))
p = p + theme(legend.position="top")
p
fn = "out/YoungAdult_Amplification_ClinAct_BarplotBinary_Count.pdf"
ggsave(fn, w=7, h=3.5, useDingbats=FALSE)

################################ RESPONSIVE SUMMARY TABLES ################################

ya_amplification_summ_df <- data.frame(matrix(ncol = 7, nrow = 0))
lo_amplification_summ_df <- data.frame(matrix(ncol = 7, nrow = 0))

for (cancer in unique(clin_amplification_therapeutics_ggplot$acronym)) {
  
  temp = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$acronym==cancer,]
  
  ya_sample_size = sum(temp$age_binary==TRUE)
  lo_sample_size = sum(temp$age_binary==FALSE)
  
  ya_AonLab = sum(temp$ab_drug=="A on-label" & temp$age_binary==TRUE)
  ya_AonLab_perc = round((ya_AonLab/ya_sample_size)*100, digits=1)
  
  ya_AoffLab = sum(temp$ab_drug=="A off-label" & temp$age_binary==TRUE)
  ya_AoffLab_perc = round((ya_AoffLab/ya_sample_size)*100, digits=1)
  
  ya_BonLab = sum(temp$ab_drug=="B on-label" & temp$age_binary==TRUE)
  ya_BonLab_perc = round((ya_BonLab/ya_sample_size)*100, digits=1)
  
  ya_BoffLab = sum(temp$ab_drug=="B off-label" & temp$age_binary==TRUE)
  ya_BoffLab_perc = round((ya_BoffLab/ya_sample_size)*100, digits=1)
  
  ya_None = sum(temp$ab_drug=="None" & temp$age_binary==TRUE)
  ya_None_perc = round((ya_None/ya_sample_size)*100, digits=1)
  
  lo_AonLab = sum(temp$ab_drug=="A on-label" & temp$age_binary==FALSE)
  lo_AonLab_perc = round((lo_AonLab/lo_sample_size)*100, digits=1)
  
  lo_AoffLab = sum(temp$ab_drug=="A off-label" & temp$age_binary==FALSE)
  lo_AoffLab_perc = round((lo_AoffLab/lo_sample_size)*100, digits=1)
  
  lo_BonLab = sum(temp$ab_drug=="B on-label" & temp$age_binary==FALSE)
  lo_BonLab_perc = round((lo_BonLab/lo_sample_size)*100, digits=1)
  
  lo_BoffLab = sum(temp$ab_drug=="B off-label" & temp$age_binary==FALSE)
  lo_BoffLab_perc = round((lo_BoffLab/lo_sample_size)*100, digits=1)
  
  lo_None = sum(temp$ab_drug=="None" & temp$age_binary==FALSE)
  lo_None_perc = round((lo_None/lo_sample_size)*100, digits=1)
  
  ya_row = cbind(cancer, ya_sample_size, ya_AonLab_perc, ya_AoffLab_perc, ya_BonLab_perc, ya_BoffLab_perc, ya_None_perc)
  lo_row = cbind(cancer, lo_sample_size, lo_AonLab_perc, lo_AoffLab_perc, lo_BonLab_perc, lo_BoffLab_perc, lo_None_perc)
  
  ya_amplification_summ_df <- rbind(ya_amplification_summ_df, ya_row)
  lo_amplification_summ_df <- rbind(lo_amplification_summ_df, lo_row)
  
}

##############################################################################################################################

ya_amplification_summ_df$cancer = as.character(ya_amplification_summ_df$cancer)
ya_amplification_summ_df = ya_amplification_summ_df[order(ya_amplification_summ_df$cancer),]
rownames(ya_amplification_summ_df) <- ya_amplification_summ_df$cancer
ya_amplification_summ_df$cancer = NULL

colnames(ya_amplification_summ_df) <- c("Young adult cases","% A on-label","% A off-label","% B on-label","% B off-label", "% None")

##############################################################################################################################

lo_amplification_summ_df$cancer = as.character(lo_amplification_summ_df$cancer)
lo_amplification_summ_df = lo_amplification_summ_df[order(lo_amplification_summ_df$cancer),]
rownames(lo_amplification_summ_df) <- lo_amplification_summ_df$cancer
lo_amplification_summ_df$cancer = NULL

colnames(lo_amplification_summ_df) <- c("Later onset cases","% A on-label","% A off-label","% B on-label","% B off-label", "% None")

################################ RESPONSIVE AMPEVNT TABLES ################################

ampEvnt_list = list(); aeCount_list = list()

for (ae in unique(clin_amplification_actionable$ampEvnt)) {
  ampEvnt_list = append(ampEvnt_list, ae)
  aeCount_list = append(aeCount_list, sum(clin_amplification_actionable$ampEvnt==ae))
}

ampEvnt_vec <- data.frame(matrix(unlist(ampEvnt_list),byrow=T),stringsAsFactors=FALSE)
aeCount_vec <- data.frame(matrix(unlist(aeCount_list),byrow=T),stringsAsFactors=FALSE)

temp = cbind(ampEvnt_vec, aeCount_vec)
colnames(temp)[1] = "ampEvnt"
colnames(temp)[2] = "aeCount"
temp = temp[order(temp$aeCount, decreasing=T),]

temp = temp[1:10,]

###############################################

summ_ampEvnt_df <- data.frame(matrix(ncol = 12, nrow = 0))

young_adult_n = c()
later_onset_n = c()

amp_evnt_vec = c("CCND1:AMP", "FGFR1:AMP", "ERBB2:AMP", "EGFR:AMP", "MDM2:AMP",
                 "CCND2:AMP", "KRAS:AMP", "CDK4:AMP", "CDK6:AMP", "MET:AMP")

clin_amplification_therapeutics_ggplot[is.na(clin_amplification_therapeutics_ggplot)] <- ""

for (cancer in unique(clin_amplification_therapeutics_ggplot$acronym)) {
  
  current_cancer_df = clin_amplification_therapeutics_ggplot[clin_amplification_therapeutics_ggplot$acronym==cancer,]
  
  unique_cases = current_cancer_df[!duplicated(current_cancer_df$bcr_patient_barcode),]
  ya_sample_size = sum(unique_cases$age_binary==TRUE)
  lo_sample_size = sum(unique_cases$age_binary==FALSE)
  
  ampEvnt_row <- data.frame(matrix(ncol = 0, nrow = 1))
  
  for (ampEvnt in amp_evnt_vec) {
    
    num_current_ampEvnt_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$ampEvnt==ampEvnt)
    num_current_ampEvnt_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$ampEvnt==ampEvnt)
    
    ya_current_ampEvnt_perc = round((num_current_ampEvnt_yng/ya_sample_size)*100, digits=1)
    lo_current_ampEvnt_perc = round((num_current_ampEvnt_ltr/lo_sample_size)*100, digits=1)
    
    current_ampEvnt_perc = paste(ya_current_ampEvnt_perc, lo_current_ampEvnt_perc, sep=" | ")
    
    ampEvnt_row = cbind(ampEvnt_row, current_ampEvnt_perc)
    
  }
  
  young_adult_n = c(young_adult_n, ya_sample_size)
  later_onset_n = c(later_onset_n, lo_sample_size)
  
  ampEvnt_row = cbind(ampEvnt_row, cancer)
  
  summ_ampEvnt_df <- rbind(summ_ampEvnt_df, ampEvnt_row)
  
}

##############################################################################################################################

colnames(summ_ampEvnt_df) <- c("% CCND1:AMP", "% FGFR1:AMP", "% ERBB2:AMP", "% EGFR:AMP", "% MDM2:AMP",
                               "% CCND2:AMP", "% KRAS:AMP", "% CDK4:AMP", "% CDK6:AMP", "% MET:AMP", "cancer")

young_later_n = paste(young_adult_n, later_onset_n, sep=" | ")

summ_ampEvnt_df = add_column(summ_ampEvnt_df, young_later_n, .before ="% CCND1:AMP")

summ_ampEvnt_df$cancer = as.character(summ_ampEvnt_df$cancer)
summ_ampEvnt_df = summ_ampEvnt_df[order(summ_ampEvnt_df$cancer),]

colnames(summ_ampEvnt_df)[1] = "Young adult cases | Later onset cases"

col_idx <- grep("cancer", names(summ_ampEvnt_df))
summ_ampEvnt_df <- summ_ampEvnt_df[,c(col_idx,(1:ncol(summ_ampEvnt_df))[-col_idx])]