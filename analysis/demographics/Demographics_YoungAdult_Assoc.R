##### Demographics_YoungAdult_Assoc.R #####
# William Lee @ June 2020, updated October 2021

library(dplyr)
library(tidyr)

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/Demographics"
setwd(bdir)

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

# only include entries where acronym, age at diagnosis, and gender are not missing values 
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

clin_merge_demo = clin_merge_subtype_avail_complete
clin_merge_demo = clin_merge_demo[complete.cases(clin_merge_demo),] 

clin_merge_demo$age_at_initial_pathologic_diagnosis = as.numeric(clin_merge_demo$age_at_initial_pathologic_diagnosis)

selected_14_cancers <- data.frame(matrix(ncol = 28, nrow = 0))
colnames(selected_14_cancers) = colnames(clin_merge_demo)
for (cancer in unique(clin_merge_demo$acronym)) {
  
  temp = clin_merge_demo[clin_merge_demo$acronym==cancer,]
  
  ya_cases = sum(temp$age_binary==TRUE)
  lo_cases = sum(temp$age_binary==FALSE)
  
  if (ya_cases >= 40 & lo_cases >= 40) {
    selected_14_cancers = rbind(selected_14_cancers,temp)
  }
}

demo_table_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(demo_table_df) <- c("cancer","Sample_size","Female_ratio","Avg_onset_age","young_ct","reglr_ct","Young_adult_ratio")

for (cancer in unique(clin_merge_demo$acronym)) {
  temp = clin_merge_demo[clin_merge_demo$acronym==cancer,]
  
  Sample_size = nrow(temp)
  
  male_ct = sum(temp$gender=="MALE")
  fmle_ct = sum(temp$gender=="FEMALE")
  Female_ratio = round(fmle_ct/(male_ct+fmle_ct), digits=2)
  
  Avg_age = round(mean(temp$age_at_initial_pathologic_diagnosis), digits=1)
  Std_dev = round(sd(temp$age_at_initial_pathologic_diagnosis), digits=1)
  Avg_onset_age = paste(Avg_age, "±", Std_dev)
  
  young_ct = sum(temp$age_binary=="TRUE")
  reglr_ct = sum(temp$age_binary=="FALSE")
  
  Young_adult_ratio = round(young_ct/(reglr_ct+young_ct), digits=2)
  
  demo_table_row = cbind(cancer, Sample_size, Female_ratio, Avg_onset_age, young_ct, reglr_ct, Young_adult_ratio)
  
  demo_table_df <- rbind(demo_table_df, demo_table_row)
}

colnames(demo_table_df) <- c("temp","Sample size","Female ratio","Age at diagnosis (mean ± SD)","Young adult cases (N)","Later-onset cases (N)","Young adult ratio")

demo_table_df$temp = as.character(demo_table_df$temp)
demo_table_df = demo_table_df[order(demo_table_df$temp),]

demo_table_df$`Young adult cases (N)` = as.character(demo_table_df$`Young adult cases (N)`)
demo_table_df$`Young adult cases (N)` = as.numeric(demo_table_df$`Young adult cases (N)`)

demo_table_df$`Later-onset cases (N)` = as.character(demo_table_df$`Later-onset cases (N)`)
demo_table_df$`Later-onset cases (N)` = as.numeric(demo_table_df$`Later-onset cases (N)`)

demo_table_fin = demo_table_df[demo_table_df$`Young adult cases (N)` >= 40  & demo_table_df$`Later-onset cases (N)` >= 40,]

Cancer = c("Breast invasive carcinoma","Cervical squamous cell carcinoma and endocervical adenocarcinoma",
            "Colon adenocarcinoma","Head and Neck squamous cell carcinoma","Kidney renal clear cell carcinoma",
            "Kidney renal papillary cell carcinoma","Brain Lower Grade Glioma","Liver hepatocellular carcinoma",
            "Ovarian serous cystadenocarcinoma","Pheochromocytoma and Paraganglioma","Sarcoma","Skin Cutaneous Melanoma",
            "Thyroid carcinoma","Uterine Corpus Endometrial Carcinoma")

demo_table_fin = add_column(demo_table_fin, Cancer, .before ="Sample size")

############################### SUBTYPE PLOTS ###############################

library(ggplot2)
library(scales)

clin_subtype_df = clin_merge_demo

clin_subtype_df$plot_age[clin_subtype_df$age_binary==TRUE] = "50"
clin_subtype_df$plot_age[clin_subtype_df$age_binary==FALSE] = "> 50"

clin_subtype_female = clin_subtype_df[clin_subtype_df$acronym=="BRCA" | clin_subtype_df$acronym=="CESC" | clin_subtype_df$acronym=="UCEC",]

clin_subtype_allgen = clin_subtype_df[clin_subtype_df$acronym=="COAD" | clin_subtype_df$acronym=="HNSC" |
                                      clin_subtype_df$acronym=="LGG"  | clin_subtype_df$acronym=="SARC",]

##############################################################################################################################

for (cancer in unique(clin_subtype_female$acronym)) {
  for (subt in unique(clin_subtype_female$SUBTYPE)) {
    
    clin_subtype_female$subt_prop[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==TRUE] =
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==TRUE,])/
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$age_binary==TRUE,])
    
    clin_subtype_female$subt_prop[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==FALSE] =
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$SUBTYPE==subt & clin_subtype_female$age_binary==FALSE,])/
      nrow(clin_subtype_female[clin_subtype_female$acronym==cancer & clin_subtype_female$age_binary==FALSE,])
  }
}

clin_subtype_female$subt_perc = round(clin_subtype_female$subt_prop*100, digits=1)

clin_subtype_female$SUBTYPE[clin_subtype_female$SUBTYPE=="CN_HIGH"] = "CN high"
clin_subtype_female$SUBTYPE[clin_subtype_female$SUBTYPE=="CN_LOW"] = "CN low"

p = ggplot(data=clin_subtype_female, aes(x = plot_age, fill = SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scale="free", nrow=3)
p = p + geom_bar(position = "fill")
p = p + scale_fill_brewer(palette = "Set3")
p = p + labs(x = "age at pathologic diagnosis (yrs.)", #y = "Subtype (%)",
             fill = "TCGA Subtype") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + coord_flip()
p = p + ylab(element_blank()) ###
p
fn = "out/Subtype_Female_hBarplotBinary_Percent.pdf"
ggsave(fn, w=5, h=4, useDingbats=FALSE)

##############################################################################################################################

for (cancer in unique(clin_subtype_allgen$acronym)) {
  for (subt in unique(clin_subtype_allgen$SUBTYPE)) {
    
    clin_subtype_allgen$subt_prop[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==TRUE] =
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==TRUE,])/
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$age_binary==TRUE,])
    
    clin_subtype_allgen$subt_prop[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==FALSE] =
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$SUBTYPE==subt & clin_subtype_allgen$age_binary==FALSE,])/
      nrow(clin_subtype_allgen[clin_subtype_allgen$acronym==cancer & clin_subtype_allgen$age_binary==FALSE,])
  }
}

clin_subtype_allgen$subt_perc = round(clin_subtype_allgen$subt_prop*100, digits=1)

manualFillGeneVar = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(clin_subtype_allgen$SUBTYPE))) 
p = ggplot(data=clin_subtype_allgen, aes(x = plot_age, fill = SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scale="free", nrow=4)
p = p + geom_bar(position = "fill")
p = p + scale_fill_manual(values=manualFillGeneVar)
p = p + labs(x = "age at pathologic diagnosis (yrs.)", #y = "Subtype (%)",
             fill = "TCGA Subtype") + theme_bw()
p = p + scale_y_continuous(labels=percent)
p = p + coord_flip()
p = p + ylab(element_blank()) ###
p
fn = "out/Subtype_hBarplotBinary_Percent.pdf"
ggsave(fn, w=5, h=5, useDingbats=FALSE)

######################  FISHER'S EXACT SUBTYPES ######################

clin_subtype_avail_df = clin_subtype_df[clin_subtype_df$acronym=="BRCA" | clin_subtype_df$acronym=="CESC" |
                                        clin_subtype_df$acronym=="COAD" | clin_subtype_df$acronym=="HNSC" |
                                        clin_subtype_df$acronym=="LGG"  | clin_subtype_df$acronym=="SARC" |
                                        clin_subtype_df$acronym=="UCEC",]

fisher_all_subtype_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(fisher_all_subtype_df) <- c("cancer", "subtype", "p_value", "FDR")

for (cancer in unique(clin_subtype_avail_df$acronym)) {

  fisher_subtype_df <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(fisher_subtype_df) <- c("cancer", "subtype", "p_value")
  
  current_cancer_df = clin_subtype_avail_df[clin_subtype_avail_df$acronym==cancer,]
  
  for (subtype in unique(current_cancer_df$SUBTYPE)) {
    
    num_current_subtype_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$SUBTYPE==subtype)
    num_current_subtype_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$SUBTYPE==subtype)
    
    num_not_current_subtype_ltr = sum(current_cancer_df$age_binary==FALSE & current_cancer_df$SUBTYPE!=subtype)
    num_not_current_subtype_yng = sum(current_cancer_df$age_binary==TRUE & current_cancer_df$SUBTYPE!=subtype)
    
    yng_set = c(num_not_current_subtype_yng, num_current_subtype_yng)
    late_set = c(num_not_current_subtype_ltr, num_current_subtype_ltr)
    
    subtype_fisher = cbind(yng_set, late_set)
    
    subtype_fisher_pval = fisher.test(subtype_fisher)$p.val
    
    temp = cbind(cancer, subtype, subtype_fisher_pval)
    colnames(temp) <- c("cancer", "subtype", "p_value")
    
    fisher_subtype_df = rbind(fisher_subtype_df, temp)

  }
  
  fisher_subtype_df$p_value = as.numeric(as.character(fisher_subtype_df$p_value))
  fisher_subtype_df$FDR = p.adjust(fisher_subtype_df[,"p_value"], method="fdr") 
  
  fisher_all_subtype_df = rbind(fisher_all_subtype_df, fisher_subtype_df)
  
}

fisher_all_subtype_df=fisher_all_subtype_df[order(fisher_all_subtype_df$FDR, decreasing=FALSE),]

###################### CNV segments ###################### 

# source the files with plotting related functions
source("../global_aes_out.R")

cnv_seg = "TCGA_mastercalls.abs_segtabs.fixed.CNVperSample.tsv"
CNV_data = read.table(header=T, quote = "", sep="\t", file = cnv_seg, stringsAsFactors=FALSE)
colnames(CNV_data) = c("CNV_segs","SAMPLE_BARCODE")

CNV_seg_merge_df = merge(selected_14_cancers,CNV_data, by="SAMPLE_BARCODE")

CNV_seg_merge_df$plot_age[CNV_seg_merge_df$age_binary==TRUE] = "  50"
CNV_seg_merge_df$plot_age[CNV_seg_merge_df$age_binary==FALSE] = "> 50"

p1 = ggplot(data=CNV_seg_merge_df,aes(x=plot_age,y=log10(CNV_segs),fill=acronym))
p1 = p1 + geom_violin(trim = TRUE)
p1 = p1 + facet_grid( .~acronym, drop=T, space ="free", scale = "free")
p1 = p1 + scale_x_discrete(drop=FALSE) 
p1 = p1 + ylab("log10(copy number segments)") + xlab("age at pathologic diagnosis (yrs.)")
p1 = p1 + theme_bw()
p1 = p1 + getPCACancerFill() 
p1 = p1 + theme(legend.position="none")
p1 = p1 + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.75)
p1
fn = 'out/YoungAdult_CNV_Segment_Violin.pdf'
ggsave(fn,h=3, w=11.5,useDingbats=F)

### wilcox cnv test ###
wilcox_cnv_df <- data.frame(matrix(ncol = 5, nrow = 0))
for (cancer in unique(CNV_seg_merge_df$acronym)) {
  
  temp = CNV_seg_merge_df[CNV_seg_merge_df$acronym==cancer,]
  
  temp_young = temp[temp$age_binary==TRUE,]
  temp_later = temp[temp$age_binary==FALSE,]
  
  Avg_young_segs = round(mean(temp_young$CNV_segs), digits=1)
  Std_young_segs = round(sd(temp_young$CNV_segs), digits=1)
  young_adult_segments = paste(Avg_young_segs, "±", Std_young_segs)
  
  Avg_later_segs = round(mean(temp_later$CNV_segs), digits=1)
  Std_later_segs = round(sd(temp_later$CNV_segs), digits=1)
  later_onset_segments = paste(Avg_later_segs, "±", Std_later_segs)
  
  wilcox = wilcox.test(temp_young$CNV_segs, temp_later$CNV_segs, paired=FALSE)
  
  cur_row = cbind(cancer, nrow(temp_young), nrow(temp_later), young_adult_segments, later_onset_segments, wilcox$statistic, wilcox$p.value)
  
  wilcox_cnv_df = rbind(wilcox_cnv_df, cur_row)
  
}

colnames(wilcox_cnv_df) <- c("Cancer", "Young adult cases (N)", "Later-onset cases (N)", "Young adult CN segments (mean ± SD)", "Later-onset CN segments (mean ± SD)", "U statistic", "p-value")
wilcox_cnv_df$`p-value` = as.numeric(wilcox_cnv_df$`p-value`)
wilcox_cnv_df$FDR = p.adjust(wilcox_cnv_df[,"p-value"], method="fdr")
wilcox_cnv_df = wilcox_cnv_df[order(wilcox_cnv_df$`p-value`, decreasing=FALSE),]

###################### Aneuploidy score ###################### 

aneuploidy_df = select(immune_feat_df, bcr_patient_barcode, Aneuploidy_Score)
aneuploidy_df = aneuploidy_df[complete.cases(aneuploidy_df),]

aneuploidy_merge_df = merge(selected_14_cancers, aneuploidy_df, by= "bcr_patient_barcode")

### wilcox aneuploidy test ###
wilcox_aneuploidy_df <- data.frame(matrix(ncol = 5, nrow = 0))
for (cancer in unique(aneuploidy_merge_df$acronym)) {
  
  temp = aneuploidy_merge_df[aneuploidy_merge_df$acronym==cancer,]
  
  temp_young = temp[temp$age_binary==TRUE,]
  temp_later = temp[temp$age_binary==FALSE,]
  
  Avg_young_aneuploidy = round(mean(temp_young$Aneuploidy_Score), digits=1)
  Std_young_aneuploidy = round(sd(temp_young$Aneuploidy_Score), digits=1)
  young_adult_aneuploidy = paste(Avg_young_aneuploidy, "±", Std_young_aneuploidy)
  
  Avg_later_aneuploidy = round(mean(temp_later$Aneuploidy_Score), digits=1)
  Std_later_aneuploidy = round(sd(temp_later$Aneuploidy_Score), digits=1)
  later_onset_aneuploidy = paste(Avg_later_aneuploidy, "±", Std_later_aneuploidy)
  
  wilcox = wilcox.test(temp_young$Aneuploidy_Score, temp_later$Aneuploidy_Score, paired=FALSE)
  
  cur_row = cbind(cancer, nrow(temp_young), nrow(temp_later), young_adult_aneuploidy, later_onset_aneuploidy, wilcox$statistic, wilcox$p.value)
  
  wilcox_aneuploidy_df = rbind(wilcox_aneuploidy_df, cur_row)
  
}

colnames(wilcox_aneuploidy_df) <- c("Cancer", "Young adult cases (N)", "Later-onset cases (N)", "Young adult aneuploidy score (mean ± SD)", "Later-onset aneuploidy score (mean ± SD)", "U statistic", "p-value")
wilcox_aneuploidy_df$`p-value` = as.numeric(wilcox_aneuploidy_df$`p-value`)
wilcox_aneuploidy_df$FDR = p.adjust(wilcox_aneuploidy_df[,"p-value"], method="fdr")
wilcox_aneuploidy_df = wilcox_aneuploidy_df[order(wilcox_aneuploidy_df$`p-value`, decreasing=FALSE),]

###################### Leukocyte fraction ###################### 

leukocyte_df = select(immune_feat_df, bcr_patient_barcode, Leukocyte_Fraction)
leukocyte_df = leukocyte_df[complete.cases(leukocyte_df),]

leukocyte_merge_df = merge(selected_14_cancers, leukocyte_df, by= "bcr_patient_barcode")

### wilcox leukocyte test ###
wilcox_leukocyte_df <- data.frame(matrix(ncol = 5, nrow = 0))
for (cancer in unique(leukocyte_merge_df$acronym)) {
  
  temp = leukocyte_merge_df[leukocyte_merge_df$acronym==cancer,]
  
  temp_young = temp[temp$age_binary==TRUE,]
  temp_later = temp[temp$age_binary==FALSE,]
  
  Avg_young_leukocyte = round(mean(temp_young$Leukocyte_Fraction), digits=2)
  Std_young_leukocyte = round(sd(temp_young$Leukocyte_Fraction), digits=2)
  young_adult_leukocyte = paste(Avg_young_leukocyte, "±", Std_young_leukocyte)
  
  Avg_later_leukocyte = round(mean(temp_later$Leukocyte_Fraction), digits=2)
  Std_later_leukocyte = round(sd(temp_later$Leukocyte_Fraction), digits=2)
  later_onset_leukocyte = paste(Avg_later_leukocyte, "±", Std_later_leukocyte)
  
  wilcox = wilcox.test(temp_young$Leukocyte_Fraction, temp_later$Leukocyte_Fraction, paired=FALSE)
  
  cur_row = cbind(cancer, nrow(temp_young), nrow(temp_later), young_adult_leukocyte, later_onset_leukocyte, wilcox$statistic, wilcox$p.value)
  
  wilcox_leukocyte_df = rbind(wilcox_leukocyte_df, cur_row)
  
}

colnames(wilcox_leukocyte_df) <- c("Cancer", "Young adult cases (N)", "Later-onset cases (N)", "Young adult leukocyte fraction (mean ± SD)", "Later-onset leukocyte fraction (mean ± SD)", "U statistic", "p-value")
wilcox_leukocyte_df$`p-value` = as.numeric(wilcox_leukocyte_df$`p-value`)
wilcox_leukocyte_df$FDR = p.adjust(wilcox_leukocyte_df[,"p-value"], method="fdr")
wilcox_leukocyte_df = wilcox_leukocyte_df[order(wilcox_leukocyte_df$`p-value`, decreasing=FALSE),]