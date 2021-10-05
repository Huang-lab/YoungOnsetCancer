##### ICGC_YoungAdultCancer.R #####
# Will Lee @ January 2021, updated October 2021

library(readxl) 
library(stringr) 
library(dplyr) 

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/icgc_validation"
setwd(bdir)

sample_path = "pcawg_sample_sheet.tsv"
sample_df = read.table(header=T, quote = "", sep="\t", file = sample_path, stringsAsFactors=FALSE)
colnames(sample_df)[6] = "sample_id"
sample_df = sample_df %>% select(donor_unique_id, sample_id) # 7255 obs.

drivers_path = "TableS3_panorama_driver_mutations_ICGC_samples.public.tsv.gz" 
drivers_df = read.table(header=T, quote = "", sep="\t", file = drivers_path, stringsAsFactors=FALSE)
drivers_df = drivers_df %>% select(sample_id, ttype, gene, top_category)
drivers_df = merge(drivers_df, sample_df, by="sample_id")

drivers_df = drivers_df[!str_detect(drivers_df$gene, 'lncrna'),]
drivers_df = drivers_df[!str_detect(drivers_df$gene, 'enhancers'),]
drivers_df = drivers_df[!str_detect(drivers_df$gene, 'telomere'),]

somatic_drivers = drivers_df[drivers_df$top_category=='mutational',] # 3661 obs.
cnv_drivers = drivers_df[drivers_df$top_category=='CNA',] # 3695 obs.

histology_path = "pcawg_specimen_histology_August2016_v9.xlsx"
histology_df = read_excel(histology_path)
colnames(histology_df)[1] = "icgc_specimen_id"
colnames(histology_df)[13] = 'SUBTYPE'
histology_df = histology_df %>% select(donor_unique_id, SUBTYPE)
histology_df = histology_df[!duplicated(histology_df$donor_unique_id),] # 2834 obs.

clinical_path = "pcawg_donor_clinical_August2016_v9.xlsx"
clinical_df = read_excel(clinical_path) 
colnames(clinical_df)[1] = "donor_unique_id"
clinical_df = clinical_df[is.na(clinical_df$tcga_donor_uuid),] # remove 885 TCGA cases (2834 -> 1949 obs.)
temp = strsplit(clinical_df$project_code, '-') ##
temp_df = as.data.frame(do.call(rbind,temp)) ##
clinical_df$acronym = temp_df$V1
clinical_df$country = temp_df$V2
colnames(clinical_df)[6] = 'gender'
colnames(clinical_df)[11] = 'age_at_initial_pathologic_diagnosis'
clinical_df$age_binary = clinical_df$age_at_initial_pathologic_diagnosis <= 50
clinical_df = clinical_df %>% select(donor_unique_id, acronym, country, gender, age_at_initial_pathologic_diagnosis, age_binary) # 1949 obs.

clin_hist_merge = merge(clinical_df, histology_df, by="donor_unique_id", all.x=T, all.y=F) # 1949 obs.

icgc_complete_df = clin_hist_merge[complete.cases(clin_hist_merge),] # 1820 obs.
icgc_complete_df = icgc_complete_df[!duplicated(icgc_complete_df$donor_unique_id),] # 1820 obs.

icgc_complete_df$acronym[icgc_complete_df$acronym=="LICA"] = "LIHC"
icgc_complete_df$acronym[icgc_complete_df$acronym=="LINC"] = "LIHC"
icgc_complete_df$acronym[icgc_complete_df$acronym=="LIRI"] = "LIHC"

icgc_complete_df$acronym[icgc_complete_df$acronym=="MELA"] = "SKCM"

###############################################################################################

source("../stat_functions.R")

# initiate empty list and index = 1 to store results
results_list = as.list(NULL)
results_list_binary = as.list(NULL)
i=1
# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(icgc_complete_df$acronym)){
  
  # manually subset analysis to cancers of interest
  if (cancer=='BRCA' | cancer=='LIHC' | cancer=='OV' | cancer=='RECA' | cancer=='SKCM') {
  
    icgc_complete_df_c = icgc_complete_df[icgc_complete_df$acronym==cancer,]
    icgc_complete_df_c = icgc_complete_df_c[complete.cases(icgc_complete_df_c),]
    
    for (gene in unique(somatic_drivers$gene)){ # run for somatic driver validation
    # for (gene in unique(cnv_drivers$gene)){ # run for CNV driver validation
      
      driver_g = somatic_drivers[somatic_drivers$gene==gene,] # run for somatic driver validation
      # driver_g = cnv_drivers[cnv_drivers$gene==gene,] # run for CNV driver validation
      
      driver_g_uniq = driver_g[!duplicated(driver_g$donor_unique_id),]
      icgc_complete_df_merge_driver = merge(icgc_complete_df_c,driver_g_uniq, by="donor_unique_id", all.x=T, all.y=F)
      
      icgc_complete_df_merge_driver$gene[!is.na(icgc_complete_df_merge_driver$gene)] = 1
      icgc_complete_df_merge_driver$gene[is.na(icgc_complete_df_merge_driver$gene)] = 0
      icgc_complete_df_merge_driver$gene = as.numeric(icgc_complete_df_merge_driver$gene)
      
      if (sum(icgc_complete_df_merge_driver$gene) < 3) {next} # only test for cancer types with at least 3 mutation carriers
      
      # model onset age as a binary variable
      model_results = run_glm(data = icgc_complete_df_merge_driver, yi = "gene", xi = "age_binary", ytype = "Binary", covi = c("SUBTYPE","gender","country"))
      cancer_stat = data.frame(cbind(cancer, gene, model_results))
      results_list_binary[[i]] = cancer_stat
      
      # increment index
      i = i + 1
    }
  }
}

# store binary result
tt = do.call(rbind,results_list_binary)
colnames(tt) = c("cancer","gene","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="fdr") 
icgc_binary_plot=tt[order(tt$P_Chi, decreasing=FALSE),]
tn = "out/somatic_icgc_onsetAgeBinary_assoc.txt"
# tn = "out/cnv_icgc_onsetAgeBinary_assoc.txt"
write.table(icgc_binary_plot, quote=F, sep="\t", file = tn, row.names = F)

###########################################################################################################################################################

library(ggplot2)
library(ggrepel)
source("../global_aes_out.R")

###########################################################################################################################################################

somatic_log_adj = read.table(header=T, quote = "", sep="\t", file='out/somatic_icgc_onsetAgeBinary_assoc.txt', stringsAsFactors=FALSE)

somatic_log_adj$neg_vals = somatic_log_adj$coefficient < 0
somatic_log_adj$coefficient = abs(somatic_log_adj$coefficient)
somatic_log_adj$coefficient = log10(somatic_log_adj$coefficient)

somatic_log_adj$coefficient[somatic_log_adj$coefficient > 0 & somatic_log_adj$neg_vals==TRUE] = 
somatic_log_adj$coefficient[somatic_log_adj$coefficient > 0 & somatic_log_adj$neg_vals==TRUE]*(-1)

somatic_log_adj$coefficient[somatic_log_adj$coefficient < 0 & somatic_log_adj$neg_vals==FALSE] = 
somatic_log_adj$coefficient[somatic_log_adj$coefficient < 0 & somatic_log_adj$neg_vals==FALSE]*(-1)

# bubble plot, somatic binary #
p = ggplot(data=somatic_log_adj,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + scale_size(name = expression(-log10(FDR)), breaks = c(0.1,0.25,0.50,0.75))
p = p + geom_text_repel(aes(label=ifelse(P_Chi<0.05,gene,NA)), size=2.5, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "log10(Coefficient)")  # + ylim()
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/ICGC_Somatic_Bubble_Binary.pdf'
ggsave(fn,h=5, w=4,useDingbat=F)

###########################################################################################################################################################

cnv_log_adj = read.table(header=T, quote = "", sep="\t", file='out/cnv_icgc_onsetAgeBinary_assoc.txt', stringsAsFactors=FALSE)

cnv_log_adj$neg_vals = cnv_log_adj$coefficient < 0
cnv_log_adj$coefficient = abs(cnv_log_adj$coefficient)
cnv_log_adj$coefficient = log10(cnv_log_adj$coefficient)

cnv_log_adj$coefficient[cnv_log_adj$coefficient > 0 & cnv_log_adj$neg_vals==TRUE] = 
cnv_log_adj$coefficient[cnv_log_adj$coefficient > 0 & cnv_log_adj$neg_vals==TRUE]*(-1)

cnv_log_adj$coefficient[cnv_log_adj$coefficient < 0 & cnv_log_adj$neg_vals==FALSE] = 
cnv_log_adj$coefficient[cnv_log_adj$coefficient < 0 & cnv_log_adj$neg_vals==FALSE]*(-1)

# bubble plot, cnv binary #
p = ggplot(data=cnv_log_adj,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + scale_size(name = expression(-log10(FDR)), breaks = c(0.1,0.25,0.50,0.75,1.0))
p = p + geom_text_repel(aes(label=ifelse(P_Chi<0.05,gene,NA)), size=2.5, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "log10(Coefficient)")  # + ylim()
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/ICGC_CNV_Bubble_Binary.pdf'
ggsave(fn,h=5, w=4,useDingbat=F)