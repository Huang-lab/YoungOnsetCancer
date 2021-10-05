##### Somatic_NonHyperMut.R #####
# Will Lee @ January 2021, updated October 2021

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/somatic_mutation"
setwd(bdir)

# read in statistical functions for glm wrapper
source("../stat_functions.R")

### MAIN ###

### MC3 mutation call file (only include the likely driver for the first-pass analysis) ###
somatic_likelyfunctional_driver_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriver.tsv"
somatic_likelyfunctional_driver = read.table(header=T, quote = "", sep="\t", file = somatic_likelyfunctional_driver_f, stringsAsFactors=FALSE)

somatic_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
unique_somatic_samples = unique(substr(somatic$Tumor_Sample_Barcode,1,12))

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
clin_merge_subtype_avail = clin_merge_subtype[clin_merge_subtype$bcr_patient_barcode %in% unique_somatic_samples,]
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

hypermut_path = "~/HypermutatorSamples.txt"
hypermut_df = read.table(header=T, quote = "", sep="\t", file = hypermut_path, stringsAsFactors=FALSE)
hypermut_ids = hypermut_df$TCGA_Barcode # 312 hypermutated samples; 8943-312=8631

clin_merge_subtype_avail_complete$hypermut_status = clin_merge_subtype_avail_complete$bcr_patient_barcode %in% hypermut_ids
clin_merge_subtype_avail_complete=clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$hypermut_status==FALSE,]

# initiate empty list and index = 1 to store results
results_list = as.list(NULL)
results_list_binary = as.list(NULL)
i=1
# conduct tests by each cancer (denoted as acronym)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  if (cancer=='BRCA' | cancer=='CESC' | cancer=='COAD' | cancer=='HNSC' | cancer=='KIRC' | cancer=='KIRP' | cancer=='LGG' |
      cancer=='LIHC' | cancer=='OV' | cancer=='PCPG' | cancer=='SARC' | cancer=='SKCM' | cancer=='THCA' | cancer=='UCEC') {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    # conduct the test by gene
    for (gene in unique(somatic_likelyfunctional_driver$Hugo_Symbol)){
      
      # subset mutation to that gene and get unique mutation
      somatic_likelyfunctional_driver_g = somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$Hugo_Symbol==gene,]
      somatic_likelyfunctional_driver_g_uniq = somatic_likelyfunctional_driver_g[!duplicated(somatic_likelyfunctional_driver_g$bcr_patient_barcode),]
      clin_merge_subtype_avail_complete_c_merge_somatic = merge(clin_merge_subtype_avail_complete_c,somatic_likelyfunctional_driver_g_uniq, by= "bcr_patient_barcode", all.x=T, all.y = F)
      
      # model whether the sample carry somatic mutation in the gene (Hugo Symbol is a standardized gene name)
      # samples with mutations considered as 1, no mutations considered as 0
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol[!is.na(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)] = 1
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol[is.na(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)] = 0
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol = as.numeric(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)
      
      if (sum(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol) < 5){next} # only test for cancer types with at least 5 mutation carriers
      
      # model onset age as a linear variable
      model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_somatic, yi = "Hugo_Symbol", xi = "age_at_initial_pathologic_diagnosis", ytype = "Binary", covi = c("SUBTYPE","gender","PC1","PC2"))
      cancer_stat = data.frame(cbind(cancer, gene, model_results))
      results_list[[i]] = cancer_stat
      # 
      # model onset age as a binary variable
      model_results = run_glm(data = clin_merge_subtype_avail_complete_c_merge_somatic, yi = "Hugo_Symbol", xi = "age_binary", ytype = "Binary", covi = c("SUBTYPE","gender","PC1","PC2"))
      cancer_stat = data.frame(cbind(cancer, gene, model_results))
      results_list_binary[[i]] = cancer_stat
      
      # increment index
      i = i + 1
    }
  }
}

##### Compile and store results #####
# store linear result in hypermut_linear_plot
tt = do.call(rbind,results_list)
colnames(tt) = c("cancer","gene","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="fdr") 
hypermut_linear_plot=tt[order(tt$P_Chi, decreasing=FALSE),]

# store binary result in hypermut_binary_plot
tt = do.call(rbind,results_list_binary)
colnames(tt) = c("cancer","gene","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="fdr") 
hypermut_binary_plot=tt[order(tt$P_Chi, decreasing=FALSE),]

##### plotting #####

library(ggplot2)
library(ggrepel)
source("../global_aes_out.R")

# bubble plot, somatic linear #
p = ggplot(data=hypermut_linear_plot,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.15,gene,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "Coefficient") # + ylim(-2.5,2.6)
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/SomaticHypermutBubble_Linear.pdf'
ggsave(fn,h=5.75, w=5,useDingbat=F)

# bubble plot, somatic binary #
p = ggplot(data=hypermut_binary_plot,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + scale_size(name = expression(-log10(FDR)), breaks = c(1, 2.5, 5.0, 7.5))
p = p + geom_text_repel(aes(label=ifelse(FDR<0.15,gene,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "Coefficient")  + ylim(-2.5,2.6)
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/SomaticHypermutBubble_Binary.pdf'
ggsave(fn,h=6, w=5,useDingbat=F)