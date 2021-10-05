# set working directory
indir_somatic<-"C:/Users/wangz21/Box Sync/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/"
indir_exp<-"C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/YoungOnsetCancer/analysis/gene_expression/data/gene_expression"
indir_function_glm<-"C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/YoungOnsetCancer/analysis"
indir_clin<-"C:/Users/wangz21/Box Sync/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/"
indir_pc<-"C:/Users/wangz21/Box Sync/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018"
indir_subtype<-"C:/Users/wangz21/Box Sync/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/"

outdir = "C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/YoungOnsetCancer/analysis/gene_expression/out"

# read in statistical functions for glm wrapper
setwd(indir_function_glm)
source("stat_functions.R")

### Gene expression file ###
setwd(indir_exp)
exp<-read.table(unz("cancer_exp.txt.zip","cancer_exp.txt"),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
colnames(exp)<-gsub("\\.","-",colnames(exp))
ma<-read.table("entrezid_symbol_in_expression_profile.txt",header = T,sep = "\t",stringsAsFactors = F)
##### clinical files #####
# clinical file #
setwd(indir_clin)
clin_complete = read.table(header=T, quote = "", sep="\t", fill =T, file = "clinical_PANCAN_patient_with_followup.tsv", stringsAsFactors=FALSE)
clin_complete<-clin_complete[clin_complete$age_at_initial_pathologic_diagnosis!="[Not Available]",]
clin_complete$age_at_initial_pathologic_diagnosis<-as.numeric(clin_complete$age_at_initial_pathologic_diagnosis)
clin_brief = clin_complete[,c("bcr_patient_barcode","acronym","age_at_initial_pathologic_diagnosis","gender")]

# Principal Component (PC) file #
setwd(indir_pc)
PCs<-read.table("2017-04-24_GDAN_AIM_PCA_ethnicity_assigned_WashU.tsv",header=T, quote = "", sep="\t", fill =T, stringsAsFactors=FALSE)[,c("Case","PC1","PC2")]
colnames(PCs)[1] = "bcr_patient_barcode"
PCs = PCs[!duplicated(PCs$bcr_patient_barcode),]
clin_merge = merge(clin_brief,PCs, by= "bcr_patient_barcode",all.x=T,all.y=F)

# subtype file #
setwd(indir_subtype)
subtype = data.frame(readxl::read_xlsx("Subtype Assignments.xlsx"),stringsAsFactors = F)
colnames(subtype)[1] = "bcr_patient_barcode"
clin_merge_subtype = merge(clin_merge,subtype, by= "bcr_patient_barcode")

# create binary age variable; typically <=50 considered as young onset cancer
clin_merge_subtype$age_binary = clin_merge_subtype$age_at_initial_pathologic_diagnosis <= 50

# only include the ones without NA
clin_merge_subtype_avail = clin_merge_subtype
clin_merge_subtype_avail[clin_merge_subtype_avail=="NA"]<-NA
apply(clin_merge_subtype_avail,2,function(x){sum(is.na(x))})
clin_merge_subtype_avail_complete = clin_merge_subtype_avail[complete.cases(clin_merge_subtype_avail) & clin_merge_subtype_avail$gender != "",]

# only include the ones with expression profile
inter_sample<-intersect(colnames(exp),clin_merge_subtype_avail_complete$bcr_patient_barcode)
clin_merge_subtype_avail_complete<-clin_merge_subtype_avail_complete[match(inter_sample,clin_merge_subtype_avail_complete$bcr_patient_barcode),]
exp<-exp[,inter_sample]

setwd(outdir)
write.table(exp,"Gene_exp_with_sample_having_cli_info.txt",row.names = T,col.names = T,quote = F,sep = "\t")
write.table(clin_merge_subtype_avail_complete,"Gene_exp_with_sample_having_cli_info_cli.txt",row.names = F,col.names = T,quote = F,sep = "\t")

##### statistical testing: pan-gene analysis within each cancer type #####
# somatic mutation in driver genes ~ onset age + covariates ("gender","PC1","PC2")
# initiate empty list and index = 1 to store results
results_list = as.list(NULL)
results_list_binary = as.list(NULL)
i=1
# conduct tests by each cancer (denoted as acronym)
gene_name<-rownames(exp)
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  tmpind<-clin_merge_subtype_avail_complete$acronym==cancer
  
  # only analyze cancer with more than 5 patients
  if(sum(tmpind)>5){
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[tmpind,]
    exp_c<-exp[,tmpind]
    exp_c1<-t(exp_c)
    # conduct the test by gene
    for (tmpgene in gene_name){
      tmpexp_c1<-exp_c1[,tmpgene]
      if(sum(is.na(tmpexp_c1))==0){
        Fn<-ecdf(tmpexp_c1)
        tmpexp_c2<-Fn(tmpexp_c1)
        clin_merge_subtype_avail_complete_c_gene<-data.frame(tmpexp_c2,clin_merge_subtype_avail_complete_c,stringsAsFactors = F)
        colnames(clin_merge_subtype_avail_complete_c_gene)[1]<-"Expression"
        
        # model onset age as a linear variable: TODO: add in gender as a covariate although it breaks the code now...
        model_results = run_glm(data = clin_merge_subtype_avail_complete_c_gene, yi = "Expression", xi = "age_at_initial_pathologic_diagnosis", ytype = "Continuous", covi = c("SUBTYPE","gender","PC1","PC2"))
        cancer_stat = data.frame(cbind(cancer, tmpgene, model_results))
        results_list[[i]] = cancer_stat
        
        # model onset age as a binary variable
        model_results = run_glm(data = clin_merge_subtype_avail_complete_c_gene, yi = "Expression", xi = "age_binary", ytype = "Continuous", covi = c("SUBTYPE","gender","PC1","PC2"))
        cancer_stat = data.frame(cbind(cancer, tmpgene, model_results))
        results_list_binary[[i]] = cancer_stat
        
        # increment index
        i = i + 1
      }#for if
    }#for tmpgene
  }#for if
  
  cat(cancer)
  cat(" ")
}#for cancer

##### Compile and store results #####
# compile linear result
tt = do.call(rbind,results_list)
colnames(tt) = c("cancer","Entrez","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F","P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="BH") 
tt=tt[order(tt$FDR, decreasing=FALSE),]
tt$Entrez<-as.integer(as.character(tt$Entrez))
tt$Symbol<-ma[match(tt$Entrez,ma$Entrez),"Symbol"]
setwd(outdir)
write.table(tt, quote=F, sep="\t", file = "Cancer_gene_exp_age_assoc.txt", row.names = F,col.names = T)

# compile binary result
tt = do.call(rbind,results_list_binary)
colnames(tt) = c("cancer","Entrez","y","y_type","x","degrees_freedom","deviance","residual_degrees_freedom","residual_deviance",
                 "F","P_Chi","coefficient","covariates");
# multiple-testing correction and sort by significance
tt$FDR = p.adjust(tt[,"P_Chi"], method="BH") 
tt=tt[order(tt$FDR, decreasing=FALSE),]
tt$Entrez<-as.integer(as.character(tt$Entrez))
tt$Symbol<-ma[match(tt$Entrez,ma$Entrez),"Symbol"]
setwd(outdir)
write.table(tt, quote=F, sep="\t", file = "Cancer_gene_exp_agebinary_assoc.txt", row.names = F)
