##### plotEpisil_hmap.R #####
# William Lee @ May 2020, updated October 2021

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

episil_df <- cbind(oncosig_df[,1], select(oncosig_df, starts_with("EPISIL"))) # 17
names(episil_df)[1]='bcr_patient_barcode'

## clinical files
clin_complete_f = "~/Box/Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv" #for analysis locally
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

# cleaning of merged df
clin_merge_subtype_avail_complete$cancer = NULL
clin_merge_subtype_avail_complete$ethnicity = NULL
clin_merge_subtype_avail_complete$washu_assigned_ethnicity = NULL
clin_merge_subtype_avail_complete$Sample = NULL
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "NA" ] <- NA
clin_merge_subtype_avail_complete[ clin_merge_subtype_avail_complete == "[Not Available]"] <- NA

clin_episil_sig_sub <- data.frame(matrix(ncol = 28, nrow = 0))
colnames(clin_episil_sig_sub) <- colnames(clin_merge_subtype_avail_complete)

clin_episil_ggplot <- data.frame(matrix(ncol = 32, nrow = 0))
colnames(clin_episil_ggplot) <- colnames(clin_merge_subtype_avail_complete)
colnames(clin_episil_ggplot)[29:32] <- c("status", "plot_age", "event", "event_stat")

for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  if (cancer == "BRCA" | cancer == "LGG") {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    clin_episil_sig_sub = rbind(clin_episil_sig_sub, clin_merge_subtype_avail_complete_c)
    
    if (cancer == "BRCA") {sig_event_vec = c("EPISIL.TCF7","EPISIL.CDKN2A","EPISIL.MGA")}
    
    if (cancer == "LGG") {sig_event_vec = c("EPISIL.LATS2","EPISIL.TCF7","EPISIL.HES4","EPISIL.NOV")}
    
    # conduct the test by each alteration event
    for (event in sig_event_vec){
      
      mf_df_current = select(episil_df, bcr_patient_barcode, event)
      mf_df_current = mf_df_current[mf_df_current[,event] != 'NA',]
      if (!is.logical(mf_df_current[,event])) {
        mf_df_current[,mf] = as.logical(mf_df_current[,event])
      }
      
      clin_merge_subtype_avail_complete_c_merge_mf = merge(clin_merge_subtype_avail_complete_c, mf_df_current, by= "bcr_patient_barcode", all.x=T, all.y = F)
      clin_merge_subtype_avail_complete_c_merge_mf = clin_merge_subtype_avail_complete_c_merge_mf[clin_merge_subtype_avail_complete_c_merge_mf[,event] != 'NA',]
      
      clin_merge_subtype_avail_complete_c_merge_mf[,29][(clin_merge_subtype_avail_complete_c_merge_mf[,29])==TRUE] = "+"
      clin_merge_subtype_avail_complete_c_merge_mf[,29][(clin_merge_subtype_avail_complete_c_merge_mf[,29])==FALSE]= "-"
      
      colnames(clin_merge_subtype_avail_complete_c_merge_mf)[29] = "status"
      
      clin_merge_subtype_avail_complete_c_merge_mf$age_at_initial_pathologic_diagnosis = as.numeric(clin_merge_subtype_avail_complete_c_merge_mf$age_at_initial_pathologic_diagnosis)
      
      clin_merge_subtype_avail_complete_c_merge_mf$plot_age[clin_merge_subtype_avail_complete_c_merge_mf$age_binary==TRUE] = "<= 50"
      clin_merge_subtype_avail_complete_c_merge_mf$plot_age[clin_merge_subtype_avail_complete_c_merge_mf$age_binary==FALSE] = "> 50"
      
      clin_merge_subtype_avail_complete_c_merge_mf$event = event
      
      clin_merge_subtype_avail_complete_c_merge_mf$event_stat = paste(clin_merge_subtype_avail_complete_c_merge_mf$event,
                                                                      clin_merge_subtype_avail_complete_c_merge_mf$status,sep=": ")
      
      clin_episil_ggplot <- rbind(clin_episil_ggplot, clin_merge_subtype_avail_complete_c_merge_mf)
    }  
  }
}

clin_episil_BRCA = clin_episil_sig_sub[clin_episil_sig_sub$acronym=="BRCA",]
clin_episil_LGG = clin_episil_sig_sub[clin_episil_sig_sub$acronym=="LGG",]

clin_episil_ggplot = clin_episil_ggplot[!is.na(clin_episil_ggplot$bcr_patient_barcode),]

clin_episil_ggplot_BRCA = clin_episil_ggplot[clin_episil_ggplot$acronym=="BRCA",]
clin_episil_ggplot_LGG = clin_episil_ggplot[clin_episil_ggplot$acronym=="LGG",]


################################# BRCA HEATMAP #################################

BRCA_hmap = clin_episil_ggplot_BRCA[clin_episil_ggplot_BRCA$status=="+",]
total_yng = nrow(clin_episil_BRCA[clin_episil_BRCA$age_binary==TRUE,])
total_late = nrow(clin_episil_BRCA[clin_episil_BRCA$age_binary==FALSE,])

BRCA_hmap$hmap_label[BRCA_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
BRCA_hmap$hmap_label[BRCA_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (evnt in unique(BRCA_hmap$event_stat)) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    BRCA_hmap$hmap_data[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE] =
      nrow(BRCA_hmap[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])/total_yng
    
    BRCA_hmap$hmap_data[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE] =
      nrow(BRCA_hmap[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])/total_late
    
  }
}

BRCA_hmap$subt_perc = round(BRCA_hmap$hmap_data*100, digits=1)

p = ggplot(data=BRCA_hmap, mapping = aes(x = SUBTYPE, y = reorder(event_stat,desc(event_stat)), fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65))
p = p + guides(fill=guide_colorbar(title="% BRCA"))
#p = p + theme(legend.position = "top")
p
fn = 'out/BRCA_episil_subtype_heatmap.pdf'
ggsave(fn,h=1.75, w=5.75,useDingbats=F)

################################# LGG HEATMAP #################################

LGG_hmap = clin_episil_ggplot_LGG[clin_episil_ggplot_LGG$status=="+",]
total_yng = nrow(clin_episil_LGG[clin_episil_LGG$age_binary==TRUE,])
total_late = nrow(clin_episil_LGG[clin_episil_LGG$age_binary==FALSE,])

LGG_hmap$hmap_label[LGG_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
LGG_hmap$hmap_label[LGG_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (evnt in unique(LGG_hmap$event_stat)) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    LGG_hmap$hmap_data[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE] =
      nrow(LGG_hmap[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])/total_yng
    
    LGG_hmap$hmap_data[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE] =
      nrow(LGG_hmap[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])/total_late
    
  }
}

LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-non-codel"] = "IDHmut-non\n-codel"
LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-codel"] = "IDHmut\n-codel"

LGG_hmap$subt_perc = round(LGG_hmap$hmap_data*100, digits=1)

p = ggplot(data=LGG_hmap, mapping = aes(x = SUBTYPE, y = reorder(event_stat,desc(event_stat)), fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=4.25)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.y = element_text(size=10))
p = p + guides(fill=guide_colorbar(title="% LGG"))
#p = p + theme(legend.position = "top")
p
fn = 'out/LGG_episil_subtype_heatmap.pdf'
ggsave(fn,h=3, w=6.25,useDingbats=F)


################################# BRCA COUNT HEATMAP #################################

BRCA_hmap = clin_episil_ggplot_BRCA[clin_episil_ggplot_BRCA$status=="+",]
total_yng = nrow(clin_episil_BRCA[clin_episil_BRCA$age_binary==TRUE,])
total_late = nrow(clin_episil_BRCA[clin_episil_BRCA$age_binary==FALSE,])

BRCA_hmap$hmap_label[BRCA_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
BRCA_hmap$hmap_label[BRCA_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (evnt in unique(BRCA_hmap$event_stat)) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    BRCA_hmap$hmap_data[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE] =
      nrow(BRCA_hmap[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])
    
    BRCA_hmap$hmap_data[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE] =
      nrow(BRCA_hmap[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])
    
  }
}

p = ggplot(data=BRCA_hmap, mapping = aes(x = SUBTYPE, y = reorder(event_stat,desc(event_stat)), fill = hmap_data))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=hmap_data),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65))
p = p + guides(fill=guide_colorbar(title="# BRCA"))
#p = p + theme(legend.position = "top")
p
fn = 'out/BRCA_episil_subtype_count_heatmap.pdf'
ggsave(fn,h=1.75, w=5.75,useDingbats=F)

################################# LGG COUNT HEATMAP #################################

LGG_hmap = clin_episil_ggplot_LGG[clin_episil_ggplot_LGG$status=="+",]
total_yng = nrow(clin_episil_LGG[clin_episil_LGG$age_binary==TRUE,])
total_late = nrow(clin_episil_LGG[clin_episil_LGG$age_binary==FALSE,])

LGG_hmap$hmap_label[LGG_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
LGG_hmap$hmap_label[LGG_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (evnt in unique(LGG_hmap$event_stat)) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    LGG_hmap$hmap_data[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE] =
      nrow(LGG_hmap[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])
    
    LGG_hmap$hmap_data[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE] =
      nrow(LGG_hmap[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])
    
  }
}

LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-non-codel"] = "IDHmut-non\n-codel"
LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-codel"] = "IDHmut\n-codel"

p = ggplot(data=LGG_hmap, mapping = aes(x = SUBTYPE, y = reorder(event_stat,desc(event_stat)), fill = hmap_data))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=hmap_data),color="orangered",size=4.25)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.y = element_text(size=10))
p = p + guides(fill=guide_colorbar(title="# LGG"))
#p = p + theme(legend.position = "top")
p
fn = 'out/LGG_episil_subtype_count_heatmap.pdf'
ggsave(fn,h=3, w=6.25,useDingbats=F)

#
#
#
#
#

################################# FISHER'S EXACT BRCA ANALYSIS #################################

# BRCA event-subtype analysis #
BRCA_hmap = clin_episil_ggplot_BRCA[clin_episil_ggplot_BRCA$status=="+",]
total_yng = nrow(clin_episil_BRCA[clin_episil_BRCA$age_binary==TRUE,])
total_late = nrow(clin_episil_BRCA[clin_episil_BRCA$age_binary==FALSE,])

BRCA_fisher_episil_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(BRCA_fisher_episil_df) <- c("event", "subtype", "p_value")

for (evnt in unique(BRCA_hmap$event_stat)) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    fisher_yng_subevnt = nrow(BRCA_hmap[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])
    fisher_yng_n_subevnt = total_yng - fisher_yng_subevnt
    
    fisher_late_subevnt = nrow(BRCA_hmap[BRCA_hmap$event_stat==evnt & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])
    fisher_late_n_subevnt = total_late - fisher_late_subevnt
    
    yng_set = c(fisher_yng_n_subevnt, fisher_yng_subevnt)
    late_set = c(fisher_late_n_subevnt, fisher_late_subevnt)
    
    BRCA_fisher = cbind(yng_set, late_set)
    
    BRCA_fisher_pval = fisher.test(BRCA_fisher)$p.val
    
    temp = cbind(evnt, subt, BRCA_fisher_pval)
    colnames(temp) <- c("event", "subtype", "p_value")
    
    BRCA_fisher_episil_df = rbind(BRCA_fisher_episil_df, temp)
    
  }
}

BRCA_fisher_episil_df$p_value = as.numeric(as.character(BRCA_fisher_episil_df$p_value))

BRCA_fisher_episil_df$FDR = p.adjust(BRCA_fisher_episil_df[,"p_value"], method="fdr") 
BRCA_fisher_episil_df=BRCA_fisher_episil_df[order(BRCA_fisher_episil_df$p_value, decreasing=FALSE),]

################################# FISHER'S EXACT LGG ANALYSIS #################################

# LGG event-subtype analysis #
LGG_hmap = clin_episil_ggplot_LGG[clin_episil_ggplot_LGG$status=="+",]
total_yng = nrow(clin_episil_LGG[clin_episil_LGG$age_binary==TRUE,])
total_late = nrow(clin_episil_LGG[clin_episil_LGG$age_binary==FALSE,])

LGG_fisher_episil_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(LGG_fisher_episil_df) <- c("event", "subtype", "p_value")

for (evnt in unique(LGG_hmap$event_stat)) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    fisher_yng_subevnt = nrow(LGG_hmap[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])
    fisher_yng_n_subevnt = total_yng - fisher_yng_subevnt
    
    fisher_late_subevnt = nrow(LGG_hmap[LGG_hmap$event_stat==evnt & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])
    fisher_late_n_subevnt = total_late - fisher_late_subevnt
    
    yng_set = c(fisher_yng_n_subevnt, fisher_yng_subevnt)
    late_set = c(fisher_late_n_subevnt, fisher_late_subevnt)
    
    LGG_fisher = cbind(yng_set, late_set)
    
    LGG_fisher_pval = fisher.test(LGG_fisher)$p.val
    
    temp = cbind(evnt, subt, LGG_fisher_pval)
    colnames(temp) <- c("event", "subtype", "p_value")
    
    LGG_fisher_episil_df = rbind(LGG_fisher_episil_df, temp)
    
  }
}

LGG_fisher_episil_df$p_value = as.numeric(as.character(LGG_fisher_episil_df$p_value))

LGG_fisher_episil_df$FDR = p.adjust(LGG_fisher_episil_df[,"p_value"], method="fdr") 
LGG_fisher_episil_df=LGG_fisher_episil_df[order(LGG_fisher_episil_df$p_value, decreasing=FALSE),]