##### Immune_Mut_Scatter.R #####
# William Lee @ May 2020, updated October 2021

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

nonsil_mut_df <- data.frame(matrix(ncol = 29, nrow = 0))
colnames(nonsil_mut_df) <- colnames(clin_merge_subtype_avail_complete)
colnames(nonsil_mut_df)[29] <- "Nonsilent_Mutation_Rate"

silent_mut_df <- data.frame(matrix(ncol = 29, nrow = 0))
colnames(silent_mut_df) <- colnames(clin_merge_subtype_avail_complete)
colnames(silent_mut_df)[29] <- "Silent_Mutation_Rate"

leuk_frac_df <- data.frame(matrix(ncol = 29, nrow = 0))
colnames(leuk_frac_df) <- colnames(clin_merge_subtype_avail_complete)
colnames(leuk_frac_df)[29] <- "Leukocyte_Fraction"

aneuploidy_df <- data.frame(matrix(ncol = 29, nrow = 0))
colnames(aneuploidy_df) <- colnames(clin_merge_subtype_avail_complete)
colnames(aneuploidy_df)[29] <- "Aneuploidy_Score"

for (cancer in unique(clin_merge_subtype_avail_complete$acronym)) {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    if ((sum(clin_merge_subtype_avail_complete_c$age_binary==TRUE) >= 40) &
       (sum(clin_merge_subtype_avail_complete_c$age_binary==FALSE) >= 40)) {
    
    # conduct the test by immune feature
    for (immune_feat in immune_features) {
      
      if (immune_feat == "Nonsilent_Mutation_Rate") {
      
      # if (immune_feat == "Silent_Mutation_Rate") {
      
      # if (immune_feat == "Leukocyte_Fraction") {
        
      # if (immune_feat == "Aneuploidy_Score") {
      
      # if (immune_feat == "Homologous_Recombination_Defects") {
        
        immune_feat_df_current = select(immune_feat_df, bcr_patient_barcode, immune_feat)
        immune_feat_df_current = na.omit(immune_feat_df_current)
        
        clin_merge_subtype_avail_complete_c_merge_if = merge(clin_merge_subtype_avail_complete_c,immune_feat_df_current, by= "bcr_patient_barcode", all.x = T, all.y = F)
        clin_merge_subtype_avail_complete_c_merge_if = clin_merge_subtype_avail_complete_c_merge_if[complete.cases(clin_merge_subtype_avail_complete_c_merge_if[29]),]
        
        nonsil_mut_df = rbind(nonsil_mut_df, clin_merge_subtype_avail_complete_c_merge_if)
        # silent_mut_df = rbind(silent_mut_df, clin_merge_subtype_avail_complete_c_merge_if)
        # leuk_frac_df = rbind(leuk_frac_df, clin_merge_subtype_avail_complete_c_merge_if)
        # aneuploidy_df = rbind(aneuploidy_df, clin_merge_subtype_avail_complete_c_merge_if)
        # hr_defects_df = rbind(hr_defects_df, clin_merge_subtype_avail_complete_c_merge_if)
       
      }
    }
  }
}

nonsil_mut_noSubt = nonsil_mut_df[nonsil_mut_df$SUBTYPE=="Not_Applicable",]

silent_mut_noSubt = silent_mut_df[silent_mut_df$SUBTYPE=="Not_Applicable",]

nonsil_mut_allgen = nonsil_mut_df[nonsil_mut_df$acronym=="COAD" | nonsil_mut_df$acronym=="HNSC" |
                                  nonsil_mut_df$acronym=="LGG"  | nonsil_mut_df$acronym=="SARC",]

silent_mut_allgen = silent_mut_df[silent_mut_df$acronym=="COAD" | silent_mut_df$acronym=="HNSC" |
                                  silent_mut_df$acronym=="LGG"  | silent_mut_df$acronym=="SARC",]

nonsil_mut_female = nonsil_mut_df[nonsil_mut_df$acronym=="BRCA" | nonsil_mut_df$acronym=="CESC" | nonsil_mut_df$acronym=="UCEC",]

silent_mut_female = silent_mut_df[silent_mut_df$acronym=="BRCA" | silent_mut_df$acronym=="CESC" | silent_mut_df$acronym=="UCEC",]

# nonsilent scatter, no subtype #
p = ggplot(data=nonsil_mut_noSubt,aes(x=log10(Nonsilent_Mutation_Rate), y=as.numeric(age_at_initial_pathologic_diagnosis),color=acronym))
p = p + facet_wrap( .~acronym, drop=T,scales="free_y",nrow=7, strip.position = "right")
p = p + geom_point(alpha=0.8,show.legend = FALSE)
p = p + getPCACancerColor() #+ xlim(-2.5,2.5)
p = p + scale_y_continuous(breaks=c(50,75))
p = p + labs(x="log10(Nonsilent Mutation Rate)", y="age at initial pathologic diagnosis (yrs.)")
p = p + geom_vline(xintercept = 0, alpha=0.2)
p = p + geom_hline(yintercept = 50, alpha=0.8)
p = p + theme_bw() +
theme(axis.title.x = element_text(size=8), axis.title.y = element_text(size=8)) +
theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))
p
fn = 'out/Nonsilent_Scatter_NoSub.pdf'
ggsave(fn,h=6.5, w=2.5,useDingbat=F)

# nonsilent scatter, all gender subtype #
manualFillSubtype = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(nonsil_mut_allgen$SUBTYPE))) 
p = ggplot(data=nonsil_mut_allgen,aes(x=log10(Nonsilent_Mutation_Rate), y=as.numeric(age_at_initial_pathologic_diagnosis),color=SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scales="free_y",nrow=4, strip.position = "right")
p = p + geom_point(alpha=0.8, show.legend = FALSE)
p = p + scale_color_manual(values=manualFillSubtype)
p = p + scale_y_continuous(breaks=c(50,75))
p = p + labs(x="log10(Nonsilent Mutation Rate)", y="age at initial pathologic diagnosis (yrs.)", color="subtype")
p = p + geom_vline(xintercept = 0, alpha=0.2)
p = p + geom_hline(yintercept = 50, alpha=0.8)
p = p  + theme_bw() 
p
fn = 'out/Nonsilent_Scatter_allGen.pdf'
ggsave(fn,h=5.25, w=3.75,useDingbat=F)

# nonsilent scatter, female subtype #
p = ggplot(data=nonsil_mut_female,aes(x=log10(Nonsilent_Mutation_Rate), y=as.numeric(age_at_initial_pathologic_diagnosis),color=SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scales="free_y",nrow=3, strip.position = "right")
p = p + geom_point(alpha=0.8, show.legend = FALSE)
p = p + scale_color_brewer(palette = "Set3")
p = p + scale_y_continuous(breaks=c(50,75))
p = p + labs(x="log10(Nonsilent Mutation Rate)", y="age at initial pathologic diagnosis (yrs.)", color="subtype")
p = p + geom_vline(xintercept = 0, alpha=0.2)
p = p + geom_hline(yintercept = 50, alpha=0.8)
p = p  + theme_bw() 
p
fn = 'out/Nonsilent_Scatter_Female.pdf'
ggsave(fn,h=4.25, w=3.75,useDingbat=F)

####################################################################################################################

# silent scatter, no subtype #
p = ggplot(data=silent_mut_noSubt,aes(x=log10(Silent_Mutation_Rate), y=as.numeric(age_at_initial_pathologic_diagnosis),color=acronym))
p = p + facet_wrap( .~acronym, drop=T,scales="free_y",nrow=7, strip.position = "right")
p = p + geom_point(alpha=0.8,show.legend = FALSE)
p = p + getPCACancerColor() #+ xlim(-2.5,2.5)
p = p + scale_y_continuous(breaks=c(50,75))
p = p + labs(x="log10(Silent Mutation Rate)",y="age at initial pathologic diagnosis (yrs.)")
p = p + geom_vline(xintercept = 0, alpha=0.2)
p = p + geom_hline(yintercept = 50, alpha=0.8)
p = p  + theme_bw() +
theme(axis.title.x = element_text(size=8), axis.title.y = element_text(size=8)) +
theme(axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))
p
fn = 'out/Silent_Scatter_NoSub.pdf'
ggsave(fn,h=6.5, w=2.5,useDingbat=F)

# silent scatter, all gender subtype #
manualFillSubtype = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(silent_mut_allgen$SUBTYPE))) 
p = ggplot(data=silent_mut_allgen,aes(x=log10(Silent_Mutation_Rate), y=as.numeric(age_at_initial_pathologic_diagnosis),color=SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scales="free_y",nrow=4, strip.position = "right")
p = p + geom_point(alpha=0.8, show.legend = FALSE)
p = p + scale_color_manual(values=manualFillSubtype)
p = p + scale_y_continuous(breaks=c(50,75))
p = p + labs(x="log10(Silent Mutation Rate)", y="age at initial pathologic diagnosis (yrs.)" ,color="subtype")
p = p + geom_vline(xintercept = 0, alpha=0.2)
p = p + geom_hline(yintercept = 50, alpha=0.8)
p = p  + theme_bw() 
p
fn = 'out/Silent_Scatter_allGen.pdf'
ggsave(fn,h=5.25, w=3.75,useDingbat=F)

# silent scatter, female subtype #
p = ggplot(data=silent_mut_female,aes(x=log10(Silent_Mutation_Rate), y=as.numeric(age_at_initial_pathologic_diagnosis),color=SUBTYPE))
p = p + facet_wrap( .~acronym, drop=T,scales="free_y",nrow=3, strip.position = "right")
p = p + geom_point(alpha=0.8, show.legend = FALSE)
p = p + scale_color_brewer(palette = "Set3")
p = p + scale_y_continuous(breaks=c(50,75))
p = p + labs(x="log10(Silent Mutation Rate)", y="age at initial pathologic diagnosis (yrs.)", color="subtype")
p = p + geom_vline(xintercept = 0, alpha=0.2)
p = p + geom_hline(yintercept = 50, alpha=0.8)
p = p  + theme_bw() 
p
fn = 'out/Silent_Scatter_Female.pdf'
ggsave(fn,h=4.25, w=3.75,useDingbat=F)


################################ Leukocyte Fraction & Aneuploidy Score ################################

library(plyr)
library(dplyr)
library(reshape2)
library(UpSetR)
library(ggrepel)

source("../global_aes_out.R")

# Leuk Frac
leuk_frac_df$plot_age[leuk_frac_df$age_binary==TRUE] = "  50"
leuk_frac_df$plot_age[leuk_frac_df$age_binary==FALSE] = "> 50"
#
p1 = ggplot(data=leuk_frac_df,aes(x=plot_age,y=Leukocyte_Fraction,fill=acronym))
p1 = p1 + geom_violin(trim = TRUE)
p1 = p1 + facet_grid( .~acronym, drop=T, space ="free", scale = "free")
p1 = p1 + scale_x_discrete(drop=FALSE) #+ scale_y_discrete(drop=FALSE)
p1 = p1 + ylab("Leukocyte Fraction") + xlab("age at pathologic diagnosis (yrs.)")
p1 = p1 + theme_bw()
p1 = p1 + getPCACancerFill()
p1 = p1 + theme(legend.position="none")
p1 = p1 + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.75)
p1
fn = 'out/YoungAdult_LeukFrac_Violin.pdf'
ggsave(fn,h=2.75, w=11,useDingbats=F)

# Aneuploidy
aneuploidy_df$plot_age[aneuploidy_df$age_binary==TRUE] = "  50"
aneuploidy_df$plot_age[aneuploidy_df$age_binary==FALSE] = "> 50"
#
p1 = ggplot(data=aneuploidy_df,aes(x=plot_age,y=Aneuploidy_Score,fill=acronym))
p1 = p1 + geom_violin(trim = TRUE)
p1 = p1 + facet_grid( .~acronym, drop=T, space ="free", scale = "free")
p1 = p1 + scale_x_discrete(drop=FALSE) #+ scale_y_discrete(drop=FALSE)
p1 = p1 + ylab("Aneuploidy Score") + xlab("age at pathologic diagnosis (yrs.)")
p1 = p1 + theme_bw()
p1 = p1 + getPCACancerFill()
p1 = p1 + theme(legend.position="none")
p1 = p1 + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.75)
p1
fn = 'out/YoungAdult_Aneuploidy_Violin.pdf'
ggsave(fn,h=2.75, w=11,useDingbats=F)