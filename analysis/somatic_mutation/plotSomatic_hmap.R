##### plotSomatic_hmap.R #####
# William Lee @ May 2020, updated October 2021

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/somatic_mutation"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

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

clin_somatic_sig_sub <- data.frame(matrix(ncol = 28, nrow = 0))
colnames(clin_somatic_sig_sub) <- colnames(clin_merge_subtype_avail_complete)

clin_somatic_mut_ggplot <- data.frame(matrix(ncol = 31, nrow = 0))
colnames(clin_somatic_mut_ggplot) <- colnames(clin_merge_subtype_avail_complete)
colnames(clin_somatic_mut_ggplot)[29:31] <-c("plot_age", "gene", "gene_m_stat")

# for somatic analysis, we are just interested in BRCA, UCEC, & LGG
for (cancer in unique(clin_merge_subtype_avail_complete$acronym)){
  
  if (cancer == "BRCA" | cancer == "UCEC" | cancer == "LGG") {
    
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete[clin_merge_subtype_avail_complete$acronym==cancer,]
    clin_merge_subtype_avail_complete_c = clin_merge_subtype_avail_complete_c[complete.cases(clin_merge_subtype_avail_complete_c),]
    
    clin_somatic_sig_sub = rbind(clin_somatic_sig_sub, clin_merge_subtype_avail_complete_c)
    
    # conduct the test by gene
    
    if (cancer == "BRCA") {sig_gene_vec = c("CDH1","GATA3","KMT2C")}
    
    if (cancer == "LGG") {sig_gene_vec = c("ATRX","EGFR","IDH1","TP53")}
    
    if (cancer == "UCEC") {sig_gene_vec = c("ATRX","BRD7","CNBD1","CTNNB1","FLT3","LATS1","PTEN","RPS6KA3","SIN3A")}
    
    for (gene in sig_gene_vec) {
      
      # subset mutation to that gene and get unique mutation
      somatic_likelyfunctional_driver_g = somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$Hugo_Symbol==gene,]
      somatic_likelyfunctional_driver_g_uniq = somatic_likelyfunctional_driver_g[!duplicated(somatic_likelyfunctional_driver_g$bcr_patient_barcode),]
      clin_merge_subtype_avail_complete_c_merge_somatic = merge(clin_merge_subtype_avail_complete_c,somatic_likelyfunctional_driver_g_uniq, by= "bcr_patient_barcode", all.x=T, all.y = F)
      
      # check whether the sample carry somatic mutation in the gene (Hugo Symbol is a standardized gene name)
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol[!is.na(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)] = "MUT"
      clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol[is.na(clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol)] = "WT"
      
      clin_merge_subtype_avail_complete_c_merge_somatic$age_at_initial_pathologic_diagnosis = as.numeric(clin_merge_subtype_avail_complete_c_merge_somatic$age_at_initial_pathologic_diagnosis)
      
      clin_merge_subtype_avail_complete_c_merge_somatic$plot_age[clin_merge_subtype_avail_complete_c_merge_somatic$age_binary==TRUE] = "<= 50"
      clin_merge_subtype_avail_complete_c_merge_somatic$plot_age[clin_merge_subtype_avail_complete_c_merge_somatic$age_binary==FALSE] = "> 50"
      
      clin_merge_subtype_avail_complete_c_merge_somatic$gene = gene
      
      clin_merge_subtype_avail_complete_c_merge_somatic$gene_m_stat = paste(clin_merge_subtype_avail_complete_c_merge_somatic$gene,
                                                                            clin_merge_subtype_avail_complete_c_merge_somatic$Hugo_Symbol,sep=":\n")
      
      clin_somatic_mut_ggplot <- rbind(clin_somatic_mut_ggplot, clin_merge_subtype_avail_complete_c_merge_somatic)
    }
  }
}

clin_somatic_UCEC = clin_somatic_sig_sub[clin_somatic_sig_sub$acronym=="UCEC",]
clin_somatic_BRCA = clin_somatic_sig_sub[clin_somatic_sig_sub$acronym=="BRCA",]
clin_somatic_LGG = clin_somatic_sig_sub[clin_somatic_sig_sub$acronym=="LGG",]

clin_somatic_mut_ggplot_UCEC = clin_somatic_mut_ggplot[clin_somatic_mut_ggplot$acronym=="UCEC",]
clin_somatic_mut_ggplot_BRCA = clin_somatic_mut_ggplot[clin_somatic_mut_ggplot$acronym=="BRCA",]
clin_somatic_mut_ggplot_LGG = clin_somatic_mut_ggplot[clin_somatic_mut_ggplot$acronym=="LGG",]

#
#
#
#
#

################################# BRCA HEATMAP #################################

BRCA_hmap = clin_somatic_mut_ggplot_BRCA[clin_somatic_mut_ggplot_BRCA$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_BRCA[clin_somatic_BRCA$age_binary==TRUE,])
total_late = nrow(clin_somatic_BRCA[clin_somatic_BRCA$age_binary==FALSE,])

BRCA_hmap$hmap_label[BRCA_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
BRCA_hmap$hmap_label[BRCA_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (gmut in unique(BRCA_hmap$gene_m_stat)) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])/total_yng
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])/total_late
    
  }
}

BRCA_hmap$subt_perc = round(BRCA_hmap$hmap_data*100, digits=1)

# BRCA heatmap #
p = ggplot(data=BRCA_hmap, mapping = aes(x = SUBTYPE, y = reorder(gene_m_stat,desc(gene_m_stat)), fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="% BRCA"))
#p = p + theme(legend.position = "top")
p
fn = 'out/BRCA_mut_subtype_heatmap.pdf'
ggsave(fn,h=1.75, w=5.25,useDingbats=F)

################################# UCEC HEATMAP #################################

UCEC_hmap = clin_somatic_mut_ggplot_UCEC[clin_somatic_mut_ggplot_UCEC$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_UCEC[clin_somatic_UCEC$age_binary==TRUE,])
total_late = nrow(clin_somatic_UCEC[clin_somatic_UCEC$age_binary==FALSE,])

UCEC_hmap$hmap_label[UCEC_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
UCEC_hmap$hmap_label[UCEC_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (gmut in unique(UCEC_hmap$gene_m_stat)) {
  for (subt in unique(UCEC_hmap$SUBTYPE)) {
    
    UCEC_hmap$hmap_data[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE] =
      nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE,])/total_yng
    
    UCEC_hmap$hmap_data[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE] =
      nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE,])/total_late
    
  }
}

UCEC_hmap$subt_perc = round(UCEC_hmap$hmap_data*100, digits=1)

# UCEC heatmap #
p = ggplot(data=UCEC_hmap, mapping = aes(x = SUBTYPE, y = reorder(gene_m_stat,desc(gene_m_stat)), fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65,size=8), axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="% UCEC"))
#p = p + theme(legend.position = "top")
p
fn = 'out/UCEC_mut_subtype_heatmap.pdf'
ggsave(fn,h=3.25, w=5.25,useDingbats=F)

################################# LGG HEATMAP #################################

LGG_hmap = clin_somatic_mut_ggplot_LGG[clin_somatic_mut_ggplot_LGG$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_LGG[clin_somatic_LGG$age_binary==TRUE,])
total_late = nrow(clin_somatic_LGG[clin_somatic_LGG$age_binary==FALSE,])

LGG_hmap$hmap_label[LGG_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
LGG_hmap$hmap_label[LGG_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (gmut in unique(LGG_hmap$gene_m_stat)) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    LGG_hmap$hmap_data[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE] =
      nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])/total_yng
    
    LGG_hmap$hmap_data[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE] =
      nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])/total_late
    
  }
}

LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-non-codel"] = "IDHmut-non\n-codel"
LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-codel"] = "IDHmut\n-codel"

LGG_hmap$subt_perc = round(LGG_hmap$hmap_data*100, digits=1)

# LGG heatmap #
p = ggplot(data=LGG_hmap, mapping = aes(x = SUBTYPE, y = reorder(gene_m_stat,desc(gene_m_stat)), fill = subt_perc))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=subt_perc),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(size=8))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="% LGG"))
#p = p + theme(legend.position = "top")
p
fn = 'out/LGG_mut_subtype_heatmap.pdf'
ggsave(fn,h=2.25, w=5.25,useDingbats=F)

#
#
#
#
#

################################# BRCA COUNT HEATMAP #################################

BRCA_hmap = clin_somatic_mut_ggplot_BRCA[clin_somatic_mut_ggplot_BRCA$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_BRCA[clin_somatic_BRCA$age_binary==TRUE,])
total_late = nrow(clin_somatic_BRCA[clin_somatic_BRCA$age_binary==FALSE,])

BRCA_hmap$hmap_label[BRCA_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
BRCA_hmap$hmap_label[BRCA_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (gmut in unique(BRCA_hmap$gene_m_stat)) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])
    
    BRCA_hmap$hmap_data[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE] =
      nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])
    
  }
}

# BRCA count heatmap #
p = ggplot(data=BRCA_hmap, mapping = aes(x = SUBTYPE, y = reorder(gene_m_stat,desc(gene_m_stat)), fill = hmap_data))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=hmap_data),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="# BRCA"))
#p = p + theme(legend.position = "top")
p
fn = 'out/BRCA_mut_subtype_count_heatmap.pdf'
ggsave(fn,h=1.75, w=5.25,useDingbats=F)

################################# UCEC COUNT HEATMAP #################################

UCEC_hmap = clin_somatic_mut_ggplot_UCEC[clin_somatic_mut_ggplot_UCEC$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_UCEC[clin_somatic_UCEC$age_binary==TRUE,])
total_late = nrow(clin_somatic_UCEC[clin_somatic_UCEC$age_binary==FALSE,])

UCEC_hmap$hmap_label[UCEC_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
UCEC_hmap$hmap_label[UCEC_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (gmut in unique(UCEC_hmap$gene_m_stat)) {
  for (subt in unique(UCEC_hmap$SUBTYPE)) {
    
    UCEC_hmap$hmap_data[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE] =
      nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE,])
    
    UCEC_hmap$hmap_data[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE] =
      nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE,])
    
  }
}

# UCEC count heatmap #
p = ggplot(data=UCEC_hmap, mapping = aes(x = SUBTYPE, y = reorder(gene_m_stat,desc(gene_m_stat)), fill = hmap_data))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=hmap_data),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(angle=22.5,vjust=0.65,size=8), axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="# UCEC"))
#p = p + theme(legend.position = "top")
p
fn = 'out/UCEC_mut_subtype_count_heatmap.pdf'
ggsave(fn,h=3.25, w=5.25,useDingbats=F)

################################# LGG COUNT HEATMAP #################################

LGG_hmap = clin_somatic_mut_ggplot_LGG[clin_somatic_mut_ggplot_LGG$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_LGG[clin_somatic_LGG$age_binary==TRUE,])
total_late = nrow(clin_somatic_LGG[clin_somatic_LGG$age_binary==FALSE,])

LGG_hmap$hmap_label[LGG_hmap$age_binary==TRUE] = paste("young adult, n = ", total_yng, sep="")
LGG_hmap$hmap_label[LGG_hmap$age_binary==FALSE] = paste("later onset, n = ", total_late, sep="")

for (gmut in unique(LGG_hmap$gene_m_stat)) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    LGG_hmap$hmap_data[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE] =
      nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])
    
    LGG_hmap$hmap_data[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE] =
      nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])
    
  }
}

LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-non-codel"] = "IDHmut-non\n-codel"
LGG_hmap$SUBTYPE[LGG_hmap$SUBTYPE=="IDHmut-codel"] = "IDHmut\n-codel"

# LGG count heatmap #
p = ggplot(data=LGG_hmap, mapping = aes(x = SUBTYPE, y = reorder(gene_m_stat,desc(gene_m_stat)), fill = hmap_data))
p = p + facet_wrap( .~hmap_label, drop=T,scale="fixed", nrow=1)
p = p + geom_tile()
p = p + geom_text(aes(label=hmap_data),color="orangered",size=3.5)
p = p + theme_bw()
p = p + scale_fill_distiller(direction = 1)
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + theme(axis.text.x = element_text(size=8))#, axis.text.y = element_text(size=8))
p = p + guides(fill=guide_colorbar(title="# LGG"))
#p = p + theme(legend.position = "top")
p
fn = 'out/LGG_mut_subtype_count_heatmap.pdf'
ggsave(fn,h=2.25, w=5.25,useDingbats=F)

#
#
#
#
#

################################# FISHER'S EXACT BRCA ANALYSIS #################################

# BRCA mutation-subtype analysis #
BRCA_hmap = clin_somatic_mut_ggplot_BRCA[clin_somatic_mut_ggplot_BRCA$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_BRCA[clin_somatic_BRCA$age_binary==TRUE,])
total_late = nrow(clin_somatic_BRCA[clin_somatic_BRCA$age_binary==FALSE,])

BRCA_fisher_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(BRCA_fisher_df) <- c("gene:mutation", "subtype", "p_value")

for (gmut in unique(BRCA_hmap$gene_m_stat)) {
  for (subt in unique(BRCA_hmap$SUBTYPE)) {
    
    fisher_yng_submut = nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==TRUE,])
    fisher_yng_n_submut = total_yng - fisher_yng_submut
    
    fisher_late_submut = nrow(BRCA_hmap[BRCA_hmap$gene_m_stat==gmut & BRCA_hmap$SUBTYPE==subt & BRCA_hmap$age_binary==FALSE,])
    fisher_late_n_submut = total_late - fisher_late_submut
    
    yng_set = c(fisher_yng_n_submut, fisher_yng_submut)
    late_set = c(fisher_late_n_submut, fisher_late_submut)
    
    BRCA_fisher = cbind(yng_set, late_set)
    
    BRCA_fisher_pval = fisher.test(BRCA_fisher)$p.val
    
    temp = cbind(gmut, subt, BRCA_fisher_pval)
    colnames(temp) <- c("gene:mutation", "subtype", "p_value")

    BRCA_fisher_df = rbind(BRCA_fisher_df, temp)
    
  }
}

BRCA_fisher_df$p_value = as.numeric(as.character(BRCA_fisher_df$p_value))

BRCA_fisher_df$FDR = p.adjust(BRCA_fisher_df[,"p_value"], method="fdr") 
BRCA_fisher_df=BRCA_fisher_df[order(BRCA_fisher_df$p_value, decreasing=FALSE),]

################################# FISHER'S EXACT UCEC ANALYSIS #################################

# UCEC mutation-subtype analysis #
UCEC_hmap = clin_somatic_mut_ggplot_UCEC[clin_somatic_mut_ggplot_UCEC$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_UCEC[clin_somatic_UCEC$age_binary==TRUE,])
total_late = nrow(clin_somatic_UCEC[clin_somatic_UCEC$age_binary==FALSE,])

UCEC_fisher_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(UCEC_fisher_df) <- c("gene:mutation", "subtype", "p_value")

for (gmut in unique(UCEC_hmap$gene_m_stat)) {
  for (subt in unique(UCEC_hmap$SUBTYPE)) {
    
    fisher_yng_submut = nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==TRUE,])
    fisher_yng_n_submut = total_yng - fisher_yng_submut
    
    fisher_late_submut = nrow(UCEC_hmap[UCEC_hmap$gene_m_stat==gmut & UCEC_hmap$SUBTYPE==subt & UCEC_hmap$age_binary==FALSE,])
    fisher_late_n_submut = total_late - fisher_late_submut
    
    yng_set = c(fisher_yng_n_submut, fisher_yng_submut)
    late_set = c(fisher_late_n_submut, fisher_late_submut)
    
    UCEC_fisher = cbind(yng_set, late_set)
    
    UCEC_fisher_pval = fisher.test(UCEC_fisher)$p.val
    
    temp = cbind(gmut, subt, UCEC_fisher_pval)
    colnames(temp) <- c("gene:mutation", "subtype", "p_value")
    
    UCEC_fisher_df = rbind(UCEC_fisher_df, temp)
    
  }
}

UCEC_fisher_df$p_value = as.numeric(as.character(UCEC_fisher_df$p_value))

UCEC_fisher_df$FDR = p.adjust(UCEC_fisher_df[,"p_value"], method="fdr") 
UCEC_fisher_df=UCEC_fisher_df[order(UCEC_fisher_df$p_value, decreasing=FALSE),]

################################# FISHER'S EXACT LGG ANALYSIS #################################

# LGG mutation-subtype analysis #
LGG_hmap = clin_somatic_mut_ggplot_LGG[clin_somatic_mut_ggplot_LGG$Hugo_Symbol=="MUT",]
total_yng = nrow(clin_somatic_LGG[clin_somatic_LGG$age_binary==TRUE,])
total_late = nrow(clin_somatic_LGG[clin_somatic_LGG$age_binary==FALSE,])

LGG_fisher_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(LGG_fisher_df) <- c("gene:mutation", "subtype", "p_value")

for (gmut in unique(LGG_hmap$gene_m_stat)) {
  for (subt in unique(LGG_hmap$SUBTYPE)) {
    
    fisher_yng_submut = nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==TRUE,])
    fisher_yng_n_submut = total_yng - fisher_yng_submut
    
    fisher_late_submut = nrow(LGG_hmap[LGG_hmap$gene_m_stat==gmut & LGG_hmap$SUBTYPE==subt & LGG_hmap$age_binary==FALSE,])
    fisher_late_n_submut = total_late - fisher_late_submut
    
    yng_set = c(fisher_yng_n_submut, fisher_yng_submut)
    late_set = c(fisher_late_n_submut, fisher_late_submut)
    
    LGG_fisher = cbind(yng_set, late_set)
    
    LGG_fisher_pval = fisher.test(LGG_fisher)$p.val
    
    temp = cbind(gmut, subt, LGG_fisher_pval)
    colnames(temp) <- c("gene:mutation", "subtype", "p_value")
    
    LGG_fisher_df = rbind(LGG_fisher_df, temp)
    
  }
}

LGG_fisher_df$p_value = as.numeric(as.character(LGG_fisher_df$p_value))

LGG_fisher_df$FDR = p.adjust(LGG_fisher_df[,"p_value"], method="fdr") 
LGG_fisher_df=LGG_fisher_df[order(LGG_fisher_df$p_value, decreasing=FALSE),]