##### plotImmuneAssoc_Continuous.R #####
# William Lee @ September 2019, updated October 2021

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/immune_features"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

#################### immune heatmaps #################### 

assoc_IGS_n = "out/IGS_onsetAge_assoc.txt"
assoc_IGS = read.table(sep="\t",header=T,file=assoc_IGS_n, stringsAsFactors=FALSE)

assoc_IGS$ID = "Immune Gene Signature"

assoc_IGS$immune_feature[assoc_IGS$immune_feature=="TGF_beta_Response"] = "TGF-b Response"
assoc_IGS$immune_feature[assoc_IGS$immune_feature=="IFN_gamma_Response"] = "IFN-y Response"
assoc_IGS$immune_feature[assoc_IGS$immune_feature=="Lymphocyte_Infiltration_Signature_Score"] = "Lymphocyte Infiltration"

assoc_IGS$immune_feature <- gsub('_', ' ', assoc_IGS$immune_feature)

p = ggplot(data=assoc_IGS, mapping = aes(x = cancer, y = reorder(immune_feature,desc(immune_feature)), fill = coefficient))
p = p + facet_wrap( .~ID, drop=T)#,scales="free_y",nrow=2)
p = p + geom_tile()
p = p + geom_text(aes(label=ifelse(FDR<0.05,"***",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.05 & FDR<0.10,"**",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.10 & FDR<0.15,"*",NA)))
p = p + theme_bw()
p = p + scale_fill_gradient2(low = "dodgerblue4",mid = "white",high = "red",midpoint = 0,guide = "colourbar",aesthetics = "fill")
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(element_blank())
p = p + guides(fill=guide_colorbar(title="coefficient",barheight=5.25))
p = p + theme(axis.text.x = element_blank())
p = p + theme(axis.ticks.x = element_blank())
#p = p + theme(legend.position = "top")
p
fn = 'out/IGS_heatmap_continuous.pdf'
ggsave(fn,h=1.75, w=6.75,useDingbats=F)

#
#
#
#
#

assoc_IIC_n = "out/IIC_onsetAge_assoc.txt"
assoc_IIC = read.table(sep="\t",header=T,file=assoc_IIC_n, stringsAsFactors=FALSE)

assoc_IIC$ID = "Immune Infiltrate Composition"

assoc_IIC$immune_feature[assoc_IIC$immune_feature=="Macrophages_M1"] = "M1 Macrophages"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="Macrophages_M2"] = "M2 Macrophages"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="NK_Cells_Resting"] = "NK Cells, Resting"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="NK_Cells_Activated"] = "NK Cells, Activated"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="B_Cells_Naive"] = "B Cells, Naive"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="T_Cells_CD8"] = "CD8+ T Cells"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="T_Cells_CD4_Naive"] = "CD4+ T Cells, Naive"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="T_Cells_Regulatory_Tregs"] = "Regulatory T Cells"

assoc_IIC$immune_feature <- gsub('_', ' ', assoc_IIC$immune_feature)

p = ggplot(data=assoc_IIC, mapping = aes(x = cancer, y = reorder(immune_feature,desc(immune_feature)), fill = coefficient))
p = p + facet_wrap( .~ID, drop=T)#,scales="free_y",nrow=2)
p = p + geom_tile()
p = p + geom_text(aes(label=ifelse(FDR<0.05,"***",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.05 & FDR<0.10,"**",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.10 & FDR<0.15,"*",NA)))
p = p + theme_bw()
p = p + scale_fill_gradient2(low = "dodgerblue4",mid = "white",high = "red",midpoint = 0,guide = "colourbar",aesthetics = "fill")
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(label = "cancer") + ylab(element_blank())
p = p + guides(fill=guide_colorbar(title="coefficient",barheight=7))
p = p + theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5))
#p = p + theme(legend.position = "top")
p
fn = 'out/IIC_heatmap_continuous.pdf'
ggsave(fn,h=3.25, w=7,useDingbats=F)

#
#
#
#
#

assoc_Th_Cells_n = "out/Th_Cells_onsetAge_assoc.txt"
assoc_Th_Cells = read.table(sep="\t",header=T,file=assoc_Th_Cells_n, stringsAsFactors=FALSE)

assoc_Th_Cells$ID = "T-helper Cells"

assoc_Th_Cells$immune_feature[assoc_Th_Cells$immune_feature=="Th1_Cells"] = "Th1"
assoc_Th_Cells$immune_feature[assoc_Th_Cells$immune_feature=="Th2_Cells"] = "Th2"
assoc_Th_Cells$immune_feature[assoc_Th_Cells$immune_feature=="Th17_Cells"] = "Th17"

assoc_Th_Cells$immune_feature <- factor(assoc_Th_Cells$immune_feature, levels = c("Th1","Th2","Th17"))

p = ggplot(data=assoc_Th_Cells, mapping = aes(x = immune_feature, y = reorder(cancer,desc(cancer)), fill = coefficient))
p = p + facet_wrap( .~ID, drop=T)#,scales="free_y",nrow=2)
p = p + geom_tile()
p = p + geom_text(aes(label=ifelse(FDR<0.05,"***",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.05 & FDR<0.10,"**",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.10 & FDR<0.15,"*",NA)))
p = p + theme_bw()
p = p + scale_fill_gradient2(low = "dodgerblue4",mid = "white",high = "red",midpoint = 0,guide = "colourbar",aesthetics = "fill")
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(label = "cancer")
p = p + guides(fill=guide_colorbar(title="coefficient",barheight=7))
p = p + theme(axis.text.y = element_text(colour="black", size=10, vjust=0.5))
#p = p + theme(legend.position = "top")
p
fn = 'out/ThCells_heatmap_continuous.pdf'
ggsave(fn,h=4, w=3.25,useDingbats=F)

#
#
#
#
#

assoc_neo_n = "out/Neoantigens_onsetAge_assoc.txt"
assoc_neo = read.table(sep="\t",header=T,file=assoc_neo_n, stringsAsFactors=FALSE)

assoc_neo$ID = "Neoantigens"

assoc_neo$immune_feature[assoc_neo$immune_feature=="SNV_Neoantigens"] = "SNV"
assoc_neo$immune_feature[assoc_neo$immune_feature=="Indel_Neoantigens"] = "Indel"

p = ggplot(data=assoc_neo, mapping = aes(x = immune_feature, y = reorder(cancer,desc(cancer)), fill = coefficient))
p = p + facet_wrap( .~ID, drop=T)#,scales="free_y",nrow=2)
p = p + geom_tile()
p = p + geom_text(aes(label=ifelse(FDR<0.05,"***",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.05 & FDR<0.10,"**",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.10 & FDR<0.15,"*",NA)))
p = p + theme_bw()
p = p + scale_fill_gradient2(low = "dodgerblue4",mid = "white",high = "red",midpoint = 0,guide = "colourbar",aesthetics = "fill")
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(label = "cancer")
p = p + guides(fill=guide_colorbar(title="coefficient",barheight=7))
p = p + theme(axis.text.y = element_text(colour="black", size=10, vjust=0.5))
#p = p + theme(legend.position = "top")
p
fn = 'out/Neoantigens_heatmap_continuous.pdf'
ggsave(fn,h=4, w=2.75,useDingbats=F)

#
#
#
#
#

assoc_mut_n = "out/Mut_Rate_onsetAge_assoc.txt"
assoc_mut = read.table(sep="\t",header=T,file=assoc_mut_n, stringsAsFactors=FALSE)

assoc_mut$ID = "Mutation Rate"

assoc_mut$immune_feature[assoc_mut$immune_feature=="Nonsilent_Mutation_Rate"] = "Nonsil."
assoc_mut$immune_feature[assoc_mut$immune_feature=="Silent_Mutation_Rate"] = "Silent"

p = ggplot(data=assoc_mut, mapping = aes(x = immune_feature, y = reorder(cancer,desc(cancer)), fill = coefficient))
p = p + facet_wrap( .~ID, drop=T)#,scales="free_y",nrow=2)
p = p + geom_tile()
p = p + geom_text(aes(label=ifelse(FDR<0.05,"***",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.05 & FDR<0.10,"**",NA)))
p = p + geom_text(aes(label=ifelse(FDR>=0.10 & FDR<0.15,"*",NA)))
p = p + theme_bw()
p = p + scale_fill_gradient2(low = "dodgerblue4",mid = "white",high = "red",midpoint = 0,guide = "colourbar",aesthetics = "fill")
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(label = "cancer")
p = p + guides(fill=guide_colorbar(title="coefficient",barheight=7))
p = p + theme(axis.text.y = element_text(colour="black", size=10, vjust=0.5))
#p = p + theme(legend.position = "top")
p
fn = 'out/Mutation_heatmap_continuous.pdf'
ggsave(fn,h=4, w=2.75,useDingbats=F)