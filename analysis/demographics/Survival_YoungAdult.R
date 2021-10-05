### Survival_YoungAdult.R ###
# Will Lee @ October 2020, updated October 2021

library(ggplot2)
library(ggpubr)
library(ggpmisc)

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/immune_features"
setwd(bdir)

#############################################################################

assoc_IGS_n = "out/IGS_onsetAgeBinary_assoc.txt"
assoc_IGS = read.table(sep="\t",header=T,file=assoc_IGS_n, stringsAsFactors=FALSE)

assoc_IGS$immune_feature[assoc_IGS$immune_feature=="Lymphocyte_Infiltration_Signature_Score"] = "Lymphocyte Infiltration"

lymphocyte_infil = assoc_IGS[assoc_IGS$immune_feature=="Lymphocyte Infiltration",]

#############################################################################

assoc_IIC_n = "out/IIC_onsetAgeBinary_assoc.txt"
assoc_IIC = read.table(sep="\t",header=T,file=assoc_IIC_n, stringsAsFactors=FALSE)

assoc_IIC$ID = "Immune Infiltrate Composition"

assoc_IIC$immune_feature[assoc_IIC$immune_feature=="Macrophages_M1"] = "M1 Macrophages"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="Macrophages_M2"] = "M2 Macrophages"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="B_Cells_Naive"] = "B Cells, Naive"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="T_Cells_CD8"] = "CD8+ T Cells"
assoc_IIC$immune_feature[assoc_IIC$immune_feature=="T_Cells_CD4_Naive"] = "CD4+ T Cells, Naive"

m1_macro = assoc_IIC[assoc_IIC$immune_feature=="M1 Macrophages",]
m2_macro = assoc_IIC[assoc_IIC$immune_feature=="M2 Macrophages",]
b_cell_n = assoc_IIC[assoc_IIC$immune_feature=="B Cells, Naive",]
cd8_cell = assoc_IIC[assoc_IIC$immune_feature=="CD8+ T Cells",]
cd4_cell = assoc_IIC[assoc_IIC$immune_feature=="CD4+ T Cells, Naive",]

#############################################################################

# set working directory for dev and debug
bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/Demographics"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

survival_multivar = read.table(sep="\t", header=T, file="../../data/Table_SamsteinCox_YoungOnset_Multivar.txt", stringsAsFactors=FALSE)
survival_multivar$var = "Multivar."
survival_multivar$cancer[survival_multivar$cancer=="COADREAD"] = "COAD"

survival_univar = read.table(sep="\t", header=T, file="../../data/Table_SamsteinCox_YoungOnset_Univar.txt", stringsAsFactors=FALSE)
survival_univar$var = "Univar."
survival_univar$cancer[survival_univar$cancer=="COADREAD"] = "COAD"

survival_combo <- rbind(survival_multivar, survival_univar)
survival_combo$ID = "Hazard Ratio"

p = ggplot(data=survival_combo, mapping = aes(x = var, y = reorder(cancer,desc(cancer)), fill = exp_coef))
p = p + facet_wrap( .~ID, drop=T)#,scales="free_y",nrow=2)
p = p + geom_tile()
p = p + geom_text(aes(label=ifelse(p.value<0.05,"***",NA)), size=6.5, vjust=0.75)
p = p + geom_text(aes(label=ifelse(p.value>=0.05 & p.value<0.10,"**",NA)), size=6.5, vjust=0.75)
p = p + geom_text(aes(label=ifelse(p.value>=0.10 & p.value<0.15,"*",NA)), size=6.5, vjust=0.75)
p = p + theme_bw()
p = p + scale_fill_gradient2(low = "dodgerblue4",mid = "white",high = "red",midpoint = 1,guide = "colourbar",aesthetics = "fill")
p = p + theme(panel.grid.major = element_blank())
p = p + xlab(element_blank()) + ylab(label = "cancer")
p = p + guides(fill=guide_colorbar(title="hazard\nratio",barheight=7))
p = p + theme(axis.text.y = element_text(colour="black", size=10, vjust=0.5))
#p = p + theme(legend.position = "top")
p
fn = 'out/Survival_heatmap.pdf'
ggsave(fn,h=3.75, w=3.25,useDingbats=F)

#############################################################################

lymphocyte_survival = merge(lymphocyte_infil, survival_multivar, by="cancer")
temp = round(summary(lm(coef ~ coefficient, data=lymphocyte_survival))$r.squared, 2)

formula = y ~ x
p = ggplot(data = lymphocyte_survival, aes(x = coefficient, y = coef, color=cancer))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = formula)
p = p + geom_point(size=4, alpha=0.7)
p = p + annotate("text", x=0.1, y=0.5, label=paste("R2 =", temp, sep = " "))
p = p + labs(x="Lymphoctye infiltration", y="Survival (coefficient)")
p = p + getPCACancerColor()
p = p + theme_bw()
p = p + guides(color=FALSE)
p
fn = 'out/LymphocyteSurvival.pdf'
ggsave(fn,h=2.25, w=2.25,useDingbat=F)

#############################################################################

m1_macro_survival = merge(m1_macro, survival_multivar, by="cancer")
temp = round(summary(lm(coef ~ coefficient, data=m1_macro_survival))$r.squared, 2)

formula = y ~ x
p = ggplot(data = m1_macro_survival, aes(x = coefficient, y = coef, color=cancer))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = formula)
p = p + geom_point(size=4, alpha=0.7)
p = p + annotate("text", x=0.002, y=0.5, label=paste("R2 =", temp, sep = " "))
p = p + labs(x="M1 Macrophages", y="Survival (coefficient)")
p = p + getPCACancerColor()
p = p + theme_bw()
p
fn = 'out/M1MacroSurvival.pdf'
ggsave(fn,h=2.25, w=3.25,useDingbat=F)

#############################################################################

m2_macro_survival = merge(m2_macro, survival_multivar, by="cancer")
temp = round(summary(lm(coef ~ coefficient, data=m2_macro_survival))$r.squared, 2)

formula = y ~ x
p = ggplot(data = m2_macro_survival, aes(x = coefficient, y = coef, color=cancer))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = formula)
p = p + geom_point(size=4, alpha=0.7)
p = p + annotate("text", x=0, y=0.45, label=paste("R2 =", temp, sep = " "))
p = p + labs(x="M2 Macrophages", y="Survival (coefficient)")
p = p + getPCACancerColor()
p = p + theme_bw()
p
fn = 'out/M2MacroSurvival.pdf'
ggsave(fn,h=2.25, w=3.25,useDingbat=F)

#############################################################################

b_cell_survival = merge(b_cell_n, survival_multivar, by="cancer")
temp = round(summary(lm(coef ~ coefficient, data=b_cell_survival))$r.squared, 2)

formula = y ~ x
p = ggplot(data = b_cell_survival, aes(x = coefficient, y = coef, color=cancer))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = formula)
p = p + geom_point(size=4, alpha=0.7)
p = p + annotate("text", x=0.0075, y=0.5, label=paste("R2 =", temp, sep = " "))
p = p + labs(x="Naive B Cells", y="Survival (coefficient)")
p = p + getPCACancerColor()
p = p + theme_bw()
p
fn = 'out/BcellSurvival.pdf'
ggsave(fn,h=2.25, w=3.25,useDingbat=F)

#############################################################################

cd4_cell_survival = merge(cd4_cell, survival_multivar, by="cancer")
temp = round(summary(lm(coef ~ coefficient, data=cd4_cell_survival))$r.squared, 2)

formula = y ~ x
p = ggplot(data = cd4_cell_survival, aes(x = coefficient, y = coef, color=cancer))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = formula)
p = p + geom_point(size=4, alpha=0.7)
p = p + annotate("text", x=0.001, y=0.5, label=paste("R2 =", temp, sep = " "))
#p = p + labs(x="Naive CD4+ T Cells", y="Survival (coefficient)")
p = p + labs(x="Naive CD4+ T Cells", y=element_blank())
p = p + getPCACancerColor()
p = p + theme_bw()
p = p + guides(color=FALSE)
p
fn = 'out/CD4cellSurvival.pdf'
ggsave(fn,h=2.25, w=2.125,useDingbat=F)

#############################################################################

cd8_cell_survival = merge(cd8_cell, survival_multivar, by="cancer")
temp = round(summary(lm(coef ~ coefficient, data=cd8_cell_survival))$r.squared, 2)

formula = y ~ x
p = ggplot(data = cd8_cell_survival, aes(x = coefficient, y = coef, color=cancer))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = formula)
p = p + geom_point(size=4, alpha=0.7)
p = p + annotate("text", x=0, y=0.5, label=paste("R2 =", temp, sep = " "))
#p = p + labs(x="CD8+ T Cells", y="Survival (coefficient)")
p = p + labs(x="CD8+ T Cells", y=element_blank())
p = p + getPCACancerColor()
p = p + theme_bw()
p
fn = 'out/CD8cellSurvival.pdf'
ggsave(fn,h=2.25, w=3.125,useDingbat=F)