##### plotSomaticAssoc_Updated.R #####
# Updated by Will Lee, October 2021

bdir = "~/Box/Huang_lab/manuscripts/YoungOnsetCancer/analysis/somatic_mutation"
setwd(bdir)

# source the files with plotting related functions
source("../global_aes_out.R")

### associations - updated by Will Lee ###

assoc_somatic_n = "out/somatic_likelyfunctional_driver_onsetAgeBinary_assoc_Updated.txt"
assoc_somatic = read.table(sep="\t",header=T,file=assoc_somatic_n, stringsAsFactors=FALSE)

assoc_somatic_n_continuous = "out/somatic_likelyfunctional_loose_onsetAge_assoc_Updated.txt"
assoc_somatic_continuous = read.table(sep="\t",header=T,file=assoc_somatic_n_continuous, stringsAsFactors=FALSE)

assoc_somatic_c = "out/somatic_likelyfunctional_driver_onsetAgeBinary_Clinical.txt"
assoc_somatic_clinical = read.table(sep="\t",header=T,file=assoc_somatic_c, stringsAsFactors=FALSE)

assoc_somatic_c_continuous = "out/somatic_likelyfunctional_driver_onsetAge_Clinical.txt"
assoc_somatic_clin_cont = read.table(sep="\t",header=T,file=assoc_somatic_c_continuous, stringsAsFactors=FALSE)

##################################################################################################################################

# bubble plot, somatic continuous #
p = ggplot(data=assoc_somatic_continuous,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,gene,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "Coefficient") #+ ylim(-2.5,2.6)
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/SomaticAssocOnsetBubble_Updated.pdf'
ggsave(fn,h=5.75, w=5,useDingbat=F)

### volcano plot, somatic continuous ###
p = ggplot(data=assoc_somatic_continuous,aes(x=coefficient,y = -log10(FDR),color = cancer))
p = p + geom_point(alpha=0.7)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,gene,NA)), show.legend = FALSE)
p = p + getPCACancerColor() + xlim(-35,5)
#p = p + labs(x="Cancer",y= "Coefficient")
p = p + geom_vline(xintercept = 0, alpha=0.5)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/SomaticAssocOnsetVolcano.pdf'
ggsave(fn,h=4.25, w=3.5,useDingbat=F)

### volcano plot, somatic continuous clinical ###
p = ggplot(data=assoc_somatic_clin_cont,aes(x=coefficient,y = -log10(FDR),color = cancer))
p = p + geom_point(alpha=0.7)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,geneVar,NA)), colour="black", show.legend = FALSE)
p = p + getPCACancerColor() + xlim(-5,5)
#p = p + labs(x="Cancer",y= "Coefficient")
p = p + geom_vline(xintercept = 0, alpha=0.5)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = 'out/SomaticAssocOnsetVolcano_Clinical.pdf'
ggsave(fn,h=4.25, w=3.5,useDingbat=F)

##################################################################################################################################

# bubble plot, somatic binary #
p = ggplot(data=assoc_somatic,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + scale_size(name = expression(-log10(FDR)), breaks = c(1, 2.5, 5.0, 7.5))
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,gene,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "Coefficient") + ylim(-2.5,2.6)
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/SomaticAssocOnsetBinaryBubble_Updated.pdf'
ggsave(fn,h=6, w=5,useDingbat=F)

# bubble plot, somatic binary clinical #
p = ggplot(data=assoc_somatic_clinical,aes(x=cancer,y = coefficient,color = cancer))
p = p + geom_point(aes(size = -log10(FDR)),alpha=0.8)
p = p + geom_text_repel(aes(label=ifelse(FDR<0.05,geneVar,NA)), size=3, colour="black", show.legend = FALSE)
p = p + getPCACancerColor()
p = p + labs(x="Cancer",y= "Coefficient") + ylim(-5,25)
p = p + theme(legend.position = "top")
p = p + geom_hline(yintercept = 0, alpha=0.1)
p = p + theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=12, angle= 90, vjust=0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())
p
fn = 'out/SomaticAssocOnsetBinaryBubble_Clinical.pdf'
ggsave(fn,h=5.75, w=5,useDingbat=F)