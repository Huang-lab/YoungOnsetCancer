library(ggplot2)
library(reshape2)

###setup work directory
indir_data<-"C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/YoungOnsetCancer/analysis/gene_expression/out"
indir_plot<-"C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/YoungOnsetCancer/analysis/"

outdir_plot<-indir_data

###source some function for ggplot
setwd(indir_plot)
source("global_aes_out.R")

###read in data
setwd(indir_data)
gene_age<-read.table("Cancer_gene_exp_age_assoc1.txt",sep = "\t",stringsAsFactors = F,header = T)
gene_age_binary<-read.table("Cancer_gene_exp_agebinary_assoc1.txt",sep = "\t",stringsAsFactors = F,header = T)

gene_age<-gene_age[!is.na(gene_age$FDR),]
gene_age_binary<-gene_age_binary[!is.na(gene_age_binary$FDR),]

cs<-sort(unique(gene_age$cancer))
length(cs)
cs1<-sort(unique(gene_age_binary$cancer))
length(cs1)
length(intersect(cs,cs1))

fdr_lev<-c(0.01,0.05,0.1,0.15)
sta_ori<-matrix(0,nrow = length(cs),ncol = length(fdr_lev))
rownames(sta_ori)<-cs
colnames(sta_ori)<-fdr_lev
sta_binary<-sta_ori
for (i in 1:length(cs)) {
  tmpgene_age<-gene_age[gene_age$cancer==cs[i],]
  tmpgene_age_binary<-gene_age_binary[gene_age_binary$cancer==cs[i],]
  for (j in 1:length(fdr_lev)) {
    sta_ori[i,j]<-sum(tmpgene_age$FDR<fdr_lev[j])
    sta_binary[i,j]<-sum(tmpgene_age_binary$FDR<fdr_lev[j])
  }#for j
}#for i
sta<-cbind(sta_ori,sta_binary)
colnames(sta)<-c(paste("ori_fdr",fdr_lev,sep = "_"),paste("binary_fdr",fdr_lev,sep = ""))
sta<-rbind(sta,colSums(sta))
rownames(sta)[nrow(sta)]<-"TOTAL"
setwd(outdir_plot)
write.table(sta,"cancer_gene_different_fdr.txt",row.names = T,col.names = T,quote = F,sep = "\t")

setwd(outdir_plot)
for (i in 1:length(fdr_lev)) {
  tmpgene_age<-gene_age[gene_age$FDR<fdr_lev[i],]
  tmpgene_age_binary<-gene_age_binary[gene_age_binary$FDR<fdr_lev[i],]

  p <- ggplot(data=tmpgene_age,aes(y=coefficient,x=cancer,color = cancer))
  p <- p + geom_point(aes(size=-log10(FDR)),alpha=0.8)
  p <- p + getPCACancerColor()
  p <- p + labs(x="cancer",y="coefficient")
  p <- p + geom_hline(yintercept = 0, alpha=0.5)
  p <- p + theme_bw() + 
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 90, vjust = 0.5),axis.text.y = element_text(colour = "black",size = 8),axis.ticks = element_blank())
  p <- p + guides(fill = guide_legend(override.aes = aes(label="")))
  p
  
  ggsave(paste("cancer_gene_age_original_",fdr_lev[i],".pdf",sep = ""), h=7, w=7, useDingbat=F)
  
  p <- ggplot(data=tmpgene_age_binary,aes(y=coefficient,x=cancer,color = cancer))
  p <- p + geom_point(aes(size=-log10(FDR)),alpha=0.8)
  p <- p + getPCACancerColor()
  p <- p + labs(x="cancer",y="coefficient")
  p <- p + geom_hline(yintercept = 0, alpha=0.5)
  p <- p + theme_bw() + 
    theme(axis.text.x = element_text(colour = "black", size = 8, angle = 90, vjust = 0.5),axis.text.y = element_text(colour = "black",size = 8),axis.ticks = element_blank())
  p
  
  ggsave(paste("cancer_gene_age_binary_",fdr_lev[i],".pdf",sep = ""), h=7, w=8, useDingbat=F)
  
  cat(i)
  cat("\n")
}#for i

setwd(outdir_plot)
sta<-read.table("cancer_gene_different_fdr.txt",row.names = 1,header = T,sep = "\t",stringsAsFactors = F)
sta1<-sta[-nrow(sta),]

sta2<-data.frame(rownames(sta1),sta1,stringsAsFactors = F)
colnames(sta2)[1]<-"cancer"
sta3<-melt(sta2,id.vars="cancer")
colnames(sta3)<-c("cancer","FDRlev","count")

p<-ggplot(data = sta3,aes(x=FDRlev,y=cancer,fill=count))+geom_tile(color="grey")
p
p<-p+theme_minimal()+
  theme(axis.text.x = element_text(vjust = 1,hjust = 1,size = 10,angle = 45),axis.title.x = element_blank(),axis.title.y = element_blank(),
        panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
p
p<-p+coord_fixed(ratio = 0.25)
p
p<-p+scale_fill_gradient2(low = "white",high = "red")
p
p<-p+geom_text(aes(x=FDRlev,y=cancer,label=ifelse(count!=0,count,NA)),color="black",size=2.5)
p

ggsave("cancer_gene_different_fdr.pdf", h=7, w=8, useDingbat=F)
