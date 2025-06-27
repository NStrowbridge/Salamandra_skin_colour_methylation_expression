#DGE analyses script
# by NS
#intended for use with output from 04_edgeR_methylation_transcript_level script

####~load libraries~~~~~~~~~~####
library(rstudioapi)
library(ggplot2)
library(ggrepel)
library(amap)
library(reshape2)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/01_scripts") #set wd
count = 00

###~~output dir~~~~~~~~~~~~~~####
output = "../06_dm_analysis_skin_transcript_level_morph/"
dir.create(output) #make folder for output
setwd(output)

###~~specify data~~~~~~~~~~~~####
ss_file="../02_reference_data/sample_sheet.csv" #sample sheet
#fa_file="../02_reference_data/functional_annotation.csv" #load eggnog mappings
em_file="../04_edgeR_methylation_transcript_level_morph/cpm_skin.csv" #methylation matrix with CPM
de_file_start="../04_edgeR_methylation_transcript_level_morph/results_" #differential methylation output from edgeR
analyses = c("BrovStri", "YelvBroM", "YelvStri")
custom_colors <- c("Yellow" = "yellow", "Striped" = "grey", "Brown" = "brown")

###~~logfile~~~~~~~~~~~~~~~~~####
#log_file=file(paste("05_dge_analysis_",Sys.Date(),".log",sep=""))
#sink(log_file,append=TRUE,type="output")
#sink(log_file,append=TRUE,type="message")
#Sys.time()

####~big loop~~~~~~~~~~~~~~~~####
for (analysis in analyses) {

####~analysis start~~~~~~~~~~####
  count = count + 1
  if (count < 10) {
    anal_dir=paste("0",count,"_",analysis,sep = "") #if the analyses' index is less than 10, add a trailing 0
  } else {
    anal_dir=paste(count,"_",analysis,sep = "") 
  } 
  dir.create(anal_dir)

####~load data~~~~~~~~~~~~~~~####
ss=read.csv(ss_file,row.names=1) #loads sample sheet
em=read.csv(em_file,row.names=1) #loads methylation matrix
de=read.csv(paste(de_file_start,analysis,".csv",sep = ""),row.names=1) #loads differential methylation
setwd(anal_dir) #output will now go to the directory for the current analyses

####~parse data~~~~~~~~~~~~~~####
sortByFDR = function(df, x = "FDR") {
  order_of_x=order(df[,x],decreasing=FALSE)
  df[order_of_x,]
}

###~~methylation matrix~~~~~~~####
em = em[,row.names(ss)]#select and reorder columns in em based on row names in ss
em_scaled=data.frame(scale(data.frame(em)))
em_scaled=na.omit(em_scaled)

##~~~write out em~~~~~~~~~~~~###
write.csv(em, file = "em.csv")
write.csv(em_scaled, file = "em_scaled.csv")
###~~master~~~~~~~~~~~~~~~~~~####
master=merge(em,de,by.x=0,by.y=0) #combine DGE results with CPM to make master
names(master)[1]="SYMBOL"
master$mean=rowMeans(master[,2:(nrow(ss)+1)])
master$mlog10p=-log10(master$PValue)
master$mlog10FDR=-log10(master$FDR)
master$sig=as.factor(master$PValue<0.05&abs(master$logFC)>1.0)
master$sigFDR=as.factor(master$FDR<0.1&abs(master$logFC)>1.0)
row.names(master)=master[,"SYMBOL"]

##~~~sig genes~~~~~~~~~~~~~~~####
master_sig=subset(master,sig==TRUE)
master_sig=sortByFDR(df = master_sig, x = "FDR")
write.csv(master_sig, file = "sig.csv") #write out master sig
sig_genes=master_sig$SYMBOL
em_sig=em[sig_genes,]
em_scaled_sig=em_scaled[sig_genes,]

##~~~sig genes (FDR)~~~~~~~~~####
master_fdr=subset(master,sigFDR==TRUE)
fdr_genes=master_fdr$SYMBOL
em_fdr=em[fdr_genes,]
em_scaled_fdr=em_scaled[fdr_genes,]

##~~~sig up and sig down~~~~~####
master_sig_up=subset(master_sig,logFC>0)
master_sig_down=subset(master_sig,logFC<0)
master_non_sig=subset(master,sig==FALSE)
master_fdr_up=subset(master_fdr, logFC>0)
master_fdr_down=subset(master_fdr, logFC<0)
master_fdr_non_sig=subset(master,sigFDR==FALSE)

# Check if there are any significant upregulated genes
if (nrow(master_sig_up) > 0) {
  master_sig_up$direction="up"
} else {
  master_sig_up <- data.frame() # Create an empty data frame if no significant upregulated genes
}

# Check if there are any significant downregulated genes
if (nrow(master_sig_down) > 0) {
  master_sig_down$direction="down"
} else {
  master_sig_down <- data.frame() # Create an empty data frame if no significant downregulated genes
}

# Combine significant up and downregulated genes
master_sig=rbind(master_sig_up,master_sig_down)

# Add non-significant genes
master_non_sig$direction="ns"
master=rbind(master_sig,master_non_sig)
master$direction=factor(master$direction,levels=c("up","down","ns"))

# Check if there are any significant upregulated genes (FDR)
if (nrow(master_fdr_up) > 0) {
  master_fdr_up$direction="up"
} else {
  master_fdr_up <- data.frame() # Create an empty data frame if no significant upregulated genes (FDR)
}

# Check if there are any significant downregulated genes (FDR)
if (nrow(master_fdr_down) > 0) {
  master_fdr_down$direction="down"
} else {
  master_fdr_down <- data.frame() # Create an empty data frame if no significant downregulated genes (FDR)
}

# Combine significant up and downregulated genes (FDR)
master_fdr_sig=rbind(master_fdr_up,master_fdr_down)

# Add non-significant genes (FDR)
master_fdr_non_sig$direction="ns"
master=rbind(master_fdr_sig,master_fdr_non_sig)
master$direction=factor(master$direction,levels=c("up","down","ns"))
##~~~write out masters~~~~~~~####
write.csv(master, file = "master.csv")
write.csv(master_sig, file = "sig.csv")
write.csv(master_fdr, file = "sig_fdr.csv")
write.csv(master_fdr_non_sig, file = "non_sig_fdr.csv")
write.csv(master_fdr_up, file = "fdr_up.csv")
write.csv(master_fdr_down, file = "fdr_down.csv")
write.csv(master_sig_up, file = "sig_up.csv")
write.csv(master_sig_down, file = "sig_down.csv")

###~~top up and down genes~~~####
order_of_p=order(master_sig_up[,"PValue"],decreasing=FALSE)
master_sig_up=master_sig_up[order_of_p,]
top5_sig_up=master_sig_up[1:5,]
top10_sig_up=master_sig_up[1:10,]
top20_sig_up=master_sig_up[1:20,]
order_of_p=order(master_sig_down[,"PValue"],decreasing=FALSE)
master_sig_down=master_sig_down[order_of_p,]
top5_sig_down=master_sig_down[1:5,]
top10_sig_down=master_sig_down[1:10,]
top20_sig_down=master_sig_down[1:20,]

###~~re-sort master~~~~~~~~~~####
master$direction=factor(master$direction,levels=c("up","down","ns"))
master$sig=factor(master$sig,levels=c("TRUE","FALSE"))
order_of_p=order(master[,"PValue"],decreasing=FALSE)
master=master[order_of_p,]

###~~gene lists~~~~~~~~~~~~~~####
all_genes = row.names(master)
write(all_genes, "gene_universe.txt")
write(sig_genes, "genes_sig.txt")
genes_non_sig = row.names(master_non_sig)
write(genes_non_sig, "genes_non_sig.txt")
genes_sig_up = row.names(master_sig_up)
write(genes_sig_up, "genes_sig_up.txt")
genes_sig_down = row.names(master_sig_down)
write(genes_sig_down, "genes_sig_down.txt")
genes_sig_fdr = row.names(master_fdr)
write(genes_sig_fdr, "genes_sig_fdr.txt")
genes_fdr_up = row.names(master_fdr_up)
write(genes_fdr_up, "genes_fdr_up.txt")
genes_fdr_down = row.names(master_fdr_down)
write(genes_fdr_down, "genes_fdr_down.txt")



####~theme~~~~~~~~~~~~~~~~~~~####

ns_theme=theme(
  plot.title=element_text(size=14),
  axis.text.x=element_text(size=10),
  axis.text.y=element_text(size=10),
  axis.title.x=element_text(size=18),
  axis.title.y=element_text(size=18)
)

####~make plots~~~~~~~~~~~~~~####

####Prepare data for venn diagram
fdr_threshold <- 0.10
fdr.degs <- row.names(master_fdr[which(master_fdr$FDR <= fdr_threshold), ])
write(fdr.degs, "fdr_degs.csv")

###~~MA~~~~~~~~~~~~~~~~~~~~~~####
ma_plot=ggplot(master,aes(x=log10(mean),y=logFC,colour=direction))+
  geom_point(size=0.9)+
  labs(title="MA plot",x="P",y="Log2FC")+
  theme_bw()+
  geom_vline(xintercept=2,linetype="dashed",colour="grey",linewidth=0.5)+
  geom_hline(yintercept=0,linetype="dashed",colour="grey",linewidth=0.5)+
  scale_colour_manual(values=c("red","blue","black"),labels=c("Up","Down","Non-sig"),name="")+
  geom_label_repel(data=top5_sig_up, aes(label=SYMBOL),show.legend=FALSE)+
  geom_label_repel(data=top5_sig_down, aes(label=SYMBOL),show.legend=FALSE)
ggsave("ma.svg",plot = ma_plot)

###~~volcano~~~~~~~~~~~~~~~~~####
volcano_plot=ggplot(master,aes(x=logFC,y=mlog10FDR,colour=direction))+
  geom_point()+
  labs(title="Volcano plot",x="Log Fold Change",y="Log FDR")+
  theme_bw()+
  geom_vline(xintercept=-1,linetype="dashed",colour="grey",size=0.5)+
  geom_vline(xintercept=1,linetype="dashed",colour="grey",size=0.5)+
  geom_hline(yintercept=-log10(0.1),linetype="dashed",colour="grey")+
  scale_colour_manual(values=c("red","blue","black"),labels=c("Up","Down","Non-sig"),name="")+
  geom_label_repel(data=top5_sig_up, aes(label=SYMBOL),show.legend=FALSE)+
  geom_label_repel(data=top5_sig_down, aes(label=SYMBOL),show.legend=FALSE)
ggsave("volcano.svg",plot = volcano_plot)


###~~heatmap~~~~~~~~~~~~~~~~~####

##~~~allsig~~~~~~~~~~~~~~~~~~####

#~~~~make matrix~~~~~~~~~~~~~####
hm.matrix=as.matrix(em_scaled_sig) #convert df to matrix to keep gene names
y.dist=amap::Dist(hm.matrix,method="spearman")
y.cluster=hclust(y.dist,method="average")
y.dd=as.dendrogram(y.cluster)
y.dd.reorder=reorder(y.dd,0,FUN="average")
y.order=order.dendrogram(y.dd.reorder)
hm.matrix_clustered=hm.matrix[y.order,]
hm.matrix_clustered=melt(hm.matrix_clustered) #melt the matrix
names(hm.matrix_clustered)=c("gene","sample","methylation")

#~~~~make the heatmap~~~~~~~~####
colours=c("blue","grey","red")
colorRampPalette(colours)(200)
hm=ggplot(hm.matrix_clustered,aes(x=sample,y=gene,fill=methylation))+
  geom_tile()+
  scale_fill_gradientn(colours=colorRampPalette(colours)(200))+
  ylab("")+
  xlab("")+
  theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
          legend.title=element_blank(),legend.spacing = unit(0.25,"cm"),
          axis.text.x=element_text(angle=45,hjust=1))
ggsave("hm.svg",plot = hm)

###~~boxplots~~~~~~~~~~~~~~~~####

##~~~top 10~~~~~~~~~~~~~~~~~~####
top10=master[1:10,]
candidate_genes=as.vector(row.names(top10)) #get top 10 genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot_candidate_genes=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("boxplot_top10.svg", boxplot_candidate_genes)

#~~~~make faceted boxplot~~~~####
faceted_boxplot_top10=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_boxplot_top10.svg", faceted_boxplot_top10)

##~~~top 20~~~~~~~~~~~~~~~~~~####
top20=master[1:20,]
candidate_genes=as.vector(row.names(top20)) #get top 20 genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(log2(em[candidate_genes,])))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot_candidate_genes_20=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("boxplot_top20.svg",plot = boxplot_candidate_genes_20)

#~~~~make faceted boxplot~~~~####
faceted_boxplot_top20=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_boxplot_top20.svg", plot = faceted_boxplot_top20)

##~~~top 5 up~~~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top5_sig_up)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("5up.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_5up.svg", plot = boxplot)

##~~~top 10 up~~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top10_sig_up)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("10up.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_10up.svg", plot = boxplot)

##~~~top 20 up~~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top20_sig_up)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("20up.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_20up.svg", plot = boxplot)

##~~~top 5 down~~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top5_sig_down)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("5down.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_5down.svg", plot = boxplot)

##~~~top 10 down~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top10_sig_down)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("10down.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_10down.svg", plot = boxplot)

##~~~top 20 down~~~~~~~~~~~~~####
candidate_genes=as.vector(row.names(top20_sig_down)) #get top genes as a vector

#~~~~make gene table~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot~~~~~~~~~~~~####
boxplot=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("20down.svg",plot = boxplot)

#~~~~make faceted boxplot~~~~####
boxplot=ggplot(gene_data.m,aes(y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+   
  scale_fill_manual(values = custom_colors)+
  facet_wrap(~variable,ncol=5)
ggsave("faceted_20down.svg", plot = boxplot)

####~close loop~~~~~~~~~~~~~~####
setwd("..")
}

##~~~PMEL~~~~~~~~~~~~~~~~~~####
PMEL=master["rb_78195",]
candidate_genes=as.vector(row.names(PMEL)) #get top 10 genes as a vector

#~~~~make gene table for PMEL~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot of PMEL~~~~~~~~~~~~####
boxplot_candidate_genes=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme_classic()+
  scale_fill_manual(values=c("grey3", "brown4", "yellow3")) +
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("boxplot_PMEL.svg", boxplot_candidate_genes)

##~~~GPNMB~~~~~~~~~~~~~~~~~~####
GPNMB=master["rb_27044",]
candidate_genes=as.vector(row.names(GPNMB)) #get top 10 genes as a vector

#~~~~make gene table for GPNMB~~~~~~~~~####
gene_data=data.frame(t(em_scaled[candidate_genes,]))
gene_data$sample_group=ss$Condition
gene_data.m=melt(gene_data,id.vars = "sample_group")

#~~~~make boxplot of GPNMB~~~~~~~~~~~~####
boxplot_candidate_genes=ggplot(gene_data.m,aes(x=variable,y=value,fill=sample_group))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  theme_classic()+
  scale_fill_manual(values=c("grey3", "brown4", "yellow3")) +
  theme(axis.text.x=element_text(angle=45,hjust=1))#must be placed after all other theme, rotates x axis text
ggsave("boxplot_GPNMB.svg", boxplot_candidate_genes)

####~end of script~~~~~~~~~~~####
closeAllConnections()

