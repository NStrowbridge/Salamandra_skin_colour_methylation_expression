#DM analysis script
# by Nic
#intended for use with output from 04_edgeR_methylation_transcript_level

####~load libraries~~~~~~~~~~####
library(rstudioapi)
library(reshape2)
library(ggplot2)
library(ggrepel)

####~housekeeping~~~~~~~~~~~~####

rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/01_scripts") #set wd to Scripts folderdir.create("../05_density_pca") #make folder for output
dir.create("../05_density_pca_skin_colour_methylation_site_level") #make folder for output
setwd("../05_density_pca_skin_colour_methylation_site_level") #set output folder as wd

###~~specify data~~~~~~~~~~~~####

ss_file=("../02_reference_data/sample_sheet.csv") #sample sheet
em_file_skin=("../04_edger_methylation_site_level/cpm_skin.csv") #expression matrix with CPM
log2cpm_file_skin=("../04_edger_methylation_site_level/log2cpm_skin.csv")


###~~logfile~~~~~~~~~~~~~~~~~####

log_file=file(paste("02_density_and_pca_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")
Sys.time()

####~load data~~~~~~~~~~~~~~~####

em_skin=read.csv(em_file_skin,row.names=1) #loads expression matrix
ss=read.csv(ss_file,row.names=1) #loads sample sheet
log2cpm_skin=read.csv(log2cpm_file_skin,row.names = 1)

####~parse data~~~~~~~~~~~~~~####

#select and reorder columns in em based on row names in ss
em_skin = em_skin[,row.names(ss)]

#make a log10 expression matrix
log10cpm_skin = log10(em_skin)

#make a scaled expression matrix
em_scaled_skin=data.frame(t(scale(data.frame(t(em_skin)))))
em_scaled_skin=na.omit(em_scaled_skin)

#melt the em matrices
em.m_skin=melt(em_skin)
em_scaled.m_skin=melt(em_scaled_skin)
log2cpm.m_skin=melt(log2cpm_skin)
log10cpm.m_skin=melt(log10cpm_skin)

###~~save data files~~~~~~~~~####
write.csv(ss,file="ss.csv")
write.csv(em_skin,file="em_skin.csv")
write.csv(em_scaled_skin,file="em_scaled_skin.csv")

####~theme~~~~~~~~~~~~~~~~~~~####

js_theme=theme(
  plot.title=element_text(size=14),
  axis.text.x=element_text(size=10),
  axis.text.y=element_text(size=10),
  axis.title.x=element_text(size=18),
  axis.title.y=element_text(size=18)
)

####~make plots~~~~~~~~~~~~~~####

###~~density~~~~~~~~~~~~~~~~~####

##~~~faceted~~~~~~~~~~~~~~~~~####
density_plot_skin=ggplot(em.m_skin,aes(x=log10(value),colour=variable))+
  geom_density(alpha=0.75)+
  facet_wrap(~variable,ncol=6)+ #make ncol a variable determined by sample no for generalisability
  js_theme+
  theme(strip.background=element_rect(fill="transparent",linewidth=0),
        legend.position="none")+
  labs(x="Log10(CPM)",y="Density")
#~~~~save plot~~~~~~~~~~~~~~~####
svglite("density_faceted_skin.svg")
print(density_plot_skin)
dev.off()

##~~~combined~~~~~~~~~~~~~~~~####
density_plot2_skin=ggplot(em.m_skin,aes(x=log10(value),colour=variable))+
  geom_density(alpha=0.75)+
  js_theme+
  theme(strip.background=element_rect(fill="transparent",linewidth=0),
        legend.position="right")+
  labs(x="Log10(CPM)",y="Density")

#~~~~save plot~~~~~~~~~~~~~~~####
svglite("density_overlapping_skin.svg", width = 4, height = 4)
print(density_plot2_skin)
dev.off()

###~~boxplots~~~~~~~~~~~~~~~~####

##~~~scaled CPM~~~~~~~~~~~~~~####
boxplot_skin_cpm=ggplot(em_scaled.m_skin,aes(y=value,fill=variable))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  labs(title = "Expression (All Samples)", y = "Scaled CPM") +
  ylim(-5,5)
ggsave("density_boxplot_scaled_cpm_skin.svg", plot = boxplot_skin_cpm)

##~~~log2 CPM~~~~~~~~~~~~~~~~####
boxplot_skin_cpm2=ggplot(log2cpm.m_skin,aes(y=value,fill=variable))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  labs(title = "Expression (All Samples)", y = "Log2 CPM")
ggsave("density_boxplot_log2cpm_skin.svg", plot = boxplot_skin_cpm2)

##~~~log10 CPM~~~~~~~~~~~~~~~####
boxplot_skin_cpm10=ggplot(log10cpm.m_skin,aes(y=value,fill=variable))+ 
  geom_boxplot(outlier.size=0,show.legend=TRUE)+
  labs(title = "Expression (All Samples)", y = "Log10 CPM")
ggsave("density_boxplot_log10cpm_skin.svg", plot = boxplot_skin_cpm10, width = 10, height = 8)

###~~PCAs~~~~~~~~~~~~~~~~~~~~####
em_matrix_skin=t(as.matrix(sapply(em_scaled_skin,as.numeric)))
pca_skin=prcomp(em_matrix_skin)
pca_coord_skin=data.frame(pca_skin$x)

##~~~PCs 1 & 2~~~~~~~~~~~~~~~####

#~~~~make axis labels~~~~~~~~####
vars=apply(pca_skin$x,2,var)
prop_x=round(vars["PC1"]/sum(vars),4)*100
prop_y=round(vars["PC2"]/sum(vars),4)*100
x_axis_label=paste("PC1","(",prop_x,"%)",sep="")
y_axis_label=paste("PC2","(",prop_y,"%)",sep="")

#~~~~make plot~~~~~~~~~~~~~~~####
pca_plot_1_2_skin=ggplot(pca_coord_skin,aes(x=PC1,y=PC2,colour=ss$Condition))+
  geom_point(aes(color=ss$Condition))+
  scale_color_manual(values=c("grey3", "brown4", "yellow3"))+
  geom_text_repel(aes(label=row.names(ss)),show.legend=FALSE)+
  labs(title="PCA",x=x_axis_label,y=y_axis_label)
ggsave(file="pca_1_2_skin.svg", plot = pca_plot_1_2_skin, width = 10, height = 8)
pca_plot_1_2_skin

##~~~PCs 3 & 4~~~~~~~~~~~~~~~####

#~~~~make axis labels~~~~~~~~####
vars=apply(pca_skin$x,2,var)
prop_x=round(vars["PC3"]/sum(vars),4)*100
prop_y=round(vars["PC4"]/sum(vars),4)*100
x_axis_label=paste("PC3","(",prop_x,"%)",sep="")
y_axis_label=paste("PC4","(",prop_y,"%)",sep="")

#~~~~make plot~~~~~~~~~~~~~~~####
pca_plot_3_4_skin=ggplot(pca_coord_skin,aes(x=PC3,y=PC4,colour=ss$Condition))+
  geom_point(aes(color=ss$Condition))+
  scale_color_manual(values=c("grey3", "brown4", "yellow3"))+
  geom_text_repel(aes(label=row.names(ss)),show.legend=FALSE)+
  labs(title="PCA",x=x_axis_label,y=y_axis_label)
ggsave(file="pca_3_4_skin.svg", plot = pca_plot_3_4_skin, width = 10, height = 8)
pca_plot_3_4_skin

##~~~PCs 5 & 6~~~~~~~~~~~~~~~####

#~~~~make axis labels~~~~~~~~####
vars=apply(pca_skin$x,2,var)
prop_x=round(vars["PC5"]/sum(vars),4)*100
prop_y=round(vars["PC6"]/sum(vars),4)*100
x_axis_label=paste("PC5","(",prop_x,"%)",sep="")
y_axis_label=paste("PC6","(",prop_y,"%)",sep="")

#~~~~make plot~~~~~~~~~~~~~~~####
pca_plot_5_6_skin=ggplot(pca_coord_skin,aes(x=PC5,y=PC6,colour=ss$Condition))+
  geom_point(aes(color=ss$Condition))+
  scale_color_manual(values=c("grey3", "brown4", "yellow3"))+
  geom_text_repel(aes(label=row.names(ss)),show.legend=FALSE)+
  labs(title="PCA",x=x_axis_label,y=y_axis_label)
ggsave(file="pca_5_6_skin.svg", plot = pca_plot_5_6_skin, width = 10, height = 8)
pca_plot_5_6_skin


####~end of script~~~~~~~~~~~####
closeAllConnections()