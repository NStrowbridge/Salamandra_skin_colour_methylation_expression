# Variance partitioning analysis
# by Nic
# for analysis of methylation in yellow, brown and striped morphs

####~load libraries~~~~~~~~~~####
library(variancePartition)
library(lme4)
library(tximport)
library(car)
library("limma")
library("edgeR")


####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Local_Desktop/Molecular_mechanisms_of_colouration/Differential_methylation/01_scripts") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../03.1_variance_partitioning_methylation/" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept
saldata = "../04_n_mod_reads/" #specify where the modified read data is

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 1)

##~~~import txname by geneid~~~~~~~~~~~~#### Y
txname_geneid = read.csv("../02_reference_data/txname_geneid.csv")
txname_geneid <- txname_geneid[-c(1)]

##~~~load filtered_n_mod_reads files~~~~~~~~####
filtered_n_mod_reads_m6anet = read.csv("../04_n_mod_reads/filtered_n_mod_reads_df.csv", header = TRUE)
filtered_n_mod_reads_m6anet <- filtered_n_mod_reads_m6anet[,-1]
rownames(filtered_n_mod_reads_m6anet) <- filtered_n_mod_reads_m6anet[,1]

####~~~~identify genes that pass expression cutoff~~~~#####
isexpr <- rowSums(cpm(filtered_n_mod_reads_m6anet) > 1) >= 0.5 * ncol(filtered_n_mod_reads_m6anet)

####~~~~create data structure with only expressed genes~~~~#####
gExpr <- DGEList(counts = filtered_n_mod_reads_m6anet[isexpr, ])

####~~~~Perform TMM normalization~~~~#####
gExpr <- calcNormFactors(gExpr)

####~~~~Specify variables to be included in the voom() estimates of uncertainty~~~~#####
design <- model.matrix(~Sequencing.kit, ss)

####~~~~Estimate precision weights for each gene and sampley~~~~#####
vobjGenes <- voom(gExpr, design)

####~~~~Fit the variance partitioning model~~~~####
form <- ~ (1|Condition) + (1|Sequencing.kit) + (1|Individual) + (1|Morph.landmark) + (1|Morph)
varPart <- fitExtractVarPartModel(vobjGenes, form, ss)
plotVarPart(varPart)
print(varPart)
png("variance_partitioning_plot.png")
plotVarPart(varPart)
dev.off()


