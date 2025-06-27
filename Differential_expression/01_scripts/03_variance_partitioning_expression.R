# Variance partitioning analysis
# by Nic
# for analysis of DGE in yellow, brown and striped morphs

####~load libraries~~~~~~~~~~####
library(variancePartition)
library(lme4)
library(tximport)
library(car)
library("limma")
library("edgeR")


####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../03_variance_partitioning_with_txname_by_geneid/" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept
saldata = "../03_salmon_data/" #specify where the salmon quant files are

###~~logfile~~~~~~~~~~~~~~~~~####
log_file=file(paste("01_edgeR_",Sys.Date(),".log",sep=""))
sink(log_file,append=TRUE,type="output")
sink(log_file,append=TRUE,type="message")

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 1)
ss$landmark <- sapply(strsplit(as.character(ss$Morph.landmark), "_"), `[`, 2)
###~~quant files~~~~~~~~~~~~~####

##~~~import txname by geneid~~~~~~~~~~~~#### Y
txname_geneid = read.csv("../02_reference_data/txname_geneid.csv")
txname_geneid <- txname_geneid[-c(1)]

##~~~load quant files~~~~~~~~####
salmonQuantFiles = file.path("../03_salmon_data",paste(ss$RunID),"quant.sf") #makes a list of filepaths to the quant data
names(salmonQuantFiles) = row.names(ss) #associate filepaths with sampleIDs from ss 
txi = tximport(salmonQuantFiles, type = "salmon", tx2gene = txname_geneid) #import salmon quant files for DGE, specificy txout to import without gene level summarization
cts = txi$counts #get gene counts for DGE
write.csv(cts, "salmon_counts.csv") #save gene counts

####~~~~identify genes that pass expression cutoff~~~~#####
isexpr <- rowSums(cpm(cts) > 1) >= 0.5 * ncol(cts)

####~~~~create data structure with only expressed genes~~~~#####
gExpr <- DGEList(counts = cts[isexpr, ])

####~~~~Perform TMM normalization~~~~#####
gExpr <- calcNormFactors(gExpr)

####~~~~Specify variables to be included in the voom() estimates of uncertainty~~~~#####
design <- model.matrix(~Batch, ss)

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


