# EdgeR analysis script 
# by Nic
# for analysis of DGM in yellow, brown and striped morphs

####~load libraries~~~~~~~~~~####
library(BiocManager)
library(GenomicFeatures)
library(tximport)
library(edgeR)
library(svglite)
#install.packages("locfit")
library(locfit)
#BiocManager::install("sva")
library(sva)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("script_folder") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../04_edgeR_methylation_transcript_level_stridorvstrilat" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 1)

##~~~load filtered_n_mod_reads files~~~~~~~~####
n_methylated_reads = read.csv("../04_n_methylated_reads/n_methylated_reads.csv", header = TRUE)
# Set the first column as row names
rownames(n_methylated_reads) <- n_methylated_reads$X
# Remove the first column
n_methylated_reads <- n_methylated_reads[ , -1]

# Ensure column names match sample identifiers
sample_ids <- colnames(n_methylated_reads)
ss <- ss[match(sample_ids, ss$RunID), ]

write.csv(n_methylated_reads, "n_methylated_reads_transcript_level.csv") #save gene counts

# Create DGEList object
n_methylated_reads <- as.matrix(n_methylated_reads)
y <- DGEList(n_methylated_reads, group = ss$Morph.landmark)
y$samples$batch <- ss$Sequencing.kit
y$samples$Individual <- ss$Individual
design = model.matrix(~ group, data = y$samples) #design for filtering
keep1 = filterByExpr(y, design, min.count = 15) 
y = y[keep1, ] #should keep only genes with ~15+ reads in at least one group for Direct RNA methylation data
y$counts <- ComBat_seq(y$counts, batch = y$samples$batch )
y <- normLibSizes(y)
y <- estimateDisp(y, design)
sqrt(y$common.dispersion)

###~~get CPM (skin colour)~~~~~~~~~~~~~~~~~####
cpms = edgeR::cpm(y, offset = y$offset, log = FALSE)
write.csv(cpms, "cpm_skin.csv")

###~~get log2 CPM~~~~~~~~~~~~####
logcpm = cpm(y, log = TRUE)
write.csv(logcpm, "log2cpm_skin.csv")

####~some basic plots~~~~~~~~####
svglite("mds.svg", width = 4, height = 4)
plotMDS(y) #visualise variation between samples
dev.off()


svglite("bcv.svg", width = 4, height = 4)
plotBCV(y) #visualise dispersion estimates
dev.off()

####~GLM analysis of DM~~~~~####
###~~fit GLM (morphxlandmark)~~~~~~~~~~~~~~~~~####
# Create the design matrix
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(factor(make.names(y$samples$group)))

# Convert counts to a matrix
counts_matrix <- as.matrix(y$counts)

# Estimate the correlation between duplicate measurements
corfit <- duplicateCorrelation(counts_matrix, design, block = y$samples$Individual)

# Fit the GLM with the estimated correlation
fit = glmQLFit(y, design, correlation = corfit$consensus)

###~~compare groups (morph.landmark)~~~~~~~~~~####
my.contrasts = makeContrasts(
  Stri_dorvStri_lat = Striped_dor-Striped_lat,
  levels = design)
qlf.Stri_dorvStri_lat = glmQLFTest(fit, contrast=my.contrasts[,"Stri_dorvStri_lat"])
###~~get DGE results (Morph,landmark)~~~~~~~~~####
res.Stri_dorvStri_lat = topTags(qlf.Stri_dorvStri_lat, n=nrow(y), sort.by = "PValue")
##~~~write out DGE tables (morph.landmark~~~~####
write.csv(res.Stri_dorvStri_lat, "results_Stri_dorvStri_lat.csv")
####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()
