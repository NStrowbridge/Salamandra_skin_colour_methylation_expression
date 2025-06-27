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
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_methylation/01_scripts") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../04_edgeR_methylation_transcript_level_RNA002" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept
saldata = "../04_n_mod_reads/" #specify where the modified read data is

###~~logfile~~~~~~~~~~~~~~~~~####
#log_file=file(paste("01_edgeR_",Sys.Date(),".log",sep=""))
#sink(log_file,append=TRUE,type="output")
#sink(log_file,append=TRUE,type="message")

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 1)

##~~~load filtered_n_mod_reads files~~~~~~~~####
filtered_n_mod_reads_m6anet = read.csv("../03_m6anet_data/n_mod_reads_summedxtranscript_df.csv", header = TRUE)
# Set the first column as row names
rownames(filtered_n_mod_reads_m6anet) <- filtered_n_mod_reads_m6anet$X
# Remove the first column
filtered_n_mod_reads_m6anet <- filtered_n_mod_reads_m6anet[ , -1]

# Function to rename columns
rename_columns <- function(colname) {
  # Extract the sample name
  sample_name <- gsub(".*_m6a_|_total_n_mod_reads", "", colname)
  
  # Determine the suffix
  suffix <- ifelse(grepl("_dor_", colname), "dor", "lat")
  
  # Check if the sample name consists only of numbers
  if (grepl("^[0-9]+$", sample_name)) {
    sample_name <- paste("ELT", sample_name, sep = "_")
  }
  
  # Combine the sample name and suffix
  paste(sample_name, suffix, sep = "_")
}

# Apply the function to all column names
colnames(filtered_n_mod_reads_m6anet) <- sapply(colnames(filtered_n_mod_reads_m6anet), rename_columns)
# Ensure column names match sample identifiers
sample_ids <- colnames(filtered_n_mod_reads_m6anet)
ss <- ss[match(sample_ids, ss$RunID), ]

# Filter the sample sheet to keep only samples with Sequencing.Kit as RNA002
filtered_ss <- ss[ss$Sequencing.kit == "RNA002", ]

# Get the sample IDs that match the filtered sample sheet
filtered_sample_ids <- filtered_ss$RunID
write.csv(filtered_ss, "../02_reference_data/filtered_sample_sheet_RNA002.csv")

# Filter the filtered_n_mod_reads_m6anet dataframe to keep only the columns that match the filtered sample IDs
filtered_n_mod_reads_m6anet <- filtered_n_mod_reads_m6anet[, filtered_sample_ids]

# Print the filtered dataframe
print(filtered_n_mod_reads_m6anet)

write.csv(filtered_n_mod_reads_m6anet, "n_mod_reads_transcript_level_RNA002.csv") #save gene counts

# Create DGEList object
filtered_n_mod_reads_m6anet <- as.matrix(filtered_n_mod_reads_m6anet)
y <- DGEList(filtered_n_mod_reads_m6anet, group = filtered_ss$Condition)
y$samples$Individual <- filtered_ss$Individual
design = model.matrix(~ group, data = y$samples) #design for filtering
keep1 = filterByExpr(y, design, min.count = 2) 
y = y[keep1, ] #should keep only genes with ~2+ reads in at least one group for Direct RNA
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

####~GLM analysis of DGE~~~~~####

###~~fit GLM (skin colour)~~~~~~~~~~~~~~~~~####
# Create the design matrix
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(factor(make.names(y$samples$group)))

# Convert counts to a matrix
counts_matrix <- as.matrix(y$counts)

# Estimate the correlation between duplicate measurements
corfit <- duplicateCorrelation(counts_matrix, design, block = y$samples$Individual)

# Fit the GLM with the estimated correlation
fit = glmQLFit(y, design, correlation = corfit$consensus)

###~~compare groups (skin)~~~~~~~~~~####
my.contrasts = makeContrasts(
  YelvBla = Yellow-Black, #compare yellow skin vs black
  YelvBro = Yellow-Brown,#compare yellow vs brown
  BrovBla = Brown-Black, #compare brown vs black
  levels = design)
qlf.YelvBla = glmQLFTest(fit, contrast=my.contrasts[,"YelvBla"])
qlf.YelvBro = glmQLFTest(fit, contrast=my.contrasts[,"YelvBro"])
qlf.BrovBla = glmQLFTest(fit, contrast=my.contrasts[,"BrovBla"])

###~~get DGE results (skin)~~~~~~~~~####
res.YelvBla = topTags(qlf.YelvBla , n=nrow(y), sort.by = "PValue")
res.YelvBro = topTags(qlf.YelvBro, n=nrow(y), sort.by = "PValue")
res.BrovBla = topTags(qlf.BrovBla, n=nrow(y), sort.by = "PValue")

##~~~write out DGE tables (skin)~~~~####
write.csv(res.YelvBla, "results_YelvBla.csv")
write.csv(res.YelvBro, "results_YelvBro.csv")
write.csv(res.BrovBla, "results_BrovBla.csv")

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()

