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
output = "../04_edgeR_methylation_transcript_level_striyelvyelyel" #specify where the output should go
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

# Create a new column in ss based on Morph.landmark conditions
ss$skin_colour <- ifelse(ss$Morph.landmark %in% c("Striped_dor", "Striped_yellow"), "Striped_yellow",
                         ifelse(ss$Morph.landmark %in% c("Yellow_lat", "Yellow_dor"), "Yellow_Yellow", NA))

# Filter the dataframe to include only the relevant categories
filtered_ss <- ss[ss$skin_colour %in% c("Striped_yellow", "Yellow_Yellow"), ]


# Filter the cts dataframe based on the filtered_ss samples and match columns
filtered_n_methylated_reads <- n_methylated_reads[, filtered_ss$RunID]
filtered_n_methylated_reads <- filtered_n_methylated_reads[, match(filtered_ss$RunID, colnames(filtered_n_methylated_reads))]

# Prepare DGEList with the filtered data
filtered_n_methylated_reads <- as.matrix(filtered_n_methylated_reads)
y <- DGEList(filtered_n_methylated_reads, group = filtered_ss$skin_colour)
y$samples$batch <- filtered_ss$Sequencing.kit
y$samples$Individual <- filtered_ss$Individual
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

####~GLM analysis of DGE~~~~~####

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
  Yel_strivYel_yel = Striped_yellow-Yellow_Yellow, 
  levels = design)
qlf.Yel_strivYel_yel = glmQLFTest(fit, contrast=my.contrasts[,"Yel_strivYel_yel"])

###~~get DGE results (Morph,landmark)~~~~~~~~~####
res.Yel_strivYel_yel = topTags(qlf.Yel_strivYel_yel, n=nrow(y), sort.by = "PValue")


##~~~write out DGE tables (morph.landmark~~~~####
write.csv(res.Yel_strivYel_yel, "results_Yel_strivYel_yel.csv")


####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()

