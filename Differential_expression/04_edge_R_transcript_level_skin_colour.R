# EdgeR analysis script 
# by Nic
# for analysis of DGE in yellow, brown and striped morphs

####~load libraries~~~~~~~~~~####
library(BiocManager)
library(GenomicFeatures)
library(tximport)
library(edgeR)
library(svglite)
library(sva)
library(limma)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../04_edge_R_transcript_level_skin_colour/" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept
saldata = "../03_salmon_data/" #specify where the salmon quant files are

###~~logfile~~~~~~~~~~~~~~~~~####
#log_file=file(paste("01_edgeR_",Sys.Date(),".log",sep=""))
#sink(log_file,append=TRUE,type="output")
#sink(log_file,append=TRUE,type="message")

####~load data~~~~~~~~~~~~~~~####

###~~sample sheet~~~~~~~~~~~~####
ss = read.csv(paste(refdata,"sample_sheet.csv",sep = ""), row.names = 1)

###~~quant files~~~~~~~~~~~~~####

##~~~import txname by geneid~~~~~~~~~~~~#### Y
txname_geneid = read.csv("../02_reference_data/txname_geneid.csv")
txname_geneid <- txname_geneid[-c(1)]

##~~~load quant files~~~~~~~~####
salmonQuantFiles = file.path("../03_salmon_data",paste(ss$RunID),"quant.sf") #makes a list of filepaths to the quant data
names(salmonQuantFiles) = row.names(ss) #associate filepaths with sampleIDs from ss 
txi = tximport(salmonQuantFiles, type = "salmon", txOut = TRUE) #import salmon quant files for DGE, specificy txout to import without gene level summarization
cts = txi$counts #get gene counts for DGE
write.csv(cts, "salmon_counts.csv") #save gene counts

# Ensure column names match sample identifiers
sample_ids <- colnames(cts)
ss <- ss[match(sample_ids, ss$RunID), ]

####~prepare DGEList (condiition = skin colour) ~~~~~~~~~####
cts <- as.matrix(cts)
y = DGEList(cts, group = ss$Condition)
y$samples$batch <- ss$Sequencing.kit
y$samples$Individual <- ss$Individual
design = model.matrix(~ group + batch, data = y$samples) #design for filtering
keep1 = filterByExpr(y, design, min.count = 2) 
y = y[keep1, ] #should keep only genes with ~2+ reads in at least one group for Direct RNA
y$counts <- ComBat_seq(y$counts, batch = y$samples$batch )
y = normLibSizes(y) #normalise DGEList for library size
y = estimateDisp(y, design)
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

