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
output = "../04_edgeR_methylation_transcript_level_morphxlandmark" #specify where the output should go
dir.create(output) #create directory for output
setwd(output) #set the new output directory as the working directory

###~~specify data~~~~~~~~~~~~####
refdata = "../02_reference_data/" #specify where the reference data is kept

###~~logfile~~~~~~~~~~~~~~~~~####
#log_file=file(paste("01_edgeR_",Sys.Date(),".log",sep=""))
#sink(log_file,append=TRUE,type="output")
#sink(log_file,append=TRUE,type="message")

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

####~GLM analysis of DGE~~~~~####

###~~fit GLM (morphxlandmark)~~~~~~~~~~~~~~~~~####
# Create the design matrix
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(factor(make.names(y$samples$group)))

# Convert counts to a matrix
counts_matrix <- as.matrix(y$counts)

# Estimate the correlation between duplicate measurements
# Fit the GLM with the estimated correlation
fit = glmQLFit(y, design)

###~~compare groups (morph.landmark)~~~~~~~~~~####
my.contrasts = makeContrasts(
  Yel_dorvYel_lat = Yellow_dor-Yellow_lat, #compare yellow dorsal vs yellow lateral
  Yel_dorvStri_dor = Yellow_dor-Striped_dor,#compare yellow dorsal vs Striped dorsal
  Yel_dorvStri_lat = Yellow_dor-Striped_lat, #compare yellow dorsal vs Striped lateral
  Yel_latvStri_dor = Yellow_lat-Striped_dor, #compare yellow lateral vs striped dorsal
  Yel_latvStri_lat = Yellow_lat-Striped_lat, #compare yellow lateral vs striped lateral
  Yel_dorvBro_dor = Yellow_dor-Brown_dor, #compare yellow dorsal vs Brown dorsal
  Yel_dorvBro_lat = Yellow_dor-Brown_lat, #compare yellow dorsal vs Brown lateral
  Yel_latvBro_dor = Yellow_lat-Brown_dor, #compare yellow lateral vs Brown dorsal
  Yel_latvBro_lat = Yellow_lat-Brown_lat, #compare yellow lateral vs Brown lateral
  Stri_dortvBro_dor = Striped_dor-Brown_dor, #compare Striped dorsal vs Brown dorsal
  Stri_dortvBro_lat = Striped_dor-Brown_lat, #compare Striped dorsal vs Brown lateral
  Stri_lattvBro_dor = Striped_lat-Brown_dor, #compare Striped lateral vs Brown dorsal
  Stri_lattvBro_lat = Striped_lat-Brown_lat, #compare Striped lateral vs Brown lateral
  Stri_dorvStri_lat = Striped_dor-Striped_lat, #compare Striped dorsal vs Striped lateral
  Bro_dorvBro_lat = Brown_dor-Brown_lat, #compare Brown dorsal vs Brown lateral
  levels = design)
qlf.Yel_dorvYel_lat = glmQLFTest(fit, contrast=my.contrasts[,"Yel_dorvYel_lat"])
qlf.Yel_dorvStri_dor = glmQLFTest(fit, contrast=my.contrasts[,"Yel_dorvStri_dor"])
qlf.Yel_dorvStri_lat = glmQLFTest(fit, contrast=my.contrasts[,"Yel_dorvStri_lat"])
qlf.Yel_latvStri_dor = glmQLFTest(fit, contrast=my.contrasts[,"Yel_latvStri_dor"])
qlf.Yel_latvStri_lat = glmQLFTest(fit, contrast=my.contrasts[,"Yel_latvStri_lat"])
qlf.Yel_dorvBro_dor = glmQLFTest(fit, contrast=my.contrasts[,"Yel_dorvBro_dor"])
qlf.Yel_dorvBro_lat = glmQLFTest(fit, contrast=my.contrasts[,"Yel_dorvBro_lat"])
qlf.Yel_latvBro_dor = glmQLFTest(fit, contrast=my.contrasts[,"Yel_latvBro_dor"])
qlf.Yel_latvBro_lat = glmQLFTest(fit, contrast=my.contrasts[,"Yel_latvBro_lat"])
qlf.Stri_dortvBro_dor= glmQLFTest(fit, contrast=my.contrasts[,"Stri_dortvBro_dor"])
qlf.Stri_dortvBro_lat= glmQLFTest(fit, contrast=my.contrasts[,"Stri_dortvBro_lat"])
qlf.Stri_lattvBro_dor  = glmQLFTest(fit, contrast=my.contrasts[,"Stri_lattvBro_dor"])
qlf.Stri_lattvBro_lat= glmQLFTest(fit, contrast=my.contrasts[,"Stri_lattvBro_lat"])
qlf.Stri_dorvStri_lat = glmQLFTest(fit, contrast=my.contrasts[,"Stri_dorvStri_lat"])
qlf.Bro_dorvBro_lat = glmQLFTest(fit, contrast=my.contrasts[,"Bro_dorvBro_lat"])

###~~get DGE results (Morph,landmark)~~~~~~~~~####
res.Yel_dorvYel_lat = topTags(qlf.Yel_dorvYel_lat, n=nrow(y), sort.by = "PValue")
res.Yel_dorvStri_dor = topTags(qlf.Yel_dorvStri_dor, n=nrow(y), sort.by = "PValue")
res.Yel_dorvStri_lat = topTags(qlf.Yel_dorvStri_lat, n=nrow(y), sort.by = "PValue")
res.Yel_latvStri_dor = topTags(qlf.Yel_latvStri_dor, n=nrow(y), sort.by = "PValue")
res.Yel_latvStri_lat = topTags(qlf.Yel_latvStri_lat, n=nrow(y), sort.by = "PValue")
res.Yel_dorvBro_dor = topTags(qlf.Yel_dorvBro_dor, n=nrow(y), sort.by = "PValue")
res.Yel_dorvBro_lat = topTags(qlf.Yel_dorvBro_lat, n=nrow(y), sort.by = "PValue")
res.Yel_latvBro_dor = topTags(qlf.Yel_latvBro_dor, n=nrow(y), sort.by = "PValue")
res.Yel_latvBro_lat = topTags(qlf.Yel_latvBro_lat, n=nrow(y), sort.by = "PValue")
res.Stri_dortvBro_dor = topTags(qlf.Stri_dortvBro_dor, n=nrow(y), sort.by = "PValue")
res.Stri_dortvBro_lat = topTags(qlf.Stri_dortvBro_lat, n=nrow(y), sort.by = "PValue")
res.Stri_lattvBro_dor = topTags(qlf.Stri_lattvBro_dor, n=nrow(y), sort.by = "PValue")
res.Stri_lattvBro_lat = topTags(qlf.Stri_lattvBro_lat, n=nrow(y), sort.by = "PValue")
res.Stri_dorvStri_lat = topTags(qlf.Stri_dorvStri_lat, n=nrow(y), sort.by = "PValue")
res.Bro_dorvBro_lat = topTags(qlf.Bro_dorvBro_lat, n=nrow(y), sort.by = "PValue")

##~~~write out DGE tables (morph.landmark~~~~####
write.csv(res.Yel_dorvYel_lat, "results_Yel_dorvYel_lat.csv")
write.csv(res.Yel_dorvStri_dor, "results_Yel_dorvStri_dor.csv")
write.csv(res.Yel_dorvStri_lat, "results_Yel_dorvStri_lat.csv")
write.csv(res.Yel_latvStri_dor, "results_Yel_latvStri_dor.csv")
write.csv(res.Yel_latvStri_lat, "results_Yel_latvStri_lat.csv")
write.csv(res.Yel_dorvBro_dor, "results_Yel_dorvBro_dor.csv")
write.csv(res.Yel_dorvBro_lat, "results_Yel_dorvBro_lat.csv")
write.csv(res.Yel_latvBro_dor, "results_Yel_latvBro_dor.csv")
write.csv(res.Yel_latvBro_lat, "results_Yel_latvBro_lat.csv")
write.csv(res.Stri_dortvBro_dor, "results_Stri_dortvBro_dor.csv")
write.csv(res.Stri_dortvBro_lat, "results_Stri_dortvBro_lat.csv")
write.csv(res.Stri_lattvBro_dor, "results_Stri_lattvBro_dor.csv")
write.csv(res.Stri_lattvBro_lat, "results_Stri_lattvBro_lat.csv")
write.csv(res.Stri_dorvStri_lat, "results_Stri_dorvStri_lat.csv")
write.csv(res.Bro_dorvBro_lat, "results_Bro_dorvBro_lat.csv")

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()

