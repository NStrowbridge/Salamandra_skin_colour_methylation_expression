####~load libraries~~~~~~~~~~####
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(svglite)
library(pheatmap)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/01_scripts") #set wd to Scripts folder

###~~output directory~~~~~~~~####
output = "../04_DESeq2_transcript_level_skin_colour/" #specify where the output should go
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
salmonQuantFiles = file.path(saldata, paste(ss$RunID), "quant.sf")
names(salmonQuantFiles) = row.names(ss)

txi = tximport(salmonQuantFiles, type = "salmon", txOut = TRUE)


# Extract raw counts
raw_counts <- txi$counts

# Filter low-count genes
keep <- rowSums(raw_counts) > 10
filtered_counts <- raw_counts[keep, ]

# Apply ComBat-Seq
batch <- ss$Sequencing.kit
condition <- ss$Condition
combat_counts <- ComBat_seq(counts = filtered_counts, batch = batch, group = condition)


combat_counts_rounded <- round(combat_counts)

dds <- DESeqDataSetFromMatrix(countData = combat_counts_rounded,
                              colData = ss,
                              design = ~ Condition)


# Run DESeq2
dds <- DESeq(dds)


###~~get normalized counts~~~~~~~~####
norm_counts = counts(dds, normalized = TRUE)
write.csv(norm_counts, "normalized_counts.csv")

log_norm_counts = assay(rlog(dds))
write.csv(log_norm_counts, "log2_normalized_counts.csv")

####~basic plots~~~~~~~~~~~~~~####
svglite("pca_plot.svg", width = 5, height = 5)
plotPCA(rlog(dds), intgroup = "Condition")
dev.off()

svglite("sample_dists.svg", width = 5, height = 5)
sampleDists <- dist(t(log_norm_counts))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)
dev.off()

####~DGE analysis~~~~~~~~~~~~~####
res.YelvBla = results(dds, contrast = c("Condition", "Yellow", "Black"))
res.YelvBro = results(dds, contrast = c("Condition", "Yellow", "Brown"))
res.BrovBla = results(dds, contrast = c("Condition", "Brown", "Black"))

write.csv(as.data.frame(res.YelvBla), "results_YelvBla.csv")
write.csv(as.data.frame(res.YelvBro), "results_YelvBro.csv")
write.csv(as.data.frame(res.BrovBla), "results_BrovBla.csv")

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()

