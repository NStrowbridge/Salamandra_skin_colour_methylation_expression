# Get functional gene annotation using eggnogg protein mappings
# R script
# by Nic

####~load libraries~~~~~~~~~~####
library(phylotools)
library(stringr)
library(readr)
library(rstudioapi)

####~housekeeping~~~~~~~~~~~~####
rm(list=ls()) #clear the environment
setwd("/Users/nicstrowbridge/Desktop/Nic_PhD_files_2/DirectRNA_Colour_bernardezi/Differential_expression/02_reference_data") #set wd to reference data folder

###~~logfile~~~~~~~~~~~~~~~~~####
#log_file=file(paste("02_get_functional_annotation_",Sys.Date(),".log",sep=""))
#sink(log_file,append=TRUE,type="output")
#sink(log_file,append=TRUE,type="message")
#Sys.time()

####~get eggnogg annotation~~####
eggnog = read.delim("eggnogmapper_diamond_output_kmer19.emapper.annotations.tsv", header = FALSE)
colnames(eggnog) = eggnog[5,] #get correct column headers
eggnog = eggnog[-c(1:5,19108,19109,19110),] #delete empty rows and header row
colnames(eggnog)[1] = "query" #remove hashtag from 1st column name
 
###~~make txname*geneid~~~~~~~~~~####
library("dplyr") #load dplyr to change name of specific columns
#First make new column with txname (currently txname ORF data)
eggnog$TXNAME = substr(eggnog$query,1,nchar(eggnog$query)-3)
# Rename seed_ortholog to gene name
eggnog <- eggnog %>% rename("GENEID" = "seed_ortholog")
print(head(eggnog))
txname_geneid = eggnog[,c("TXNAME","GENEID")]

###~~make txname*preferred_name+GO~~~~~~~~~~####
txname_preferredname_GOcat = eggnog[,c("TXNAME","Preferred_name", "GOs")]
txname_preferredname = eggnog[,c("TXNAME","Preferred_name")]

####~map file for topGO~~~~~~####
gene2GO = eggnog[,c("GENEID","GOs")] #make gene2GO annotation list

####~write out annotation~~~~####
write.csv(txname_geneid, file = "../02_reference_data/txname_geneid.csv")
write.csv(txname_preferredname_GOcat, file = "../02_reference_data/txname_preferredname_GOcat.csv")
write.csv(txname_preferredname, file = "../02_reference_data/txname_preferredname.csv")
write_tsv(gene2GO, file = "../02_reference_data/gene2GO.map", col_names = FALSE)

####~fin~~~~~~~~~~~~~~~~~~~~~####
closeAllConnections()

