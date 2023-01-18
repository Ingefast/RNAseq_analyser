#Juan Santos, August 2021
## ---------------------------
##
## sRNA.size_distribution_plotter.r
##
## Plot sRNA size distribution of reads mapping into genes and TEs
##
## Author: Juan Santos, SLU, 2020 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## set working directory
setwd("/media/diskc/project_RNAseq_tmp/summary_output")


library("DESeq2")


## reads count data file
countData <- read.table("total.gene.counts.txt", header=TRUE)

#Adds pseudocount if desired
#countData <- countData + 1

head(countData)


## reads metadata file

metaData <- read.table("/example/genomic_reference_files/metadata.txt", header=TRUE)
head(metaData)


######################################################################################
####   Comparison 1: control VS treated; wt_VS_mutant
######################################################################################

## Selecting conditions to be compared

countData_sel <- countData[ , c("id", "wt_rep1", "wt_rep2", "wt_rep3", "mutant_rep1", "mutant_rep2", "mutant_rep3") ]
head(countData_sel)
dim(countData_sel)


metaData_sel <- metaData[ metaData$dex == "wt" |  metaData$dex == "mutant" , ]
head(metaData_sel)
dim(metaData_sel)


## Construct DESEQ DataSet Object
## converting counts to integer mode

dds <- DESeqDataSetFromMatrix(countData=countData_sel, colData=metaData_sel, design=~dex, tidy = TRUE)


## Now we’re ready to run DESEQ function

dds <- DESeq(dds)

## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing


#Take a look at the results table

res <- results(dds)
head(res)
dim(res)

## Summary of differential gene expression
summary(res)


write.table(res, "DESeq2.result.wt_VS_mutant.txt", sep = "\t", quote = FALSE)
