### Juan Santos, 2022_Mars
### This script takes as input files processed by RNA pipeline for nicolas_RNAseq_arf22

#Instructions following partly :
#https://lashlock.github.io/compbio/R_presentation.html

rm(list = ls())
options(scipen=999)

library("DESeq2")

setwd("/Users/jsantos/Documents/PROJECTS/nicolas_RNAseq_arf22_Mac/analysis/Q0_differential_analysis/wt_VS_arf22_hash6")


## reads count data file
countData <- read.table("/Users/jsantos/Documents/PROJECTS/nicolas_RNAseq_arf22_Mac/normalized_output/htseq.REPLICATES_TOTAL/total.gene.counts.txt", header=TRUE)

#countData <- countData + 1

head(countData)



## reads metadata file

metaData <- read.table("/Users/jsantos/Documents/PROJECTS/nicolas_RNAseq_arf22_Mac/analysis/Q0_differential_analysis/metadata.txt", header=TRUE)

head(metaData)



######################################################################################
####   Comparison 1: control VS treated; wt_VS_arf22_hash6     
######################################################################################

## Selecting conditions to be compared

countData_sel <- countData[ , c("id", "wt_rep1", "wt_rep2", "wt_rep3", "arf22_hash6_rep1", "arf22_hash6_rep2", "arf22_hash6_rep3") ]
head(countData_sel)
dim(countData_sel)


metaData_sel <- metaData[ metaData$dex == "wt" |  metaData$dex == "arf22_hash6" , ]
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


write.table(res, "DESeq2.result.wt_VS_arf22_hash6.txt", sep = "\t", quote = FALSE)
