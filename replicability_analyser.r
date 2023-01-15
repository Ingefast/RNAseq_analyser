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

library(corrgram)
library(lattice)

##########################################################################################
##########################################################################################
#######################    MULTIVARIATE ANALYSIS       ###################################
##########################################################################################
##########################################################################################


data1 <- read.table("total.gene.TPM.txt", header = TRUE, row.names = 1)

head(data1)
dim(data1)
summary(data1)

data1 <- t(data1)

dim(data1)


#reduces the n of bins for experimental temporal purposes
#data1 <- data1[ , 1:10000]
dim(data1)

#This removes columns with only zeroes
data1 <- data1[, colSums(data1 != 0, na.rm = TRUE) > 0]


dim(data1)


##########################################################################
#########################    THE PLOTTING      ###########################
##########################################################################


#making the ordination with PCA
#ord <- rda(data1)

#making the ordination with DCA
#ord <- decorana(data1)

#making the ordination with NMDS
ord <- metaMDS(data1,autotransform = FALSE, zerodist = "ignore")



plot(ord, display = "sites", type="n", main = "NMDS - Genes")

#Defines a series of symbols
#vpch<- c('M', "1", "2", "3", "1", "2", "3" )
vpch<- c(15, 17, 19)

#Defines a series of colors
vcols <- c("gold", "gold", "gold", "red", "red", "red")

#Defines the names of the samples for plot identity
ftypes<-c( "wt_rep1", "wt_rep2", "wt_rep3", "mutant_rep1", "mutant_rep2", "mutant_rep3")

#Shows the points labeled after Plot identity
points(ord, col = rep(vcols,each=nrow(data1)/length(unique(ftypes))), cex = 1.5, pch = rep(vpch,each=nrow(data1)/length(unique(ftypes))))

#creates a vector for the category names to be shown in the legend
unique(ftypes)->leg.txt

#creates a legend inside the diagram
legend("bottomright", cex=0.8, ncol=1, inset=.02, title="Plots", legend=leg.txt, col=vcols, pch=vpch, horiz=FALSE)



##########################################################################################
##########################################################################################
#######################    CORRELOGRAM                 ###################################
##########################################################################################
##########################################################################################


total <- read.table("total.gene.TPM.txt", header = TRUE, row.names = 1)

total <- total[c(-1)] 
head(total)
total[total == 0] <- NA
head(total)

total <- log(total)
head(total)

#this removes outliers
#total <- ifelse(total <= -5, NA, total)
#total <- ifelse(total <= -5, NA, total)


sample_name <- c("gene_expression_corr")

png(paste(sample_name,".png", sep=""), width=13, height=13, units="in", res=150)


#corrgram(total[1:3], order = TRUE, lower.panel = panel.ellipse, upper.panel=panel.pts, text.panel=panel.txt, diag.panel=panel.minmax, main = "Correlation") 

#adds the color scale
Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")


#this tells the panel.cor function to consider only pairwise complete obs otherwise you #get no corr coeff

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


#this only plots the histone marks                  
pairs(main = sample_name, total[ , c(1:6)], xlim=c(0, 10), ylim=c(0, 10), lower.panel = function(...) smoothScatter(..., colramp = Lab.palette, nrpoints = 0, add = TRUE), upper.panel = panel.cor, diag.panel = panel.minmax)

dev.off()
