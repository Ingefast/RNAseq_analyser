![This is an image](/images/rnaseq_analyser_title.png)


# INTRODUCTION

This is a general pipeline for analysis of gene expression RNA sequencing data from Illumina. It has mainly been developed with focus on the *Arabidopsis* TAIR10 genome and some related species, but it should be fully functional with other organisms. The scripts do not require to pass command-line arguments; settings like input data and reference genomic files have to be specified by editing the script in a text editor. Therefore, some very basic knowledge of linux and R is required. The scripts are generously commented in hashes (#) with complementary suggestions and hints (worth to read as a complement to this instruction). It produces basic background output for customised downstream analysis.

# WORKFLOW
![This is an image](/images/flowchart.png)

# SUPPORTED PLATFORMS

Linux, Mac OS

Shell scripts (*.sh) of this software were developed and tested using GNU bash (v4.4.20) in a Ubuntu 18.04 linux system. R scripts were developed using the R console (v4.1.1) in macOS Monterey.

# DEPENDENCIES

The following tools need to be installed and ideally available in the PATH environment. The pipeline is fully functional and tested with the following versions of the packages listed below. Other versions are very likely functional as well, but a detailed compatibility review of older and newer versions has not been done here. 

``fastqc`` (v0.11.9)

``trim_galore`` (v0.6.7)

``hisat2`` (v2.1.0)

``stringtie`` (v2.1.2)

``samtools`` (v1.3.1)

``DESeq2`` (v1.38.2)

# SETTING UP THE WORKING DIRECTORY AND THE GENOMIC REFERENCE FILES

The raw data (Illumina single-end fastq files) should be allocated in sample folders under a parent directory (/**REPLICATES_TOTAL**) following the file structure below. Single replicates are the basic units of this pipeline. Intermediary and final output files will be generated in respective sample folders. In the following example six samples are used: two conditions (*mutant* and *wt*) with three replicates each (*rep1*, *rep2* and *rep3*). 


```
└── REPLICATES_TOTAL
    ├── mutant_rep1
    │   ├── sample_R1.fastq
    │   └── sample_R2.fastq
    ├── mutant_rep2
    │   ├── sample_R1.fastq
    │   └── sample_R2.fastq
    ├── mutant_rep3
    │   ├── sample_R1.fastq
    │   └── sample_R2.fastq
    ├── wt_rep1
    │   ├── sample_R1.fastq
    │   └── sample_R2.fastq
    ├── wt_rep2
    │   ├── sample_R1.fastq
    │   └── sample_R2.fastq
    └── wt_rep3
        ├── sample_R1.fastq
        └── sample_R2.fastq
```

The paths to the working and sample directories, and reference genomic files have to be specified by editing the header of the respective shell script (*.sh).

```
##################################################################################
###########  Reference files to be used and working dir        ###################
##################################################################################

#hisat2 indexes to be used
genome_hisat2_idx="/home/jsantos/genome_ref/TAIR10/hisat2/TAIR10_chr_all"

#annotation GTF file
genes_annot_gtf="/home/jsantos/genome_ref/TAIR10/TAIR10_GFF3_genes.protein_coding.FIVE_CHR.gtf"

##################################################################################

#working directory where all the sample directories are allocated
working_dir="/media/diskc/project_RNAseq_tmp/REPLICATES_TOTAL"

#this defines the samples to be analysed (dir names)
sample_list="mutant_rep1  mutant_rep2  mutant_rep3  wt_rep1  wt_rep2  wt_rep3";
#sample_list="test_rep1";

##################################################################################

```


Ordinary assembly fasta files and ``hisat2`` indexed fasta files (hisat2-build) should be available for the relevant genome and are usually downloadable from general or organism-specific genome repositories like ([TAIR10](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release)).

Annotation files in gtf format can be downloaded/adapted from the TAIR10 **Arabidopsis** genome repository above, but a functional gtf file for only the five chromosomes is provided under the examples folder. To prepare bed files out of gtf or gff3 files is not straightforward. The [`gff2bed`](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html) tool from ``BEDOPS`` suit is an option. Another possibility, often more pragmatic, is to process it with a combination of linux regular expressions and/or manual editing in a text editor.


# INSTALLATION

Shell scripts can be cloned and run directly on a linux server.

```
git clone https://github.com/Ingefast/RNAseq_analyser.git
cd RNAseq_analyser
```


## 1. Quality Control and Size Trimming.

Read quality (``fastqc``) is assessed before and after size trimming (``trim_galore``) if desired. 

Usage:
```
nohup bash RNAseq.QC_Trimmer.sh
```
Two pair-end files with trimmed sequences in fastq format are generated per sample (**trimmed_1.fastq**, **trimmed_2.fastq**).


## 2. Mapping of trimmed reads to the reference genome

Once the pair-end reads are quality trimmed, mapping to the reference genome is performed using ``hisat2``

Usage:
```
nohup bash RNAseq.mapper.sh
```

Once the alignment is sorted and indexed, the number of reads mapping to each gene is counted using the packages ``htseq-count`` and ``stringtie``.

``htseq-count`` produces a table with the number of raw reads for each feature (**counts.htseq.txt**)

``stringtie`` produces tables with expression values normalised in several ways (**abund.stringtie.txt**).


## 3. Summarising 

Usage:
```
nohup bash RNAseq.total_table_maker.sh
```
Produces two contingency tables of expression (genes x samples):

**total.gene.counts.txt**
Table with raw read counts per gene to be used down the line to run a differential analysis of gene expression with ``DESeq2``, and to assess the replicability of the experiment.

**total.gene.TPM.txt**
Table with normalised gene expression values as Transcripts Per Million (TPM) appropriate to diverse downstream analysis and for submission to final NGS data repositories as the processed final output of the study.


## 4. Differential analysis of gene expression 

The ``R`` script **RNAseq.differential_analyser.DESEQ2.r** uses the above indicated table **total.gene.counts.txt** as input to perform an standard differential analysis. The output text table can be processed further in ``R`` or ``excel``.


## 5. Replicability assessment 

With the ``R`` script **replicability_analyser.r** is possible to analyse the tables **total.gene.counts.txt** or **total.gene.TPM.txt** to plot the correlation between particular gene expression values across all the samples simultaneously. This is a good way to check for deviant samples and assess replicability. A scatterplot in the lower diagonal panel is presented and Pearson correlation coefficients in the upper panel (Figure 1B). This script performs also a multivariate analysis using the same input table to evaluate the similarities between samples. By default the samples are ordinated with Nonmetric multidimensional scaling (NMDS) as implemented in the R library vegan, but other alternative popular ordination methods such as Principal component analysis (PCA) (Figure 1B) or redundancy analysis (RDA) are also available.

*Figure 1*. Replicability assessment of gene expression. (A) Correlogram with gene expression values (TPM) in four conditions with three replicates each. (C) PCA ordinaton diagram of the same dataset.

 
# REFERENCES

Modified versions of this pipeline have been used to process the RNA-seq datasets in e.g. the following papers:

1. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.

2. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.

# CONTACT
juan.santos at slu.se
