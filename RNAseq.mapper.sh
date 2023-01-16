#!/bin/bash

#Juan Santos, 2022 November

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


### Step 1 . Selecting sRNA size_1830 range and mapping

for sample in $sample_list  ; do

	cd $working_dir/${sample};
	
	hisat2 -p 20 --dta -x $genome_hisat2_idx -1 trimmed_1.fastq -2 trimmed_2.fastq -S trimmed_pe.sam 2>align_summary.txt;
	samtools sort -@ 20 trimmed_pe.sam > trimmed_pe.sorted.bam ;
	samtools index trimmed_pe.sorted.bam;

done
wait


### Step 2 . Counting expression in different ways

for sample in $sample_list  ; do

	cd $working_dir/${sample};
	
#for stringtie see https://bioinformatics.uconn.edu/rnaseq-arabidopsis/#mapping
stringtie -p 20 trimmed_pe.sorted.bam -A abund.stringtie.txt -G $genes_annot_gtf -o counts.stringtie.txt;

#with htseq after https://scilifelab.github.io/courses/rnaseq/labs/
htseq-count -s no -q trimmed_pe.sorted.bam $genes_annot_gtf > counts.htseq.txt;

	#for featureCounts see http://bioinf.wehi.edu.au/featureCounts/
#https://digibio.blogspot.com/2017/11/rna-seq-analysis-hisat2-featurecounts.html
#featureCounts -T 10 -t exon -g gene_id -a $genes_annot_gtf -o counts.featureCounts.txt trimmed_pe.sorted.bam;

done
wait

