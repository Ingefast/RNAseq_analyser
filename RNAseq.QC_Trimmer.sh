#!/bin/bash

#Juan Santos, Nov 2022
#This script runs Quality Control of Raw Reads, Trimming and Quality Control of Trimmed Reads.

##################################################################################

#working directory where all the sample directories are allocated
working_dir="/media/diskc/project_RNAseq_tmp/REPLICATES_TOTAL"

#this defines the samples to be analysed (dir names)
sample_list="mutant_rep1  mutant_rep2  mutant_rep3  wt_rep1  wt_rep2  wt_rep3";
#sample_list="test_rep1";

##################################################################################


####### Step 1. Quality Control of Raw Sequences

for sample in $sample_list ; do

	cd $working_dir/${sample};
	
	mkdir QC_R1;
	mkdir QC_R2;

	fastqc *R1.fastq -o QC_R1&
	fastqc *R2.fastq -o QC_R2&

done
wait



####### Step 2. Quailty trimming of sequences

for sample in $sample_list ; do

	cd $working_dir/${sample};

	trim_galore --paired  --fastqc --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 20 --three_prime_clip_R2 20 *R1.fastq *R2.fastq&

done
wait



####### Step 3. Renaming

for sample in $sample_list ; do

	cd $working_dir/${sample};

	mv *_val_1.fq trimmed_1.fastq;
	mv *_val_2.fq trimmed_2.fastq;

done
wait




