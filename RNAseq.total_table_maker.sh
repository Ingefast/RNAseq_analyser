#Juan Santos. 2022 Oct


##################################################################################

#working directory where all the sample directories are allocated
working_dir="/media/diskc/project_RNAseq_tmp"

mkdir $working_dir/summary_output;

#this defines the samples to be analysed (dir names)
sample_list="mutant_rep1  mutant_rep2  mutant_rep3  wt_rep1  wt_rep2  wt_rep3";
#sample_list="test_rep1";

##################################################################################

####### Step 0. Makes a dir to place all the counting tables

	mkdir $working_dir/summary_output;


####### Step 1. Copying and gathering gene expression counting files together

for sample in $sample_list ; do
	
	cp $working_dir/REPLICATES_TOTAL/${sample}/abund.stringtie.txt $working_dir/summary_output/${sample}.abund.stringtie.txt;
	
	cp $working_dir/REPLICATES_TOTAL/${sample}/counts.htseq.txt $working_dir/summary_output/${sample}.counts.htseq.txt;

done
wait


####### Step 2. Copying and gathering gene expression counting files together

cd $working_dir/summary_output;

####	this makes a table with raw read count values from htseq-count

grep '^AT' mutant_rep1.counts.htseq.txt|sort| cut -f 1 > tmp_ids

grep '^AT' mutant_rep1.counts.htseq.txt|sort| cut -f 2 > tmp1
grep '^AT' mutant_rep2.counts.htseq.txt|sort| cut -f 2 > tmp2
grep '^AT' mutant_rep3.counts.htseq.txt|sort| cut -f 2 > tmp3

grep '^AT' wt_rep1.counts.htseq.txt|sort| cut -f 2 > tmp4
grep '^AT' wt_rep2.counts.htseq.txt|sort| cut -f 2 > tmp5
grep '^AT' wt_rep3.counts.htseq.txt|sort| cut -f 2 > tmp6

paste tmp_ids tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 >tmp_total;

echo -e "id\tmutant_rep1\tmutant_rep2\tmutant_rep3\twt_rep1\twt_rep2\twt_rep3">tmp_header;

cat tmp_header tmp_total>total.gene.counts.txt;

rm tmp*;


#### this makes a table with TPM values form stringtie


grep '^AT' mutant_rep1.abund.stringtie.txt|sort| cut -f 1 > tmp_ids

grep '^AT' mutant_rep1.abund.stringtie.txt|sort| cut -f 9 > tmp1
grep '^AT' mutant_rep2.abund.stringtie.txt|sort| cut -f 9 > tmp2
grep '^AT' mutant_rep3.abund.stringtie.txt|sort| cut -f 9 > tmp3

grep '^AT' wt_rep1.abund.stringtie.txt|sort| cut -f 9 > tmp4
grep '^AT' wt_rep2.abund.stringtie.txt|sort| cut -f 9 > tmp5
grep '^AT' wt_rep3.abund.stringtie.txt|sort| cut -f 9 > tmp6

paste tmp_ids tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 >tmp_total;

echo -e "id\tmutant_rep1\tmutant_rep2\tmutant_rep3\twt_rep1\twt_rep2\twt_rep3">tmp_header;

cat tmp_header tmp_total>total.gene.TPM.txt;

rm tmp*;
