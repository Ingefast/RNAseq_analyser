#Juan Santos. 2022 Oct

##################################################################################

#working directory where all the sample directories are allocated
working_dir="/home/jsantos/project_RNAseq_tmp/REPLICATES_TOTAL/summary_output"

####	this makes a table with count values form htseq-count

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

