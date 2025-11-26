#!/bin/bash

samples=$1

for file in $(cat $samples)
do
/data/leuven/338/vsc33813/miniconda3/bin/picard CollectInsertSizeMetrics I=$file'_trimmed_aligned_merged_marked_duplicates.bam' O=/lustre1/project/stg_00019/research/Kate/fragmentome/results/read_length/$file'.picard_insert_size_metrics.txt' H=/lustre1/project/stg_00019/research/Kate/fragmentome/results/read_length/$file'_histogram.pdf' M=0.00000000000000000000000000000000000000000000000000000000005 DEVIATIONS=40
done

