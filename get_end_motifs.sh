#!/bin/bash

#Explanation:

#samtools view -f 2 -F 284 -q 1 input.bam: This command filters the BAM file to keep only the properly paired reads, with both mates mapped to the same chromosome, with mapping quality at least 1, and excludes secondary and supplementary alignments. The -F 284 flag is a combination of the flags for unmapped reads (4), mates unmapped (8), not primary alignement (256), read reverse strand (16) totaling 284.

#awk '$9<=600 && $9>=-600': This awk command further filters the alignments to keep only those with an insert size less than or equal to 600 bases.

# awk && $6 !~ /S/ : This removes soft clipped reads, as we don't want adapter sequences in the SEQ field to influence the 4-mer frequencies 

#cut -f 10: The cut command is used to extract the 10th field (column) of the SAM file, which is the SEQ (sequence) column.

#cut -c 1-4: This cut command trims the sequence to keep only the first four nucleotides (4-mer).

#sort: This sorts the list of 4-mers.

#uniq -c: This command counts the unique 4-mers in the sorted list.

#The output of this command will be a list of unique 4-mers (from the first four nucleotides of the SEQ column of the properly paired forward alignments with insert size ≤ 600), along with their counts in the input BAM file.


samples=$1

for file in $(cat $samples)
do
/data/leuven/338/vsc33813/miniconda3/bin/samtools view -f 2 -F 284 -q 1 data/fastq_ZOL/$file'_trimmed_aligned_merged_marked_duplicates.bam' | awk '$9<=600 && $9>=-600 && $6 !~ /S/' | cut -f 10 | cut -c 1-4 | sort | uniq -c > results/end_motifs/$file'_end_motif_freq.txt'
done


