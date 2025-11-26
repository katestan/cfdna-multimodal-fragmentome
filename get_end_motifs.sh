#!/bin/bash

#Explanation:

#samtools view -f 2 -F 268 -q 1 input.bam: This command filters the BAM file to keep only the properly paired reads, with both mates mapped to the same chromosome, with mapping quality at least 1, and excludes secondary and supplementary alignments. The -F 268 flag is a combination of the flags for unmapped reads (4), mates unmapped (8), read is PCR or optical duplicate (256), totaling 268.

#awk '$9<=600 && $9>=-600': This awk command further filters the alignments to keep only those with an insert size less than or equal to 600 bases.

# awk && $6 !~ /S/ : This removes soft clipped reads, as we don't want adapter sequences in the SEQ field to influence the 4-mer frequencies 

#cut -f 10: The cut command is used to extract the 10th field (column) of the SAM file, which is the SEQ (sequence) column.

#cut -c 1-4: This cut command trims the sequence to keep only the first four nucleotides (4-mer).

#sort: This sorts the list of 4-mers.

#uniq -c: This command counts the unique 4-mers in the sorted list.

#The output of this command will be a list of unique 4-mers (from the first four nucleotides of the SEQ column of the properly paired alignments with insert size ≤ 600), along with their counts in the input BAM file.

#The SEQ field in the bam field is always represented in the 5' to 3' direction (note that this means that the SEQ field for the reverse strand will be the reverse complement of the original read). For our purposes this is easier, as we need 4-mers from the 5' end only, so first 4 nuceotides of SEQ field. 

samples=$1

for file in $(cat $samples)
do
/data/leuven/338/vsc33813/miniconda3/bin/samtools view -f 2 -F 268 -q 1 data/fastq_ZOL/$file'_trimmed_aligned_merged_marked_duplicates.bam' | awk '$9<=600 && $9>=-600 && $6 !~ /S/' | cut -f 10 | cut -c 1-4 | sort | uniq -c > results/end_motifs/$file'_end_motif_freq.txt'
done


