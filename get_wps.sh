#!/bin/bash

# Set working directory
WORKDIR="/path/to/your/working/directory"
cd "$WORKDIR"

# Read in tab delimited file with GC-codes listed
samples=$1

for sample in $(cat $samples)
do

mkdir body/$sample
cd body/$sample
mkdir counts
cd counts
python2.7 ../../src/extractReadStartsFromBAM_Region_WPS.py --minInsert=120 --maxInsert=180 -o 'block_%s.tsv.gz' -i ../../data/Ensemble_canonical_GRCh38.body.tsv data/$sample'_trimmed_aligned_merged_marked_duplicates.bam'

cd "$WORKDIR"
mkdir body/$sample/fft
( cd body/$sample/counts; ls block_*.tsv.gz ) | xargs -n 500 Rscript src/fft_path.R body/$sample/counts body/$sample/fft

mkdir body/$sample/fft_summaries
python2.7 src/convert_files.py -a ref/Ensemble_canonical_GRCh38.body.tsv -t "$WORKDIR" -r "$WORKDIR" -p body -i $sample

done
