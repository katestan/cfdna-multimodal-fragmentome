# cfDNA Multimodal Fragmentome

This repository contains scripts for analysing cfDNA sequencing data. Below are the requirements, folder structure, and instructions for running the scripts.

## Folder Structure

The repository expects the following folder structure:

```
cfdna-multimodal-fragmentome/
├── data/
│   ├── $sample_trimmed_aligned_merged_marked_duplicates.bam
│   ├── samples
│   └── metadata
├── results/
│   ├── end_motifs/
│   ├── read_length/
│   └── body/
├── ref/
│   └── tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_broad_cell_types_vtannotated_normbatchcorrected_expression_data.txt.gz
│   └── tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_tissue_general_cell_types_vtannotated_normbatchcorrected_expression_data.txt.gz
├── src/
│   ├── celloforigin_correlations.R
│   ├── convert_files.py
│   ├── fft_path.R
│   ├── extractReadStartsFromBAM_Region_WPS.py
│   ├── model_training_testing_functions.R
│   └── run_model_training_testing_unimodal.R
├── get_end_motifs.sh
├── get_read_lengths.sh
├── get_wps.sh
└── get_correlations.sh
```

### Required Files

1. **`data/` folder**:
   - BAM files named as `$sample_trimmed_aligned_merged_marked_duplicates.bam`.
   - A `samples` file containing a list of sample IDs, one per line.
   - A `metadata` file with the following columns:
     - `GC-AR_gc`: Sample ID.
     - `class`: Indicates whether the sample is `case` or `control`.
     - Columns for each condition to be predicted (e.g., `SAD` (yes/no), `Crohns` (yes/no)).
       - Controls should have `no` in all condition columns and `control` in the `class` column.

2. **`results/` folder**:
   - Subfolders `end_motifs/`, `read_length/`, and `body/` for storing script outputs.

3. **`ref/` folder**:
   - Contains the reference atlas file required for the `get_correlations.sh` script.

## Script Execution

### 1. Generate End Motifs
Run the following command to generate end motifs:
```bash
bash get_end_motifs.sh $samples
```
**Output**: `results/end_motifs/$sample_end_motif_freq.txt`

### 2. Calculate Read Lengths
Run the following command to calculate read lengths:
```bash
bash get_read_lengths.sh $samples
```
**Output**: `results/read_length/$sample.picard_insert_size_metrics.txt`

### 3. Compute WPS (Windowed Protection Score)
Run the following command to compute WPS:
```bash
bash get_wps.sh $samples
```
**Output**:
- `results/body/$sample/fft`
- `results/body/$sample/count`
- `results/body/$sample/fft_summaries`

### 4. Perform Correlations
Once `get_wps.sh` is complete, run the following command to perform correlations:
```bash
bash get_correlations.sh $samples ref/$reference_atlas
```
**Output**: `results/body/$sample/fft_summaries/$sample_<reference_atlas_name>_Ave193-199bp_correlation.csv`

## Notes
- The `get_correlations.sh` script depends on the output of `get_wps.sh`. Ensure `get_wps.sh` completes successfully before running `get_correlations.sh`.
- Replace `$samples` with the path to the `samples` file containing the list of sample IDs.
- Replace `$reference_atlas` with the name of the reference atlas file located in the `ref/` folder.

## Example
Assume the `samples` file contains:
```
sample1
sample2
sample3
```

To run the scripts:
```bash
bash get_end_motifs.sh data/samples
bash get_read_lengths.sh data/samples
bash get_wps.sh data/samples
bash get_correlations.sh data/samples ref/tabsap_vt2023_reference_atlas.txt
```
