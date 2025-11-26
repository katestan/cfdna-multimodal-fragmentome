#!/bin/bash

# Set working directory
WORKDIR="/path/to/your/working/directory"
cd "$WORKDIR"

# Read in tab delimited file with GC-codes listed
samples=$1

Rscript --vanilla src/celloforigin_correlations.R samplelist ref/tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_broad_cell_types_vtannotated_normbatchcorrected_expression_data.txt



