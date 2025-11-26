library(data.table)
library(rstatix)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggforce)
library(tidytext)
library(readxl)
library(viridis)
library(data.table)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(igraph)
library(FNN)
library(cowplot)
library(caret)
library(pROC)
library(RColorBrewer)
source("src/model_training_testing_functions.R")
args <- commandArgs(TRUE)

#read in metadata for each sample, must have GC-AR_gc column with sample ID and class column (case/control) and a column per condition that you would like to run. For example SAD (yes/no), Crohns (yes/no). These should be case, control should be no in all condition columns and "control" in class column.
metadata <- read.table(metadata, header=T)

#read in list of all files output by get_correlation.sh, each file should look like: $sample_tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_broad_cell_types_vtannotated_normbatchcorrected_expression_Ave193-199bp_correlation.csv
files <- read.table("tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbacth_broad_cell_types_vtannotated_normbatchcorrected_expression_Ave193-199bp_correlation_files", header=F)

#read in cell rank data
dat <- lapply(files, function(x) fread(x, header=T))

# Define cell types to exclude
exclude_patterns <- c("eye", "retina", "pigmented", "salivary", "surface.ectodermal", 
                      "duct.epithelial.cell", "acinar.cell.of.salivary.gland", 
                      "mucus.secreting.cell", "secretory.cell", "serous.cell.of.epithelium.of.trachea", 
                      "tracheal.goblet.cell", "epithelial.cell.of.lacrimal.sac", "corneal", 
                      "radial", "microglial", "Mueller", "keratocyte", 
                      "conjunctival", "keratinocyte", "Schwann", "tongue", 
                      "sperm", "prostate", "unknown")

# Remove cell types that are not relevant to study
for (i in seq_along(dat)) {
  dat[[i]]$GC_code <- files[i]
  colnames(dat[[i]]) <- c("feature", "correlation", "value", "GC-AR_gc")
  
  # Combine exclusion patterns into a single regular expression
  exclude_regex <- paste(exclude_patterns, collapse = "|")
  
  # Exclude rows with cell types matching any pattern in exclude_patterns
  dat[[i]] <- dat[[i]][!grepl(exclude_regex, dat[[i]]$feature), ]
  
  # Replace "platelet" with "platelet.megakaryocyte"
  dat[[i]]$feature <- gsub("platelet", "platelet.megakaryocyte", dat[[i]]$feature)
  
  # Recreate rank column
  dat[[i]]$value <- seq_len(nrow(dat[[i]]))

  dat[[i]] <- dat[[i]] %>% select(-correlation)
}

# Combine files
dat <- do.call(rbind, dat)

# Should join by GC-AR_gc 
dat <- left_join(dat, metadata)

dat$class <- as.factor(dat$class)

ranks <- dat

files <- read.table("tabsap_vt2023trophoblasts_censusversion2024-04-15_normalized_log1p_scale_integrated_on_sequencingbatch_tissue_cell_types_normbatchcorrected_Ave193-199bp_correlation_files", header=F)

#read in cell rank data
dat <- lapply(files, function(x) fread(x, header=T))

include_celltypes <- fread("include_celltypes", header=F)
include_celltypes <- include_celltypes$V1

for (i in seq_along(dat)) {
  dat[[i]]$GC_code <- files[i]
  colnames(dat[[i]]) <- c("feature", "correlation", "value", "GC-AR_gc")

  # Combine exclusion patterns into a single regular expression
  include_regex <- paste(include_celltypes, collapse = "|")

  # Exclude rows with cell types matching any pattern in exclude_patterns
  dat[[i]] <- dat[[i]][grepl(include_regex, dat[[i]]$feature), ]

  # Replace "platelet" with "platelet.megakaryocyte"
  dat[[i]]$feature <- gsub("platelet", "platelet.megakaryocyte", dat[[i]]$feature)

  # Recreate rank column
  dat[[i]]$value <- seq_len(nrow(dat[[i]]))

  dat[[i]] <- dat[[i]] %>% select(-correlation)
}

# Combine files
dat <- do.call(rbind, dat)
print(length(unique(dat$feature)))

dat <- left_join(dat, metadata)

dat$class <- as.factor(dat$class)

tissueranks <- dat

files <- read.table("insert_size_metrics.txt", header=F)

#read in read length data
dat <- lapply(files$V1, function(x) fread(x, skip=12, nrows=600, header=T))

for (i in seq_along(dat)) {

	colnames(dat[[i]]) <- c("feature", "All_Reads.fr_count", "All_Reads.rf_count", "All_Reads.tandem_count")

	dat[[i]] <- dat[[i]] %>%
                mutate(value = All_Reads.fr_count / sum(All_Reads.fr_count))

	dat[[i]] <- dat[[i]] %>% select(feature, value)

	dat[[i]] <- subset(dat[[i]], feature <= 600)

	dat[[i]]$`GC-AR_gc` <- files$V1[i]

}

dat <- do.call(rbind, dat)

dat <- left_join(dat, metadata)

dat$feature <- paste0("insert_size_", dat$feature)

dat$class <- as.factor(dat$class)

fragmentsize <- dat

files <- read.table("end_motif_metrics.txt", header=F)

dat <- lapply(files$V1, function(x) fread(x))

for (i in seq_along(dat)) {

        colnames(dat[[i]]) <- c("value", "feature")

	dat[[i]] <- subset(dat[[i]], !feature %like% "N")

        dat[[i]] <- dat[[i]] %>%
                mutate(value = value / sum(value))

        dat[[i]]$`GC-AR_gc` <- files$V1[i]

}

dat <- do.call(rbind, dat)
print(length(unique(dat$feature)))

dat <- left_join(dat, metadata)
dat$feature <- paste0("endmotif_", dat$feature)

dat$class <- as.factor(dat$class)

dat <- dat %>% select(feature, value, everything())

endmotifs <- dat

condition=args[1]
feature=args[2] ## has to be either fragmentsize, endmotifs, tissueranks, or ranks

dat <- feature

df <- dat[dat[[condition]] == "yes" | dat$class == "control", ]

##Format data wider

df.wider <- df %>%
  filter(!is.na(feature)) %>%
  select(feature, value, `GC-AR_gc`, class) %>%
  pivot_wider(names_from = feature, values_from = value) %>%
  as.data.frame() %>%
  mutate(rownames = paste(`GC-AR_gc`, class, sep = "_")) %>%  # Create a new column for row names
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  column_to_rownames("rownames") %>% select(-`GC-AR_gc`) %>% # Set the new row names
  data.matrix() %>%
  as.data.frame()

#Feature Selection and Model Tuning/Training/Testing
df.wider <- df.wider[,c(2:ncol(df.wider), 1)]
df.wider$class <- gsub("1", "case", df.wider$class)
df.wider$class <- gsub("2", "control", df.wider$class)
df.wider$class <- factor(df.wider$class, levels=c("case", "control"))
print(levels(df.wider$class))  # Should output: "case" "control"

result <- tuneMLmodels(df.wider, rfe=TRUE, rfeModel = "rf", feature=feature, split=0.8)
saveRDS(result, paste0("RF_", feature, condition, "_RFEmarkerselection_tuneMLmodels.rds"))

output <- getDirectionality(df.wider, prefix=condition)
write.table(output, paste0(feature, "_", condition, "_directionality.txt", row.names=F, col.names=T, quote=F, sep="\t")

