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
library(grid)
library(cowplot)
library(ggpubr)

conditions <- c(
  "DCDA",
  "MCDA",
  "CMV",
  "SAD",
  "TAD",
  "Crohns",
  "Psoriasis",
  "PretermBirth_32"
)

anno <- fread("../ref/broadcelltypes142_annotations")

dat = readRDS("models/highrisk_figure_2g_4d_figure_S2_unimodal_models.rds")

direction <- fread("../ref/directionality_allphenos.txt")
colnames(direction) <- c("foldchange", "Feature", "phenotype")

plots=list()
legends=list()
datalist=list()
# Loop through each dataset and print the modified file name and the plot
lapply(seq_along(dat), function(i) {
    # Get the original file name
    clean_name <- names(dat)[i]
    
    #Number of total selected features
    nsf=nrow(dat[[i]]$ordered_features)

    # Create the frequency table as a string
    freq_table <- as.data.frame(table(dat[[i]]$model_rf$pred$obs))
    freq_text <- paste(apply(freq_table, 1, function(x) paste(x, collapse = ": ")), collapse = ", ")
    
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features[order(-dat[[i]]$ordered_features$Importance.MeanDecreaseGini),]

    
    num_rows <- if (grepl("tissue", names(dat)[i])) { 
	    min(nrow(dat[[i]]$ordered_features), 50)
    } else {
	    min(nrow(dat[[i]]$ordered_features), 50)
    }
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features[1:num_rows, ]
    
    print(clean_name)
    print(nrow(dat[[i]]$ordered_features))
    print(rownames(dat[[i]]$ordered_features))
    print(num_rows)

    fill_color <- if (grepl("tissue", names(dat)[i])) {
  viridis(20)
} else if (grepl("endmotif", names(dat)[i])) {
  colorRampPalette(brewer.pal(9, "Greens")[4:8])(4)
  #viridis(4)
} else if (grepl("fragment", names(dat)[i])) {
  colorRampPalette(brewer.pal(9, "Purples")[4:7])(5)
  #viridis(3)
} else {
  "#404788FF"  # Default color if none of the conditions match
}

    clean_name <- gsub("tissuecelltypes204|broadcelltypes140|endmotifs256|fragmentlength600|fragmentlength", "", names(dat)[i])
    
    tissue_types <- c("adipose.tissue", "bladder.organ", "blood", "bone.marrow", 
                  "endocrine.gland", "exocrine.gland", "heart", "large.intestine",
                  "liver", "lymph.node", "pancreas", "musculature", "placenta", 
                  "respiratory.system", "skin.of.body", "small.intestine", 
                  "uterus", "vasculature", "lung", "kidney")

    insert_types <- c("ultrashort (<100bp)", "sub-nucleosomal", "mono-nucleosomal", "di-nucleosomal", "long (>450bp)")
    endmotif_types <- c("A", "T", "G", "C")
    
    if (grepl("tissue", names(dat)[i])) {
    direction_condition <- subset(direction, phenotype == clean_name)
    dat[[i]]$ordered_features <- left_join(dat[[i]]$ordered_features, direction_condition)
    dat[[i]]$ordered_features$broad <- gsub("adipose.tissue.", "", dat[[i]]$ordered_features$Feature)
    dat[[i]]$ordered_features$broad <- gsub("bronchial.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("cardiac.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("vein.endothelial.cell", "vascular.endothelial.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("vessel.endothelial.cell", "vascular.endothelial.cell", dat[[i]]$ordered_features$broad) 
    dat[[i]]$ordered_features$broad <- gsub("endothelial.cell.of.artery", "vascular.endothelial.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("capillary.endothelial.cell", "vascular.endothelial.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("endothelial.cell.of.vascular.tree", "vascular.endothelial.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("adipose.tissue.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("of.tissue", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("adipose.tissue.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("placenta.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("VCT_p", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("associated.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("VCT_CCC", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("VCT_fusing", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("VCT", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("SCT", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("EVT_1", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("EVT_2", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("iEVT", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("eEVT", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("GC", "trophoblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("bladder.organ.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.mammary.gland", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("lung.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("lung.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("gut.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.salivary.gland", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.small.intestine", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.cardiac.tissue", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.trachea", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("pulmonary.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.breast", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.hepatic.sinusoid", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("kidney.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("blood.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("bone.marrow.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("endocrine.gland.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("exocrine.gland.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("heart.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("large.intestine.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.large.intestine.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("liver.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("lymph.node.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("pancreas.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("musculature.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("placenta.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("respiratory.system.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("respiratory.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("skin.of.body.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("small.intestine.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("uterus.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("vasculature.", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("positive..alpha.beta", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.uterus", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.colon", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("neutrophil", "granulocyte", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("enterocyte", "intestinal.absorptive.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("intestinal.enterendocrine.cell", "intestinal.endocrine.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("intestinal.tuft.cell", "intestinal.secretory.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("enterocyte.of.epithelium.of.large.intestine", "intestinal.absorptive.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("enterocyte.of.epithelium", "intestinal.absorptive.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("intestinal.enteroendocrine.cell", "intestinal.endocrine.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("goblet.cell", "intestinal.secretory.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("paneth.cell.of.epithelium", "intestinal.secretory.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("paneth.cell", "intestinal.secretory.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("transit.amplifying.cell", "intestinal.absorptive.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("transit.amplifying.cell.of.colon", "intestinal.absorptive.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("duodenum.glandular.cell", "intestinal.secretory.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("intestinal.crypt.stem.cell.of.large.intestine", "intestinal.stem.cell", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub(".of.bronchus", "", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features$broad <- gsub("fibroblast.", "fibroblast", dat[[i]]$ordered_features$broad)
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  rowwise() %>%
  mutate(tissue = str_extract(Feature, paste(tissue_types, collapse = "|"))) %>%
  ungroup()

    } else if (grepl("fragment", names(dat)[i])) {
    direction_condition <- subset(direction, phenotype == clean_name)
    dat[[i]]$ordered_features <- left_join(dat[[i]]$ordered_features, direction_condition)
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  mutate(
    # Remove "insert_size_" from Feature and convert it to numeric
    Feature_numeric = as.numeric(str_replace(Feature, "insert_size_", "")),
    
    # Create the 'broad' category based on Feature_numeric value
    broad = case_when(
      Feature_numeric < 100 ~ "ultrashort (<100bp)",
      Feature_numeric >= 100 & Feature_numeric < 160 ~ "sub-nucleosomal",
      Feature_numeric >= 160 & Feature_numeric < 250 ~ "mono-nucleosomal",
      Feature_numeric >= 250 & Feature_numeric < 450 ~ "di-nucleosomal",
      Feature_numeric >= 450 ~ "long (>450bp)",
      TRUE ~ NA_character_  # Catch any NA values if they exist
    )
  )
  dat[[i]]$ordered_features$tissue <- dat[[i]]$ordered_features$broad
} else if (grepl("endmotif", names(dat)[i])) {
  direction_condition <- subset(direction, phenotype == clean_name)
  dat[[i]]$ordered_features <- left_join(dat[[i]]$ordered_features, direction_condition)
  dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  mutate(broad=substr(Feature, 10, 11))
  dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  mutate(tissue=substr(Feature, 10, 10))
} else {
  dat[[i]]$ordered_features$broad <- NA  # Default color if none of the conditions match
}
    dat[[i]]$ordered_features$broad <- as.factor(dat[[i]]$ordered_features$broad)
        dat[[i]]$ordered_features$tissue <- as.factor(dat[[i]]$ordered_features$tissue)
    dat[[i]]$ordered_features$Feature <- gsub("small.intestine.", "", dat[[i]]$ordered_features$Feature)

    if (grepl("tissue", names(dat)[i])) {
    #tissue_color_map <- setNames(fill_color, tissue_types)
    tissue_color_map <- c(
  "adipose.tissue" = viridis(20)[1],
  "bladder.organ" = viridis(20)[2],
  "blood" = viridis(20)[3],
  "bone.marrow" = viridis(20)[4],
  "endocrine.gland" = viridis(20)[5],
  "exocrine.gland" = viridis(20)[6],
  "heart" = viridis(20)[7],
  "large.intestine" = viridis(20)[8],
  "liver" = viridis(20)[9],
  "lymph.node" = viridis(20)[10],
  "pancreas" = viridis(20)[11],
  "musculature" = viridis(20)[12],
  "placenta" = viridis(20)[13],
  "respiratory.system" = viridis(20)[14],
  "skin.of.body" = viridis(20)[15],
  "small.intestine" = viridis(20)[16],
  "uterus" = viridis(20)[17],
  "vasculature" = viridis(20)[18],
  "lung" = viridis(20)[19],
  "kidney" = viridis(20)[20]
)
    #dat[[i]]$ordered_features$Importance.MeanDecreaseGini <- abs(dat[[i]]$ordered_features$Importance.MeanDecreaseGini)
    dat[[i]]$ordered_features$Importance.MeanDecreaseGini[dat[[i]]$ordered_features$foldchange < 1] <- dat[[i]]$ordered_features$Importance.MeanDecreaseGini[dat[[i]]$ordered_features$foldchange < 1] * -1
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  group_by(broad) %>%
  mutate(total_importance = sum(Importance.MeanDecreaseGini)) %>%
  ungroup()
  lim=max(abs(dat[[i]]$ordered_features$total_importance))+2
    plot_data <- dat[[i]]$ordered_features %>% select(broad, Feature, tissue, Importance.MeanDecreaseGini, total_importance) %>% arrange(desc(total_importance))
    plot_data$model <- names(dat)[i]
    plot_data <- plot_data %>% select(model, broad, Feature, tissue, Importance.MeanDecreaseGini, total_importance)
    feature_importance_plot=ggplot(dat[[i]]$ordered_features, aes(x = reorder(broad, total_importance), y = Importance.MeanDecreaseGini, group=Feature, fill=tissue)) +
      geom_bar(position="stack", stat = "identity", color="black") + scale_fill_manual(values=tissue_color_map)+
      labs(x = "", y = "") + coord_flip() +
      ggtitle(clean_name) +
      theme_light()+theme(axis.text.y=element_text(size=24), axis.text.x=element_text(size=34), legend.position = "none", plot.title = element_text(hjust = 0.5, size=34))+
      ylim(-lim, lim)+
    #annotate("text", x = -Inf, y = Inf, label = clean_name,
             #hjust = 1, vjust = -2, size = 10, color = "black")+
    annotate("text", x = -Inf, y = Inf, label = paste0("N.total = : ", nsf),
             hjust = 1, vjust = -0.5, size = 10, color = "black")
    #annotate("text", x = -Inf, y = Inf, label = freq_text,
             #hjust = 1, vjust = -2, size = 10, color = "black")  # Add condition text
    } else if (grepl("endmotif", names(dat)[i])){
    tissue_color_map <- setNames(fill_color, endmotif_types)
    dat[[i]]$ordered_features$Importance.MeanDecreaseGini[dat[[i]]$ordered_features$foldchange < 1] <- dat[[i]]$ordered_features$Importance.MeanDecreaseGini[dat[[i]]$ordered_features$foldchange < 1] * -1
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  group_by(broad) %>%
  mutate(total_importance = sum(Importance.MeanDecreaseGini)) %>%
  ungroup()
    lim=max(abs(dat[[i]]$ordered_features$total_importance))+14
    plot_data <- dat[[i]]$ordered_features %>% select(broad, total_importance, Importance.MeanDecreaseGini, Feature, tissue) %>% arrange(desc(total_importance))
    plot_data$model <- names(dat)[i]
    plot_data <- plot_data %>% select(model, broad, Feature, tissue, Importance.MeanDecreaseGini, total_importance)
    feature_importance_plot=ggplot(dat[[i]]$ordered_features, aes(x = reorder(broad, total_importance), y = Importance.MeanDecreaseGini, group=Feature, fill=tissue)) +
      geom_bar(position="stack", stat = "identity", color="black") + scale_fill_manual(values=tissue_color_map)+
      labs(x = "", y = "") + coord_flip() +
      theme_light()+theme(axis.text.y=element_text(size=24), axis.text.x=element_text(size=34), legend.position = "none", plot.title = element_text(hjust = 0.5, size=34))+
      ylim(-lim, lim)+
    #annotate("text", x = -Inf, y = Inf, label = clean_name,
             #hjust = 1, vjust = -2, size = 5, color = "black")+
    annotate("text", x = -Inf, y = Inf, label = paste0("N.total: ", nsf),
             hjust = 1, vjust = -0.5, size = 10, color = "black")
    #annotate("text", x = -Inf, y = Inf, label = freq_text,
             #hjust = 1, vjust = -2, size = 10, color = "black")  # Add condition text
    } else if (grepl("fragment", names(dat)[i])){
    tissue_color_map <- setNames(fill_color, insert_types)
    dat[[i]]$ordered_features$Importance.MeanDecreaseGini[dat[[i]]$ordered_features$foldchange < 1] <- dat[[i]]$ordered_features$Importance.MeanDecreaseGini[dat[[i]]$ordered_features$foldchange < 1] * -1
    dat[[i]]$ordered_features <- dat[[i]]$ordered_features %>%
  group_by(broad) %>%
  mutate(total_importance = sum(Importance.MeanDecreaseGini)) %>%
  ungroup()
    lim=max(abs(dat[[i]]$ordered_features$total_importance))+14
    plot_data <- dat[[i]]$ordered_features %>% select(broad, total_importance, Importance.MeanDecreaseGini, Feature, tissue) %>% arrange(desc(total_importance))
    plot_data$model <- names(dat)[i]
    plot_data <- plot_data %>% select(model, broad, Feature, tissue, Importance.MeanDecreaseGini, total_importance)
    feature_importance_plot=ggplot(dat[[i]]$ordered_features, aes(x = reorder(broad, total_importance), y = Importance.MeanDecreaseGini, group=Feature, fill=tissue)) +
      geom_bar(position="stack", stat = "identity", , color="black") + scale_fill_manual(values=tissue_color_map)+
      labs(x = "", y = "") + coord_flip() +
      theme_light()+theme(axis.text.y=element_text(size=24), axis.text.x=element_text(size=34), legend.position = "none", plot.title = element_text(hjust = 0.5, size=34))+
      #ylim(-0.002, 0.016)+
      ylim(-lim, lim)+
    #annotate("text", x = -Inf, y = Inf, label = clean_name,
             #hjust = 1, vjust = -2, size = 5, color = "black")+
    annotate("text", x = -Inf, y = Inf, label = paste0("N.total: ", nsf),
             hjust = 1, vjust = -0.5, size = 10, color = "black")
    #annotate("text", x = -Inf, y = Inf, label = freq_text,
             #hjust = 1, vjust = -2, size = 10, color = "black")  # Add condition text
    }
    # Print the combined plot to the PDF
    if (grepl("tissue", names(dat)[i])) {
    legend.plot <- ggplot(dat[[i]]$ordered_features, aes(x = reorder(broad, total_importance), y = Importance.MeanDecreaseGini, group=Feature, fill=tissue)) +
      geom_bar(position="stack", stat = "identity") + scale_fill_manual(values=tissue_color_map)+
      labs(x = "", y = "") + coord_flip() +
      theme_light()+theme(axis.text.y=element_text(size=14), legend.position = "right")
    } else {
    legend.plot <- ggplot(dat[[i]]$ordered_features, aes(x = reorder(broad, total_importance), y = Importance.MeanDecreaseGini, group=Feature, fill=tissue)) +
      geom_bar(position="stack", stat = "identity") + scale_fill_manual(values=fill_color)+
      labs(x = "", y = "") + coord_flip() +
      theme_light()+theme(axis.text.y=element_text(size=14), legend.position = "right")
    }
    legends[[i]] <<- cowplot::get_legend(legend.plot)
    plots[[i]] <<- feature_importance_plot
    datalist[[i]] <<- plot_data

})

print(warnings())

sub <- plots[c(1,2,3,4,9,10,11,12,17,18,19,20, 5,6,7,8,13,14,15,16,21,22,23,24)]

out <- plot_grid(plotlist = sub, ncol = 4, align = "hv", rel_heights=c(3,1,2,3,1,2))

out_data <- rbindlist(datalist, use.names = TRUE)

write.table(out_data, "plot_RFE_celloforigin.txt", sep="\t", row.names=F, col.names=T, quote=F)

pdf("plot_RFE_celloforigin.pdf", width=46, height=38)
print(out)
dev.off()

