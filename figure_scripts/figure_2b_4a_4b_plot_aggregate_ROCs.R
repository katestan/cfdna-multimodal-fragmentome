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

conditions=c(
	     "DCDA",
	     "MCDA",
	     "CMV",
	     "Diabetes",
	     "SAD",
	     "TAD",
	     "Crohns",
	     "Psoriasis"
	     )

#read in models
dat = readRDS("figure_scripts/models/highrisk_models.rds")

#remove APO prediction models from list
dat = dat[!grepl("APO", names(dat))]

# Initialize an empty list to store plots
plots <- list()
datalist <- list()

# Loop through each condition
for (i in 1:length(conditions)) {
    # Get the original file name
    condition <- conditions[i]
    print(condition)

    model <- dat[grepl(condition, names(dat))]
    print(length(model))
    
    # Remove unwanted parts from the file name
    print(names(model))

    #train ensembl model (stacking)
    probs.train <- model[[1]]$model_rf$pred
    probs.train <- probs.train[, c("SampleID", "case")]
    probs.train$model <- "Discovery"
    probs.train$target <- as.factor(sub(".*_", "", probs.train$SampleID))

    #Get Sample ID index
    map <- as.data.frame(probs.train$SampleID)
    colnames(map) <- c("SampleID")
    map$rowIndex <- rep(1:nrow(map))

    probs.test <- model[[1]]$predictions
    probs.test$SampleID <- rownames(probs.test)
    probs.test <- probs.test[, c("SampleID", "case")]
    probs.test$model <- "Validation"
    probs.test$target <- as.factor(sub(".*_", "", probs.test$SampleID))


    roc_data_models <- data.frame()  # Initialize as a dataframe

    probs.train$target <- factor(probs.train$target, levels = c("control", "case"))
    roc_obj <- roc(probs.train$target, probs.train$case, direction = "<", levels = c("control", "case"))
    auc_value <- auc(roc_obj)
    auc_ci <- ci.auc(roc_obj)
    sensitivity_at_5_fpr <- coords(roc_obj, x = 0.95, input = "specificity", ret = "sensitivity")

        roc_data <- data.frame(
            TPR = rev(roc_obj$sensitivities),
            FPR = rev(1 - roc_obj$specificities),
            model = "Discovery",
            auc = round(auc_value, 3),
	    auc_ci = paste0("(", round(auc_ci[1], 3), "-", round(auc_ci[3],3), ")"),
            sensitivity_at_5_fpr = sensitivity_at_5_fpr$sensitivity
        )

    roc_data_models <- rbind(roc_data_models, roc_data)

    # Create the ROC data for ensemble test set
    probs.test$target <- factor(probs.test$target, levels = c("control", "case"))
    roc_obj <- roc(probs.test$target, probs.test$case, direction = "<", levels = c("control", "case"))
    auc_value <- auc(roc_obj)
    auc_ci <- ci.auc(roc_obj)
    sensitivity_at_5_fpr <- coords(roc_obj, x = 0.95, input = "specificity", ret = "sensitivity")

        roc_data <- data.frame(
            TPR = rev(roc_obj$sensitivities),
            FPR = rev(1 - roc_obj$specificities),
            model = "Validation",
            auc = round(auc_value, 3),
	    auc_ci = paste0("(", round(auc_ci[1], 3), "-", round(auc_ci[3],3), ")"),
            sensitivity_at_5_fpr = sensitivity_at_5_fpr$sensitivity
        )

    roc_data_models <- rbind(roc_data_models, roc_data)


    # Frequency table and text for annotations
    freq_table_train_case <- as.data.frame(table(probs.train$target))[2,2]
    freq_table_train_ctrl <- as.data.frame(table(probs.train$target))[1,2]
    
    freq_table_test_case <- as.data.frame(table(probs.test$target))[2,2]
    freq_table_test_ctrl <- as.data.frame(table(probs.test$target))[1,2]

    # Unique AUC values for each model
    perf <- roc_data_models %>% select(model, auc, auc_ci, sensitivity_at_5_fpr) %>% unique()

    # Create the plot
    plot <- ggplot(roc_data_models, aes(x = FPR, y = TPR, color = model)) +
        geom_line(linewidth = 2.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
        xlab("False Positive Rate") +
        ylab("True Positive Rate") +
        scale_color_manual(values = c("#b9b9b9", "#354b86")) +
	#scale_color_manual(values = c("#440154FF", "#35B779FF", "#31688EFF")) +
        annotate("text", x = 1, y = 0.38, label = condition, color = "black", hjust = 1, vjust = 1.5, size = 12) +
	annotate("text", x = 1, y = 0.29, label = paste("AUC =", round(perf$auc[1], 3), perf$auc_ci[1]), color = "#b9b9b9", hjust = 1, vjust = 1.5, size = 10) +

	annotate("text", x = 1, y = 0.21, label = paste("Sens 5% FPR =", round(perf$sensitivity_at_5_fpr[1], 3)), color = "#b9b9b9", hjust = 1, vjust = 1.5, size = 10) +
	annotate("text", x = 1, y = 0.13, label = paste("AUC =", round(perf$auc[2], 3), perf$auc_ci[2]), color = "#354b86", hjust = 1, vjust = 1.5, size = 10) +
	annotate("text", x = 1, y = 0.05, label = paste("Sens 5% FPR =", round(perf$sensitivity_at_5_fpr[2], 3)), color = "#354b86", hjust = 1, vjust = 1.5, size = 10) +
        annotate("text", x = 0, y = 1, label = paste("n. cases: ", freq_table_train_case), hjust = 0, color = "#b9b9b9", size = 10, parse = FALSE) +
	annotate("text", x = 0, y = 0.94, label = paste("n. ctrls: ",freq_table_train_ctrl), hjust = 0, color = "#b9b9b9", size = 10, parse = FALSE) +
	annotate("text", x = 0, y = 0.86, label = paste("n. cases: ",freq_table_test_case), hjust = 0, color = "#354b86", size = 10, parse = FALSE) +
	annotate("text", x = 0, y = 0.80, label = paste("n. ctrls: ",freq_table_test_ctrl), hjust = 0, color = "#354b86", size = 10, parse = FALSE) +
        theme_light() + theme(legend.position = "none")

    # Store the plot in the list
    plots[[i]] <- plot
    roc_data_models$condition <- condition
    datalist[[i]] <- roc_data_models
}

pdf("ROCs_aggregated_matconditions_onevsctrls.pdf", width=30, height=16)
plot_grid(plotlist = plots, ncol = 4, align = "hv")
dev.off()

out_data <- do.call(rbind, datalist)
write.table(out_data, "ROCs_aggregated_matconditions_onevsctrls.txt", sep="\t", row.names=F, col.names=T, quote=F)

