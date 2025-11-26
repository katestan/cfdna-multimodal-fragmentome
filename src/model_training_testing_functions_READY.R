##function to calculate foldchange and wilcox p-value for ranks of each cell type between cases and controls
getDirectionality <- function(dat.anno, prefix) {
  alltruegroup <- c()
  cell <- unique(colnames(dat.anno))
  cell <- setdiff(cell, "class")
  for(i in 1:length(cell)) {
    value=cell[i]	  
    df.mean <- dat.anno %>% group_by(class) %>% summarize(Mean = mean(!!sym(value), na.rm=TRUE)) %>% as.data.frame()
    print(df.mean)
    if (value %like% "endmotif" | value %like% "insert_size"){
    df.mean.FC <- df.mean$Mean[1]/df.mean$Mean[2] %>% as.data.frame()
    } else {
    df.mean.FC <- df.mean$Mean[2]/df.mean$Mean[1] %>% as.data.frame()
    }
    print(paste0("feature: ", value))
    print(paste0("foldchange control: ", df.mean[2,1]))
    colnames(df.mean.FC) <- c("foldchange")
    df.mean.FC$cell_type <- value
    alltruegroup <- rbind(alltruegroup, df.mean.FC)
  }
   alltruegroup$phenotype <- prefix
   return(alltruegroup)
}

##function to perform interative feature selection using random forest, and then tune models and apply best performing model on the test set
tuneMLmodels <- function(df, rfe=TRUE, rfeModel = "rf", feature, split){
# Split the data into training (80%) and testing (20%) sets with stratified sampling
set.seed(200)
df$class <- factor(df$class, levels = c("case", "control"))
trainIndex <- createDataPartition(df$class, p = split, list = FALSE)
trainData <- df[trainIndex, ]
print(levels(trainData$class))
testData <- df[-trainIndex, ]
print(levels(testData$class))

classDistribution <- table(trainData$class)
classWeights <- max(classDistribution) / classDistribution

# Map class weights to the correct classes
weightVector <- ifelse(trainData$class == names(classWeights)[1], classWeights[1], classWeights[2])

map <- as.data.frame(rownames(trainData))
colnames(map) <- c("SampleID")
map$rowIndex <- rep(1:nrow(map))

# Define control function for training with 10-fold cross-validation
control <- trainControl(method = "LOOCV", classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")

# Define hyperparameter grids for each model
svmGrid <- expand.grid(C = 2^(-2:2), sigma = 2^(-2:2))
svmGridLinear <- expand.grid(cost = 2^(-2:2), weight=c(1,5))
rfGrid <- expand.grid(mtry = 2:4)
enetGrid <- expand.grid(alpha = seq(0, 1, length = 10), lambda = seq(0.001, 0.1, length = 10))
rpartGrid <- expand.grid(cp = seq(0.001, 0.1, length = 10))

#recursive elimination of features
    # Feature selection using RFE
    set.seed(200)

    if (rfeModel == "rf") {
      rfeControl <- rfeControl(functions = rfFuncs, method = "LOOCV")
    } else if (rfeModel == "lm") {
      rfeControl <- rfeControl(functions = lmFuncs, method = "LOOCV")
    } else if (rfeModel == "treebag") {
      rfeControl <- rfeControl(functions = treebagFuncs, method = "LOOCV")
    }

    if (feature == "ranks") {
    rfeResults <- rfe(trainData[, -ncol(trainData)], trainData$class, sizes = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 140), rfeControl = rfeControl)
    print(rfeResults)
    } else if (feature == "endmotifs") {
    rfeResults <- rfe(trainData[, -ncol(trainData)], trainData$class, sizes = c(10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 190, 220, 250, 256), rfeControl = rfeControl)
    print(rfeResults)
    } else if (feature == "fragmentlength") {
    rfeResults <- rfe(trainData[, -ncol(trainData)], trainData$class, sizes = c(10, 25, 50, 75, 100, 200, 300, 400, 500, 600), rfeControl = rfeControl)
    } else if (feature == "ensemble") {
    print(ncol(trainData)-1)
    rfeResults <- rfe(trainData[, -ncol(trainData)], trainData$class, sizes = seq(10, ncol(trainData)-1, by = 40), rfeControl = rfeControl)
    } else if (feature == "tissueranks") {
    print(ncol(trainData)-1)
    rfeResults <- rfe(trainData[, -ncol(trainData)], trainData$class, sizes = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200), rfeControl = rfeControl)
    }

    # Get the selected features
    selectedFeatures <- predictors(rfeResults)
    print(selectedFeatures)
    trainData <- trainData[, c(selectedFeatures, "class")]
    testData <- testData[, c(selectedFeatures, "class")]

    feature_data <- data.frame(
      Feature = selectedFeatures,
      Importance = rfeResults$fit$importance[selectedFeatures,]
    )
    ordered_features <- feature_data[order(rfeResults$optVariables), ]

    # Plotting the barplot
    feature_importance_plot=ggplot(ordered_features, aes(x = reorder(Feature, Importance.MeanDecreaseAccuracy), y = Importance.MeanDecreaseAccuracy)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(x = "Feature", y = "Importance") +
      ggtitle("Feature Importance Selected by RFE") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Train models with hyperparameter tuning
set.seed(200)
model_rf <- train(class ~ ., data = trainData, method = "rf", metric = "ROC", tuneGrid = rfGrid, trControl = control, weights = weightVector)
saveRDS(model_rf, "master_detailedphenotypes_loocv_onevsctrls_selectedFeatures_model.RDS")
model_rf$pred <- left_join(model_rf$pred, map)

# Collect results
print(model_rf)
best_auc_rf <- max(model_rf$results$ROC)
print(cat("Best AUC for Random Forest:", best_auc_rf, "\n"))

#Apply models on test set
predictions = predict(model_rf, newdata = testData, type = "prob")

  # Calculate ROC for the best Random Forest model using cross-validation predictions

  cv_predictions <- model_rf$pred
  cv_predictions$obs <- factor(cv_predictions$obs, levels=c("control", "case")) 
  roc_obj <- roc(cv_predictions$obs, cv_predictions$case, direction = "<", levels = c("control", "case"))
  auc_value <- auc(roc_obj)
  coords_optimal_f1 <- coords(roc_obj, "best", best.method = "closest.topleft", ret = c("sensitivity", "specificity"))

  # Create a data frame for ggplot
  roc_data <- data.frame(
    TPR = rev(roc_obj$sensitivities),
    FPR = rev(1 - roc_obj$specificities)
  )

  # Plot the ROC curve with AUC, sensitivity, and specificity
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggtitle("ROC Curve for Best Random Forest Model (CV Predictions)") +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_value, 2))) +
    annotate("text", x = 0.8, y = 0.15, label = paste("Sensitivity =", round(coords_optimal_f1$sensitivity, 2))) +
    annotate("text", x = 0.8, y = 0.1, label = paste("Specificity =", round(coords_optimal_f1$specificity, 2)))

return(list(model_rf = model_rf,
    ordered_features = ordered_features,
    predictions = predictions,
    roc_plot = roc_plot,
    feature_importance_plot = feature_importance_plot
    ))
}

finalMLmodels <- function(df, rfe=TRUE, rfeModel = "rf", feature){

map <- as.data.frame(rownames(df))
colnames(map) <- c("SampleID")
map$rowIndex <- rep(1:nrow(map))

# Define control function for training with 10-fold cross-validation
control <- trainControl(method = "LOOCV", classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")

# Define hyperparameter grids for each model
rfGrid <- expand.grid(mtry = 2:4)

#recursive elimination of features
    # Feature selection using RFE
    set.seed(200)

    if (rfeModel == "rf") {
      rfeControl <- rfeControl(functions = rfFuncs, method = "LOOCV")
    } else if (rfeModel == "lm") {
      rfeControl <- rfeControl(functions = lmFuncs, method = "LOOCV")
    } else if (rfeModel == "treebag") {
      rfeControl <- rfeControl(functions = treebagFuncs, method = "LOOCV")
    }

    if (feature == "ranks") {
    rfeResults <- rfe(df[, -ncol(df)], df$class, sizes = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380), rfeControl = rfeControl)
print(rfeResults)
    } else if (feature == "endmotifs") {
    rfeResults <- rfe(df[, -ncol(df)], df$class, sizes = c(10, 30, 50, 70, 90, 110, 130, 150, 170, 190, 190, 220, 250, 256), rfeControl = rfeControl)
    print(rfeResults)
    } else if (feature == "fragmentlength") {
    rfeResults <- rfe(df[, -ncol(df)], df$class, sizes = c(10, 25, 50, 75, 100, 200, 300, 400, 500, 600), rfeControl = rfeControl)
    } else if (feature == "tissueranks") {
    rfeResults <- rfe(df[, -ncol(df)], df$class, sizes = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200), rfeControl = rfeControl)
    } else if (feature == "ensemble") {
    print(ncol(df)-1)
    rfeResults <- rfe(df[, -ncol(df)], df$class, sizes = seq(5, ncol(df)-1, by = 10), rfeControl = rfeControl)
    }

    # Get the selected features
    selectedFeatures <- predictors(rfeResults)
    print(selectedFeatures)
    df <- df[, c(selectedFeatures, "class")]

    feature_data <- data.frame(
      Feature = selectedFeatures,
      Importance = rfeResults$fit$importance[selectedFeatures,]
    )
    ordered_features <- feature_data[order(rfeResults$optVariables), ]

    # Plotting the barplot
    feature_importance_plot=ggplot(ordered_features, aes(x = reorder(Feature, Importance.MeanDecreaseAccuracy), y = Importance.MeanDecreaseAccuracy)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(x = "Feature", y = "Importance") +
      ggtitle("Feature Importance Selected by RFE") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Train models with hyperparameter tuning
set.seed(200)
model_rf <- train(class ~ ., data = df, method = "rf", metric = "ROC", tuneGrid = rfGrid, trControl = control)
model_rf$pred <- left_join(model_rf$pred, map)

print(model_rf)

# Extract the best model AUCs
best_auc_rf <- max(model_rf$results$ROC)

# Print the best AUCs
print(cat("Best AUC for Random Forest:", best_auc_rf, "\n"))

#Calculate ROC for the best Random Forest model using cross-validation predictions

  cv_predictions <- model_rf$pred
  cv_predictions$obs <- factor(cv_predictions$obs, levels=c("control", "case"))
  roc_obj <- roc(cv_predictions$obs, cv_predictions$case, direction = "<", levels = c("control", "case"))
  auc_value <- auc(roc_obj)
  coords_optimal_f1 <- coords(roc_obj, "best", best.method = "closest.topleft", ret = c("sensitivity", "specificity"))

  # Create a data frame for ggplot
  roc_data <- data.frame(
    TPR = rev(roc_obj$sensitivities),
    FPR = rev(1 - roc_obj$specificities)
  )

  # Plot the ROC curve with AUC, sensitivity, and specificity
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "#404788FF") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    ggtitle("ROC Curve for Final Model") +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    annotate("text", x = 0.8, y = 0.1, label = paste("AUC =", round(auc_value, 2))) +
    theme_light()

return(list(model_rf = model_rf,
    ordered_features = ordered_features,
    roc_plot = roc_plot,
    roc_data = roc_data,
    feature_importance_plot = feature_importance_plot
    ))
}}
