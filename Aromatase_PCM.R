library(readxl)
library(prospectr)
library(Rcpi)
library(pls)
library(caret)
library(genefilter)

data <- read_excel("aromatase_data.xlsx")
### preparing data for 13 models 
compound <- data[, 5:311]
protein <- data[, 313:351]
pIC50 <- data[2]
compound = compound[, -nearZeroVar(compound)]
compound[compound == "0"] <- -1
CxP <- getCPI(compound, protein, type = "tensorprod")
CxP <- as.data.frame(CxP)
dfcompound <- names(data.frame(compound[,1:37]))
dfprotein <- names(data.frame(protein[,1:39]))
compoundNamecross <- rep(dfcompound, each = 39)
proteinNamecross <- rep(dfprotein, times = 37)
label <- paste(compoundNamecross, proteinNamecross, sep="_")
colnames(CxP) <- label
PxP <- getCPI(protein, protein, type = "tensorprod")
proteinName2 <- rep(dfprotein, times = 39)
proteinName1 <- rep(dfprotein, each = 39)
label_protein <- paste(proteinName1, proteinName2, sep = "_")
colnames(PxP) <- label_protein
index <- seq(1, 1521, by = 40)
protein_selfcross <- PxP[, -index]
transposedIndexed_protein <- t(protein_selfcross)
index1 <- which(duplicated(transposedIndexed_protein))
removed_duplicated_protein <- transposedIndexed_protein[-index1, ]
PxP <- t(removed_duplicated_protein)

CxC <- getCPI(compound, compound, type = "tensorprod")
compoundName2 <- rep(dfcompound, times = 37)
compoundName1 <- rep(dfcompound, each = 37)
label <- paste(compoundName1, compoundName2, sep = "_")
colnames(CxC) <- label
index3 <- seq(1, 1369, by = 38)
compound_selfcross <- CxC[, -index3]
transposedIndexed_compound <- t(compound_selfcross)
index4 <- which(duplicated(transposedIndexed_compound))
removed_compound <- transposedIndexed_compound[-index4, ]
compound_finalcrossterms <- t(removed_compound)
CxC <- compound_finalcrossterms

#compound
C <- compound
#protein
P <- protein
#CxP
#CxC
#PxP
C_P <- cbind(C, P)
C_P_CxP_data_block_scale <- cbind(C, P, CxP) * (1/sqrt(length(C)+length(P)+length(CxP)))
#A_B_AxB_data <- cbind(affinity, A_B_AxB_data_block_scale)
C_P_CxC_data_block_scale <- cbind(C, P, 
                                  CxC) * (1/sqrt(length(C)+length(P)+length(CxC)))
#A_B_AxA_data <- cbind(affinity, A_B_AxA_data_block_scale)
C_P_PxP_data_block_scale <- cbind(C, P,
                                  PxP) * (1/sqrt(length(C)+length(P)+length(PxP)))
#A_B_BxB_data <- cbind(affinity, A_B_BxB_data_block_scale)
C_P_CxP_CxC_data_block_scale <- cbind(C, P, CxP,
                                      CxC) * (1/sqrt(length(C)+length(P)+length(CxP)+length(CxC)))
#A_B_AxB_AxA_data <- cbind(affinity, A_B_AxB_AxA_data_block_scale)
C_P_CxP_PxP_data_block_scale <- cbind(C, P, CxP,
                                      PxP) * (1/sqrt(length(C)+length(P)+length(CxP)+length(PxP)))
#A_B_AxB_BxB_data <- cbind(affinity, A_B_AxB_BxB_data_block_scale)
C_P_CxC_PxP_data_block_scale <- cbind(C, P, CxC,
                                      PxP) * (1/sqrt(length(C)+length(P)+length(CxC)+length(PxP)))
#A_B_AxA_BxB_data <- cbind(affinity, A_B_AxA_BxB_data_block_scale)
C_P_CxP_CxC_PxP_data_block_scale <- cbind(C, P, CxP, CxC, PxP) * (1/sqrt(length(C)+length(P)+
                                                    length(CxC)+length(PxP)))
#A_B_AxB_AxA_BxB_data <- cbind(affinity, A_B_AxB_AxA_BxB_data_block_scale)

C <- cbind(pIC50, C)
P <- cbind(pIC50, P)
CxP <- cbind(pIC50, CxP)
CxC <- cbind(pIC50, CxC)
PxP <- cbind(pIC50, PxP)
C_P <- cbind(pIC50, C_P)
C_P_CxP <- cbind(pIC50, C_P_CxP_data_block_scale)
C_P_CxC <- cbind(pIC50, C_P_CxC_data_block_scale)
C_P_PxP <- cbind(pIC50, C_P_PxP_data_block_scale)
C_P_CxP_CxC <- cbind(pIC50, C_P_CxP_CxC_data_block_scale)
C_P_CxP_PxP <- cbind(pIC50, C_P_CxP_PxP_data_block_scale)
C_P_CxC_PxP <- cbind(pIC50, C_P_CxC_PxP_data_block_scale)
C_P_CxP_CxC_PxP <- cbind(pIC50, C_P_CxP_CxC_PxP_data_block_scale)


pls_training <- function(x){
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
  tune <- train(pIC50 ~., data = x,  method = "pls", tuneLength = 10,
               trControl = ctrl)
  results <- list(100)
  for (i in 1:100) {
    sel <- naes(x, k = 90, pc = 5, iter.max = 100)
    Train <- x[sel$model, ]
    Test <- x[sel$test, ]
    model <- plsr(pIC50~., data = Train, ncomp = tune$bestTune[[1]])
    prediction <- predict(model, Train, ncomp = tune$bestTune[[1]])
    value <- data.frame(obs = Train$pIC50, pred = prediction)
    labeling <- c("obs", "pred")
    colnames(value) <- labeling
    results[[i]] <- defaultSummary(value)
  }
  return(results)
}

mean_and_sd <- function(x) {
  c(round(rowMeans(x, na.rm = TRUE), digits = 4),
    round(genefilter::rowSds(x, na.rm = TRUE), digits = 4))
}

pls_train <- function(x) {
  ok <- pls_training(x)
  data <- data.frame(ok)
  result <- mean_and_sd(data)
  df <- data.frame(result)
  R2_and_RMSE <- t(df)
  label <- c("RMSE_Mean", "Rsquared_Mean", "RMSE_SD", "Rsquared_SD")
  colnames(R2_and_RMSE) <- label
  return(R2_and_RMSE)
}

### name of the data C P CxP CxC PxP C_P C_P_CxP C_P_CxC C_P_PxP C_P_CxP_CxC C_P_CxP_PxP C_P_CxC_PxP C_P_CxP_CxC_PxP

input <- list(C = C, P = P, CxP = CxP, CxC = CxC, PxP = PxP, C_P = C_P,
              C_P_CxP = C_P_CxP, C_P_CxC = C_P_CxC, C_P_PxP = C_P_PxP, 
              C_P_CxP_CxC = C_P_CxP_CxC, C_P_CxP_PxP = C_P_CxP_PxP,
              C_P_CxC_PxP = C_P_CxC_PxP, C_P_CxP_CxC_PxP = C_P_CxP_CxC_PxP)
### Training for PLS
set.seed(3009)
results_PLS_training <- lapply(input, function(x) {
  models <- pls_train(x)
  return(models)
})
#### Cross Validation for PLS

pls_cross_validation <- function(x){
  tune <- train(pIC50 ~., data = x,  method = "pls", tuneLength = 10,
                trControl = ctrl)
  results <- list(100)
  for (i in 1:100) {
    sel <- naes(x, k = 90, pc = 5, iter.max = 100)
    myData <- x[sel$model, ]
    Test <- x[sel$test, ]
    k = 10
    index <- sample(1:k, nrow(myData), replace = TRUE)
    folds <- 1:k
    myRes <- data.frame()
    for (j in 1:k)
      training <- subset(myData, index %in% folds[-j])
    testing <- subset(myData, index %in% c(j))
    ctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 1)
    tune <- train(pIC50 ~., data = training, method = "pls", tunLength = 10, trControl = ctrl)
    model <- plsr(pIC5~., data = training, ncomp = tune$bestTune[[1]])
    prediction <- predict(model, testing, ncomp = tune$bestTune[[1]])
    value <- data.frame(obs = testing$pIC50, pred = prediction)
    labeling <- c("obs", "pred")
    colnames(value) <- labeling
    results[[i]] <- defaultSumary(value)
  }
  return(results)
}

pls_10_CV <- function(x) {
  ok <- pls_cross_validation(x)
  data <- data.frame(ok)
  result <- mean_and_sd(data)
  df <- data.frame(result)
  R2_and_RMSE <- t(df)
  label <- c("RMSE_Mean", "Rsquared_Mean", "RMSE_SD", "Rsquared_SD")
  colnames(R2_and_RMSE) <- label
  return(R2_and_RMSE)
}
### results for 10 fold CV
set.seed(39)
results_PLS_10_CV <- lapply(input, function(x) {
  models <- pls_10_CV(x)
  return(models)
})

#### PLS Testing

pls_testing <- function(x){
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
  tune <- train(pIC50 ~., data = x,  method = "pls", tuneLength = 10,
                trControl = ctrl)
  results <- list(100)
  for (i in 1:100) {
    sel <- naes(x, k = 90, pc = 5, iter.max = 100)
    Train <- x[sel$model, ]
    Test <- x[sel$test, ]
    model <- plsr(pIC50~., data = Train, ncomp = tune$bestTune[[1]])
    prediction <- predict(model, Test, ncomp = tune$bestTune[[1]])
    value <- data.frame(obs = Test$pIC50, pred = prediction)
    labeling <- c("obs", "pred")
    colnames(value) <- labeling
    results[[i]] <- defaultSummary(value)
  }
  return(results)
}

pls_test <- function(x) {
  ok <- pls_testing(x)
  data <- data.frame(ok)
  result <- mean_and_sd(data)
  df <- data.frame(result)
  R2_and_RMSE <- t(df)
  label <- c("RMSE_Mean", "Rsquared_Mean", "RMSE_SD", "Rsquared_SD")
  colnames(R2_and_RMSE) <- label
  return(R2_and_RMSE)
}

set.seed(300)
results_PLS_testing <- lapply(input, function(x) {
  models <- pls_test(x)
  return(models)
})

##### leave compound out

data <- read_excel("aromatase_data.xlsx")
compound_index <- data$Compound
leave_one_compound <- function(x) {
dat <- cbind(compound_index, x)
subs <- unique(dat$compound_index)
model_these <- vector(mode = "list", length = length(subs))
for (i in seq_along(subs))
  model_these[[i]] <- which(dat$compound_index != subs[i])
names(model_these) <- paste0("compound_index", subs)
LOCO <- train(x = dat[, 3:length(dat)],
              y = dat[, 2],
              method = "pls",
              trControl = ctrl)
prediction <- predict(LOCO, dat)
value <- data.frame(obs = dat$pIC50, pred = prediction)
labeling <- c("obs", "pred")
colnames(value) <- labeling
results <- defaultSummary(value)
return(results)
}

set.seed(100)
results_PLS_LOCO <- lapply(input, function(x) {
  models <- leave_one_compound(x)
  return(models)
})

##### leave protein out

data <- read_excel("aromatase_data.xlsx")
protein_index <- data$Protein
leave_one_protein <- function(x) {
  dat <- cbind(protein_index, x)
  subs <- unique(dat$protein_index)
  model_these <- vector(mode = "list", length = length(subs))
  for (i in seq_along(subs))
    model_these[[i]] <- which(dat$protein_index != subs[i])
  names(model_these) <- paste0("protein_index", subs)
  LOPO <- train(x = dat[, 3:length(dat)],
                y = dat[, 2],
                method = "pls",
                trControl = ctrl)
  prediction <- predict(LOPO, dat)
  value <- data.frame(obs = dat$pIC50, pred = prediction)
  labeling <- c("obs", "pred")
  colnames(value) <- labeling
  results <- defaultSummary(value)
  return(results)
}

set.seed(90)
results_PLS_LOPO <- lapply(input, function(x) {
  models <- leave_one_protein(x)
  return(models)
})

#### Models for each compounds
data <- read_excel("aromatase_data.xlsx")
compound_index <- data$Compound
compound <- data[, 5:311]
QSAR_data <- cbind(compound_index, pIC50, compound)
subs <- unique(QSAR_data$compound_index)

QSAR_input <- list(OHA = "4_OHA", alphaAPTADD = "7alpha-APTADD", AG = "AG", 
                     CGS20267 = "CGS20267", ICID1033 = "ICID1033", MDL101003 = "MDL101003", 
                   MR20492 = "MR20492", MR20494 = "MR20494", MR20814 = "MR20814", Vorozole= "Vorozole")


results_QSAR_PLS_Training <- lapply(QSAR_input, function(x) {
  data <- subset(QSAR_data, compound_index = "x")
  data <- data[, 2:length(data)]
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
  tune <- suppressWarnings(train(x = data[, 2:length(data)], y = data[, 1],
                                 method = "pls", tuneLength = 10, trControl = ctrl))
  ok <- list(100)
  for (i in 1:100) {
    sel <- naes(data, k = 7, pc = 5, iter.max = 100)
  Train <- data[sel$model, ]
  Test <- data[sel$test, ]
  model <- plsr(pIC50~., data = Train, ncomp = 1)
  prediction <- predict(model, Train, ncomp = 1)
  value <- data.frame(obs = Train$pIC50, pred = prediction)
  labeling <- c("obs", "pred")
  colnames(value) <- labeling
  ok[[i]] <- defaultSummary(value)
  }
  data <- data.frame(ok)
  result <- mean_and_sd(data)
  df <- data.frame(result)
  R2_and_RMSE <- t(df)
  label <- c("RMSE_Mean", "Rsquared_Mean", "RMSE_SD", "Rsquared_SD")
  colnames(R2_and_RMSE) <- label
  return(R2_and_RMSE)
})

### Testing
results_QSAR_PLS_testing <- lapply(QSAR_input, function(x) {
  data <- subset(QSAR_data, compound_index = "x")
  data <- data[, 2:length(data)]
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 1)
  tune <- suppressWarnings(train(x = data[, 2:length(data)], y = data[, 1],
                                 method = "pls", tuneLength = 10, trControl = ctrl))
  ok <- list(100)
  for (i in 1:100) {
    sel <- naes(data, k = 7, pc = 5, iter.max = 100)
    Train <- data[sel$model, ]
    Test <- data[sel$test, ]
    model <- plsr(pIC50~., data = Train, ncomp = 1)
    prediction <- predict(model, Test, ncomp = 1)
    value <- data.frame(obs = Test$pIC50, pred = prediction)
    labeling <- c("obs", "pred")
    colnames(value) <- labeling
    ok[[i]] <- defaultSummary(value)
  }
    data <- data.frame(ok)
    result <- mean_and_sd(data)
    df <- data.frame(result)
    R2_and_RMSE <- t(df)
    label <- c("RMSE_Mean", "Rsquared_Mean", "RMSE_SD", "Rsquared_SD")
    colnames(R2_and_RMSE) <- label
    return(R2_and_RMSE)
})


results_QSAR_PLS_10_fold_CV <- lapply(QSAR_input, function(x) {
  data <- subset(QSAR_data, compound_index = "x")
  data <- data[, 2:length(data)]
  ok <- list(100)
  for (i in 1:100) {
    sel <- naes(data, k = 7, pc = 5, iter.max = 100)
    myData <- data[sel$model, ]
    Test <- data[sel$test, ]
    k = 10
    index <- sample(1:k, nrow(myData), replace = TRUE)
    folds <- 1:k
    myRes <- data.frame()
    for (j in 1:k)
      training <- subset(myData, index %in% folds[-j])
    testing <- subset(myData, index %in% c(j))
    ctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 1)
    tune <- suppressWarnings(train(pIC50 ~., data = training,
                                   method = "pls", tunLength = 10, trControl = ctrl))
    model <- plsr(pIC5~., data = training, ncomp = 1)
    prediction <- predict(model, testing, ncomp = 1)
    value <- data.frame(obs = testing$pIC50, pred = prediction)
    labeling <- c("obs", "pred")
    colnames(value) <- labeling
    ok[[i]] <- defaultSumary(value)
  }
  data <- data.frame(ok)
  result <- mean_and_sd(data)
  df <- data.frame(result)
  R2_and_RMSE <- t(df)
  label <- c("RMSE_Mean", "Rsquared_Mean", "RMSE_SD", "Rsquared_SD")
  colnames(R2_and_RMSE) <- label
  return(R2_and_RMSE)
})
























