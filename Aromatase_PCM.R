library(readxl)
library(prospectr)
library(Rcpi)
library(pls)
library(caret)
library(dplyr)
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
dfprotein <- colnames(protein)
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

### PLS Coefficient
pls_Coefficient <- function(x){
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
    tune <- train(pIC50 ~., data = training, method = "pls", tuneLength = 10, trControl = ctrl)
    model <- plsr(pIC50~., data = training, ncomp = tune$bestTune[[1]])
    ok <- model$coefficients
    yay <- data.frame(ok)
    myRes <- rbind(myRes, yay)
    mean <- apply(data.frame(myRes), 1, mean)
    sd <- apply(data.frame(myRes), 1, sd)
    results[[i]] <- cbind(mean, sd)
  }
  return(results)
}


coefficient_input <- list(C_P = C_P, C_P_CxP_CxC = C_P_CxP_CxC, C_P_CxP_CxC_PxP = C_P_CxP_CxC_PxP)
coefficient_input <- list(C_P = C_P)


set.seed(39)
results_PLS_coefficient <- lapply(coefficient_input, function(x) {
  models <- pls_Coefficient(x)
  return(models)
})






##### chemical space
data <- read_excel("aromatase_data.xlsx")
compound <- data[, 5:311]
set.seed(122)
compound = compound[, -nearZeroVar(compound)]
compound[compound == "0"] <- -1

compoundname <- as.character(unique(data$Compound))
dataframe <- cbind(Compound = data$Compound, compound)
sample_compound <- function(x) {
  myRef  <- data.frame()
  for (i in compoundname) {
    index_data <- which(x$Compound == i)
    ok <- x[index_data, ]
    sample <- sample_n(ok, size = 1, replace = TRUE)
    myRef <- rbind(myRef, sample)
  }
  return(myRef)
}
df <- sample_compound(dataframe)[, 2:length(dataframe)]
set.seed(2000)
df = df[, -nearZeroVar(df)]
### optium PC with Horn's Parallel Analysis
paran(df, iterations = 5000)
df.sc <- scale(df)
df.svd <- svd(df.sc)
df.scores <- df.svd$u %*% diag(df.svd$d)
df.loadings <- df.svd$v
df.vars <- df.svd$d^2 / (nrow(df) -1)
df.totalvar <- sum(df.vars)
df.relvars <- df.vars / df.totalvar
variances <- 100 * round(df.relvars, digits = 3)
variances[1:5]
par(mfrow = c(2,2))
barplot(df.vars[1:10], main = "Variances",
        names.arg = paste("PC", 1:10))
barplot(log(variances[1:10]), main = "log(Variances)", ylim = c(0, 4),
        names.arg = paste("PC", 1:10))
barplot(100*df.relvars[1:10], main = "Relative variances (%)",
        names.arg = paste("PC", 1:10))
barplot(cumsum(100*df.relvars[1:10]),
        main = "Cumulative variances (%)",
        names.arg = paste("PC", 1:10), ylim = c(0, 100))
        
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
screeplot(pca, type = "line", npcs = 10)
library(nFactors)
cv = eigen(cor(df))
ap = parallel(subject= nrow(df), var = ncol(swiss), rep = 100, cent = .05)
nS = nScree(x = ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)
scores <- pca$x[,1:5]
loadings <- pca$rotation[,1:5]
km <- kmeans(scores, center=2, nstart=5)
ggdata <- data.frame(scores, Cluster=km$cluster)
compoundnumber <- c(4, 5, 6, 10, 2, 1, 3, 8, 9, 7)
#ggdata <- cbind(compoundname, ggdata)
ggdata <- cbind(compoundnumber, ggdata)
### paper numbering
library(grid)
set.seed(23)
x <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle(" ") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=compoundnumber), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    #plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-15, 15), xlim = c(-15, 15))
  #scale_y_continuous(limits = c(-15,15), labels = comma) +
  #scale_x_continuous(limits = c(-15,15), labels = comma) 
### loading plot
data <- read_excel("aromatase_data.xlsx")
compound <- data[, 5:311]
set.seed(122)
compound = compound[, -nearZeroVar(compound)]
compound[compound == "0"] <- -1

compoundname <- as.character(unique(data$Compound))
dataframe <- cbind(Compound = data$Compound, compound)
sample_compound <- function(x) {
  myRef  <- data.frame()
  for (i in compoundname) {
    index_data <- which(x$Compound == i)
    ok <- x[index_data, ]
    sample <- sample_n(ok, size = 1, replace = TRUE)
    myRef <- rbind(myRef, sample)
  }
  return(myRef)
}
df <- sample_compound(dataframe)[, 2:length(dataframe)]
set.seed(2000)
df = df[, -nearZeroVar(df)]
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
scores <- pca$x[,1:5]
set.seed(300)
loadings <- pca$rotation[,1:5]
km <- kmeans(loadings, center=3, nstart=5)
ggdata <- data.frame(loadings, Cluster=km$cluster)
labelcompound <- rownames(loadings)
ggdata <- cbind(labelcompound, ggdata)
### paper numbering
library(grid)
set.seed(23)
a <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle(" ") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=labelcompound), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    #plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-0.8, 0.8), xlim = c(-0.8, 0.8))

## Protein  space
data <- read_excel("aromatase_data.xlsx")
protein <- data[, 313:351]
set.seed(200)
proteinname <- as.character(unique(data$Protein))
dataframe <- cbind(Protein = data$Protein, protein)
sample_protein <- function(x) {
  myRef  <- data.frame()
  for (i in proteinname) {
    index_data <- which(x$Protein == i)
    ok <- x[index_data, ]
    sample <- sample_n(ok, size = 1, replace = TRUE)
    myRef <- rbind(myRef, sample)
  }
  return(myRef)
}
df <- sample_protein(dataframe)[, 2:length(dataframe)]
set.seed(100)
df = df[, -nearZeroVar(df)]
paran <- paran(df, iterations = 5000)
df.sc <- scale(df)
df.svd <- svd(df.sc)
df.scores <- df.svd$u %*% diag(df.svd$d)
df.loadings <- df.svd$v
df.vars <- df.svd$d^2 / (nrow(df) -1)
df.totalvar <- sum(df.vars)
df.relvars <- df.vars / df.totalvar
variances <- 100 * round(df.relvars, digits = 3)
variances[1:5]
par(mfrow = c(2,2))
barplot(df.vars[1:10], main = "Variances",
      names.arg = paste("PC", 1:10))
barplot(log(variances[1:10]), main = "log(Variances)",
        names.arg = paste("PC", 1:10))
barplot(100*df.relvars[1:10], main = "Relative variances (%)",
     names.arg = paste("PC", 1:10))
barplot(cumsum(100*df.relvars[1:10]),
 main = "Cumulative variances (%)",
names.arg = paste("PC", 1:10), ylim = c(0, 100))

pca <- prcomp(df, retx=TRUE,scale.=TRUE)
scores <- pca$x[,1:5]
loadings <- pca$rotation[,1:5]
km <- kmeans(scores, center=1, nstart=5)
ggdata <- data.frame(scores, Cluster=km$cluster)
ggdata <- cbind(proteinname, ggdata)
### paper numbering
library(grid)
set.seed(2300)
y <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle(" ") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=proteinname), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
   # plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 5))
#scale_y_continuous(limits = c(-15,15), labels = comma) +
#scale_x_continuous(limits = c(-15,15), labels = comma) 
#grid.arrange(x, y, ncol = 1, main = " ")

### loading plot
protein <- data[, 313:351]
set.seed(200)
proteinname <- as.character(unique(data$Protein))
dataframe <- cbind(Protein = data$Protein, protein)
sample_protein <- function(x) {
  myRef  <- data.frame()
  for (i in proteinname) {
    index_data <- which(x$Protein == i)
    ok <- x[index_data, ]
    sample <- sample_n(ok, size = 1, replace = TRUE)
    myRef <- rbind(myRef, sample)
  }
  return(myRef)
}
df <- sample_protein(dataframe)[, 2:length(dataframe)]
set.seed(100)
df = df[, -nearZeroVar(df)]
column_names <- colnames(df)
label_names2 <- gsub("zscl", "z", column_names)
colnames(df) <- label_names2
pca <- prcomp(df, retx=TRUE,scale.=TRUE)
scores <- pca$x[,1:5]
loadings <- pca$rotation[,1:5]
km <- kmeans(loadings, center=2, nstart=5)
ggdata <- data.frame(loadings, Cluster=km$cluster)
labelprotein <- rownames(loadings)
ggdata <- cbind(labelprotein, ggdata)
### paper numbering
library(grid)
set.seed(23)
b <- ggplot(ggdata, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(fill=factor(Cluster)), size=5, shape=20, pch = 21, alpha = 0.8) +
  ggtitle(" ") +
  stat_ellipse(aes(fill=factor(Cluster)), colour = "black", 
               geom="polygon", level=0.95, alpha=0.2) +
  guides(color=guide_legend("Cluster"),fill=guide_legend("Cluster")) +
  #geom_text(aes(label=labelprotein), size=7, hjust=0.5, vjust= 1.5, alpha=0.45) +
  theme(
    legend.position=("none"),
    #plot.title = element_text(size=20, face="bold", colour="black", vjust = 2, hjust=-0.07),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text.x = element_text(size = 15),
    legend.title=element_blank(),
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20)) +
  coord_cartesian(ylim = c(-1.5, 1.5), xlim = c(-1.5, 1.5))

plot_grid(x, y,a, b,  labels = c("A", "B", "C", "D"), ncol = 2,  label_size = 20)

##### PLS Coefficents Feature Importance Anaplysis 100 times 
input <- list(C = C, P = P, CxP = CxP, CxC = CxC, PxP = PxP, C_P = C_P,
              C_P_CxP = C_P_CxP, C_P_CxC = C_P_CxC, C_P_PxP = C_P_PxP, 
              C_P_CxP_CxC = C_P_CxP_CxC, C_P_CxP_PxP = C_P_CxP_PxP,
              C_P_CxC_PxP = C_P_CxC_PxP, C_P_CxP_CxC_PxP = C_P_CxP_CxC_PxP)


pls_Coefficient <- function(x){
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
    tune <- train(pIC50 ~., data = training, method = "pls", tuneLength = 10, trControl = ctrl)
    model <- plsr(pIC50~., data = training, ncomp = tune$bestTune[[1]])
    ok <- model$coefficients
    yay <- data.frame(ok)
    myRes <- rbind(myRes, yay)
    mean <- apply(data.frame(myRes), 1, mean)
    sd <- apply(data.frame(myRes), 1, sd)
    results[[i]] <- cbind(mean, sd)
  }
  return(results)
}


coefficient_input <- list(C_P = C_P, C_P_CxP_CxC = C_P_CxP_CxC, C_P_CxP_CxC_PxP = C_P_CxP_CxC_PxP)
coefficient_input <- list(C_P = C_P)


set.seed(39)
results_PLS_coefficient <- lapply(coefficient_input, function(x) {
  models <- pls_Coefficient(x)
  return(models)
})

### PLS coefficient PLot

##Model6
Model6 <- read.csv("Model6_coefficient.csv", header = TRUE)
labels <- c("descriptors", "mean", "sd")
colnames(Model6) <- labels
descriptors <- Model6$descriptors
descriptors <- gsub("zscl", "z", descriptors, perl = TRUE)
Model6 <- apply(Model6[, 2:3], 2, round, digits = 4)
Model6 <- data.frame(descriptors, Model6)
Model6$colour <- ifelse(Model6$mean < 0, "negative", "positive")
a <- ggplot(Model6, aes(x = descriptors, y = mean, fill = colour)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(.9)) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), colour = "black", position = position_dodge(.9)) +
  theme( 
    # plot.title = element_text(size=70, face="bold", colour="black", vjust = 2, hjust = -0.145),
    legend.position=("none"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x= element_blank(),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  ylab("PLS Coefficient") +
  scale_fill_manual(values = c(positive = "firebrick1", negative = "steelblue"))


# Model 10
Model10 <- read.csv("Model10_coefficient.csv", header = TRUE)
labels <- c("descriptors", "mean", "sd")
colnames(Model10) <- labels
descriptors <- Model10$descriptors
descriptors <- gsub("zscl", "z", descriptors, perl = TRUE)
Model10 <- apply(Model10[, 2:3], 2, round, digits = 3)
Model10 <- data.frame(descriptors, Model10)
Model10$colour <- ifelse(Model10$mean < 0, "negative", "positive")
b <- ggplot(Model10, aes(x = descriptors, y = mean, fill = colour)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(.9)) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), colour = "black", position = position_dodge(.9)) +
  theme( 
    # plot.title = element_text(size=70, face="bold", colour="black", vjust = 2, hjust = -0.145),
    legend.position=("none"),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x= element_blank(),
    axis.title.y = element_text(size = 20),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1)) +
  ylab("PLS Coefficient") +
  scale_fill_manual(values = c(positive = "firebrick1", negative = "steelblue"))



# Model13
Model13 <- read.csv("Model13_coefficient.csv", header = TRUE)
labels <- c("descriptors", "mean", "sd")
colnames(Model13) <- labels
descriptors <- Model13$descriptors
Model13 <- apply(Model13[, 2:3], 2, round, digits = 3)
Model13 <- data.frame(descriptors, Model13)
Model13$colour <- ifelse(Model13$mean < 0, "negative", "positive")


### predictive models with 7.0 cutoff PCM models 

correlation_cutoff <- function(x) {
  library(caret)
  data <- x[-1]
  corr <- cor(data)
  data_corr <- corr[1:as.numeric(ncol(corr)), 1:as.numeric(ncol(corr))]
  high <- findCorrelation(data_corr, cutoff = .70)
  processed_data <- data[, -high]
  return(processed_data)
}

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

#Preparing the input for the filtering of descriptors higher than .70


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


set.seed(300)
output_descriptors_filtered <- lapply(input, function(x) {
  models <- correlation_cutoff(x)
  return(models)
})

C <- cbind(pIC50, output_descriptors_filtered$C)
P <- cbind(pIC50, output_descriptors_filtered$P)
CxP <- cbind(pIC50, output_descriptors_filtered$CxP)
CxC <- cbind(pIC50, output_descriptors_filtered$CxC)
PxP <- cbind(pIC50, output_descriptors_filtered$PxP)
C_P <- cbind(pIC50, output_descriptors_filtered$C_P)
C_P_CxP <- cbind(pIC50, output_descriptors_filtered$C_P_CxP_data_block_scale)
C_P_CxC <- cbind(pIC50, output_descriptors_filtered$C_P_CxC_data_block_scale)
C_P_PxP <- cbind(pIC50, output_descriptors_filtered$C_P_PxP_data_block_scale)
C_P_CxP_CxC <- cbind(pIC50, output_descriptors_filtered$C_P_CxP_CxC_data_block_scale)
C_P_CxP_PxP <- cbind(pIC50, output_descriptors_filtered$C_P_CxP_PxP_data_block_scale)
C_P_CxC_PxP <- cbind(pIC50, output_descriptors_filtered$C_P_CxC_PxP_data_block_scale)
C_P_CxP_CxC_PxP <- cbind(pIC50, output_descriptors_filtered$C_P_CxP_CxC_PxP_data_block_scale)

