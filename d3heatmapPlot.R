##aromatase d3heatmap plot
## For compound 
library(caret)
library(readxl)
library(d3heatmap)
library(dplyr)

data <- read_excel("aromatase_data.xlsx")
compound <- data[, 5:311]
set.seed(112)
compound <- compound[, -nearZeroVar(compound)]
compoundname <- as.character(unique(data$Compound))
dataframe <- cbind(Compound = data$Compound, compound)
sample_compound <- function(x) {
  myRef <- data.frame()
  for (i in compoundname) {
    index_data <- which(x$Compound == i)
    ok <- x[index_data, ]
    sample <- sample_n(ok, size = 1, replace = TRUE)
    myRef <- rbind(myRef, sample)
  }
  return(myRef)
}
df <- sample_compound(dataframe)
names <- df$Compound
compound <- df[, 2:ncol(df)]
rownames(compound) <- names
d3heatmap(compound, scale = "column", colors = "Spectral")
#### D3 plot for protein
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
df <- sample_protein(dataframe)
set.seed(100)
names <- df$Protein
protein <- df[, 2:ncol(df)]
rownames(protein) <- names
d3heatmap(protein, scale = "column", colors = "Spectral")



