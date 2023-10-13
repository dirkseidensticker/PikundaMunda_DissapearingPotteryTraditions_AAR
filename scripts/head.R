
library(ggplot2)
library(tidyr)

library(FactoMineR)
library(factoextra)

#if(!require('devtools')) install.packages('devtools')
#devtools::install_github('ISAAKiel/quantAAR')

library(quantAAR)

obj <- data.table::fread("input/TransGenTN_Samples_Excerpt.csv", encoding = "UTF-8")


