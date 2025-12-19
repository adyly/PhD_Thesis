

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
#library(ggvenn)
library(ggrepel)
#library(writexl)


# function for normalizing each channel by total intensity of each channel
NormalizeTotalIntensity <- function(data, match){
  
  # print the number of matches to see that it works as expected
  print('normalized intensity columns:')
  print(sum(grepl(match, colnames(data))))
  
  # calculate normalization factors
  normalizationFactors <- colSums(data[grepl(match, colnames(data))], na.rm = TRUE)
  
  # apply the normalization factors over the approapriate columns
  normalized_data <- sweep(data[grepl(match, colnames(data))], 2, normalizationFactors, FUN = '/')
  
  # change the colnames
  colnames(normalized_data) <- paste0(colnames(data[grepl(match, colnames(data))]), '_n')
  
  return(cbind(data, normalized_data))
}