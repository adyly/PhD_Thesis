# Script with functions for the project

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(tidyr)

# theme for plotting
theme_bw_extended <- function(){theme_bw() + theme(
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text = element_text(color = 'black'),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = 'black'),
  axis.ticks = element_line(color = 'black')
)
}


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




raw_data<-read.csv('DIAJodiplusAth-first-pass.pg_matrix.csv',header = T)

colnames(raw_data[4]) <- 'Gene.names'


# normalize data and pca plot ---------------------------------------------

normalized_data <- NormalizeTotalIntensity(raw_data, 'REP')


pca_data <- prcomp(x= na.omit(normalized_data[grepl('_n', colnames(normalized_data))]), scale. = TRUE)
View(pca_data$rotation)

pca_plot_data <- data.frame(pca_data$rotation)

# version without batch correction
pca_plot_data$names <- rownames(pca_plot_data)
no_correction <- pca_plot_data %>% 
  mutate(operator = case_when(str_detect(names, 'jodi') ~ 'Jodi',
                              TRUE ~ 'Athanasia'),
         cell_line = substr(names, 1,3),
         cell_compartment = substr(names, 5, 5)) %>% 
  ggplot(aes(x= PC1, y = PC2, color = cell_line, shape = operator, size = cell_compartment))+
  geom_point()+
  ggtitle('Not batch corrected')+
  guides(size=guide_legend(title="Cell compartment"),
         shape=guide_legend(title="Operator"),
         color=guide_legend(title="Cell line"))+
  theme_bw()







# redistributing data -----------------------------------------------------

normalized_data_2 <- normalized_data %>% 
  dplyr::select(Accession, Protein.Ids, Protein.Names, Gene.names, First.Protein.Description, ends_with('_n')) %>% 
  pivot_longer(cols = ends_with('_n'), names_to = 'sample', values_to = 'normalized_intensities') %>% 
  mutate(sample = substr(sample, 1, nchar(sample)-2),
         cell_line = substr(sample, 1,3),
         cell_compartment = case_when(str_detect(sample, 'N') ~ 'Nucleus',
                                      TRUE ~ 'Whole_lysate'),
         operator = case_when(str_detect(sample, 'jodi') ~ 'Jodi',
                              TRUE ~ 'Athanasia'),
         replicate = case_when(str_detect(sample, 'WL') ~ substr(sample, 12, 12),
                               TRUE ~ substr(sample, 11,11))) %>% 
  group_by(cell_line, cell_compartment, operator, Gene.names, Accession) %>% # keep accession since some gene.names are double, and some has no gena name
  summarise(mean_normalized_intensity = mean(normalized_intensities, na.rm = TRUE)) %>% # does not do anything since Accession is unique
  pivot_wider(id_cols = c(cell_line, operator, Gene.names, Accession), names_from = cell_compartment, values_from = mean_normalized_intensity) %>% 
  mutate(ratio = Nucleus/Whole_lysate,
         log2fc = log2(ratio))



# histogram of log2fc
log2fc_histogram <- normalized_data_2 %>% 
  mutate(cell_line = factor(cell_line, levels = c('E13', 'E28', 'E21', 'E57', 'E31', 'E34')),
         cell_line_operator = paste(cell_line, operator, sep = ' ')) %>% 
  filter(cell_line_operator != "E13 Jodi") %>% 
  ggplot(aes(x = log2fc))+
  #facet_grid(operator~cell_line)+
  facet_wrap(cell_line~operator, ncol = 2)+
  geom_histogram(bins = 100, color = 'black', fill = 'white')


log2fc_histogram  
#ggsave('log2fc_histogram.png', log2fc_histogram)



