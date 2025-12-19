# script for making volcanos and other plots for the thesis
source('Normalisation.R')
library(plotly)

# read set of adhesion proteins
library(readxl)
adhesion_proteins <- read_xlsx('adhesion.xlsx', sheet = 'adhesion')
adhesion_proteins <- adhesion_proteins %>% 
  dplyr::select(Gene, Category) %>% 
  rename(Gene.names = Gene)


# make list with ensembl gene name convertion -----------------------------

library(biomaRt)

markers <- read.csv('markers.csv')
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

IDs <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
             filters = "ensembl_gene_id", values = markers$Gene,
             mart = mart) %>% 
  rename(Gene.names = hgnc_symbol,
         Gene = ensembl_gene_id)

markers <- markers %>% 
  left_join(IDs, by = 'Gene')

#write.csv(markers, file = 'markers_location_gene.csv')

uncertain_markers <- read.csv('uncertain_markers.tsv', sep = '\t') %>% 
  rename(Gene.names = Gene.name,
         Location_Approved = Main.location)

normalized_data_markers <- normalized_data %>% 
  left_join(markers, by = 'Gene.names')

completeData_markers <- completeData %>% 
  left_join(markers, by = "Gene.names")



# Volcano  -------------------------------------------

colnames(normalized_data)[1] <- "Accession"

nonimputed_volcano <- normalized_data %>% 
  select(Accession, Protein.Ids, Protein.Names, Gene.names, First.Protein.Description, ends_with('_n')) %>% 
  na.omit() %>% 
  pivot_longer(cols = ends_with('_n'), names_to = 'sample', values_to = 'normalized_intensities') %>% 
  mutate(sample = substr(sample, 1, nchar(sample)-2),
         cell_line = substr(sample, 1,3),
         cell_compartment = case_when(str_detect(sample, 'N') ~ 'Nucleus',
                                      TRUE ~ 'Whole_lysate'),
         operator = case_when(str_detect(sample, 'jodi') ~ 'Jodi',
                              TRUE ~ 'Athanasia'),
         replicate = case_when(str_detect(sample, 'WL') ~ substr(sample, 12, 12),
                               TRUE ~ substr(sample, 11,11))) %>% 
  filter(sample %in% c("E13_WL_REP_1", "E13_WL_REP_2", "E13_WL_REP_3", "E13_N_REP_1", "E13_N_REP_2", "E13_N_REP_3")) %>% 
  pivot_wider(id_cols = c("Accession", "Gene.names", "cell_line"), names_from = "cell_compartment", values_from = c("normalized_intensities")) %>% 
  group_by(Accession, Gene.names) %>% 
  mutate(pvalue = t.test(unlist(Nucleus), unlist(Whole_lysate))$p.value,
         ratio = mean(unlist(Nucleus))/mean(unlist(Whole_lysate)),
         log2ratio = log2(ratio)) %>% 
  left_join(adhesion_proteins, by = 'Gene.names') %>% 
  mutate(significant_adhesion = if_else(abs(log2ratio) > 1 & pvalue < 0.05, true = Category, false = NA)) %>% 
  ggplot(aes(x = log2ratio, y = -log10(pvalue), fill = significant_adhesion, label = Gene.names))+
  geom_point(pch = 21, color = 'black', size = 3)+
  theme_bw_extended()+
  geom_text_repel()+
  scale_color_manual(values = c('red', 'blue', 'green'), breaks = c('consensus', 'Meta', 'LitCur'), na.value = 'grey')

ggplotly(nonimputed_volcano)
# Save at .png
# save the widget at .html format
library(htmlwidgets)
saveWidget(ggplotly(nonimputed_volcano), file="myFile.html")

library(webshot)
webshot::install_phantomjs()
webshot("file:///Users/ady/Desktop/GSEA_Thesis/myFile.html" , "output.png", delay = 0.2)

ggsave(ggplotly(nonimputed_volcano), filename = 'test.png', width = 2000, height = 2000, units = 'px')


volcano_plot <- function(normalized_data, imputed_data = NA, op, cell, title){
  new_list <- scan(file = "nucleoadhesome.txt", what = character(), sep = "\n")
  new_list <- toupper(new_list)
  
  data <- normalized_data %>% 
    dplyr::select(Accession, Protein.Ids, Protein.Names, Gene.names, First.Protein.Description, ends_with('_n')) %>% 
    #na.omit() %>% 
    pivot_longer(cols = ends_with('_n'), names_to = 'sample', values_to = 'normalized_intensities') %>% 
    mutate(sample = substr(sample, 1, nchar(sample)-2),
           cell_line = substr(sample, 1,3),
           cell_compartment = case_when(str_detect(sample, 'N') ~ 'Nucleus',
                                        TRUE ~ 'Whole_lysate'),
           operator = case_when(str_detect(sample, 'jodi') ~ 'Jodi',
                                TRUE ~ 'Athanasia'),
           replicate = case_when(str_detect(sample, 'WL') ~ substr(sample, 12, 12),
                                 TRUE ~ substr(sample, 11,11))) %>% 
    filter(operator == op, cell_line == cell)  %>% 
    pivot_wider(id_cols = c("Accession", "Gene.names", "cell_line"), names_from = "cell_compartment", values_from = c("normalized_intensities")) %>% 
    rowwise() %>% 
    mutate(nucleus_NA = sum(is.na(Nucleus))>0,
           wl_NA = sum(is.na(Whole_lysate))>0) %>% 
    filter(nucleus_NA == FALSE,
           wl_NA == FALSE) %>% 
    group_by(Accession, Gene.names) %>% 
    mutate(pvalue = t.test(unlist(Nucleus), unlist(Whole_lysate))$p.value,
           ratio = mean(unlist(Nucleus))/mean(unlist(Whole_lysate)),
           log2ratio = log2(ratio)) %>% 
    left_join(adhesion_proteins, by = 'Gene.names') %>% 
    left_join(uncertain_markers, by = 'Gene.names') %>% 
    mutate(significant_adhesion = if_else(abs(log2ratio) > 1 & pvalue < 0.05, true = Category, false = NA),
           nuclear_proteins = if_else(grepl('Nuc', Location_Approved), true = "Nuclear", false = "Non nuclear"),
           significant = if_else(pvalue < 0.05, "Significant", "Non significant"),
           data = 'Non imputed') %>% 
    filter(!is.na(Location_Approved)) %>% 
    mutate(significant_adhesion = case_when(significant_adhesion == 'consensus' ~ 'Consensus',
                                            significant_adhesion == 'LitCur' ~ 'Literature Curated',
                                            grepl('^H2A', Gene.names) ~ 'Histone marker',
                                            grepl('^H2B', Gene.names) ~ 'Histone marker',
                                            grepl('^H3', Gene.names) ~ 'Histone marker',
                                            grepl('^H4', Gene.names) ~ 'Histone marker',
                                            grepl('^H5', Gene.names) ~ 'Histone marker',
                                            TRUE ~ significant_adhesion))
  
   
  all_data <- data %>%
    mutate(non_nuclear_significant_adhesion_protein = if_else(nuclear_proteins == 'Non nuclear' &
                                                                significant == 'Significant' &
                                                                !is.na(significant_adhesion) &
                                                                log2ratio > 1,
                                                              'Significant non nuclear adhesion protein',
                                                              'Other'))
  
  plot <- all_data %>% # all_data %>% 
    ggplot(aes(x = log2ratio, y = -log10(pvalue), fill = significant_adhesion, label = Gene.names, shape = nuclear_proteins, alpha = significant, size = non_nuclear_significant_adhesion_protein))+
    #facet_wrap(~data) +
    geom_point(data = all_data %>% filter(is.na(significant_adhesion)))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Meta', nuclear_proteins == 'Nuclear'))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Consensus', nuclear_proteins == 'Nuclear'))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Literature Curated', nuclear_proteins == 'Nuclear'))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Meta', nuclear_proteins == 'Non nuclear'))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Consensus', nuclear_proteins == 'Non nuclear'))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Literature Curated', nuclear_proteins == 'Non nuclear'))+
    geom_point(data = all_data %>% filter(significant_adhesion == 'Histone marker'))+
    theme_bw_extended()+
    #geom_text_repel()+
    scale_fill_manual(values = c('#003f5c', '#ffa600', '#bc5090', 'red'), breaks = c('Consensus', 'Meta', 'Literature Curated', 'Histone marker'), na.value = 'grey')+
    #scale_color_manual(values = c('RED', 'black'), breaks = c("Nuclear", "Non nuclear"))+
    scale_alpha_manual(values = c(0.5, 1), breaks = c("Non significant", "Significant"))+
    scale_shape_manual(values = c(21, 22), breaks = c("Nuclear", "Non nuclear"))+
    ggtitle(paste(cell, op, title, sep = ' '))+
    scale_size_manual(values = c(3, 3), breaks = c('Significant non nuclear adhesion protein', 'Other'), guide = "none")+
    geom_label_repel(data = all_data %>% filter(non_nuclear_significant_adhesion_protein == "Significant non nuclear adhesion protein"), max.overlaps = Inf, size = 3)+
    guides(fill=guide_legend(title="Significant adhesion proteins", override.aes = aes(label = "")),
           shape = guide_legend(title = 'Nucelar or non nuclear protein'),
           alpha = guide_legend(title = 'Significant protein', override.aes = aes(label = "")))
  #scale_color_manual(values = c('#003f5c', '#ffa600', '#bc5090'), breaks = c('consensus', 'Meta', 'LitCur'), na.value = 'grey')
  
  plot
  
  
  return(list(plot, all_data))
  
}
library(htmlwidgets)
library(plotly)
library(xlsx)


for(cell in c('E34', 'E31')){
  stuff <- volcano_plot(normalized_data, imputed_data = NA, 'Jodi', cell, '')
  
  plot <- stuff[[1]]
  
  saveWidget(ggplotly(plot), file= paste0('non-imputed-final-figures-forthesis/', cell, '_Jodi_uncertainmarkers',".html"))
  ggsave(plot, filename = paste0('non-imputed-final-figures-forthesis/', cell, '_Jodi_uncertainmarkers', '.png'), width = 2500, height = 2000, units = 'px')
  
  data_to_save <- stuff[[2]] %>% 
    filter(log2ratio > 0,
           !is.na(Category),
           pvalue < 0.05) %>% 
    dplyr::select(-Nucleus, -Whole_lysate)
  
  full_set_data <- stuff[[2]] %>% 
    dplyr::select(-Nucleus, -Whole_lysate)
  
  write.csv(x = data_to_save, file =  paste0('non-imputed-final-figures-forthesis/', cell, '_Jodi_uncertainmarkers', '.csv'))
  write.csv(x = full_set_data, file = paste0('non-imputed-final-figures-forthesis/', cell, '_full_set', '.csv'))
}

for(cell in c('E13', 'E21', 'E57', 'E28')){
  stuff <- volcano_plot(normalized_data, completeData, 'Athanasia', cell, '')
  
  plot <- stuff[[1]]
  
  saveWidget(ggplotly(plot), file= paste0('non-imputed-final-figures-forthesis/', cell, '_Athanasia_uncertainmarkers',".html"))
  ggsave(plot, filename = paste0('non-imputed-final-figures-forthesis/', cell, '_Athanasia_uncertainmarkers', '.png'), width = 4000, height = 2000, units = 'px')
  
  # For Ath: remember to filter out the nuclear when you are planning to point out the adhesion proteins on the volcano plots
  data_to_save <- stuff[[2]] %>% 
    filter(log2ratio > 0,
           !is.na(Category),
           pvalue < 0.05) %>% 
    dplyr::select(-Nucleus, -Whole_lysate)
  
  full_set_data <- stuff[[2]] %>% 
    dplyr::select(-Nucleus, -Whole_lysate)
  
  write.csv(x = data_to_save, file =  paste0('non-imputed-final-figures-forthesis/', cell, '_Athanasia_uncertainmarkers', '.csv'))
  write.csv(x = full_set_data, file = paste0('non-imputed-final-figures-forthesis/', cell, '_full_set', '.csv'))
}


# volcanoplot for the thesis
e13_volcano <- volcano_plot(normalized_data, completeData, 'Athanasia', 'E13', '')[[1]]
e21_volcano <- volcano_plot(normalized_data, completeData, 'Athanasia', 'E21', '')[[1]]
e57_volcano <- volcano_plot(normalized_data, completeData, 'Athanasia', 'E57', '')[[1]]
e28_volcano <- volcano_plot(normalized_data, completeData, 'Athanasia', 'E28', '')[[1]]
e31_volcano <- volcano_plot(normalized_data, completeData, 'Jodi', 'E31', '')[[1]]
e34_volcano <- volcano_plot(normalized_data, completeData, 'Jodi', 'E34', '')[[1]]

volcano_for_thesis = ggarrange(e13_volcano, 
                               e28_volcano,
                               e21_volcano, 
                               e57_volcano,  
                               e31_volcano,  
                               e34_volcano,
                               ncol = 2,
                               nrow = 3,
                               common.legend = TRUE,
                               legend = 'right')

volcano_for_thesis

ggsave(volcano_for_thesis, filename = paste0('non-imputed-final-figures-forthesis/', 'volcano_for_thesis.png'), width = 11, height = 11.7)


############################################## Barplots, piecharts, etc #############################################################################
combined_volcano_data <- rbind(e13_volcano$data, e28_volcano$data, e21_volcano$data, e57_volcano$data, e31_volcano$data, e34_volcano$data)

# barplot showing number of proteins in volcano
# Summarize counts
combined_volcano_data_counts2 <- combined_volcano_data %>%
  group_by(cell_line) %>%
  summarise(
    total = n(),
    non_na = sum(!is.na(significant_adhesion))
  )
# Plot
BarplotFinalProteins <- ggplot(combined_volcano_data_counts2, aes(x = reorder(cell_line, -total))) +
  # Background total
  geom_col(aes(y = total, fill = "Total Proteins"), width = 0.7, color = "white") +
  
  # Foreground non-NA portion
  geom_col(aes(y = non_na, fill = "Significant Proteins"), width = 0.7, color = "white") +
  
  # Label: total count on top
  geom_text(aes(y = total, label = total),
            vjust = -0.4, size = 4.5, fontface = "bold") +
  
  # Label: non-NA count inside blue section
  geom_text(aes(y = non_na / 2, label = non_na),
            color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(
    name = NULL,
    values = c("Total Proteins" = "#A8A8A8", "Significant Proteins" = "#FFC107")
  ) +
  labs(
    title = "DIA-MS number of identified proteins",
    subtitle = "Both in whole lysate and nucleus",
    x = "Glioma Stem Cells",
    y = "Number of proteins"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  expand_limits(y = max(combined_volcano_data_counts2$total) * 1.1)

BarplotFinalProteins

ggsave(BarplotFinalProteins, filename = paste0('non-imputed-final-figures-forthesis/', 'BarplotFinalProteins.png'), width = 6, height = 4)

# barplot for whole lysate and nucleus separetley
data <- normalized_data %>% 
  dplyr::select(Accession, Protein.Ids, Protein.Names, Gene.names, First.Protein.Description, ends_with('_n')) %>% 
  #na.omit() %>% 
  pivot_longer(cols = ends_with('_n'), names_to = 'sample', values_to = 'normalized_intensities') %>% 
  mutate(sample = substr(sample, 1, nchar(sample)-2),
         cell_line = substr(sample, 1,3),
         cell_compartment = case_when(str_detect(sample, 'N') ~ 'Nucleus',
                                      TRUE ~ 'Whole_lysate'),
         operator = case_when(str_detect(sample, 'jodi') ~ 'Jodi',
                              TRUE ~ 'Athanasia'),
         replicate = case_when(str_detect(sample, 'WL') ~ substr(sample, 12, 12),
                               TRUE ~ substr(sample, 11,11))) %>% 
  #filter(operator == op, cell_line == cell)  %>% 
  pivot_wider(id_cols = c("Accession", "Gene.names", "cell_line"), names_from = "cell_compartment", values_from = c("normalized_intensities")) %>% 
  rowwise() %>% 
  mutate(nucleus_NA = sum(is.na(Nucleus))>0,
         wl_NA = sum(is.na(Whole_lysate))>0) 

data_counts <- data %>%
  group_by(cell_line) %>%
  summarise(
    nucleus_false = sum(nucleus_NA == FALSE, na.rm = TRUE),
    wl_false = sum(wl_NA == FALSE, na.rm = TRUE)
  ) %>%
  # Convert to long format for easier plotting with ggplot2
  pivot_longer(cols = c(nucleus_false, wl_false), 
               names_to = "variable", values_to = "count")

BarplotWLNuc <- ggplot(data_counts, aes(x = cell_line, y = count, fill = variable)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.8), 
            vjust = -0.3, size = 4) +
  scale_fill_manual(
    values = c("nucleus_false" = "#FF7F0E", "wl_false" = "#1F77B4"),
    labels = c("Number of proteins identified in Nucleus", "Number of proteins identified in Whole Lysate")
  ) +
  labs(
    title = "Number of proteins identified in Whole Lysate and Nucleus",
    x = "Cell Line",
    y = "Number of proteins",
    fill = "Variable"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    #legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  expand_limits(y = max(combined_volcano_data_counts2$total) * 1.1)

BarplotWLNuc

ggsave(BarplotWLNuc, filename = paste0('non-imputed-final-figures-forthesis/', 'BarplotWLNuc.png'), width = 10, height = 5)




# Heatmap 
heatmap_data <- combined_volcano_data %>% 
  select(Accession, Gene.names, cell_line, ratio, log2ratio) %>% #category
  mutate(Gene_Accession = paste(Accession, Gene.names))

ggplot(heatmap_data, aes(x = cell_line, y = Gene.names, fill = log2ratio)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of log2 Fold Changes", fill = "log2FC")

library(tibble)
log2FC_matrix <- heatmap_data %>%
  select(Gene_Accession, cell_line, log2ratio) %>% #Category
  pivot_wider(names_from = cell_line, values_from = log2ratio) %>%
  column_to_rownames("Gene_Accession") %>% 
  drop_na()

#use gene names as row names (not necessary)
rownames(log2FC_matrix) <- make.unique(log2FC_matrix$Gene.names)

log2FC_matrix_num <- log2FC_matrix %>%
  select(E13, E28, E21, E57, E31, E34) %>% 
  mutate(across(where(is.list), ~ sapply(., function(x) as.numeric(x[1])))) %>% 
  #mutate(across(everything(), as.numeric)) %>%   # convert all columns to numeric
  as.matrix()

annotation_row <- log2FC_matrix %>% 
  select(Category)

library(pheatmap)
pheatmap(as.matrix(log2FC_matrix_num),
         color = colorRampPalette(c("blue", "white", "red"))(50), # blue = down, red = up
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         main = "Log2 Fold Change Heatmap")

heatmap <- pheatmap(log2FC_matrix_num,
                    color = colorRampPalette(c("blue", "white", "red"))(100),
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    cutree_rows = 5,   # 🔹 divides rows into 4 clusters
                    cutree_cols = 1,   # 🔹 divides columns into 3 clusters
                    border_color = "gray70",
                    show_rownames = FALSE,
                    main = "Clustered Heatmap of Log2 Fold Changes")

heatmap

pdf("non-imputed-final-figures-forthesis/HeatmapLog2FC.pdf", width = 6, height = 6)
heatmap
dev.off()

#test to add row anotations for meta, lit-cur, and consensus



View(heatmap$tree_row)
row_clusters <- cutree(heatmap$tree_row, k = 6)
cluster_df <- data.frame(
  Protein = names(row_clusters),
  Cluster = factor(row_clusters)
)

# Get row order as plotted in the heatmap
row_order <- heatmap$tree_row$order
rownames(log2FC_matrix_num)[row_order]

# Combine with cluster assignments
cluster_df <- data.frame(
  Protein = rownames(log2FC_matrix_num)[row_order],
  Cluster = row_clusters[rownames(log2FC_matrix_num)[row_order]]
)
head(cluster_df)

#order of clusters top to bottom
unique(cluster_df$Cluster)

library(writexl)
write_xlsx(cluster_df, "non-imputed-final-figures-forthesis/HeatmapDataOrder315624.xlsx")


# pie chart of proportion of adhesion proteins
library(dplyr)
library(ggplot2)

piechart_counts <- combined_volcano_data %>%
  mutate(Category_grouped = ifelse(is.na(Category), "NA", "Non-NA")) %>%
  group_by(cell_line, Category_grouped) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n))

piechart_adhesion <- ggplot(piechart_counts, aes(x = "", y = prop, fill = Category_grouped)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cell_line) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(
    title = "Proportion of Adhesion and Non Adhesion proteins",
    subtitle = "For each cell line",
    fill = "Category"
  ) +
  #theme_minimal(base_size = 14) #+
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_fill_manual(
    values = c(
      "NA" = "#FF5252",       # red
      "Non-NA" = "#990240"   # green    # blue (example if you have more groups)
    ),
    labels = c("Non adhesion protein", "Adhesion protein"),
    name = NULL
  ) +
  geom_text(
    aes(label = paste0(round(prop * 1000)/10, "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  )

piechart_adhesion

ggsave(piechart_adhesion, filename = paste0('non-imputed-final-figures-forthesis/', 'piechart.png'), width = 10, height = 5)

# piechart for proteins with posiive log2 fold change
piechart_counts_positive <- combined_volcano_data %>%
  filter(log2ratio > 0) %>% 
  mutate(Category_grouped = ifelse(is.na(Category), "NA", "Non-NA")) %>%
  group_by(cell_line, Category_grouped) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n))

piechart_adhesion_positive <- ggplot(piechart_counts_positive, aes(x = "", y = prop, fill = Category_grouped)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ cell_line) +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(
    title = "Proportion of Nuclear enriched Adhesion and Non Adhesion proteins",
    subtitle = "For each cell line",
    fill = "Category"
  ) +
  #theme_minimal(base_size = 14) #+
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_fill_manual(
    values = c(
      "NA" = "#265252",       # red
      "Non-NA" = "#9F3240"   # green    # blue (example if you have more groups)
    ),
    labels = c("Non adhesion protein", "Adhesion protein"),
    name = NULL
  ) +
  geom_text(
    aes(label = paste0(round(prop * 1000)/10, "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  )


piechart_adhesion_positive

ggsave(piechart_adhesion_positive, filename = paste0('non-imputed-final-figures-forthesis/', 'piechart_positivelog2fc.png'), width = 10, height = 5)



# finding proteins
final_adhesion_proteins <- combined_volcano_data %>% 
  filter(log2ratio > 0,
         !is.na(significant_adhesion),
         nuclear_proteins == 'Non nuclear')%>%
  mutate(Cell_Type = case_when(
    cell_line %in% c("E57", "E21") ~ "Mesenchymal",
    cell_line %in% c("E31", "E34") ~ "Proneural",
    cell_line %in% c("E13", "E28") ~ "Classical",
    TRUE ~ NA_character_  # default for any unexpected values
  ))

# making venn diagram

# Install if not already installed
# install.packages("ggvenn")
library(ggvenn)

# Make sure Gene.names are unique per Cell_Type
gene_lists <- split(final_adhesion_proteins$Gene.names,final_adhesion_proteins$Cell_Type)
gene_lists <- lapply(gene_lists, unique)

ggvenn(
  gene_lists,
  fill_color = c("red", "green", "blue"),  # circle colors
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 3
)


venn_diagram <- ggvenn(
  gene_lists,
  fill_color = c("#009E73", "#FC7D0B", "#A3ACB9"),  # soft, visually distinct colors
  stroke_size = 0.8,                                  # outline thickness
  set_name_size = 6,                                  # size of set labels
  text_size = 4,                                      # size of overlap numbers
  show_percentage = FALSE                              # show counts instead of percentages
) +
  ggtitle("Adhesion Protein Overlaps Across Cell Types") +        # add a title
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

venn_diagram

ggsave(venn_diagram, filename = paste0('non-imputed-final-figures-forthesis/', 'VennDiagram.png'), width = 10, height = 5)

genes_in_all_three <- Reduce(intersect, gene_lists)
print(genes_in_all_three)

# Get all pairwise intersections
pairwise_intersections <- combn(names(gene_lists), 2, function(x) {
  intersect(gene_lists[[x[1]]], gene_lists[[x[2]]])
}, simplify = FALSE)

names(pairwise_intersections) <- combn(names(gene_lists), 2, paste, collapse = "_")

print(pairwise_intersections)

# present in one cell type
genes_unique <- lapply(names(gene_lists), function(type) {
  others <- setdiff(names(gene_lists), type)
  unique_genes <- setdiff(gene_lists[[type]], unlist(gene_lists[others]))
  return(unique_genes)
})

names(genes_unique) <- names(gene_lists)
print(genes_unique)
















# test pca
normalized_data_test <- normalized_data %>% 
  left_join(uncertain_markers, by = c('Gene.names')) %>% 
  mutate(nuclear_proteins = if_else(grepl('Nuc', Location_Approved), true = "Nuclear", false = Location_Approved)) %>% 
  distinct() %>% 
  mutate(Accession = make.unique(Accession))

rownames(normalized_data_test) <- normalized_data_test$Accession

pca_data <- prcomp(x= na.omit(normalized_data_test[grepl('_n', colnames(normalized_data))]), scale. = TRUE)
View(pca_data$rotation)

pca_plot_data <- data.frame(pca_data$x)

pca_plot_data$Accession <- row.names(pca_plot_data)

pca_plot_data <- pca_plot_data %>% 
  left_join(normalized_data_test %>% dplyr::select(Accession, nuclear_proteins),
            by = 'Accession')

# version without batch correction
#pca_plot_data$names <- rownames(pca_plot_data)
no_correction <- pca_plot_data %>% 
  # mutate(operator = case_when(str_detect(names, 'jodi') ~ 'Jodi',
  #                             TRUE ~ 'Athanasia'),
  #        cell_line = substr(names, 1,3),
  #        cell_compartment = substr(names, 5, 5)) %>% 
  ggplot(aes(x= PC1, y = PC2, color = nuclear_proteins))+
  geom_point()

no_correction

# test pca with ratio
normalized_data_test2 <- normalized_data_2 %>% 
  filter(cell_line != 'E13') %>% 
  pivot_wider(values_from = ratio, names_from = cell_line, id_cols = c('Accession')) %>% 
  na.omit() %>% 
  left_join(normalized_data %>% dplyr::select(Gene.names, Accession),
            by = 'Accession') %>% 
  left_join(uncertain_markers, by = c('Gene.names')) %>% 
  mutate(nuclear_proteins = if_else(grepl('Nuc', Location_Approved), true = "Nuclear", false = Location_Approved)) %>% 
  distinct() %>% 
  mutate(Accession = make.unique(Accession))

normalized_data_test2 <- data.frame(normalized_data_test2)

rownames(normalized_data_test2) <- normalized_data_test2$Accession

pca_data <- prcomp(x= na.omit(normalized_data_test2[grepl('E', colnames(normalized_data_test2))]), scale. = TRUE)
View(pca_data$rotation)

pca_plot_data <- data.frame(pca_data$x)

pca_plot_data$Accession <- row.names(pca_plot_data)

pca_plot_data <- pca_plot_data %>% 
  left_join(normalized_data_test %>% dplyr::select(Accession, nuclear_proteins),
            by = 'Accession')

# version without batch correction
#pca_plot_data$names <- rownames(pca_plot_data)
no_correction <- pca_plot_data %>% 
  # mutate(operator = case_when(str_detect(names, 'jodi') ~ 'Jodi',
  #                             TRUE ~ 'Athanasia'),
  #        cell_line = substr(names, 1,3),
  #        cell_compartment = substr(names, 5, 5)) %>% 
  ggplot(aes(x= PC1, y = PC2, color = nuclear_proteins))+
  geom_point()

no_correction

no_correction <- pca_plot_data %>% 
  mutate(test_label = if_else(nuclear_proteins == 'Nuclear', 'Nuclear', 'Other')) %>% 
  # mutate(operator = case_when(str_detect(names, 'jodi') ~ 'Jodi',
  #                             TRUE ~ 'Athanasia'),
  #        cell_line = substr(names, 1,3),
  #        cell_compartment = substr(names, 5, 5)) %>% 
  ggplot(aes(x= PC1, y = PC2, color = test_label))+
  geom_point()

no_correction

