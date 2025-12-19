# GSEA 
# script for dotplots and ridgeplots

source('Normalisation.R')

library(DOSE)
enrichDGN
library("GSEABase")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)


normalized_data_2 <- normalized_data_2 %>% 
  filter(!is.nan(log2fc),
         !is.na(log2fc))

gapdh<- normalized_data_2 %>% 
  ungroup() %>% 
  filter(Gene.names == 'GAPDH') %>% 
  na.omit() %>% 
  dplyr::select(cell_line, operator, log2fc) %>% 
  dplyr::rename(gapdh_cutoff = log2fc)

normalized_data_2 <- normalized_data_2 %>% 
  left_join(gapdh, by = c('cell_line', 'operator')) %>% 
  filter(log2fc > gapdh_cutoff)

# drop gene.names after semicolons
normalized_data_2$Gene.names <- str_split(normalized_data_2$Gene.names, pattern = ';', simplify = TRUE)[,1]

normalized_data_2$Entrez_id <- mapIds(org.Hs.eg.db, normalized_data_2$Gene.names, 'ENTREZID', 'SYMBOL')

# export data per cell line
library(writexl)

write_xlsx(normalized_data_2 %>% filter(cell_line == 'E13') %>% ungroup() %>% dplyr::select(Gene.names, log2fc), "Dotplot data/E13.xlsx")
write_xlsx(normalized_data_2 %>% filter(cell_line == 'E21') %>% ungroup() %>% dplyr::select(Gene.names, log2fc), "Dotplot data/E21.xlsx")
write_xlsx(normalized_data_2 %>% filter(cell_line == 'E34') %>% ungroup() %>% dplyr::select(Gene.names, log2fc), "Dotplot data/E34.xlsx")
write_xlsx(normalized_data_2 %>% filter(cell_line == 'E31') %>% ungroup() %>% dplyr::select(Gene.names, log2fc), "Dotplot data/E31.xlsx")
write_xlsx(normalized_data_2 %>% filter(cell_line == 'E57') %>% ungroup() %>% dplyr::select(Gene.names, log2fc), "Dotplot data/E57.xlsx")
write_xlsx(normalized_data_2 %>% filter(cell_line == 'E28') %>% ungroup() %>% dplyr::select(Gene.names, log2fc), "Dotplot data/E28.xlsx")

#CHANGE CELL LINE HERE AND DO
cell_lines <- c('E13', 'E21', 'E34', 'E31', 'E57', 'E28')

for(cell_line in cell_lines){
  genelist <- normalized_data_2$log2fc[normalized_data_2$cell_line == cell_line]
  names(genelist) <- normalized_data_2$Entrez_id[normalized_data_2$cell_line == cell_line]
  genelist <- sort(genelist, decreasing = TRUE)
  
  ego_BP <- gseGO(geneList      = genelist,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  minGSSize     = 100,
                  maxGSSize     = 500,
                  eps           = 0,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  verbose       = FALSE)
  
  ego_BP@result <- ego_BP@result %>%
    filter(!grepl('vacuol', Description)) # vacuol to get both vacuole and vacuolar
  
  BP <- dotplot(ego_BP)+ggtitle('BP')
  BP_ridge <- ridgeplot(ego_BP) + ggtitle('BP')+theme(axis.text = element_text(size=10))
  BP_ridge
  
  ego_MF <- gseGO(geneList      = genelist,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  minGSSize     = 100,
                  maxGSSize     = 500,
                  eps           = 0,
                  pAdjustMethod = "BH",   
                  pvalueCutoff  = 0.05,
                  verbose       = FALSE)
  
  ego_MF@result <- ego_MF@result %>% 
    filter(!grepl('vacuol', Description)) # vacuol to get both vacuole and vacuolar
  
  
  
  MF <- dotplot(ego_MF)+ggtitle('MF')
  MF_ridge <- ridgeplot(ego_MF) + ggtitle('MF')
  
  ego_CC <- gseGO(geneList      = genelist,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  minGSSize     = 100,
                  maxGSSize     = 500,
                  eps           = 0,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  verbose       = FALSE)
  
  ego_CC@result <- ego_CC@result %>%
    filter(!grepl('vacuol', Description)) # vacuol to get both vacuole and vacuolar
  
  CC <- dotplot(ego_CC)+ggtitle('CC')
  CC_ridge <- ridgeplot(ego_CC)+ggtitle('CC')
  
  plot1<-ggarrange(BP, MF, CC, ncol = 3, nrow = 1)
  plot2<- ggarrange(BP_ridge, MF_ridge, CC_ridge, ncol = 3, nrow = 1)
  ggsave(paste0(cell_line, '_dotplots.png'), plot1, width = 5000, height = 2000, units = 'px')
  ggsave(paste0(cell_line, '_ridgeplot.png'), plot2, width = 7000, height = 6000, units = 'px')
}





genelist_e13 <- normalized_data_2$log2fc[normalized_data_2$cell_line == 'E13']
names(genelist_e13) <- normalized_data_2$Entrez_id[normalized_data_2$cell_line == 'E13']
genelist_e13 <- sort(genelist_e13, decreasing = TRUE)


ego <- gseGO(geneList      = genelist_e13,
             OrgDb         = org.Hs.eg.db,
             ont           = "CC",
             minGSSize     = 100,
             maxGSSize     = 500,
             eps           = 0,
             pAdjustMethod = "BH",   
             pvalueCutoff  = 0.05,
             verbose       = FALSE)

dotplot(ego)+ggtitle('CC')
ridgeplot(ego)

# 
GOcats4vis <- data.frame(ego)$Description[ sample(1:dim(data.frame(ego))[1], size = 8) ]
dotplot(ego, showCategory=GOcats4vis, orderBy = "setSize")
dotplot(ego, showCategory=GOcats4vis, orderBy = "setSize", split=".sign") +  facet_grid(.~.sign)


# 03-10-2023

#dotplot
edo2 <- gseDO(geneList)
dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

#barplot
library(DOSE)
data(geneList)
ad_gene_list <- read.csv('adhesion.csv')
ad_gene_list <- ad_gene_list %>% 
  select(Gene, Category) %>% 
  rename(name = Gene) %>% 
  mutate(Category = case_when(Category == 'consensus' ~ 'Consensus',
                              Category == 'LitCur' ~ 'Literature-Curated',
                              TRUE ~ Category)) %>% 
  mutate(cl = case_when(Category == 'Meta' ~ 0.1,
                        Category == 'Consensus' ~ 0.4,
                        Category == 'Literature-Curated' ~ 0.7,
                        TRUE ~ 1))
# mutate(cl = case_when(Category == 'Meta' ~ '#7a5195',
#                       Category == 'Consensus' ~ '#ef5675',
#                       Category == 'Literature-Curated' ~ '#ffa600',
#                       TRUE ~ '#000000'))

categorys <- c('Peritoneal adhesion', 'Intrauterine adhesions', 'Leukocyte Adhesion Deficiency Type 3', 'Leukocyte adhesion deficiency type 1')

for(cell_line in cell_lines){
  #print(cell_line)
  genelist <- normalized_data_2$log2fc[normalized_data_2$cell_line == cell_line]
  names(genelist) <- normalized_data_2$Entrez_id[normalized_data_2$cell_line == cell_line]
  genelist <- sort(genelist, decreasing = TRUE)
  
  de <- names(genelist)[abs(genelist) > 2]
  
  edo <- enrichDGN(de)
  
  ## convert gene ID to Symbol
  edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
  
  edox_csv <- edox@result
  write.csv(edox_csv, file = paste0(cell_line, '_edox.csv'))
  
  #p1 <- cnetplot(edox, foldChange=geneList, showCategory = c(edox_csv$Description[1:2], edox_csv$Description[grepl('adh', edox_csv$Description)]))
  
  # nuclear p1
  # change edox_csv$Description[1] to edox_csv$Description[1:2] for original networks with two cluster
  p1 <- cnetplot(edox, foldChange=geneList, showCategory = c(edox_csv$Description[1], edox_csv$Description[grepl('nucl', edox_csv$Description)]))
  
  #ggsave(paste(cell_line,  'Network plot of enriched terms_p1 log2fc.png',  sep = " "), p1, width = 5000, height = 4000, units = 'px')
  
  # p1_data <- p1_data %>% 
  #   left_join(ad_gene_list, by = 'name')
  p1_data <- p1$data
  p1_data$color <- NA
  p1_data$color[p1_data$name %in% ad_gene_list$name[ad_gene_list$Category == 'Meta']] <- 0
  p1_data$color[p1_data$name %in% ad_gene_list$name[ad_gene_list$Category == 'Consensus']] <- 1
  p1_data$color[p1_data$name %in% ad_gene_list$name[ad_gene_list$Category == 'Literature-Curated']] <- -1
  p1_data$color <- factor(p1_data$color)
  p1$data <- p1_data
  # p1_data$color[p1_data$Category == 'Meta'] <- 2
  # p1_data$color[p1_data$Category == 'Literature-Curated'] <- '#bc5090'
  # p1_data$color[p1_data$Category == 'Consensus'] <- '#ffa600'
  # p1_data <- p1_data %>% 
  #   mutate(color = case_when(Category == 'Meta' ~ '#003f5c',
  #                            Category == 'Literature-Curated' ~ '#bc5090',
  #                            Category == 'Consensus' ~ '#ffa600',
  #                            TRUE ~ color))
  
  p1 <- p1 + scale_color_gradient2(low = '#bc5090', mid = '#ffa600', high = '#003f5c', midpoint = 0, na.value = 'grey')
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
  p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
  #output <- cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
  
  ggsave(paste(cell_line,  'Network plot of enriched terms_p1 nuclear.png',  sep = " "), p1, width = 5000, height = 4000, units = 'px')
  ggsave(paste(cell_line, 'Network plot of enriched terms_p2.png', sep = " "), p2, width = 5000, height = 4000, units = 'px')
  ggsave(paste(cell_line, 'Network plot of enriched terms_p3.png', sep = " "), p3, width = 5000, height = 4000, units = 'px')
  
  p1 <- heatplot(edox, showCategory=5)
  p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
  heatplot_output <- cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  
  ggsave(paste(cell_line, 'heatplot.png', sep = " "), heatplot_output, width = 5000, height = 4000, units = 'px')
  
  edox2 <- pairwise_termsim(edox)
  p1 <- treeplot(edox2)
  p2 <- treeplot(edox2, hclust_method = "average")
  treeplot_output <- aplot::plot_list(p1, p2, tag_levels='A')
  
  ggsave(paste(cell_line, 'treeplot.png', sep = " "), treeplot_output, width = 7000, height = 4000, units = 'px')
  
  #enrichment map make letters smaller
  edo <- pairwise_termsim(edo)
  p1 <- emapplot(edo)
  p2 <- emapplot(edo, cex_category=1.5)
  p3 <- emapplot(edo, layout="kk")
  p4 <- emapplot(edo, cex_category=1.5,layout="kk")
  #cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
  
  ggsave(paste(cell_line, 'enrichment_map_p1.png', sep = " "), p1, width = 5000, height = 4000, units = 'px')
  ggsave(paste(cell_line, 'enrichment_map_p2.png', sep = " "), p2, width = 5000, height = 4000, units = 'px')
  ggsave(paste(cell_line, 'enrichment_map_p3.png', sep = " "), p3, width = 5000, height = 4000, units = 'px')
  ggsave(paste(cell_line, 'enrichment_map_p4.png', sep = " "), p4, width = 5000, height = 4000, units = 'px')
  
  upsetp = upsetplot(edo)
  ggsave(paste(cell_line, "upsetplot.png", sep = " "), upsetp, width = 3000, height = 2000, units = "px")
  #articles
  terms <- edo$Description[1:5]
  p <- pmcplot(terms, 2010:2020)
  p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
  articleplot = ggarrange(p, p2, ncol=2)
  
  ggsave(paste(cell_line, "articleplot.png", sep = " "), articleplot, width = 4000, height = 2000, units = "px")
}
de <- names(genelist_e13)[abs(genelist_e13) > 2]

edo <- enrichDGN(de)

barplot(edo, showCategory=20) 
barplot(edo)
mutate(edo, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

#gene concept network size warning error and no gene symbols in plot

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
#output <- cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

ggsave(paste(cell_line,  'Network plot of enriched terms_p1.png',  sep = " "), p1, width = 5000, height = 4000, units = 'px')
ggsave(paste(cell_line, 'Network plot of enriched terms_p2.png', sep = " "), p2, width = 5000, height = 4000, units = 'px')
ggsave(paste(cell_line, 'Network plot of enriched terms_p3.png', sep = " "), p3, width = 5000, height = 4000, units = 'px')


p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
heatplot_output <- cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

ggsave('heatplot.png', heatplot_output, width = 5000, height = 4000, units = 'px')

#treeplot make letters small

edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
treeplot_output <- aplot::plot_list(p1, p2, tag_levels='A')

ggsave('treeplot.png', treeplot_output, width = 7000, height = 4000, units = 'px')

#enrichment map make letters smaller
edo <- pairwise_termsim(edo)
p1 <- emapplot(edo)
p2 <- emapplot(edo, cex_category=1.5)
p3 <- emapplot(edo, layout="kk")
p4 <- emapplot(edo, cex_category=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

ggsave('enrichment_map_p1.png', p1, width = 5000, height = 4000, units = 'px')
ggsave('enrichment_map_p2.png', p2, width = 5000, height = 4000, units = 'px')
ggsave('enrichment_map_p3.png', p3, width = 5000, height = 4000, units = 'px')
ggsave('enrichment_map_p4.png', p4, width = 5000, height = 4000, units = 'px')

