# PCA for Lopit experiment

# Libraries ---------------------------------------------------------------
source('Normalisation.R')
library(readxl)
library(tidyr)
library(dplyr)
library("MSnbase")
library("pRoloc")
library("pRolocdata")
library("pRolocGUI")
library(mclust)

#E13 CELL LINE

# Import dataset -----------------------------------------------------------------
df <- readxl::read_xlsx("20230907_nU5_Ex2_Athanasia_E13_24F_proteins.xlsx")

# Normalise dataset -----------------------------------------------------------------
df <- df %>% 
  filter(Contaminant == FALSE) %>% 
  mutate(unique_and_razor = `# Unique Peptides` + `# Razor Peptides`) %>% 
  filter(unique_and_razor > 1,
         !is.na(`Gene Symbol`)) 

#make unique genes MSN  set  does not allow duplicates

df$`Gene Symbol` <- make.unique(df$`Gene Symbol`)

markers <- read.csv('markers.csv')

markers$Location_Approved[grepl('Nuc', markers$Location_Approved)] <- 'Nucleus'
markers <- markers %>% 
  rename(`Ensembl Gene ID` = Gene) %>%
  filter(Location_Approved %in% c("Cytosol", "Focal adhesion sites", "Nucleus", "unknown", "Endoplasmic reticulum", "Mitochondria", "Plasma membrane", "Golgi apparatus", "Actin filaments"))

df <- df %>% 
  left_join(markers, by = c('Ensembl Gene ID')) %>% 
  group_by(Location_Approved) %>% 
  mutate(n = n()) %>% 
  filter(n > 6)

# normalize to total intensity and ratio to the common sample
df_norm <- NormalizeTotalIntensity(df, 'E13') %>% 
  #NormalizeTotalIntensity('Lysate') %>% 
  mutate(across(55:64, ~ .x / WL_Mix_E13_E21_SVG_n))

# change colnames to ratio
colnames(df_norm)[57:66] <- paste0(colnames(df_norm[57:66]), '_ratio')

f0 <- df_norm %>% 
  dplyr::select(Accession, `Gene Symbol`, Location_Approved, 58:66) %>% 
  na.omit() %>% 
  filter_if(~is.numeric(.), all_vars(!is.infinite(.)))

colnames(f0)[4:12] <- str_split(colnames(f0)[4:12], n = 2, pattern = '_n', simplify = TRUE)[,1]

f1 <- f0 %>% 
  ungroup() %>% 
  dplyr::select(`Gene Symbol`, c(4:12)) %>% 
  arrange(desc(`Gene Symbol`))

f2 <- f0 %>% 
  ungroup() %>% 
  dplyr::select(Accession, `Gene Symbol`, Location_Approved) %>% 
  arrange(desc(`Gene Symbol`)) %>% 
  dplyr::select(`Gene Symbol`, Accession, Location_Approved)

write.csv(f1, file = 'f1_best.csv', quote = FALSE, row.names = FALSE)
write.csv(f2, file = 'f2_best.csv', quote = FALSE, row.names = FALSE)

dataset <- readMSnSet(exprsFile = "f1_best.csv",
                      featureDataFile = "f2_best.csv",
                      phenoDataFile = "f3.csv",
                      sep = ","
                      , fill=T, row.names = 1L)

res <- phenoDisco(dataset, fcol = 'Location_Approved')



res1 <- svmClassification(res, p, fcol="Location_Approved")

plot2D(res1, fcol="svm")
addLegend(res1, fcol = "svm",
          where = "bottomright",
          cex = .5)

ptsze <- exp(fData(res1)$svm.scores)
plot2D(res1, fcol = "svm", fpch = "svm",
       cex = ptsze)

pRolocVis(res1)


# PCA of samples (f1) -----------------------------------------------------
library(factoextra)
f1_pca <- data.frame(f1)
row.names(f1_pca) <- f1_pca$Gene.Symbol
f1_pca <- f1_pca[-1]

pca <- prcomp(f1_pca)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

sample_pca <- prcomp(t(f1_pca))

fviz_pca_ind(sample_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
