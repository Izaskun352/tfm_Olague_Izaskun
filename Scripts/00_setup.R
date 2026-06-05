
# Cargar librerias


library(tidyverse)
library(vegan)
library(reshape2) # pasar la matriz a formato de columnas
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(circlize)
library(RColorBrewer)   # paleta específica: colorRampPalette(brewer.pal(9, "YlGnBu"))
library(igraph)
library(forcats)
library(ggVennDiagram)
library(rlang)
library(stringr)
library(tibble)
library(pals)
library(arrow)
library(stringr) # Para manipular el texto de los IDs fácilmente
library(purrr)   # Para manejar la columna de tipo lista
library(qgraph)
library(tools)
library(rvest)


# Definir opciones globañes

# Definir funciones auxiliares (opcional)


#Distancia_Traits_GWAS_RARE_Comun <- readRDS("./Data./Matrices_Distancias./Distancia_Traits_GWAS_RARE_Comun.rds")  
#Similitud_Traits_GWAS_RARE_Comun <- readRDS("./Data./Matrices_Distancias./Similitud_Traits_GWAS_RARE_Comun.rds")
#matriz_indice_jaccard_comun <- readRDS("./Data./Matrices_Distancias/Matriz_indice_jaccard_comun.rds")

#Distancia_Traits_GWAS_RARE_Rara <- readRDS("./Data./Matrices_Distancias./Distancia_Traits_GWAS_RARE_Rara.rds")
#Similitud_Traits_GWAS_RARE_Rara <- readRDS("./Data./Matrices_Distancias./Similitud_Traits_GWAS_RARE_Rara.rds")
#matriz_indice_jaccard_rara <- readRDS("./Data./Matrices_Distancias/Matriz_indice_jaccard_rara.rds")




