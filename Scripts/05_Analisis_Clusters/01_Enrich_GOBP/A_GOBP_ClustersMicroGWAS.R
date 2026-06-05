
# Cargar scripr necesario

source("Scripts/04_Estudio_Clusters/03_Enrich_GOBP/0_Funciones_GOBP.R")
library(clusterProfiler)
library(org.Hs.eg.db)

# INPUTS
interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
universo_genes <- unique(na.omit(interactoma[,1]))
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# CALCULAR ENRIQUECIMIENTO GOBP DE LOS CLUSTERS MICROGWAS

analizar_GO_carpeta(
  carpeta_clusters = "./Output/Piloto_Microbiota/Clusters_MicroGWAS",
  carpeta_salida   = "./Output/Piloto_Microbiota/GO_Clusters_MicroGWAS",
  universo_genes   = universo_genes,
  tabla_traits     = traits_MicroGWAS_areas
)

# --- REPRESENTAR HEATMAP

# Cargamos scrips con funciones para cambiar nombre clusters + colores establecidos

source("Scripts/Renombrar_Clusters.R")

heatmap_GO_clusters(
  archivo_maestro = "./Output/Piloto_Microbiota/GO_Clusters_MicroGWAS/GO_Resumen_Simplificado.csv",
  archivo_salida  = "./Output/Gráficos/MicroGWAS/Heatmap_Enrich_GOBP_ClustersMicroGWAS.pdf",
  tabla_traits    = traits_MicroGWAS_areas,
  colores_traits  = colores_6_traits,
  n_top           = 3,
  fun_acortar     = acortar_nombres_microbioma,
  mapeo_clusters  = dic_clusters_MicroGWAS,
  titulo = "GO Biological Process enrichment by cluster (-log10 p.adjust)"
  )

# --- REPRESENTAR DOTPLOT

dotplot_GO_clusters(
  go_maestro     = "./Output/Piloto_Microbiota/GO_Clusters_MicroGWAS/GO_Resumen_Simplificado.csv",
  archivo_salida = "./Output/Gráficos/MicroGwas/Dotplot_Enrich_GOBP_ClustersMicroGWAS.pdf",
  tabla_traits   = traits_MicroGWAS_areas,
  n_top          = 2,
  fun_acortar    = acortar_nombres_microbioma,
  mapeo_clusters = dic_clusters_MicroGWAS,
  titulo = "GO-BP enrichment by MicroGWAS Cluster",
  ancho = 23,
  alto = 18
)




