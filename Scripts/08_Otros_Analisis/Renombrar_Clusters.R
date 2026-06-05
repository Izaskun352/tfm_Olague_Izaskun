

library(tidyverse)

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# FUNCION PARA CREAR DICCIONARIO Y RENOMBRAR CLUSTERS

crear_mapeo_clusters <- function(carpeta) {
  
  archivos <- list.files(carpeta, pattern = "\\.csv$", full.names = FALSE)
  nombres_originales <- tools::file_path_sans_ext(archivos)
  
  df <- tibble(original = nombres_originales) %>%
    mutate(
      # Rasgo = todo antes de "_Cluster_"
      rasgo = str_extract(original, "^(.+?)(?=_Cluster_)"),
      # Número de cluster = todo después de "_Cluster_"
      num_cluster = str_extract(original, "(?<=_Cluster_)[\\d\\.]+$")
    )
  
  # Función para ordenación jerárquica: 1 < 1.1 < 1.1.7 < 1.2 < 2
  ordenar_jerarquico <- function(x) {
    partes <- str_split(x, "\\.")[[1]]
    paste(str_pad(partes, 4, pad = "0"), collapse = ".")
  }
  
  df <- df %>%
    mutate(orden = sapply(num_cluster, ordenar_jerarquico)) %>%
    arrange(rasgo, orden) %>%
    group_by(rasgo) %>%
    mutate(
      indice = row_number(),                                    # reinicia en 1 por rasgo
      nuevo_nombre = paste0(rasgo, "_Cluster_", indice)
    ) %>%
    ungroup()
  
  mapeo <- setNames(df$nuevo_nombre, df$original)
  
  return(mapeo)
}

# CLUSTERS MICROGWAS

carpeta_clusters_microGWAS <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
dic_clusters_MicroGWAS <- crear_mapeo_clusters(carpeta_clusters_microGWAS)

# CLUSTERS NP

carpeta_clusters_NP <- "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric"
dic_clusters_NP <- crear_mapeo_clusters(carpeta_clusters_NP)

# CLUSTERS IM

carpeta_clusters_IM <- "./Output/Piloto_Microbiota/Clusters_Immune"
dic_clusters_IM <- crear_mapeo_clusters(carpeta_clusters_IM)

# DICCIONARIO UNIFICADO
cluster_names <- c(dic_clusters_MicroGWAS, dic_clusters_NP, dic_clusters_IM)


# ESTANDARIZAR NOMBRES DE UNA TABLA
estandarizar_tabla_NP <- function(tabla) {
  tabla %>%
    dplyr::rename(
      Cluster_MicroGWAS  = Cluster_1,       # ajusta si tus nombres finales son distintos
      Cluster_NP_IM      = Cluster_2,
      Trait_Micro_Nombre = Trait_Nombre_1,
      Trait_NP_IM_Nombre = Trait_Nombre_2
    )
}

estandarizar_tabla_IM <- function(tabla) {
  tabla %>%
    dplyr::rename(
      Cluster_NP_IM      = Cluster_IM,
      Trait_NP_IM_Nombre = Trait_IM_Nombre
    )
}

# ── Colores MicroGWAS ya establecidos ----
acortar_nombres_microbioma <- function(nombres) {
  n <- gsub(" microbiome measurement", " Micr.", nombres, ignore.case = TRUE)
  return(n)}

colores_6_traits <- c(
  "EFO_0007753" = "#FFF2AE",
  "EFO_0007874" = "pink",
  "EFO_0007883" = "#BFEFFF",
  "EFO_0011013" = "#FDBF6F",
  "EFO_0801228" = "#D2B4DE",
  "EFO_0801229" = "#C1E1C1"
)

efo_a_nombre <- traits_MicroGWAS_areas %>%
  dplyr::filter(Rasgo %in% names(colores_6_traits)) %>%
  dplyr::select(Rasgo, name) %>%
  dplyr::distinct() %>%
  tibble::deframe()

colores_micro_por_nombre <- setNames(
  colores_6_traits[names(efo_a_nombre)],
  acortar_nombres_microbioma(efo_a_nombre)   # mismo acortamiento que se aplica en Par
)


# --- Colores Traits NP - IM

todos_traits_NP_IM <- c(
  "atopic eczema", "psoriasis", "inflammatory bowel disease", 
  "ankylosing spondylitis", "seasonal allergic rhinitis", 
  "systemic lupus erythematosus", "attention deficit hyperactivity disorder", 
  "obsessive-compulsive disorder", "major depressive disorder", 
  "substance abuse", "bipolar disorder", "anorexia nervosa", 
  "peripheral neuropathy", "opioid dependence", 
  "autism spectrum disorder", "insomnia"
)

paleta_np_im <- c(
  "orangered", "palegreen4", "purple3", "peru", "peachpuff4", "red2", 
  "royalblue", "seagreen1", "seashell2", "sienna1", "#67001F", 
  "#142157", "#E08214", "palevioletred1", "#ffff00", "slategray"
)

# Crear el mapeo fijo
colores_traits_NP_IM <- setNames(paleta_np_im, todos_traits_NP_IM)



