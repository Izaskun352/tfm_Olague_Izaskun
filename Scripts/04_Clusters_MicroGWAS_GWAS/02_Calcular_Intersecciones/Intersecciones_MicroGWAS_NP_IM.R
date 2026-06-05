

# ======================================================================================
# CALCULAR INTERSECCIONES ENTRE CLUSTERS MICROGWAS CON RASGOS NP Y RASGOS IM
# ======================================================================================

# Cargar scripts necesario

source("Scripts/04_Estudio_Clusters/00_Funciones_Calculo_Interseccion.R")
source("Scripts/Renombrar_Clusters.R")

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")
interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
universo_genes <- unique(na.omit(interactoma[,1]))

# OUTPUTS
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
resultado_jaccard_np <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterNervousPsychiatric/Jaccard_Completa_Interseccion_NP_MicroGWAS.csv")

resultado_jaccard_im <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterImmune/Jaccard_Completa_Interseccion_IM_MicroGWAS.csv")
lista_intersecciones_im <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")

# -------------------------------------------------------------------
# RASGOS NP  ----

# 1--- calcular solapamiento con los clusters MicroGWAS

get_nombre_trait <- function(nombre_cluster) {
  trait_id  <- sub("_Cluster_.*$", "", nombre_cluster)
  id_limpio <- sub("^ZSCO\\.", "", trait_id)
  nombre    <- dicc_traits %>%
    dplyr::filter(Rasgo == id_limpio) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  if (length(nombre) == 0 || is.na(nombre)) return(trait_id)
  return(nombre)
}

clusters_Micro_NP <- cargar_clusters(
  carpetas = c(
    "./Output/Piloto_Microbiota/Clusters_MicroGWAS",
    "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric"
  )
)

# 2--- Calcular jaccard

resultados_jaccard_np <- calcular_jaccard_pares(
  lista_clusters   = clusters_Micro_NP,
  carpeta_salida   = "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterNervousPsychiatric",
  nombre_archivo   = "Jaccard_Completa_Interseccion_NP_MicroGWAS.csv",
  fun_nombre_trait = get_nombre_trait
)
#resultado_jaccard_np <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterNervousPsychiatric/Jaccard_Completa_Interseccion_NP_MicroGWAS.csv")

# 3--- Sacamos la lista con las intersecciones

calcular_intersecciones(resultados_jaccard = resultado_jaccard_np,
                        lista_clusters = clusters_Micro_NP,
                        archivo_salida = "./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
#lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")

# 4--- Barplot

pares_mixtos_np <- lapply(names(lista_intersecciones_np), function(nombre_par) {
  inter <- lista_intersecciones_np[[nombre_par]]
  
  data.frame(
    etiqueta               = nombre_par,
    Cluster_1              = inter$c1,
    Cluster_2              = inter$c2,
    Trait_Nombre_1         = inter$trait_1,
    Trait_Nombre_2         = inter$trait_2,
    Indice_Jaccard         = inter$jaccard,
    Genes_Compartidos      = length(inter$genes_ensembl),
    Senal_Inicial_Compartida = sum(inter$genes_tabla$es_semilla_en_alguno),
    Senal_Solo_Cluster1    = sum(inter$genes_tabla$es_semilla_en_c1 & !inter$genes_tabla$es_semilla_en_c2),
    Senal_Solo_Cluster2    = sum(inter$genes_tabla$es_semilla_en_c2 & !inter$genes_tabla$es_semilla_en_c1),
    Senal_Ambos_Clusters   = sum(inter$genes_tabla$es_semilla_en_c1 & inter$genes_tabla$es_semilla_en_c2)
  )
}) %>% dplyr::bind_rows()

barplot_intersecciones(
  pares_mixtos         = pares_mixtos_np,
  colores_traits       = colores_6_traits,
  tabla_traits         = traits_MicroGWAS_areas,
  archivo_salida       = "./Output/Gráficos/MicroGWAS/Barplot_semilla_intersecciones_NP_MicroGWAS.pdf",
  mapeo_clusters_micro = dic_clusters_MicroGWAS,
  mapeo_clusters_np    = dic_clusters_NP,
  fun_acortar          = acortar_nombres_microbioma
)

# -------------------------------------------------------------------

# -------------------------------------------------------------------
# RASGOS IM ----

# 1--- calcular solapamiento con los clusters MicroGWAS

clusters_Micro_IM <- cargar_clusters(
  carpetas = c(
    "./Output/Piloto_Microbiota/Clusters_MicroGWAS",
    "./Output/Piloto_Microbiota/Clusters_Immune"
  )
)

# 2--- Calcular jaccard

resultados_jaccard_im <- calcular_jaccard_pares(
  lista_clusters   = clusters_Micro_IM,
  carpeta_salida   = "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterImmune",
  nombre_archivo   = "Jaccard_Completa_Interseccion_IM_MicroGWAS.csv",
  fun_nombre_trait = get_nombre_trait
)

resultado_jaccard_im <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterImmune/Jaccard_Completa_Interseccion_IM_MicroGWAS.csv")

# 3--- Sacamos la lista con las intersecciones

calcular_intersecciones(resultados_jaccard = resultado_jaccard_im,
                        lista_clusters = clusters_Micro_IM,
                        archivo_salida = "./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")
#lista_intersecciones_im <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")

# 4--- Barplot

pares_mixtos_im <- lapply(names(lista_intersecciones_im), function(nombre_par) {
  inter <- lista_intersecciones_im[[nombre_par]]
  
  data.frame(
    etiqueta               = nombre_par,
    Cluster_1              = inter$c2,
    Cluster_2              = inter$c1,
    Trait_Nombre_1         = inter$trait_2,
    Trait_Nombre_2         = inter$trait_1,
    Indice_Jaccard         = inter$jaccard,
    Genes_Compartidos      = length(inter$genes_ensembl),
    Senal_Inicial_Compartida = sum(inter$genes_tabla$es_semilla_en_alguno),
    Senal_Solo_Cluster1    = sum(inter$genes_tabla$es_semilla_en_c2 & !inter$genes_tabla$es_semilla_en_c1),
    Senal_Solo_Cluster2    = sum(inter$genes_tabla$es_semilla_en_c1 & !inter$genes_tabla$es_semilla_en_c2),
    Senal_Ambos_Clusters   = sum(inter$genes_tabla$es_semilla_en_c1 & inter$genes_tabla$es_semilla_en_c2)
  )
}) %>% dplyr::bind_rows()

barplot_intersecciones(
  pares_mixtos         = pares_mixtos_im,
  colores_traits       = colores_6_traits,
  tabla_traits         = traits_MicroGWAS_areas,
  archivo_salida       = "./Output/Gráficos/MicroGWAS/Barplot_semilla_intersecciones_IM_MicroGWAS.pdf",
  mapeo_clusters_micro = dic_clusters_MicroGWAS,
  mapeo_clusters_np    = dic_clusters_IM,
  fun_acortar          = acortar_nombres_microbioma
)

# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BARPLOT NP + IM
# -------------------------------------------------------------------

pares_mixtos <- dplyr::bind_rows(pares_mixtos_im, pares_mixtos_np)

barplot_intersecciones(
  pares_mixtos         = pares_mixtos,
  colores_traits       = colores_6_traits,
  tabla_traits         = traits_MicroGWAS_areas,
  archivo_salida       = "./Output/Gráficos/MicroGWAS/Barplot_semilla_intersecciones_NP_IM_MicroGWAS.pdf",
  mapeo_clusters_micro = dic_clusters_MicroGWAS,
  mapeo_clusters_np    = c(dic_clusters_NP,dic_clusters_IM),
  fun_acortar          = acortar_nombres_microbioma
)





