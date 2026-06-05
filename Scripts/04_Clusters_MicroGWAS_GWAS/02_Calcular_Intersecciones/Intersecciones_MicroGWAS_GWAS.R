

# ======================================================================================
# CALCULAR INTERSECCIONES ENTRE CLUSTERS MICROGWAS CON RASGOS GWAS 
# ======================================================================================

# LIBRERIAS
source("scripts/00_setup.R")

# Cargar scripts necesarios

source("Scripts/04_Estudio_Clusters/00_Funciones_Calculo_Interseccion.R")
source("Scripts/Renombrar_Clusters.R")

# INPUTS

interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
universo_genes <- unique(na.omit(interactoma[,1]))
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

carpeta_clusters_microGWAS <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
carpeta_clusters_VarComun <- "./Output/Piloto_Microbiota/Clusters_Traits_VarComun"

# OUTPUTS
resultado_jaccard_vc <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_MicroGWAS/Jaccard_Completa_VarComun_MicroGWAS.csv")
lista_intersecciones_vc <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_VarComun_MicroGWAS.rds")
lista_intersecciones_vc_vc <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_VarComun_VarComun.rds")

resultados_jaccard_vc_np <- readRDS("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_NervousPsychiatric_VarComun/Pleiotropia_Jaccard_Alta_Significativa.rds")
resultados_jaccard_vc_np <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_NervousPsychiatric_VarComun/Pleiotropia_Jaccard_Alta_Significativa.csv")

resultados_jaccard_vc_im <- readRDS("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_VarComun/Pleiotropia_Jaccard_Alta_Significativa.rds")
resultados_jaccard_vc_im <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_VarComun/Pleiotropia_Jaccard_Alta_Significativa.csv")

resultados_jaccard_im_np <- readRDS("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_NervousPsychiatric/Pleiotropia_Jaccard_Alta_Significativa.rds")
resultados_jaccard_im_np <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_NervousPsychiatric/Pleiotropia_Jaccard_Alta_Significativa.csv")


# FUNCION

get_nombre_trait <- function(nombre_cluster) {
  trait_id  <- sub("_Cluster_.*$", "", nombre_cluster)
  id_limpio <- sub("^ZSCO\\.", "", trait_id)
  nombre    <- dicc_traits %>%
    dplyr::filter(Rasgo == id_limpio) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  if (length(nombre) == 0 || is.na(nombre)) return(trait_id)
  return(nombre)}

# ----------------------------------------------------------------------
# 1- SOLAPAMIENTO MICROGWAS vs GWAS (sin NP / IM) ----

carpeta_pleiotropia <- "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_MicroGWAS"
dir.create(carpeta_pleiotropia, showWarnings = FALSE)

clusters_Micro_VC <- cargar_clusters(
  carpetas = c(carpeta_clusters_microGWAS, carpeta_clusters_VarComun))

resultados_jaccard_vc <- calcular_jaccard_pares(
  lista_clusters   = clusters_Micro_VC,
  carpeta_salida   = carpeta_pleiotropia,
  nombre_archivo   = "Jaccard_Completa_VarComun_MicroGWAS.csv",
  fun_nombre_trait = get_nombre_trait)

# --- Sacamos la lista con las intersecciones

calcular_intersecciones(resultados_jaccard = resultado_jaccard_vc,
                        lista_clusters = clusters_Micro_VC,
                        archivo_salida = "./Output/Piloto_Microbiota/lista_intersecciones_VarComun_MicroGWAS.rds")
lista_intersecciones_vc <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_VarComun_MicroGWAS.rds")

# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# 2- SOLAPAMIENTO GWAS vs GWAS (sin NP/IM) ----

# # Forma más optimizada para calcular jaccard de todos vs todos ----

# 1. Construir "all"
message("Leyendo clusters de VarComun...")
carpeta_VarComun <- "./Output/Piloto_Microbiota/Clusters_Traits_VarComun"
archivos_vc <- list.files(carpeta_VarComun, pattern = "\\.csv$", full.names = TRUE)

all <- do.call(rbind, lapply(archivos_vc, function(f) {
  nombre <- gsub("\\.csv$", "", basename(f))
  message("  Leyendo: ", nombre)
  df <- read.csv2(f)
  data.frame(
    cluster.walktrap = nombre,
    ENSG             = gsub(";.*$", "", df$gene.ENSG),
    stringsAsFactors = FALSE
  )
}))
message("✔️ Total filas en 'all': ", nrow(all), " | Clusters únicos: ", length(unique(all$cluster.walktrap)))

# 2. Construir jacc_tab
message("Generando pares de clusters...")
clusters_vc <- unique(all$cluster.walktrap)
jacc_tab    <- as.data.frame(t(combn(clusters_vc, 2)), stringsAsFactors = FALSE)
colnames(jacc_tab) <- c("g1", "g2")
message("  Pares antes de filtrar mismo trait: ", nrow(jacc_tab))

jacc_tab$trait_g1 <- sub("_Cluster_.*$", "", jacc_tab$g1)
jacc_tab$trait_g2 <- sub("_Cluster_.*$", "", jacc_tab$g2)
jacc_tab <- jacc_tab %>%
  dplyr::filter(trait_g1 != trait_g2) %>%
  dplyr::select(-trait_g1, -trait_g2)
message("✔️ Pares tras eliminar mismo trait: ", nrow(jacc_tab))

# 3. Preparar estructura optimizada
message("Preparando estructura optimizada...")
groups        <- split(all$ENSG, all$cluster.walktrap)
groups        <- lapply(groups, unique)
gene_universe <- unique(all$ENSG)
gene_id       <- seq_along(gene_universe)
names(gene_id) <- gene_universe
int_groups    <- lapply(groups, function(g) sort(gene_id[g]))
group_index   <- seq_along(int_groups)
names(group_index) <- names(int_groups)
pairs_i <- group_index[jacc_tab$g1]
pairs_j <- group_index[jacc_tab$g2]
message("✔️ Estructura lista. Comenzando cálculo de Jaccard...")

# 4. Calcular Jaccard con progress
total_pares <- nrow(jacc_tab)
contador    <- 0L

resultados_jaccard <- mapply(function(i, j) {
  contador <<- contador + 1L
  if (contador %% 1000 == 0 || contador == total_pares) {
    message(sprintf("  Procesando par %d / %d (%.1f%%)...",
                    contador, total_pares, 100 * contador / total_pares))
  }
  a     <- int_groups[[i]]
  b     <- int_groups[[j]]
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) 0 else inter / uni
}, pairs_i, pairs_j)

# 5. Montar tabla resultado
jacc_tab$Indice_Jaccard <- round(resultados_jaccard, 3)
jacc_tab <- jacc_tab %>% dplyr::arrange(desc(Indice_Jaccard))

# 6. Filtrar >= 0.5
jacc_significativa <- jacc_tab %>% dplyr::filter(Indice_Jaccard >= 0.5)
jacc_positivo <- jacc_tab %>% dplyr::filter(Indice_Jaccard > 0)
message("Pares con Jaccard > 0:     ", nrow(dplyr::filter(jacc_tab, Indice_Jaccard > 0)))
message("Pares con Jaccard >= 0.50: ", nrow(jacc_significativa))

# 7. Guardar
carpeta_out <- "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_VarComun"
dir.create(carpeta_out, showWarnings = FALSE)

write.csv2(jacc_positivo,          file.path(carpeta_out, "Matriz_Jaccard_Positiva.csv"),              row.names = FALSE)
write.csv2(jacc_significativa, file.path(carpeta_out, "Pleiotropia_Jaccard_Alta_Significativa.csv"), row.names = FALSE)
saveRDS(jacc_tab,             file.path(carpeta_out, "Matriz_Jaccard_Completa.rds"))
saveRDS(jacc_significativa,   file.path(carpeta_out, "Pleiotropia_Jaccard_Alta_Significativa.rds"))
pares_varcomun_varcomun <- readRDS("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_VarComun/Pleiotropia_Jaccard_Alta_Significativa.rds")

# Creamos lista de intersecciones (no usamos funcion pq estructura de resultados_jaccard es diferente) ----

resultados_jaccard_vc <- read.csv2("./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_VarComun/Pleiotropia_Jaccard_Alta_Significativa.csv")
message("Pares MicroGWAS vs VarComun con Jaccard >= 0.5: ", nrow(resultados_jaccard_vc))

# Cargar clusters con Ensembl IDs
carpeta_clusters <-"./Output/Piloto_Microbiota/Clusters_Traits_VarComun"
archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)

lista_genes_ensembl  <- list()
lista_datos_clusters <- list()

for (f in archivos_clusters) {
  nombre <- gsub("\\.csv$", "", basename(f))
  df     <- read.csv2(f)
  lista_genes_ensembl[[nombre]]  <- gsub(";.*$", "", df$gene.ENSG)
  lista_datos_clusters[[nombre]] <- df
}
head(resultados_jaccard_vc)
# Bucle principal
lista_intersecciones_vc_vc <- list()

for (i in 1:nrow(resultados_jaccard_vc)) {
  
  # Adaptado a tus columnas reales g1 y g2
  c1 <- resultados_jaccard_vc$g1[i]   
  c2 <- resultados_jaccard_vc$g2[i]        
  nombre_par <- paste0(c1, "_vs_", c2)
  
  genes_c1 <- lista_genes_ensembl[[c1]]
  genes_c2 <- lista_genes_ensembl[[c2]]
  
  if (is.null(genes_c1) || is.null(genes_c2)) {
    message("⚠️ No se encontraron genes en las carpetas para: ", nombre_par)
    next
  }
  
  genes_interseccion <- intersect(genes_c1, genes_c2)
  
  if (length(genes_interseccion) < 5) {
    message("⚠️ Menos de 5 genes en la intersección: ", nombre_par)
    next}
  
  datos_c1 <- lista_datos_clusters[[c1]]
  datos_c2 <- lista_datos_clusters[[c2]]
  
  # Eliminados los 'trait_origen' antiguos porque no existen en este dataframe
  filas_c1 <- datos_c1 %>%
    dplyr::filter(gsub(";.*$", "", gene.ENSG) %in% genes_interseccion) %>%
    dplyr::mutate(cluster_origen = c1)
  
  filas_c2 <- datos_c2 %>%
    dplyr::filter(gsub(";.*$", "", gene.ENSG) %in% genes_interseccion) %>%
    dplyr::mutate(cluster_origen = c2)
  
  genes_tabla <- dplyr::bind_rows(filas_c1, filas_c2) %>%
    dplyr::select(gene.ENSG, gene.gene, gene.senal_inicial, cluster_origen) %>%
    dplyr::group_by(gene.gene) %>%
    dplyr::summarise(
      gene.ENSG            = dplyr::first(gene.ENSG),
      es_semilla_en_c1     = any(gene.senal_inicial == "Si" & cluster_origen == c1),
      es_semilla_en_c2     = any(gene.senal_inicial == "Si" & cluster_origen == c2),
      es_semilla_en_alguno = any(gene.senal_inicial == "Si"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      nombre_par    = nombre_par,
      Cluster_g1    = c1,                                                                    
      Cluster_g2    = c2)
  
  simbolos <- unique(genes_tabla$gene.gene)
  
  # Guardado en la lista correcta: lista_intersecciones_vc_vc
  lista_intersecciones_vc_vc[[nombre_par]] <- list(
    c1             = c1,
    c2             = c2,
    jaccard        = resultados_jaccard_vc$Indice_Jaccard[i],
    genes_ensembl  = genes_interseccion,
    genes_simbolo  = simbolos,
    genes_tabla    = genes_tabla
  )
  message("✔️ ", nombre_par, " — ", length(genes_interseccion), " genes en intersección")
}

saveRDS(lista_intersecciones_vc_vc, file = "./Output/Piloto_Microbiota/lista_intersecciones_VarComun_VarComun.rds")

# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# 3- SOLAPAMIENTO GWAS vs NP / IM ----

# Creamos funcion para calcular jaccard de manera mas dirigida y optimizada

calcular_jaccard_dirigido <- function(carpeta_A, carpeta_B, 
                                      nombre_A, nombre_B,
                                      carpeta_output,
                                      umbral = 0.5) {
  
  # 1. Leer archivos
  archivos_A <- list.files(carpeta_A, pattern = "\\.csv$", full.names = TRUE)
  archivos_B <- list.files(carpeta_B, pattern = "\\.csv$", full.names = TRUE)
  message(nombre_A, ": ", length(archivos_A), " clusters")
  message(nombre_B, ": ", length(archivos_B), " clusters")
  
  # 2. Leer grupos
  leer_grupo <- function(archivos) {
    do.call(rbind, lapply(archivos, function(f) {
      nom <- gsub("\\.csv$", "", basename(f))
      df  <- read.csv2(f)
      data.frame(
        cluster = nom,
        ENSG    = gsub(";.*$", "", df$gene.ENSG),
        stringsAsFactors = FALSE
      )
    }))
  }
  
  all_A <- leer_grupo(archivos_A)
  all_B <- leer_grupo(archivos_B)
  
  # 3. Estructura optimizada con enteros
  gene_universe  <- unique(c(all_A$ENSG, all_B$ENSG))
  gene_id        <- seq_along(gene_universe)
  names(gene_id) <- gene_universe
  
  grupos_A <- lapply(split(all_A$ENSG, all_A$cluster), unique)
  grupos_B <- lapply(split(all_B$ENSG, all_B$cluster), unique)
  
  int_A <- lapply(grupos_A, function(g) sort(gene_id[g]))
  int_B <- lapply(grupos_B, function(g) sort(gene_id[g]))
  
  nombres_A <- names(int_A)   # <-- ahora sí está definido después de int_A
  nombres_B <- names(int_B)
  
  total_pares <- length(nombres_A) * length(nombres_B)
  message("Calculando Jaccard para ", total_pares, " pares ", 
          nombre_A, " x ", nombre_B, "...")
  
  # 4. Bucle guardando SOLO pares >= umbral
  contador         <- 0L
  lista_resultados <- list()
  
  for (g1 in nombres_A) {
    for (g2 in nombres_B) {
      
      contador <- contador + 1L
      if (contador %% 10000 == 0 || contador == total_pares) {
        message(sprintf("  Par %d / %d (%.1f%%) | Pares significativos: %d",
                        contador, total_pares,
                        100 * contador / total_pares,
                        length(lista_resultados)))
      }
      
      a     <- int_A[[g1]]
      b     <- int_B[[g2]]
      inter <- length(intersect(a, b))
      if (inter == 0) next
      
      uni     <- length(a) + length(b) - inter
      jaccard <- inter / uni
      
      if (jaccard >= umbral) {
        lista_resultados[[length(lista_resultados) + 1]] <- data.frame(
          g1             = g1,
          g2             = g2,
          Indice_Jaccard = round(jaccard, 3),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # 5. Montar tabla
  if (length(lista_resultados) == 0) {
    message("⚠️ No se encontraron pares con Jaccard >= ", umbral)
    return(invisible(NULL))
  }
  
  pares <- dplyr::bind_rows(lista_resultados)
  colnames(pares) <- c(paste0("Cluster_", nombre_A),
                       paste0("Cluster_", nombre_B),
                       "Indice_Jaccard")
  
  # 6. Añadir nombres con join
  pares <- pares %>%
    dplyr::mutate(
      trait_id_A = sub("^ZSCO\\.", "", sub("_Cluster_.*$", "",
                                           .[[paste0("Cluster_", nombre_A)]])),
      trait_id_B = sub("^ZSCO\\.", "", sub("_Cluster_.*$", "",
                                           .[[paste0("Cluster_", nombre_B)]]))
    ) %>%
    dplyr::left_join(dicc_nombres, by = c("trait_id_A" = "Rasgo")) %>%
    dplyr::rename(!!paste0("Trait_Nombre_", nombre_A) := name) %>%
    dplyr::left_join(dicc_nombres, by = c("trait_id_B" = "Rasgo")) %>%
    dplyr::rename(!!paste0("Trait_Nombre_", nombre_B) := name) %>%
    dplyr::select(-trait_id_A, -trait_id_B) %>%
    dplyr::arrange(desc(Indice_Jaccard))
  
  # 7. Guardar
  dir.create(carpeta_output, showWarnings = FALSE)
  write.csv2(pares,
             file.path(carpeta_output, "Pleiotropia_Jaccard_Alta_Significativa.csv"),
             row.names = FALSE)
  saveRDS(pares,
          file.path(carpeta_output, "Pleiotropia_Jaccard_Alta_Significativa.rds"))
  
  message("✔️ Completado! Pares con Jaccard >= ", umbral, ": ", nrow(pares))
  
  return(invisible(pares))
}

# Ejecutar los 3 cruces
np_vc <- calcular_jaccard_dirigido(
  carpeta_A      = "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric",
  carpeta_B      = "./Output/Piloto_Microbiota/Clusters_Traits_VarComun",
  nombre_A       = "NervousPsychiatric",
  nombre_B       = "VarComun",
  carpeta_output = "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_NervousPsychiatric_VarComun")

immune_vc <- calcular_jaccard_dirigido(
  carpeta_A      = "./Output/Piloto_Microbiota/Clusters_Immune",
  carpeta_B      = "./Output/Piloto_Microbiota/Clusters_Traits_VarComun",
  nombre_A       = "Immune",
  nombre_B       = "VarComun",
  carpeta_output = "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_VarComun")

immune_np <- calcular_jaccard_dirigido(
  carpeta_A      = "./Output/Piloto_Microbiota/Clusters_Immune",
  carpeta_B      = "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric",
  nombre_A       = "Immune",
  nombre_B       = "NervousPsychiatric",
  carpeta_output = "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_NervousPsychiatric")

# ----------------------------------------------------------------------