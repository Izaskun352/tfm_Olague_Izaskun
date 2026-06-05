
source("scripts/00_setup.R")  ## Abrimos script con las librerias 
library(clusterProfiler)

# =========================================================================
# ANÁLISIS DE EXPRESIÓN EN TEJIDOS (Gene Ontology)  

# INTERSECCION MICORGWAS - NERVOUS / PSYCHIATRIC
# =========================================================================

# ----------------------------------------------------------------------------------------------------
# INPUTS
# ----------------------------------------------------------------------------------------------------

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
genes_interactoma <- interactoma[,1]
message("Genes en el interactoma: ", length(genes_interactoma))
# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# FUNCIONES
# ----------------------------------------------------------------------------------------------------

heatmap_expresion <- function(datos_cl, tipo = "log2") {
  
  trait1      <- unique(datos_cl$trait_1)
  trait2      <- unique(datos_cl$trait_2)
  jaccard_val <- unique(datos_cl$jaccard)
  
  # --- Construir matriz log2 ---
  matriz_base <- datos_cl %>%
    dplyr::select(Gene.name, Tissue, log2_nTPM) %>%
    tidyr::pivot_wider(names_from  = Tissue,
                       values_from = log2_nTPM) %>%
    tibble::column_to_rownames("Gene.name") %>%
    as.matrix()
  
  # --- Aplicar transformación según tipo ---
  if (tipo == "log2") {
    matriz_plot <- matriz_base
    titulo_tipo <- "[log2 nTPM]"
    
  } else if (tipo == "zscore") {
    matriz_plot <- t(apply(matriz_base, 1, zscore_robusto))
    colnames(matriz_plot) <- colnames(matriz_base)
    titulo_tipo <- "[Z-score robusto]"
  }
  
  paleta      <- colorRampPalette(c("#EDF8B1", "#7FCDBB", "#2C7FB8", "#081D58"))(100)
  # --- Distancias (siempre desde matriz sin NAs) ---
  matriz_calculo <- matriz_plot
  matriz_calculo[is.na(matriz_calculo)] <- 0
  dist_tejidos <- as.dist(1 - cor(matriz_calculo, method = "spearman"))
  dist_genes   <- dist(matriz_calculo, method = "euclidean")
  
  # --- Tamaños dinámicos ---
  n_genes <- nrow(matriz_plot)
  cellheight_din   <- dplyr::case_when(
    n_genes <= 20  ~ 20, n_genes <= 50  ~ 12,
    n_genes <= 100 ~ 8,  n_genes <= 200 ~ 5,
    TRUE           ~ 3
  )
  fontsize_row_din <- dplyr::case_when(
    n_genes <= 20  ~ 9, n_genes <= 50  ~ 7,
    n_genes <= 100 ~ 6, n_genes <= 200 ~ 5,
    TRUE           ~ 4
  )
  
  # --- Heatmap ---
  tryCatch({
    pheatmap::pheatmap(
      mat                      = matriz_plot,
      color                    = paleta,
      cluster_rows             = TRUE,
      cluster_cols             = TRUE,
      clustering_distance_rows = dist_genes,
      clustering_distance_cols = dist_tejidos,
      clustering_method        = "average",
      treeheight_row           = 0,
      treeheight_col           = 30,
      show_rownames            = TRUE,
      show_colnames            = TRUE,
      fontsize_col             = 8,
      fontsize_row             = fontsize_row_din,
      angle_col                = 45,
      cellwidth                = 15,
      cellheight               = cellheight_din,
      main                     = paste0(titulo_tipo, " ", trait1, " vs ", trait2,
                                        "\n(Jaccard = ", jaccard_val,
                                        " | ", n_genes, " genes | varianza > 1)"),
      na_col                   = "whitesmoke",
      border_color             = "#DCDCDC"
    )
  }, error = function(e) {
    message("  ⚠️ Error: ", e$message)
  })
}


# ----------------------------------------------------------------------------------------------------
# BLOQUE 1- ==== DESCARGA DEL DATASET DE EXPRESIÓN  HPA ====

# --- Descomprimir el zip descargado manualmente   -- de HPA -->  rna_tissue_consensus.tsv.zip: ---
unzip("./Data/Diccionarios/rna_tissue_consensus.tsv.zip", 
      exdir = "./Data/Diccionarios")

# --- Leer y explorar ---
hpa_data <- read.delim("./Data/Diccionarios/rna_tissue_consensus.tsv")

message("Tejidos únicos: ", n_distinct(hpa_data$Tissue))
message("Genes únicos: ", n_distinct(hpa_data$Gene))

# Transformar nTPM a log2 poniendo NA donde nTPM == 0
hpa_data <- hpa_data %>%
  dplyr::mutate(
    log2_nTPM = ifelse(nTPM == 0, NA, log2(nTPM))
  )

# --- Filtrar solo genes del interactoma  ---
hpa_universo <- hpa_data %>%
  dplyr::filter(Gene %in% genes_interactoma)

tejidos <- unique(hpa_universo$Tissue)
message("Tejidos a analizar: ", length(tejidos))

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 2: ===== Expresion global del cluster por tejido =====

resultados_jaccard_np <- read.csv2(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterNervousPsychiatric/Matriz_Jaccard_Completa.csv"
)


# Identificar pares mixtos: un MicroGWAS y un Nervous/Psychiatric
pares_mixtos <- resultados_jaccard_np %>%
  dplyr::filter(Indice_Jaccard >= 0.5) %>%
  dplyr::filter(
    # Cluster_1 es MicroGWAS y Cluster_2 es Nervous/Psychiatric
    (!grepl("^ZSCO\\.", Cluster_1) & grepl("^ZSCO\\.", Cluster_2)) |
      # O al revés
      (grepl("^ZSCO\\.", Cluster_1) & !grepl("^ZSCO\\.", Cluster_2))
  )

message("Pares mixtos MicroGWAS vs Nervous/Psychiatric: ", nrow(pares_mixtos))

pares_mixtos %>%
  dplyr::select(Cluster_1, Cluster_2, Indice_Jaccard, Genes_Compartidos) %>%
  print()

#  Recargar los clusters con Ensembl IDs
carpeta_clusters <- c("./Output/Piloto_Microbiota/Clusters_MicroGWAS", "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric")
archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)

lista_genes_ensembl <- lapply(archivos_clusters, function(f) {
  df <- read.csv2(f)
  gsub(";.*$", "", df$gene.ENSG)  # Ensembl IDs limpios
})
names(lista_genes_ensembl) <- gsub("\\.csv$", "", basename(archivos_clusters))

#

# ---- Función: test KS para un cluster en un tejido ----
test_ks_tejido <- function(genes_cluster, genes_universo, hpa_data, tejido) {
  
  # Expresión en este tejido
  expr_tejido <- hpa_data %>%
    dplyr::filter(Tissue == tejido)
  
  # log2_nTPM de los genes del cluster en este tejido (sin NAs)
  expr_cluster <- expr_tejido %>%
    dplyr::filter(Gene %in% genes_cluster) %>%
    dplyr::pull(log2_nTPM) %>%
    na.omit()
  
  # log2_nTPM del resto del universo en este tejido (sin NAs)
  expr_resto <- expr_tejido %>%
    dplyr::filter(!Gene %in% genes_cluster,
                  Gene %in% genes_universo) %>%
    dplyr::pull(log2_nTPM) %>%
    na.omit()
  
  # Necesitamos suficientes genes en ambos grupos
  if (length(expr_cluster) < 3 | length(expr_resto) < 3) return(NULL)
  
  # Test KS de dos colas
  ks <- ks.test(expr_cluster, expr_resto, alternative = "two.sided")
  
  # Dirección del enriquecimiento: media cluster vs media resto
  direccion <- ifelse(mean(expr_cluster) > mean(expr_resto), "higher", "lower")
  
  return(data.frame(
    tissue         = tejido,
    n_genes_cluster = length(expr_cluster),
    n_genes_resto  = length(expr_resto),
    media_cluster  = mean(expr_cluster),
    media_resto    = mean(expr_resto),
    fold_change    = mean(expr_cluster) / mean(expr_resto),
    direccion      = direccion,
    statistic      = ks$statistic,
    pvalue         = ks$p.value
  ))
}

# ---- Bucle para calcular expresion de todos los clusters en todos los tejidos ----

carpeta_expresion_pleiotropia <- "./Output/Piloto_Microbiota/Expresion_Pleiotropia_MicroGWAS_NervousPsychiatric"
dir.create(carpeta_expresion_pleiotropia, showWarnings = FALSE)

lista_resultados_pleiotropia_tejidos <- list()

for (i in 1:nrow(pares_mixtos)) {
  
  c1 <- pares_mixtos$Cluster_1[i]
  c2 <- pares_mixtos$Cluster_2[i]
  jaccard <- pares_mixtos$Indice_Jaccard[i]
  
  message("=== Par ", i, ": ", c1, " vs ", c2, " (Jaccard = ", jaccard, ") ===")
  
  # Extraer genes Ensembl de cada cluster
  genes_c1 <- lista_genes_ensembl[[c1]]
  genes_c2 <- lista_genes_ensembl[[c2]]
  
  if (is.null(genes_c1) || is.null(genes_c2)) {
    message("  ⚠️ No se encontraron genes para alguno de los clusters")
    next
  }
  
  # Genes en la intersección
  genes_interseccion <- intersect(genes_c1, genes_c2)
  message("  Genes en intersección: ", length(genes_interseccion))
  
  if (length(genes_interseccion) < 5) {
    message("  ⚠️ Menos de 5 genes en la intersección, saltando...")
    next
  }
  # --- 2. APLICAR TEST KS A LOS GENES COMPARTIDOS ---
  
  # Aplicar test KS para cada tejido usando SOLO los genes de la intersección
  res_tejidos <- lapply(tejidos, function(tej) {
    test_ks_tejido(
      genes_cluster  = genes_interseccion, # ¡El cambio clave está aquí!
      genes_universo = genes_interactoma,
      hpa_data       = hpa_universo,
      tejido         = tej
    )
  })
  
  # Combinar resultados y limpiar nulos (por si algún tejido no tenía datos)
  res_df <- dplyr::bind_rows(res_tejidos)
  
  # Si el test se pudo hacer al menos en algún tejido, ajustamos p-valores
  if (nrow(res_df) > 0) {
    res_df <- res_df %>%
      dplyr::mutate(
        Cluster_1 = c1,
        Cluster_2 = c2,
        Par = paste0(c1, "_vs_", c2),
        Trait_Nombre_1 = pares_mixtos$Trait_Nombre_1[i], # Opcional: añadir los nombres legibles
        Trait_Nombre_2 = pares_mixtos$Trait_Nombre_2[i],
        p.adjust = p.adjust(pvalue, method = "BH")
      ) %>%
      dplyr::arrange(p.adjust)
    
    message("  ✔️ Tejidos significativos (p.adj < 0.05, higher): ",
            sum(res_df$p.adjust < 0.05 & res_df$direccion == "higher"))
    
    # Guardar CSV individual
    nombre_par <- paste0(c1, "_vs_", c2)
    write.csv2(res_df,
               file = file.path(carpeta_expresion_pleiotropia,
                                paste0("Expresion_Pleiotropia_", nombre_par, ".csv")),
               row.names = FALSE)
    
    # Añadir a la lista para luego hacer una tabla global
    lista_resultados_pleiotropia_tejidos[[nombre_par]] <- res_df
    
  } else {
    message("  ⚠️ No hubo suficientes datos de expresión para evaluar los tejidos en este par.")
  }
}

# --- 3. CREAR TABLA MAESTRA DE RESULTADOS ---
# Unimos todos los resultados en un solo dataframe para facilitar la graficación posterior
if(length(lista_resultados_pleiotropia_tejidos) > 0) {
  tabla_maestra_expresion <- dplyr::bind_rows(lista_resultados_pleiotropia_tejidos)
  
  write.csv2(tabla_maestra_expresion,
             file = file.path(carpeta_expresion_pleiotropia, "Tabla_Maestra_Expresion_Pleiotropia.csv"),
             row.names = FALSE)
  message("\n¡Análisis completado! Tabla maestra guardada.")
}

# ---- Exportar resultados simplificados ====

carpeta_resultados_pleio <- "./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_NervousPsychiatric_MicroGWAS"
dir.create(carpeta_resultados_pleio, showWarnings = FALSE, recursive = TRUE)

# --- 1. Resumen global (Solo los pares que tienen algún tejido significativo) ---
resumen_ks_pleio <- tabla_maestra_expresion %>%
  dplyr::filter(p.adjust < 0.05, direccion == "higher") %>%
  dplyr::group_by(Par, Trait_Nombre_1, Trait_Nombre_2) %>% 
  dplyr::summarise(
    n_tejidos_sig  = n(),
    tejidos        = paste(tissue, collapse = ", "),
    mejor_p.adjust = min(p.adjust),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mejor_p.adjust)

write.csv2(resumen_ks_pleio,
           file = file.path(carpeta_resultados_pleio, "resumen_global_KS_Pleiotropia.csv"),
           row.names = FALSE)

# --- 2. Clasificar cada Par en una categoría interpretativa ---
# Extraemos todos los pares que llegaron a evaluarse (los que tenían >= 5 genes)
pares_todos_df <- pares_mixtos %>%
  dplyr::mutate(Par = paste0(Cluster_1, "_vs_", Cluster_2)) %>%
  dplyr::select(Par, Trait_Nombre_1, Trait_Nombre_2) %>%
  dplyr::distinct()

# Identificamos cuáles sí llegaron a evaluarse en la tabla maestra
pares_evaluados <- unique(tabla_maestra_expresion$Par)

resumen_interpretado_pleio <- pares_todos_df %>%
  dplyr::left_join(resumen_ks_pleio, by = c("Par", "Trait_Nombre_1", "Trait_Nombre_2")) %>%
  dplyr::mutate(
    categoria = dplyr::case_when(
      # Si el par no está en los evaluados, es que se saltó por falta de genes
      !Par %in% pares_evaluados ~ "Skipped (< 5 shared genes or no data)",
      !is.na(n_tejidos_sig) & n_tejidos_sig >= 50 ~ "Ubiquitous expression",
      !is.na(n_tejidos_sig) & n_tejidos_sig > 0   ~ "Tissue-specific enrichment",
      TRUE                                        ~ "No significant enrichment / Lower"
    ),
    n_tejidos_sig  = ifelse(is.na(n_tejidos_sig), 0, n_tejidos_sig),
    tejidos        = ifelse(is.na(tejidos), "-", tejidos),
    mejor_p.adjust = ifelse(is.na(mejor_p.adjust), NA, mejor_p.adjust)
  ) %>%
  dplyr::arrange(categoria, desc(n_tejidos_sig))

print(head(resumen_interpretado_pleio))

write.csv2(resumen_interpretado_pleio,
           file = file.path(carpeta_resultados_pleio, "resumen_interpretado_todos_pares_pleiotropia.csv"),
           row.names = FALSE)

message("\n¡Resúmenes generados y clasificados con éxito! Se incluyeron todas las intersecciones originales.")


# ---- VISUALIZACIÓN DE RESULTADOS ----

acortar_nombres_microbioma <- function(nombres) {
  n <- gsub(" microbiome measurement", " Micr.", nombres, ignore.case = TRUE)
  return(n)
}
# ── Colores MicroGWAS ya establecidos 

colores_6_traits <- c(
  "EFO_0007753" = "#FFFF99",
  "EFO_0007874" = "#FFDAC1",
  "EFO_0007883" = "#AEC6CF",
  "EFO_0011013" = "#B5EAD7",
  "EFO_0801228" = "#C9B1FF",
  "EFO_0801229" = "#FFD1DC"
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
tabla_maestra_expresion <- read.csv2(
  "./Output/Piloto_Microbiota/Expresion_Pleiotropia_MicroGWAS_NervousPsychiatric/Tabla_Maestra_Expresion_Pleiotropia.csv"
)

# =====> BOXPLOT EXPRESIÓN CLUSTER POR TEJIDO ====

# Identificar pares mixtos: un MicroGWAS y un Nervous/Psychiatric (esto está ya en el Bloque 2, pero para ejecutar este bloque solo directamente)

resultados_jaccard_np <- read.csv2(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterNervousPsychiatric/Matriz_Jaccard_Completa.csv"
)
pares_mixtos <- resultados_jaccard_np %>%
  dplyr::filter(Indice_Jaccard >= 0.5) %>%
  dplyr::filter(
    # Cluster_1 es MicroGWAS y Cluster_2 es Nervous/Psychiatric
    (!grepl("^ZSCO\\.", Cluster_1) & grepl("^ZSCO\\.", Cluster_2)) |
      # O al revés
      (grepl("^ZSCO\\.", Cluster_1) & !grepl("^ZSCO\\.", Cluster_2))
  )

message("Pares mixtos MicroGWAS vs Nervous/Psychiatric: ", nrow(pares_mixtos))

pares_mixtos %>%
  dplyr::select(Cluster_1, Cluster_2, Indice_Jaccard, Genes_Compartidos) %>%
  print()

#  Recargar los clusters con Ensembl IDs  (esto está ya en el Bloque 2, pero para ejecutar este bloque solo directamente)
carpeta_clusters <- c("./Output/Piloto_Microbiota/Clusters_MicroGWAS", "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric")
archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)

lista_genes_ensembl <- lapply(archivos_clusters, function(f) {
  df <- read.csv2(f)
  gsub(";.*$", "", df$gene.ENSG)  # Ensembl IDs limpios
})
names(lista_genes_ensembl) <- gsub("\\.csv$", "", basename(archivos_clusters))

# Crear la carpeta

carpeta_boxplots_pleio <- "./Output/Piloto_Microbiota/Boxplots_Expresion_Intersecciones_NervousPsychiatric_MicroGWAS"
dir.create(carpeta_boxplots_pleio, showWarnings = FALSE)

# --- 1. Construir dataframe con genes de cada intersección ---
# Para cada par mixto, extraer genes de la intersección y etiquetarlos

lookup_intersecciones <- lapply(1:nrow(pares_mixtos), function(i) {
  c1 <- pares_mixtos$Cluster_1[i]
  c2 <- pares_mixtos$Cluster_2[i]
  
  genes_c1 <- lista_genes_ensembl[[c1]]
  genes_c2 <- lista_genes_ensembl[[c2]]
  
  if (is.null(genes_c1) || is.null(genes_c2)) return(NULL)
  
  genes_int <- intersect(genes_c1, genes_c2)
  if (length(genes_int) < 5) return(NULL)
  
  # Extraer trait MicroGWAS (el que no tiene ZSCO)
  trait_micro <- ifelse(!grepl("^ZSCO\\.", c1), c1, c2)
  trait_np    <- ifelse(grepl("^ZSCO\\.", c1), c1, c2)
  # Nombre legible del trait MicroGWAS
  trait_micro_id <- gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", trait_micro)
  
  # Quitar prefijo ZSCO. y extraer ID del trait NP
  trait_np_id <- gsub("^ZSCO\\.", "", trait_np)           # quitar ZSCO.
  trait_np_id <- gsub("^([^_]+_[0-9]+)_Cluster.*$", "\\1", trait_np_id)  # extraer ID
  
  # Nombres legibles
  trait_micro_nombre <- traits_MicroGWAS_areas %>%
    dplyr::filter(Rasgo == trait_micro_id) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  
  trait_np_nombre <- traits_MicroGWAS_areas %>%
    dplyr::filter(Rasgo == trait_np_id) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  # Si no encuentra el nombre NP usar el ID como fallback
  if (is.na(trait_np_nombre) || length(trait_np_nombre) == 0) {
    trait_np_nombre <- trait_np_id
  }
  
  # Nombre corto del par para la etiqueta
  cluster_micro_corto <- gsub("^EFO_[0-9]+_", "", trait_micro)
  cluster_np_corto    <- gsub("^ZSCO\\.[^_]+_", "", trait_np)
  
  etiqueta <- paste0(cluster_micro_corto, " vs ", cluster_np_corto,
                     "\n(", trait_np_nombre, ")")
  
  data.frame(
    Gene               = genes_int,
    etiqueta           = etiqueta,
    trait_micro_id     = trait_micro_id,
    trait_micro_nombre = trait_micro_nombre,
    trait_np_nombre    = trait_np_nombre,
    Par                = paste0(c1, "_vs_", c2),
    stringsAsFactors   = FALSE
  )
}) %>% dplyr::bind_rows()

message("Total genes en intersecciones: ", nrow(lookup_intersecciones))
message("Pares representados: ", n_distinct(lookup_intersecciones$Par))

# --- 2. Unir con hpa_universo ---
datos_intersecciones <- hpa_universo %>%
  dplyr::left_join(lookup_intersecciones, by = "Gene") %>%
  dplyr::mutate(
    etiqueta = ifelse(is.na(etiqueta), "Universe", etiqueta)
  )

# Asignar color a cada etiqueta de intersección
colores_intersecciones <- lookup_intersecciones %>%
  dplyr::select(etiqueta, trait_micro_id) %>%
  dplyr::distinct() %>%
  dplyr::mutate(color = colores_6_traits[trait_micro_id])

colores_completos <- c(
  "Universe" = "#4D4D4D",
  setNames(colores_intersecciones$color, colores_intersecciones$etiqueta)
)

# --- 4. Generar un PDF por trait MicroGWAS ---
traits_micro_unicos <- unique(lookup_intersecciones$trait_micro_id)

for (trait_id in traits_micro_unicos) {
  
  trait_nombre <- traits_MicroGWAS_areas %>%
    dplyr::filter(Rasgo == trait_id) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  
  message("Procesando trait: ", trait_nombre)
  
  # Etiquetas de las intersecciones de este trait
  etiquetas_trait <- lookup_intersecciones %>%
    dplyr::filter(trait_micro_id == trait_id) %>%
    dplyr::pull(etiqueta) %>%
    unique()
  
  # Filtrar datos
  datos_trait <- datos_intersecciones %>%
    dplyr::filter(etiqueta %in% c("Universe", etiquetas_trait))
  
  # Orden: Universe arriba, intersecciones abajo
  orden_trait <- c(rev(etiquetas_trait), "Universe")
  datos_trait$etiqueta <- factor(datos_trait$etiqueta, levels = orden_trait)
  
  # Colores para este trait
  colores_trait <- c(
    "Universe" = "#4D4D4D",
    setNames(rep(colores_6_traits[trait_id], length(etiquetas_trait)),
             etiquetas_trait)
  )
  
  nombre_pdf <- gsub("[^a-zA-Z0-9_]", "_", trait_nombre)
  
  pdf(file.path(carpeta_boxplots_pleio,
                paste0(trait_id, "_", nombre_pdf, "_intersecciones.pdf")),
      width = 10, height = 8)
  
  for (tej in tejidos) {
    
    datos_tej <- datos_trait %>%
      dplyr::filter(Tissue == tej)
    
    p <- ggplot(datos_tej, aes(y = etiqueta, x = log1p(nTPM), fill = etiqueta)) +
      geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
      scale_fill_manual(values = colores_trait) +
      theme_bw() +
      theme(
        axis.text.y     = element_text(size = 8),
        axis.title.y    = element_blank(),
        legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 12, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, size = 8, color = "grey50")
      ) +
      labs(
        title    = tej,
        subtitle = paste0(trait_nombre, " — Intersecciones con Nervous/Psychiatric"),
        x        = "log1p(nTPM)"
      )
    
    print(p)
  }
  
  dev.off()
  message("  ✔️ PDF generado para: ", trait_nombre)
}

message("¡Completado! ", length(traits_micro_unicos), " PDFs generados")


# --- PDF global con todas las intersecciones ---

# Orden fijo: Universe arriba, luego todas las intersecciones
todas_etiquetas <- lookup_intersecciones %>%
  dplyr::select(etiqueta, trait_micro_id) %>%
  dplyr::distinct() %>%
  dplyr::arrange(trait_micro_id) %>%  # agrupar por trait MicroGWAS
  dplyr::pull(etiqueta)

orden_global <- c(rev(todas_etiquetas), "Universe")

datos_intersecciones$etiqueta <- factor(datos_intersecciones$etiqueta, 
                                        levels = orden_global)

# Colores globales
colores_global <- c(
  "Universe" = "#4D4D4D",
  setNames(colores_intersecciones$color, colores_intersecciones$etiqueta)
)

pdf(file.path(carpeta_boxplots_pleio, "TODAS_las_intersecciones_global.pdf"),
    width = 25, height = 30)

for (tej in tejidos) {
  
  datos_tej <- datos_intersecciones %>%
    dplyr::filter(Tissue == tej)
  
  p <- ggplot(datos_tej, aes(y = etiqueta, x = log1p(nTPM), fill = etiqueta)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    scale_fill_manual(values = colores_global) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 16),
      axis.title.y    = element_blank(),
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5, size = 12, color = "grey50")
    ) +
    labs(
      title    = paste0("Expression distribution — ", tej),
      subtitle = "All MicroGWAS vs Nervous/Psychiatric intersections",
      x        = "log1p(nTPM)"
    )
  
  print(p)
}

dev.off()

message("¡PDF global generado con ", length(tejidos), " páginas!")


# =====> HEATMAP ====

# ── Colores NP: distintos de los microbioma 

paleta_np <- c(
  "#142157", "plum", "#E08214", "#2166AC", "#D6604D",
  "#08A045", "#8C510A", "#01665E", "#C51B7D", "#35978F",
  "#a7a7a7", "#40004B", "#ffff00", "#053061", "#67001F"
)

# --- 1. PREPARAR DATOS ---

datos_heatmap_completos <- tabla_maestra_expresion %>%
  dplyr::mutate(
    # Aplicamos limpieza de nombres
    Trait_1_Corto = acortar_nombres_microbioma(Trait_Nombre_1),
    Cluster_1_corto = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_2_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    
    # Nombre final para la columna del heatmap
    Par_Limpio = paste0(Trait_1_Corto, " (", Cluster_1_corto, ")",
                        " vs ",
                        Trait_Nombre_2, " (", Cluster_2_corto, ")"),
    
    # Solo nos quedamos con el log_pval si es significativo y higher
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05, -log10(p.adjust), 0)
  )

# --- 2. CONSTRUIR MATRIZ ---
matriz_pval_completa <- datos_heatmap_completos %>%
  dplyr::select(tissue, Par_Limpio, log_pval) %>%
  tidyr::pivot_wider(names_from = Par_Limpio, 
                     values_from = log_pval,
                     values_fill = 0) %>%
  tibble::column_to_rownames("tissue")

# Limpieza básica: Quitar solo los tejidos (filas) que se han quedado a 0 en absolutamente TODAS las columnas
matriz_pval_completa <- matriz_pval_completa[rowSums(matriz_pval_completa) > 0, , drop = FALSE]

message("Dimensiones de la matriz total: ", nrow(matriz_pval_completa), " tejidos x ", 
        ncol(matriz_pval_completa), " pares pleiotrópicos")

# --- 3. ANOTACIÓN DE COLUMNAS (Microbioma y SNC) 

traits_por_par <- tabla_maestra_expresion %>%
  dplyr::mutate(
    # ¡IMPORTANTE! Hacemos exactamente la misma transformación que arriba para que cuadren
    Trait_1_Corto = acortar_nombres_microbioma(Trait_Nombre_1),
    Cluster_1_corto = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_2_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    Par = paste0(Trait_1_Corto, " (", Cluster_1_corto, ") vs ",
                 Trait_Nombre_2, " (", Cluster_2_corto, ")")
  ) %>%
  dplyr::select(Par, Trait_1_Corto, Trait_Nombre_2) %>%
  dplyr::distinct()

# Dos columnas: una con los traits MicroGWAS cortos y otra con los traits Nervous/Psychiatric
annotation_col <- data.frame(
  Microbioma = traits_por_par$Trait_1_Corto, 
  Enfermedad_SNC = traits_por_par$Trait_Nombre_2,
  row.names = traits_por_par$Par 
)

# Colores Microbioma:
niveles_micro <- sort(unique(traits_por_par$Trait_1_Corto))
colores_micro_c <- colores_micro_por_nombre[niveles_micro]
colores_micro_c <- colores_micro_c[!is.na(colores_micro_c)]   # quitar si algún nombre no matchea

# Colores NP
niveles_snc <- sort(unique(traits_por_par$Trait_Nombre_2))
colores_snc_c <- setNames(paleta_np[seq_along(niveles_snc)], niveles_snc)

annotation_colors_completa <- list(
  Microbioma     = colores_micro_c,
  Enfermedad_SNC = colores_snc_c
)

# --- 4. PALETA Y HEATMAP ---
paleta <- colorRampPalette(c("#EDF8B1", "#7FCDBB", "#2C7FB8", "#081D58"))(100)

# Al tener todas las columnas, necesitamos un PDF bastante grande
pdf("./Output/Gráficos/MicroGWAS/Heatmap_Expr_Tejidos_Interseccion_NervousPsychiatric.pdf", height = 22, width = 30)
print(pheatmap(
  mat               = as.matrix(matriz_pval_completa),
  color             = paleta,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors_completa,
  
  # Agrupamos filas y columnas para que la estructura sea visible
  cluster_rows      = TRUE,
  cluster_cols      = TRUE, 
  
  treeheight_row = 0, 
  treeheight_col = 0,
  
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 8,
  fontsize_col      = 7, # Letra un poco más pequeña para que quepan todos los pares
  angle_col         = 45,
  
  cellheight = 10, 
  cellwidth  = 25,
  
  main              = "Global Tissue Expression of ALL Pleiotropic Intersections (-log10 FDR)",
  border_color      = "white"
))
dev.off()


# =====> DOTPLOT ====
# --- 1. Preparar datos ---
datos_dotplot_expr <- tabla_maestra_expresion %>%
  dplyr::mutate(
    Trait_1_Corto   = acortar_nombres_microbioma(Trait_Nombre_1),
    Cluster_1_corto = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_2_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    
    # Eje X: enfermedad SNC + cluster
    Eje_X = paste0(Trait_Nombre_2, "\n(", Cluster_2_corto, ")"),
    
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05, -log10(p.adjust), NA)
  ) %>%
  dplyr::filter(!is.na(log_pval)) %>%
  dplyr::group_by(Cluster_1, Cluster_2) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

# Ordenar tejidos por p.adjust global (más significativos arriba)
datos_dotplot_expr$tissue <- factor(
  datos_dotplot_expr$tissue,
  levels = unique(datos_dotplot_expr$tissue[order(datos_dotplot_expr$p.adjust, decreasing = TRUE)])
)

# --- 2. Dotplot ---
ggplot(datos_dotplot_expr, aes(x = Eje_X, y = tissue)) +
  geom_point(aes(size = fold_change, color = p.adjust)) +
  scale_color_viridis_c(
    option    = "magma",
    direction = 1,
    trans     = "log10",
    name      = "p.adjust"
  ) +
  scale_size_continuous(name = "Fold Change", range = c(3, 9)) +
  theme_bw() +
  # Facets por rasgo de microbiota (igual que el dotplot de referencia)
  facet_grid(~ Trait_1_Corto, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y        = element_text(size = 14, color = "black"),
    axis.title         = element_blank(),
    plot.title         = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5, size = 10, color = "grey40"),
    
    strip.background   = element_blank(),
    strip.text         = element_text(size = 10, face = "bold", color = "black"),
    
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.spacing      = unit(0.05, "lines"),
    
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  labs(
    title    = "Tissue Expression of Pleiotropic Intersections (higher, FDR < 0.05)",
    subtitle = "Microbiota traits vs Nervous System / Psychiatric Disorders"
  )

# --- 3. Guardar ---
ggsave(
  filename  = "./Output/Gráficos/MicroGWAS/Dotplot_Expr_Tejidos_Interseccion_NervousPsychiatric.pdf",
  width     = 20,
  height    = 15,
  dpi       = 300,
  limitsize = FALSE
)
# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# BLOQUE 3: EXPRESION POR CLUSTER   ----

# Cargar datos de intersección
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds") 
nombres_pares <- names(lista_intersecciones_np)

# ---- Iteramos directamente sobre los nombres de la lista  ----
resultados_ks_log2 <- lapply(nombres_pares, function(nombre_par) {
  
  # Extraemos la información de este par
  datos_par <- lista_intersecciones_np[[nombre_par]]
  genes_cl  <- datos_par$genes_ensembl
  
  message("Procesando par: ", nombre_par)
  
  # --- Test ks por tejido ---
  res_tejidos <- lapply(tejidos, function(tej) {
    test_ks_tejido(
      genes_cluster  = genes_cl,
      genes_universo = genes_interactoma,
      hpa_data       = hpa_universo,
      tejido         = tej
    )
  })
  
  res_df <- dplyr::bind_rows(res_tejidos) 
  
  # Solo ajustamos p-valores si el test devolvió resultados
  if (nrow(res_df) > 0) {
    res_df <- res_df %>%
      dplyr::mutate(
        cluster  = nombre_par,
        trait_1  = datos_par$trait_1,
        trait_2  = datos_par$trait_2,
        jaccard  = datos_par$jaccard,
        p.adjust = p.adjust(pvalue, method = "BH")
      ) %>%
      dplyr::arrange(p.adjust)
    
    message("  — Tejidos significativos (p.adjust < 0.05): ",
            sum(res_df$p.adjust < 0.05))
  }
  
  # --- Expresión gen a gen para este cluster ---
  expr_gen <- hpa_universo %>%
    dplyr::filter(Gene %in% genes_cl) %>%
    dplyr::select(Gene, Gene.name, Tissue, log2_nTPM) %>%
    dplyr::mutate(
      cluster  = nombre_par,
      trait_1  = datos_par$trait_1,
      trait_2  = datos_par$trait_2,
      jaccard  = datos_par$jaccard
    )
  
  # Devolver lista con ambos resultados
  return(list(
    test = res_df,
    expr = expr_gen
  ))
})

names(resultados_ks_log2) <- nombres_pares

# --- Extraer los dos componentes ---
tests_ks_completo <- lapply(resultados_ks_log2, function(x) x$test) %>%
  dplyr::bind_rows()

expresion_completa_ks <- lapply(resultados_ks_log2, function(x) x$expr) %>%
  dplyr::bind_rows()

message("Dimensiones expresión completa: ", nrow(expresion_completa_ks),
        " filas x ", ncol(expresion_completa_ks), " columnas")


# ---- EXPORTAR RESULTADOS   =====

carpeta_resultados_ks <- "./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_NP_MicroGWAS"
dir.create(carpeta_resultados_ks, showWarnings = FALSE, recursive = TRUE)

# Guardar como RDS
saveRDS(expresion_completa_ks,
        file = file.path(carpeta_resultados_ks, 
                         "expresion_genes_por_tejido_all.rds"))
# Resumen global
resumen_ks <- dplyr::bind_rows(tests_ks_completo) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_tejidos_sig  = n(),
    tejidos        = paste(tissue, collapse = ", "),
    mejor_p.adjust = min(p.adjust)
  ) %>%
  dplyr::arrange(mejor_p.adjust)
write.csv2(resumen_ks,
           file = file.path(carpeta_resultados_ks, "resumen_global_KS.csv"),
           row.names = FALSE)

# Guardar objeto para no recalcular
saveRDS(resultados_ks_log2,
        file = file.path(carpeta_resultados_ks, "resultados_KS.rds"))


expresion_completa <- readRDS("./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_NP_MicroGWAS/expresion_genes_por_tejido_all.rds")


# ---- VISUALIZACIÓN DE LOS RESULTADOS ====

# Clusters únicos
clusters_unicos <- unique(expresion_completa_ks$cluster)
message("Clusters a representar: ", length(clusters_unicos))


pdf("./Output/Gráficos/MicroGWAS/Heatmap_Exp_Tejidos_por_Cluster_KS_NP_MicroGWAS.pdf",
    height = 20, width = 20)

for (cl in clusters_unicos) {
  message("Procesando (log2): ", cl)
  datos_cl <- expresion_completa %>% dplyr::filter(cluster == cl)
  heatmap_expresion(datos_cl, tipo = "log2")
}

dev.off()
message("PDF generado con ", length(clusters_unicos), " heatmaps!")

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 4: ===== CALCULAMOS VARIANZA DE LA EXPRESIÓN DE CADA GEN EN TODOS LOS TEJIDOS ====

## Por cada intersección
# 1- Calculamos varianza de cada gen en todos los tejidos   + Z-score
# 2- Ordenamos de mayor a menor varianza
# 3- Heatmap con los genes con mayor varianza

# --- 1. Cargar datos ---
expresion_completa <- readRDS("./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_NP_MicroGWAS/expresion_genes_por_tejido_all.rds")

# --- 2. Calcular varianza por gen y cluster ---
varianza_genes <- expresion_completa %>%
  dplyr::group_by(cluster, Gene, Gene.name) %>%
  dplyr::summarise(
    varianza  = var(log2_nTPM, na.rm = TRUE),
    mediana   = median(log2_nTPM, na.rm = TRUE),
    n_tejidos = sum(!is.na(log2_nTPM)),  # tejidos con expresión no nula
    .groups = "drop"
  ) %>%
  dplyr::arrange(cluster, desc(varianza))

# --- 3. Ver distribución global de varianzas ---
summary(varianza_genes$varianza)

# --- 4. Visualizar distribución por cluster ---
ggplot(varianza_genes, aes(x = varianza)) +
  geom_histogram(bins = 50, fill = "#AEC6CF", color = "white") +
  geom_vline(xintercept = quantile(varianza_genes$varianza, 0.75, na.rm = TRUE),
             color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = quantile(varianza_genes$varianza, 0.90, na.rm = TRUE),
             color = "darkred", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = quantile(varianza_genes$varianza, 0.75, na.rm = TRUE),
           y = Inf, label = "P75", vjust = 2, color = "red", size = 3) +
  annotate("text", x = quantile(varianza_genes$varianza, 0.90, na.rm = TRUE),
           y = Inf, label = "P90", vjust = 2, color = "darkred", size = 3) +
  facet_wrap(~ cluster, scales = "free") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 5),
    axis.text  = element_text(size = 6)
  ) +
  labs(title = "Distribución de varianza por gen y cluster",
       x = "Varianza (log2_nTPM)", y = "Número de genes")


# Empezamos eligiendo varianza >1

genes_alta_varianza <- varianza_genes %>%
  dplyr::filter(varianza > 1.5) %>%
  dplyr::arrange(cluster, desc(varianza))

genes_alta_varianza %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_genes = n()) %>%
  dplyr::arrange(desc(n_genes)) %>%
  print()

expresion_filtrada <- expresion_completa %>%
  dplyr::inner_join(
    genes_alta_varianza %>% dplyr::select(cluster, Gene, varianza),
    by = c("cluster", "Gene")
  )

message("Genes únicos tras filtro: ", n_distinct(expresion_filtrada$Gene))
message("Dimensiones: ", nrow(expresion_filtrada), " filas")

# --- 5. Función para calcular Zscore robusto

zscore_robusto <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  if (mad_val == 0) return(rep(0, length(x)))  # evitar división por cero
  return((x - med) / mad_val)
}

# --- 6. REPRESENTACION GRAFICA - HEATMAP
clusters_unicos <- unique(expresion_filtrada$cluster)

# ---- log2 npTM ----

pdf("./Output/Gráficos/MicroGWAS/Heatmap_Expresion_Varianza_log2_NP_Micro.pdf",
    height = 15, width = 20)

for (cl in clusters_unicos) {
  message("Procesando (log2): ", cl)
  datos_cl <- expresion_filtrada %>% dplyr::filter(cluster == cl)
  heatmap_expresion(datos_cl, tipo = "log2")
}

dev.off()
message("PDF log2 generado!")

# ---- Z-score robusto

pdf("./Output/Gráficos/MicroGWAS/Heatmap_Expresion_varianza_Zscore_NP_Micro_.pdf",
    height = 15, width = 20)

for (cl in clusters_unicos) {
  message("Procesando (Z-score): ", cl)
  datos_cl <- expresion_filtrada %>% dplyr::filter(cluster == cl)
  heatmap_expresion(datos_cl, tipo = "zscore")
}

dev.off()
message("PDF Z-score robusto generado!")




