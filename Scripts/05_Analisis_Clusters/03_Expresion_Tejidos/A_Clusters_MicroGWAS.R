
source("scripts/00_setup.R")  ## Abrimos script con las librerias 
library(clusterProfiler)

# =========================================================================
# ANÁLISIS DE EXPRESIÓN EN TEJIDOS (Gene Ontology)  - en ClustersMicroGWAS
# =========================================================================

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
genes_interactoma <- interactoma[,1]
message("Genes en el interactoma: ", length(genes_interactoma))

# ----------------------------------------------------------------------------------------------------
# BLOQUE 1- ==== DESCARGA DEL DATASET DE EXPRESIÓN  HPA ====

# --- Descomprimir el zip descargado manualmente   -- de HPA -->  rna_tissue_consensus.tsv.zip: ---
unzip("./Data/Diccionarios/rna_tissue_consensus.tsv.zip", 
      exdir = "./Data/Diccionarios")

# --- Leer y explorar ---
hpa_data <- read.delim("./Data/Diccionarios/rna_tissue_consensus.tsv")

# Transformar nTPM a log2 poniendo NA donde nTPM == 0
hpa_data <- hpa_data %>%
  dplyr::mutate(
    log2_nTPM = ifelse(nTPM == 0, NA, log2(nTPM))
  )


message("Tejidos únicos: ", n_distinct(hpa_data$Tissue))
message("Genes únicos: ", n_distinct(hpa_data$Gene))

# --- Filtrar solo genes del interactoma  ---
hpa_universo <- hpa_data %>%
  dplyr::filter(Gene %in% genes_interactoma)

tejidos <- unique(hpa_universo$Tissue)
message("Tejidos a analizar: ", length(tejidos))

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 2: ===== LECTURA DE LOS CLUSTERS =====

## leemos todos los archivos de la carpeta
carpeta_clusters <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
archivos_clusters <- list.files(carpeta_clusters, 
                                pattern = "\\.csv$", 
                                full.names = TRUE)

## Extraer el nombre de cada cluster desde el nombre del archivo
nombres_clusters <- gsub("\\.csv$", "", basename(archivos_clusters))

##  Leer cada CSV y extraer la columna de Ensembl IDs 
cluster_genes <- lapply(archivos_clusters, function(f) {
  df <- read.csv(f)
  # Extraer solo la parte antes del primer ";"
  gsub(";.*$", "", df$gene.ENSG)
})
names(cluster_genes) <- nombres_clusters   # genera una lista llamada `cluster_genes` donde cada elemento es un vector de Ensembl IDs
message("Ejemplo gen en cluster tras corrección: ", cluster_genes[[1]][1])

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 3: ==== GSEA-like CON TEST KS POR CLUSTER Y TEJIDO ====

# --- Función que aplica el test KS para un cluster en un tejido ---
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

# --- Aplicar a todos los clusters y tejidos ---
resultados_ks <- lapply(seq_along(cluster_genes), function(i) {
  genes_cl <- cluster_genes[[i]]
  cl <- nombres_clusters[i]
  message("Procesando cluster ", i, ": ", cl)
  
  # Aplicar test KS para cada tejido
  res_tejidos <- lapply(tejidos, function(tej) {
    test_ks_tejido(
      genes_cluster  = genes_cl,
      genes_universo = genes_interactoma,
      hpa_data       = hpa_universo,
      tejido         = tej
    )
  })
  
  # Combinar resultados y corregir p-valores
  res_df <- bind_rows(res_tejidos) %>%
    dplyr::mutate(
      cluster    = cl,
      p.adjust   = p.adjust(pvalue, method = "BH")
    ) %>%
    dplyr::arrange(p.adjust)
  
  message("  — Tejidos significativos (p.adjust < 0.05, higher): ",
          sum(res_df$p.adjust < 0.05 & res_df$direccion == "higher"))
  
  return(res_df)
})
names(resultados_ks) <- nombres_clusters

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 4: ==== EXPORTAR LOS RESULTADOS ====

carpeta_resultados <- "./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_ClustersMicroGWAS"
dir.create(carpeta_resultados, showWarnings = FALSE, recursive = TRUE)

# Un CSV por cluster con tejidos significativos
invisible(lapply(seq_along(resultados_ks), function(i) {
  res <- resultados_ks[[i]]
  cl <- nombres_clusters[i]
  
  res_sig <- res %>%
    dplyr::filter(p.adjust < 0.05, direccion == "higher") %>%
    dplyr::arrange(p.adjust)
  
  if (nrow(res_sig) > 0) {
    write.csv2(res_sig,
               file = file.path(carpeta_resultados,
                                paste0(cl, "_KS_HPA.csv")),
               row.names = FALSE)
  }
}))

# --- Resumen global ---
resumen_ks <- bind_rows(resultados_ks) %>%
  dplyr::filter(p.adjust < 0.05, direccion == "higher") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_tejidos_sig  = n(),
    tejidos        = paste(tissue, collapse = ", "),
    mejor_p.adjust = min(p.adjust)
  ) %>%
  dplyr::arrange(mejor_p.adjust)

print(resumen_ks)

write.csv2(resumen_ks,
           file = file.path(carpeta_resultados, "resumen_global_KS.csv"),
           row.names = FALSE)

# --- Resumen incluyendo todos los clusters ---

todos_clusters_df <- data.frame(cluster = nombres_clusters)

resumen_completo <- todos_clusters_df %>%
  dplyr::left_join(resumen_ks, by = "cluster") %>%
  dplyr::mutate(
    n_tejidos_sig  = ifelse(is.na(n_tejidos_sig), 0, n_tejidos_sig),
    tejidos        = ifelse(is.na(tejidos), "No significant tissue enrichment", tejidos),
    mejor_p.adjust = ifelse(is.na(mejor_p.adjust), NA, mejor_p.adjust)
  ) %>%
  dplyr::arrange(desc(n_tejidos_sig))

print(resumen_completo)

write.csv2(resumen_completo,
           file = file.path(carpeta_resultados, "resumen_completo_todos_clusters.csv"),
           row.names = FALSE)


# --- Ver resultados exploratorios para clusters sin significancia ---
clusters_sin_sig <- c("EFO_0007883_Cluster_1.3",
                      "EFO_0007883_Cluster_1.4.3",
                      "EFO_0007883_Cluster_3.1.4",
                      "EFO_0801228_Cluster_1.1.7",
                      "EFO_0801228_Cluster_1.4",
                      "EFO_0801228_Cluster_2",
                      "EFO_0801229_Cluster_9",
                      "EFO_0011013_Cluster_1.2.2")

lapply(clusters_sin_sig, function(cl) {
  res <- resultados_ks[[cl]]
  
  message("=== ", cl, " ===")
  res %>%
    dplyr::filter(direccion == "higher") %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(tissue, fold_change, pvalue, p.adjust) %>%
    head(5) %>%
    print()
})

# Ver distribución de direcciones para estos clusters
lapply(clusters_sin_sig, function(cl) {
  res <- resultados_ks[[cl]]
  message("=== ", cl, " ===")
  message("  higher: ", sum(res$direccion == "higher"))
  message("  lower:  ", sum(res$direccion == "lower"))
  message("  mejor pvalue (cualquier dirección): ", 
          round(min(res$pvalue), 4))
  message("  fold_change rango: ", 
          round(min(res$fold_change), 3), " - ", 
          round(max(res$fold_change), 3))
})

# Clasificar cada cluster en una categoría interpretativa
resumen_interpretado <- data.frame(cluster = nombres_clusters) %>%
  dplyr::left_join(resumen_ks, by = "cluster") %>%
  dplyr::mutate(
    categoria = dplyr::case_when(
      !is.na(n_tejidos_sig) & n_tejidos_sig == 51 ~ "Ubiquitous expression",
      !is.na(n_tejidos_sig) & n_tejidos_sig > 0   ~ "Tissue-specific enrichment",
      cluster %in% clusters_sin_sig               ~ "Lower expression than background",
      TRUE                                         ~ "No significant enrichment"
    ),
    n_tejidos_sig  = ifelse(is.na(n_tejidos_sig), 0, n_tejidos_sig),
    tejidos        = ifelse(is.na(tejidos), "-", tejidos),
    mejor_p.adjust = ifelse(is.na(mejor_p.adjust), NA, mejor_p.adjust)
  ) %>%
  dplyr::arrange(categoria, desc(n_tejidos_sig))

print(resumen_interpretado)

write.csv2(resumen_interpretado,
           file = file.path(carpeta_resultados, 
                            "resumen_interpretado_todos_clusters.csv"),
           row.names = FALSE)

# ----------------------------------------------------------------------------------------------------

 #### ==== -------------------------------------------------------------------------------------------

# Guardar resultados_ks para no tener que recalcular
saveRDS(resultados_ks, 
        file = "./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_ClustersMicroGWAS/resultados_ks.rds")

resultados_ks <- readRDS("./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_ClustersMicroGWAS/resultados_ks.rds")

#### ====  -------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# BLOQUE 5: ===== VISUALIZAR DISTRIBUCIÓN EXPRESIÓN POR TEJIDO DE CADA CLUSTER  ----



# Si falta cluster_genes
carpeta_clusters <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)
nombres_clusters <- gsub("\\.csv$", "", basename(archivos_clusters))
cluster_genes <- lapply(archivos_clusters, function(f) {
  df <- read.csv(f)
  gsub(";.*$", "", df$gene.ENSG)
})
names(cluster_genes) <- nombres_clusters


carpeta_boxplots <- "./Output/Piloto_Microbiota/Boxplots_Expresion_Tejidos_MicroGWAS"
dir.create(carpeta_boxplots, showWarnings = FALSE)

dicc_cluster_trait <- data.frame(
  cluster_full  = nombres_clusters,
  cluster_corto = gsub("^EFO_[0-9]+_", "", nombres_clusters),
  trait_id      = gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", nombres_clusters)
) %>%
  dplyr::left_join(
    traits_MicroGWAS_areas %>%
      dplyr::select(Rasgo, name) %>%
      dplyr::distinct(),
    by = c("trait_id" = "Rasgo")
  ) %>%
  dplyr::mutate(
    etiqueta = paste0(name, "\n(", cluster_corto, ")")
  )
lookup_gen_cluster <- lapply(nombres_clusters, function(cl) {
  data.frame(
    Gene     = cluster_genes[[cl]],
    etiqueta = dicc_cluster_trait$etiqueta[dicc_cluster_trait$cluster_full == cl],
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()

datos_todos <- hpa_universo %>%
  dplyr::left_join(lookup_gen_cluster, by = "Gene") %>%
  dplyr::mutate(
    etiqueta = ifelse(is.na(etiqueta), "Universe", etiqueta)
  )

# --- 1. Obtener traits únicos ---
traits_unicos <- unique(dicc_cluster_trait$trait_id)
message("Traits únicos: ", length(traits_unicos))

# Asignarle un color a cada cluster

colores_6_traits <- c(
  "EFO_0007753" = "#FFFF99",  # amarillo
  "EFO_0007874" = "#FFDAC1",  # naranja
  "EFO_0007883" = "#AEC6CF",  # azul
  "EFO_0011013" = "#B5EAD7",  # verde
  "EFO_0801228" = "#C9B1FF",  # morado
  "EFO_0801229" = "#FFD1DC"   # rosa
)

# Asignar color de trait a cada cluster
dicc_cluster_trait <- dicc_cluster_trait %>%
  dplyr::mutate(color = colores_6_traits[trait_id])

colores_clusters <- setNames(dicc_cluster_trait$color, dicc_cluster_trait$etiqueta)
colores_todos <- c("Universe" = "#4D4D4D", colores_clusters)


# --- 2. Generar un PDF por trait ---
for (trait_id in traits_unicos) {
  
  # Nombre legible del trait
  trait_nombre <- traits_MicroGWAS_areas %>%
    dplyr::filter(Rasgo == trait_id) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  
  message("Procesando trait: ", trait_nombre)
  
  # Clusters de este trait
  clusters_trait <- dicc_cluster_trait %>%
    dplyr::filter(trait_id == !!trait_id) %>%
    dplyr::pull(etiqueta)
  
  # Filtrar datos para este trait + Universe
  datos_trait <- datos_todos %>%
    dplyr::filter(etiqueta %in% c("Universe", clusters_trait))
  
  # Orden fijo: Universe primero, luego clusters del trait
  orden_trait <- c(rev(clusters_trait), "Universe")
  datos_trait$etiqueta <- factor(datos_trait$etiqueta, levels = orden_trait)
  
  # Color del trait + Universe oscuro
  color_trait <- colores_6_traits[trait_id]
  colores_trait <- c("Universe" = "#4D4D4D",
                     setNames(rep(color_trait, length(clusters_trait)), 
                              clusters_trait))
  
  # Nombre del archivo PDF (sin caracteres especiales)
  nombre_pdf <- gsub("[^a-zA-Z0-9_]", "_", trait_nombre)
  
  pdf(file.path(carpeta_boxplots, paste0(trait_id, "_", nombre_pdf, ".pdf")),
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
        subtitle = trait_nombre,
        x        = "log1p(nTPM)"
      )
    
    print(p)
  }
  
  dev.off()
  message("  ✔️ PDF generado para: ", trait_nombre)
}

message("¡Completado! ", length(traits_unicos), " PDFs generados")





# --- PDF con todos los clusters de todos los traits ---

# Orden fijo: Universe arriba, luego todos los clusters agrupados por trait
orden_todos <- c(rev(dicc_cluster_trait$etiqueta), "Universe")

datos_todos$etiqueta <- factor(datos_todos$etiqueta, levels = orden_todos)

# Colores: todos los clusters con su color de trait + Universe oscuro
colores_completos <- c("Universe" = "#4D4D4D", colores_clusters)

pdf(file.path(carpeta_boxplots, "TODOS_los_clusters_todos_traits.pdf"),
    width = 20, height = 30)

for (tej in tejidos) {
  
  datos_tej <- datos_todos %>%
    dplyr::filter(Tissue == tej)
  
  p <- ggplot(datos_tej, aes(y = etiqueta, x = log1p(nTPM), fill = etiqueta)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    scale_fill_manual(values = colores_completos) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 7),
      axis.title.y    = element_blank(),
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    labs(
      title = paste0("Expression distribution — ", tej, " — All clusters"),
      x     = "log1p(nTPM)"
    )
  
  print(p)
}

dev.off()

message("¡PDF global generado con ", length(tejidos), " páginas!")
# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 6: ==== VISUALIZACIÓN RESULTADOS - HEATMAP 
acortar_nombres_microbioma <- function(nombres) {
  n <- gsub(" microbiome measurement", " Micr.", nombres, ignore.case = TRUE)
  return(n)}

colores_6_traits <- c(
  "EFO_0007753" = "#FFFF99",
  "EFO_0007874" = "#FFDAC1",
  "EFO_0007883" = "#AEC6CF",
  "EFO_0011013" = "#B5EAD7",
  "EFO_0801228" = "#C9B1FF",
  "EFO_0801229" = "#FFD1DC")

efo_a_nombre <- traits_MicroGWAS_areas %>%
  dplyr::filter(Rasgo %in% names(colores_6_traits)) %>%
  dplyr::select(Rasgo, name) %>%
  dplyr::distinct() %>%
  tibble::deframe()

colores_micro_por_nombre <- setNames(
  colores_6_traits[names(efo_a_nombre)],
  acortar_nombres_microbioma(efo_a_nombre))

# ==== HEATMAP ====

# 1. --- CONSTRUIR MATRIZ DE -log10(p.adjust)  ---

# Extraer todos los resultados en un único dataframe
todos_resultados <- bind_rows(resultados_ks) %>%
  dplyr::mutate(cluster = dplyr::recode(cluster, !!!dic_clusters_MicroGWAS))

# 1. --- MATRIZ ---
matriz_pval <- todos_resultados %>%
  dplyr::mutate(
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05,
                      -log10(p.adjust), 0)
  ) %>%
  dplyr::select(tissue, cluster, log_pval) %>%
  tidyr::pivot_wider(names_from  = cluster,
                     values_from = log_pval,
                     values_fill = 0) %>%
  tibble::column_to_rownames("tissue")

message("Dimensiones de la matriz: ", nrow(matriz_pval), " tejidos x ",
        ncol(matriz_pval), " clusters")

# 2. --- ANOTACIÓN DE COLUMNAS ---
efo_por_cluster      <- gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", colnames(matriz_pval))
efo_a_nombre_completo <- setNames(traits_MicroGWAS_areas$name,
                                  traits_MicroGWAS_areas$Rasgo)
nombres_reales       <- efo_a_nombre_completo[efo_por_cluster]
nombres_reales_cortos <- acortar_nombres_microbioma(nombres_reales)

annotation_col <- data.frame(
  Trait     = nombres_reales_cortos,
  row.names = colnames(matriz_pval)
)

mis_colores_heatmap <- list(Trait = colores_micro_por_nombre)

colnames_cortos <- gsub("^EFO_[0-9]+_", "", colnames(matriz_pval))

# Distancias (reemplazando NAs por 0 para el cálculo)
matriz_calculo <- as.matrix(matriz_pval)
matriz_calculo[is.na(matriz_calculo)] <- 0

# Distancias desde la matriz completa con ceros (sin filtrar filas/cols)
dist_tejidos  <- dist(matriz_calculo, method = "euclidean")
dist_clusters <- dist(t(matriz_calculo), method = "euclidean")

# Matriz para plotear (ceros → NA para que salgan grises)
matriz_plot <- as.matrix(matriz_pval)
matriz_plot[matriz_plot == 0] <- NA

#ponemos breaks
breaks <- seq(0, 30, length.out = 101)
paleta <- colorRampPalette(
  c("#EEF9C4", "#7FCDBB", "#2C7FB8", "#1A4E88"))(100)

pdf("./Output/Gráficos/MicroGWAS/Heatmap_Enrich_Tejidos_ClustersMicroGWAS.pdf",
    height = 20, width = 20)
print(pheatmap(
  mat                      = matriz_plot,
  color                    = paleta,
  breaks = breaks,
  annotation_col           = annotation_col,
  annotation_colors        = mis_colores_heatmap,
  labels_col               = colnames_cortos,
  cluster_rows             = TRUE,
  cluster_cols             = TRUE,
  clustering_distance_rows = dist_tejidos,
  clustering_distance_cols = dist_clusters,
  clustering_method        = "average",
  treeheight_row           = 0,
  treeheight_col           = 0,
  show_rownames            = TRUE,
  show_colnames            = TRUE,
  fontsize_row             = 9,
  fontsize_col             = 8,
  angle_col                = 45,
  cellheight               = 10,
  cellwidth                = 15,
  main                     = "Tissue enrichment by cluster (-log10 p.adjust)",
  na_col                   = "whitesmoke",
  border_color             = "#DCDCDC"
))
dev.off()


# ==== DOTPLOT ====

# --- 1. Preparar datos ---
datos_dotplot_ks <- todos_resultados %>%
  dplyr::mutate(
    # Nombres cortos de cluster (igual que el heatmap)
    Cluster_corto = gsub("^EFO_[0-9]+_", "", cluster),
    
    # EFO y nombre del trait (igual que el heatmap)
    efo_id      = gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", cluster),
    Trait_Nombre = efo_a_nombre[efo_id],
    Trait_Corto   = acortar_nombres_microbioma(Trait_Nombre), 
    Cluster_label = paste0(Trait_Corto, "\n(", Cluster_corto, ")"),
    
    # Eje X: cluster corto
    Eje_X = Cluster_corto,
    
    # Solo significativos y higher
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05, -log10(p.adjust), NA)
  ) %>%
  dplyr::filter(!is.na(log_pval)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

# Ordenar tejidos por p.adjust global (más significativos arriba)
datos_dotplot_ks$tissue <- factor(
  datos_dotplot_ks$tissue,
  levels = unique(datos_dotplot_ks$tissue[order(datos_dotplot_ks$p.adjust, decreasing = TRUE)])
)

# --- 2. Dotplot ---
ggplot(datos_dotplot_ks, aes(x = Eje_X, y = tissue)) +
  geom_point(aes(size = fold_change, color = p.adjust)) +
  scale_color_viridis_c(
    option    = "magma",
    direction = 1,
    trans     = "log10",
    name      = "p.adjust"
  ) +
  scale_size_continuous(name = "Fold Change", range = c(3, 8)) +
  theme_bw() +
  # Facets por nombre del trait (igual que el dotplot de referencia)
  facet_grid(~ Trait_Corto, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y        = element_text(size = 12, color = "black"),
    axis.title         = element_blank(),
    plot.title         = element_text(hjust = 0.5, size = 14, face = "bold"),
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
    title    = "Tissue enrichment by MicroGWAS Cluster (higher, FDR < 0.05)",
    subtitle = "Kolmogorov-Smirnov test vs interactome background"
  )

# --- 3. Guardar ---
ggsave(
  filename  = "./Output/Gráficos/MicroGWAS/Dotplot_Enrich_Tejidos_ClustersMicroGWAS.pdf",
  width     = 18,
  height    = 14,
  dpi       = 300,
  limitsize = FALSE
)

# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# BLOQUE 7: ==== EXPRESIÓN TEJIDOS EN LA INTERSECCIÓN ====


# --- 1. Cargar matriz de Jaccard y filtrar pares ---
resultados_jaccard <- read.csv2(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterGWAS/Matriz_Jaccard_Completa.csv"
)

pares_pleiotropia <- resultados_jaccard %>%
  dplyr::filter(Indice_Jaccard >= 0.5)

message("Pares con Jaccard >= 0.5: ", nrow(pares_pleiotropia))

# --- 2. Recargar clusters con Ensembl IDs ---
carpeta_clusters <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)

lista_genes_ensembl <- list()

for (f in archivos_clusters) {
  nombre <- gsub("\\.csv$", "", basename(f))
  df     <- read.csv2(f)
  lista_genes_ensembl[[nombre]] <- gsub(";.*$", "", df$gene.ENSG)
}

# --- 3. Aplicar test KS a la intersección de cada par ---
carpeta_expresion_pleio <- "./Output/Piloto_Microbiota/Expresion_Pleiotropia_MicroGWAS"
dir.create(carpeta_expresion_pleio, showWarnings = FALSE, recursive = TRUE)

lista_resultados_expr_pleio <- list()

for (i in 1:nrow(pares_pleiotropia)) {
  
  c1      <- pares_pleiotropia$Cluster_1[i]
  c2      <- pares_pleiotropia$Cluster_2[i]
  jaccard <- pares_pleiotropia$Indice_Jaccard[i]
  
  message("=== Par ", i, ": ", c1, " vs ", c2, " (Jaccard = ", jaccard, ") ===")
  
  # Genes de cada cluster
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
  
  # Aplicar test KS para cada tejido usando los genes de la intersección
  res_tejidos <- lapply(tejidos, function(tej) {
    test_ks_tejido(
      genes_cluster  = genes_interseccion,
      genes_universo = genes_interactoma,
      hpa_data       = hpa_universo,
      tejido         = tej
    )
  })
  
  # Combinar y ajustar p-valores
  res_df <- dplyr::bind_rows(res_tejidos)
  
  if (nrow(res_df) == 0) {
    message("  ⚠️ No hubo suficientes datos de expresión para este par")
    next
  }
  
  res_df <- res_df %>%
    dplyr::mutate(
      Cluster_1      = c1,
      Cluster_2      = c2,
      Par            = paste0(c1, "_vs_", c2),
      Trait_Nombre_1 = pares_pleiotropia$Trait_Nombre_1[i],
      Trait_Nombre_2 = pares_pleiotropia$Trait_Nombre_2[i],
      Indice_Jaccard = jaccard,
      N_genes_interseccion = length(genes_interseccion),
      p.adjust       = p.adjust(pvalue, method = "BH")
    ) %>%
    dplyr::arrange(p.adjust)
  
  message("  ✔️ Tejidos significativos (p.adjust < 0.05, higher): ",
          sum(res_df$p.adjust < 0.05 & res_df$direccion == "higher"))
  
  # Guardar CSV individual
  nombre_par <- paste0(c1, "_vs_", c2)
  write.csv2(res_df,
             file = file.path(carpeta_expresion_pleio,
                              paste0("Expresion_KS_", nombre_par, ".csv")),
             row.names = FALSE)
  
  lista_resultados_expr_pleio[[nombre_par]] <- res_df
}

# --- 4. Tabla maestra ---
if (length(lista_resultados_expr_pleio) > 0) {
  
  tabla_maestra_expr_pleio <- dplyr::bind_rows(lista_resultados_expr_pleio)
  
  write.csv2(tabla_maestra_expr_pleio,
             file = file.path(carpeta_expresion_pleio, 
                              "Tabla_Maestra_Expresion_Pleiotropia_MicroGWAS.csv"),
             row.names = FALSE)
  
  message("\n¡Análisis completado! Tabla maestra guardada con ",
          nrow(tabla_maestra_expr_pleio), " filas.")
}

# --- 5. Resumen interpretado ---
pares_evaluados <- unique(tabla_maestra_expr_pleio$Par)

resumen_expr_pleio <- tabla_maestra_expr_pleio %>%
  dplyr::filter(p.adjust < 0.05, direccion == "higher") %>%
  dplyr::group_by(Par, Trait_Nombre_1, Trait_Nombre_2, Indice_Jaccard) %>%
  dplyr::summarise(
    n_tejidos_sig  = n(),
    tejidos        = paste(tissue, collapse = ", "),
    mejor_p.adjust = min(p.adjust),
    .groups        = "drop"
  ) %>%
  dplyr::arrange(mejor_p.adjust)

# Incluir pares sin significancia
pares_todos_df <- pares_pleiotropia %>%
  dplyr::mutate(Par = paste0(Cluster_1, "_vs_", Cluster_2)) %>%
  dplyr::select(Par, Trait_Nombre_1, Trait_Nombre_2, Indice_Jaccard)

resumen_interpretado_expr <- pares_todos_df %>%
  dplyr::left_join(resumen_expr_pleio, 
                   by = c("Par", "Trait_Nombre_1", "Trait_Nombre_2", "Indice_Jaccard")) %>%
  dplyr::mutate(
    categoria = dplyr::case_when(
      !Par %in% pares_evaluados                   ~ "Skipped (< 5 shared genes or no data)",
      !is.na(n_tejidos_sig) & n_tejidos_sig >= 50 ~ "Ubiquitous expression",
      !is.na(n_tejidos_sig) & n_tejidos_sig > 0   ~ "Tissue-specific enrichment",
      TRUE                                         ~ "No significant enrichment / Lower"
    ),
    n_tejidos_sig  = ifelse(is.na(n_tejidos_sig), 0, n_tejidos_sig),
    tejidos        = ifelse(is.na(tejidos), "-", tejidos)
  ) %>%
  dplyr::arrange(categoria, desc(n_tejidos_sig))

print(resumen_interpretado_expr)

write.csv2(resumen_interpretado_expr,
           file = file.path(carpeta_expresion_pleio,
                            "Resumen_Interpretado_Expresion_Pleiotropia_MicroGWAS.csv"),
           row.names = FALSE)


# ===== HEATMAP - Expresión tejidos intersecciones pleiotrópicas MicroGWAS =====

# --- 1. Preparar datos ---

datos_heatmap <- tabla_maestra_expr_pleio %>%
  dplyr::mutate(
    Trait_1_Corto   = acortar_nombres_microbioma(Trait_Nombre_1),
    Trait_2_Corto   = acortar_nombres_microbioma(Trait_Nombre_2),
    Cluster_1_corto = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_2_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    
    Par_Limpio = paste0(Trait_1_Corto, " (", Cluster_1_corto, ")",
                        " vs ",
                        Trait_2_Corto, " (", Cluster_2_corto, ")"),
    
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05, -log10(p.adjust), 0)
  )

# --- 2. Construir matriz: filas = tejidos, columnas = pares ---
matriz_heatmap <- datos_heatmap %>%
  dplyr::select(tissue, Par_Limpio, log_pval) %>%
  tidyr::pivot_wider(names_from  = Par_Limpio,
                     values_from = log_pval,
                     values_fill = 0) %>%
  tibble::column_to_rownames("tissue")

# Quitar tejidos con todo a 0
matriz_heatmap <- matriz_heatmap[rowSums(matriz_heatmap) > 0, , drop = FALSE]

message("Dimensiones matriz: ", nrow(matriz_heatmap), " tejidos x ", 
        ncol(matriz_heatmap), " pares")

# --- 3. Anotación de columnas ---
anotaciones_df <- datos_heatmap %>%
  dplyr::select(Par_Limpio, Trait_1_Corto, Trait_2_Corto) %>%
  dplyr::distinct()

annotation_col <- data.frame(
  Microbioma_1 = anotaciones_df$Trait_1_Corto,
  Microbioma_2 = anotaciones_df$Trait_2_Corto,
  row.names    = anotaciones_df$Par_Limpio
)

# Paletas anotación
niveles_micro1 <- unique(anotaciones_df$Trait_1_Corto)
niveles_micro2 <- unique(anotaciones_df$Trait_2_Corto)

paleta_1 <- c("#264653","#2A9D8F","#E9C46A","#F4A261","#E76F51","#6D597A",
              "#3D405B","#81B29A","#F2CC8F","#E07A5F")
paleta_2 <- c("#AEC6CF","#FFD1DC","#B5EAD7","#C9B1FF","#FFDAC1","#FF8B94",
              "#B5C0D0","#D8A7B1","#EFD3D7","#DCEDC1")

colores_micro1 <- setNames(rep(paleta_1, length.out = length(niveles_micro1)), niveles_micro1)
colores_micro2 <- setNames(rep(paleta_2, length.out = length(niveles_micro2)), niveles_micro2)

annotation_colors <- list(
  Microbioma_1 = colores_micro1,
  Microbioma_2 = colores_micro2
)

# --- 4. Paleta y heatmap ---
paleta <- colorRampPalette(c("#EDF8B1","#7FCDBB","#2C7FB8","#081D58"))(100)

pdf("./Output/Gráficos/MicroGWAS/Heatmap_Expr_Tejidos_Interseccion_MicroGWAS.pdf", 
    height = 20, width = 22)
print(pheatmap(
  mat               = as.matrix(matriz_heatmap),
  color             = paleta,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 8,
  fontsize_col      = 7,
  angle_col         = 45,
  cellheight        = 10,
  cellwidth         = 20,
  main              = "Tissue Expression of Pleiotropic Intersections MicroGWAS (-log10 FDR)",
  border_color      = "white"
))
dev.off()

# ==== DOTPLOT - Expresión tejidos intersecciones pleiotrópicas MicroGWAS ====

# --- 1. Preparar datos ---
datos_dotplot <- tabla_maestra_expr_pleio %>%
  dplyr::mutate(
    Trait_1_Corto   = acortar_nombres_microbioma(Trait_Nombre_1),
    Trait_2_Corto   = acortar_nombres_microbioma(Trait_Nombre_2),
    Cluster_1_corto = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_2_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    
    # Eje X: cluster 2 dentro de cada facet de microbioma 1
    Eje_X = paste0(Trait_2_Corto, "\n(", Cluster_2_corto, ")"),
    
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05, -log10(p.adjust), NA)
  ) %>%
  dplyr::filter(!is.na(log_pval)) %>%
  dplyr::group_by(Cluster_1, Cluster_2) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

# Ordenar tejidos por significancia global
datos_dotplot$tissue <- factor(
  datos_dotplot$tissue,
  levels = unique(datos_dotplot$tissue[order(datos_dotplot$p.adjust, decreasing = TRUE)])
)

# --- 2. Dotplot ---
ggplot(datos_dotplot, aes(x = Eje_X, y = tissue)) +
  geom_point(aes(size = fold_change, color = p.adjust)) +
  scale_color_viridis_c(
    option    = "magma",
    direction = 1,
    trans     = "log10",
    name      = "p.adjust"
  ) +
  scale_size_continuous(name = "Fold Change", range = c(4, 12)) +
  theme_bw() +
  facet_grid(~ Trait_1_Corto, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y        = element_text(size = 9, color = "black"),
    axis.title         = element_blank(),
    plot.title         = element_text(hjust = 0.5, size = 14, face = "bold"),
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
    title    = "Tissue Expression of Pleiotropic Intersections MicroGWAS (higher, FDR < 0.05)",
    #subtitle = "Kolmogorov-Smirnov test vs interactome background"
  )

# --- 3. Guardar ---
ggsave(
  filename  = "./Output/Gráficos/MicroGWAS/Dotplot_Expr_Tejidos_Interseccion_MicroGWAS.pdf",
  width     = 12,
  height    = 10,
  dpi       = 300,
  limitsize = FALSE
)




# ----------------------------------------------------------------------------------------------------






