
source("scripts/00_setup.R")  ## Abrimos script con las librerias 
library(clusterProfiler)

# =========================================================================
# ANÁLISIS DE EXPRESIÓN EN TEJIDOS (Gene Ontology)  

# INTERSECCION MICORGWAS - NERVOUS / PSYCHIATRIC
# =========================================================================

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
genes_interactoma <- interactoma[,1]
message("Genes en el interactoma: ", length(genes_interactoma))

lista_intersecciones_im <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds") 

# ----- Datos de HPA
hpa_data <- read.delim("./Data/Diccionarios/rna_tissue_consensus.tsv")
# Transformar nTPM a log2 poniendo NA donde nTPM == 0
hpa_data <- hpa_data %>%
  dplyr::mutate(
    log2_nTPM = ifelse(nTPM == 0, NA, log2(nTPM))
  )

hpa_universo <- hpa_data %>%
  dplyr::filter(Gene %in% genes_interactoma)

tejidos <- unique(hpa_universo$Tissue)
message("Tejidos a analizar: ", length(tejidos))

# FUNCION

acortar_nombres_microbioma <- function(nombres) {
  n <- gsub(" microbiome measurement", " Micr.", nombres, ignore.case = TRUE)
  return(n)
}

# ── Colores MicroGWAS ya establecidos ────────────────────────────────────────

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

# ========================================================================
# BLOQUE 1: ANÁLISIS EXPRESIÓN HPA ====

# ---- Función que aplica el test KS para un cluster en un tejido 
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


carpeta_expresion_im <- "./Output/Piloto_Microbiota/Expresion_Tejidos_IM_MicroGWAS"
dir.create(carpeta_expresion_im, showWarnings = FALSE)

lista_resultados_pleiotropia_tejidos <- list()

for (nombre_par in names(lista_intersecciones_im)) {
  
  par <- lista_intersecciones_im[[nombre_par]]
  
  message("=== Par: ", nombre_par, " (Jaccard = ", round(par$jaccard, 3), ") ===")
  
  genes_interseccion <- par$genes_tabla$gene.ENSG 
  
  # Limpiamos los Ensembl IDs por si vienen múltiples separados por ";" o con versión (".1")
  genes_interseccion <- gsub(";.*$", "", genes_interseccion)
  genes_interseccion <- unique(genes_interseccion)
  
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
        Cluster_1 = par$c1,
        Cluster_2 = par$c2,
        Par = nombre_par,
        Trait_Nombre_1 = par$trait_1, 
        Trait_Nombre_2 = par$trait_2,
        p.adjust = p.adjust(pvalue, method = "BH")
      ) %>%
      dplyr::arrange(p.adjust)
    
    message("  ✔️ Tejidos significativos (p.adj < 0.05, higher): ",
            sum(res_df$p.adjust < 0.05 & res_df$direccion == "higher"))
    
    # Guardar CSV individual
    write.csv2(res_df,
               file = file.path(carpeta_expresion_im,
                                paste0("Expresion_Pleiotropia_", nombre_par, ".csv")),
               row.names = FALSE)
    
    # Añadir a la lista para luego hacer una tabla global
    lista_resultados_pleiotropia_tejidos[[nombre_par]] <- res_df
    
  } else {
    message("  ⚠️ No hubo suficientes datos de expresión para evaluar los tejidos en este par.")
  }
}


# ---- Tabla Maestra

if(length(lista_resultados_pleiotropia_tejidos) > 0) {
  tabla_maestra_expresion <- dplyr::bind_rows(lista_resultados_pleiotropia_tejidos)
  
  write.csv2(tabla_maestra_expresion,
             file = file.path(carpeta_expresion_im, "Tabla_Maestra_Expresion_Pleiotropia.csv"),
             row.names = FALSE)
  message("\n¡Análisis completado! Tabla maestra guardada.")
}

# ========================================================================

# ========================================================================
# BLOQUE 2: ===== EXPORTAR LOS RESULTADOS SIMPLIFICADOS ====

# --- 1. Resumen global (Solo los pares que tienen algún tejido significativo) ---
resumen_ks <- tabla_maestra_expresion %>%
  dplyr::filter(p.adjust < 0.05, direccion == "higher") %>%
  dplyr::group_by(Par, Trait_Nombre_1, Trait_Nombre_2) %>% 
  dplyr::summarise(
    n_tejidos_sig  = n(),
    tejidos        = paste(tissue, collapse = ", "),
    mejor_p.adjust = min(p.adjust),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mejor_p.adjust)

write.csv2(resumen_ks,
           file = file.path(carpeta_expresion_im, "resumen_global_KS.csv"),
           row.names = FALSE)

# --- 2. Clasificar cada Par en una categoría interpretativa ---

pares_todos_df <- data.frame(
  Par = names(lista_intersecciones_im),
  stringsAsFactors = FALSE
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Trait_Nombre_1 = lista_intersecciones_im[[Par]]$trait_1,
    Trait_Nombre_2 = lista_intersecciones_im[[Par]]$trait_2
  ) %>%
  dplyr::ungroup()

# Identificamos cuáles sí llegaron a evaluarse en la tabla maestra (los que tenían >= 5 genes)
pares_evaluados <- unique(tabla_maestra_expresion$Par)

# Cruzamos la información y asignamos las categorías (¡CORREGIDO EL LEFT_JOIN!)
resumen_interpretado <- pares_todos_df %>%
  dplyr::left_join(resumen_ks, by = c("Par", "Trait_Nombre_1", "Trait_Nombre_2")) %>%
  dplyr::mutate(
    categoria = dplyr::case_when(
      # Si el par no está en los evaluados, es que se saltó por falta de genes compartidos
      !Par %in% pares_evaluados ~ "Skipped (< 5 shared genes or no data)",
      # Si es significativo en muchísimos tejidos, es expresión basal/ubicua
      !is.na(n_tejidos_sig) & n_tejidos_sig >= 50 ~ "Ubiquitous expression",
      # Si es significativo en algunos tejidos concretos (ej. solo cerebro)
      !is.na(n_tejidos_sig) & n_tejidos_sig > 0   ~ "Tissue-specific enrichment",
      # Si se evaluó pero no dio significativo en ningún sitio
      TRUE                                        ~ "No significant enrichment / Lower"
    ),
    n_tejidos_sig  = ifelse(is.na(n_tejidos_sig), 0, n_tejidos_sig),
    tejidos        = ifelse(is.na(tejidos), "-", tejidos),
    mejor_p.adjust = ifelse(is.na(mejor_p.adjust), NA, mejor_p.adjust)
  ) %>%
  dplyr::arrange(categoria, desc(n_tejidos_sig))

write.csv2(resumen_interpretado,
           file = file.path(carpeta_expresion_im, "resumen_interpretado.csv"),
           row.names = FALSE)

message("\n¡Resúmenes generados y clasificados con éxito! Se incluyeron todas las intersecciones originales.")








# ========================================================================

# ========================================================================
# BLOQUE 3: ===== VISUALIZACIÖN DE LOS RESULTADOS ====

tabla_maestra_expresion <- read.csv2(
  "./Output/Piloto_Microbiota/Expresion_Tejidos_IM_MicroGWAS/Tabla_Maestra_Expresion_Pleiotropia.csv"
)

# ── Colores NP: distintos de los microbioma ───────────────────────────────────

paleta_im <- c(
  "#142157", "plum", "#E08214", "#2166AC", "#D6604D",
  "#08A045", "#8C510A", "#01665E", "#C51B7D", "#35978F",
  "#a7a7a7", "#40004B", "#ffff00", "#053061", "#67001F"
)

# --- 1. PREPARAR DATOS ---

datos_heatmap_completos <- tabla_maestra_expresion %>%
  dplyr::mutate(
    # Aplicamos limpieza de nombres
    Trait_Micro_Corto = acortar_nombres_microbioma(Trait_Nombre_2),
    Trait_IM_Nombre   = Trait_Nombre_1,
    
    Cluster_IM_corto    = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_Micro_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    
    # Nombre final para la columna del heatmap
    Par_Limpio = paste0(Trait_Micro_Corto, " (", Cluster_Micro_corto, ")",
                        " vs ",
                        Trait_IM_Nombre, " (", Cluster_IM_corto, ")"),
    
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

# --- 3. ANOTACIÓN DE COLUMNAS 

traits_por_par <- datos_heatmap_completos %>%
  dplyr::select(Par_Limpio, Trait_Micro_Corto, Trait_IM_Nombre) %>%
  dplyr::distinct()

# Dos columnas: una con los traits MicroGWAS cortos y otra con los traits Nervous/Psychiatric
annotation_col <- data.frame(
  Microbioma    = traits_por_par$Trait_Micro_Corto, 
  Enfermedad_IM = traits_por_par$Trait_IM_Nombre,
  row.names     = traits_por_par$Par_Limpio 
)

# Colores Microbioma:
niveles_micro <- sort(unique(traits_por_par$Trait_Micro_Corto))
colores_micro_c <- colores_micro_por_nombre[niveles_micro]

# Colores IM
niveles_IM <- sort(unique(traits_por_par$Trait_IM_Nombre))
colores_IM_c <- setNames(paleta_im, niveles_IM)
colores_IM_c <- colores_IM_c[!is.na(names(colores_IM_c))]

annotation_colors_completa <- list(
  Microbioma    = colores_micro_c,
  Enfermedad_IM = colores_IM_c
)


# --- 4. PALETA Y HEATMAP ---
paleta <- colorRampPalette(c("#EDF8B1", "#7FCDBB", "#2C7FB8", "#081D58"))(100)

# Al tener todas las columnas, necesitamos un PDF bastante grande
pdf("./Output/Gráficos/MicroGWAS/Heatmap_Expr_Tejidos_Interseccion_IM.pdf", height = 22, width = 30)
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

# ==== DOTPLOT ====
# --- 1. Preparar datos ---
datos_dotplot_expr <- tabla_maestra_expresion %>%
  dplyr::mutate(
    Trait_Micro_Corto = acortar_nombres_microbioma(Trait_Nombre_2),
    Trait_IM_Nombre   = Trait_Nombre_1,
    
    # 2. Acortar el texto "Cluster" a "Cl_" para que ocupe menos
    Cluster_IM_corto    = gsub(".*_Cluster_", "Cl_", Cluster_1),
    Cluster_Micro_corto = gsub(".*_Cluster_", "Cl_", Cluster_2),
    
    # 3. Eje X: enfermedad IM + cluster
    # (Si quisieras que salieran los dos, podrías usar la variable Par_Limpio del heatmap)
    Eje_X = paste0(Trait_IM_Nombre, "\n(", Cluster_IM_corto, ")"),
    
    # 4. Calcular el logaritmo del p-valor solo para los significativos y mayores
    log_pval = ifelse(direccion == "higher" & p.adjust < 0.05, -log10(p.adjust), NA)
  ) %>%
  # Quitar los que no son significativos
  dplyr::filter(!is.na(log_pval)) %>%
  # Agrupar por cada par de clusters
  dplyr::group_by(Cluster_1, Cluster_2) %>%
  # Ordenar de menor a mayor p-valor (los más significativos primero)
  dplyr::arrange(p.adjust) %>%
  # Quedarnos solo con los 10 tejidos más significativos por cada par
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
  facet_grid(~ Trait_Micro_Corto, scales = "free_x", space = "free_x") +
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
    title    = "Tissue Expression (higher, FDR < 0.05)",
    subtitle = "Microbiota traits vs Immune System Disease"
  )

# --- 3. Guardar ---
ggsave(
  filename  = "./Output/Gráficos/MicroGWAS/Dotplot_Expr_Tejidos_Interseccion_IM.pdf",
  width     = 20,
  height    = 15,
  dpi       = 300,
  limitsize = FALSE
)
# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# ESTUDIAMOS EXPRESION DE CADA GEN DEL CLUSTER  ----

# Cargar datos de intersección
lista_intersecciones_im <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds") 
nombres_pares <- names(lista_intersecciones_im)

# Iteramos directamente sobre los nombres de la lista
resultados_ks_log2 <- lapply(nombres_pares, function(nombre_par) {
  
  # Extraemos la información de este par
  datos_par <- lista_intersecciones_im[[nombre_par]]
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

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 4: EXPORTAR RESULTADOS   =====

carpeta_resultados_ks <- "./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_IM_MicroGWAS"
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


expresion_completa_ks <- readRDS("./Output/Piloto_Microbiota/Resultados_Exp_Tejidos_KS_IM_MicroGWAS/expresion_genes_por_tejido_all.rds")


# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# BLOQUE 5: VISUALIZACIÓN DE LOS RESULTADOS

##### EXPRESIÓN POR GEN ====

# Clusters únicos
clusters_unicos <- unique(expresion_completa_ks$cluster)
message("Clusters a representar: ", length(clusters_unicos))

# Paleta
paleta_heatmap <- colorRampPalette(c("#EDF8B1", "#7FCDBB", "#2C7FB8", "#081D58"))(100)

pdf("./Output/Gráficos/MicroGWAS/Heatmap_Exp_Tejidos_por_Cluster_KS_IM_MicroGWAS.pdf",
    height = 20, width = 20)

for (cl in clusters_unicos) {
  
  message("Procesando: ", cl)
  
  # Filtrar datos de este cluster
  datos_cl <- expresion_completa_ks %>%
    dplyr::filter(cluster == cl)
  
  # Etiqueta para el título
  trait1 <- unique(datos_cl$trait_1)
  trait2 <- unique(datos_cl$trait_2)
  jaccard_val <- unique(datos_cl$jaccard)
  
  # Construir matriz original (con NAs para visualización)
  matriz_cl <- datos_cl %>%
    dplyr::select(Gene.name, Tissue, log2_nTPM) %>%
    tidyr::pivot_wider(names_from  = Tissue,
                       values_from = log2_nTPM) %>%
    tibble::column_to_rownames("Gene.name") %>%
    as.matrix()
  
  # Matriz sin NAs para calcular distancias
  matriz_calculo <- matriz_cl
  matriz_calculo[is.na(matriz_calculo)] <- 0
  
  
  # Distancia de tejidos
  dist_tejidos <- as.dist(1 - cor(matriz_calculo, method = "spearman"))
  
  # Distancias de genes (filas) — usando la matriz sin NAs
  dist_genes <- dist(matriz_calculo, method = "euclidean")
  
  # Dentro del bucle, antes del pheatmap:
  n_genes <- nrow(matriz_cl)
  
  # Altura de celda dinámica: más genes = celdas más pequeñas
  cellheight_dinamico <- dplyr::case_when(
    n_genes <= 20  ~ 20,
    n_genes <= 50  ~ 12,
    n_genes <= 100 ~ 8,
    n_genes <= 200 ~ 5,
    TRUE           ~ 4
  )
  
  # Tamaño de fuente dinámico
  fontsize_row_dinamico <- dplyr::case_when(
    n_genes <= 20  ~ 9,
    n_genes <= 50  ~ 7,
    n_genes <= 100 ~ 6,
    n_genes <= 200 ~ 5,
    TRUE           ~ 4
  )
  
  message("  Genes: ", n_genes, 
          " — cellheight: ", cellheight_dinamico,
          " — fontsize: ", fontsize_row_dinamico)
  
  # Heatmap
  tryCatch({
    pheatmap::pheatmap(
      mat                      = matriz_cl,
      color                    = paleta_heatmap,
      cluster_rows             = TRUE,
      cluster_cols             = TRUE,
      clustering_distance_rows = dist_genes,
      clustering_distance_cols = dist_tejidos,
      clustering_method        = "average",
      treeheight_row           = 0,
      treeheight_col           = 0,
      show_rownames            = TRUE, # ifelse(nrow(matriz_cl) <= 50, TRUE, FALSE)
      show_colnames            = TRUE,
      fontsize_col             = 12,
      cellwidth                = 20,
      cellheight  = cellheight_dinamico,
      fontsize_row = fontsize_row_dinamico,
      main                     = paste0(trait1, " vs ", trait2,
                                        "\n(Jaccard = ", jaccard_val, 
                                        " | ", nrow(matriz_cl), " genes) --> log2(nTPM), KS test"),
      na_col                   = "whitesmoke",
      border_color             = "#F5F5F5"
    )
  }, error = function(e) {
    message("  ⚠️ Error en cluster: ", cl, " — ", e$message)
  })
}

dev.off()
message("PDF generado con ", length(clusters_unicos), " heatmaps!")












