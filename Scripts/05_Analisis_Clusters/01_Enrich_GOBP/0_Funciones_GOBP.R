
# LIBRERIAS

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# -------------------------------------
# CALCULO ENRIQUECIMIENTO GOBP
# -------------------------------------

# Funcion para calcular enriquecimiento GOBP de un solo cluster

analizar_GO_cluster <- function(
    archivo,
    carpeta_salida,
    universo_genes,
    tabla_traits    = NULL,       # data.frame con Rasgo y name para obtener el nombre del trait
    columna_rasgo   = "Rasgo",
    columna_nombre  = "name",
    org_db          = org.Hs.eg.db,
    ontologia       = "BP",
    p_cutoff        = 0.05,
    q_cutoff        = 0.05,
    simplify_cutoff = 0.7
) {
  
  nombre_archivo   <- tools::file_path_sans_ext(basename(archivo))
  datos_cluster    <- read.csv2(archivo)
  
  # Extraer genes
  genes <- as.character(datos_cluster$gene.ENSG)
  genes <- genes[!is.na(genes) & genes != ""]
  n_genes <- length(genes)
  
  # Señal inicial
  n_senal_inicial <- sum(datos_cluster$gene.senal_inicial == "Si", na.rm = TRUE)
  
  # Nombre del trait
  nombre_trait_id    <- sub("_Cluster_.*$", "", nombre_archivo)
  nombre_trait_label <- nombre_trait_id  # fallback por defecto
  
  if (!is.null(tabla_traits)) {
    label <- tabla_traits %>%
      filter(.data[[columna_rasgo]] == nombre_trait_id) %>%
      pull(.data[[columna_nombre]]) %>%
      dplyr::first()
    
    if (!is.na(label) && length(label) > 0) nombre_trait_label <- label
  }
  
  message("Analizando: ", nombre_archivo,
          " (", n_genes, " genes | trait: ", nombre_trait_label, ")")
  
  # Análisis GO 
  resultado_go <- enrichGO(
    gene          = genes,
    universe      = universo_genes,
    OrgDb         = org_db,
    keyType       = "ENSEMBL",
    ont           = ontologia,
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff
  )
  
  # Sin resultados significativos
  if (is.null(resultado_go) ||
      nrow(resultado_go@result %>% filter(p.adjust < p_cutoff)) == 0) {
    message("  -> ⚠️  Sin términos GO significativos.")
    return(NULL)
  }
  
  # Columnas extra comunes
  cols_extra <- list(
    Cluster_Origen  = nombre_archivo,
    Trait_ID        = nombre_trait_id,
    Trait_Nombre    = nombre_trait_label,
    N_genes_cluster = n_genes,
    N_senal_inicial = n_senal_inicial
  )
  
  #  Tabla completa 
  tabla_completa <- as.data.frame(resultado_go) %>%
    mutate(!!!cols_extra)
  
  write.csv2(tabla_completa,
             file.path(carpeta_salida, paste0("GO_", nombre_archivo, ".csv")),
             row.names = FALSE)
  
  #  Tabla simplificada 
  resultado_simplificado <- tryCatch(
    simplify(resultado_go, cutoff = simplify_cutoff,
             by = "p.adjust", select_fun = min),
    error = function(e) resultado_go
  )
  
  tabla_simplificada <- as.data.frame(resultado_simplificado) %>%
    mutate(!!!cols_extra)
  
  write.csv2(tabla_simplificada,
             file.path(carpeta_salida, paste0("GO_Simplificado_", nombre_archivo, ".csv")),
             row.names = FALSE)
  
  message("  -> ✔️  Términos GO: ", nrow(tabla_completa),
          " originales | ", nrow(tabla_simplificada), " simplificados.")
  
  # Devuelve ambas tablas
  list(completa = tabla_completa, simplificada = tabla_simplificada)
}


# Funcion para calcular enriquecimiento GOBP de todos los clusters de una carpeta

analizar_GO_carpeta <- function(
    carpeta_clusters,
    carpeta_salida  = NULL,       # si es NULL, crea subcarpeta automáticamente
    universo_genes,
    tabla_traits    = NULL,
    columna_rasgo   = "Rasgo",
    columna_nombre  = "name",
    org_db          = org.Hs.eg.db,
    ontologia       = "BP",
    p_cutoff        = 0.05,
    q_cutoff        = 0.05,
    simplify_cutoff = 0.7
) {
  #  Validaciones 
  if (!dir.exists(carpeta_clusters)) {
    stop("❌ La carpeta de clusters no existe: ", carpeta_clusters)
  }
  
  archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(archivos_clusters) == 0) {
    stop("❌ No se encontraron archivos .csv en: ", carpeta_clusters)
  }
  
  # Carpeta de salida automática si no se especifica
  if (is.null(carpeta_salida)) {
    carpeta_salida <- file.path(dirname(carpeta_clusters), 
                                paste0("GO_", basename(carpeta_clusters)))
  }
  dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
  
  message("Universo de genes: ",  length(universo_genes), " genes")
  message("Clusters a analizar: ", length(archivos_clusters))
  
  #  Bucle principal 
  lista_completas     <- list()
  lista_simplificadas <- list()
  
  for (archivo in archivos_clusters) {
    
    resultado <- tryCatch(
      analizar_GO_cluster(
        archivo         = archivo,
        carpeta_salida  = carpeta_salida,
        universo_genes  = universo_genes,
        tabla_traits    = tabla_traits,
        columna_rasgo   = columna_rasgo,
        columna_nombre  = columna_nombre,
        org_db          = org_db,
        ontologia       = ontologia,
        p_cutoff        = p_cutoff,
        q_cutoff        = q_cutoff,
        simplify_cutoff = simplify_cutoff
      ),
      error = function(e) {
        message("  -> ❌ Error en ", basename(archivo), ": ", e$message)
        NULL
      }
    )
    
    if (!is.null(resultado)) {
      nombre <- tools::file_path_sans_ext(basename(archivo))
      lista_completas[[nombre]]     <- resultado$completa
      lista_simplificadas[[nombre]] <- resultado$simplificada
    }
  }
  
  #  Archivos maestros
  message("\n¡Proceso terminado!")
  message("Clusters con GO significativo: ", length(lista_completas))
  
  if (length(lista_completas) > 0) {
    
    tabla_maestra_completa <- bind_rows(lista_completas)
    write.csv2(tabla_maestra_completa,
               file.path(carpeta_salida, "GO_Resumen_Completo.csv"),
               row.names = FALSE)
    
    tabla_maestra_simplificada <- bind_rows(lista_simplificadas)
    write.csv2(tabla_maestra_simplificada,
               file.path(carpeta_salida, "GO_Resumen_Simplificado.csv"),
               row.names = FALSE)
    
    message("✔️  Resúmenes globales guardados en: ", carpeta_salida)
    
  } else {
    message("⚠️  Ningún cluster arrojó resultados GO significativos.")
  }
  
  invisible(list(completa = lista_completas, simplificada = lista_simplificadas))
}

# Funcion para calcular GOBP de la interseccion

analizar_GO_intersecciones <- function(
    lista_intersecciones,
    carpeta_salida,
    universo_genes,
    tabla_traits    = NULL,
    columna_rasgo   = "Rasgo",
    columna_nombre  = "name",
    org_db          = org.Hs.eg.db,
    ontologia       = "BP",
    p_cutoff        = 0.05,
    q_cutoff        = 0.05,
    simplify_cutoff = 0.7
) {
  dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
  
  lista_completas     <- list()
  lista_simplificadas <- list()
  
  for (nombre_par in names(lista_intersecciones)) {
    inter <- lista_intersecciones[[nombre_par]]
    genes <- inter$genes_ensembl
    
    message("Analizando: ", nombre_par, " (", length(genes), " genes)")
    
    resultado_go <- enrichGO(
      gene          = genes,
      universe      = universo_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENSEMBL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
    
    if (is.null(resultado_go) ||
        nrow(resultado_go@result %>% filter(p.adjust < 0.05)) == 0) {
      message("  -> ⚠️  Sin términos GO significativos.")
      next
    }
    
    cols_extra <- list(
      Par_Origen     = nombre_par,
      Cluster_1      = inter$c1,
      Cluster_2      = inter$c2,
      Trait_Nombre_1 = inter$trait_1,
      Trait_Nombre_2 = inter$trait_2,
      Jaccard        = inter$jaccard,
      N_genes        = length(genes)
    )
    
    tabla_completa <- as.data.frame(resultado_go) %>% mutate(!!!cols_extra)
    write.csv2(tabla_completa,
               file.path(carpeta_salida, paste0("GO_", nombre_par, ".csv")),
               row.names = FALSE)
    
    resultado_simplificado <- tryCatch(
      simplify(resultado_go, cutoff = 0.7, by = "p.adjust", select_fun = min),
      error = function(e) resultado_go
    )
    tabla_simplificada <- as.data.frame(resultado_simplificado) %>% mutate(!!!cols_extra)
    write.csv2(tabla_simplificada,
               file.path(carpeta_salida, paste0("GO_Simplificado_", nombre_par, ".csv")),
               row.names = FALSE)
    
    lista_completas[[nombre_par]]     <- tabla_completa
    lista_simplificadas[[nombre_par]] <- tabla_simplificada
    
    message("  -> ✔️  ", nrow(tabla_completa), " términos | ",
            nrow(tabla_simplificada), " simplificados")
  }
  
  # Maestros globales
  if (length(lista_completas) > 0) {
    bind_rows(lista_completas) %>%
      write.csv2(file.path(carpeta_salida, "GO_Resumen_Completo.csv"), row.names = FALSE)
    bind_rows(lista_simplificadas) %>%
      write.csv2(file.path(carpeta_salida, "GO_Resumen_Simplificado.csv"), row.names = FALSE)
    message("✔️  Resúmenes globales guardados en: ", carpeta_salida)
  }
  
  invisible(list(completa = lista_completas, simplificada = lista_simplificadas))
}

# -------------------------------------
# REPRESENTACION GRAFICA
# -------------------------------------

# HEATMAP DE GOBP

heatmap_GO_clusters <- function(
    archivo_maestro,                    # ruta al CSV maestro simplificado
    archivo_salida,                     # ruta del PDF de salida
    tabla_traits,                       # data.frame con Rasgo y name
    colores_traits,                     # named vector: c("EFO_xxx" = "#color", ...)
    n_top,                # top N términos GO por cluster
    columna_rasgo = "Rasgo",
    columna_nombre  = "name",
    fun_acortar     = NULL,             # función opcional para acortar nombres de traits
    patron_efo = "^(EFO_[0-9]+)_Cluster.*$",  # patrón para extraer el ID del trait
    mapeo_clusters   = NULL,  
    altura_pdf     = 20,
    anchura_pdf  = 20,
    fontsize_row  = 12,
    fontsize_col = 10,
    cellheight = 12,
    cellwidth  = 17,
    titulo
) {
  # Validaciones 
  if (!file.exists(archivo_maestro)) {
    stop("❌ No se encuentra el archivo maestro: ", archivo_maestro)
  }
  
  cols_necesarias <- c(columna_rasgo, columna_nombre)
  cols_faltantes  <- setdiff(cols_necesarias, colnames(tabla_traits))
  if (length(cols_faltantes) > 0) {
    stop("❌ Columnas faltantes en 'tabla_traits': ", paste(cols_faltantes, collapse = ", "))
  }
  
  dir.create(dirname(archivo_salida), showWarnings = FALSE, recursive = TRUE)
  
  # Cargar datos 
  go_maestro <- read.csv2(archivo_maestro)
  message("Archivo maestro cargado: ", nrow(go_maestro), " filas")
  if (!is.null(mapeo_clusters)) {
    go_maestro <- go_maestro %>%
      dplyr::mutate(Cluster_Origen = dplyr::recode(Cluster_Origen, !!!mapeo_clusters))
  }
  
  #  Paleta de colores de traits
  efo_a_nombre <- tabla_traits %>%
    dplyr::filter(.data[[columna_rasgo]] %in% names(colores_traits)) %>%
    dplyr::select(all_of(c(columna_rasgo, columna_nombre))) %>%
    dplyr::distinct() %>%
    tibble::deframe()
  
  nombres_colores <- efo_a_nombre
  if (!is.null(fun_acortar)) nombres_colores <- fun_acortar(nombres_colores)
  
  colores_globales     <- setNames(colores_traits[names(efo_a_nombre)], nombres_colores)
  annotation_colors    <- list(Trait = colores_globales)
  
  # Top N términos GO por cluster 
  top_terms <- go_maestro %>%
    dplyr::group_by(Cluster_Origen) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup() %>%
    dplyr::pull(Description) %>%
    unique()
  
  message("Términos GO únicos seleccionados: ", length(top_terms))
  
  # Construir matriz 
  matriz_go <- go_maestro %>%
    dplyr::filter(Description %in% top_terms) %>%
    dplyr::mutate(log_pval = -log10(p.adjust)) %>%
    dplyr::select(Description, Cluster_Origen, log_pval) %>%
    tidyr::pivot_wider(names_from  = Cluster_Origen,
                       values_from = log_pval,
                       values_fill = 0) %>%
    tibble::column_to_rownames("Description")
  
  message("Dimensiones matriz: ", nrow(matriz_go), " términos x ", ncol(matriz_go), " clusters")
  
  
  #  Anotación de columnas 
  efo_por_cluster <- gsub(patron_efo, "\\1", colnames(matriz_go))
  efo_a_nombre_full <- setNames(tabla_traits[[columna_nombre]], tabla_traits[[columna_rasgo]])
  nombres_reales <- efo_a_nombre_full[efo_por_cluster]
  
  if (!is.null(fun_acortar)) nombres_reales <- fun_acortar(nombres_reales)
  
  annotation_col <- data.frame(
    Trait     = nombres_reales,
    row.names = colnames(matriz_go)
  )
  
  # Extraer n_genes por cluster desde el maestro
  n_genes_por_cluster <- go_maestro %>%
    dplyr::select(Cluster_Origen, N_genes_cluster) %>%
    dplyr::distinct()
  
  colnames_cortos <- sapply(colnames(matriz_go), function(cl) {
    n <- n_genes_por_cluster$N_genes_cluster[n_genes_por_cluster$Cluster_Origen == cl]
    n <- if (length(n) == 0) "?" else n[1]
    paste0(gsub("^EFO_[0-9]+_", "", cl), " (n=", n, ")")
  })
  
  # Paleta del heatmap 
  paleta <- colorRampPalette(c("#EEF9C4", "#7FCDBB", "#2C7FB8", "#1A4E88"))(100)
  
  breaks <- seq(0, 60, length.out = 101)
  
  #  Generar heatmap 
  pdf(archivo_salida, height = altura_pdf, width = anchura_pdf)
  print(pheatmap(
    mat               = as.matrix(matriz_go),
    color             = paleta,
    breaks = breaks,
    annotation_col    = annotation_col,
    annotation_colors = annotation_colors,
    labels_col        = colnames_cortos,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    treeheight_row    = 0,
    treeheight_col    = 0,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = fontsize_row,
    fontsize_col      = fontsize_col,
    cellheight        = cellheight,
    cellwidth         = cellwidth,
    angle_col         = 45,
    main              = titulo,
    border_color      = "white"
  ))
  dev.off()
  
  message("✔️  Heatmap guardado en: ", archivo_salida)
  invisible(matriz_go)
}

# DOTPLOT DE GOBP

dotplot_GO_clusters <- function(
    go_maestro,                         # data.frame maestro (ya cargado) o ruta a CSV
    archivo_salida,                     # ruta del PDF de salida
    tabla_traits,                       # data.frame con Rasgo y name
    n_top  = 2,                # top N términos GO por cluster
    columna_rasgo   = "Rasgo",
    columna_nombre  = "name",
    fun_acortar     = NULL,             # función opcional para acortar nombres de traits
    patron_efo      = "^(EFO_[0-9]+)_Cluster.*$",
    titulo,
    mapeo_clusters  = NULL,
    ancho = 30,
    alto  = 18,
    dpi = 300,
    size_texto_x  = 22,
    size_texto_y  = 22,
    size_strip  = 14,
    rango_puntos = c(6, 12)
) {
  # Cargar datos si se pasa una ruta 
  if (is.character(go_maestro)) {
    if (!file.exists(go_maestro)) stop("❌ No se encuentra el archivo: ", go_maestro)
    go_maestro <- read.csv2(go_maestro)
    message("Archivo maestro cargado: ", nrow(go_maestro), " filas")
  }
  
  dir.create(dirname(archivo_salida), showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(mapeo_clusters)) {
    go_maestro <- go_maestro %>%
      dplyr::mutate(Cluster_Origen = dplyr::recode(Cluster_Origen, !!!mapeo_clusters))
  }
  
  # Mapeo EFO → nombre 
  efo_a_nombre <- setNames(tabla_traits[[columna_nombre]], tabla_traits[[columna_rasgo]])
  
  #  Preparar datos 
  datos_dotplot <- go_maestro %>%
    dplyr::mutate(
      Cluster_corto = gsub("^EFO_[0-9]+_", "", Cluster_Origen),
      
      GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)),
      GeneRatio_den = as.numeric(sub(".*/", "", GeneRatio)),
      GeneRatio_val = GeneRatio_num / GeneRatio_den,
      
      efo_id       = gsub(patron_efo, "\\1", Cluster_Origen),
      Trait_Nombre  = efo_a_nombre[efo_id],
      Trait_Corto   = if (!is.null(fun_acortar)) fun_acortar(Trait_Nombre) else Trait_Nombre,
      Cluster_label = paste0("(", Cluster_corto, ", n = ", N_genes_cluster, ")")
    ) %>%
    dplyr::group_by(Cluster_Origen) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup()
  
  # Ordenar términos GO 
  datos_dotplot$Description <- factor(
    datos_dotplot$Description,
    levels = unique(datos_dotplot$Description[order(datos_dotplot$p.adjust, decreasing = TRUE)])
  )
  
  # Dotplot 
  p <- ggplot(datos_dotplot, aes(x = Cluster_label, y = Description)) +
    geom_point(aes(size = GeneRatio_val, color = p.adjust)) +
    scale_color_viridis_c(
      option    = "magma",
      direction = 1,
      trans     = "log10",
      name      = "p.adjust"
    ) +
    scale_size_continuous(name = "GeneRatio", range = rango_puntos) +
    facet_grid(~ Trait_Corto, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = size_texto_x),
      axis.text.y        = element_text(size = size_texto_y),
      axis.title         = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 20, face = "bold"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = size_strip, face = "bold", color = "black"),
      panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      panel.spacing      = unit(0.05, "lines"),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    ) +
    labs(title = titulo)
  
  # Guardar 
  ggsave(filename = archivo_salida, plot = p, width = ancho, height = alto, dpi = dpi)
  message("✔️  Dotplot guardado en: ", archivo_salida)
  
  invisible(p)
}







