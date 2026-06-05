
#source("scripts/00_setup.R")

# -------------------------------------
# CALCULO ENRIQUECIMIENTO GOBP
# -------------------------------------

# ==============================================================================
# F1: Análisis InterPro para un único cluster
# ==============================================================================

analizar_InterPro_cluster <- function(
    archivo,
    carpeta_salida,
    universo_genes,
    TERM2GENE,                    
    TERM2NAME,                    
    tabla_traits    = NULL,
    columna_rasgo   = "Rasgo",
    columna_nombre  = "name",
    p_cutoff        = 0.05,
    q_cutoff        = 0.05,
    n_top_simplificado = 30       # top N dominios para la versión simplificada
) {
  nombre_archivo <- tools::file_path_sans_ext(basename(archivo))
  datos_cluster  <- read.csv2(archivo)
  
  # Extraer genes
  genes   <- as.character(datos_cluster$gene.ENSG)
  genes   <- genes[!is.na(genes) & genes != ""]
  n_genes <- length(genes)
  
  # Señal inicial
  n_senal_inicial <- sum(datos_cluster$gene.senal_inicial == "Si", na.rm = TRUE)
  
  # Nombre del trait
  nombre_trait_id    <- sub("_Cluster_.*$", "", nombre_archivo)
  nombre_trait_label <- nombre_trait_id
  
  if (!is.null(tabla_traits)) {
    label <- tabla_traits %>%
      dplyr::filter(.data[[columna_rasgo]] == nombre_trait_id) %>%
      dplyr::pull(.data[[columna_nombre]]) %>%
      dplyr::first()
    if (!is.na(label) && length(label) > 0) nombre_trait_label <- label
  }
  
  message("Analizando: ", nombre_archivo,
          " (", n_genes, " genes | trait: ", nombre_trait_label, ")")
  
  # Enriquecimiento InterPro 
  resultado <- tryCatch(
    clusterProfiler::enricher(
      gene          = genes,
      universe      = universo_genes,
      TERM2GENE     = TERM2GENE,
      TERM2NAME     = TERM2NAME,
      pAdjustMethod = "BH",
      pvalueCutoff  = p_cutoff,
      qvalueCutoff  = q_cutoff
    ),
    error = function(e) {
      message("  -> ❌ Error en enricher(): ", e$message)
      NULL
    }
  )
  
  if (is.null(resultado) ||
      nrow(resultado@result %>% dplyr::filter(p.adjust < p_cutoff)) == 0) {
    message("  -> ⚠️  Sin dominios InterPro significativos.")
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
  
  # Tabla completa 
  tabla_completa <- as.data.frame(resultado) %>%
    dplyr::mutate(!!!cols_extra)
  
  write.csv2(tabla_completa,
             file.path(carpeta_salida, paste0("InterPro_", nombre_archivo, ".csv")),
             row.names = FALSE)
  
  # Tabla simplificada (top N por p.adjust y GeneRatio) 
  tabla_simplificada <- tabla_completa %>%
    dplyr::mutate(
      GeneRatio_val = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      })
    ) %>%
    dplyr::arrange(p.adjust, dplyr::desc(GeneRatio_val)) %>%
    dplyr::slice_head(n = n_top_simplificado) %>%
    dplyr::select(-GeneRatio_val)
  
  write.csv2(tabla_simplificada,
             file.path(carpeta_salida, paste0("InterPro_Simplificado_", nombre_archivo, ".csv")),
             row.names = FALSE)
  
  message("  -> ✔️  Dominios InterPro: ", nrow(tabla_completa),
          " totales | ", nrow(tabla_simplificada), " simplificados.")
  
  list(completa = tabla_completa, simplificada = tabla_simplificada)
}


# ==============================================================================
# F2: Análisis InterPro para todos los clusters de una carpeta
# ==============================================================================

analizar_InterPro_carpeta <- function(
    carpeta_clusters,
    carpeta_salida     = NULL,
    universo_genes,
    TERM2GENE,
    TERM2NAME,
    tabla_traits       = NULL,
    columna_rasgo      = "Rasgo",
    columna_nombre     = "name",
    p_cutoff           = 0.05,
    q_cutoff           = 0.05,
    n_top_simplificado = 30
) {
  # Validaciones 
  if (!dir.exists(carpeta_clusters)) {
    stop("❌ La carpeta de clusters no existe: ", carpeta_clusters)
  }
  
  archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)
  if (length(archivos_clusters) == 0) {
    stop("❌ No se encontraron archivos .csv en: ", carpeta_clusters)
  }
  
  if (is.null(carpeta_salida)) {
    carpeta_salida <- file.path(dirname(carpeta_clusters),
                                paste0("InterPro_", basename(carpeta_clusters)))
  }
  dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
  
  message("Universo de genes: ",   length(universo_genes), " genes")
  message("Clusters a analizar: ", length(archivos_clusters))
  
  #  Bucle principal
  lista_completas     <- list()
  lista_simplificadas <- list()
  
  for (archivo in archivos_clusters) {
    
    resultado <- tryCatch(
      analizar_InterPro_cluster(
        archivo            = archivo,
        carpeta_salida     = carpeta_salida,
        universo_genes     = universo_genes,
        TERM2GENE          = TERM2GENE,
        TERM2NAME          = TERM2NAME,
        tabla_traits       = tabla_traits,
        columna_rasgo      = columna_rasgo,
        columna_nombre     = columna_nombre,
        p_cutoff           = p_cutoff,
        q_cutoff           = q_cutoff,
        n_top_simplificado = n_top_simplificado
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
  
  # Archivos maestros 
  message("\n¡Proceso terminado!")
  message("Clusters con enriquecimiento InterPro significativo: ", length(lista_completas))
  
  if (length(lista_completas) > 0) {
    
    dplyr::bind_rows(lista_completas) %>%
      write.csv2(file.path(carpeta_salida, "InterPro_Resumen_Completo.csv"),
                 row.names = FALSE)
    
    dplyr::bind_rows(lista_simplificadas) %>%
      write.csv2(file.path(carpeta_salida, "InterPro_Resumen_Simplificado.csv"),
                 row.names = FALSE)
    
    message("Resúmenes globales guardados en: ", carpeta_salida)
    message("Nombre de los archivos, : InterPro_Resumen_Completo.csv ; InterPro_Resumen_Simplificado.csv", carpeta_salida)
    
  } else {
    message("Ningún cluster arrojó resultados InterPro significativos.")
  }
  
  invisible(list(completa = lista_completas, simplificada = lista_simplificadas))
}

# ==============================================================================
# F3: Análisis InterPro para una intersección de clusters
# ==============================================================================

analizar_InterPro_intersecciones <- function(
    lista_intersecciones,
    carpeta_salida,
    universo_genes,
    TERM2GENE,
    TERM2NAME,
    p_cutoff           = 0.05,
    q_cutoff           = 0.05,
    n_top_simplificado = 30
) {
  dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
  
  lista_completas     <- list()
  lista_simplificadas <- list()
  
  for (nombre_par in names(lista_intersecciones)) {
    inter <- lista_intersecciones[[nombre_par]]
    genes <- inter$genes_ensembl
    
    message("Analizando: ", nombre_par, " (", length(genes), " genes)")
    
    # Enriquecimiento InterPro 
    resultado <- tryCatch(
      clusterProfiler::enricher(
        gene   = genes,
        universe  = universo_genes,
        TERM2GENE  = TERM2GENE,
        TERM2NAME  = TERM2NAME,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_cutoff,
        qvalueCutoff  = q_cutoff
      ),
      error = function(e) {
        message("  -> ❌ Error en enricher(): ", e$message)
        NULL
      }
    )
    
    if (is.null(resultado) ||
        nrow(resultado@result %>% dplyr::filter(p.adjust < p_cutoff)) == 0) {
      message("  -> ⚠️  Sin dominios InterPro significativos.")
      next
    }
    
    #  Columnas extra 
    cols_extra <- list(
      Par     = nombre_par,
      Cluster_1  = inter$c1,
      Cluster_2  = inter$c2,
      Trait_Nombre_1 = inter$trait_1,
      Trait_Nombre_2 = inter$trait_2,
      Jaccard  = inter$jaccard,
      N_genes = length(genes)
    )
    
    #  Tabla completa
    tabla_completa <- as.data.frame(resultado) %>%
      dplyr::mutate(!!!cols_extra)
    
    write.csv2(tabla_completa,
               file.path(carpeta_salida, paste0("InterPro_", nombre_par, ".csv")),
               row.names = FALSE)
    
    # Tabla simplificada (top N por p.adjust y GeneRatio) 
    tabla_simplificada <- tabla_completa %>%
      dplyr::mutate(
        GeneRatio_val = sapply(GeneRatio, function(x) {
          parts <- strsplit(x, "/")[[1]]
          as.numeric(parts[1]) / as.numeric(parts[2])
        })
      ) %>%
      dplyr::arrange(p.adjust, dplyr::desc(GeneRatio_val)) %>%
      dplyr::slice_head(n = n_top_simplificado) %>%
      dplyr::select(-GeneRatio_val)
    
    write.csv2(tabla_simplificada,
               file.path(carpeta_salida, paste0("InterPro_Simplificado_", nombre_par, ".csv")),
               row.names = FALSE)
    
    lista_completas[[nombre_par]]     <- tabla_completa
    lista_simplificadas[[nombre_par]] <- tabla_simplificada
    
    message("  -> ✔️  ", nrow(tabla_completa), " dominios totales | ",
            nrow(tabla_simplificada), " simplificados.")
  }
  
  #  Archivos maestros 
  message("\n¡Proceso terminado!")
  message("Pares con enriquecimiento InterPro significativo: ", length(lista_completas))
  
  if (length(lista_completas) > 0) {
    dplyr::bind_rows(lista_completas) %>%
      write.csv2(file.path(carpeta_salida, "InterPro_Resumen_Completo.csv"),
                 row.names = FALSE)
    dplyr::bind_rows(lista_simplificadas) %>%
      write.csv2(file.path(carpeta_salida, "InterPro_Resumen_Simplificado.csv"),
                 row.names = FALSE)
    message("✔️  Resúmenes globales guardados en: ", carpeta_salida)
  } else {
    message("⚠️  Ningún par arrojó resultados InterPro significativos.")
  }
  
  invisible(list(completa = lista_completas, simplificada = lista_simplificadas))
}

# -------------------------------------
# REPRESENTACION GRÄFICA
# -------------------------------------

