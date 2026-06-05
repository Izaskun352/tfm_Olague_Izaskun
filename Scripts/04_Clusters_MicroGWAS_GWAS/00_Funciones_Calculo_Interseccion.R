
# LIBRERIAS

# INPUTS

# FUNCIONES PARA SACAR INTERSECCION ENTRE CLUSTERS

#------------------------------------------------------------------------------
# F1: Calcular índice de Jaccard entre dos vectores de genes
# ------------------------------------------------------------------------------

calc_jaccard <- function(grupo_A, grupo_B) {
  interseccion <- length(intersect(grupo_A, grupo_B))
  union        <- length(union(grupo_A, grupo_B))
  if (union == 0) return(0)
  return(interseccion / union)
}

# ------------------------------------------------------------------------------
# F2: Leer clusters desde una o varias carpetas
# ------------------------------------------------------------------------------

cargar_clusters <- function(
    carpetas,
    columna_gen  = "gene.gene",
    columna_ensg = "gene.ENSG",     # para extraer Ensembl IDs
    patron       = "\\.csv$"
) {
  archivos <- list.files(path = carpetas, pattern = patron, full.names = TRUE)
  if (length(archivos) == 0) stop("❌ No se encontraron archivos .csv en las carpetas indicadas.")
  message("Cargando ", length(archivos), " clusters...")
  
  lista_genes  <- list()
  lista_ensembl <- list()
  lista_datos  <- list()
  
  for (archivo in archivos) {
    nombre                <- tools::file_path_sans_ext(basename(archivo))
    datos                 <- read.csv2(archivo)
    lista_genes[[nombre]]   <- as.character(datos[[columna_gen]])
    lista_ensembl[[nombre]] <- gsub(";.*$", "", as.character(datos[[columna_ensg]]))
    lista_datos[[nombre]]   <- datos
  }
  
  list(genes = lista_genes, ensembl = lista_ensembl, datos = lista_datos)
}


# ------------------------------------------------------------------------------
# F3: Calcular solapamiento Jaccard para todos los pares de clusters
# ------------------------------------------------------------------------------

calcular_jaccard_pares <- function(
    lista_clusters,
    carpeta_salida   = NULL,
    nombre_archivo   = "Matriz_Jaccard_Completa.csv",
    umbral           = 0,             # solo guardar pares con Jaccard > umbral
    columna_gen      = "gene.gene",
    columna_senal    = "gene.senal_inicial",
    valor_senal      = "Si",
    fun_nombre_trait = NULL           # función opcional: cluster -> nombre legible
) {
  lista_genes <- lista_clusters$genes
  lista_datos <- lista_clusters$datos
  
  nombres_clusters <- names(lista_genes)
  pares_posibles   <- combn(nombres_clusters, 2, simplify = FALSE)
  
  message("Cruzando ", length(pares_posibles), " pares posibles...")
  
  lista_resultados <- list()
  
  for (par in pares_posibles) {
    c1 <- par[1]; c2 <- par[2]
    
    genes1  <- lista_genes[[c1]]
    genes2  <- lista_genes[[c2]]
    j_index <- calc_jaccard(genes1, genes2)
    
    if (j_index > umbral) {
      
      genes_comunes <- intersect(genes1, genes2)
      datos1        <- lista_datos[[c1]]
      datos2        <- lista_datos[[c2]]
      
      senal_c1 <- datos1 %>%
        dplyr::filter(.data[[columna_gen]] %in% genes_comunes,
                      .data[[columna_senal]] == valor_senal) %>%
        dplyr::pull(.data[[columna_gen]])
      
      senal_c2 <- datos2 %>%
        dplyr::filter(.data[[columna_gen]] %in% genes_comunes,
                      .data[[columna_senal]] == valor_senal) %>%
        dplyr::pull(.data[[columna_gen]])
      
      ambos   <- intersect(senal_c1, senal_c2)
      solo_c1 <- setdiff(senal_c1, senal_c2)
      solo_c2 <- setdiff(senal_c2, senal_c1)
      
      lista_resultados[[length(lista_resultados) + 1]] <- data.frame(
        Cluster_1                = c1,
        Trait_Nombre_1           = if (!is.null(fun_nombre_trait)) fun_nombre_trait(c1) else c1,
        Cluster_2                = c2,
        Trait_Nombre_2           = if (!is.null(fun_nombre_trait)) fun_nombre_trait(c2) else c2,
        Indice_Jaccard           = round(j_index, 3),
        Genes_Compartidos        = length(genes_comunes),
        Total_Genes_Union        = length(union(genes1, genes2)),
        Senal_Inicial_Compartida = length(union(senal_c1, senal_c2)),
        Senal_Cluster1           = length(senal_c1),
        Senal_Cluster2           = length(senal_c2),
        Senal_Ambos_Clusters     = length(ambos),
        Senal_Solo_Cluster1      = length(solo_c1),
        Senal_Solo_Cluster2      = length(solo_c2),
        Genes_Senal_C1           = paste(senal_c1, collapse = ";"),
        Genes_Senal_C2           = paste(senal_c2, collapse = ";"),
        Genes_Senal_Ambos        = paste(ambos,    collapse = ";"),
        stringsAsFactors         = FALSE
      )
    }
  }
  
  resultados <- dplyr::bind_rows(lista_resultados)
  
  message("✔️  Pares con Jaccard > ", umbral, ": ", nrow(resultados))
  
  if (!is.null(carpeta_salida)) {
    dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
    write.csv2(resultados,
               file.path(carpeta_salida, nombre_archivo),
               row.names = FALSE)
    message("✔️  Guardado en: ", file.path(carpeta_salida, nombre_archivo))
  }
  
  invisible(resultados)
}



# ------------------------------------------------------------------------------
# F4: Calcular lista de intersecciones para pares mixtos
# ------------------------------------------------------------------------------

calcular_intersecciones <- function(
    resultados_jaccard,              # data.frame output de calcular_jaccard_pares()
    lista_clusters,                  # output de cargar_clusters()
    patron_grupo_A   = "^ZSCO\\.",   # patrón que identifica un grupo (ej: MicroGWAS)
    umbral_jaccard   = 0.5,
    min_genes        = 5,            # mínimo de genes en la intersección para incluir el par
    archivo_salida   = NULL,         # ruta .rds donde guardar la lista (opcional)
    columna_gen      = "gene.gene",
    columna_ensg     = "gene.ENSG",
    columna_senal    = "gene.senal_inicial",
    valor_senal      = "Si"
) {
  lista_ensembl <- lista_clusters$ensembl
  lista_datos   <- lista_clusters$datos
  
  #  Filtrar pares mixtos (un cluster de cada grupo) 
  pares_mixtos <- resultados_jaccard %>%
    dplyr::filter(Indice_Jaccard >= umbral_jaccard) %>%
    dplyr::filter(
      (!grepl(patron_grupo_A, Cluster_1) &  grepl(patron_grupo_A, Cluster_2)) |
        ( grepl(patron_grupo_A, Cluster_1) & !grepl(patron_grupo_A, Cluster_2))
    )
  
  message("Pares mixtos encontrados: ", nrow(pares_mixtos))
  
  if (nrow(pares_mixtos) == 0) {
    message("⚠️  Ningún par mixto supera el umbral de Jaccard >= ", umbral_jaccard)
    return(invisible(list()))
  }
  
  lista_intersecciones <- list()
  
  for (i in 1:nrow(pares_mixtos)) {
    
    c1         <- pares_mixtos$Cluster_1[i]
    c2         <- pares_mixtos$Cluster_2[i]
    nombre_par <- paste0(c1, "_vs_", c2)
    
    genes_c1 <- lista_ensembl[[c1]]
    genes_c2 <- lista_ensembl[[c2]]
    
    if (is.null(genes_c1) || is.null(genes_c2)) {
      message("⚠️  No se encontraron genes para: ", nombre_par)
      next
    }
    
    genes_interseccion <- intersect(genes_c1, genes_c2)
    
    if (length(genes_interseccion) < min_genes) {
      message("⚠️  Menos de ", min_genes, " genes en la intersección: ", nombre_par)
      next
    }
    
    datos_c1 <- lista_datos[[c1]]
    datos_c2 <- lista_datos[[c2]]
    
    filas_c1 <- datos_c1 %>%
      dplyr::filter(gsub(";.*$", "", .data[[columna_ensg]]) %in% genes_interseccion) %>%
      dplyr::mutate(cluster_origen = c1, trait_origen = pares_mixtos$Trait_Nombre_1[i])
    
    filas_c2 <- datos_c2 %>%
      dplyr::filter(gsub(";.*$", "", .data[[columna_ensg]]) %in% genes_interseccion) %>%
      dplyr::mutate(cluster_origen = c2, trait_origen = pares_mixtos$Trait_Nombre_2[i])
    
    genes_tabla <- dplyr::bind_rows(filas_c1, filas_c2) %>%
      dplyr::select(
        all_of(c(columna_ensg, columna_gen, columna_senal)),
        cluster_origen, trait_origen
      ) %>%
      dplyr::group_by(.data[[columna_gen]]) %>%
      dplyr::summarise(
        gene.ENSG            = dplyr::first(.data[[columna_ensg]]),
        es_semilla_en_c1     = any(.data[[columna_senal]] == valor_senal & cluster_origen == c1),
        es_semilla_en_c2     = any(.data[[columna_senal]] == valor_senal & cluster_origen == c2),
        es_semilla_en_alguno = any(.data[[columna_senal]] == valor_senal),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        nombre_par     = nombre_par,
        Cluster_1      = c1,
        Cluster_2      = c2,
        Trait_Nombre_1 = pares_mixtos$Trait_Nombre_1[i],
        Trait_Nombre_2 = pares_mixtos$Trait_Nombre_2[i]
      )
    
    lista_intersecciones[[nombre_par]] <- list(
      c1            = c1,
      c2            = c2,
      jaccard       = pares_mixtos$Indice_Jaccard[i],
      trait_1       = pares_mixtos$Trait_Nombre_1[i],
      trait_2       = pares_mixtos$Trait_Nombre_2[i],
      genes_ensembl = genes_interseccion,
      genes_simbolo = unique(genes_tabla[[columna_gen]]),
      genes_tabla   = genes_tabla
    )
    
    message("✔️  ", nombre_par, " — ", length(genes_interseccion), " genes en intersección")
  }
  
  message("lista_intersecciones generada: ", length(lista_intersecciones), " pares")
  
  if (!is.null(archivo_salida)) {
    dir.create(dirname(archivo_salida), showWarnings = FALSE, recursive = TRUE)
    saveRDS(lista_intersecciones, file = archivo_salida)
    message("✔️  Guardada en: ", archivo_salida)
  }
  
  invisible(lista_intersecciones)
}


# ------------------------------------------------------------------------------
# F5: Barplot intersecciones: numero de genes + genes semilla
# ------------------------------------------------------------------------------

barplot_intersecciones <- function(
    pares_mixtos,
    colores_traits,              # named vector: EFO_id -> color (para Non-seed genes)
    tabla_traits,                # para extraer nombres legibles de traits
    archivo_salida,
    mapeo_clusters_micro = NULL, # diccionario renombrado clusters MicroGWAS
    mapeo_clusters_np    = NULL, # diccionario renombrado clusters NP
    columna_rasgo        = "Rasgo",
    columna_nombre       = "name",
    fun_acortar          = NULL,
    patron_efo           = "^(EFO_[0-9]+)_Cluster.*$",
    color_semilla_micro  = "#FF8C42",
    color_semilla_np     = "#4A90B8",
    color_semilla_ambos  = "#9B2335",
    titulo               = "Number of genes per pleiotropic intersection",
    ancho                = 14
  ) {
    efo_a_nombre <- setNames(tabla_traits[[columna_nombre]], tabla_traits[[columna_rasgo]])
    
    #  Aplicar mapeo de nombres y construir etiquetas 
    datos <- pares_mixtos %>%
      dplyr::mutate(
        # Renombrar clusters si hay diccionario
        Cluster_1_label = if (!is.null(mapeo_clusters_micro)) 
          dplyr::recode(Cluster_1, !!!mapeo_clusters_micro) 
        else Cluster_1,
        Cluster_2_label = if (!is.null(mapeo_clusters_np))    
          dplyr::recode(Cluster_2, !!!mapeo_clusters_np)    
        else Cluster_2,
        
        # Nombres cortos (sin prefijo EFO o ZSCO)
        cluster_1_corto = gsub("^EFO_[0-9]+_", "", Cluster_1_label),
        cluster_2_corto = gsub("^.*(?=Cluster_)", "", Cluster_2_label, perl = TRUE),
        
        # EFO del cluster MicroGWAS para los colores
        trait_micro_efo = gsub(patron_efo, "\\1", Cluster_1),
        
        # Nombre legible del trait NP
        trait_np_nombre = {
          np_id <- gsub("^ZSCO\\.", "", Cluster_2)
          np_id <- gsub("^([^_]+_[0-9]+)_Cluster.*$", "\\1", np_id)
          n <- efo_a_nombre[np_id]
          ifelse(is.na(n), np_id, n)
        },
        trait_np_nombre = if (!is.null(fun_acortar)) fun_acortar(trait_np_nombre) 
        else trait_np_nombre,
        
        trait_np_nombre = stringr::str_to_sentence(trait_np_nombre),
        
        # Nombre legible del trait MicroGWAS
        trait_micro_nombre = {
          n <- efo_a_nombre[trait_micro_efo]
          ifelse(is.na(n), trait_micro_efo, n)
        },
        trait_micro_nombre = if (!is.null(fun_acortar)) fun_acortar(trait_micro_nombre)
        else trait_micro_nombre,
        
        # Etiqueta final
        etiqueta = paste0(cluster_1_corto, " vs \n", trait_np_nombre, " (",cluster_2_corto, ")"
                          ),
        n_no_semilla = Genes_Compartidos - Senal_Inicial_Compartida
      )
    
    # Formato largo 
    datos_largo <- datos %>%
      dplyr::select(etiqueta, trait_micro_efo, trait_micro_nombre,  
                    Genes_Compartidos, Senal_Inicial_Compartida,
                    n_no_semilla, Senal_Solo_Cluster1,
                    Senal_Solo_Cluster2, Senal_Ambos_Clusters) %>%
      tidyr::pivot_longer(
        cols = c(n_no_semilla, Senal_Solo_Cluster1,
                      Senal_Solo_Cluster2, Senal_Ambos_Clusters),
        names_to  = "tipo",
        values_to = "n_genes"
      ) %>%
      dplyr::mutate(
        tipo = dplyr::case_when(
          tipo == "n_no_semilla"         ~ "Non-seed genes",
          tipo == "Senal_Solo_Cluster1"  ~ "Seed only MicroGWAS",
          tipo == "Senal_Solo_Cluster2"  ~ "Seed only NP",
          tipo == "Senal_Ambos_Clusters" ~ "Seed in both"
        ),
        tipo = factor(tipo, levels = c("Non-seed genes", "Seed only MicroGWAS",
                                       "Seed only NP", "Seed in both"))
      )
    
    #  Orden de clusters 
    orden <- datos %>%
      dplyr::arrange(trait_micro_efo, dplyr::desc(Genes_Compartidos)) %>%
      dplyr::pull(etiqueta)
    
    datos_largo$etiqueta <- factor(datos_largo$etiqueta, levels = rev(orden))
    
    #  Colores por etiqueta + tipo 
    ref_efo <- datos %>% dplyr::select(etiqueta, trait_micro_efo) %>% dplyr::distinct()
    
    colores_fill <- datos_largo %>%
      dplyr::mutate(
        color = dplyr::case_when(
          tipo == "Non-seed genes"      ~ colores_traits[trait_micro_efo],
          tipo == "Seed only MicroGWAS" ~ color_semilla_micro,
          tipo == "Seed only NP"        ~ color_semilla_np,
          tipo == "Seed in both"        ~ color_semilla_ambos
        ),
        fill_key = paste0(etiqueta, "_", tipo)
      ) %>%
      dplyr::select(fill_key, color) %>%
      dplyr::distinct()
    
    colores_vector <- setNames(colores_fill$color, colores_fill$fill_key)
    datos_largo    <- datos_largo %>%
      dplyr::mutate(fill_key = paste0(etiqueta, "_", tipo))
    
    # Leyenda 
    traits_micro <- datos %>%
      dplyr::select(trait_micro_efo, trait_micro_nombre) %>%
      dplyr::distinct()
    
    leyenda_colores <- c(
      setNames(colores_traits[traits_micro$trait_micro_efo], traits_micro$trait_micro_nombre),
      "Non-seed genes"      = "#D3D3D3",
      "Seed only MicroGWAS" = color_semilla_micro,
      "Seed only NP"        = color_semilla_np,
      "Seed in both"        = color_semilla_ambos
    )
    
    # Plot 
    n_pares <- dplyr::n_distinct(datos_largo$etiqueta)
    
    p <- ggplot(datos_largo, aes(y = etiqueta, x = n_genes, fill = fill_key)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      geom_text(
        data        = datos,
        aes(y       = etiqueta,
            x       = Genes_Compartidos,
            label   = paste0(Genes_Compartidos, " (", Senal_Inicial_Compartida, ")")),
        hjust       = -0.15,
        size        = 6,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = colores_vector, guide = "none") +
      ggnewscale::new_scale_fill() +
      geom_point(
        data        = data.frame(x = NA, y = NA, tipo = names(leyenda_colores)),
        aes(x = x, y = y, fill = tipo),
        shape       = 22,
        size        = 4,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = leyenda_colores, name = "") +
      scale_x_continuous(
        position = "top",
        limits   = c(0, max(datos$Genes_Compartidos) * 1.25),
        expand   = c(0.02, 0)
      ) +
      labs(title = titulo, x = "Number of genes", y = "") +
      theme_minimal() +
      theme(
        legend.position  = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(hjust = 0.5, face = "bold",
                                        margin = margin(t = 10, b = 10), size = 12),
        axis.line.x.top  = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x.top = element_line(color = "black"),
        axis.title.x.top = element_text(margin = margin(b = 10, t = 10), size = 11),
        axis.text.x.top  = element_text(size = 14,margin = margin(b = 5)),
        axis.text.y      = element_text(size = 12),
        plot.margin      = margin(t = 10, r = 80, b = 10, l = 10)
      )
    
    pdf(archivo_salida, width = ancho, height = n_pares * 0.5 + 3)
    print(p)
    dev.off()
    
    message("✔️  PDF generado con ", n_pares, " intersecciones: ", archivo_salida)
    invisible(p)
  }


