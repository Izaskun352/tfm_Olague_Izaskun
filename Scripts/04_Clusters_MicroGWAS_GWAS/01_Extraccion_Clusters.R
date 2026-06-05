
# LIBRERIAS

source("scripts/00_setup.R")
source("Scripts/Renombrar_Clusters.R")

## EXTRACCION CLUSTERS

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# FUNCIONES

# ---- Funcion para extraer los clusters

extraer_clusters <- function(
    carpeta_entrada,
    carpeta_salida,
    columna_filtro   = "Selected.KS",
    valor_filtro     = "1",
    columna_cluster  = "cluster.walktrap",
    columnas_genes   = c("ENSG", "gene", "padj")
) {
  # Validaciones 
  if (!dir.exists(carpeta_entrada)) {
    stop("❌ La carpeta de entrada no existe: ", carpeta_entrada)
  }
  
  archivos_traits <- list.files(
    path       = carpeta_entrada,
    pattern    = "\\.[Rr][Dd][Ss]$",
    full.names = TRUE
  )
  
  if (length(archivos_traits) == 0) {
    stop("❌ No se encontraron archivos .rds en: ", carpeta_entrada)
  }
  
  #  Crear carpeta de salida 
  dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
  
  total_clusters <- 0
  
  # Bucle principal 
  for (archivo in archivos_traits) {
    
    nombre_trait  <- tools::file_path_sans_ext(basename(archivo))
    datos_trait   <- as.data.frame(readRDS(archivo))
    
    # Verificar que las columnas necesarias existen
    cols_necesarias <- c(columna_filtro, columna_cluster, columnas_genes)
    cols_faltantes  <- setdiff(cols_necesarias, colnames(datos_trait))
    
    if (length(cols_faltantes) > 0) {
      message("⚠️ Columnas faltantes en '", nombre_trait, "': ",
              paste(cols_faltantes, collapse = ", "), " — saltando.")
      next
    }
    
    # Filtrar clusters significativos
    datos_filtrados <- datos_trait %>%
      filter(.data[[columna_filtro]] == valor_filtro)
    
    if (nrow(datos_filtrados) == 0) {
      message("⚠️  Ningún cluster significativo para: ", nombre_trait)
      next
    }
    
    # Dividir por cluster
    lista_clusters <- split(
      datos_filtrados[, columnas_genes, drop = FALSE],
      datos_filtrados[[columna_cluster]]
    )
    
    # Guardar cada cluster
    for (id_cluster in names(lista_clusters)) {
      
      id_limpio        <- gsub(";", ".", id_cluster)
      genes_del_cluster <- lista_clusters[[id_cluster]]
      
      # Añadir columna señal inicial si 'padj' está disponible
      if ("padj" %in% columnas_genes) {
        genes_del_cluster <- genes_del_cluster %>%
          mutate(senal_inicial = ifelse(padj == 100, "Si", "No"))
      }
      
      nombre_archivo <- file.path(
        carpeta_salida,
        paste0(nombre_trait, "_Cluster_", id_limpio, ".csv")
      )
      
      write.csv2(data.frame(gene = genes_del_cluster),
                 file = nombre_archivo, row.names = FALSE)
      
      total_clusters <- total_clusters + 1
    }
    
    message("✔️  Archivos generados para: ", nombre_trait)
  }
  
  # Resumen final 
  message("¡Proceso terminado con éxito!")
  message("¡Se han generado un total de ", total_clusters, " clusters significativos!")
  
  invisible(total_clusters)
}

# --- Funcion para extraer resultados propagacion de traits de un área específica

copiar_traits_por_area <- function(
    carpeta_origen,
    carpeta_destino,
    tabla_traits,                                        # data.frame con los traits y sus áreas
    columna_rasgo   = "Rasgo",                          # columna con el nombre del trait
    columna_areas   = "therapeuticAreas",               # columna con las áreas terapéuticas
    areas_filtro, # IDs de las áreas a filtrar
    sobreescribir   = FALSE
) {
  # Validaciones
  if (!dir.exists(carpeta_origen)) {
    stop("❌ La carpeta de origen no existe: ", carpeta_origen)
  }
  
  if (!is.data.frame(tabla_traits)) {
    stop("❌ 'tabla_traits' debe ser un data.frame.")
  }
  
  cols_necesarias <- c(columna_rasgo, columna_areas)
  cols_faltantes  <- setdiff(cols_necesarias, colnames(tabla_traits))
  if (length(cols_faltantes) > 0) {
    stop("❌ Columnas faltantes en 'tabla_traits': ", paste(cols_faltantes, collapse = ", "))
  }
  
  archivos_GWAS <- list.files(carpeta_origen, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE)
  
  if (length(archivos_GWAS) == 0) {
    stop("❌ No se encontraron archivos .rds en: ", carpeta_origen)
  }
  
  #  Crear carpeta de destino
  dir.create(carpeta_destino, showWarnings = FALSE, recursive = TRUE)
  
  # Filtrar traits por área terapéutica 
  patron_areas <- paste(areas_filtro, collapse = "|")
  
  traits_objetivo <- tabla_traits %>%
    filter(grepl(patron_areas, .data[[columna_areas]], ignore.case = TRUE)) %>%
    pull(.data[[columna_rasgo]]) %>%
    trimws() %>%
    unique()
  
  message("Traits encontrados: ", length(traits_objetivo))
  
  if (length(traits_objetivo) == 0) {
    message("⚠️  Ningún trait coincide con las áreas especificadas.")
    return(invisible(0L))
  }
  
  # Buscar archivos que correspondan a esos traits
  message("Buscando archivos...")
  
  archivos_a_copiar <- unique(unlist(lapply(traits_objetivo, function(trait) {
    grep(trait, archivos_GWAS, value = TRUE, ignore.case = TRUE, fixed = TRUE)
  })))
  
  #  Copiar archivos 
  if (length(archivos_a_copiar) > 0) {
    file.copy(from = archivos_a_copiar, to = carpeta_destino, overwrite = sobreescribir)
    message("✔️  Hecho. Traits: ",         length(traits_objetivo),
            " | Archivos copiados: ", length(archivos_a_copiar),
            " → ",                    carpeta_destino)
  } else {
    message("⚠️  No se encontró ningún archivo que coincida con los traits seleccionados.")
  }
  
  invisible(length(archivos_a_copiar))
}

# -------------------------------------------
# MICROGWAS ----

# Crear carpeta para guardar los clusters
carpeta_clusters <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
dir.create(carpeta_clusters, showWarnings = FALSE)

extraer_clusters (
  carpeta_entrada = "./Output/Piloto_Microbiota/MicroGWAS",
  carpeta_salida = carpeta_clusters
)

# ==== Barplot: genes / cluster ====

# --- 1. Leer todos los CSVs y construir dataframe ---
archivos_clusters <- list.files(carpeta_clusters, pattern = "\\.csv$", full.names = TRUE)

datos_barplot <- lapply(archivos_clusters, function(f) {
  df <- read.csv2(f)
  nombre <- gsub("\\.csv$", "", basename(f))
  data.frame(
    cluster      = nombre,
    n_total      = nrow(df),
    n_semilla    = sum(df$gene.senal_inicial == "Si", na.rm = TRUE),
    n_no_semilla = sum(df$gene.senal_inicial == "No", na.rm = TRUE)
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(
    cluster = dplyr::recode(cluster, !!!dic_clusters_MicroGWAS)  
  )

# --- 2. Añadir info de trait ---
datos_barplot <- datos_barplot %>%
  dplyr::mutate(
    trait_id      = gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", cluster),
    cluster_corto = gsub("^EFO_[0-9]+_", "", cluster)
  ) %>%
  dplyr::left_join(
    traits_MicroGWAS_areas %>% dplyr::select(Rasgo, name) %>% dplyr::distinct(),
    by = c("trait_id" = "Rasgo")
  ) %>%
  dplyr::mutate(
    etiqueta = paste0(acortar_nombres_microbioma(name), "\n(", cluster_corto, ")")
  )

# --- 3. Pasar a formato largo para barras apiladas ---
datos_largo <- datos_barplot %>%
  dplyr::select(etiqueta, trait_id, n_semilla, n_no_semilla) %>%
  tidyr::pivot_longer(
    cols      = c(n_semilla, n_no_semilla),
    names_to  = "tipo",
    values_to = "n_genes"
  ) %>%
  dplyr::mutate(
    tipo = ifelse(tipo == "n_semilla", "Seed genes", "Non-seed genes")
  )

# --- 4. Orden: agrupar por trait y dentro por tamaño descendente ---
orden_clusters <- datos_barplot %>%
  dplyr::arrange(trait_id, desc(n_total)) %>%
  dplyr::pull(etiqueta)

datos_largo$etiqueta <- factor(datos_largo$etiqueta, levels = rev(orden_clusters))

# --- 5. Colores por trait para los genes NO semilla
#         y versión más oscura para los genes semilla ---
colores_6_traits <- c(
  "EFO_0007753" = "#FFF2AE",
  "EFO_0007874" = "pink",
  "EFO_0007883" = "#BFEFFF",
  "EFO_0011013" = "#FDBF6F",
  "EFO_0801228" = "#D2B4DE",
  "EFO_0801229" = "#C1E1C1"
)

# Colores oscuros para genes semilla (mismo tono pero más saturado)
colores_semilla <- c(
  "EFO_0007753" = "#E5C494",
  "EFO_0007874" = "magenta4",
  "EFO_0007883" = "skyblue4",
  "EFO_0011013" = "tan4",
  "EFO_0801228" = "#6A3D9A",
  "EFO_0801229" = "#2F4F4F"
)

efo_a_nombre <- traits_MicroGWAS_areas %>%
  dplyr::filter(Rasgo %in% names(colores_6_traits)) %>%    
  dplyr::select(Rasgo, name) %>%
  dplyr::distinct() %>%
  tibble::deframe()                     
colores_globales <- setNames(
  colores_6_traits[names(efo_a_nombre)],  
  efo_a_nombre                           
)
# Construir vector de colores para scale_fill
# cada etiqueta tiene dos colores (semilla y no semilla)
colores_fill <- datos_largo %>%
  dplyr::mutate(
    color = ifelse(tipo == "Seed genes",
                   colores_semilla[trait_id],
                   colores_6_traits[trait_id]),
    fill_key = paste0(etiqueta, "_", tipo)
  ) %>%
  dplyr::select(fill_key, color) %>%
  dplyr::distinct()

colores_vector <- setNames(colores_fill$color, colores_fill$fill_key)

datos_largo <- datos_largo %>%
  dplyr::mutate(fill_key = paste0(etiqueta, "_", tipo))

# --- 6. Barplot ---
n_clusters <- n_distinct(datos_largo$etiqueta)

p <- ggplot(datos_largo, aes(y = etiqueta, x = n_genes, fill = fill_key)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(
    data = datos_barplot,
    aes(y = etiqueta, x = n_total, label = paste0(n_total, " (", n_semilla, ")")),
    hjust   = -0.15,
    size    = 6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = colores_vector,
    guide  = guide_legend(
      override.aes = list(fill = c("#4D4D4D", "#AEC6CF")),  # ejemplo: oscuro=semilla, claro=no semilla
      title = ""
    )) + 
  scale_x_continuous(
    position = "top",
    limits   = c(0, max(datos_barplot$n_total) * 1.2),
    expand   = c(0.02, 0)
  ) +
  labs(
    title = "Number of genes per cluster",
    x     = "Number of genes",
    y     = "",
    fill  = ""
  ) +
  theme_minimal() +
  theme(
    legend.position  = "none",    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold",
                                    margin = margin(t = 10, b = 10), size = 16),
    axis.line.x.top  = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x.top = element_line(color = "black"),
    axis.text.x.top  = element_text(size = 14,margin = margin(b = 5)),
    axis.title.x.top = element_text(margin = margin(b = 10, t = 10), size = 16),
    axis.text.y      = element_text(size = 11),
    plot.margin      = margin(t = 10, r = 80, b = 10, l = 10)
  )

# --- 7. Guardar ---
pdf("./Output/Gráficos/MicroGWAS/Barplot_genes_por_cluster.pdf",
    width  = 12,
    height = n_clusters * 0.4 + 0.5)

print(p)

dev.off()


# -------------------------------------------

# -------------------------------------------
# TRAITS NERVOUS / PSYCHIATRIC ----

# --- 1- Extraemos matrices de propagacion

copiar_traits_por_area(
  carpeta_origen = "./Output/Piloto_Microbiota/MicroGWAS_GWAS",
  carpeta_destino = "./Output/Piloto_Microbiota/Traits_Nervous_Psychiatric",
  tabla_traits = traits_MicroGWAS_areas,
  areas_filtro    = c("EFO_0000319", "MONDO_0004995")
)

# --- 2- Extraemos clusters

extraer_clusters(
  carpeta_entrada = "./Output/Piloto_Microbiota/Traits_Nervous_Psychiatric",
  carpeta_salida = "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric"
)
# -------------------------------------------

# -------------------------------------------
# TRAITS IMMUNE SYSTEM ----

# --- 1- Extraemos matrices de propagacion

copiar_traits_por_area(
  carpeta_origen = "./Output/Piloto_Microbiota/MicroGWAS_GWAS",
  carpeta_destino = "./Output/Piloto_Microbiota/Traits_Immune",
  tabla_traits = traits_MicroGWAS_areas,
  areas_filtro    = "EFO_0000540"
)

# --- 2- Extraemos clusters

extraer_clusters(
  carpeta_entrada = "./Output/Piloto_Microbiota/Traits_Immune",
  carpeta_salida = "./Output/Piloto_Microbiota/Clusters_Immune"
)
# -------------------------------------------

# -------------------------------------------
# RESTO DE TRAITS VAR COMUN ----

# --- 1- Extraemos matrices de propagacion

carpeta_origen  <- "./Output/Piloto_Microbiota/MicroGWAS_GWAS"
carpeta_destino <- "./Output/Piloto_Microbiota/Traits_VarComun"

archivos_GWAS <- list.files(carpeta_origen, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE)
dir.create(carpeta_destino, showWarnings = FALSE)

ID_MicroGWAS <- c(
  "EFO_0007753", 
  "EFO_0007874", 
  "EFO_0007883", 
  "EFO_0011013", 
  "EFO_0801228", 
  "EFO_0801229"
)
traits_objetivo <- traits_MicroGWAS_areas %>%   # tenemos 1628 traits
  filter(
    !grepl("EFO_0000618", therapeuticAreas, ignore.case = TRUE) &
      !grepl("MONDO_0002025", therapeuticAreas, ignore.case = TRUE) &
      !grepl("EFO_0000540", therapeuticAreas, ignore.case = TRUE) &
      !(Rasgo %in% ID_MicroGWAS)
  ) %>%
  pull(Rasgo) %>%
  trimws() %>%
  unique()   # <-- evita duplicados


message("Traits encontrados: ", length(traits_objetivo))  # tenemos los traits que queremos

archivos_a_copiar <- character(0)

message("Buscando archivos...")

# # Buscar los archivos que correspondan a esos traits
for (trait in traits_objetivo) {  
  coincidencias <- grep(
    trait, archivos_GWAS,
    value      = TRUE,
    ignore.case = TRUE,
  )
  archivos_a_copiar <- c(archivos_a_copiar, coincidencias)
}

## Copiar achivos en la carpeta final
if (length(archivos_a_copiar) > 0) {
  file.copy(from = archivos_a_copiar, to = carpeta_destino, overwrite = FALSE)
  message("Hecho. Traits: ", length(traits_objetivo),
          " | Archivos copiados: ", length(archivos_a_copiar),
          " → ", carpeta_destino)
} else {
  message("No se encontró ningún archivo que coincida con los traits seleccionados.")
}

# --- 2- Extraemos clusters

extraer_clusters(
  carpeta_entrada = "./Output/Piloto_Microbiota/Traits_VarComun",
  carpeta_salida =  "./Output/Piloto_Microbiota/Clusters_Traits_VarComun"
)

# -------------------------------------------

