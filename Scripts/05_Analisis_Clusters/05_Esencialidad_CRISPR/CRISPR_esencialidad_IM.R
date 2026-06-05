
source("scripts/00_setup.R")  ## Abrimos script con las librerias 

library("depmap")
library("ExperimentHub")


# ESENCIALIDAD GENÉTICA CON CRISPR KNOCK-OUT -- GENES INTERSECCIÓN MICROGWAS + NERVOUS / PSYCHIATRIC
# -------------------------------------------------------------------------------------------------

# INPUTS
#
#


lista_intersecciones <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds") # ----


# -------------------------------------------------------------------------------------------------
# Comprueba los nombres exactos de tus columnas
colnames(crispr_data)
colnames(metadata)

# =================================================================================================
# BLOQUE 1: DESCARGAR depmap::crispr + ver líneas celulares SNC / SNP  ====

# ---- ExperimentHub descarga los datos bajo demanda y los cachea localmente 
## - R crea carpeta en disco duro (la caché)  y al descargar dataset de DepMap se guarda ahí
eh <- ExperimentHub()
query(eh, "depmap")  # Ver los datasets disponibles de depmap
#crispr <- eh[["EH3081"]]


# ---- Cargar el dataset CRISPR (Chronos scores)
crispr_data <- depmap::depmap_crispr()  # descarga el release más reciente disponible - con esto no haria falta lo anterior
glimpse(crispr_data)
# Estructura del objeto:
# - gene_name  : símbolo del gen
# - entrez_id  : Entrez ID
# - depmap_id  : identificador de línea celular
# - cell_line  : nombre de la línea celular
# - lineage    : tejido de origen
# - dependency : Chronos score (< 0 = dependencia/esencialidad)


# ---- Estudiar líneas celulares relacionadas con el sistema nervioso

# 1. Cargar la metadata (si no lo tienes en memoria)
metadata <- depmap_metadata()

# 2. Ver TODOS los linajes principales (tejidos generales)
linajes_unicos <- unique(metadata$lineage)
print(linajes_unicos)

# 3. Contar cuántas líneas celulares hay por tejido
conteo_linajes <- table(metadata$lineage)
conteo_linajes <- sort(conteo_linajes, decreasing = TRUE)
print(conteo_linajes)


# =================================================================================================


# =================================================================================================
# BLOQUE 2: DEPMAP-CRISPR ====

# Análisis por par
carpeta_depmap <- "./Output/Piloto_Microbiota/DepMap_im_MicroGWAS"
dir.create(carpeta_depmap, showWarnings = FALSE)

lista_resultados_depmap <- list()
lista_matrices_celulas        <- list()
lista_matrices_tejidos        <- list()

for (nombre_par in names(lista_intersecciones)) {
  
  par <- lista_intersecciones[[nombre_par]]
  message("=== DepMap: ", nombre_par, " (Jaccard = ", par$jaccard, ") ===")
  
  # --------- DATOS CRUDOS - una fila por gen x línea celular
  
  crispr_par <- crispr_data %>%
    dplyr::filter(gene_name %in% par$genes_simbolo)
  
  if (nrow(crispr_par) == 0) {
    message("  ⚠️ Ningún gen encontrado en DepMap")
    next
  }
  
  # Añadir info de linaje
  crispr_par_anotado <- crispr_par %>%
    dplyr::left_join(
      metadata %>% dplyr::select(depmap_id, lineage,
                                 cell_line_meta = cell_line),
      by = "depmap_id"
    ) %>%
    # Añadir info de cluster y semilla
    dplyr::left_join(
      par$genes_tabla,
      by = c("gene_name" = "gene.gene")
    )
  
  # ---------- MATRIZ A: UNA COLUMNA POR LÍNEA CELULAR
  
  matriz_depmap_celulas <- crispr_par_anotado %>%
    dplyr::mutate(col_name = paste0(lineage, "__", cell_line_meta)) %>%
    dplyr::select(gene_name, nombre_par, Cluster_1, Cluster_2,
                  Trait_Nombre_1, Trait_Nombre_2,
                  es_semilla_en_c1, es_semilla_en_c2, es_semilla_en_alguno,
                  col_name, dependency) %>%
    tidyr::pivot_wider(
      names_from  = col_name,
      values_from = dependency,
      values_fn   = mean
    ) %>%
    dplyr::arrange(gene_name)
  
  # ---------- MATRIZ B: UNA COLUMNA POR LINAJE (tejido)
  
  matriz_depmap_tejidos <- crispr_par_anotado %>%
    dplyr::group_by(gene_name, nombre_par, Cluster_1, Cluster_2,
                    Trait_Nombre_1, Trait_Nombre_2,
                    es_semilla_en_c1, es_semilla_en_c2,
                    es_semilla_en_alguno, lineage) %>%
    dplyr::summarise(
      mean_dependency = mean(dependency, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from  = lineage,
      values_from = mean_dependency
    ) %>%
    dplyr::arrange(gene_name)
  
  # ----------- ESENCIALIDAD
  
  # Agregar scores por gen y unir con info de semilla
  esencialidad <- crispr_par %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarise(
      mean_dependency       = mean(dependency, na.rm = TRUE),
      median_dependency     = median(dependency, na.rm = TRUE),
      sd_dependency         = sd(dependency, na.rm = TRUE),
      n_cell_lines          = n(),
      prop_essential        = mean(dependency < -0.5, na.rm = TRUE),
      prop_highly_essential = mean(dependency <= -1,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      esencialidad = dplyr::case_when(
        mean_dependency < -1   ~ "Esencial fuerte",
        mean_dependency < -0.5 ~ "Esencial moderado",
        mean_dependency < 0    ~ "Efecto leve",
        TRUE                   ~ "No esencial"
      )
    ) %>%
    dplyr::left_join(
      par$genes_tabla,
      by = c("gene_name" = "gene.gene")
    ) %>%
    dplyr::arrange(mean_dependency)
  
  # ------------- GUARDAR CSV 
  
  write.csv2(matriz_depmap_celulas,
             file = file.path(carpeta_depmap,
                              paste0("DepMap_Matriz_LineaCelular_", nombre_par, ".csv")),
             row.names = FALSE)
  
  write.csv2(matriz_depmap_tejidos,
             file = file.path(carpeta_depmap,
                              paste0("DepMap_Matriz_Linaje_", nombre_par, ".csv")),
             row.names = FALSE)
  
  # Acumular en listas
  lista_resultados_depmap[[nombre_par]] <- esencialidad
  lista_matrices_celulas[[nombre_par]]        <- matriz_depmap_celulas
  lista_matrices_tejidos[[nombre_par]]        <- matriz_depmap_tejidos
}

# ===== TABLAS MAESTRAS ====

if (length(lista_resultados_depmap) > 0) {
  
  # Guardamos la esencialidad- si es esencial como media de todos los tejidos
  
  tabla_maestra_depmap <- dplyr::bind_rows(lista_resultados_depmap)
  
  write.csv2(tabla_maestra_depmap,
             file = file.path(carpeta_depmap, "DepMap_Esencialidad_Resumen.csv"),
             row.names = FALSE)
  
  saveRDS(tabla_maestra_depmap,
          file = file.path(carpeta_depmap, "DepMap_Esencialidad_Resumen.rds"))
  
  # MATRIZ LINEAS CELULARES
  
  matriz_celulas_maestra <- dplyr::bind_rows(lista_matrices_celulas) %>%
    dplyr::arrange(nombre_par, gene_name)
  
  write.csv2(matriz_celulas_maestra,
             file = file.path(carpeta_depmap, "DepMap_Matriz_Celulas_Maestra.csv"),
             row.names = FALSE)
  
  saveRDS(matriz_celulas_maestra,
          file = file.path(carpeta_depmap, "DepMap_Matriz_Celulas_Maestra.rds"))
  
  # MATRIZ LINAJES
  
  matriz_tejidos_maestra <- dplyr::bind_rows(lista_matrices_tejidos) %>%
    dplyr::arrange(nombre_par, gene_name)
  
  write.csv2(matriz_tejidos_maestra,
             file = file.path(carpeta_depmap, "DepMap_Matriz_Tejidos_Maestra.csv"),
             row.names = FALSE)
  
  saveRDS(matriz_tejidos_maestra,
          file = file.path(carpeta_depmap, "DepMap_Matriz_Tejidos_Maestra.rds"))
  
  message("Pares procesados: ",    length(lista_resultados_depmap))
  message("Total genes (esenc.): ", nrow(tabla_maestra_depmap))
  message("Total genes (mat. A): ", nrow(matriz_celulas_maestra))
  message("Total genes (mat. B): ", nrow(matriz_tejidos_maestra))
}






# ===== REPRESENTACIÖN GRÄFICA ---- PENDIENTE HEATMAP 
colores_6_traits <- c(
  "EFO_0007753" = "#FFFF99",
  "EFO_0007874" = "#FFDAC1",
  "EFO_0007883" = "#AEC6CF",
  "EFO_0011013" = "#B5EAD7",
  "EFO_0801228" = "#C9B1FF",
  "EFO_0801229" = "#FFD1DC"
)
# ===== BOXPLOT ====

carpeta_boxplots_depmap <- "./Output/Gráficos/MicroGWAS/Boxplots_Esencialidad_DepMap"
dir.create(carpeta_boxplots_depmap, showWarnings = FALSE)

message("Cargando datos...")

matriz_tejidos <- readRDS("./Output/Piloto_Microbiota/DepMap_im_MicroGWAS/DepMap_Matriz_Tejidos_Maestra.rds")

# 3. Transformar la Matriz (De ancha a larga) para ggplot
datos_plot <- matriz_tejidos %>%
  tidyr::pivot_longer(
    # Excluimos las columnas de metadatos, así todo el resto (los tejidos) se agrupan
    cols = -c(gene_name, nombre_par, Cluster_1, Cluster_2, 
              Trait_Nombre_1, Trait_Nombre_2, 
              es_semilla_en_c1, es_semilla_en_c2, es_semilla_en_alguno),
    names_to = "lineage",
    values_to = "dependency_score"
  ) %>%
  # Quitar los NA (genes que no se midieron en un tejido concreto)
  dplyr::filter(!is.na(dependency_score)) %>%
  # Crear las etiquetas y extraer el ID de MicroGWAS para los colores
  dplyr::mutate(
    trait_micro_id = ifelse(grepl("EFO", Cluster_1), 
                            sub("_Cluster.*", "", Cluster_1), 
                            sub("_Cluster.*", "", Cluster_2)),
    etiqueta = paste0(Trait_Nombre_1, " vs \n", Trait_Nombre_2)
  )

# 4. Generar los PDFs
traits_micro_unicos <- unique(datos_plot$trait_micro_id)
linajes_unicos <- unique(datos_plot$lineage)

for (trait_id in traits_micro_unicos) {
  
  message("Generando PDF para: ", trait_id)
  
  # Filtrar datos solo para este Trait
  datos_trait <- datos_plot %>%
    dplyr::filter(trait_micro_id == trait_id)
  
  # Obtener las etiquetas únicas para este trait
  etiquetas_trait <- unique(datos_trait$etiqueta)
  
  # Asignar los colores (todos del mismo color base del EFO)
  colores_trait <- setNames(rep(colores_6_traits[trait_id], length(etiquetas_trait)), etiquetas_trait)
  
  # Crear el PDF
  pdf(file.path(carpeta_boxplots_depmap, paste0("DepMap_Esencialidad_", trait_id, ".pdf")),
      width = 10, height = 8)
  
  # Un gráfico por cada página (Linaje)
  for (lin in linajes_unicos) {
    
    datos_linaje <- datos_trait %>% dplyr::filter(lineage == lin)
    
    # Evitar gráficos vacíos
    if(nrow(datos_linaje) == 0) next
    
    p <- ggplot(datos_linaje, aes(y = etiqueta, x = dependency_score, fill = etiqueta)) +
      geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
      # La línea roja absoluta que marca la muerte celular (-1)
      geom_vline(xintercept = -1, linetype = "dashed", color = "red", size = 1) +
      scale_fill_manual(values = colores_trait) +
      theme_classic() +
      theme(
        axis.text.y     = element_text(size = 10),
        axis.title.y    = element_blank(),
        legend.position = "none",
        plot.title      = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey50"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
      ) +
      labs(
        title    = toupper(gsub("_", " ", lin)),
        subtitle = "Chronos Dependency Score (Línea roja = Esencialidad fuerte)",
        x        = "Score de Dependencia"
      )
    
    print(p)
  }
  
  dev.off()
  message("  ✔️ PDF completado.")
}

message("¡Proceso finalizado! PDFs guardados en la carpeta de salida.")

# ===== PDF GLOBAL CON TODAS LAS INTERSECCIONES ===

message("Generando PDF global con todas las intersecciones...")

# 1. Orden fijo: agrupar las etiquetas por el ID de MicroGWAS para que salgan juntas
todas_etiquetas <- datos_plot %>%
  dplyr::select(etiqueta, trait_micro_id) %>%
  dplyr::distinct() %>%
  dplyr::arrange(trait_micro_id) %>% 
  dplyr::pull(etiqueta)

# Aplicar el orden al factor (usamos rev() para que el primero de la lista salga arriba en el gráfico)
datos_plot$etiqueta <- factor(datos_plot$etiqueta, levels = rev(todas_etiquetas))

# 2. Preparar el diccionario de colores globales
colores_global <- datos_plot %>%
  dplyr::select(etiqueta, trait_micro_id) %>%
  dplyr::distinct() %>%
  # Asignar el color base a cada etiqueta según su EFO
  dplyr::mutate(color = colores_6_traits[trait_micro_id]) %>%
  dplyr::select(etiqueta, color) %>%
  tibble::deframe() # Convierte el dataframe en un vector con nombres

# 3. Crear el PDF gigante
# Usamos un tamaño grande (height = 15) para que los nombres de las 25 intersecciones no se solapen
pdf(file.path(carpeta_boxplots_depmap, "DepMap_Esencialidad_TODAS_intersecciones_global.pdf"),
    width = 15, height = 15)

for (lin in linajes_unicos) {
  
  datos_linaje <- datos_plot %>% dplyr::filter(lineage == lin)
  
  if(nrow(datos_linaje) == 0) next
  
  p <- ggplot(datos_linaje, aes(y = etiqueta, x = dependency_score, fill = etiqueta)) +
    geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
    geom_vline(xintercept = -1, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = colores_global) +
    theme_classic() + 
    theme(
      # Hacemos la letra del eje Y un poco más grande para que se lea bien en el PDF global
      axis.text.y      = element_text(size = 12),
      axis.title.y     = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5, size = 14, color = "grey50"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title    = paste0("DISTRIBUCIÓN GLOBAL — ", toupper(gsub("_", " ", lin))),
      subtitle = "MicroGWAS vs Nervous/Psychiatric intersections (Red line = Strong Essentiality)",
      x        = "Chronos Dependency Score"
    )
  
  print(p)
}

dev.off()

message("¡PDF global generado con éxito!")


# =================================================================================================

# =================================================================================================
# BLOQUE 3: CRISPRBrain

# Descargamos dataset
## Web oficial descargamos archivo: 'Glutameric Neuron-Survival- CRISPR'

crisprbrain_data <- read.csv("./Data/CRIPRbrain_iCRISPR_Survival_GlutamericNeurons.csv", stringsAsFactors = FALSE)
head(crisprbrain_data)

carpeta_crisprbrain <- "./Output/Piloto_Microbiota/CRISPRbrain_Esencialidad_im_MicroGWAS"
dir.create(carpeta_crisprbrain, showWarnings = FALSE)

lista_resultados_crisprbrain <- list()

for (nombre_par in names(lista_intersecciones)) {
  
  par <- lista_intersecciones[[nombre_par]]
  message("=== CRISPRbrain: ", nombre_par, " ===")
  # 1. Filtramos usando la columna TSS en lugar de Gene
  crispr_par <- crisprbrain_data %>%
    dplyr::filter(TSS %in% par$genes_simbolo)
  
  if (nrow(crispr_par) == 0) {
    message("  ⚠️ Ningún gen encontrado en CRISPRbrain")
    next
  }
  
  # 2. Preparar la tabla de resultados adaptada a tus columnas
  esencialidad_cb <- crispr_par %>%
    # Seleccionamos las columnas
    dplyr::select(gene_name = TSS, Phenotype, P.Value, Gene.Score, Hit.Class) %>%
    dplyr::mutate(
      # Usamos la clasificación oficial de CRISPRbrain
      esencialidad = dplyr::case_when(
        Hit.Class == "Negative Hit" ~ "Esencial (muerte celular)",
        Hit.Class == "Positive Hit" ~ "Ventaja de crecimiento",
        TRUE ~ "No esencial / Ruido"
      )
    ) %>%
    # Unimos usando "Gene_Symbol" (que era TSS)
    dplyr::left_join(
      par$genes_tabla %>%
        dplyr::distinct(gene.gene, es_semilla_en_c1 , es_semilla_en_c2,es_semilla_en_alguno),
      by = c("gene_name" = "gene.gene")
    ) %>%
    dplyr::mutate(
      nombre_par     = nombre_par,
      Cluster_1      = par$c1,
      Cluster_2      = par$c2,
      Trait_Nombre_1 = par$trait_1,
      Trait_Nombre_2 = par$trait_2
    ) %>%
    dplyr::arrange(Phenotype) # Ordenamos por el score (los más negativos arriba)
  
  # Imprimir recuento por consola
  esenciales_count <- sum(esencialidad_cb$esencialidad == "Esencial (Muerte neuronal)", na.rm = TRUE)
  message("  ✔️ Genes letales en neuronas sanas: ", esenciales_count)
  
  # Exportar el CSV de este par concreto
  write.csv2(esencialidad_cb,
             file = file.path(carpeta_crisprbrain, paste0("CRISPRbrain_", nombre_par, ".csv")),
             row.names = FALSE)
  
  # Guardar en la lista
  lista_resultados_crisprbrain[[nombre_par]] <- esencialidad_cb
}


# -------- TABLA MAESTRA

if (length(lista_resultados_crisprbrain) > 0) {
  
  # Combinar todas las intersecciones en una sola tabla
  tabla_maestra_cb <- dplyr::bind_rows(lista_resultados_crisprbrain)
  
  # Exportar CSV
  write.csv2(tabla_maestra_cb,
             file = file.path(carpeta_crisprbrain, "CRISPRbrain_Esencialidad_Resumen.csv"),
             row.names = FALSE)
  
  # Exportar RDS
  saveRDS(tabla_maestra_cb,
          file = file.path(carpeta_crisprbrain, "CRISPRbrain_Esencialidad_Resumen.rds"))
  
  message("\n¡Análisis completado!")
  message("Total pares procesados: ", length(lista_resultados_crisprbrain))
  message("Total genes evaluados en neuronas: ", nrow(tabla_maestra_cb))
} else {
  message("\nNo se encontraron coincidencias en ningún par.")
}


# ==== REPRESENTACIÖN GRAFICA ====
# ==== BOXPLOT ====

tabla_maestra_cb <- readRDS("./Output/Piloto_Microbiota/CRISPRbrain_Esencialidad_im_MicroGWAS/CRISPRbrain_Esencialidad_Resumen.rds")

carpeta_boxplots_cb <- "./Output/Gráficos/MicroGWAS/Boxplots_Esencialidad_CRISPRbrain_im_MicroGWAS"
dir.create(carpeta_boxplots_cb, showWarnings = FALSE)

# 3. Preparar los datos y crear etiquetas
datos_plot_cb <- tabla_maestra_cb %>%
  # Quitar por si acaso algún gen no tuviera score
  dplyr::filter(!is.na(Phenotype)) %>%
  dplyr::mutate(
    # Extraer el ID de MicroGWAS para que funcione el diccionario de colores
    trait_micro_id = ifelse(grepl("EFO", Cluster_1), 
                            sub("_Cluster.*", "", Cluster_1), 
                            sub("_Cluster.*", "", Cluster_2)),
    etiqueta = paste0(Trait_Nombre_1, " vs \n", Trait_Nombre_2)
  )

# ===== pdf individual
traits_micro_unicos <- unique(datos_plot_cb$trait_micro_id)

for (trait_id in traits_micro_unicos) {
  
  message("Generando PDF para: ", trait_id)
  
  # Filtrar los datos para este Trait específico
  datos_trait <- datos_plot_cb %>% dplyr::filter(trait_micro_id == trait_id)
  
  # Asignar los colores
  etiquetas_trait <- unique(datos_trait$etiqueta)
  colores_trait <- setNames(rep(colores_6_traits[trait_id], length(etiquetas_trait)), etiquetas_trait)
  
  # Crear el PDF
  pdf(file.path(carpeta_boxplots_cb, paste0("CRISPRbrain_Esencialidad_", trait_id, ".pdf")),
      width = 10, height = 8)
  
  p <- ggplot(datos_trait, aes(y = etiqueta, x = Phenotype, fill = etiqueta)) +
    geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
    # Línea en el 0: Hacia la izquierda es toxicidad neuronal
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 1) +
    scale_fill_manual(values = colores_trait) +
    theme_classic() +
    theme(
      axis.text.y      = element_text(size = 10),
      axis.title.y     = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey50"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title    = "CRISPRbrain Phenotype Score ",
      subtitle = "(Valores < 0 indican muerte neuronal)",
      x        = "Phenotype Score"
    )
  
  print(p)
  dev.off()
}

# ===== pdf global
message("Generando PDF global...")

# Ordenar los datos por ID para que los colores salgan agrupados en bloque
datos_plot_cb <- datos_plot_cb %>%
  dplyr::arrange(desc(trait_micro_id)) %>%
  dplyr::mutate(etiqueta = factor(etiqueta, levels = unique(etiqueta)))

# Extraer colores para el global
colores_global <- setNames(colores_6_traits[datos_plot_cb$trait_micro_id], datos_plot_cb$etiqueta)

# Crear el PDF global grande
pdf(file.path(carpeta_boxplots_cb, "CRISPRbrain_Esencialidad_TODAS_intersecciones.pdf"), 
    width = 15, height = 15)

p_global <- ggplot(datos_plot_cb, aes(y = etiqueta, x = Phenotype, fill = etiqueta)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 1) +
  scale_fill_manual(values = colores_global) +
  theme_classic() +
  theme(
    axis.text.y      = element_text(size = 12),
    axis.title.y     = element_blank(),
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 14, color = "grey50"),
  ) +
  labs(
    title    = "ESENCIALIDAD GENÉTICA - CRISPRbrain",
    subtitle = "All MicroGWAS vs Nervous/Psychiatric intersections (Values < 0 = Neuronal toxicity)",
    x        = "Phenotype Score"
  )

print(p_global)
dev.off()

message("¡Todos los PDFs generados con éxito!")











# =================================================================================================



