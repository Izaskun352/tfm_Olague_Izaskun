
# Cargar scripts necesario

source("Scripts/04_Estudio_Clusters/B_Enrich_GOBP/0_Funciones_GOBP.R")
source("Scripts/Renombrar_Clusters.R")

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")
interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
universo_genes <- unique(na.omit(interactoma[,1]))
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
lista_intersecciones_im <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")

# OUTPUTS

gobp_completo_np <-read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_NP_MicroGWAS/GO_Resumen_Completo.csv")
gobp_simplificado_np <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_NP_MicroGWAS/GO_Resumen_Simplificado.csv")

gobp_completo_im <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_Immune_MicroGWAS/GO_Resumen_Completo.csv")
gobp_simplificado_im <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_Immune_MicroGWAS/GO_Resumen_Interseccion_Immune_MicroGWAS.csv")

# -------------------------------------------------------------------
# GOBP DE LAS INTERSECCIONES
# -------------------------------------------------------------------
# Traits NP + MicroGWAS

analizar_GO_intersecciones(
  lista_intersecciones = lista_intersecciones_np,
  carpeta_salida       = "./Output/Piloto_Microbiota/GO_Intersecciones_MicroGWAS_NervousPsychiatric",
  universo_genes       = universo_genes
)

gobp_completo_np <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_MicroGWAS_NervousPsychiatric/GO_Resumen_Completo.csv")
gobp_simplificado_np <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_MicroGWAS_NervousPsychiatric/GO_Pleiotropia_Resumen.csv")

# Traits IM + MicroGWAS

analizar_GO_intersecciones(
  lista_intersecciones = lista_intersecciones_im,
  carpeta_salida       = "./Output/Piloto_Microbiota/GO_Intersecciones_Immune_MicroGWAS",
  universo_genes       = universo_genes
)

gobp_completo_im <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_Immune_MicroGWAS/GO_Resumen_Completo.csv")
gobp_simplificado_im <- read.csv2("./Output/Piloto_Microbiota/GO_Intersecciones_Immune_MicroGWAS/GO_Resumen_Interseccion_Immune_MicroGWAS.csv")


# -------------------------------------------------------------------
# UNIMOS RESULTADOS
# -------------------------------------------------------------------

tabla_GO_intersecciones <- bind_rows(
  gobp_simplificado_im %>% mutate(Fuente = "IM"),
  gobp_simplificado_np %>% mutate(Fuente = "NP")
)
tabla_GO_intersecciones <- bind_rows(
  estandarizar_tabla_IM(gobp_simplificado_im) %>% mutate(Fuente = "IM"),
  estandarizar_tabla_NP(gobp_simplificado_np) %>% mutate(Fuente = "NP")
)

# -------------------------------------------------------------------
# REPRESENTACION GRAFICA
# -------------------------------------------------------------------
# --- Top N términos GO por par ---
n_top <- 5
top_terms_por_par <- tabla_GO_intersecciones %>%
  dplyr::group_by(Cluster_MicroGWAS, Cluster_NP_IM) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = n_top) %>%
  dplyr::ungroup() %>%
  dplyr::pull(Description) %>%
  unique()

# --- Construir matriz y anotación ---
tabla_procesada <- tabla_GO_intersecciones %>%
  dplyr::mutate(
    Trait_1_Corto       = acortar_nombres_microbioma(Trait_Micro_Nombre),
    Trait_1_Corto       = stringr::str_to_sentence(Trait_1_Corto),
    Trait_NP_IM_Nombre  = stringr::str_to_sentence(Trait_NP_IM_Nombre),
    Cluster_Micro_label = gsub(".*_Cluster_", "Cluster_", cluster_names[Cluster_MicroGWAS]),
    Cluster_NP_IM_label = gsub(".*_Cluster_", "Cluster_", cluster_names[Cluster_NP_IM]),
    Par = paste0(Cluster_Micro_label, " vs ", Cluster_NP_IM_label, " (n = ", N_genes_interseccion, ")\n", Trait_NP_IM_Nombre),
    log_pval = -log10(p.adjust) 
  )

# -------------------------------------------------------------------
# HEATMAP ----

# Matriz
matriz_go <- tabla_procesada %>%
  dplyr::filter(Description %in% top_terms_por_par) %>%
  dplyr::select(Description, Par, log_pval) %>%
  tidyr::pivot_wider(names_from  = Par,
                     values_from = log_pval,
                     values_fill = 0) %>%
  tibble::column_to_rownames("Description")


# --- 4. Anotación de columnas ---
traits_por_par <- tabla_procesada %>%
  dplyr::select(Par, Trait_1_Corto, Trait_NP_IM_Nombre, Fuente) %>%
  dplyr::distinct()

# Subsets de colores
niveles_micro  <- sort(unique(traits_por_par$Trait_1_Corto))
niveles_np_im  <- sort(unique(traits_por_par$Trait_NP_IM_Nombre))

names(colores_micro_por_nombre) <- stringr::str_to_sentence(names(colores_micro_por_nombre))
names(colores_traits_NP_IM)     <- stringr::str_to_sentence(names(colores_traits_NP_IM))

colores_micro_c  <- colores_micro_por_nombre[niveles_micro]
colores_np_im_c  <- colores_traits_NP_IM[niveles_np_im]  

annotation_col <- data.frame(
  Microbioma = factor(traits_por_par$Trait_1_Corto,    levels = names(colores_micro_c)),
  Enfermedad = factor(traits_por_par$Trait_NP_IM_Nombre, levels = names(colores_np_im_c)),
  #Fuente     = factor(traits_por_par$Fuente,           levels = c("IM", "NP")),
  row.names  = traits_por_par$Par
)

annotation_colors_completa <- list(
  Microbioma = colores_micro_c,
  Enfermedad = colores_np_im_c,
  Fuente     = c("IM" = "#E8A838", "NP" = "#6A5ACD")
)

# Plot

paleta <- colorRampPalette(c("#EEF9C4", "#7FCDBB", "#2C7FB8", "#1A4E88"))(100)

pdf("./Output/Gráficos/MicroGWAS/Heatmap_Enrich_GOBP_IM_NP.pdf", height = 30, width = 40)
print(pheatmap(
  mat               = as.matrix(matriz_go),
  color             = paleta,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors_completa,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  cellwidth         = 45,
  cellheight        = 15,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 11,
  fontsize_col      = 9,
  angle_col         = 45,
  main              = "GO-BP enrichment in pleiotropic gene intersections (MicroGWAS vs NP & IM)",
  border_color      = "white"
))
dev.off()

# -------------------------------------------------------------------

# -------------------------------------------------------------------
# DOTPLOT ----

datos_dotplot <- tabla_procesada %>%
  # Mantener solo los términos top seleccionados para el heatmap
  dplyr::filter(Description %in% top_terms_por_par) %>%
  dplyr::mutate(
    # Parseo del GeneRatio (idéntico a la función)
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatio_den = as.numeric(sub(".*/", "", GeneRatio)),
    GeneRatio_val = GeneRatio_num / GeneRatio_den,
    
    # Eje X: Nombre del rasgo NP/IM + los dos sub-clusters debajo
    Cluster_Pair_label = paste0(Cluster_Micro_label, " vs ", Cluster_NP_IM_label, " (n = ", N_genes_interseccion, ")\n", Trait_NP_IM_Nombre)
  )

# --- 2. Ordenar los términos GO de forma decreciente por p-valor ---
datos_dotplot$Description <- factor(
  datos_dotplot$Description,
  levels = unique(datos_dotplot$Description[order(datos_dotplot$p.adjust, decreasing = TRUE)])
)

# --- 3. Construir el Gráfico con la estructura exacta de tu función ---
dotplot <- ggplot(datos_dotplot, aes(x = Cluster_Pair_label, y = Description)) +
  # Puntos mapeados por tamaño (GeneRatio) y color (p.adjust)
  geom_point(aes(size = GeneRatio_val, color = p.adjust)) +
  
  # Escala de colores "magma" con transformación logarítmica
  scale_color_viridis_c(
    option    = "magma",
    direction = 1,
    trans     = "log10",
    name      = "p.adjust"
  ) +
  
  # Rango de tamaño de los puntos de la función c(6, 12)
  scale_size_continuous(name = "GeneRatio", range = c(6, 12)) +
  
  # CORRECCIÓN EN PANEL: Facetado únicamente por el Trait MicroGWAS
  facet_grid(~ Trait_1_Corto, scales = "free_x", space = "free_x") +
  
  # Estilo visual de tu función original
  theme_bw() +
  theme(
    # Ángulo 45° con ajuste horizontal como pedía tu función para leer los rasgos de abajo
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 18), 
    axis.text.y        = element_text(size = 22),
    axis.title         = element_blank(),
    plot.title         = element_text(hjust = 0.5, size = 20, face = "bold"),
    strip.background   = element_blank(), 
    strip.text         = element_text(size = 14, face = "bold", color = "black"), 
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.spacing      = unit(0.05, "lines"), # El espaciado ultra-estrecho de tu función original
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  labs(title = "GO-BP Enrichment in Gene Intersections")

# --- Guardar en PDF ---
pdf("./Output/Gráficos/MicroGWAS/Dotplot_Enrich_GOBP_IM_NP.pdf", width = 40, height = 18)
print(dotplot)
dev.off()
# -------------------------------------------------------------------