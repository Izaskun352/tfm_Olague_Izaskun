
# Librerias
source("scripts/00_setup.R")

# Cargar script necesario

source("Scripts/04_Estudio_Clusters/04_Enrich_InterPro/01_Funciones_Enrich_InterPro.R")
source("Scripts/Renombrar_Clusters.R")

# INPUTS
interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
universo_genes <- unique(na.omit(interactoma[,1]))
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

TERM2GENE_interpro <- readRDS("./Data/Diccionarios/Anotaciones_InterPro/TERM2GENE_InterPro.rds")
TERM2NAME_interpro <- readRDS("./Data/Diccionarios/Anotaciones_InterPro/TERM2NAME_InterPro.rds")

# OUTPUTS

carpeta_InterPro_clusters <- "./Output/Piloto_Microbiota/InterPro_Clusters_MicroGWAS"

# CALCULAR ENRIQUECIMIENTO INTERPRO DE LOS CLUSTERS MICROGWAS

carpeta_InterPro_clusters <- "./Output/Piloto_Microbiota/InterPro_Clusters_MicroGWAS"
dir.create(carpeta_InterPro_clusters, showWarnings = FALSE)

analizar_InterPro_carpeta(
  carpeta_clusters   = "./Output/Piloto_Microbiota/Clusters_MicroGWAS",
  carpeta_salida     = carpeta_InterPro_clusters,
  universo_genes     = universo_genes,
  TERM2GENE          = TERM2GENE_interpro,
  TERM2NAME          = TERM2NAME_interpro,
  tabla_traits       = traits_MicroGWAS_areas,
  n_top_simplificado = 30)

#resultados_interpro_microGWAS <- "./Output/Piloto_Microbiota/InterPro_Clusters_MicroGWAS/InterPro_Resumen_Simplificado.csv"

# --- REPRESENTAR GRAFICAMENTE

tabla_interpro <- read.csv2(
  file.path(carpeta_InterPro_clusters, "InterPro_Resumen_Simplificado.csv")) %>%
  dplyr::mutate(
    Cluster_Origen = dplyr::recode(Cluster_Origen, !!!dic_clusters_MicroGWAS)) %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)),
    GeneRatio_den = as.numeric(sub(".*/", "", GeneRatio)),
    GeneRatio_val = GeneRatio_num / GeneRatio_den,
    efo_id        = gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", Cluster_Origen),
    Trait_Nombre  = efo_a_nombre[efo_id],
    Cluster_corto = gsub("^EFO_[0-9]+_", "", Cluster_Origen),
    log_padj      = -log10(p.adjust + 1e-10)) # evitar log(0)

top2_por_cluster <- tabla_interpro %>%
  group_by(Cluster_Origen) %>%
  arrange(p.adjust) %>%
  slice_head(n = 2) %>%
  ungroup()

# IDs únicos resultantes (pueden ser más de 3*n_clusters por solapamiento)
ids_top2 <- unique(top2_por_cluster$ID)
message("Términos InterPro únicos tras top2 por cluster: ", length(ids_top2))

# Filtrar la tabla maestra a esos términos  

tabla_top2 <- tabla_interpro %>%
  dplyr::filter(ID %in% ids_top2) %>%
  dplyr::mutate(
    Termino_label = substr(Description, 1, 50), 
    Trait_Nombre  = efo_a_nombre_full[efo_id],   # log_padj_raw   = -log10(p.adjust),
                                                # log_padj       = pmin(log_padj_raw, 8)
    Trait_Corto   = stringr::str_to_title(acortar_nombres_microbioma(Trait_Nombre))  
  )


#  HEATMAP ----

matriz_top2 <- tabla_top2 %>%
  dplyr::select(Termino_label, Cluster_Origen, log_padj) %>%
  group_by(Termino_label, Cluster_Origen) %>%
  summarise(log_padj = min(max(log_padj), 8), .groups = "drop") %>%
  pivot_wider(
    names_from  = Cluster_Origen,
    values_from = log_padj,
    values_fill = 0
  ) %>%
  column_to_rownames("Termino_label")
cat("Rango de la matriz:", range(as.matrix(matriz_top2)), "\n")

# Anotación columnas
efo_por_cluster <- gsub("^(EFO_[0-9]+)_Cluster.*$", "\\1", colnames(matriz_top2))
efo_a_nombre_full <- setNames(traits_MicroGWAS_areas$name, traits_MicroGWAS_areas$Rasgo)
nombres_reales <- stringr::str_to_title(acortar_nombres_microbioma(efo_a_nombre_full[efo_por_cluster]))

annotation_col <- data.frame(
  Trait     = nombres_reales,
  row.names = colnames(matriz_top2)
)

# Colores de anotación
efo_a_nombre_corto <- setNames(
  stringr::str_to_title(acortar_nombres_microbioma(efo_a_nombre_full[names(colores_6_traits)])),
  names(colores_6_traits))

colores_globales <- setNames(
  colores_6_traits[names(efo_a_nombre_corto)],
  efo_a_nombre_corto)
annotation_colors_interpro <- list(Trait = colores_globales)

colnames_cortos <- gsub("^EFO_[0-9]+_", "", colnames(matriz_top2))

paleta <- colorRampPalette(c("#EEF9C4", "#7FCDBB", "#2C7FB8", "#1A4E88"))(100)
# Plot
pdf("./Output/Gráficos/MicroGWAS/Heatmap_Enrich_InterPro_top2_ClustersMicroGWAS.pdf",
    height = 20, width = 20)
print(pheatmap(
  mat               = as.matrix(matriz_top2),
  color             = paleta,
  breaks            = seq(0, 8, length.out = 101),
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors_interpro,
  labels_col        = colnames_cortos,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 10,
  fontsize_col      = 9,
  cellheight        = 10,
  cellwidth         = 15,
  angle_col         = 45,
  main              = "InterPro family enrichment by cluster — Top 3 per cluster (-log10 p.adjust, capped at 8)",
  border_color      = "white"
))
dev.off()
message("✔️ Heatmap top2 guardado.")

#  DOTPLOT  ----

datos_dotplot_top2 <- tabla_top2 %>%
  dplyr::mutate(Cluster_label = paste0("(", Cluster_corto, ")")) %>%
  tidyr::complete(
    Termino_label,
    tidyr::nesting(Cluster_Origen, Cluster_label, Cluster_corto,
                   Trait_Corto, Trait_Nombre, efo_id) )

# Ordenar términos: más significativos arriba
orden_terminos <- tabla_top2 %>%
  group_by(Termino_label) %>%
  summarise(padj_medio = mean(p.adjust, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(padj_medio)) %>%
  pull(Termino_label)

datos_dotplot_top2$Termino_label <- factor(
  datos_dotplot_top2$Termino_label,
  levels = orden_terminos)

ggplot(datos_dotplot_top2, aes(x = Cluster_label, y = Termino_label)) +
  geom_point(aes(size = GeneRatio_val, color = p.adjust)) +
  scale_color_viridis_c(
    option    = "magma",
    direction = 1,
    trans     = "log10",
    name      = "p.adjust",
    na.value  = "grey90"     # celdas sin valor = gris claro
  ) +
  scale_size_continuous(
    name  = "GeneRatio",
    range = c(4, 10)
  ) +
  facet_grid(~ Trait_Corto, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
    axis.text.y        = element_text(size = 18, color = "black"),
    axis.title         = element_blank(),
    plot.title         = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5, size = 10, color = "grey40"),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 16, face = "bold", color = "black"),
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    panel.spacing      = unit(0.05, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  ) +
  labs(
    title    = "InterPro family enrichment by MicroGWAS Cluster",
    subtitle = "Top 2 terms per cluster | Simplified | p.adjust < 0.05 (BH)"
  )

ggsave(
  filename  = "./Output/Gráficos/MicroGWAS/Dotplot_Enrich_InterPro_top2_ClustersMicroGWAS.pdf",
  width     = 18,
  height    = 12,
  dpi       = 300,
  limitsize = FALSE
)
message("✔️ Dotplot top2 guardado.")
