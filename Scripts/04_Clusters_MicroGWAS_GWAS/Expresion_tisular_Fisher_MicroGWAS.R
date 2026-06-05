ruta_clusters <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
archivos_clusters <- list.files(path = ruta_clusters, pattern = "\\.csv$", 
                                full.names = TRUE, recursive = FALSE)

cluster_all <- map_dfr(archivos_clusters, function(ruta) {
  df <- read.csv2(ruta)
  nombre <- file_path_sans_ext(basename(ruta))
  df %>%
    select(gene.ENSG, gene.gene) %>%
    mutate(Cluster = nombre)
})

# ── 2. Cargar atlas (igual que antes) ────────────────────────────────────────
ruta_atlas <- "./Data/Tissue_Expression_Atlas/association_scores/"
archivos_atlas <- list.files(path = ruta_atlas, pattern = "\\.csv$", full.names = TRUE)
nombres_tejidos <- basename(archivos_atlas) %>%
  str_replace("cohorts_combined_", "") %>%
  str_replace("_avg_outer_prob.csv", "")

df_conversion <- readRDS("./Data/Tissue_Expression_Atlas/conversion_df_ENSG.rds")

atlas_lista <- map2(archivos_atlas, nombres_tejidos, function(archivo, tejido) {
  read_csv(archivo, show_col_types = FALSE) %>%
    rename(prot1 = 1, prot2 = 2, score = 3) %>%
    left_join(df_conversion, by = c("prot1" = "to_id")) %>% rename(ENSG1 = from_id) %>%
    left_join(df_conversion, by = c("prot2" = "to_id")) %>% rename(ENSG2 = from_id) %>%
    filter(!is.na(ENSG1), !is.na(ENSG2)) %>%
    mutate(tejido = tejido)
})
atlas_df <- bind_rows(atlas_lista)

# ── 3. Pares por cluster ──────────────────────────────────────────────────────
genes_por_cluster <- cluster_all %>%
  group_by(Cluster) %>%
  summarise(genes = list(unique(gene.ENSG)), .groups = "drop")

pares_por_cluster <- genes_por_cluster %>%
  mutate(pares = map2(genes, Cluster, function(gs, cl) {
    if (length(gs) < 2) return(tibble(ENSG1 = character(), ENSG2 = character()))
    combn(gs, 2, simplify = FALSE) %>%
      map_dfr(~tibble(ENSG1 = .x[1], ENSG2 = .x[2])) %>%
      mutate(Cluster = cl)
  })) %>%
  select(pares) %>%
  unnest(pares) %>%
  mutate(
    g1 = if_else(ENSG1 < ENSG2, ENSG1, ENSG2),
    g2 = if_else(ENSG1 < ENSG2, ENSG2, ENSG1)
  ) %>%
  distinct(Cluster, g1, g2)

# ── 4. Universo y pares por tejido (igual que antes) ─────────────────────────
universo_pares <- atlas_df %>%
  mutate(
    g1 = if_else(ENSG1 < ENSG2, ENSG1, ENSG2),
    g2 = if_else(ENSG1 < ENSG2, ENSG2, ENSG1)
  ) %>%
  distinct(g1, g2)

N_pares <- nrow(universo_pares)

UMBRAL <- 0.8
pares_por_tejido <- atlas_df %>%
  filter(score >= UMBRAL) %>%
  mutate(
    g1 = if_else(ENSG1 < ENSG2, ENSG1, ENSG2),
    g2 = if_else(ENSG1 < ENSG2, ENSG2, ENSG1)
  ) %>%
  distinct(tejido, g1, g2)

# ── 5. Test de Fisher: 27 clusters × 11 tejidos = 297 tests ──────────────────
n_clusters <- unique(pares_por_cluster$Cluster)
n_tejidos  <- unique(pares_por_tejido$tejido)

resultados_fisher_clusters <- expand_grid(
  Cluster = n_clusters,
  tejido  = n_tejidos
) %>%
  mutate(fisher = map2(Cluster, tejido, function(cl, tej) {
    
    pares_comp <- pares_por_cluster %>%
      dplyr::filter(Cluster == cl) %>%
      dplyr::select(g1, g2)
    
    pares_tej <- pares_por_tejido %>%
      dplyr::filter(tejido == tej) %>%
      dplyr::select(g1, g2)
    
    a <- nrow(inner_join(pares_comp, pares_tej, by = c("g1", "g2")))
    
    if (a >= 10) {
      b <- nrow(anti_join(pares_tej,  pares_comp, by = c("g1", "g2")))
      c <- nrow(anti_join(pares_comp, pares_tej,  by = c("g1", "g2")))
      d <- N_pares - a - b - c
      
      mat <- matrix(c(a, c, b, d), nrow = 2,
                    dimnames = list(c("En_cluster", "Fuera_cluster"),
                                    c("Expresado",  "No_expresado")))
      ft     <- fisher.test(mat, alternative = "greater")
      or_val <- ft$estimate
      p_val  <- ft$p.value
    } else {
      or_val <- NA_real_
      p_val  <- NA_real_
    }
    
    tibble(
      n_pares_cluster = nrow(pares_comp),
      n_pares_tejido  = nrow(pares_tej),
      overlap         = a,
      odds_ratio      = or_val,
      p_value         = p_val
    )
  })) %>%
  unnest(fisher)

# ── 6. Corrección múltiple ────────────────────────────────────────────────────
resultados_fisher_clusters <- resultados_fisher_clusters %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Cluster, p_adj_BH)

saveRDS(resultados_fisher_clusters, 
        "./Output/Piloto_Microbiota/Clusters_MicroGWAS/Resultados_Fisher_Expresion_Tejidos.rds")

# ── 7. Heatmap ────────────────────────────────────────────────────────────────
# Número de genes por cluster para la etiqueta

# cabiar nombres

crear_mapeo_clusters <- function(carpeta) {
  
  archivos <- list.files(carpeta, pattern = "\\.csv$", full.names = FALSE)
  nombres_originales <- tools::file_path_sans_ext(archivos)
  
  df <- tibble(original = nombres_originales) %>%
    mutate(
      # Rasgo = todo antes de "_Cluster_"
      rasgo = str_extract(original, "^(.+?)(?=_Cluster_)"),
      # Número de cluster = todo después de "_Cluster_"
      num_cluster = str_extract(original, "(?<=_Cluster_)[\\d\\.]+$")
    )
  
  # Función para ordenación jerárquica: 1 < 1.1 < 1.1.7 < 1.2 < 2
  ordenar_jerarquico <- function(x) {
    partes <- str_split(x, "\\.")[[1]]
    paste(str_pad(partes, 4, pad = "0"), collapse = ".")
  }
  
  df <- df %>%
    mutate(orden = sapply(num_cluster, ordenar_jerarquico)) %>%
    arrange(rasgo, orden) %>%
    group_by(rasgo) %>%
    mutate(
      indice = row_number(),                                    # reinicia en 1 por rasgo
      nuevo_nombre = paste0(rasgo, "_Cluster_", indice)
    ) %>%
    ungroup()
  
  mapeo <- setNames(df$nuevo_nombre, df$original)
  
  return(mapeo)
}

carpeta_clusters_microGWAS <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
dic_clusters_MicroGWAS <- crear_mapeo_clusters(carpeta_clusters_microGWAS)
acortar_nombres_microbioma <- function(nombres) {
  n <- gsub(" microbiome measurement", " Micr.", nombres, ignore.case = TRUE)
  return(n)}

n_genes_cluster <- cluster_all %>%
  group_by(Cluster) %>%
  summarise(n_genes = n_distinct(gene.ENSG), .groups = "drop")

df_heatmap_clusters <- resultados_fisher_clusters %>%
  left_join(n_genes_cluster, by = "Cluster") %>%
  mutate(
    # Renombrar con el diccionario jerárquico
    Cluster_renombrado = dic_clusters_MicroGWAS[Cluster],
    # Extraer el trait del nombre renombrado
    Rasgo      = str_extract(Cluster_renombrado, "^(.+?)(?=_Cluster_)"),
    # Extraer el número ya reordenado
    num_cluster = str_extract(Cluster_renombrado, "(?<=_Cluster_)[\\d]+$")
  ) %>%
  left_join(traits_MicroGWAS_areas %>% select(Rasgo, name), by = "Rasgo") %>%
  mutate(
    nombre_trait   = if_else(is.na(name), Rasgo, name),
    nombre_trait   = acortar_nombres_microbioma(nombre_trait),
    cluster_label  = paste0(nombre_trait, " (", num_cluster, ")\n(", n_genes, " genes)"),
    # Ordenar el eje x por rasgo y número de cluster
    cluster_label  = fct_reorder(cluster_label, as.numeric(num_cluster)),
    odds_ratio_cap = pmin(odds_ratio, 4),
    odds_ratio_cap = ifelse(p_adj_BH < 0.05, odds_ratio_cap, NA)
  )

ggplot(df_heatmap_clusters, aes(x = factor(cluster_label), y = tejido, fill = odds_ratio_cap)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low      = "peachpuff",
    mid      = "peachpuff3",
    high     = "skyblue4",
    midpoint = 2,
    limits   = c(1, 4),
    name     = "Odds Ratio",
    na.value = "grey80"
  ) +
  labs(
    title = "Enriquecimiento de tejido por cluster",
    x     = "Cluster",
    y     = "Tejido"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 12),
    panel.grid      = element_blank(),
    plot.title      = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("./Output/Piloto_Microbiota/Heatmap_Fisher_Clusters_MicroGWAS.pdf",
       width = 16, height = 7, dpi = 300)
