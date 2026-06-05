



# INPUTS

ruta_atlas <- "./Data/Tissue_Expression_Atlas/association_scores/"
archivos_atlas <- list.files(path = ruta_atlas, pattern = "\\.csv$", full.names = TRUE)
nombres_tejidos <- basename(archivos_atlas) %>%
  str_replace("cohorts_combined_", "") %>%
  str_replace("_avg_outer_prob.csv", "")
df_conversion <- readRDS("./Data/Tissue_Expression_Atlas/conversion_df_ENSG.rds")

red_pleiotropia_varComun_MicroGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Red_Pleiotropia.rds")
comp <- components(red_pleiotropia_varComun_MicroGWAS)

pares_todos <- igraph::as_data_frame(red_pleiotropia_varComun_MicroGWAS, what = "edges")
pares_todos <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Todos_Los_Pares_Jaccard.rds")


Clusters_Immune <- read.csv2("./Output/Piloto_Microbiota/Clusters_Immune/ZSCO.EFO_0000094_Cluster_1.2.csv")
Clusters_MicroGWAS <- "./Output/Piloto_Microbiota/Clusters_MicroGWAS"
Clusters_NervousPsychiatric <- "./Output/Piloto_Microbiota/Clusters_NervousPsychiatric"
Clusters_Traits_VarComun <- "./Output/Piloto_Microbiota/Clusters_Traits_VarComun"

# hacer df con los genes de todos los clusters

base_dir <- "./Output/Piloto_Microbiota"
carpetas <- c(
  "Clusters_Immune",
  "Clusters_MicroGWAS",
  "Clusters_NervousPsychiatric",
  "Clusters_Traits_VarComun"
)

rutas_carpetas <- file.path(base_dir, carpetas)
archivos_clusters <- list.files(path = rutas_carpetas, pattern = "\\.csv$", 
                                full.names = TRUE, recursive = FALSE)

#funcion para editar cada archivo
procesar_archivo <- function(ruta) {
  df <- read.csv2(ruta)
  nombre_archivo <- file_path_sans_ext(basename(ruta))
  
  df_limpio <- df %>%
    select(gene.gene, gene.ENSG) %>%      
    mutate(Cluster = nombre_archivo)    

  return(df_limpio)
}

# Aplicar la función a todos los archivos y unirlos en un único dataframe
# map_dfr aplica la función a la lista y hace un bind_rows() automáticamente
cluster_all <- map_dfr(archivos_csv, procesar_archivo)

# Leer todos los tejidos y convertir a ENSG
atlas_lista <- map2(archivos_atlas, nombres_tejidos, function(archivo, tejido) {
  read_csv(archivo, show_col_types = FALSE) %>%
    rename(prot1 = 1, prot2 = 2, score = 3) %>%
    left_join(df_conversion, by = c("prot1" = "to_id")) %>% rename(ENSG1 = from_id) %>%
    left_join(df_conversion, by = c("prot2" = "to_id")) %>% rename(ENSG2 = from_id) %>%
    filter(!is.na(ENSG1), !is.na(ENSG2)) %>%
    mutate(tejido = tejido)})
atlas_df <- bind_rows(atlas_lista)

# pares posibles por componente

df_componentes <- tibble(
  Cluster    = names(comp$membership),
  componente = comp$membership)
df_genes_comp <- cluster_all %>%
  left_join(df_componentes, by = c("Cluster" = "Cluster")) %>%
  filter(!is.na(componente)) %>%
  distinct(gene.ENSG, componente)

genes_por_componente <- df_genes_comp %>%
  group_by(componente) %>%
  summarise(genes = list(unique(gene.ENSG)), .groups = "drop")

pares_por_componente <- genes_por_componente %>%
  mutate(pares = map2(genes, componente, function(gs, comp) {
    if (length(gs) < 2) return(tibble(ENSG1 = character(), ENSG2 = character()))
    combn(gs, 2, simplify = FALSE) %>%
      map_dfr(~tibble(ENSG1 = .x[1], ENSG2 = .x[2])) %>%
      mutate(componente = comp)
  })) %>%
  select(pares) %>%
  unnest(pares)
# definimos niverso - todas las intersecciones (únicas) de todos los tejidos

universo_pares <- atlas_df %>%
  distinct(ENSG1, ENSG2)
universo_pares <- universo_pares %>%   # normalizar
  mutate(
    g1 = if_else(ENSG1 < ENSG2, ENSG1, ENSG2),
    g2 = if_else(ENSG1 < ENSG2, ENSG2, ENSG1)) %>%
  distinct(g1, g2)

N_pares <- nrow(universo_pares)


#ponemos las parejas siempre en el mismo orden (normalizar)

UMBRAL = 0.8
pares_por_tejido <- atlas_df %>%
  filter(score >= UMBRAL) %>%
  mutate(
    g1 = if_else(ENSG1 < ENSG2, ENSG1, ENSG2),
    g2 = if_else(ENSG1 < ENSG2, ENSG2, ENSG1)
  ) %>%
  distinct(tejido, g1, g2)

pares_por_componente <- pares_por_componente %>%
  mutate(
    g1 = if_else(ENSG1 < ENSG2, ENSG1, ENSG2),
    g2 = if_else(ENSG1 < ENSG2, ENSG2, ENSG1)
  ) %>%
  distinct(componente, g1, g2)



# test de fisher 

n_componentes <- sort(unique(pares_por_componente$componente))
n_tejidos     <- unique(pares_por_tejido$tejido)

resultados_fisher <- expand_grid(  #crear todas las combinaciones
  componente = n_componentes,
  tejido     = n_tejidos
) %>%
  mutate(fisher = map2(componente, tejido, function(comp, tej) {
    
    pares_comp   <- pares_por_componente %>% 
      dplyr::filter(componente == comp) %>%
      dplyr::select(g1, g2)
    
    pares_tej    <- pares_por_tejido %>% 
      dplyr::filter(tejido == tej) %>%
      dplyr::select(g1, g2)
    
    # pares SI componente y SI tejido
    a <- nrow(inner_join(pares_comp, pares_tej, by = c("g1", "g2")))
    
    if (a >= 10) {  #mínimo 10 interacciones
      # pares NO componente y SI tejido
      b <- nrow(anti_join(pares_tej, pares_comp, by = c("g1", "g2")))
      # c: pares SI componente NO tejido
      c <- nrow(anti_join(pares_comp, pares_tej, by = c("g1", "g2")))
      # d: pares NO componente NO tejido
      d <- N_pares - a - b - c
      
      mat <- matrix(c(a, c, b, d), nrow = 2,
                    dimnames = list(c("En_comp", "Fuera_comp"),
                                    c("Expresado", "No_expresado")))
      ft <- fisher.test(mat, alternative = "greater")
      or_val <- ft$estimate
      p_val  <- ft$p.value
    } else {  # si hay menos de 10 no calculamos fisher y devolvemos NA
      or_val <- NA_real_
      p_val  <- NA_real_
    }
    
    tibble(
      n_pares_comp   = nrow(pares_comp),
      n_pares_tejido = nrow(pares_tej),
      overlap        = a,
      odds_ratio     = or_val,   
      p_value        = p_val
    )
  })) %>%
  unnest(fisher)

resultados_fisher <- resultados_fisher %>%
  dplyr::group_by(componente) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(componente, p_adj_BH)
saveRDS(resultados_fisher, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Resultados_Fisher_Expresion_Tejidos.rds")
#heatmap

df_componentes_info <- tibble(
  componente = 1:comp$no,
  csize      = comp$csize)

df_heatmap <- resultados_fisher %>%
  left_join(df_componentes_info, by = "componente") %>%
  mutate(
    comp_label    = paste0("Comp. ", componente, "\n(", csize, " nodos)"),
    odds_ratio_cap = pmin(odds_ratio, 4),
    odds_ratio_cap = ifelse(p_adj_BH < 0.05, odds_ratio_cap, NA)
  )

ggplot(df_heatmap, aes(x = factor(comp_label), y = tejido, fill = odds_ratio_cap)) +
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
    title    = "Enriquecimiento de tejido por componente",
    x        = "Componente",
    y        = "Tejido"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y      = element_text(size = 14),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 14),
    legend.position  = "right"
  )

ggsave("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Heatmap_Fisher_Componentes_Tejidos.pdf", 
       width = 6.5, height = 7, dpi = 300)


