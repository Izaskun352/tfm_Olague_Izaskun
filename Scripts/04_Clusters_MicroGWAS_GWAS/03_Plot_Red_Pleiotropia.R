

# CONSTRUIR RED DE PLEIOTROPIA (Clusters MicroGWAS vs Clusters GWAS) 

# Clusters MicroGWAS - verdes
# Clusters Nervous / Psychiatric - marron
# Clusters Immune - azul
# Clusters GWAS - negro

# LIBRERIAS

source("scripts/00_setup.R")
library(ggraph)

# OUTPUTS
red_pleiotropia_varComun_MicroGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Red_Pleiotropia.rds")

# 1. Construir diccionario ID -> nombre legible
dicc_nombres <- traits_MicroGWAS_areas %>%
  dplyr::select(Rasgo, name) %>%
  dplyr::distinct()

get_nombre <- function(nombre_cluster) {
  trait_id  <- sub("_Cluster_.*$", "", nombre_cluster)  # quitar _Cluster_X
  id_limpio <- sub("^ZSCO\\.", "", trait_id)             # quitar ZSCO.
  nombre    <- dicc_nombres %>%
    dplyr::filter(Rasgo == id_limpio) %>%
    dplyr::pull(name) %>%
    dplyr::first()
  if (length(nombre) == 0 || is.na(nombre)) return(id_limpio)
  return(nombre)
}

# 2. Cargar pares (Jaccard >= 0.5)

# --- MicroGWAS vs VarComun ---
pares_micro_vc <- read.csv2(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_MicroGWAS/Jaccard_Alta_Significativa_VarComun_MicroGWAS.csv"
) %>%
  dplyr::transmute(g1 = Cluster_MicroGWAS, g2 = Cluster_Comun,
                   Indice_Jaccard, tipo_par = "MicroGWAS_VarComun")

# --- MicroGWAS vs NervousPsychiatric + NP vs NP ---
pares_np_raw <- read.csv2(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterNervousPsychiatric/Jaccard_Completa_Interseccion_NP_MicroGWAS.csv"
) %>% dplyr::filter(Indice_Jaccard >= 0.5)

pares_micro_np <- pares_np_raw %>%
  dplyr::filter(
    (!grepl("^ZSCO\\.", Cluster_1) &  grepl("^ZSCO\\.", Cluster_2)) |
      ( grepl("^ZSCO\\.", Cluster_1) & !grepl("^ZSCO\\.", Cluster_2))
  ) %>%
  dplyr::transmute(g1 = Cluster_1, g2 = Cluster_2,
                   Indice_Jaccard, tipo_par = "MicroGWAS_NervousPsychiatric")

pares_np_np <- pares_np_raw %>%
  dplyr::filter(grepl("^ZSCO\\.", Cluster_1) & grepl("^ZSCO\\.", Cluster_2)) %>%
  dplyr::transmute(g1 = Cluster_1, g2 = Cluster_2,
                   Indice_Jaccard, tipo_par = "NervousPsychiatric_NervousPsychiatric")

# --- MicroGWAS vs Immune + Immune vs Immune ---
pares_immune_raw <- read.csv2(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_ClusterImmune/Jaccard_Completa_Interseccion_Immune_MicroGWAS.csv"
) %>% dplyr::filter(Indice_Jaccard >= 0.5)

pares_micro_immune <- pares_immune_raw %>%
  dplyr::filter(
    (!grepl("^ZSCO\\.", Cluster_1) &  grepl("^ZSCO\\.", Cluster_2)) |
      ( grepl("^ZSCO\\.", Cluster_1) & !grepl("^ZSCO\\.", Cluster_2))
  ) %>%
  dplyr::transmute(g1 = Cluster_1, g2 = Cluster_2,
                   Indice_Jaccard, tipo_par = "MicroGWAS_Immune")

pares_immune_immune <- pares_immune_raw %>%
  dplyr::filter(grepl("^ZSCO\\.", Cluster_1) & grepl("^ZSCO\\.", Cluster_2)) %>%
  dplyr::transmute(g1 = Cluster_1, g2 = Cluster_2,
                   Indice_Jaccard, tipo_par = "Immune_Immune")

# --- VarComun vs VarComun ---
pares_vc_vc <- readRDS(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_VarComun_VarComun/Pleiotropia_Jaccard_Alta_Significativa.rds") %>%
  dplyr::transmute(g1, g2, Indice_Jaccard, tipo_par = "VarComun_VarComun")

# --- NervousPsychiatric vs VarComun ---
pares_np_vc <- readRDS(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_NervousPsychiatric_VarComun/Pleiotropia_Jaccard_Alta_Significativa.rds") %>%
  dplyr::transmute(g1 = Cluster_NervousPsychiatric, g2 = Cluster_VarComun,
                   Indice_Jaccard, tipo_par = "NervousPsychiatric_VarComun")

# --- Immune vs VarComun ---
pares_immune_vc <- readRDS(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_VarComun/Pleiotropia_Jaccard_Alta_Significativa.rds") %>%
  dplyr::transmute(g1 = Cluster_Immune, g2 = Cluster_VarComun,
                   Indice_Jaccard, tipo_par = "Immune_VarComun")

# --- Immune vs NervousPsychiatric ---
pares_immune_np <- readRDS(
  "./Output/Piloto_Microbiota/Modulos_Pleiotropicos_Immune_NervousPsychiatric/Pleiotropia_Jaccard_Alta_Significativa.rds") %>%
  dplyr::transmute(g1 = Cluster_Immune, g2 = Cluster_NervousPsychiatric,
                   Indice_Jaccard, tipo_par = "Immune_NervousPsychiatric")
# 3. Unir todos los pares
pares_todos <- dplyr::bind_rows(
  pares_micro_vc,
  pares_micro_np,
  pares_np_np,
  pares_micro_immune,
  pares_immune_immune,
  pares_vc_vc,
  pares_np_vc,
  pares_immune_vc,
  pares_immune_np)

saveRDS(pares_todos, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Todos_Los_Pares_Jaccard.rds")

# 4. Identificar nodos por tipo
nodos_microGWAS <- unique(c(pares_micro_vc$g1, pares_micro_np$g1, pares_micro_immune$g2))
nodos_np        <- unique(c(pares_micro_np$g2, pares_np_np$g1, pares_np_np$g2,
                            pares_np_vc$g1, pares_immune_np$g2))
nodos_immune    <- unique(c(pares_micro_immune$g1, pares_immune_immune$g1,
                            pares_immune_immune$g2, pares_immune_vc$g1, pares_immune_np$g1))


# 5. Construir red inicial
red <- igraph::graph_from_data_frame(
  d        = pares_todos %>% dplyr::select(g1, g2, Indice_Jaccard, tipo_par),
  directed = FALSE
)

# 6. Filtrar: quedarse solo con componentes que tengan >= 1 nodo MicroGWAS
componentes <- igraph::components(red)
comp_ids    <- componentes$membership

# Identificar qué componentes tienen al menos un nodo MicroGWAS
comps_con_micro <- unique(comp_ids[names(comp_ids) %in% nodos_microGWAS])

message("Componentes con nodo MicroGWAS: ", length(comps_con_micro))
message("Nodos MicroGWAS encontrados en la red: ", 
        sum(names(comp_ids) %in% nodos_microGWAS))

nodos_validos <- names(comp_ids[comp_ids %in% comps_con_micro])
red_filtrada  <- igraph::induced_subgraph(red, vids = nodos_validos)
message("Componentes con MicroGWAS: ", length(comps_con_micro))

# 7. Filtro adicional: nodos VarComun deben conectar con MicroGWAS Y con >= 1 VarComun
nodos_a_eliminar <- c()
for (nodo in igraph::V(red_filtrada)$name) {
  if (!nodo %in% nodos_microGWAS) {  # es VarComun
    vecinos <- names(igraph::neighbors(red_filtrada, nodo))
    tiene_micro <- any(vecinos %in% nodos_microGWAS)
    tiene_otros    <- any(!vecinos %in% nodos_microGWAS)
    if (!tiene_micro || !tiene_otros) {
      nodos_a_eliminar <- c(nodos_a_eliminar, nodo)
    }
  }
}

if (length(nodos_a_eliminar) > 0) {
  red_filtrada <- igraph::delete_vertices(red_filtrada, nodos_a_eliminar)
  message("Nodos VarComun eliminados por no cumplir condición: ", length(nodos_a_eliminar))
}

message("Nodos en la red final: ", igraph::vcount(red_filtrada))
message("Aristas en la red final: ", igraph::ecount(red_filtrada))

# 8. Añadir atributos de nodos
tipo_nodo <- ifelse(igraph::V(red_filtrada)$name %in% nodos_microGWAS, "MicroGWAS",
                    ifelse(igraph::V(red_filtrada)$name %in% nodos_np,         "NervousPsychiatric",
                           ifelse(igraph::V(red_filtrada)$name %in% nodos_immune,     "Immune",
                                  "VarComun")))

igraph::V(red_filtrada)$tipo_nodo <- tipo_nodo
igraph::V(red_filtrada)$label     <- sapply(igraph::V(red_filtrada)$name, get_nombre)

# 9. Añadir atributos de aristas
igraph::E(red_filtrada)$categoria_jaccard <- ifelse(
  igraph::E(red_filtrada)$Indice_Jaccard >= 0.7, "Jaccard >= 0.7",
  "Jaccard 0.5-0.7")

# 10. Añadimos layout

e <- as_edgelist(red_filtrada, names = F)
layout_red <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(red_filtrada),
                                                area = 4 * (vcount(red_filtrada)^2),
                                                repulse.rad = vcount(red_filtrada)^3)
# Convertir el layout a dataframe con nombres de columna x e y
layout_df <- data.frame(
  x = layout_red[, 1],
  y = layout_red[, 2])

# Crear el tbl_graph
tbl_red <- tidygraph::as_tbl_graph(red_filtrada)

# 11. Visualizar con ggraph

plot_red <- ggraph(tbl_red, layout = "manual", x = layout_df$x, y = layout_df$y) +
  geom_edge_link(
    aes(color = categoria_jaccard, width = Indice_Jaccard),
    alpha = 0.7
  ) +
  geom_node_point(
    aes(color = tipo_nodo),  # color del perímetro
    fill   = "white",          # interior blanco
    shape  = 21,               # círculo con borde
    size   = 6,
    stroke = 2                 # grosor del borde
  ) +
  
  geom_node_text(
    aes(label = label, color = tipo_nodo),
    repel        = TRUE,
    size         = 3,
    max.overlaps = 20,
    family       = "sans"
  ) +
  scale_edge_color_manual(
    values = c("Jaccard >= 0.7"  = "#E63946",
               "Jaccard 0.5-0.7" = "#457B9D"),
    name   = "Fuerza de solapamiento"
  ) +
  scale_edge_width_continuous(range = c(0.5, 2.5), guide = "none") +
  scale_color_manual(
    values = c("MicroGWAS"          = "palegreen4",
               "NervousPsychiatric" = "peachpuff4",
               "Immune"             = "slateblue4",
               "VarComun"           = "black"),
    name   = "Tipo de trait"
  ) +
  theme_graph() +
  theme(
    legend.position   = "right",
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(size = 10)
  ) +
  labs(title = "Red de pleiotropía MicroGWAS × VarComun",
       subtitle = "Jaccard >= 0.5 | Componentes con al menos un trait MicroGWAS")

# 12. Exportar
carpeta_red <- "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun"
dir.create(carpeta_red, showWarnings = FALSE)

pdf(file.path(carpeta_red, "Red_Pleiotropia_2.pdf"), 
    width = 30, height = 30)
print(plot_red)  
dev.off()

png(file.path(carpeta_red, "Red_Pleiotropia_2.png"),
    width = 6000, height = 6000, res = 300)
print(plot_red)
dev.off()

# 13. Guardar
saveRDS(red_filtrada, file.path(carpeta_red, "Red_Pleiotropia.rds"))
red_pleiotropia_varComun_MicroGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Red_Pleiotropia.rds")
igraph::components(red_pleiotropia_varComun_MicroGWAS)
vertex_attr_names(red_pleiotropia_varComun_MicroGWAS)
str(red_pleiotropia_varComun_MicroGWAS)