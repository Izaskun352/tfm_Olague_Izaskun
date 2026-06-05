
library(ggraph)
library(tidygraph)

# FILTRAR COMPONENTES RED PLEIOTROPIA MICROGWAS - VARIACION COMUN

# Aislar los componentes de la red de pleiotropia MicroGWAS - GWAS

# INPUTS

red_pleiotropia_varComun_MicroGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Red_Pleiotropia.rds")
comp <- components(red_pleiotropia_varComun_MicroGWAS)
carpeta_red <- "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun"
comp$csize

# FUNCION CREAR RED

plot_red_pleiotropia <- function(
    tbl_red,
    layout_df,
    umbral_jaccard = c("0.5", "0.7"),   # "0.5" = ambas categorías, "0.7" = solo >= 0.7
    titulo    = "Red de pleiotropía MicroGWAS × VarComun",
    subtitulo = "Jaccard >= 0.5 | Componentes con al menos un trait MicroGWAS",
    colores_nodo = c(
      "MicroGWAS"          = "palegreen4",
      "NervousPsychiatric" = "peachpuff4",
      "Immune"             = "slateblue4",
      "VarComun"           = "black"
    ),
    colores_arista = c(
      "Jaccard >= 0.7"  = "#E63946",
      "Jaccard 0.5-0.7" = "#457B9D"
    ),
    rango_grosor = c(0.5, 2.5),
    alpha_arista = 0.7,
    tam_nodo     = 6,
    tam_texto    = 3,
    max_overlaps = 20
) {
  umbral_jaccard <- match.arg(umbral_jaccard)
  
  # Filtrar aristas según umbral
  if (umbral_jaccard == "0.7") {
    tbl_red <- tbl_red %>%
      activate(edges) %>%
      filter(categoria_jaccard == "Jaccard >= 0.7")
    
    colores_arista <- colores_arista["Jaccard >= 0.7"]
  }
  
  ggraph(tbl_red, layout = "manual", x = layout_df$x, y = layout_df$y) +
    geom_edge_link(
      aes(color = categoria_jaccard, width = Indice_Jaccard),
      alpha = alpha_arista
    ) +
    geom_node_point(
      aes(color = color_nodo),
      fill   = "white",
      shape  = 21,
      size   = tam_nodo,
      stroke = 2
    ) +
    geom_node_text(
      aes(label = label, color = color_nodo),
      repel        = TRUE,
      size         = tam_texto,
      max.overlaps = max_overlaps,
      family       = "sans"
    ) +
    scale_edge_color_manual(
      values = colores_arista,
      name   = "Fuerza de solapamiento"
    ) +
    scale_edge_width_continuous(range = rango_grosor, guide = "none") +
    scale_color_manual(
      values = colores_nodo,
      name   = "Tipo de trait"
    ) +
    theme_graph() +
    theme(
      legend.position   = "right",
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.title      = element_text(face = "bold"),
      legend.text       = element_text(size = 10)
    ) +
    labs(title = titulo, subtitle = subtitulo)
}

# ===============================================================
# 1- COMPONENTE GOBP = TRADUCCION ====

# Aislar el componente que queremos

red_traduccion <- induced_subgraph(red_pleiotropia_varComun_MicroGWAS, vids = which(comp$membership == 1))

# layout 
e <- as_edgelist(red_traduccion, names = F)
layout_red_traduccion <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(red_traduccion),
                                                area = 4 * (vcount(red_traduccion)^2),
                                                repulse.rad = vcount(red_traduccion)^3)

layout_red_traduccion <- data.frame(
  x = layout_red_traduccion[, 1],
  y = layout_red_traduccion[, 2]
)

# Crear el tbl_graph
tbl_red_traduccion <- tidygraph::as_tbl_graph(red_traduccion)

# Hacer red
plot_red_traduccion <- plot_red_pleiotropia(tbl_red_traduccion, layout_red_traduccion, umbral_jaccard = "0.7", titulo = "Red Traduccion, Jaccard > 0,7")

# Guardar
pdf(file.path(carpeta_red, "Red_Traduccion.pdf"), 
    width = 15, height = 15)
print(plot_red_traduccion)  
dev.off()

vertex_attr(red_pleiotropia_varComun_MicroGWAS)

# ===============================================================







