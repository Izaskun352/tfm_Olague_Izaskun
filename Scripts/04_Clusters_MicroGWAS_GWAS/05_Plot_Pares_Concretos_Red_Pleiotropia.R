
library(ggraph)
library(tidygraph)

# CONSTRUIR RED DE LA INTERSECCION DE DOS CLUSTERS - DE UN PAR CONCRETO

# -- A partir de la red pleiotropia MicroGWAS - GWAS

# INPUTS

red_pleiotropia_varComun_MicroGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Red_Pleiotropia.rds")
carpeta_red <- "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun"
pares_todos <- igraph::as_data_frame(red_pleiotropia_varComun_MicroGWAS, what = "edges")

# FUNCION

generar_red_par <- function(par_nombre, interactoma, lista_intersecciones, excluir_db = NULL) {
  
  # 1. Procesar el interactoma
  interactoma_multidb <- as.data.frame(interactoma) %>%
    dplyr::rename(ENSG_A = targetA, ENSG_B = targetB) %>%
    dplyr::filter(ENSG_A != ENSG_B) %>%
    { if (!is.null(excluir_db)) dplyr::filter(., !sourceDatabase %in% excluir_db) else . } %>%
    dplyr::mutate(
      par_id = paste(pmin(ENSG_A, ENSG_B), pmax(ENSG_A, ENSG_B), sep = "_")
    ) %>%
    dplyr::distinct(par_id, sourceDatabase, .keep_all = TRUE)  # único por par + DB
  
  prioridad_db <- c("intact" = 1, "reactome" = 2, "signor" = 3, "string" = 4)
  
  interactoma_multidb <- interactoma_multidb %>%
    dplyr::mutate(prioridad = prioridad_db[sourceDatabase]) %>%
    dplyr::group_by(par_id) %>%
    dplyr::mutate(
      dbs_presentes  = paste(sort(unique(sourceDatabase)), collapse = ";"),
      n_dbs          = n_distinct(sourceDatabase)
    ) %>%
    dplyr::slice_min(prioridad, n = 1, with_ties = FALSE) %>%  # quédate con la DB prioritaria
    dplyr::ungroup()
  
  
  # 2. Extraer los datos del par introducido
  par <- lista_intersecciones[[par_nombre]]
  c1      <- par$c1
  c2      <- par$c2
  jaccard <- par$jaccard
  
  tabla_genes <- par$genes_tabla
  
  # Limpiar Ensembl IDs
  genes_interseccion <- unique(gsub(";.*$", "", tabla_genes$gene.ENSG))
  
  if (length(genes_interseccion) == 0) {
    message("  ⚠️ Intersección vacía, saltando...")
    return(NULL) # Cambiado 'next' por 'return(NULL)' para salir de la función
  }
  
  message("  Genes en intersección: ", length(genes_interseccion))
  
  # 3. Filtrar aristas del interactoma solo para estos genes
  aristas <- interactoma_multidb %>%
    dplyr::filter(ENSG_A %in% genes_interseccion & ENSG_B %in% genes_interseccion)
  
  if (nrow(aristas) > 0) {
    red <- graph_from_data_frame(
      d        = aristas[, c("ENSG_A", "ENSG_B","sourceDatabase", "dbs_presentes", "n_dbs")],
      directed = FALSE,
      vertices = data.frame(name = genes_interseccion)
    )
  } else {
    red <- make_empty_graph(n = 0, directed = FALSE) %>%
      add_vertices(length(genes_interseccion), name = genes_interseccion)
    message("  ⚠️ Sin aristas entre genes — red solo con nodos")
  }
  
  # 4. Tabla con ENSG - nombre gen + categoria: semilla / propagado
  anotacion <- tabla_genes %>%
    dplyr::mutate(ENSG_limpio = gsub(";.*$", "", gene.ENSG)) %>%
    dplyr::select(ENSG_limpio, gene.gene, es_semilla_en_c1, es_semilla_en_c2) %>%
    dplyr::distinct(ENSG_limpio, .keep_all = TRUE)  %>% 
    dplyr::mutate(
      categoria_nodo = dplyr::case_when(
        es_semilla_en_c1 & es_semilla_en_c2  ~ "Ambos",
        es_semilla_en_c1 & !es_semilla_en_c2 ~ "Solo_C1",
        !es_semilla_en_c1 & es_semilla_en_c2 ~ "Solo_C2",
        TRUE                                  ~ "Propagado"
      )
    )
  
  # Ponemos a cada nodo el nobre del gen y su categoria (semilla / propagado)
  idx <- match(V(red)$name, anotacion$ENSG_limpio)
  V(red)$gene_name      <- anotacion$gene.gene[idx]
  V(red)$categoria_nodo <- anotacion$categoria_nodo[idx]
  
  # 5. layout 
  e <- as_edgelist(red, names = F)
  layout_red <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(red),
                                                  area = 4 * (vcount(red)^2),
                                                  repulse.rad = vcount(red)^3)
  
  # Convertir el layout a dataframe con nombres de columna x e y
  layout_df <- data.frame(
    x = layout_red[, 1],
    y = layout_red[, 2]
  )
  
  # Crear el tbl_graph
  tbl_red <- tidygraph::as_tbl_graph(red)
  
  # 6. Visualizar red
  plot_red <- ggraph(tbl_red, layout = "manual", x = layout_df$x, y = layout_df$y) +
    geom_edge_link(
      aes(color = sourceDatabase), 
      width = 1,                 
      alpha = 0.8
    ) +
    geom_node_point(
      aes(color = categoria_nodo),  # color del perímetro
      fill   = "white",          # interior blanco
      shape  = 21,               # círculo con borde
      size   = 8,
      stroke = 2                 # grosor del borde
    ) +
    
    geom_node_text(
      aes(label = gene_name, color = categoria_nodo),
      repel        = TRUE,
      size         = 5,
      max.overlaps = 20,
      family       = "sans"
    ) +
    scale_edge_color_manual(
      values = c("intact"   = "#E69F00",  # naranja
                 "signor"   = "#56B4E9",  # azul claro
                 "reactome" = "#009E73",  # verde
                 "string"   = "#CC79A7"),   # rosa/morado
      name   = "Base de datos"
    ) +
    scale_color_manual(
      values = c("Ambos"          = "slateblue4",
                 "Solo_C1" = "palegreen4",
                 "Solo_C2"             = "peachpuff4",
                 "Propagado"            = "black"),
      name   = "Tipo de trait"
    ) +
    theme_graph(base_family = "sans") +
    theme(
      legend.position   = "right",
      legend.display    = "block",
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.title      = element_text(face = "bold", size = 10),
      legend.text       = element_text(size = 9)
    ) +
    labs(title = par_nombre,
         subtitle = "Solamente Intact")
  
  # MUY IMPORTANTE: Devolver el objeto del gráfico al final
  return(plot_red)
}


# -----------------------------------------------------------------------
#  ---- oral vs hyperuricemia ----

par_nombre <- names(lista_intersecciones_vc)[
  sapply(lista_intersecciones_vc, function(x) {
    grepl("oral microbiome", x$trait_microGWAS, ignore.case = TRUE) &
      grepl("hyperuricemia",   x$trait_comun,      ignore.case = TRUE)
  })
]
print(par_nombre)

genes_interseccion <- par$genes_ensembl
writeLines(genes_interseccion, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/genes_interseccion_Oral_Hyperuricemia")

carpeta_red <- "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun"

plot_red <- generar_red_par(par_nombre = par_nombre, interactoma = interactoma_edges_withSource, lista_intersecciones = lista_intersecciones_vc)
pdf(file.path(carpeta_red, "Red_Oral_Hyperuricemia.pdf"), 
    width = 30, height = 30)
print(plot_red)  
dev.off()
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# ---- skin vs haplogobina ----

skin_haplogobina_name <- names(lista_intersecciones_vc)[
  sapply(lista_intersecciones_vc, function(x) {
    grepl("skin microbiome", x$trait_microGWAS, ignore.case = TRUE) &
      grepl("haptoglobin measurement",   x$trait_comun,      ignore.case = TRUE)
  })
]
print(skin_haplogobina_name)
skin_haplogobina <- lista_intersecciones_vc[[skin_haplogobina_name]]

genes_interseccion <- skin_haplogobina$genes_ensembl
writeLines(genes_interseccion, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/genes_interseccion_Skin_Haplogobina")

plot_red <- generar_red_par(par_nombre = skin_haplogobina_name, interactoma = interactoma_edges_withSource, lista_intersecciones = lista_intersecciones_vc)
pdf(file.path(carpeta_red, "Red_Skin_Haplogobina.pdf"), 
    width = 30, height = 30)
print(plot_red)  
dev.off()
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# ---- Gut vs neuroticsm measurement  ----

gut_neurotricsm_name <- names(lista_intersecciones_vc)[
  sapply(lista_intersecciones_vc, function(x) {
    grepl("gut microbiome", x$trait_microGWAS, ignore.case = TRUE) &
      grepl("neuroticsm measurement",   x$trait_comun,      ignore.case = TRUE)
  })
]
print(gut_neurotricsm_name)
gut_neurotricsm <- lista_intersecciones_vc[["EFO_0007874_Cluster_1.3.4.2_vs_ZSCO.EFO_0007660_Cluster_1.5.2"]]

genes_interseccion <- gut_neurotricsm$genes_ensembl
writeLines(genes_interseccion, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/genes_interseccion_Gut_Neuroticsm")

plot_red <- generar_red_par(par_nombre = "EFO_0007874_Cluster_1.3.4.2_vs_ZSCO.EFO_0007660_Cluster_1.5.2", interactoma = interactoma_edges_withSource, 
                            lista_intersecciones = lista_intersecciones_vc,
                            excluir_db = "string")
pdf(file.path(carpeta_red, "Red_Gut_Neuroticsm.pdf"), 
    width = 20, height = 20)
print(plot_red)  
dev.off()

# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# ---- Gut vs neuroticsm measurement vs intraocular pressure measurement ----

cluster_A <- "EFO_0007874_Cluster_1.3.4.2"  
cluster_B <- "ZSCO.EFO_0007660_Cluster_1.5.2"
cluster_C <- "ZSCO.EFO_0004695_Cluster_1.6.1"

# Par A-B
par_AB_name <- names(lista_intersecciones_vc)[
  sapply(lista_intersecciones_vc, function(x) {
    (x$c1 == cluster_A & x$c2 == cluster_B) |
      (x$c1 == cluster_B & x$c2 == cluster_A)})]
# Par A-C
par_AC_name <- names(lista_intersecciones_vc)[
  sapply(lista_intersecciones_vc, function(x) {
    (x$c1 == cluster_A & x$c2 == cluster_C) |
      (x$c1 == cluster_C & x$c2 == cluster_A)})]

# Par B-C
par_BC_name <- names(lista_intersecciones_vc_vc)[
  sapply(lista_intersecciones_vc_vc, function(x) {
    (x$c1 == cluster_B & x$c2 == cluster_C) |
      (x$c1 == cluster_C & x$c2 == cluster_B)})]

cat("Par AB:", par_AB_name, "\n")
cat("Par AC:", par_AC_name, "\n")
cat("Par BC:", par_BC_name, "\n")

# Extraer genes en interseccion
genes_AB <- lista_intersecciones_vc[[par_AB_name]]$genes_ensembl
genes_AC <- lista_intersecciones_vc[[par_AC_name]]$genes_ensembl
genes_BC <- lista_intersecciones_vc_vc[[par_BC_name]]$genes_ensembl

genes_interseccion_3 <- Reduce(intersect, list(genes_AB, genes_AC, genes_BC))
cat("Genes en la intersección de los tres clusters:", length(genes_interseccion_3), "\n")


writeLines(genes_interseccion_3, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/genes_interseccion_3clusters")

# -----------------------------------------------------------------------
