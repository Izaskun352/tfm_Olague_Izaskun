
source("scripts/00_setup.R")

#------------------------------------------------------------------------------------------------------
# SKIN x HAPLOGOBINA --> CILLIUM
#------------------------------------------------------------------------------------------------------

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
lista_intersecciones_im  <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")
lista_intersecciones_vc  <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_VarComun_MicroGWAS.rds")
lista_intersecciones_vc_vc  <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_VarComun_VarComun.rds")
all_pairs <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Todos_Los_Pares_Jaccard.rds")

interactoma_edges_withSource <- readRDS ("./Data/nasertic/input/interactome_withSource_nodeFilter.rds")
interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")

variantes_raras <- readRDS("./Data/Variantes/rare_filter_Select_Col.rds")
variantes_raras_codificantes <- readRDS("./Data/Variantes/variantes_raras_codificantes.rds")

# OUTPUTS
genes_all_clusters <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes_all_clusters.rds")
genes_all <- readLines( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes_all.txt")
genes_all_interactoma_intact <- readLines( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_intact.txt")

genes_all_pathway <- readLines( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_pathway.txt")
genes_all_intact_pathway <- readLines( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_intact_pathway.txt")

genes_all_cluster_string <- readRDS( file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes_all_clusters_sin_string.rds")
genes_all_interactoma_string <- readLines( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_string.txt")
genes_all_string_pathway <- readLines( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_string_pathway.txt")

# CaRGAR LO COMÚN

load("Data/Redes/config_red.RData")               # colores_db, prioridad_db, add_layer
aristas_color <- readRDS("Data/Redes/aristas_color.rds")
source("Scripts/Scripts_Limpios/Red_Capas_Finales/01_Funciones_Redes.R")

#------------------------------------------------------------------------------------------------------
# ESTUDIAR SOLAPAMIENTO
#------------------------------------------------------------------------------------------------------

skin_haplogobina_name <- names(lista_intersecciones_vc)[
  sapply(lista_intersecciones_vc, function(x) {
    grepl("skin microbiome", x$trait_microGWAS, ignore.case = TRUE) &
      grepl("adiponectin",   x$trait_comun,      ignore.case = TRUE)
  })
]
print(skin_haplogobina_name)
skin_haplogobina <- lista_intersecciones_vc[[skin_haplogobina_name]]

# Construir grafo base

g <- build_network(skin_haplogobina, aristas_color, colores_db)

# Añadimos seed genes
g <- add_layer_seed(g, skin_haplogobina$genes_tabla)

.#--------------------------------------------------------
# AÑADIMOS CAPAS / ATRIBUTOS
#--------------------------------------------------------
# ---- GENES CON VARIANTES RARAS DE INTERES ----

# DF con variantes raras en lugares codificantes
variantes_raras <- readRDS("./Data/Variantes/rare_filter_Select_Col.rds")
variantes_raras_codificantes <- readRDS("./Data/Variantes/variantes_raras_codificantes.rds")
g <- add_layer_rare_variants(g, variantes_raras_codificantes)

df_nodos <- igraph::as_data_frame(g, what = "vertices")

df_enfermedades <- df_nodos %>%
  filter(has_rare_variant == TRUE) %>%
  select(name, symbol, diseases_rare_variant)  

df_enfermedades_limpio <- df_enfermedades %>%
  tidyr::separate_rows(diseases_rare_variant, sep = ";") %>%
  left_join(traits_MicroGWAS_areas, by = c("diseases_rare_variant" = "Rasgo")) %>%
  select(
    name_ensembl = name.x,
    symbol,
    disease_id = diseases_rare_variant,
    disease_name = name.y,
    Nombre_area,
    therapeuticAreas
  )

ids_sin_anotar <- df_enfermedades_limpio %>%
  filter(is.na(disease_name)) %>%
  pull(disease_id) %>%
  unique()





#----------------------------------------------------------
# PLOT
#----------------------------------------------------------

# filtrar por score de string
string_score_min <- 0.4
g_plot <- igraph::delete_edges(g, which(
  E(g)$sourceDatabase == "string" & 
    (is.na(E(g)$scoring) | E(g)$scoring < string_score_min)
))

g_plot <- igraph::delete_vertices(g_plot, which(igraph::degree(g_plot) == 0))

# Asignar colores si tiene variante rara o no
V(g_plot)$color_nodo <- ifelse(V(g_plot)$has_rare_variant, "slateblue", "grey85")


set.seed(42)
plot <- ggraph(g_plot, layout = "fr") +
  geom_edge_link(color = "grey70", alpha = 0.6, width = 0.8) +
  
  geom_node_point(aes(color = color_nodo), size = 8) +
  scale_color_identity(
    name   = "Rare Variant",
    labels = c("slateblue" = "Associated", "grey85" = "None"),
    breaks = c("slateblue", "grey85"),
    guide  = guide_legend(override.aes = list(size = 5))
  ) +
  
  geom_node_text(
    aes(label = symbol), 
    repel = TRUE, 
    size = 3.5, 
    fontface = "bold",
    color = "black"
  ) +
  
  # Estética final
  ggtitle(paste("Gene Network (STRING score >", string_score_min, ")")) +
  theme_graph() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold")
  )

cairo_pdf("./Output/Redes_Capas/Skin_Haplogobina.pdf", width = 25, height = 15)
print(plot)
dev.off()


# Estudiar atributos nodos

df_atr_nodos <- as_data_frame(g, what = "vertices")
print(df_atr_nodos["ENSG00000157423",])


#---------------------------------------------------------
# MEGACLUSTER CON TODOS LOS TRAITS VARCOMUN
#---------------------------------------------------------

# hacemos red con los solapamientos
library(ggraph)
buscar_par <- function(lista, cA, cB) {
  nombre <- names(lista)[
    sapply(lista, function(x) {
      (x$c1 == cA & x$c2 == cB) | (x$c1 == cB & x$c2 == cA)
    })
  ]
  if (length(nombre) == 0) stop(paste("Par no encontrado:", cA, "vs", cB))
  nombre
}
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
  
  return(list(plot = plot_red, red = red))
}

cluster_A <- "EFO_0801228_Cluster_2"
cluster_B <- "ZSCO.EFO_0004640_Cluster_1"
cluster_C <- "ZSCO.EFO_0004884_Cluster_2"
cluster_D <- "ZSCO.EFO_0004502_Cluster_1"

# ── 1. Buscar pares en su lista correcta ─────────────────────────────────────
par_AB_name <- buscar_par(lista_intersecciones_vc,    cluster_A, cluster_B)
par_AC_name <- buscar_par(lista_intersecciones_vc,    cluster_A, cluster_C)
par_AD_name <- buscar_par(lista_intersecciones_vc,    cluster_A, cluster_D)
par_BC_name <- buscar_par(lista_intersecciones_vc_vc, cluster_B, cluster_C)
#par_BD_name <- buscar_par(lista_intersecciones_vc_vc, cluster_B, cluster_D)
par_CD_name <- buscar_par(lista_intersecciones_vc_vc, cluster_C, cluster_D)

# ── 2. Extraer genes de su lista correcta ────────────────────────────────────
genes_AB <- lista_intersecciones_vc[[par_AB_name]]$genes_ensembl
genes_AC <- lista_intersecciones_vc[[par_AC_name]]$genes_ensembl
genes_AD <- lista_intersecciones_vc[[par_AD_name]]$genes_ensembl
genes_BC <- lista_intersecciones_vc_vc[[par_BC_name]]$genes_ensembl
#genes_BD <- lista_intersecciones_vc_vc[[par_BD_name]]$genes_ensembl
genes_CD <- lista_intersecciones_vc_vc[[par_CD_name]]$genes_ensembl

# ── 3. Intersección de los 4 clusters ────────────────────────────────────────
genes_interseccion_4 <- Reduce(intersect, list(genes_AB, genes_AC, genes_AD,
                                               genes_BC,  genes_CD))
cat("Genes en la intersección de los 4 clusters:", length(genes_interseccion_4), "\n")

genes_extra_ensembl <- c("ENSG00000197653", "ENSG00000188596")  # <- tus genes
genes_extra_simbolo <- c("DNAH10", "CFAP54")                         # <- sus símbolos

filas_extra <- tibble::tibble(
  gene.ENSG        = genes_extra_ensembl,
  gene.gene        = genes_extra_simbolo,
  es_semilla_en_c1 = TRUE,
  es_semilla_en_c2 = FALSE
)

# ── 6. Construir genes_tabla de referencia ────────────────────────────────────
tabla_ref <- lista_intersecciones_vc[[par_AB_name]]$genes_tabla %>%
  dplyr::mutate(ENSG_limpio = gsub(";.*$", "", gene.ENSG)) %>%
  dplyr::filter(ENSG_limpio %in% genes_interseccion_4) %>%
  dplyr::bind_rows(filas_extra)

# ── 7. Construir entry sintético ──────────────────────────────────────────────
entry_sintetico <- list(
  c1              = cluster_A,
  c2              = paste(cluster_B, cluster_C, cluster_D, sep = " & "),
  jaccard         = NA,
  trait_microGWAS = lista_intersecciones_vc[[par_AB_name]]$trait_microGWAS,
  trait_comun     = paste(
    lista_intersecciones_vc[[par_AB_name]]$trait_comun,
    lista_intersecciones_vc[[par_AC_name]]$trait_comun,
    lista_intersecciones_vc[[par_AD_name]]$trait_comun,
    sep = " / "
  ),
  genes_ensembl   = c(genes_interseccion_4, genes_extra_ensembl),
  genes_simbolo   = c(
    lista_intersecciones_vc[[par_AB_name]]$genes_simbolo[
      lista_intersecciones_vc[[par_AB_name]]$genes_ensembl %in% genes_interseccion_4
    ],
    genes_extra_simbolo
  ),
  genes_tabla     = tabla_ref
)


nombre_sintetico <- paste(cluster_A, cluster_B, cluster_C, cluster_D, sep = "_vs_")
lista_sintetica  <- setNames(list(entry_sintetico), nombre_sintetico)

# ── 5. Generar y guardar la red ───────────────────────────────────────────────
resultado <- generar_red_par(
  par_nombre           = nombre_sintetico,
  interactoma          = interactoma_edges_withSource,
  lista_intersecciones = lista_sintetica,
  #excluir_db           = "string"   # opcional
)
plot <- resultado$plot
red <- resultado$red
saveRDS(resultado$red, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/red_todos_interseccion.rds")
pdf(file.path( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/red_todos_interseccion.pdf"), width = 12, height = 12)
print(plot)
dev.off()

gen <- "DNAH11"
V(red)$gene_name[neighbors(red, v = which(V(red)$gene_name == gen))]

#sacamos todos los genes----
genes_skin <- as.data.frame(read.csv2("./Output/Piloto_Microbiota/Clusters_MicroGWAS/EFO_0801228_Cluster_2.csv")) %>%
  select(ENSG = gene.ENSG, p.adj = gene.padj, gene = gene.gene) %>%
  mutate(cluster = "skin") %>%
  mutate(n_seed_genes = sum(p.adj == 100, na.rm = TRUE))
genes_haptoglobin <- as.data.frame(read.csv2("./Output/Piloto_Microbiota/Clusters_Traits_VarComun/ZSCO.EFO_0004640_Cluster_1.csv")) %>%
  select(ENSG = gene.ENSG, p.adj = gene.padj, gene = gene.gene) %>%
  mutate(cluster = "haptoglobin") %>%
  mutate(n_seed_genes = sum(p.adj == 100, na.rm = TRUE))
genes_breast <- as.data.frame(read.csv2("./Output/Piloto_Microbiota/Clusters_Traits_VarComun/ZSCO.EFO_0004884_Cluster_2.csv")) %>%
  select(ENSG = gene.ENSG, p.adj = gene.padj, gene = gene.gene) %>%
  mutate(cluster = "breast size") %>%
  mutate(n_seed_genes = sum(p.adj == 100, na.rm = TRUE))
genes_adiponectin <- as.data.frame(read.csv2("./Output/Piloto_Microbiota/Clusters_Traits_VarComun/ZSCO.EFO_0004502_Cluster_1.csv")) %>%
  select(ENSG = gene.ENSG, p.adj = gene.padj, gene = gene.gene) %>%
  mutate(cluster = "adiponectin") %>%
  mutate(n_seed_genes = sum(p.adj == 100, na.rm = TRUE))

genes_especificos <- Reduce(intersect,list(genes_skin$gene, genes_adiponectin$gene, genes_breast$gene, genes_haptoglobin$gene))
genes_all <- unique(c(genes_skin$ENSG, genes_adiponectin$ENSG, genes_breast$ENSG, genes_haptoglobin$ENSG))
genes_all_cluster <- bind_rows(genes_skin, genes_adiponectin, genes_breast, genes_haptoglobin)
writeLines(genes_all, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes_all.txt")

# añadimos los genes directamente unidos del interactoma

#solo IntAct
interactoma_intact <- as.data.frame(interactoma_edges_withSource) %>%   
  dplyr::filter(sourceDatabase == "intact") %>%
  dplyr::rename(ENSG_A = targetA, ENSG_B = targetB) %>%
  dplyr::filter(ENSG_A != ENSG_B) %>%
  dplyr::mutate(
    par_id = paste(pmin(ENSG_A, ENSG_B), pmax(ENSG_A, ENSG_B), sep = "_")
  ) %>%
  dplyr::distinct(par_id, .keep_all = TRUE) %>%
  dplyr::select(ENSG_A, ENSG_B)

genes_all_cluster <- genes_all_cluster %>%
  left_join(
    bind_rows(
      data.frame(ENSG = interactoma_intact[,1], interactor = interactoma_intact[,2]),
      data.frame(ENSG = interactoma_intact[,2], interactor = interactoma_intact[,1])
    ) %>%
      group_by(ENSG) %>%
      summarise(genes_unidos = paste(unique(interactor), collapse = ", "), .groups = 'drop'),
    by = "ENSG"
  )
saveRDS(genes_all_cluster, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes_all_clusters.rds")

genes_all_interactoma_intact <- unique(c(
  genes_all, 
  unique(unlist(strsplit(na.omit(genes_all_cluster$genes_unidos), ", ")))
))
writeLines(genes_all_interactoma_intact, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_intact.txt")

  #sacamos los genes de DAVID del pathway KEGG: Motor proteins
library(KEGGREST)

info_pathway <- keggGet("hsa04814")

genes_pathway <- info_pathway[[1]]$GENE
#entrez_ids <- genes_brutos[c(TRUE, FALSE)]
genes_pathway   <- as.data.frame(genes_pathway[c(FALSE, TRUE)]) %>%
  rename(genes = 1) %>%
  mutate(gene_limpio = sub(";.*", "", genes)) %>%
  pull(gene_limpio)

# interseccion con genes del pathway Motor proteins ----

genes_all_pathway <- as.data.frame(interactoma) %>%
  filter(gene %in% genes_pathway) %>%
  pull(ENSG) %>%
  na.omit() %>%
  intersect(genes_all)
writeLines(genes_all_pathway, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_pathway.txt")
  
genes_all_intact_pathway <- as.data.frame(interactoma) %>%
  filter(gene %in% genes_pathway) %>%
  pull(ENSG) %>%
  na.omit() %>%                                                                                         
  intersect(genes_all_interactoma_intact)
writeLines(genes_all_interactoma_pathway, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_intact_pathway.txt")

#sin string
interactoma_string <- as.data.frame(interactoma_edges_withSource) %>%   
  dplyr::filter(sourceDatabase != "string") %>%
  dplyr::rename(ENSG_A = targetA, ENSG_B = targetB) %>%
  dplyr::filter(ENSG_A != ENSG_B) %>%
  dplyr::mutate(
    par_id = paste(pmin(ENSG_A, ENSG_B), pmax(ENSG_A, ENSG_B), sep = "_")
  ) %>%
  dplyr::distinct(par_id, .keep_all = TRUE) %>%
  dplyr::select(ENSG_A, ENSG_B)

genes_all_cluster_string <- genes_all_cluster %>%
  left_join(
    bind_rows(
      data.frame(ENSG = interactoma_string[,1], interactor = interactoma_string[,2]),
      data.frame(ENSG = interactoma_string[,2], interactor = interactoma_string[,1])
    ) %>%
      group_by(ENSG) %>%
      summarise(genes_unidos = paste(unique(interactor), collapse = ", "), .groups = 'drop'),
    by = "ENSG"
  )
saveRDS(genes_all_cluster_string, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes_all_clusters_sin_string.rds")

genes_all_interactoma_string <- unique(c(
  genes_all, 
  unique(unlist(strsplit(na.omit(genes_all_cluster_string$genes_unidos), ", ")))
))
writeLines(genes_all_interactoma_string, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_string.txt")

genes_all_string_pathway <- as.data.frame(interactoma) %>%
  filter(gene %in% genes_pathway) %>%
  pull(ENSG) %>%
  na.omit() %>%                                                                                         
  intersect(genes_all_interactoma_string)
writeLines(genes_all_string_pathway, "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/genes-all_string_pathway.txt")

# vemos si son DEG y en qué enfermedades----

DEG_atlas <- open_dataset("./Data/Diccionarios/Diff_expr_OTAR") %>%
  collect()
lista_genes_DEG_np <- readRDS("./Output/Piloto_Microbiota/Expresion_Diferencial_OTAR/df_genes_NP.rds")
genes_DEG <- intersect(genes_all_interactoma_pathway, DEG_atlas$targetId)

diseases_deg <- DEG_atlas %>%
  filter(targetId %in% genes_DEG) %>%
  group_by(diseaseId, contrast) %>%
  
  summarise(
    num_genes = n_distinct(targetFromSourceId),
    genes_implicados = paste(unique(targetFromSourceId), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(num_genes)) %>%
  left_join(traits_MicroGWAS_areas, by = c("diseaseId" = "Rasgo")) %>%
  relocate(name, diseaseId, num_genes, contrast, Nombre_area, therapeuticAreas)

#vemos cuales son los genes comunes por área (seleccionamos enfermedades con al menos 10 genes asociados)
diccionario_genes <- setNames(interactoma[,3], interactoma[,1])
areas_con_genes_comunes <- diseases_deg %>%
  filter(num_genes >= 10, !is.na(therapeuticAreas)) %>%
  
  separate_rows(therapeuticAreas, sep = ",\\s*") %>%
  
  group_by(therapeuticAreas) %>%
  
  summarise(
    num_contrastes = n_distinct(contrast),
    contrastes_agrupados = paste(unique(contrast), collapse = " | "),
    nombres_areas = paste(unique(na.omit(Nombre_area)), collapse = " | "),
    genes_lista = list(Reduce(intersect, strsplit(genes_implicados, ",\\s*"))),
    .groups = "drop"
  ) %>%
  
  filter(num_contrastes > 1) %>%
  
  mutate(
    num_genes_comunes = lengths(genes_lista),
    genes_comunes_ENSG = sapply(genes_lista, paste, collapse = ", "),
    
    simbolos_comunes = sapply(genes_lista, function(x) {
      simbolos <- diccionario_genes[x]
      simbolos[is.na(simbolos)] <- x[is.na(simbolos)]
      paste(simbolos, collapse = ", ")
    })
  ) %>%
  
  filter(num_genes_comunes > 0) %>%

  select(-genes_lista) %>%
  relocate(nombres_areas, .after = therapeuticAreas) %>%
  relocate(simbolos_comunes, .after = genes_comunes_ENSG) %>%
  arrange(desc(num_genes_comunes))


# variantes raras ----

variantes_raras_codificantes <- readRDS("./Data/Variantes/variantes_raras_codificantes.rds")
all_genes_rare <- as.data.frame(variantes_raras_codificantes) %>%
  dplyr::filter(targetFromSourceId %in% genes_all_pathway)

resumen_all_genes_rare <- all_genes_rare %>%
  group_by(targetFromSourceId) %>%
  summarise(
    num_variantes_distintas = n_distinct(variantId),
    enfermedades_asociadas = paste(unique(diseaseFromSourceMappedId), collapse = " | "),
    .groups = "drop"
  ) %>%
  left_join(
    as.data.frame(interactoma) %>% select(ENSG, simbolo_gen = gene), 
    by = c("targetFromSourceId" = "ENSG")
  ) %>%
  relocate(simbolo_gen, .after = targetFromSourceId) %>%
  arrange(desc(num_variantes_distintas))
saveRDS(all_genes_rare, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/all_genes_rare.rds")
saveRDS(all_genes_rare, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/resumen_all_genes_rare.rds")


all_genes_interactoma_rare <- as.data.frame(variantes_raras_codificantes) %>%
  dplyr::filter(targetFromSourceId %in% genes_all_interactoma_pathway)
resumen_all_genes_rare_interactoma <- all_genes_interactoma_rare %>%
  group_by(targetFromSourceId) %>%
  summarise(
    num_variantes_distintas = n_distinct(variantId),
    enfermedades_asociadas = paste(unique(diseaseFromSourceMappedId), collapse = " | "),
    .groups = "drop"
  ) %>%
  left_join(
    as.data.frame(interactoma) %>% select(ENSG, simbolo_gen = gene), 
    by = c("targetFromSourceId" = "ENSG")
  ) %>%
  relocate(simbolo_gen, .after = targetFromSourceId) %>%
  arrange(desc(num_variantes_distintas))
saveRDS(all_genes_interactoma_rare, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/all_genes_interactoma_rare.rds")
saveRDS(resumen_all_genes_rare_interactoma, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/resumen_all_genes_interactoma_rare.rds")

# filtramos con los genes que forman parte del complejo de la dineina
library(stringr)
genes_dineina <- genes_all_clusters %>%
  filter(str_detect(gene, "^DYNC|^DNAH|^DNAI|^DNAL|^DNAAF")) %>%
  select(ENSG, gene) %>%
  distinct()

variantes_dineina <- as.data.frame(variantes_raras) %>%
  inner_join(genes_dineina, by = c("targetFromSourceId" = "ENSG"))

variantes_dineina_codif <- as.data.frame(variantes_raras_codificantes) %>%
  inner_join(genes_dineina, by = c("targetFromSourceId" = "ENSG"))

#cogemos las variantes de ProtVar----
zips <- list.files("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/variantes_protVar", 
                   pattern = "\\.zip$", 
                   full.names = TRUE)

df_total <- map_dfr(zips, function(zip) {
  csv_name <- unzip(zip, list = TRUE)$Name[1]
  read_csv(unz(zip, csv_name), 
           show_col_types = FALSE,
           col_types = cols(.default = "c"))  
})

df_clean <- df_total |>
  select(
    Gene,
    Chromosome,
    Coordinate,
    Reference_allele,
    Alternative_allele,
    Amino_acid_position,
    Amino_acid_change,
    Consequences,
    CADD_phred_like_score,
    `Uniprot_canonical_isoform_(non_canonical)`,
    Protein_name
  ) |>
  # Convertir numéricas
  mutate(
    Amino_acid_position  = as.numeric(Amino_acid_position),
    Chromosome           = as.numeric(Chromosome),
    Coordinate           = as.numeric(Coordinate),
    CADD_phred_like_score = as.numeric(CADD_phred_like_score)

  ) |>
  # Convertir "N/A" en texto a NA real
  mutate(across(where(is.character), ~ na_if(., "N/A")))

df_clean <- df_clean |>
  filter(
    Consequences == "stop gained" |
      (Consequences == "missense" & CADD_phred_like_score >= 20)
  )
saveRDS(df_clean, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/ProtVar_clean.rds")
#dominios dynein proteins----
dominios <- read_csv("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/dominios_genes_dynein.csv")

df <- read_csv("C:/Users/Usuario/Downloads/46dc8739-1530-4c96-a251-1306a0310ccb.csv/46dc8739-1530-4c96-a251-1306a0310ccb.csv")

#cruzamos con df de dominios
df_variantes_dominios <- df_clean |>
  mutate(Amino_acid_position = as.numeric(Amino_acid_position)) |>
  inner_join(dominios, by = c("Gene" = "gene"),
             relationship = "many-to-many") |>
  filter(Amino_acid_position >= inicio & Amino_acid_position <= fin)

#boxplot variantes en dominios vs variantes fuera de dominios----
df_variantes_out <- df_clean |>
  mutate(Amino_acid_position = as.numeric(Amino_acid_position)) |>
  anti_join(df_variantes_dominios, 
            by = c("Gene", "Chromosome", "Coordinate", 
                   "Reference_allele", "Alternative_allele"))

#plot
datos_ggplot <- bind_rows( 
  df_variantes_dominios |>
    distinct(Gene, Chromosome, Coordinate, Reference_allele, Alternative_allele, 
             CADD_phred_like_score) |>
    mutate(grupo = "En dominio"),
  
  df_variantes_out |>
    distinct(Gene, Chromosome, Coordinate, Reference_allele, Alternative_allele,
             CADD_phred_like_score) |>
    mutate(grupo = "Fuera de dominio")
) |>
  mutate(CADD_phred_like_score = as.numeric(CADD_phred_like_score))

ggplot(datos_ggplot, aes(x = grupo, y = CADD_phred_like_score, fill = grupo, color = grupo)) +
  geom_violin(trim = FALSE,alpha = 0.3, width = 0.6) +          
  geom_boxplot(width = 0.1, outlier.shape = NA,    
               alpha = 0.5, median.linewidth = 2) +
  
  stat_summary(fun = median, geom = "text",
               aes(label = round(after_stat(y), 1)),
               vjust = -0.5, hjust = -0.8, size = 4, fontface = "bold") +
  
  scale_fill_manual(values = c( "#FDBF6F", "#C1E1C1")) +
  scale_color_manual(values = c("tan4", "#2F4F4F")) +
  theme_minimal() +
  labs(
    title = "Distribución de variantes", x = "", y = " ") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(face = "bold")
  )

ggsave("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/Boxplot_variantes_dominio.pdf", width = 8, height = 6)

#---------------------------------------------------
#SELECCIONAR VARIANTES EN FUNCION DE LA ENFERMEDAD
#---------------------------------------------------

# cruzamos con df variantes raras --> estudiar que enfermedad tienen asociada
#mapear con disease
diseases <- readRDS("./Data/nasertic/input/All_diseases.rds")

variantes_dineina <- variantes_dineina %>%
  left_join(as.data.frame(diseases), by = c("diseaseFromSourceMappedId" = "ID"))

diseases_2 <- tibble(
  diseaseFromSourceMappedId = c(
    "Orphanet_93270", "Orphanet_93269", "MONDO_0015522", 
    "Orphanet_85173", "Orphanet_474", "Orphanet_388", "MONDO_0007462"),
  name_rescate = c(
    "Short rib-polydactyly syndrome, Saldino-Noonan type",
    "Short rib-polydactyly syndrome, Majewski type",
    "Nemaline myopathy 3",
    "IMAGe syndrome",
    "Jeune syndrome",
    "Hirschsprung disease",
    "Multiple sclerosis, susceptibility to"))

variantes_dineina <- variantes_dineina %>%
  left_join(diseases_2, by = "diseaseFromSourceMappedId") %>%
  mutate(name = coalesce(name, name_rescate)) %>%
  select(-name_rescate)

#cruzamos con el df de ProtVar

#filtramos por CADD >= 30

df_CADD_30 <- df_total |>
  select(
    Gene,
    Chromosome,
    Coordinate,
    Reference_allele,
    Alternative_allele,
    Amino_acid_position,
    Amino_acid_change,
    Consequences,
    CADD_phred_like_score,
    `Uniprot_canonical_isoform_(non_canonical)`,
    Protein_name
  ) |>
  # Convertir numéricas
  mutate(
    Amino_acid_position  = as.numeric(Amino_acid_position),
    Chromosome           = as.numeric(Chromosome),
    Coordinate           = as.numeric(Coordinate),
    CADD_phred_like_score = as.numeric(CADD_phred_like_score)
    
  ) |>
  mutate(across(where(is.character), ~ na_if(., "N/A"))) |>   # Convertir "N/A" en texto a NA real
  filter(Consequences == "missense" & CADD_phred_like_score >= 30)

#crear columnas para cruzar --> cromosoma + coordinate + reference_allele + alternative_allele
variantes_dineina <- variantes_dineina |>
  separate(variantId, into = c("chr_ens", "coord_ens", "ref_ens", "alt_ens"),
           sep = "_", remove = FALSE) |>
  mutate(chr_ens   = as.numeric(chr_ens),
         coord_ens = as.numeric(coord_ens))
saveRDS(variantes_dineina, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/variantes_dineina.rds")
#cruzar
df_CADD_30_enf <- df_CADD_30 |>
  left_join(
    variantes_dineina |> select(chr_ens, coord_ens, ref_ens, alt_ens, name, diseaseFromSourceMappedId, therapeuticAreas),
    by = c(
      "Chromosome"        = "chr_ens",
      "Coordinate"        = "coord_ens",
      "Reference_allele"  = "ref_ens",
      "Alternative_allele"= "alt_ens"
    )
  )  #hay 2295 variantes con enfermedad asociada
saveRDS(df_CADD_30_enf_filtrado, file = "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/ProtVar_dineina_CADD30.rds")

enfermedades_variantes <- df_CADD30_dominios |>
  filter(!is.na(name)) |>
  group_by(name) |>
  summarise(
    n_variantes = n(),
    genes       = paste(unique(Gene), collapse = ", ")
  ) |>
  arrange(desc(n_variantes))
df_CADD_30_enf_filtrado <- df_CADD_30_enf |>
  filter(!is.na(name))

df_CADD30_dominios <- df_CADD_30_enf |>
  filter(!is.na(name)) |>
  inner_join(dominios, by = c("Gene" = "gene"),
             relationship = "many-to-many") |>
  filter(Amino_acid_position >= inicio & Amino_acid_position <= fin)


#agrupar variantes en menos de 20aa + mismo gen
df_CADD30_clusters <- df_CADD_30_enf_filtrado |>
  filter(!is.na(name)) |>
  arrange(Gene, name, Amino_acid_position) |>
  group_by(Gene, name) |>
  mutate(
    dist_anterior = Amino_acid_position - lag(Amino_acid_position, default = first(Amino_acid_position)),
    cluster_id = cumsum(dist_anterior > 20)
  ) |>
  ungroup()
resumen_clusters <- df_CADD |>
  group_by(Gene, name, cluster_id) |>
  summarise(
    n_variantes = n(),
    pos_inicio  = min(Amino_acid_position),
    pos_fin     = max(Amino_acid_position),
    CADD_max    = max(CADD_phred_like_score, na.rm = TRUE),
    CADD_medio  = round(mean(CADD_phred_like_score, na.rm = TRUE), 1),
    .groups     = "drop"
  ) |>
  arrange(desc(n_variantes), desc(CADD_max))

resumen_clusters |>
  group_by(Gene) |>
  summarise(
    n_clusters        = n_distinct(cluster_id),
    n_variantes_total = sum(n_variantes),
    n_enfermedades    = n_distinct(name),
    CADD_max          = max(CADD_max)
  ) |>
  arrange(desc(n_variantes_total))
