
# RED PARA LAS INTERSECCIONES DE CLUSTERS -- MAPEAR CON INTERACTOMA

# INTERSECCION MICROGWAS + IMMUNE

# INPUTS

edge_interactoma <- readRDS("./Data/nasertic/input/Combined_STRING40_OTAR0924_FILTER.rds")
interactoma_edges_withSource <- readRDS ("./Data/nasertic/input/interactome_withSource_nodeFilter.rds")
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

lista_intersecciones_im <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")

# PALETAS DE COLORES

# ---- Paleta de nodos ---
colores_nodos <- c(
  "Solo_C1"  = "#E8593C",   # rojo    → semilla solo en trait 1
  "Solo_C2"  = "#2ECC71",   # verde   → semilla solo en trait 2
  "Ambos"    = "#9B59B6",   # morado  → semilla en los dos
  "Propagado"= "#3B8BD4"    # azul    → gen propagado
)

# INTERACTOMA  ----

# Interactoma solo Intact

interactoma_intact_filtrado <- as.data.frame(interactoma_edges_withSource) %>%   
  dplyr::filter(sourceDatabase == "intact") %>%
  dplyr::rename(ENSG_A = targetA, ENSG_B = targetB) %>%
  dplyr::filter(ENSG_A != ENSG_B) %>%
  dplyr::mutate(
    par_id = paste(pmin(ENSG_A, ENSG_B), pmax(ENSG_A, ENSG_B), sep = "_")
  ) %>%
  dplyr::distinct(par_id, .keep_all = TRUE) %>%
  dplyr::select(ENSG_A, ENSG_B)

message("Interacciones IntAct únicas: ", nrow(interactoma_intact_filtrado))

# -------------------------------------------------------------------
# BUCLE CON INTERACTOMA SOLO INTACT ----

carpeta_redes <- "./Output/Piloto_Microbiota/Redes_Interseccion_Intact_IM_MicroGWAS"
dir.create(carpeta_redes, showWarnings = FALSE)

for (nombre_par in names(lista_intersecciones_im)) {
  
  par <- lista_intersecciones_im[[nombre_par]]
  
  c1      <- par$c1
  c2      <- par$c2
  jaccard <- par$jaccard
  
  message("=== Par: ", nombre_par, " (Jaccard = ", round(jaccard, 3), ") ===")
  
  # 3.1 Extraer la tabla de genes que me has enseñado
  tabla_genes <- par$genes_tabla
  
  # Limpiar Ensembl IDs
  genes_interseccion <- unique(gsub(";.*$", "", tabla_genes$gene.ENSG))
  
  if (length(genes_interseccion) == 0) {
    message("  ⚠️ Intersección vacía, saltando...")
    next
  }
  
  message("  Genes en intersección: ", length(genes_interseccion))
  
  # 3.2 Filtrar aristas del interactoma solo para estos genes
  aristas <- interactoma_intact_filtrado %>%
    dplyr::filter(ENSG_A %in% genes_interseccion & ENSG_B %in% genes_interseccion)
  
  # 3.3 Construir la red
  if (nrow(aristas) > 0) {
    red <- graph_from_data_frame(
      d        = aristas[, c("ENSG_A", "ENSG_B")],
      directed = FALSE,
      vertices = data.frame(name = genes_interseccion)
    )
  } else {
    red <- make_empty_graph(n = 0, directed = FALSE) %>%
      add_vertices(length(genes_interseccion), name = genes_interseccion)
    message("  ⚠️ Sin aristas entre genes — red solo con nodos")
  }
  
  # 1. Preparamos una tablita limpia solo con los ENSG y el nombre del gen
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
  
  # 2. Buscamos cada nodo de la red en esa tablita
  idx <- match(V(red)$name, anotacion$ENSG_limpio)
  V(red)$gene_name     <- anotacion$gene.gene[idx]
  V(red)$categoria_nodo <- anotacion$categoria_nodo[idx]
  
  # ---- Color de nodos ----
  V(red)$color <- colores_nodos[V(red)$categoria_nodo]
  
  # Exportar CSV
  
  tabla_nodos <- data.frame(
    ENSG          = V(red)$name,
    gen           = V(red)$gene_name,
    categoria_nodo = V(red)$categoria_nodo,
    senal_inicial  = ifelse(V(red)$categoria_nodo != "Propagado", "Si", "No"),
    grado         = degree(red)
  ) %>% dplyr::arrange(desc(grado))
  
  write.csv2(tabla_nodos,
             file.path(carpeta_redes, paste0("Nodos_", nombre_par, ".csv")),
             row.names = FALSE)
  
  if (nrow(aristas) > 0) {
    write.csv2(aristas,
               file.path(carpeta_redes, paste0("Aristas_", nombre_par, ".csv")),
               row.names = FALSE)
  }
  message("  ✔️ CSVs (Nodos y Aristas) exportados.")
  
  
  # 3.6 Extraer top términos GO
  nombre_par_directo <- paste0(c1, "_vs_", c2)
  nombre_par_inverso <- paste0(c2, "_vs_", c1)
  
  ruta_go_directa <- file.path(carpeta_go_pleiotropia, paste0("GO_Pleiotropia_", nombre_par_directo, ".csv"))
  ruta_go_inversa <- file.path(carpeta_go_pleiotropia, paste0("GO_Pleiotropia_", nombre_par_inverso, ".csv"))
  
  # R comprueba cuál de los dos existe en tu ordenador
  if (file.exists(ruta_go_directa)) {
    ruta_go <- ruta_go_directa
  } else if (file.exists(ruta_go_inversa)) {
    ruta_go <- ruta_go_inversa
  } else {
    ruta_go <- NA
  }
  
  # Si encontró alguna de las dos, extrae los términos
  if (!is.na(ruta_go)) {
    message("  ✔️ GO encontrado: ", basename(ruta_go))
    top_go <- read.csv2(ruta_go) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::pull(Description) %>%
      paste(collapse = "\n")
  } else {
    message("  ⚠️ GO NO encontrado ni al derecho ni al revés.")
    top_go <- "Sin términos GO significativos"
  }
  
  
  # 3.7 DIBUJAR Y EXPORTAR PDF
  titulo <- paste0(par$trait_1, " vs ", par$trait_2, "  (Jaccard = ", round(jaccard, 3), ")")
  
  # Layout
  if(ecount(red) > 0) {
    e <- as_edgelist(red, names = F)
    layout_red <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(red),
                                                    area = 4 * (vcount(red)^2),
                                                    repulse.rad = vcount(red)^3)
  } else {
    layout_red <- layout_in_circle(red) 
  }
  
  # 4.9 Exportar PDF
  titulo   <- paste0(par$trait_1, " vs ", par$trait_2, "  (Jaccard = ", round(jaccard, 3), ")")
  ruta_pdf <- file.path(carpeta_redes, paste0("Red_", nombre_par, ".pdf"))
  
  pdf(ruta_pdf, width = 12, height = 10)
  par(mar = c(8, 1, 3, 1))
  
  plot(red,
       layout             = layout_red,
       vertex.label       = V(red)$gene_name,
       vertex.label.cex   = 0.4,
       vertex.label.color = "black",
       vertex.size        = 6,
       vertex.color       = V(red)$color,
       vertex.frame.color = "white",
       edge.color         = "#E69F00",
       edge.width         = 0.6,
       main               = titulo)
  
  mtext(paste0(c1, "  vs  ", c2), side = 3, line = 0.3, cex = 0.65, col = "gray50")
  mtext(paste0("Top GOBP: \n", top_go), side = 1, line = 5, cex = 0.75, col = "gray30")
  
  cats_presentes <- intersect(names(colores_nodos), unique(V(red)$categoria_nodo))
  etiquetas_nodos <- c(
    "Solo_C1"   = paste0("Semilla: ", par$trait_1),
    "Solo_C2"   = paste0("Semilla: ", par$trait_2),
    "Ambos"     = "Semilla: ambos traits",
    "Propagado" = "Gen propagado"
  )
  
  legend("bottomleft",
         legend = etiquetas_nodos[cats_presentes],
         fill   = colores_nodos[cats_presentes],
         bty    = "n", cex = 0.8,
         title = "Leyenda")
  
  dev.off()
  
  message("  ✔️ PDF exportado | Nodos: ", vcount(red), " | Aristas: ", ecount(red))
}

message("\n¡Proceso completado! Todas las redes generadas.")
# -------------------------------------------------------------------
