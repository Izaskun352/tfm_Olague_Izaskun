
source("scripts/00_setup.R")
# ========================================================
# FUNCIONES PARA CONSTRUIR REDES CON CApAS
# ========================================================

# --- FUNCION PARA AÑADIR CAPAS

add_layer <- function(g, gene_list, attr_name) {
  vertex_attr(g, attr_name) <- V(g)$name %in% gene_list
  g
}


# --- Funcion para construir grafo a partir de los genes

build_network <- function(interseccion,
                          aristas_color,
                          colores_db,
                          keep_isolated  = TRUE,
                          exclude_db     = NULL) {
  
  genes_interseccion <- interseccion$genes_ensembl
  
  aristas_df <- aristas_color %>%
    dplyr::filter(ENSG_A %in% genes_interseccion,
                  ENSG_B %in% genes_interseccion) %>%
    { if (!is.null(exclude_db))               # filtra si se especifica
      dplyr::filter(., !sourceDatabase %in% exclude_db)
      else . } %>%
    dplyr::select(from = ENSG_A, to = ENSG_B,
                  sourceDatabase, dbs_presentes, n_dbs,
                  scoring,
                  n_tejidos_Neural,
                  n_tejidos_Gut_microbiome,
                  n_tejidos_Immune_systemic,
                  n_tejidos_Peripheral)
  
  nodos_df <- data.frame(
    name = if (keep_isolated) genes_interseccion
    else unique(c(aristas_df$from, aristas_df$to))
  )
  
  g <- igraph::graph_from_data_frame(aristas_df, vertices = nodos_df, directed = FALSE)
  
  # --- Atributos de nodo base ---
  E(g)$color <- colores_db[E(g)$sourceDatabase]
  E(g)$width <- scales::rescale(E(g)$n_dbs, to = c(0.5, 3))
  
  V(g)$degree      <- igraph::degree(g)
  V(g)$betweenness <- igraph::betweenness(g, normalized = TRUE)
  V(g)$color       <- "grey80"
  V(g)$size        <- scales::rescale(V(g)$degree, to = c(4, 18))
  
  # --- Mapeo ENSG --> símbolo ---
  idx <- match(V(g)$name, interseccion$genes_tabla$gene.ENSG)
  V(g)$symbol <- interseccion$genes_tabla$gene.gene[idx]
  
  g
}


# -------------------------------------------------------
# AÑADIR FILTROS
# -------------------------------------------------------

#---- Añadir genes semilla ----

add_layer_seed <- function(g, genes_tabla) {
  # genes_tabla puede venir de una o varias intersecciones (con filas repetidas por gen)
  semilla_summary <- genes_tabla %>%
    dplyr::filter(gene.ENSG %in% V(g)$name,
                  es_semilla_en_alguno == TRUE) %>%
    dplyr::mutate(
      tipo_c1 = ifelse(startsWith(Cluster_1, "ZSCO"), "neuro", "microbiome"),
      tipo_c2 = ifelse(startsWith(Cluster_2, "ZSCO"), "neuro", "microbiome")
    ) %>%
    dplyr::mutate(
      seed_type_fila = dplyr::case_when(
        es_semilla_en_c1 & es_semilla_en_c2 & tipo_c1 != tipo_c2 ~ "both",
        es_semilla_en_c1 & es_semilla_en_c2 & tipo_c1 == tipo_c2 ~ tipo_c1,
        es_semilla_en_c1 ~ tipo_c1,
        es_semilla_en_c2 ~ tipo_c2,
        TRUE ~ NA_character_
      ),
      trait_semilla_contribucion = dplyr::case_when(
        es_semilla_en_c1 & es_semilla_en_c2 ~ paste(Trait_Nombre_1, Trait_Nombre_2, sep = ";"),
        es_semilla_en_c1                    ~ Trait_Nombre_1,
        es_semilla_en_c2                    ~ Trait_Nombre_2,
        TRUE                                ~ NA_character_
      )
    ) %>%
    dplyr::group_by(gene.ENSG) %>%
    dplyr::summarise(
      is_seed = TRUE,
      seed_type = dplyr::case_when(
        "both" %in% seed_type_fila ~ "both",
        "microbiome" %in% seed_type_fila & "neuro" %in% seed_type_fila ~ "both",
        "microbiome" %in% seed_type_fila ~ "microbiome",
        "neuro" %in% seed_type_fila ~ "neuro",
        TRUE ~ NA_character_
      ),
      traits_semilla = paste(sort(unique(unlist(strsplit(trait_semilla_contribucion, ";")))), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(name = gene.ENSG)
  
  # Atributo booleano
  g <- add_layer(g, semilla_summary$name, "is_seed")
  
  # Atributos cualitativos
  idx <- match(V(g)$name, semilla_summary$name)
  V(g)$seed_type      <- semilla_summary$seed_type[idx]
  V(g)$traits_semilla <- semilla_summary$traits_semilla[idx]
  
  g
}


#---- Añadir atributo DEG ----

add_layer_DEG <- function(g, lista_genes_DEG) {
  
  deg_summary <- lista_genes_DEG %>%
    dplyr::filter(targetFromSourceId %in% V(g)$name) %>%
    dplyr::group_by(targetFromSourceId) %>%
    dplyr::summarise(
      diseases_DEG = paste(sort(unique(diseaseFromSourceMappedId)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(name = targetFromSourceId)
  
  # Atributo 1: booleano
  g <- add_layer(g, deg_summary$name, "is_DEG")
  
  # Atributo 2: en qué enfermedades
  idx <- match(V(g)$name, deg_summary$name)
  V(g)$diseases_DEG <- deg_summary$diseases_DEG[idx]  # NA si no es DEG
  
  g
}


#---- Añadir atributo Drug Targets (Chembl) ----

chmbl <- arrow::open_dataset("./Data/Benchmark/chmbl") %>%
  filter( clinicalStage %in% c("PHASE_4", "APPROVAL", "PREAPPROVAL", "PHASE_3")) %>%
  dplyr::select(diseaseId, targetId)%>%
  collect() 

add_layer_drug_target <- function(g, chmbl) {
  
  drug_summary <- chmbl %>%
    dplyr::filter(targetId %in% V(g)$name) %>%
    dplyr::group_by(targetId) %>%
    dplyr::summarise(
      diseases_drug = paste(sort(unique(diseaseId)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(name = targetId)
  
  g <- add_layer(g, drug_summary$name, "is_drug_target")
  
  idx <- match(V(g)$name, drug_summary$name)
  V(g)$diseases_drug <- drug_summary$diseases_drug[idx]
  
  g
}


#---- Añadir atributo DRug Target de un área / enfermedades específicas ----

add_layer_drug_target_area <- function(g, chmbl, ids_enfermedades_interes) {
  
  drug_summary <- chmbl %>%
    dplyr::filter(
      targetId %in% V(g)$name,
      diseaseId %in% ids_enfermedades_interes
    ) %>%
    dplyr::group_by(targetId) %>%
    dplyr::summarise(
      diseases_drug_area = paste(sort(unique(diseaseId)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(name = targetId)
  
  g <- add_layer(g, drug_summary$name, "is_drug_target_area")
  idx <- match(V(g)$name, drug_summary$name)
  V(g)$diseases_drug_area <- drug_summary$diseases_drug_area[idx]
  g
}


#---- Añadir atributo variantes raras ----

add_layer_rare_variants <- function(g, variantes_raras_codificantes,
                                    ids_enfermedades_interes = NULL) {
  
  if (is.matrix(variantes_raras_codificantes)) {
    variantes_raras_codificantes <- as.data.frame(variantes_raras_codificantes,
                                                  stringsAsFactors = FALSE)
  }
  df_filtrado <- variantes_raras_codificantes %>%
    dplyr::filter(targetFromSourceId %in% V(g)$name)
  
  if (!is.null(ids_enfermedades_interes)) {
    df_filtrado <- df_filtrado %>%
      dplyr::filter(diseaseFromSourceMappedId %in% ids_enfermedades_interes)
  }
  
  variant_summary <- df_filtrado %>%
    dplyr::group_by(targetFromSourceId) %>%
    dplyr::summarise(
      diseases_rare_variant = paste(sort(unique(diseaseFromSourceMappedId)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(name = targetFromSourceId)
  
  g <- add_layer(g, variant_summary$name, "has_rare_variant")
  
  idx <- match(V(g)$name, variant_summary$name)
  V(g)$diseases_rare_variant <- variant_summary$diseases_rare_variant[idx]
  
  g
}


#---- Añadir atributo CRISPR_KO ----

add_layer_essential <- function(g, depmap_essential, crisprbrain_essential) {
  
  # --- DepMap: una capa por grupo ---
  grupos <- c("Neural", "Digestive", "Immune", "Other")
  
  for (grupo in grupos) {
    genes_grupo <- depmap_essential %>%
      dplyr::filter(grupo == !!grupo) %>%
      dplyr::pull(ENSG)
    
    attr_name <- paste0("is_essential_", tolower(grupo))
    g <- add_layer(g, genes_grupo, attr_name)
  }
  
  # --- CRISPRbrain ---
  g <- add_layer(g, crisprbrain_essential$ENSG, "is_essential_neuron")
  
  # --- Si es housekeeping
  V(g)$is_housekeeping <- V(g)$is_essential_neural &
                          V(g)$is_essential_digestive &
                          V(g)$is_essential_immune &
                          V(g)$is_essential_other
  g
}

# -----------------------------------------------
# DIBUJAR EL GRAFO
# -----------------------------------------------
library(ggnewscale)
library(ggraph)
library(ggtext)

plot_network <- function(g, titulo = "Gene network", colores_db,
                         seed = 42, exclude_db = NULL,
                         use_qgraph_layout = FALSE,
                         remove_unannotated = FALSE,
                         drug_target_layer = c("any", "area"),
                         string_score_min = NULL,
                         brain_score_min = NULL) {
  
  drug_target_layer <- match.arg(drug_target_layer)
  
  # -- Filtrar por score de STRING --
  if (!is.null(string_score_min)) {
    g <- igraph::delete_edges(g, which(
      E(g)$sourceDatabase == "string" & 
        (is.na(E(g)$scoring) | E(g)$scoring < string_score_min)
    ))
    # Elimina nodos que queden aislados sin capas
    g <- igraph::delete_vertices(g, which(
      igraph::degree(g) == 0 &
        !V(g)$is_DEG &
        !V(g)$is_drug_target &
        !V(g)$has_rare_variant &
        !V(g)$is_seed
    ))
  }
  
  # ─--Filtrar aristas si se excluye DB ---
  if (!is.null(exclude_db)) {
    
    attrs <- data.frame(
      name             = V(g)$name,
      symbol           = V(g)$symbol,
      is_DEG           = V(g)$is_DEG,
      is_drug_target   = V(g)$is_drug_target,
      has_rare_variant = V(g)$has_rare_variant,
      is_seed          = V(g)$is_seed,
      seed_type        = V(g)$seed_type,
      stringsAsFactors = FALSE
    )
    
    g <- igraph::delete_edges(g, which(E(g)$sourceDatabase %in% exclude_db))
    
    # Solo elimina nodos aislados SIN ninguna capa biológica
    g <- igraph::delete_vertices(g, which(
      igraph::degree(g) == 0 &
        !attrs$is_DEG[match(V(g)$name, attrs$name)] &
        !attrs$is_drug_target[match(V(g)$name, attrs$name)] &
        !attrs$has_rare_variant[match(V(g)$name, attrs$name)] &
        !attrs$is_seed[match(V(g)$name, attrs$name)]
    ))
    
    idx <- match(V(g)$name, attrs$name)
    V(g)$symbol           <- attrs$symbol[idx]
    V(g)$is_DEG           <- attrs$is_DEG[idx]
    V(g)$is_drug_target   <- attrs$is_drug_target[idx]
    V(g)$has_rare_variant <- attrs$has_rare_variant[idx]
    V(g)$is_seed          <- attrs$is_seed[idx]
    V(g)$seed_type        <- attrs$seed_type[idx]
  }
  
  
  # Selecciona qué atributo usar para drug target
  is_drug <- if (drug_target_layer == "any") {
    V(g)$is_drug_target
  } else {
    V(g)$is_drug_target_area
  }
  
  # --- Eliminar nodos sin ninguna capa ---
  if (remove_unannotated) {
    sin_capa <- which(
      !V(g)$is_DEG &
        !is_drug &
        !V(g)$has_rare_variant &
        !V(g)$is_seed
    )
    g <- igraph::delete_vertices(g, sin_capa)
    
    # Recalcula is_drug tras eliminar nodos
    is_drug <- if (drug_target_layer == "any") {
      V(g)$is_drug_target
    } else {
      V(g)$is_drug_target_area
    }
  }
  
  # --- Intersecciones específicas de tejido
  
  colores_tejido <- c(
    "Neural"           = "salmon2",  # rojo
    "Gut_microbiome"   = "#96ac60",  # verde
    "Immune_systemic"  = "#E9C46A",  # amarillo
    "Peripheral"       = "pink4"   # morado
  )
  
  # --- Atributos visuales ---
  V(g)$fill_cat <- dplyr::case_when(
    is_drug ~ "Drug target",
    TRUE    ~ "None"
  )
  
  V(g)$border_color <- ifelse(V(g)$has_rare_variant, "slateblue", "black")
  V(g)$shape        <- ifelse(V(g)$is_DEG, "diamond", "circle")
  
  V(g)$label_fill <- dplyr::case_when(
    V(g)$is_seed & V(g)$seed_type == "microbiome" ~ "Seed microbiome",
    V(g)$is_seed & V(g)$seed_type == "neuro"      ~ "Seed neuro",
    V(g)$is_seed & V(g)$seed_type == "both"       ~ "Seed both",
    TRUE                                           ~ "Not a seed"
  )
  
  V(g)$label_fill_color <- dplyr::recode(V(g)$label_fill,
                                         "Seed microbiome" = "#A8D8A8",
                                         "Seed neuro"      = "#A8C8E8",
                                         "Seed both"       = "#D4A8D8",
                                         "Not a seed"      = "white"
  )
  
  # Si es esencial específico
  label_text <- V(g)$symbol
  
  label_text <- ifelse(
    V(g)$is_essential_neural & !V(g)$is_housekeeping,
    paste0(label_text, " *N"),
    label_text
  )
  label_text <- ifelse(
    V(g)$is_essential_digestive & !V(g)$is_housekeeping,
    paste0(label_text, " *D"),
    label_text
  )
  label_text <- ifelse(
    V(g)$is_essential_immune & !V(g)$is_housekeeping,
    paste0(label_text, " *I"),
    label_text
  )
  label_text <- ifelse(
    V(g)$is_essential_neuron & !V(g)$is_housekeeping,
    paste0(label_text, " *Neu"),
    label_text
  )
  
  V(g)$label_text <- label_text
  
  # Paleta unificada para fill de nodos + leyenda semilla
  fill_values <- c(
    "Drug target"     = "pink2",
    "None"            = "white",
    "Seed microbiome" = "#A8D8A8",
    "Seed neuro"      = "#A8C8E8",
    "Seed both"       = "#D4A8D8",
    "Not a seed"      = "white"
  )
  
  # --- Layout ---
  if (use_qgraph_layout) {
    e          <- igraph::as_edgelist(g, names = FALSE)
    layout_mat <- qgraph::qgraph.layout.fruchtermanreingold(
      e,
      vcount      = igraph::vcount(g),
      area        = 4 * (igraph::vcount(g)^2),
      repulse.rad = igraph::vcount(g)^3
    )
  } else {
    set.seed(seed)
    layout_mat <- NULL
  }
  # --- Plot ---
  ggraph(g, layout = if (use_qgraph_layout) layout_mat else "fr") +
    
    # Arista normal : color por DB
    #geom_edge_link(aes(color = sourceDatabase, width = n_dbs)) +
    #scale_edge_color_manual(values = colores_db, name = "Source DB") +
    #scale_edge_width_continuous(range = c(1.2, 1.5), guide = "none") +
    
    # Arista: mismo color para todas
    geom_edge_link(aes(width = n_dbs), color = "grey70", alpha = 0.6) +
    scale_edge_width_continuous(range = c(0.4, 1.5), guide = "none") +
    
    # Arista paralela: solo las que tienen brain expression
    #geom_edge_parallel(
    # aes(filter = !is.na(brain_score)),
    #  color = "salmon2",
    # width = 0.3,
    #  sep   = unit(3, "pt")
    #) +
    
    # Aristas paralelas: por grupo de tejidos
    
    geom_edge_parallel(aes(filter = !is.na(n_tejidos_Neural), edge_colour = "Neural"), 
                       width = 0.8, sep = unit(3, "pt")) +
    geom_edge_parallel(aes(filter = !is.na(n_tejidos_Gut_microbiome), edge_colour = "Gut_microbiome"), 
                       width = 0.8, sep = unit(6, "pt")) +
    geom_edge_parallel(aes(filter = !is.na(n_tejidos_Immune_systemic), edge_colour = "Immune_systemic"), 
                       width = 0.8, sep = unit(9, "pt")) +
    geom_edge_parallel(aes(filter = !is.na(n_tejidos_Peripheral), edge_colour = "Peripheral"), 
                       width = 0.8, sep = unit(12, "pt")) +
    
    scale_edge_colour_manual(   # para que aparezca en la leyenda
      name   = "Tissue expression",
      breaks = c("Neural", "Gut_microbiome", "Immune_systemic", "Peripheral"),
      values = c("salmon2", "#96ac60", "#E9C46A", "pink4"), 
      labels = c("Neural", "Gut Microbiome", "Immune Systemic", "Peripheral")
    ) +
    
    # Nodos
    geom_node_point(aes(fill  = fill_cat,
                        color = border_color,
                        shape = shape,
                        size  = degree),
                    stroke = 1.2) +
    
    # Puntos fantasma para leyenda de semillas
    geom_point(
      data = data.frame(
        x        = c(Inf, Inf, Inf, Inf),
        y        = c(Inf, Inf, Inf, Inf),
        fill_cat = c("Seed microbiome", "Seed neuro", "Seed both", "Not a seed")
      ),
      aes(x = x, y = y, fill = fill_cat),
      shape = 22, size = 5, inherit.aes = FALSE
    ) +
    
    # scale_fill_manual para nodos + leyenda semilla
    scale_fill_manual(
      name   = "Node type / Seed",
      values = fill_values,
      breaks = c("Drug target", "None",
                 "Seed microbiome", "Seed neuro", "Seed both"),
      labels = c("Drug target", "Not a drug target",
                 "Seed (microbiome)", "Seed (neuro)", "Seed (both)"),
      guide  = guide_legend(override.aes = list(shape = 21, size = 5))
    ) +
    
    scale_color_identity(
      name   = "Rare variant",
      labels = c("slateblue" = "Associated", "black" = "None"),
      breaks = c("slateblue", "black"),
      guide  = guide_legend(override.aes = list(
        fill  = "white",
        shape = 21,
        size  = 5,
        color = c("slateblue", "black")
      ))
    ) +
    scale_shape_manual(
      name   = "DEG",
      values = c("diamond" = 23, "circle" = 21),
      labels = c("diamond" = "Yes", "circle" = "No")
    ) +
    scale_size_continuous(range = c(6, 12), guide = "none") +
    
    # Reset fill para geom_node_label
    ggnewscale::new_scale_fill() +
    
    # Etiquetas
    geom_node_label(
      aes(label = ifelse(is_DEG | is_drug_target | has_rare_variant | is_seed, label_text, NA),
          color = border_color,
          fill  = label_fill_color),
      size     = 5,
      repel    = TRUE,
      fontface = "bold",
      max.overlaps = Inf
    ) +
    scale_fill_identity() +  # identity para el fill de las etiquetas
    
    ggtitle(titulo) +
    labs(
      caption = "Essential Genes:\n*N: Neural   *D: Digestive   *I: Immune   *Neu: Neuron"
    ) +
    theme_graph() +
    theme(legend.position = "right",
          legend.text     = element_text(size = 9),
          legend.title    = element_text(size = 10, face = "bold"),
          plot.caption    = element_text(hjust = 1, face = "italic", size = 9, color = "grey30"))
}





