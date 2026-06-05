
library(tidyverse)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)

# INPUTS
interactoma <- readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds")
universo_genes <- unique(na.omit(interactoma[,1]))

# ESTUDIAR GOBP DE LOS MODULOS PLEIOTROPICOS DE MICROGWAS - VARIACION COMUN 

red_pleiotropia_varComun_MicroGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/Red_Pleiotropia.rds")

# --- Aislamos los componentes que queremnos

comp <- igraph::components(red_pleiotropia_varComun_MicroGWAS)
nodos <- names(comp$membership) 

# --- Sacamos genes solapados de cada megacluster

  # Carpetas base
ruta_base <- "./Output/Piloto_Microbiota/" 
carpetas <- c(
  "MicroGWAS"          = file.path(ruta_base, "Clusters_MicroGWAS"),
  "Immune"             = file.path(ruta_base, "Clusters_Immune"),
  "NervousPsychiatric" = file.path(ruta_base, "Clusters_NervousPsychiatric"),
  "VarComun"           = file.path(ruta_base, "Clusters_Traits_VarComun")
)
map_lgl(carpetas, dir.exists)

  # Funcion para leer los genes de cada cluster

leer_genes_cluster <- function(nombre_nodo) {
  
  # Intentar primero con el nombre tal cual (con ZSCO. si lo tiene)
  # y luego sin el prefijo
  nombres_a_probar <- c(
    nombre_nodo,                                    # con ZSCO. si lo tiene
    str_remove(nombre_nodo, "^ZSCO\\.")             # sin ZSCO.
  )
  
  archivo <- NULL
  for (nombre in nombres_a_probar) {
    for (carpeta in carpetas) {
      ruta <- file.path(carpeta, paste0(nombre, ".csv"))
      if (file.exists(ruta)) {
        archivo <- ruta
        break
      }
    }
    if (!is.null(archivo)) break
  }
  
  if (is.null(archivo)) {
    warning(paste("No encontrado:", nombre_nodo))
    return(NULL)
  }
  
  read.csv2(archivo) %>%
    pull(gene.ENSG) %>%
    na.omit() %>%
    unique()
}
  
  # Leer los genes de cada nodo
genes_por_nodo <- map(
  set_names(nodos),
  leer_genes_cluster
) %>% compact()
length(genes_por_nodo)

cat("Clusters leídos con éxito:", length(genes_por_nodo), "\n")

  # Agrupar los nodos por componente y solapar

componente_de   <- comp$membership[names(genes_por_nodo)]
print(componente_de)
nodos_por_comp  <- split(names(genes_por_nodo), componente_de)

interseccion_por_comp <- map(nodos_por_comp, function(nodos) {
  genes <- map(nodos, ~ genes_por_nodo[[.x]]) %>% compact()
  
  if (length(genes) == 0) return(character(0))
  if (length(genes) == 1) return(genes[[1]])
  
  Reduce(intersect, genes)
})

  # Ver cuántos genes hay en cada intersección
map_int(interseccion_por_comp, length)


# ---- ENRIQUECIMIENTO GOBP POR MEGACLUSTER / COMPONENTE ----

enrich_por_comp <- imap(interseccion_por_comp, function(genes, comp_id) {
  cat("\n── Componente", comp_id, "| genes en intersección:", length(genes), "\n")
  
  if (length(genes) < 5) {
    cat("   Pocos genes, omitiendo enriquecimiento\n")
    return(NULL)
  }
  resultado_go <- enrichGO(
    gene          = genes,
    universe      = universo_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  if (is.null(resultado_go) ||
      nrow(dplyr::filter(resultado_go@result, p.adjust < 0.05)) == 0) {
    message("  ⚠️ Sin términos GO significativos"); next
  }
  
  res <- tryCatch(
    simplify(resultado_go, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) resultado_go
  )
  
  if (!is.null(res)) cat("   Términos significativos:", nrow(res), "\n")
  res
})
  
  
# ---- RESULTADOS ----
  
  # Tabla resumen de todos los componentes juntos
resumen_enrich <- imap_dfr(enrich_por_comp, function(res, comp_id) {
    as.data.frame(res) %>%
      mutate(componente = comp_id) %>%
      dplyr::select(componente, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
  })
saveRDS(resumen_enrich, file = "./Output/resumen_enrich.rds")
resumen_enrich <- readRDS("./Output/resumen_enrich.rds")
  # Ver top 5 por componente
resumen_enrich %>%
    group_by(componente) %>%
    slice_min(p.adjust, n = 5) %>%
    print(n = Inf)

  # Hacemos tabla componente - enfermedades - top10 términos GOBP

enfermedades_por_comp <- map(nodos_por_comp, function(nodos) {
  vertex_attr(red_pleiotropia_varComun_MicroGWAS, "label", index = nodos) %>%
    unique() %>%
    paste(collapse = " | ")
})

top10_por_comp <- resumen_enrich %>%
  group_by(componente) %>%
  slice_min(p.adjust, n = 10) %>%
  summarise(top10_GOBP = paste(Description, collapse = " | "))

tabla_final <- top10_por_comp %>%
  mutate(
    enfermedades = map_chr(componente, ~ enfermedades_por_comp[[.x]] %||% "Sin datos")
  ) %>%
  dplyr::select(componente, enfermedades, top10_GOBP)


# ---- HEATMAP ----

# --- 1. Preparar datos: top 3 por componente ---
tamanio_comp <- map_int(nodos_por_comp, length)
top3_datos <- resumen_enrich %>%
  group_by(componente) %>%
  slice_min(p.adjust, n = 3) %>%
  ungroup() %>%
  mutate(
    valor      = 1,   # si o no
    componente = paste0("Comp_", componente, 
                        " (n=", tamanio_comp[componente], ")")
  ) %>%
  dplyr::select(componente, Description, valor)

# --- 2. Pasar a matriz (filas = términos GOBP, columnas = componentes) ---
mat <- top3_datos %>%
  pivot_wider(
    names_from  = componente,
    values_from = valor,
    values_fill = 0
  ) %>%
  column_to_rownames("Description") %>%
  as.matrix()

# --- 3. Heatmap ---

#pdf("./Output/Gráficos/Gráficos_Resultados_4/Heatmap_GOBP_ModulosPleiotropicos_VarComun.pdf", width = 15, height = 20)
pdf("./Output/Gráficos/Heatmap_GOBP_ModulosPleiotropicos_VarComun.pdf", width = 15, height = 20)
pheatmap(
  mat,
  color = c("cornsilk2", "skyblue4"),
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  display_numbers  = FALSE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  border_color     = "white",
  angle_col        = 45,
  cellwidth  = 25,  # anchura de cada celda en puntos
  cellheight = 20,  # altura de cada celda en puntos
  fontsize_row      = 14,
  fontsize_col      = 14,
  legend_breaks    = c(0, 10, 20, 30, 40,50,60,70,80,90,100),
  legend_labels    = c("0", "10", "20", "30", "40","50", "60", "70", "80", "90", "100"),
  name             = "-log10(p.adjust)",
  main             = "Top 3 términos GOBP por megacluster"
)
dev.off()


