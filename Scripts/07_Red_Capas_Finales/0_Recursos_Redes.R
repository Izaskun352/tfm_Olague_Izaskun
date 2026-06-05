

# RECURSOS CONSTRUIR REDES FINALES CON CAPAS

# -----------------------------------
# INTERACTOMA + ARISTAS
# -----------------------------------

# CADA DATABASE UN COLOR DIFERENte
colores_db <- c(
  "intact"   = "#84b6f4",  
  "signor"   = "pink4",  
  "reactome" = "#96ac60",  
  "string"   = "seashell3"   
)
interactoma_edges_withSource <- readRDS ("./Data/nasertic/input/interactome_withSource_nodeFilter.rds")

interactoma_multidb <- as.data.frame(interactoma_edges_withSource) %>%
  dplyr::rename(ENSG_A = targetA, ENSG_B = targetB) %>%
  dplyr::filter(ENSG_A != ENSG_B) %>%
  dplyr::mutate(
    par_id = paste(pmin(ENSG_A, ENSG_B), pmax(ENSG_A, ENSG_B), sep = "_")
  ) %>%
  dplyr::distinct(par_id, sourceDatabase, .keep_all = TRUE)  

prioridad_db <- c("intact" = 1, "reactome" = 2, "signor" = 3, "string" = 4)

interactoma_multidb <- interactoma_multidb %>%
  dplyr::mutate(prioridad = prioridad_db[sourceDatabase])

# GENERAR INTERACTOMA CON SCORE STRING

#url <- "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/interaction/"
#carpeta_entrada<-"./Data/Interactoma"

#dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
#archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
#lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))

#interactoma_string<-arrow::open_dataset("./Data/Interactoma") %>%
# select(targetA, targetB, sourceDatabase, scoring) %>%
#filter(!is.na(targetA), !is.na(targetB)) %>% #elimina las filas con vacios
#  filter(sourceDatabase == "string", scoring >= 0.4) %>% #Filtramos SOLO las de String con score>=0.4
# filter(targetA != targetB) %>% #Elimina iteraciones de un gen consigomismo
# distinct() %>% #Elimina filas duplicadas
# collect() %>%
# as.data.frame()
#saveRDS(interactoma_string, file = "./Data/Diccionarios/Interactoma_String.rds")

interactoma_string <- readRDS("./Data/Diccionarios/Interactoma_String.rds")

# --- ARISTAS GRAFO ----
aristas_color <- interactoma_multidb %>%
  dplyr::group_by(par_id) %>%
  dplyr::mutate(
    dbs_presentes  = paste(sort(unique(sourceDatabase)), collapse = ";"),
    n_dbs          = n_distinct(sourceDatabase)
  ) %>%
  dplyr::slice_min(prioridad, n = 1, with_ties = FALSE) %>%  # quédate con la DB prioritaria
  dplyr::ungroup()

aristas_color <- aristas_color %>%
  dplyr::left_join(
    interactoma_string %>% 
      dplyr::mutate(par_id = paste(pmin(targetA, targetB), 
                                   pmax(targetA, targetB), sep = "_")) %>%
      dplyr::select(par_id, scoring),
    by = "par_id"
  )


# Añadir información del atlas: intersecciones en cerebro (sin CV)
brain_expression_atlas <- readRDS("./Data/Redes/brain_expression_atlas.rds")
brain_pairs <- brain_expression_atlas %>%
  dplyr::mutate(par_id = paste(pmin(ENSG_A, ENSG_B),
                               pmax(ENSG_A, ENSG_B), sep = "_")) %>%
  dplyr::select(par_id, brain_score = brain)

aristas_color <- aristas_color %>%
  dplyr::left_join(brain_pairs, by = "par_id")

saveRDS(interactoma_multidb, "Data/Redes/interactoma_multidb.rds")      
save(colores_db, prioridad_db, 
     file = "Data/Redes/config_red.RData")

# Añadir informacion de Evidence Expression Atlas -- todos los tejidos

# intersecciones de todos los tejidos
# filtros: score > 0.8 y CV > 0.4
# tejidos agrupados en 4 grupos grandes: neural, gut-microbiome, immune_systemic, peropheral

# pivotamos atlas para que haya una fila por pareja 
atlas_wide <- atlas_merged %>%
  tidyr::pivot_wider(
    id_cols     = par_id,
    names_from  = grupo,
    values_from = n_tejidos,
    names_prefix = "n_tejidos_"
  )

# unimos a aristas_color

aristas_color <- aristas_color %>%
  dplyr::left_join(atlas_wide, by = "par_id")
saveRDS(aristas_color,       "Data/Redes/aristas_color.rds")
