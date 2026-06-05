

# INTERSECCIONES ENTRE PROTEINAS DEPENDIENTES DE TEJIDO

# INPUTS

matriz_prob <- read_csv("./Data/Tissue_Expression_Atlas/association_scores/cohorts_combined_brain_tumor_avg_outer_prob.csv")
df_conversion <- readRDS("./Data/Tissue_Expression_Atlas/conversion_df_ENSG.rds")

carpeta_data_tejidos_cv_filtrado <- "./Data/Tissue_Expression_Atlas/Data_CV_Tejidos_filtrado/"

# OUTPUTS
atlas_merged <-  readRDS("./Data/Tissue_Expression_Atlas/atlas_merged.rds")

# ---- SIN TENER EN CUENTA CV ----
matriz_prob <- read_csv("./Data/association_scores/cohorts_combined_brain_tumor_avg_outer_prob.csv")   # 3 columnas y 67847593 filas

brain_expression_atlas<- matriz_prob %>%
  # Eliminar interacciones consigo mismo y aplicar el umbral de confianza estricto
  filter(prot1 != prot2, brain > 0.8) %>%
  
  # Reordenar alfabéticamente por fila para unificar las direcciones (A-B y B-A)
  mutate(
    p1 = pmin(prot1, prot2),
    p2 = pmax(prot1, prot2)
  ) %>%
  
  # Eliminar interacciones repetidas basándonos en el nuevo orden homogéneo
  distinct(p1, p2, .keep_all = TRUE) %>%
  # seleccionamos otra vez las columnas
  dplyr::select(prot1 = p1, prot2 = p2, brain)


# ---- TENIENDO EN CUENTA CV ----
# CARGAMOS LOS ARCHIVOS FILTRADOS (score > 0.8 y CV> 0.4 para cada tejido)

archivos <- list.files(carpeta_data_tejidos_cv_filtrado, 
                       pattern = "\\.rds$",  
                       full.names = TRUE)

nombres_tejidos <- gsub("_tumor_CV_filtrado\\.rds$", "", basename(archivos))  # extraer nombre tejido

atlas_list <- lapply(archivos, readRDS) 
names(atlas_list) <- nombres_tejidos

# Tejido - grupo 
grupo_tejido <- c(
  "brain"    = "Neural",
  "colon"    = "Gut_microbiome",
  "stomach"  = "Gut_microbiome",
  "liver"    = "Gut_microbiome",
  "pancreas" = "Gut_microbiome",
  "throat"   = "Gut_microbiome",
  "blood"    = "Immune_systemic",
  "lung"     = "Immune_systemic",
  "breast"   = "Peripheral",
  "kidney"   = "Peripheral",
  "ovary"    = "Peripheral"
)

# Combinar atlas
atlas_merged <- bind_rows(atlas_list, .id = "tejido") %>%
  dplyr::mutate(
    grupo  = grupo_tejido[tejido],
    par_id = paste(pmin(ENSG1, ENSG2), pmax(ENSG1, ENSG2), sep = "_")
  ) %>%
  dplyr::group_by(par_id, grupo) %>%
  dplyr::summarise(
    tejidos_presentes = paste(sort(unique(tejido)), collapse = ";"),
    n_tejidos         = n_distinct(tejido),
    .groups           = "drop"
  )
saveRDS(atlas_merged, "./Data/Tissue_Expression_Atlas/atlas_merged.rds")

