
source("scripts/00_setup.R")

# GENES DIFERENCIALMENTE EXPRESADOS EN ENFERMEDADES

# =====================================================================================
# INPUTS
# =====================================================================================

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

traits_NP <- traits_MicroGWAS_areas %>%
  filter(
    grepl("EFO_0000618",  therapeuticAreas, ignore.case = TRUE) |
      grepl("MONDO_0002025", therapeuticAreas, ignore.case = TRUE)
  ) %>%
  pull(Rasgo) %>%
  trimws() %>%
  unique()

traits_IM <- traits_MicroGWAS_areas %>%
  filter(
    grepl("EFO_0000540",  therapeuticAreas, ignore.case = TRUE)
  ) %>%
  pull(Rasgo) %>%
  trimws() %>%
  unique()

traits_PSY <- traits_MicroGWAS_areas %>%
  filter(
    grepl("MONDO_0002025",  therapeuticAreas, ignore.case = TRUE)
  ) %>%
  pull(Rasgo) %>%
  trimws() %>%
  unique()

# CREAMOS CARPETA

carpeta_diff_expresion <- "./Output/Piloto_Microbiota/Expresion_Diferencial_OTAR"
dir.create(carpeta_diff_expresion, showWarnings = FALSE)

# 1--- Descargamos dataset --> Differential expression evidence

url <- "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_expression_atlas/"
carpeta_entrada<-"./Data/Diccionarios/Diff_expr_OTAR"
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


# 2 --- Filtramos dataset

ds <- open_dataset(carpeta_entrada)  # arrow 
ds %>% 
  dplyr::select(diseaseId, targetId, confidence, resourceScore, log2FoldChangeValue, studyId) %>%
  head(10) %>% 
  collect()


# 3 --- Filtro: Traits Nervous System / Psychiatric Disorder

genes_NP <- ds %>%
  filter(diseaseId %in% traits_NP) %>%
  collect()
genes_NP %>% head(20) %>% collect()
genes_NP_unicos <- unique(genes_NP$targetFromSourceId)  # 7109 genes - confidence. medium / high

saveRDS(genes_NP, file = file.path(carpeta_diff_expresion, "df_genes_NP.rds"))
readRDS("./Output/Piloto_Microbiota/Expresion_Diferencial_OTAR/df_genes_NP.rds")


# 4 --- Filtro: Traits Immune System

genes_IM <- ds %>%
  filter(diseaseId %in% traits_IM) %>%
  collect()

genes_IM_unicos <- unique(genes_IM$targetFromSourceId)

genes_IM_medium_high <- genes_IM %>%
  filter(confidence %in% c("medium", "high")) %>%
  collect()
genes_IM_medium_high_unicos <- unique(genes_IM_medium_high$targetFromSourceId)   # 11613 genes - confidence medium / high

saveRDS(genes_IM_medium_high, file = file.path(carpeta_diff_expresion, "df_genes_IM.rds"))


# 5 --- Filtro: Traits Psychiatric Disorder

genes_PSY <- ds %>%
  filter(diseaseId %in% traits_PSY) %>%
  collect()

genes_PSY_unicos <- unique(genes_PSY$targetFromSourceId)  # 4059 genes

genes_PSY_medium_high <- genes_PSY %>%
  filter(confidence %in% c("medium", "high")) %>%
  collect()
genes_PSY_medium_high_unicos <- unique(genes_PSY_medium_high$targetFromSourceId)   # 4059 genes - confidence medium / high




