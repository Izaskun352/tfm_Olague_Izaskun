
# LIBRERIAS

source("scripts/00_setup.R") 
library("depmap")
library("ExperimentHub")
library("BioMart")

# ===================================================
# ESENCIALIDAD GENËTICA -- CRISPRS KNOCKOUT
# ===================================================

# INPUTS

# OUTPUTS
depmap_essential <- readRDS("./Data/CRISPR_KO/depmap_essential_por_grupo.rds")
crisprbrain_essential <- readRDS("./Data/CRISPR_KO/crisprbrain_essential.rds")

#---------------------------
# DEPMAP ----
#---------------------------

# Data
crispr_data <- depmap::depmap_crispr()
metadata <- depmap_metadata()

linajes_unicos <- unique(metadata$lineage)
conteo_linajes <- table(metadata$lineage)
conteo_linajes <- sort(conteo_linajes, decreasing = TRUE)

# Agrupamos tejidos

grupo_linaje_depmap <- c(
  # Neural
  "central_nervous_system"    = "Neural",
  "peripheral_nervous_system" = "Neural",
  
  # Digestive
  "colorectal"                = "Digestive",
  "gastric"                   = "Digestive",
  "esophagus"                 = "Digestive",
  "upper_aerodigestive"       = "Digestive",
  "liver"                     = "Digestive",
  "pancreas"                  = "Digestive",
  "bile_duct"                 = "Digestive",
  
  # Immune
  "blood"                     = "Immune",
  "lymphocyte"                = "Immune",
  "plasma_cell"               = "Immune",
  
  # Other
  "lung"                      = "Other",
  "skin"                      = "Other",
  "breast"                    = "Other",
  "bone"                      = "Other",
  "ovary"                     = "Other",
  "soft_tissue"               = "Other",
  "kidney"                    = "Other",
  "uterus"                    = "Other",
  "urinary_tract"             = "Other",
  "thyroid"                   = "Other",
  "eye"                       = "Other",
  "prostate"                  = "Other",
  "cervix"                    = "Other",
  "adrenal_cortex"            = "Other",
  "fibroblast"                = "Other",
  "embryo"                    = "Other",
  "epidermoid_carcinoma"      = "Other",
  "unknown"                   = NA
)

# añadir lineaje a cripr_data
crispr_data <- crispr_data %>%
  dplyr::left_join(
    metadata %>% dplyr::select(depmap_id, lineage),
    by = "depmap_id"
  ) %>%
  dplyr::mutate(grupo = grupo_linaje_depmap[lineage]) %>%
  dplyr::filter(!is.na(grupo))
# Añadimos el grupo a cada línea celular
crispr_data <- crispr_data %>%
  dplyr::mutate(grupo = grupo_linaje_depmap[lineage]) %>%
  dplyr::filter(!is.na(grupo))

# filtro: el gen es esencial para al menos el 25% de líneas (dependency < -0.5) del grupo
depmap_essential <- crispr_data %>%
  dplyr::group_by(gene_name, grupo) %>%
  dplyr::summarise(
    n_lineas        = n(),
    n_essential     = sum(dependency < -0.5, na.rm = TRUE),
    pct_essential   = n_essential / n_lineas,
    .groups = "drop"
  ) %>%
  dplyr::filter(pct_essential >= 0.25) %>%  # esencial en al menos 25% de líneas
  dplyr::select(gene_name, grupo)
  
# mapear gene.symbol --> ENSG
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

simbolos_depmap <- unique(depmap_essential$gene_name)
mapping_depmap <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = simbolos_depmap,
  mart       = mart
) %>%
  dplyr::rename(gene_name = hgnc_symbol, ENSG = ensembl_gene_id) %>%
  dplyr::filter(ENSG != "") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

depmap_essential <- depmap_essential %>%
  dplyr::left_join(mapping_depmap, by = "gene_name") %>%
  dplyr::filter(!is.na(ENSG))
# tenemos archivo con los genes esenciales para algun grupo en al menos el 25% de las líneas
saveRDS(depmap_essential, "./Data/CRISPR_KO/depmap_essential_por_grupo.rds")

#------------------------------------
# CRIPRbrain
#------------------------------------

# Data
crisprbrain_data <- read.csv("./Data/CRIPRbrain_iCRISPR_Survival_GlutamericNeurons.csv", stringsAsFactors = FALSE)
simbolos_crisprbrain <- unique(crisprbrain_data$TSS[crisprbrain_data$Hit.Class == "Negative Hit"])

mapping_crisprbrain <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = simbolos_crisprbrain,
  mart       = mart
) %>%
  dplyr::rename(gene_name = hgnc_symbol, ENSG = ensembl_gene_id) %>%
  dplyr::filter(ENSG != "") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

crisprbrain_essential <- crisprbrain_data %>%
  dplyr::filter(Hit.Class == "Negative Hit") %>%
  dplyr::rename(gene_name = TSS) %>%
  dplyr::left_join(mapping_crisprbrain, by = "gene_name") %>%
  dplyr::filter(!is.na(ENSG)) %>%
  dplyr::select(gene_name, ENSG, Gene.Score, P.Value)

saveRDS(crisprbrain_essential, "./Data/CRISPR_KO/crisprbrain_essential.rds")

message("Genes esenciales en neuronas (CRISPRbrain): ", n_distinct(crisprbrain_essential$ENSG))




