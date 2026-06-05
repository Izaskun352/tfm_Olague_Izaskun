
# LIBRERIAS

#source("scripts/00_setup.R")
#library(clusterProfiler)
#library(biomaRt)

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# OUTPUT

TERM2GENE_interpro <- readRDS("./Data/Diccionarios/Anotaciones_InterPro/TERM2GENE_InterPro.rds")
TERM2NAME_interpro <- readRDS("./Data/Diccionarios/Anotaciones_InterPro/TERM2NAME_InterPro.rds")

# -------------------------------------
# DESCARGAR INTERPRO Y CREAR TERM2GENE (BioMart)
# -------------------------------------
  
  # ---- 1. Conexión y descarga en biomaRT
  
  mart <- useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl",
    mirror  = "useast"    # prueba: "useast", "asia", "uswest"
  )
  cromosomas <- c(1:22, "X", "Y", "MT")
  lista_interpro <- list()
  
  for (chr in cromosomas) {
    
    message("Descargando cromosoma ", chr, "...")
    
    resultado <- tryCatch(
      getBM(
        attributes = c("ensembl_gene_id",
                       "interpro",
                       "interpro_description"),
        filters    = "chromosome_name",
        values     = chr,
        mart       = mart
      ),
      error = function(e) {
        message("  -> ⚠️ Timeout en chr", chr, ", reintentando en 15s...")
        Sys.sleep(15)
        tryCatch(
          getBM(
            attributes = c("ensembl_gene_id",
                           "interpro",
                           "interpro_description"),
            filters    = "chromosome_name",
            values     = chr,
            mart       = mart
          ),
          error = function(e2) {
            message("  -> ❌ Fallo definitivo en chr", chr)
            return(NULL)
          }
        )
      }
    )
    
    if (!is.null(resultado)) {
      lista_interpro[[as.character(chr)]] <- resultado
      message("  -> ✔️ Chr", chr, ": ", nrow(resultado), " anotaciones")
    }
    
    Sys.sleep(3)  # pausa entre peticiones
  }
  
  # --- 2. Unir y limpiar
  
  interpro_anotaciones <- bind_rows(lista_interpro) %>%
    filter(interpro != "" & !is.na(interpro)) %>%
    distinct()
  
  message("Total pares gen-dominio: ", nrow(interpro_anotaciones))
  
  
  # --- 3. CONSTRUIR TERM2GENE Y TERM2NAME 
  # Estos son los dos objetos que necesita enricher()
  
  # TERM2GENE: col1 = ID del término, col2 = gen
  # ⚠ El orden de las columnas importa: primero término, luego gen
  TERM2GENE_interpro <- interpro_anotaciones %>%
    dplyr::select(interpro, ensembl_gene_id) %>%
    distinct()
  
  # TERM2NAME: col1 = ID del término, col2 = nombre legible
  TERM2NAME_interpro <- interpro_anotaciones %>%
    dplyr::select(interpro, interpro_description) %>%
    distinct()
  
  message("Términos InterPro únicos: ", n_distinct(TERM2GENE_interpro$interpro))
  message("Genes únicos anotados:    ", n_distinct(TERM2GENE_interpro$ensembl_gene_id))
  
  # --- 4. GUARDAR PARA NO REPETIR LA DESCARGA 
  
  carpeta_interpro_anotaciones <- "./Data/Diccionarios/Anotaciones_InterPro"
  dir.create(carpeta_interpro_anotaciones, showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(interpro_anotaciones, "./Data/Diccionarios/Anotaciones_InterPro/interpro_anotaciones_raw.rds")
  saveRDS(TERM2GENE_interpro,   "./Data/Diccionarios/Anotaciones_InterPro/TERM2GENE_InterPro.rds")
  saveRDS(TERM2NAME_interpro,   "./Data/Diccionarios/Anotaciones_InterPro/TERM2NAME_InterPro.rds")

