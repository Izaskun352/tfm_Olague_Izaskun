
source("scripts/00_setup.R")

# EXTRAER VARIANTES RARAS CODIFICANTES
# =====================================================================================

# =====================================================================================
# INPUTS
# =====================================================================================
Variantes_comunes <- readRDS("./Data/Variantes/common_filter_allCol.rds")

Variantes_raras <- readRDS("./Data/Variantes/rare_filter_Select_Col.rds")

SO_vector <- readRDS("./Data/Variantes/SO_vector.rds")

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
# =====================================================================================
# OUTPUTS
# =====================================================================================

variantes_raras_codificantes <- readRDS("./Data/Variantes/variantes_raras_codificantes.rds")

# =====================================================================================
# EXTRAER VARIANTES RARAS CODIFICANTES  ====

# ponemos ID en el mismo formato

names(SO_vector) <- gsub(":", "_", names(SO_vector))

# Definimos patrones biológicos que indican que se encuentra en posición codificante

coding_patterns <- c(
  "missense",
  "frameshift",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "synonymous",
  "protein_altering",
  "coding_sequence",
  "inframe_insertion",
  "inframe_deletion",
  "splice_donor",
  "splice_acceptor"
)

# Extraemos los ID
coding_SO <- names(SO_vector)[
  grepl(paste(coding_patterns, collapse="|"), SO_vector)
]

# Limpiar falsos positivos

exclude_patterns <- c(
  "intron",
  "intergenic",
  "upstream",
  "downstream",
  "UTR",
  "regulatory"
)

coding_SO <- coding_SO[
  !grepl(paste(exclude_patterns, collapse="|"), SO_vector[coding_SO])
]

unique(SO_vector[coding_SO])

# Extraemos variantes raras codificantes
Variantes_raras_codificantes <- Variantes_raras[
  Variantes_raras[,4] %in% coding_SO,
]
unique(SO_vector[unique(Variantes_raras_codificantes[,4])]) # Ver cuales hemos seleccionado
head(Variantes_codificantes)

saveRDS(Variantes_raras_codificantes, file = "./Data/Variantes/variantes_raras_codificantes.rds")

# =====================================================================================

# =====================================================================================
# EXTRAER GENES INTERSECCIN MICROGWAS + NP

genes_np_micro_todos <- unlist(
  lapply(lista_intersecciones_np, function(x) x$genes_ensembl)
)
genes_np_micro_todos <- unique(genes_np_micro_todos) # 724 genes en intersección micro - NP

# =====================================================================================


# =====================================================================================
# EXTRAEMOS GENES EN INTERSECCIÓN CON VARIANTE RARA CODIFICANTE
#=====================================================================================
# Hacemos diccionario ID-nombre_gen

mapping_list <- lapply(lista_intersecciones_np, function(x) {
  data.frame(
    ensembl = x$genes_ensembl,
    simbolo = x$genes_simbolo,
    stringsAsFactors = FALSE
  )
})
mapping_df <- unique(do.call(rbind, mapping_list))
mapping_df <- mapping_df[mapping_df$simbolo != "", ]

# Sacamos genes en inersección con variantes raras codificantes

variantes_en_genes <- as.data.frame(variantes_raras_codificantes[
  variantes_raras_codificantes[,3] %in% genes_np_micro_todos,
])
head(variantes_en_genes)

variantes_en_genes <- variantes_en_genes %>%   # 338 genes con variante rara codificante en la intersección NP-micro
  dplyr::left_join(mapping_df, 
                   by = c("targetFromSourceId" = "ensembl")) %>%
  dplyr::rename(gene_name = simbolo) %>%   
  dplyr::left_join(traits_MicroGWAS_areas, 
                   by = c("diseaseFromSourceMappedId" = "Rasgo")) %>%
  dplyr::select(-Nombre_area, -therapeuticAreas)

SO_df <- data.frame(
  SO_id = names(SO_vector),
  SO_nombre = unname(SO_vector)
)

variantes_en_genes <- variantes_en_genes %>%
  dplyr::mutate(SO_id = gsub("_", ":", variantFunctionalConsequenceId)) %>%
  dplyr::left_join(SO_df, by = "SO_id")

genes_unicos_con_variantes <- unique(variantes_en_genes[,3])

# Sacamos nombre de los genes
genes_unicos_con_variantes <- mapping_df[mapping_df$ensembl %in% genes_unicos_con_variantes, ]

# Vemos si alguno es gen semilla

genes_semilla <- unlist(
  lapply(lista_intersecciones_np, function(x) {
    tabla <- x$genes_tabla
    tabla$gene.ENSG[tabla$es_semilla_en_alguno]
  })
)
genes_semilla <- unique(genes_semilla)
genes_semilla_con_variantes <- intersect(
  genes_semilla,
  genes_unicos_con_variantes$ensembl
)
genes_semilla_variante <- genes_unicos_con_variantes[
  genes_unicos_con_variantes$ensembl %in% genes_semilla_con_variantes,
]


kmt2c <- variantes_en_genes[variantes_en_genes$gene_name == "HDAC4",]
unique(kmt2c$diseaseFromSourceMappedId)
unique(kmt2c$name)
table(kmt2c$diseaseFromSourceMappedId)



#======================================================================================
# VARIANTES EN DINEINA

library(biomaRt)
library(dplyr)

# 1. Conectar a la base de datos de Ensembl humana
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 2. Definir el término GO del complejo de dineína
termino_go_dineina <- "GO:0030286"

# 3. Extraer los IDs de Ensembl asociados a ese término GO
genes_dineina_df <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
  filters = 'go',
  values = termino_go_dineina,
  mart = ensembl
)

# 4. Aislar solo los IDs de Ensembl en un vector limpio, quitando duplicados
vector_genes_dineina <- unique(genes_dineina_df$ensembl_gene_id)

# 5. Filtrar tu dataframe original (asumiendo que se llama 'df_variantes')
# Te quedas solo con las filas donde el gen (targetFromSourceId) esté en tu vector
df_dineinas_filtrado <- df_variantes %>%
  filter(targetFromSourceId %in% vector_genes_dineina)

# Ver cuántas variantes te han quedado
nrow(df_dineinas_filtrado) 

