
source("scripts/00_setup.R")

#------------------------------------------------------------------------------------------------------
# TAXONOMIC x PERIPHERAL NEUROPATHY --> Oxidative phosphorilation, protein transmembrane transport, ATP synthesis
#------------------------------------------------------------------------------------------------------

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
lista_intersecciones_im  <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")

# CaRGAR LO COMÚN

load("Data/Redes/config_red.RData")               # colores_db, prioridad_db, add_layer
aristas_color <- readRDS("Data/Redes/aristas_color.rds")
source("Scripts/Scripts_Limpios/Red_Capas_Finales/01_Funciones_Redes.R")

#------------------------------------------------------------------------------------------------------
# ESTUDIAR SOLAPAMIENTO
#------------------------------------------------------------------------------------------------------

genes_tax_pneuropathy <- lista_intersecciones_np[[8]]$genes_ensembl  # 117 genes
genes_tabla <- lista_intersecciones_np[[8]]$genes_tabla
write.csv2(genes_tabla, file = "./Output/Redes_Capas/Tablas/tax_pneuropathy.csv" )
saveRDS(genes_tabla, file = "./Output/Redes_Capas/Tablas/tax_pneuropathy.rds" )

# Construir grafo base

interseccion <- lista_intersecciones_np[[8]]
g <- build_network(interseccion, aristas_color, colores_db)

# Añadimos seed genes
g <- add_layer_seed(g, interseccion$genes_tabla)

#--------------------------------------------------------
# AÑADIMOS CAPAS / ATRIBUTOS
#--------------------------------------------------------

# ---- GENES DIFERENCIALMENTE EXPRESADOS EN ENFERMEDADES - OTAR  ----

lista_genes_DEG_np <- readRDS("./Output/Piloto_Microbiota/Expresion_Diferencial_OTAR/df_genes_NP.rds")
genes_DEG_np <- lista_genes_DEG_np$targetId
intersect(genes_DEG_np, genes_tax_pneuropathy)

g <- add_layer_DEG(g, lista_genes_DEG_np)

# ---- GENES CON VARIANTES RARAS DE INTERES ----

# DF con variantes raras en lugares codificantes
variantes_raras_codificantes <- readRDS("./Data/Variantes/variantes_raras_codificantes.rds")

# FIltrar por variantes de enfermedades que nos interesan (áreas terapéuticas interesantes)
diseases_mapping <- readRDS("./Data/Diccionarios/diseases_mapping.rds")

areas_interes <- c(
  "EFO_0000618",  # nervous system disease
  "MONDO_0002025" # psychiatric disorder
)
# IDs de enfermedad que pertenecen a esas áreas
ids_enfermedades_interes <- diseases_mapping %>%
  dplyr::filter(therapeuticArea %in% areas_interes) %>%
  dplyr::pull(id) %>%
  unique()
#message("Enfermedades en áreas de interés: ", length(ids_enfermedades_interes))

# Cuántas variantes raras de tus genes caen en enfermedades de interés
variantes_raras_codificantes %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  dplyr::filter(
    targetFromSourceId        %in% V(g)$name,
    diseaseFromSourceMappedId %in% ids_enfermedades_interes
  ) %>%
  dplyr::summarise(
    n_variantes = n(),
    n_genes     = n_distinct(targetFromSourceId)
  ) # 1555 variantes en 10 genes de la red

# Añadimos la capa

g <- add_layer_rare_variants(g, variantes_raras_codificantes, ids_enfermedades_interes)
table(V(g)$has_rare_variant)
V(g)$name[V(g)$has_rare_variant == TRUE] # Qué genes son

# ---- GENES QUE SON DRUG TARGETS ----

intersect(chmbl$targetId, genes_tax_pneuropathy)
g <- add_layer_drug_target(g, chmbl)

# --- Genes drug targets de área específica

g <- add_layer_drug_target_area(g, chmbl, ids_enfermedades_interes)
vertex_attr_names(g)

# ---- ESENCIALIDAD GENËTICA (CRISPR KO)
depmap_essential <- readRDS("./Data/CRISPR_KO/depmap_essential_por_grupo.rds")
crisprbrain_essential <- readRDS("./Data/CRISPR_KO/crisprbrain_essential.rds")

g <- add_layer_essential(g, depmap_essential, crisprbrain_essential)

# ---- Resumen de atributos de nodo ----
data.frame(
  name              = V(g)$name,
  is_DEG            = V(g)$is_DEG,
  is_drug_target    = V(g)$is_drug_target,
  has_rare_variant  = V(g)$has_rare_variant
) %>%
  dplyr::summarise(
    n_total          = n(),
    n_DEG            = sum(is_DEG, na.rm = TRUE),
    n_drug_target    = sum(is_drug_target, na.rm = TRUE),
    n_rare_variant   = sum(has_rare_variant, na.rm = TRUE),
    n_todas_capas    = sum(is_DEG & is_drug_target & has_rare_variant, na.rm = TRUE),
    n_DEG_y_drug     = sum(is_DEG & is_drug_target, na.rm = TRUE),
    n_DEG_y_rare     = sum(is_DEG & has_rare_variant, na.rm = TRUE),
    n_drug_y_rare    = sum(is_drug_target & has_rare_variant, na.rm = TRUE),
    n_ninguna        = sum(!is_DEG & !is_drug_target & !has_rare_variant, na.rm = TRUE)
  )
# ver esencialidad genética
data.frame(
  name                   = V(g)$name,
  symbol                 = V(g)$symbol,
  is_essential_neural    = V(g)$is_essential_neural,
  is_essential_digestive = V(g)$is_essential_digestive,
  is_essential_immune    = V(g)$is_essential_immune,
  is_essential_other     = V(g)$is_essential_other,
  is_essential_neuron    = V(g)$is_essential_neuron
) %>%
  dplyr::mutate(
    n_grupos_essential = .data$is_essential_neural + 
      .data$is_essential_digestive +
      .data$is_essential_immune + 
      .data$is_essential_other,
    is_housekeeping    = n_grupos_essential == 4
  ) %>%
  dplyr::filter(is_housekeeping) %>%
  dplyr::select(symbol, is_essential_neuron, n_grupos_essential)

#n_neural n_digestive n_immune n_other n_neuron
#   30          37       43      37       39

#----------------
# PLOT
#----------------

plot <- plot_network(g, titulo = "Taxonomic x Peripheral Neuropathy", colores_db = colores_db)

cairo_pdf("./Output/Redes_Capas/D_Taxonomic_PNeuropathy.pdf", width = 15, height = 15)
print(plot)
dev.off()

# FILTRADO
plot <- plot_network(g, titulo = "Taxonomic x Peripheral Neuropathy", colores_db = colores_db, 
                     use_qgraph_layout = TRUE,
                     exclude_db = "string", 
                     remove_unannotated = TRUE,
                     drug_target_layer = "area")

cairo_pdf("./Output/Redes_Capas/D_Taxonomic_PNeuropathy_filtrado.pdf", width = 25, height = 25)
print(plot)
dev.off()

# Con Tissue Expression Atlas
plot <- plot_network(g, titulo = "Taxonomic x Peripheral Neuropathy", colores_db = colores_db, 
                     use_qgraph_layout = TRUE,
                     exclude_db = "string", 
                     remove_unannotated = TRUE,
                     drug_target_layer = "area")

cairo_pdf("./Output/Redes_Capas/D_Taxonomic_PNeuropathy_Expression_Atlas.pdf", width = 18, height = 18)
print(plot)
dev.off()


# análisis
atributos_nodos <- as.data.frame(vertex_attr(g))
map_variantes <- atributos_nodos %>%
  dplyr::filter(has_rare_variant == TRUE) %>%
  dplyr::select(name, Symbol_Red = symbol) %>% 
  dplyr::inner_join(as.data.frame(variantes_raras_codificantes), by = c("name" = "targetFromSourceId")) %>%
  dplyr::select(Gen = name, Symbol = Symbol_Red, Rasgo_ID = diseaseFromSourceMappedId) %>%
  dplyr::mutate(Capa = "Variante_Rara")

map_deg <- atributos_nodos %>%
  dplyr::filter(is_DEG == TRUE) %>% 
  dplyr::select(name, Symbol_Red = symbol) %>% 
  dplyr::inner_join(lista_genes_DEG_np, by = c("name" = "targetId")) %>%
  dplyr::select(Gen = name, Symbol = Symbol_Red, Rasgo_ID = diseaseId) %>% 
  dplyr::mutate(Capa = "DEG")

map_drugs <- atributos_nodos %>%
  dplyr::filter(is_drug_target == TRUE) %>% 
  dplyr::select(name, Symbol_Red = symbol) %>% 
  dplyr::inner_join(chmbl, by = c("name" = "targetId")) %>%
  dplyr::select(Gen = name, Symbol = Symbol_Red, Rasgo_ID = diseaseId) %>%
  dplyr::mutate(Capa = "Drug_Target")

#----------------------------------------------------------------------
# 3. Combinamos y Resumimos
#----------------------------------------------------------------------
mapeo_combinado <- dplyr::bind_rows(map_variantes, map_deg, map_drugs)

matriz_enfermedades_completa <- mapeo_combinado %>%
  dplyr::inner_join(traits_MicroGWAS_areas, by = c("Rasgo_ID" = "Rasgo")) %>%
  dplyr::group_by(Rasgo_ID, name, Nombre_area, therapeuticAreas) %>%
  dplyr::summarise(
    Genes_DEG          = paste(unique(Symbol[Capa == "DEG"]), collapse = ", "),
    Genes_DrugTarget   = paste(unique(Symbol[Capa == "Drug_Target"]), collapse = ", "),
    Genes_VarRaras     = paste(unique(Symbol[Capa == "Variante_Rara"]), collapse = ", "),
    Num_DEG            = n_distinct(Gen[Capa == "DEG"]),
    Num_DrugTarget     = n_distinct(Gen[Capa == "Drug_Target"]),
    Num_VarRaras       = n_distinct(Gen[Capa == "Variante_Rara"]),
    Total_Genes_Red    = n_distinct(Gen),
    .groups = "drop"
  ) %>%
  dplyr::rename(Enfermedad_Nombre = name) %>%
  dplyr::arrange(desc(Total_Genes_Red))

# Ver el resultado final
View(matriz_enfermedades_completa) 











