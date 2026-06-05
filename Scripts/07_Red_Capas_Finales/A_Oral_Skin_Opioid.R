
source("scripts/00_setup.R")

#------------------------------------------------------------------------------------------------------
# ORAL / SKIN x OPIOID DEPENDENCE --> SYNAPTIC TRANSMISSION
#------------------------------------------------------------------------------------------------------

# INPUTS

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
lista_intersecciones_im  <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_IM_MicroGWAS.rds")

# CaRGAR LO COMÚN

load("Data/Redes/config_red.RData")               # colores_db, prioridad_db, add_layer
aristas_color <- readRDS("Data/Redes/aristas_color.rds")
source("Scripts/Scripts_Limpios/Red_Capas_Finales/01_Funciones_Redes.R")

# OUTPUTS



#------------------------------------------------------------------------------------------------------
# ESTUDIAR SOLAPAMIENTO
#------------------------------------------------------------------------------------------------------
genes_oral_op <- lista_intersecciones_np[[25]]$genes_ensembl   # 36
write.csv2(lista_intersecciones_np[[25]]$genes_tabla, file = "./Output/Redes_Capas/Tablas/oral_opioid.csv")
saveRDS(lista_intersecciones_np[[25]]$genes_tabla, file = "./Output/Redes_Capas/Tablas/oral_opioid.rds")
genes_skin_op <- lista_intersecciones_np[[15]]$genes_ensembl   # 43 --> tiene un gen esencial que en los otros no
write.csv2(lista_intersecciones_np[[15]]$genes_tabla, file = "./Output/Redes_Capas/Tablas/skin_opioid.csv")
saveRDS(lista_intersecciones_np[[15]]$genes_tabla, file = "./Output/Redes_Capas/Tablas/skin_opioid.rds")
genes_ambos <- unique(c(genes_skin_op, genes_oral_op))   # 36 - todos los de oral_op


# Construir grafo base
genes_tabla_merged <- bind_rows(
  lista_intersecciones_np[[25]]$genes_tabla,lista_intersecciones_np[[15]]$genes_tabla) %>%
  dplyr::filter(gene.ENSG %in% genes_ambos)

# Objeto compatible con build_network
interseccion_compartida <- list(
  genes_ensembl = genes_ambos,
  genes_tabla   = genes_tabla_merged)

g <- build_network(interseccion_compartida, aristas_color, colores_db)


# Añadimos seed genes
g <- add_layer_seed(g, interseccion_compartida$genes_tabla)

#--------------------------------------------------------
# AÑADIMOS CAPAS / ATRIBUTOS
#--------------------------------------------------------

# ---- GENES DIFERENCIALMENTE EXPRESADOS EN ENFERMEDADES - OTAR  ----

lista_genes_DEG_np <- readRDS("./Output/Piloto_Microbiota/Expresion_Diferencial_OTAR/df_genes_NP.rds")
#genes_DEG_np <- lista_genes_DEG_np$targetId
#intersect(genes_DEG_np, genes_oral_op)

g <- add_layer_DEG(g, lista_genes_DEG_np)

# ---- GENES QUE SON DRUG TARGETS ----

intersect(chmbl$targetId, genes_oral_op)
g <- add_layer_drug_target(g, chmbl)

# ---- GENES CON VARIANTES RARAS DE INTERES ----

# DF con variantes raras en lugares codificantes
variantes_raras_codificantes <- readRDS("./Data/Variantes/variantes_raras_codificantes.rds")

# FIltrar por variantes de enfermedades que nos interesan (áreas terapéuticas interesantes)
diseases_OTAR <- arrow::open_dataset("./Data/Diccionarios/disease.parquet") %>%  # Abrimos diseases.parquet
  dplyr::select(id, name, therapeuticAreas) %>%   # FIltramos
  collect()
diseases_mapping <- diseases_OTAR %>%     # Quitamos la lista de therapeuticAreas y ponemos trait - disease
  tidyr::unnest(therapeuticAreas) %>%
  dplyr::rename(therapeuticArea = therapeuticAreas)
saveRDS(diseases_mapping, file = "./Data/Diccionarios/diseases_mapping.rds")

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
message("Enfermedades en áreas de interés: ", length(ids_enfermedades_interes))


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


# ---- ESENCIALIDAD GENËTICA (CRISPR KO)
depmap_essential <- readRDS("./Data/CRISPR_KO/depmap_essential_por_grupo.rds")
crisprbrain_essential <- readRDS("./Data/CRISPR_KO/crisprbrain_essential.rds")

g <- add_layer_essential(g, depmap_essential, crisprbrain_essential)

#--------------------------------------------------------
# DIBUJAR GRAFO
#--------------------------------------------------------

#---- Resumen de atributos de nodo ----
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

# VEr esencialidad genética
data.frame(
  name                   = V(g)$name,
  symbol                 = V(g)$symbol,
  is_essential_neural    = V(g)$is_essential_neural,
  is_essential_digestive = V(g)$is_essential_digestive,
  is_essential_immune    = V(g)$is_essential_immune,
  is_essential_other     = V(g)$is_essential_other,
  is_essential_neuron    = V(g)$is_essential_neuron
) #%>%
  dplyr::summarise(
    n_neural    = sum(is_essential_neural,    na.rm = TRUE),
    n_digestive = sum(is_essential_digestive, na.rm = TRUE),
    n_immune    = sum(is_essential_immune,    na.rm = TRUE),
    n_other     = sum(is_essential_other,     na.rm = TRUE),
    n_neuron    = sum(is_essential_neuron,    na.rm = TRUE)
  )

data.frame(
  gene = V(g)$name,
  seed = (ifelse(V(g)$seed!=NA, 
  DEG = V(g)$DEG,
  drug_target = V(g)$drug_target,
  rare_variant = V(g)$has_rare_variant,
  essential_depmap = V(g)$depmap_essential,
  essential_brain = V(g)$crisprbrain_essential
)

# PLOT
plot <- plot_network(g, titulo = "Oral x Opioid Dependence", colores_db = colores_db,
                     use_qgraph_layout = TRUE,
                     remove_unannotated = TRUE)

cairo_pdf("./Output/Redes_Capas/A_Oral-Skin_Opioid_definitivo.pdf", width = 10, height = 7)
print(plot)
dev.off()

# Filtro STRING > 0,7

plot_filtro_string <- plot_network(g, titulo = "Oral x Opioid Dependence (String > 0,80)", colores_db = colores_db, 
                                   use_qgraph_layout = TRUE,
                                   string_score_min = 0.8, 
                                   remove_unannotated = TRUE)

cairo_pdf("./Output/Redes_Capas/A_Oral-Skin_Opioid_Filtro_STRING.pdf", width = 20, height = 15)
print(plot_filtro_string)
dev.off()

# FILTRO STRING > 0.7 + TISSUE EXPRESSION ATLAS

plot <- plot_network(g, titulo = "Oral x Opioid Dependence (String > 0,80)", colores_db = colores_db, 
                     use_qgraph_layout = TRUE,
                     string_score_min = 0.8, 
                     remove_unannotated = TRUE)

cairo_pdf("./Output/Redes_Capas/A_Oral-Skin_Opioid_Expression_Atlas.pdf", width = 25, height = 15)
print(plot)
dev.off()



# ANALISIS
atributos_nodos <- as.data.frame(vertex_attr(g))
map_variantes <- atributos_nodos %>%
  dplyr::filter(has_rare_variant == TRUE) %>%
  dplyr::select(name) %>%
  dplyr::inner_join(as.data.frame(variantes_raras_codificantes), by = c("name" = "targetFromSourceId")) %>%
  dplyr::select(Gen = name, Rasgo_ID = diseaseFromSourceMappedId) %>%
  dplyr::mutate(Capa = "Variante_Rara")

map_deg <- atributos_nodos %>%
  dplyr::filter(is_DEG == TRUE) %>% 
  dplyr::select(name) %>%
  dplyr::inner_join(lista_genes_DEG_np, by = c("name" = "targetId")) %>%
  dplyr::select(Gen = name, Rasgo_ID = diseaseId) %>% 
  dplyr::mutate(Capa = "DEG")

map_drugs <- atributos_nodos %>%
  dplyr::filter(is_drug_target == TRUE) %>% 
  dplyr::select(name) %>%
  dplyr::inner_join(chmbl, by = c("name" = "targetId")) %>%
  dplyr::select(Gen = name, Rasgo_ID = diseaseId) %>%
  dplyr::mutate(Capa = "Drug_Target")

# Combinamos los tres mapeos
mapeo_combinado <- dplyr::bind_rows(map_variantes, map_deg, map_drugs)

# Cruzamos con el diccionario maestro
matriz_enfermedades_completa <- mapeo_combinado %>%
  dplyr::inner_join(traits_MicroGWAS_areas, by = c("Rasgo_ID" = "Rasgo")) %>%
  dplyr::group_by(Rasgo_ID, name, Nombre_area, therapeuticAreas) %>%
  dplyr::summarise(
    Genes_DEG          = paste(unique(Gen[Capa == "DEG"]), collapse = ", "),
    Genes_DrugTarget   = paste(unique(Gen[Capa == "Drug_Target"]), collapse = ", "),
    Genes_VarRaras     = paste(unique(Gen[Capa == "Variante_Rara"]), collapse = ", "),
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

enfermedades_DEG <- matriz_enfermedades_completa %>%
  dplyr::filter(Num_DEG > 0) %>%
  dplyr::pull(Enfermedad_Nombre) %>%
  unique()

# Imprimir en consola
cat("\n--- ENFERMEDADES CON GENES DEG (", length(enfermedades_DEG), ") ---\n")
print(enfermedades_DEG)

enfermedades_DrugTargets <- matriz_enfermedades_completa %>%
  dplyr::filter(Num_DrugTarget > 0) %>%
  dplyr::pull(Enfermedad_Nombre) %>%
  unique()

# Imprimir en consola
cat("\n--- ENFERMEDADES CON DRUG TARGETS (", length(enfermedades_DrugTargets), ") ---\n")
print(enfermedades_DrugTargets)


enfermedades_VarRaras <- matriz_enfermedades_completa %>%
  dplyr::filter(Num_VarRaras > 0) %>%
  dplyr::pull(Enfermedad_Nombre) %>%
  unique()

# Imprimir en consola
cat("\n--- ENFERMEDADES CON VARIANTES RARAS (", length(enfermedades_VarRaras), ") ---\n")
print(enfermedades_VarRaras)

