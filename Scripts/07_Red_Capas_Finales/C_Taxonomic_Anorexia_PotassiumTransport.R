
source("scripts/00_setup.R")

#------------------------------------------------------------------------------------------------------
# TAXONOMIC  x ANOREXIA NERVIOSA --> POTASSIUM ION TRANSPORT / ACTION POTENTIAL
#------------------------------------------------------------------------------------------------------

# INPUTS
lista_intersecciones_np <- readRDS("./Output/Piloto_Microbiota/lista_intersecciones_NP_MicroGWAS.rds")
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# CARGAR LO COMÚN

load("Data/Redes/config_red.RData")               # colores_db, prioridad_db, add_layer
aristas_color <- readRDS("Data/Redes/aristas_color.rds")
source("Scripts/Scripts_Limpios/Red_Capas_Finales/01_Funciones_Redes.R")


#------------------------------------------------------------------------------------------------------
# TAXONOMIC x ANOREXIA
#------------------------------------------------------------------------------------------------------

genes_tax_anorexia <- lista_intersecciones_np[[7]]$genes_ensembl  # 15 genes
genes_tabla <- lista_intersecciones_np[[7]]$genes_tabla
write.csv2(genes_tabla, file = "./Output/Redes_Capas/Tablas/tax_anorexia.csv" )
saveRDS(genes_tabla, file = "./Output/Redes_Capas/Tablas/tax_anorexia.rds" )

# Construir grafo base

interseccion <- lista_intersecciones_np[[7]]
g <- build_network(interseccion, aristas_color, colores_db)
g <- add_layer_seed(g, interseccion$genes_tabla)
#------------------------------------------------------------------------------------------------------
# ESTUDIAR SOLAPAMIENTO con intersecciones Traduccion
#------------------------------------------------------------------------------------------------------

indices_interes <- c(1,2,3,4,5,6, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21)

nombres_interes <- names(lista_intersecciones_np)[indices_interes]

solapamientos <- lapply(nombres_interes, function(nombre) {
  intersect(lista_intersecciones_np[[nombre]]$genes_ensembl, genes_tax_anorexia)
})
names(solapamientos) <- nombres_interes  # <-- asignar nombres explícitamente

solapamientos_top <- Filter(function(x) length(x) >= 13, solapamientos)
genes_compartidos <- Reduce(intersect, solapamientos_top)

nombres_top <- names(solapamientos_top)
tabla_compartidos <- lapply(nombres_top, function(nombre_par) {
  inter <- lista_intersecciones_np[[nombre_par]]
  inter$genes_tabla %>%
    dplyr::filter(gene.ENSG %in% genes_compartidos) %>%
    dplyr::mutate(
      par_origen = nombre_par,
      trait_1    = inter$trait_1,
      trait_2    = inter$trait_2,
      jaccard    = inter$jaccard
    )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(gene.ENSG, gene.gene) %>%
  dplyr::summarise(
    es_semilla_en_alguno = any(es_semilla_en_alguno),
    n_pares_presente     = dplyr::n(),
    traits_2             = paste(unique(trait_2), collapse = ";"),
    .groups              = "drop"
  )

write.csv2(tabla_compartidos, file = "./Output/Redes_Capas/Tablas/megacluster_potassium_transport.csv" )
saveRDS(tabla_compartidos, file = "./Output/Redes_Capas/Tablas/megacluster_potassium_transport.rds" )

genes_compartidos <- Reduce(intersect, solapamientos_top)

# Construir grafo base

nombres_top <- names(solapamientos_top)

lista_para_unir <- c(
  lapply(lista_intersecciones_np[nombres_top], \(x) x$genes_tabla),
  list(lista_intersecciones_np[[7]]$genes_tabla))

genes_tabla_merged <- bind_rows(lista_para_unir) %>%
  dplyr::filter(gene.ENSG %in% genes_compartidos)

interseccion_compartida <- list(
  genes_ensembl = genes_compartidos,
  genes_tabla   = genes_tabla_merged)

g <- build_network(interseccion_compartida, aristas_color, colores_db)

# Añadimos seed genes
g <- add_layer_seed(g, interseccion_compartida$genes_tabla)

#--------------------------------------------------------
# AÑADIMOS CAPAS / ATRIBUTOS
#--------------------------------------------------------

# ---- GENES DIFERENCIALMENTE EXPRESADOS EN ENFERMEDADES - OTAR  ----

lista_genes_DEG_np <- readRDS("./Output/Piloto_Microbiota/Expresion_Diferencial_OTAR/df_genes_NP.rds")
genes_DEG_np <- lista_genes_DEG_np$targetId
intersect(genes_DEG_np, genes_tax_anorexia)

g <- add_layer_DEG(g, lista_genes_DEG_np)

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


# ---- GENES QUE SON DRUG TARGETS ----

intersect(chmbl$targetId, genes_tax_anorexia)
g <- add_layer_drug_target(g, chmbl)

# --- Genes drug targets de área específica

g <- add_layer_drug_target_area(g, chmbl, ids_enfermedades_interes)
vertex_attr_names(g)
V(g)$is_drug_target_area
V(g)$diseases_drug_area

# ---- ESENCIALIDAD GENËTICA (CRISPR KO)
depmap_essential <- readRDS("./Data/CRISPR_KO/depmap_essential_por_grupo.rds")
crisprbrain_essential <- readRDS("./Data/CRISPR_KO/crisprbrain_essential.rds")

g <- add_layer_essential(g, depmap_essential, crisprbrain_essential)
#--------------------------------------------------------
# DIBUJAR GRAFO
#--------------------------------------------------------

# Resumen de atributos de nodo
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

# ver esencialidad genética (CRISPR KO)  --> todos 0
data.frame(
  name                   = V(g)$name,
  symbol                 = V(g)$symbol,
  is_essential_neural    = V(g)$is_essential_neural,
  is_essential_digestive = V(g)$is_essential_digestive,
  is_essential_immune    = V(g)$is_essential_immune,
  is_essential_other     = V(g)$is_essential_other,
  is_essential_neuron    = V(g)$is_essential_neuron
) %>%
  dplyr::summarise(
    n_neural    = sum(is_essential_neural,    na.rm = TRUE),
    n_digestive = sum(is_essential_digestive, na.rm = TRUE),
    n_immune    = sum(is_essential_immune,    na.rm = TRUE),
    n_other     = sum(is_essential_other,     na.rm = TRUE),
    n_neuron    = sum(is_essential_neuron,    na.rm = TRUE)
  )

# -------------------------------------
# PLOT
#--------------------------------------

# SOlo Taxonomic x Anorexia
plot <- plot_network(g, titulo = "Taxonomic x Anorexia", colores_db = colores_db)

cairo_pdf("./Output/Redes_Capas/C_Tax_Anorexia.pdf", width = 9, height = 8)
print(plot)
dev.off()

vertex_attr_names(g)

# Supercluster con intersecciones traduccion

plot <- plot_network(g, titulo = "Supercluster Potassium Transport", colores_db = colores_db)

cairo_pdf("./Output/Redes_Capas/C_Potassium_Transport.pdf", width = 15, height = 15)
print(plot)
dev.off()


# Supercluster con tissue expression atlas

plot <- plot_network(g, titulo = "Supercluster Potassium Transport", colores_db = colores_db)

cairo_pdf("./Output/Redes_Capas/C_Potassium_Transport_Expression_Atlas.pdf", width = 8, height = 8)
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









