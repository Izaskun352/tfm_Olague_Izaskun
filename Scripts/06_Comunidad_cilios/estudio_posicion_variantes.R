
# RED CLUSTER DYNEIN

red <- readRDS ("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/red_todos_interseccion.rds")
#falta poner bien el color de los genes semilla
#añadir los genes que tienen variación rara
#añadir nodos de a que enfermedad estan asociados
#añadir si son diana terapéutica (?)
#añadir DEG en enfermedades (?)
#añadir cuales estan en el KEGG pathway de Motor Proteins

#---ESTUDIO POSICION VARIANTES RARAS

df_clean = readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/ProtVar_clean.rds")
variantes_dineina <- readRDS("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/variantes_dineina.rds")

ProtVar_dineina_CADD30 <- readRDS( "./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/ProtVar_dineina_CADD30.rds")

dominios <- read_csv("./Output/Piloto_Microbiota/Red_Pleiotropia_MicroGWAS_VarComun/dynein_pathway/dominios_genes_dynein.csv")
df_CADD30_dominios <- ProtVar_dineina_CADD30 |>
  filter(!is.na(name)) |>
  inner_join(dominios, by = c("Gene" = "gene"),
             relationship = "many-to-many") |>
  filter(Amino_acid_position >= inicio & Amino_acid_position <= fin)


#agrupar variantes en menos de 20aa + mismo gen
df_CADD30_clusters <- ProtVar_dineina_CADD30 |>
  filter(!is.na(name)) |>
  arrange(Gene, name, Amino_acid_position) |>
  group_by(Gene, name) |>
  mutate(
    dist_anterior = Amino_acid_position - lag(Amino_acid_position, default = first(Amino_acid_position)),
    cluster_id = cumsum(dist_anterior > 20)
  ) |>
  ungroup()
resumen_clusters <- df_CADD30_clusters |>
  group_by(Gene, name, cluster_id) |>
  summarise(
    n_variantes = n(),
    pos_inicio  = min(Amino_acid_position),
    pos_fin     = max(Amino_acid_position),
    CADD_max    = max(CADD_phred_like_score, na.rm = TRUE),
    CADD_medio  = round(mean(CADD_phred_like_score, na.rm = TRUE), 1),
    .groups     = "drop"
  ) |>
  arrange(desc(n_variantes), desc(CADD_max))

resumen_clusters |>
  group_by(Gene) |>
  summarise(
    n_clusters        = n_distinct(cluster_id),
    n_variantes_total = sum(n_variantes),
    n_enfermedades    = n_distinct(name),
    CADD_max          = max(CADD_max)
  ) |>
  arrange(desc(n_variantes_total))
