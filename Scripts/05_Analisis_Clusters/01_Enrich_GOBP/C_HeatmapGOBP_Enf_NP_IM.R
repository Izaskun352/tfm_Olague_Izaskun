
# HEATMAP MEGACLUSTERS NP / IM x MICRO 

# Enfermedad x GOBP megacluster

datos <- tribble(
  ~Disease,      ~GOBP,
  "Attenbtion deficit hyperactivity disorder",    "Cytoplasmic translation",
  "Obsessive compulsive disorder",    "Cytoplasmic translation",
  "Major depressive disorder",      "Cytoplasmic translation",
  "Substance abuse",      "Cytoplasmic translation",
  "Bipolar Disroder",      "Cytoplasmic translation",
  "Anorexia Nerviosa",     "Cytoplasmic translation",
  #"Gut Micr Measurement", "Cytoplasmic translation",
  #"Vaginal Micr Measurement", "Cytoplasmic translation",
  #"Oral Micr Measurement", "Cytoplasmic translation",
  "Ankylosin spondytis", "Cytoplasmic translation",
  "Peripheral neuropathy", "Epigenetic regulation of gene expression",
  "Autism spectrum disorder", "Epigenetic regulation of gene expression",
  "Insomnia", "Epigenetic regulation of gene expression",
  #"Oral Micr Measurement", "Epigenetic regulation of gene expression",
  #"Gut Micr Measurement", "Epigenetic regulation of gene expression",
  "Inflammatory bowel disease", "Epigenetic regulation of gene expression",
  "Systemic lupus erythematosis", "Epigenetic regulation of gene expression",
  "Seasonal allergic rhinitis", "Epigenetic regulation of gene expression",
  "Psoriasis", "Epigenetic regulation of gene expression",
  "Atopic eczema", "Epigenetic regulation of gene expression",
  "Opioid dependence", "Synaptic transmission, glutamatergic",
  #"Oral Micr Measurement",  "Synaptic transmission, glutamatergic",
  #"Skin MIcr Measurement",  "Synaptic transmission, glutamatergic",
  "Peripheral neuropathy", "Oxidative phosphorilation",
  #"Taxonomic Micr Measurement", "Oxidative phosphorilation",
  #"Taxonomic Micr Measurement", "Potassium ion transmembrane transport",
  "Anorexia nerviosa", "Potassium ion transmembrane transport"
)
matriz <- datos %>%
  mutate(valor = 1) %>%
  pivot_wider(
    names_from = Disease,
    values_from = valor,
    values_fill = 0
  )

# Pasar nombres de fila
matriz <- as.data.frame(matriz)
rownames(matriz) <- matriz$GOBP
matriz$GOBP <- NULL
matriz <- as.matrix(matriz)


pdf("./Output/Gráficos/Gráficos_Resultados_4/Heatmap_GOBP_Megaclusters_NP_IM.pdf", width = 20, height = 20)
pheatmap(
  matriz,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  cellwidth  = 18,
  cellheight = 25, 
  color = c("cornsilk2", "skyblue4"),
  border_color      = "white",
  angle_col = 45,
  fontsize_row      = 13
)
dev.off()

pdf("./Output/Gráficos/Gráficos_Resultados_4/Heatmap_GOBP_Megaclusters_NP_IM_SIN_Micro.pdf", width = 20, height = 20)
pheatmap(
  matriz,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  cellwidth  = 18,
  cellheight = 25, 
  color = c("cornsilk2", "skyblue4"),
  border_color      = "white",
  angle_col = 45,
  fontsize_row      = 13
)
dev.off()