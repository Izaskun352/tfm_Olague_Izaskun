
source("scripts/00_setup.R") 

# VISUALIZAR LOS DATOS MICROGWAS DEL GWAS CATALOG

# INPUTS

Allgene_microGWAS <- as.data.frame(readRDS("./Data/Microbioma/all.gene.gwas_microGWAS.rds"))
genes_por_rasgo_micro <- as.data.frame(table(Allgene_microGWAS[,3]))
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

all.gene.gwas_filter_GP_nonred <- readRDS("./Data./nasertic/input/all.gene.gwas_filter_GP_nonred.rds")
genes_por_rasgo_GP <- as.data.frame(table (all.gene.gwas_filter_GP_nonred[,3]))



table(Allgene_microGWAS[,3])

ID_microGWAS <- c("EFO_0007753", "EFO_0007874", "EFO_0007883", "EFO_0011013", "EFO_0801228", "EFO_0801229")
colores_6_traits <- c(
  "EFO_0007753" = "#FFFF99",
  "EFO_0007874" = "#FFDAC1",
  "EFO_0007883" = "#AEC6CF",
  "EFO_0011013" = "#B5EAD7",
  "EFO_0801228" = "#C9B1FF",
  "EFO_0801229" = "#FFD1DC"
)

# -------------------------------
# BARPLOT ----

nombres_GWAS <- Allgene_microGWAS %>%
  left_join(traits_MicroGWAS_areas, by = c("disease" = "Rasgo")) %>%
  dplyr::select(gene, disease, name)
head(nombres_GWAS)

datos_grafico <- nombres_GWAS %>%
  group_by(name, disease) %>%
  summarise(total_genes = n(), .groups = 'drop') %>%
  # Ordenar para que el gráfico se vea organizado
  arrange(desc(total_genes))


# Dibujamos barplot

barplot <- ggplot(data = datos_grafico, aes(x = total_genes, y = reorder(name, total_genes), fill = disease)) +
  # Usamos identity y los anchos que definiste
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  
  # Etiquetas de texto con el valor exacto
  geom_text(aes(label = total_genes),
            position = position_dodge(width = 0.8),
            hjust = -0.3, # Un poco más separado para que no toque la barra
            size = 9,
            color = "black") +
  
  # Aplicar tus colores específicos para los 6 traits
  scale_fill_manual(values = colores_6_traits) +
  
  # Configuración del eje X (arriba)
  # pero puedes poner el que quieras
  scale_x_continuous(
    position = "top",
    expand = expansion(mult = c(0, 0.1)) # Da un 10% de espacio extra a la derecha para el texto
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "Genes GWAS por Rasgo",
    x = "Cantidad total de genes",
    y = "",
    fill = "ID del Rasgo"
  ) +
  
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30, color = "black"),
    axis.line.x.top = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x.top = element_line(color = "black"),
    axis.title.x.top = element_text(margin = margin(b = 15), size = 24, color = "black"),
    axis.text = element_text(color = "black", size = 20), # Nombres de rasgos legibles
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 30, color = "black"),
    plot.margin = margin(t = 10, r = 100, b = 10, l = 10) # Margen derecho para el geom_text
  )

pdf("./Output/Gráficos/MicroGWAS/Barplot_Rasgos_MicroGWAS.pdf", width = 20, height = 10)
print(barplot)
dev.off()


# -------------------------------


# -------------------------------
# VOLIN PLOT _ VARIACION COMUN vs MICROGWAS ----

datos_violin_plot = data.frame(Tipo_variante = c(rep("Rasgos variante común", length(genes_por_rasgo_GP[,1])), 
                                                 rep("Rasgos MicroGWAS", length(genes_por_rasgo_micro[,1]))),
                               Rasgos = c(genes_por_rasgo_GP[,1], genes_por_rasgo_micro[,1]),
                               distribucion_genes = as.numeric(c(genes_por_rasgo_GP[,2], genes_por_rasgo_micro[,2])))

violin <- ggplot(data = datos_violin_plot, aes(x = Tipo_variante, y = distribucion_genes, fill = Tipo_variante)) + 
  geom_violin(alpha = 0.5, color = NA) +  # trim: cortar colas , alpha: transparencia
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.4) + # añadir boxplot interno
  scale_fill_manual(values= c("Rasgos variante común" = "purple", "Rasgos MicroGWAS" = "orange")) + #
  labs (
    title = "Distribución del número de genes iniclaes por rasgo",
    x = "",
    y = "Número de genes asociados",
    fill = "Tipo de Rasgo"
  ) + 
  scale_y_log10(breaks = c(1, 5, 10, 50, 100, 200, 500, 1000, 2000, 3000)) + # escala logaritmica del eje Y pq los valores son muy bajos
  theme_bw() +  # este tema pone el borde completo (el cuadrado)
  theme(  
    panel.grid.major = element_blank(),   #quitamos las líneas del fondo
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(margin=margin(t = 0, r = 0, b = 10, l = 0),hjust = 0.5, face = "bold", size = 14),  # ponemos título en el centro
    axis.title.y = element_text(margin=margin(t = 0, r = 10, b = 0, l = 0),size = 12), # separamos el título del eje y
    axis.title.x = element_text(size = 12),
    legend.position = c(0.98, 0.98),  # ponemos la leyenda en estas coordenadas
    legend.justification = c("right", "top"),   # justificamos que la estamos poniendo arriba a la derecha
    legend.title = element_text(size = 10),
    aspect.ratio = 0.8,    # hacemos que el gráfico sea más estrecho (obliga a que sea más alto que estrecho)
    plot.margin = margin(t = 30, r = 10, b = 20, l = 10)
  )

pdf("./Output/Gráficos/Graficos_Resultados_1/ViolinPlot_Comun_Micro.pdf", width = 6, height = 7)
print(violin)
dev.off()
# -------------------------------


