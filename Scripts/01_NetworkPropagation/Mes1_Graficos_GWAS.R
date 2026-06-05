
## Abrimos script con las librerias

source("scripts/00_setup.R")

#------------------------------------------------------------
# GRÁFICO DE BARRAS HORIZONTAL  --> NÚMERO DE GENES Y RASGOS
#------------------------------------------------------------

datos_grafico_barras = data.frame(variante = c("Variante común", "Variante rara",
                                               "Variante común", "Variante rara"),
                                  medida = c("Nº de Genes", "Nº de Genes", 
                                             "Nº de Rasgos", "Nº de Rasgos"),
                                  valor = c(length(genes_comun), length(genes_raros),
                                     length(rasgos_comun), length(rasgos_raros)))


ggplot(data = datos_grafico_barras, aes(x=valor, y=medida, fill=variante)) +
  # identity: usa el número exacto de la tabla, 
  # width = 0.8 hace las barras más gordas. position_dodge al mismo valor las mantiene juntas.
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +                                                                  # position dodge : coloca las barras una al lado de otra
  labs(
    title = "Variante común vs variante rara",
    x = "Cantidad total",
    y="",
    fill = "Tipo de variante"
  )+
  geom_text(aes(label=valor),
            position = position_dodge(width = 0.8),
            hjust = -0.2,
            size = 9,
            color = "black") +
  
  coord_cartesian(clip = "off") + #  permite que el texto se dibuje fuera del área del gráfico 
  scale_x_continuous(
    position = "top",           # Pone la línea y números arriba
    limits = c(0, 11000),       # Fija el rango exacto que pediste
    expand = c(0, 0)            # espacio extra a los lados (opcional)
  ) +
  scale_fill_manual(values = c("Variante común"= "#B988b8",
                               "Variante rara" = "#FCF75E")) + 
  theme_minimal() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30, color = "black"),  # Estilo del título + posición en el centro
    axis.line.x.top = element_line(color = "black", linewidth = 0.5), # Dibuja la línea del eje X arriba
    axis.ticks.x.top = element_line(color = "black"),             # los ticks en la línea
    axis.title.x.top = element_text(margin =margin(b=15), size = 24, color = "black"), 
    axis.text = element_text(color = "black"),# Título eje x del top --> separar de abajp (b)
    # margin(t = top, r = right, b = bottom, l = left)
    legend.position = "bottom",  # ponemos la leyenda en la parte de abajo
    panel.grid.major = element_blank(),   #quitamos las líneas del fondo
    panel.grid.minor = element_blank(),
    text = element_text(size = 30), color = "black",
    
    plot.margin = margin(t = 10, r = 80, b = 10, l = 10) # damos margen a la derecha 
    # t=top, r=right, b=bottom, l=left --> en este caso aumentamo la derecha 'r'
  )

#------------------------------------------------------------
# VIOLIN PLOT --> CUANTOS GENES HAY POR RASGO
#------------------------------------------------------------

datos_violin_plot = data.frame(Tipo_variante = c(rep("Variante común", length(genes_por_rasgo_GP[,1])), 
                                                     rep("Variante rara", length(genes_por_rasgo_rare[,1]))),
                               Rasgos = c(genes_por_rasgo_GP[,1], genes_por_rasgo_rare[,1]),
                               distribucion_genes = as.numeric(c(genes_por_rasgo_GP[,2], genes_por_rasgo_rare[,2])))

ggplot(data = datos_violin_plot, aes(x = Tipo_variante, y = distribucion_genes, fill = Tipo_variante)) + 
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +  # trim: cortar colas , alpha: transparencia
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.4) + # añadir boxplot interno
  scale_fill_manual(values= c("Variante común" = "purple", "Variante rara" = "yellow")) + 
  labs (
    title = "Distribución del número de genes iniclaes por rasgo",
    x = "",
    y = "Número de genes asociados",
    fill = "Tipo de Variante"
  ) + 
  scale_y_log10(breaks = c(1, 5, 10, 50, 100, 200, 500, 1000)) + # escala logaritmica del eje Y pq los valores son muy bajos
  theme_bw() +  # este tema pone el borde completo (el cuadrado)
  theme(  
    panel.grid.major = element_blank(),   #quitamos las líneas del fondo
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(margin=margin(t = 0, r = 0, b = 10, l = 0),hjust = 0.5, face = "bold", size = 12),  # ponemos título en el centro
    axis.title.y = element_text(margin=margin(t = 0, r = 10, b = 0, l = 0),size = 10), # separamos el título del eje y
    
    legend.position = c(0.98, 0.98),  # ponemos la leyenda en estas coordenadas
    legend.justification = c("right", "top"),   # justificamos que la estamos poniendo arriba a la derecha
    legend.title = element_text(size = 10),
    aspect.ratio = 0.8,    # hacemos que el gráfico sea más estrecho (obliga a que sea más alto que estrecho)
    plot.margin = margin(t = 30, r = 10, b = 20, l = 10)
  )
    

#------------------------------------------------------------
# RASGOS CON VARIACIÓN COMÚN Y CON VARIACIÓN RARA
#------------------------------------------------------------

  ## DIAGRAMA DE VENN --> NO USAR!


datos_Venn = list(
  "Variante Común" = rasgos_comun,
  "Variante Rara" = rasgos_raros
)
length(datos_Venn[[2]])

n_interseccion <- length(rasgos_compartidos) 

ggVennDiagram(datos_Venn,
              label_alpha = 0,
              edge_size = 0.7, # tamaño de los bordes
              category.names = c("   Variante Común      ", "       Variante Rara   ")) +
  
  coord_flip() + # cambiar la orientación del gráfico
  
  scale_fill_gradient2(
    low = "darkolivegreen1", 
    high = "darkolivegreen4",
    #mid = "darkseagreen1",
    #midpoint = n_interseccion  
  ) +                          
  
  scale_x_continuous(expand = expansion(mult = 0.2)) +  # márgenes para que los nombres no se corten
  
  labs(
    title = "Rasgos asociados a variación común y variación rara",
    fill = "Nº Rasgos"
  ) + 
  
  theme(
    plot.title = element_text(margin=margin(t = 0, r = 0, b = 0, l = 0),hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom",
    aspect.ratio = 0.7,
    plot.margin = margin(t = 30, r = 10, b = 20, l = 10))

#------------------------------------------------------------
  ## DIAGRAMA DE BARRAS
#------------------------------------------------------------

datos_barras_rasgosCompartidos <- data.frame (variante = c("Variante común", "Variante rara" , "Variante común + rara"), 
                                              valor = c(length(rasgos_comun)-length(rasgos_compartidos), length(rasgos_raros)- length(rasgos_compartidos), length(rasgos_compartidos)))  

# El orden de los 'levels' se dibuja de ABAJO hacia ARRIBA en el eje Y.
datos_barras_rasgosCompartidos$variante <- factor(
  datos_barras_rasgosCompartidos$variante,
  levels = c("Variante común + rara", "Variante común", "Variante rara")
)

ggplot(data = datos_barras_rasgosCompartidos, aes(x = valor, y = variante, fill = variante)) +
  geom_bar(stat = "identity",
           position = "dodge",
           width = 0.8,
           ) +
  labs(
    title = "Solapamiento de rasgos con variante común\ny con variante rara",
     x = "Número de rasgos",
    y = "",
  ) +
  geom_text(aes(label = valor),
            hjust = -0.2,
            size = 4.5) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(
    position = "top", limits = c(0, 3500), expand = c(0.05,0)) +
  scale_fill_manual(values = c (
                                "Variante común" = "#B988b8",
                                "Variante rara" = "#FCF75E",
                                "Variante común + rara" = "grey")) + 
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),   #quitamos las líneas del fondo
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(t=10, r = 0, l = 0, b = 10), size = 16),  # Estilo del título + posición en el centro
    axis.line.x.top = element_line(color = "black", linewidth = 0.5), # Dibuja la línea del eje X arriba
    axis.ticks.x.top = element_line(color = "black"),             # los ticks en la línea
    axis.title.x.top = element_text(margin =margin(b=20, t = 10), size = 16),
    text = element_text(size = 18),
    plot.margin = margin(t = 10, r = 80, b = 10, l = 10)
  )

#------------------------------------------------------------
# DISTRIBUCIÓN ÁREAS TERAPÉUTICAS POR CADA VARIANTE
#------------------------------------------------------------
  
  #### VIOLIN PLOT -- VARIANTE COMÚN ####

    ## Preparamos datos
df_All_TherapeuticAreas_GP = as.data.frame(genes_por_rasgo_GP) %>%       # tiene que ser un dataframe
  mutate(
    therapeuticArea_counts = paste0("All (n = ", length(genes_por_rasgo_GP), ")")
  )

df_areas_GP_Violin_Plot = areas_terapeuticas_GP_ordenado %>%
  mutate(
    therapeuticArea_counts = paste0 (Nombre_area, " (n = ", n, ")")
  )

df_areas_GP_Violin_Plot = bind_rows(df_All_TherapeuticAreas_GP, df_areas_GP_Violin_Plot)# Datos necesarios para el violin plot
df_areas_GP_Violin_Plot$Numero_Genes <- 
  as.numeric(df_areas_GP_Violin_Plot$Numero_Genes)

     #unique(df_areas_GP_Violin_Plot$therapeuticArea_counts)

    ## Hacemos violin plot
ggplot(data = df_areas_GP_Violin_Plot, aes(x = Numero_Genes, y = therapeuticArea_counts)) +
  geom_violin (fill = "#B988b8", trim = FALSE, alpha = 0.8, width = 1.4, color = NA) +
  # añadir línea fina que cubra todo el rango
  stat_summary(fun.data = function(y){  # creamos la función que nos calcule los límites (rango) para cada área terapéutica (y)
    return (data.frame(
      y = median(y),  # la función necesita un punto medio de referencia a partir del cual calcular el rango
      ymin = min(y),  # mínimo
      ymax = max(y)  # máximo
    ))
  }, geom = "linerange", linewidth = 0.2, color = "black") +
  
  # añadir la línea con los cuartiles:
  stat_summary(fun.data = function(y){  # creamos la función que nos calcule los cuartiles para cada área terapéutica (y)
    return (data.frame(
      y = median(y),  # la función necesita un punto medio de referencia a partir del cual calcular los cuartiles
      ymin = quantile(y, 0.25),  # cuartil 0.25
      ymax = quantile(y, 0.75)  # cuartil 0.75
    ))
  }, geom = "linerange", linewidth = 0.5, color = "black") +   # queremos que se represente con una línea negra
  # añadir el punto rojo con la mediana
  stat_summary(fun = median,      # le decimos que la función es la mediana 
               geom = "point", size = 1, color = "red") +  # queremos que la geometria sea un punto
  labs (
    title = "Número de genes iniciales por rasgo, agrupados según área terapéutica\n (Variante común)",
    x = "Número de genes",
    y = ""
  ) +
  scale_x_log10(breaks = c(1, 5, 10, 20, 50, 200, 500, 2000)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 1, face = "bold", margin = margin(t=10, r = 0, l = 0, b = 10), size = 10),
    axis.title.x = element_text(margin = margin(t=10, r = 0, l = 0, b = 10), size = 9),
    panel.grid.major = element_blank(),   #quitamos las líneas del fondo
    panel.grid.minor = element_blank(),
    aspect.ratio = 1.4
  )


  #### VIOLIN PLOT -- VARIANTE RARA ####


    ## Preparamos datos 
df_All_TherapeuticAreas_RARE = as.data.frame(genes_por_rasgo_rare) %>%
  mutate(
    therapeuticArea_counts = paste0("All (n = ", length(genes_por_rasgo_rare), ")")
  )

df_areas_RARE_Violin_Plot = areas_terapeuticas_RARE_ordenado %>%
  mutate(
    therapeuticArea_counts = paste0 (Nombre_area, " (n = ", n, ")")
  )

df_areas_RARE_Violin_Plot = bind_rows(df_All_TherapeuticAreas_RARE, df_areas_RARE_Violin_Plot)  # Datos necesarios para el violin plot
df_areas_RARE_Violin_Plot$Numero_Genes <- 
  as.numeric(df_areas_RARE_Violin_Plot$Numero_Genes)

    ## Hacemos violon plot
ggplot(data = df_areas_RARE_Violin_Plot, aes(x = Numero_Genes, y = therapeuticArea_counts)) +
  geom_violin (fill = "#FCF75E", trim = FALSE, alpha = 0.8, width = 1.4, color = NA) +
  # añadir línea fina que cubra todo el rango
  stat_summary(fun.data = function(y){  # creamos la función que nos calcule los límites (rango) para cada área terapéutica (y)
    return (data.frame(
      y = median(y),  # la función necesita un punto medio de referencia a partir del cual calcular el rango
      ymin = min(y),  # mínimo
      ymax = max(y)  # máximo
    ))
  }, geom = "linerange", linewidth = 0.2, color = "black") +
  
  # añadir la línea con los cuartiles:
  stat_summary(fun.data = function(y){  # creamos la función que nos calcule los cuartiles para cada área terapéutica (y)
    return (data.frame(
      y = median(y),  # la función necesita un punto medio de referencia a partir del cual calcular los cuartiles
      ymin = quantile(y, 0.25),  # cuartil 0.25
      ymax = quantile(y, 0.75)  # cuartil 0.75
    ))
  }, geom = "linerange", linewidth = 0.5, color = "black") +   # queremos que se represente con una línea negra
  # añadir el punto rojo con la mediana
  stat_summary(fun = median,      # le decimos que la función es la mediana 
               geom = "point", size = 1, color = "red") +  # queremos que la geometria sea un punto
  labs (
    title = "Número de genes iniciales por rasgo, agrupados según área terapéutica\n (Variante rara)",
    x = "Número de genes",
    y = ""
  ) +
  scale_x_log10(breaks = c(1, 5, 10, 20, 50, 200, 500, 2000)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 1, face = "bold", margin = margin(t=10, r = 0, l = 0, b = 10), size = 10),
    axis.title.x = element_text(margin = margin(t=10, r = 0, l = 0, b = 10), size = 9),
    panel.grid.major = element_blank(),   #quitamos las líneas del fondo
    panel.grid.minor = element_blank(),
    aspect.ratio = 1.4
  )


#-------------------------------------------------------
#### DISTRIBUCIÓN ÁREAS TERAPÉUTICAS ---> BARPLOT - VARIANTE RARA
#-------------------------------------------------------

#### Hacemos matriz con la frecuencia de cada área ####
traits_GWAS_areas <- readRDS("./Data/Diccionarios/traits_GWAS_areas.rds")
traits_GWAS_areas <- data.frame(Rasgo = rasgos_raros) %>%
  
  left_join(        # Unimos la columna de área terapéutica a cada rasgo --> Como usamos left_join convertimos las matrices en dataframes
    as.data.frame(All_diseases),
    by = c("Rasgo" = "ID")                 
  ) %>%
  
  separate_rows(therapeuticAreas, sep = ";") %>% # Un mismo rasgo - más de un área terapéutica
  
  
  left_join(Diccionario_TherapeuticAreas, by = c("therapeuticAreas" = "ID_area")) %>%
  
  group_by(Rasgo, name) %>%
  
  summarise(
    Nombre_area = paste(unique(Nombre_area), collapse = ";"),
    therapeuticAreas = paste(unique(therapeuticAreas), collapse = ";"),
    .groups = "drop"
  )

traits_NA <- traits_GWAS_areas %>% filter(is.na(name))   %>% pull(Rasgo) 

traits_GWAS_frecuencia_areas <- traits_GWAS_areas %>%
  
  separate_rows(therapeuticAreas, Nombre_area, sep = ";") %>%
  
  filter(Nombre_area != "NA" & !is.na(Nombre_area)) %>%  # Eliminamos los NA
  
  add_count(therapeuticAreas, name = "Frecuencia") %>%
  select(therapeuticAreas, Nombre_area, Frecuencia) 

# Dibujamos barplot

pdf("./Output/Gráficos/Mes1/Traits_Rara_Frecuencia_Areas_barplot.pdf", width = 8, height = 6)

print(ggplot(data = traits_GWAS_frecuencia_areas, aes(x = Frecuencia, y = fct_reorder(Nombre_area, Frecuencia))) +
        geom_bar(stat = "identity", position = "dodge", width = 0.5, fill = "steelblue" )+     # identity: usa el número exacto de la tabla, 
        # position dodge : coloca las barras una al lado de otra
        labs(
          title = "Frecuencia de áreas terapéuticas de Enfermedades asociadas a variante rara",
          x = "Número de enfermedades",
          y="",
        )+
        geom_text(aes(label=Frecuencia),
                  position = position_dodge(width = 0.5),
                  hjust = -0.2,
                  size = 2.5) +
        
        coord_cartesian(clip = "off") + #  permite que el texto se dibuje fuera del área del gráfico 
        scale_x_continuous(
          position = "bottom",           # Pone la línea y números 
          limits = c(0, max(traits_GWAS_frecuencia_areas$Frecuencia)),       # Fija el rango exacto 
          expand = c(0, 0)            # espacio extra a los lados (opcional)
        ) +
        theme_minimal() + 
        theme(
          panel.grid.major = element_blank(),   #quitamos las líneas del fondo
          panel.grid.minor = element_blank(),
          plot.title = element_text(margin =margin(b = 10, t = 10), hjust = 0.4, face = "bold", size = 10),  # Estilo del título + posición en el centro
          axis.line.x.bottom  = element_line(color = "black", linewidth = 0.5), # Dibuja la línea del eje X arriba
          axis.ticks.x.bottom = element_line(color = "black"),             # los ticks en la línea
          axis.title.x.bottom = element_text(margin =margin(t=5), size = 8),  # Título eje x del top --> separar de abajp (b)
          legend.position = "none",  # ponemos la leyenda en la parte de abajo
          panel.grid.major.y = element_blank(), # Limpiar/quitar líneas horizontales del fondo
          text = element_text(size = 8),
          
          plot.margin = margin(t = 10, r = 40, b = 10, l = 10) # damos margen a la derecha 
        ))

dev.off()

  