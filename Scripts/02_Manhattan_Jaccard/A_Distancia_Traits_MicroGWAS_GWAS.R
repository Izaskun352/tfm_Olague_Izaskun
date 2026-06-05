

source("scripts/00_setup.R")  # Abrimos script con las librerias
source("scripts/03_Funciones_Propias.R")
source("scripts/04_Funciones_Graficos.R")# Abrimos script con las funciones
source("scripts/05_Funciones_Anotar_ID.R")

#----------------------------------------


## DISTANCIA DE MANHATTAN DE MICROGWAS CON EL RESTO DE TRAITS CON VARIACION COMUN


#--------------------------------------
# INPUTS
#--------------------------------------

carpeta_traits_GWAS <- "./Data/nasertic/output/RDS_gwas"
carpeta_solo_MicroGWAS <- "./Output/MicroGWAS"

All_diseases <- readRDS("./Data./nasertic/input/All_diseases.rds")
Diccionario_TherapeuticAreas <- readRDS("./Data./Diccionarios/Diccionario_TherapeuticAreas.rds")

ID_microGWAS <- c("EFO_0007753", "EFO_0007874", "EFO_0007883", "EFO_0011013", "EFO_0801228", "EFO_0801229")

# --------------------------------------------------------
# OUTPUTS
# -------------------------------------------------------

carpeta_MicroGWAS_GWAS <- "./Output/Piloto_Microbiota/MicroGWAS_GWAS"
Matriz_Distancia_MicroGWAS_GWAS_Sin_filtrar <- readRDS("./Data/Matrices_Distancias/Matriz_Distancia_MicroGWAS_GWAS.rds")
dim(Matriz_Distancia_MicroGWAS_GWAS)

Matriz_Distancia_MicroGWAS_GWAS <- readRDS("./Data/Matrices_Distancias/Matriz_Distancia_MicroGWAS_GWAS_1832_traits.rds")

traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# --------------------------------------------------------
# CALCULAMOS LA DISTANCIA DE MANHATTAN
# -------------------------------------------------------

## Creamos carpeta con todos los archivos

carpeta_MicroGWAS_GWAS <- "./Output/Piloto_Microbiota/MicroGWAS_GWAS"

dir.create(carpeta_MicroGWAS_GWAS, showWarnings = FALSE, recursive = TRUE)

archivos_1 <- list.files(path = carpeta_traits_GWAS, full.names = TRUE)
archivos_2 <- list.files(path = carpeta_solo_MicroGWAS, full.names = TRUE)
todos_los_archivos <- c(archivos_1, archivos_2)

resultado <- file.copy(from = todos_los_archivos, to = carpeta_MicroGWAS_GWAS, overwrite = TRUE)
if(all(resultado)) {
  print("¡Todos los archivos se copiaron con éxito!")
} else {
  print("Algunos archivos no se pudieron copiar. Verifica las rutas.")
}

## Extraemos archivos de la carpeta (matrices)
archivos_carpeta <- list.files(carpeta_MicroGWAS_GWAS,              
                               pattern = "\\.rds$", 
                               full.names = TRUE)
# Extraemos vector con los scores de todos los genes para cada trait
lista_vectores_score_comun <- extraer_vectores_scores(archivos_carpeta = archivos_carpeta,     
                                                      tipo_extraccion = "columna",
                                                      columna_extraer = 4)

## Limpieza de nombres
names(lista_vectores_score_comun) <- gsub("ZSCO.", "", names(lista_vectores_score_comun))
names(lista_vectores_score_comun) <- gsub("\\.rds$", "", names(lista_vectores_score_comun))
names(lista_vectores_score_comun) <- trimws(names(lista_vectores_score_comun))

## Creamos la matriz (Ahora las filas tendrán los nombres limpios)

matriz_vectores <- do.call(rbind, lista_vectores_score_comun)   # Creamos una matriz donde cada fila es un vector de la lista
message("Tras convertir los vectores en matriz, la matriz tiene una dimension de ", dim(matriz_vectores)[1], " x ", dim(matriz_vectores)[2])  # comprobamos las dimensiones

matriz_distancias <- as.matrix(dist(matriz_vectores, method = "manhattan"))  # creamos matriz con las distancias entre todos los traits (pairwise)
message("La matriz de distancia de manhattan tiene una dimension de ", dim(matriz_distancias)[1], " x ", dim(matriz_distancias)[2])  # comprobamos las dimensiones + vemos que la diagonal es 0

#---------------------------------------------------
# ACTUALIZAR NOMBRE/ID TRAITS + ELIMINAR DUPLICADOS
#---------------------------------------------------

traits_MicroGWAS <- rownames(Matriz_Distancia_MicroGWAS_GWAS)
duplicados_MicroGWAS <- encontrar_duplicados(traits_MicroGWAS) # hay 3 trait duplicados --> 3 TRAIT SMicroGWAS ya estaban en los traits GWAS!!!!
print(duplicados_MicroGWAS)

for (id in duplicados_MicroGWAS) {
  
  # Extraemos las dos filas gemelas de la matriz
  filas_gemelas <- Matriz_Distancia_MicroGWAS_GWAS[rownames(Matriz_Distancia_MicroGWAS_GWAS) == id, ]
  
  # Aislamos las distancias de cada versión
  distancias_micro <- as.numeric(filas_gemelas[1, ])
  distancias_gwas <- as.numeric(filas_gemelas[2, ])
  
  # Calculamos la correlación de Pearson (1 = son clones perfectos)
  correlacion <- cor(distancias_micro, distancias_gwas, use = "complete.obs")
  
  # Imprimimos el veredicto para cada ID
  print(paste("ID:", id))
  print(paste("-> Correlación entre MicroGWAS y GWAS:", round(correlacion, 4)))
  
  # Hacemos una comprobación de diferencias exactas por si acaso
  diferencia_media <- mean(abs(distancias_micro - distancias_gwas), na.rm = TRUE)
  if (diferencia_media == 0) {
    print("-> ¡Son matemáticamente idénticas (Diferencia = 0)!")
  } else {
    print(paste("-> Hay una diferencia media en las distancias de:", round(diferencia_media, 4)))
  }
  print("--------------------------------------------------")
}   # la correlacion es bastente alta pero aun así nos quedamos con los score calculados con la propagación posterior (carpeta microGWAS)

Matriz_Distancia_MicroGWAS_GWAS <- Matriz_Distancia_MicroGWAS_GWAS[!duplicated(rownames(Matriz_Distancia_MicroGWAS_GWAS)), ]
Matriz_Distancia_MicroGWAS_GWAS <- Matriz_Distancia_MicroGWAS_GWAS[, !duplicated(colnames(Matriz_Distancia_MicroGWAS_GWAS))]
dim(Matriz_Distancia_MicroGWAS_GWAS)

traits_MicroGWAS <- rownames(Matriz_Distancia_MicroGWAS_GWAS)
traits_MicroGWAS_areas <- data.frame(Rasgo = traits_MicroGWAS) %>%
  
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

traits_NA_microGWAS <- traits_MicroGWAS_areas %>% filter(is.na(name))   %>% pull(Rasgo)  # hay 107 rasgos NA


## Buscamos ID obsoletos en disease.parquet
ruta_archivo_nombre_disease <- "./Data/Diccionarios/disease.parquet"
diccionario_diseases <- read_parquet(ruta_archivo_nombre_disease)
id_parquet <- buscar_diccionario_parquet(traits_NA_microGWAS, diccionario_diseases)   # se han recuperado 53
traits_MicroGWAS_areas <- actualizar_metadatos_df(traits_MicroGWAS_areas, id_parquet)   # sustituimos por el nuevo ID del rasgo y anotamos name + área terapéutica 

traits_NA_microGWAS <- traits_MicroGWAS_areas %>% filter(is.na(name))   %>% pull(Rasgo) # ahora hay 54 NA

## Limpiar de duplicados + restos de ID antiguos

traits_MicroGWAS_areas <- traits_MicroGWAS_areas %>%
  
  filter(!Rasgo %in% id_parquet$id_antiguo) %>% # eliminamos las filas que contengan id viejos
  
  distinct(Rasgo, .keep_all = TRUE) # eliminamos filas duplicadas

print(paste("Dimensiones del dataframe limpio:", nrow(traits_MicroGWAS_areas), "filas."))   # nos quedamos con 1828 traits


# Sustituimos el Nombre_area

dict_areas <- setNames(Diccionario_TherapeuticAreas$Nombre_area, # necesario para optimizar el código y poder utilizar función ids
                       Diccionario_TherapeuticAreas$ID_area)
traits_MicroGWAS_areas <- traits_MicroGWAS_areas %>%
  mutate(
    Nombre_area = case_when(
      # Condición: Solo actuamos si Nombre_area es NA
      Nombre_area == "NA" ~ sapply(strsplit(as.character(therapeuticAreas), ";\\s*"), function(ids) {
        # ids es un vector con los IDs separados. Los buscamos en el diccionario
        nombres_encontrados <- dict_areas[ids]
        
        # Juntamos los nombres encontrados con ";" omitiendo los que no existan (NA)
        paste(unique(na.omit(nombres_encontrados)), collapse = ";")
      }),
      
      # Si no era NA, dejamos el valor que ya tenía
      TRUE ~ Nombre_area
    )
  ) %>%
  # Limpieza opcional: Si el paste dejó cadenas vacías (""), las volvemos NA reales
  mutate(Nombre_area = na_if(Nombre_area, ""))

saveRDS(traits_MicroGWAS_areas, file = "./Data/Diccionarios/traits_MicroGWAS_areas.rds")
traits_MicroGWAS_areas <- readRDS("./Data/Diccionarios/traits_MicroGWAS_areas.rds")

# -----------------------------------------------------------------
## Sustituimos nombre de los nuevos ID en la matriz
# -----------------------------------------------------------------
Matriz_Distancia_MicroGWAS_GWAS <- actualizar_dimnames_matriz(
  Matriz_Distancia_MicroGWAS_GWAS,
  vector_antiguos = id_parquet$id_antiguo,
  vector_nuevos = id_parquet$id_nuevo
)

## Limpiar la matriz de duplicaciones y de posibles residuos de ID antiguos

indices_validos <- !(rownames(Matriz_Distancia_MicroGWAS_GWAS) %in% id_parquet$id_antiguo) # Comprobamos qué nombres de fila NO son ID antiguos
Matriz_Distancia_MicroGWAS_GWAS <- Matriz_Distancia_MicroGWAS_GWAS[indices_validos, indices_validos]
print(paste("Dimensiones de la matriz con ID nuevos:", nrow(Matriz_Distancia_MicroGWAS_GWAS), "x", ncol(Matriz_Distancia_MicroGWAS_GWAS)))

indices_unicos <- !duplicated(rownames(Matriz_Distancia_MicroGWAS_GWAS)) # Comprobamos qué nombres de fila son únicos (no están duplicados) --> hay 15
Matriz_Distancia_MicroGWAS_GWAS <- Matriz_Distancia_MicroGWAS_GWAS[indices_unicos, indices_unicos]
print(paste("Dimensiones de la matriz limpia:", nrow(Matriz_Distancia_MicroGWAS_GWAS), "x", ncol(Matriz_Distancia_MicroGWAS_GWAS)))

saveRDS(Matriz_Distancia_MicroGWAS_GWAS, file = "./Data/Matrices_Distancias/Matriz_Distancia_MicroGWAS_GWAS_1832_traits.rds")

# --------------------------------------------------------------
## Sustituimos los ID nuevos en la carpeta con los archivos
# --------------------------------------------------------------

# 1. Define la ruta de la carpeta donde están tus archivos .rds
# (Asegúrate de usar barras inclinadas '/' o dobles invertidas '\\')
ruta_carpeta <- "./Output/Piloto_Microbiota/MicroGWAS_GWAS" 
archivos_actuales <- list.files(path = ruta_carpeta, pattern = "\\.rds$", full.names = TRUE)

registro_cambios <- c()
# 4. Bucle para buscar y reemplazar los nombres
for (i in 1:nrow(id_parquet)) {
  old_id <- trimws(as.character(id_parquet$id_antiguo[i]))
  new_id <- trimws(as.character(id_parquet$id_nuevo[i]))
  
  if(old_id == "" || is.na(old_id)) next
  
  # Encuentra qué archivos contienen el ID obsoleto actual
  archivos_a_cambiar <- grep(old_id, archivos_actuales, value = TRUE)
  
  # Si encuentra coincidencias, procede a renombrar
  if (length(archivos_a_cambiar) > 0) {
    for (archivo in archivos_a_cambiar) {
      
      # Genera el nuevo nombre sustituyendo el ID viejo por el nuevo
      # gsub buscará el 'old_id' en cualquier parte de la cadena y lo cambiará
      nuevo_nombre <- gsub(old_id, new_id, archivo)
      
      # Ejecuta el renombrado físico en tu disco duro
      file.rename(from = archivo, to = nuevo_nombre)
      
      mensaje <- paste(" ✔️", basename(archivo), "-->", basename(nuevo_nombre))
      registro_cambios <- c(registro_cambios, mensaje)
      
      # Actualiza nuestra lista interna para reflejar el cambio en esta iteración
      archivos_actuales[archivos_actuales == archivo] <- nuevo_nombre
    }
  }
}

if (length(registro_cambios) > 0) {
  cat("\n¡Renombrado completado con éxito! Aquí tienes el resumen:\n")
  cat(paste(registro_cambios, collapse = "\n"))
  cat("\n\nTotal de archivos modificados:", length(registro_cambios), "\n")
} else {
  cat("\nNo se encontró ningún archivo que coincidiera con los IDs obsoletos proporcionados en el dataframe.\n")
}
                                                                                          
#-------------------------------------------------------
####  DENDROGRAMA
#-------------------------------------------------------

hc_microGWAS <- hclust(as.dist(Matriz_Distancia_MicroGWAS_GWAS), method = "complete")  # as.dist: pq el primer argumento tiene que estar producido por dist
plot(hc_comun)

# elegimos method complete, también podría ser average / ward.D2 
clusters_microGWAS <- cutree(hc_microGWAS, h = 0.8) 
table(clusters_microGWAS) #View(clusters_comun)  ; plot(hc_comun) --> ilegible

## saber en qué cluster están los traits microGWAS  -------> están todos en el Cluster 1 (551 traits)

clusters_de_microGWAS <- clusters_microGWAS[ID_microGWAS]

tabla_clusters_microGWAS <- data.frame(
  ID_Enfermedad = names(clusters_de_microGWAS),
  nombre = traits_MicroGWAS_areas[id_parquet],
  Cluster_Asignado = as.integer(clusters_de_microGWAS)
)


## ponemos el nombre real de la enfermedad
hc_microGWAS_nombres <- hc_microGWAS
hc_microGWAS_nombres$labels <- traits_MicroGWAS_areas$name[match(hc_microGWAS_nombres$labels, traits_MicroGWAS_areas$Rasgo)]

pdf("./Output/Gráficos/MicroGWAS/Dendrograma_Manhatttan_GWAS_MicorGWAS.pdf", height = 10, width = 50)

plot(hc_microGWAS_nombres,
     hang = -1,        # Obliga a que todas las hojas bajen hasta la línea base (el fondo)
     cex = 0.2,        # Reduce el tamaño de la letra 
     main = "Distancia de Manhattan", 
     xlab = "ID Enfermedades", 
     ylab = "Distancia")
#rect.hclust(hc_microGWAS_nombres, h = 0.75, border = c("red", "blue", "green", "purple", "orange")) # añade rectánculos alrededor de los grupos cutree

dev.off()

## coloreamos los traits MicroGWAS

library(dendextend)
dendrograma_microGWAS <- as.dendrogram(hc_microGWAS, hang = 0.1) %>%
  color_labels(labels = ID_microGWAS, col = "red")

dendrograma_microGWAS <- branches_attr_by_labels(  # Coloreamos las ramas micoGWAS
  dendrograma_microGWAS,
  labels = ID_microGWAS,
  TF_values = c("red", "black"), # "red" si coincide , "black" para el resto
  attr = "col"                   # Atributo que queremos cambiar: el color
)

nombres_traits <- as.character(traits_MicroGWAS_areas$name[match(labels(dendrograma_microGWAS), traits_MicroGWAS_areas$Rasgo)])
dendrograma_microGWAS <- set(dendrograma_microGWAS, "labels", nombres_traits)
dendrograma_microGWAS <- set(dendrograma_microGWAS, "labels_cex", 0.2) # Aquí aplicamos tu cex = 0.2 a todo


pdf("./Output./Gráficos/MicroGWAS/Dendrograma_Manhatttan_GWAS_MicorGWAS.pdf", height = 10, width = 50)
plot(dendrograma_microGWAS,
     main = "Distancia de Manhattan", 
     xlab = "ID Enfermedades", 
     ylab = "Distancia")
# Añadir las cajas delimitadoras para los clústeres
rect.dendrogram(dendrograma_microGWAS, 
                h = 0.8, 
                border = c("dodgerblue", "darkorange", "forestgreen", "darkorchid"),   # Color de la caja 
                lty = 2,           # Estilo de línea: 2 es línea discontinua
                lwd = 1.5)         # Grosor de la línea de la caja
dev.off()


#-------------
# PHEATMAP 
#-------------

## Unimos cluster al trait + área terapeutica
df_clusters_areas_microGWAS <- generar_matriz_clusters_areas(vector_clusters = clusters_microGWAS,
                                                         df_traits_areas = traits_MicroGWAS_areas)

numero_clusters <- length(unique(clusters_microGWAS))

## Creamos matriz con la distancia entre clusters

matriz_distancia_clusters <- matrix(0, nrow = numero_clusters, ncol = numero_clusters)
rownames(matriz_distancia_clusters) <- paste("Clúster", 1:numero_clusters)
colnames(matriz_distancia_clusters) <- paste("Clúster", 1:numero_clusters)

for (i in 1:numero_clusters) {
  for (j in 1:numero_clusters) {
    # Sacamos los IDs (rasgos) que pertenecen al clúster i y al clúster j
    rasgos_i <- names(clusters_microGWAS[clusters_microGWAS == i])
    rasgos_j <- names(clusters_microGWAS[clusters_microGWAS == j])
    
    # Calculamos la distancia media de todos los cruces posibles entre esos dos grupos
    distancia_media <- mean(Matriz_Distancia_MicroGWAS_GWAS[rasgos_i, rasgos_j])
    
    # Lo guardamos en nuestra nueva minimatriz
    matriz_distancia_clusters[i, j] <- distancia_media
  }
}

## Creamos la anotación de la barra
anotaciones_df_microGWAS <- df_clusters_areas_microGWAS %>%
  separate_rows(Areas_del_Cluster, sep = ";") %>%
  
  mutate(Areas_del_Cluster = str_trim(Areas_del_Cluster)) %>%
  
  filter(Areas_del_Cluster != "" & !is.na(Areas_del_Cluster)) %>%
  
  group_by(ClusterID, Areas_del_Cluster) %>%
  summarise(Conteo = n(), .groups = "drop") %>%
  
  group_by(ClusterID) %>%
  slice_max(order_by = Conteo, n = 1, with_ties = FALSE) %>%   # cogemos el que tiene el número más alto en cada cluster
  
  ungroup() %>%
  
  select(-Conteo)

anotaciones_df_microGWAS <- as.data.frame(anotaciones_df_microGWAS)
colnames(anotaciones_df_microGWAS)[2] <- "Área Terapéutica"   # para que salga en la leyenda asi (cambiar nombre columna)
rownames(anotaciones_df_microGWAS) <- paste("Clúster", anotaciones_df_microGWAS$ClusterID)
anotaciones_df_microGWAS$ClusterID <- NULL

## Le damos un color a cada área
colores_areas <- generar_paleta_nombrada_colorRampPalette(vector_categorias = anotaciones_df_microGWAS$`Área Terapéutica`)


## Creamos el pheatmap
pdf("./Output/Gráficos/MicroGWAS/Distancia_Enfermedades_Heatmap_MicroGWAS_GWAS.pdf", height = 20, width = 20)
print(pheatmap(
  matriz_distancia_clusters, 
  color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  # breaks = seq(0.5, 1, length.out = 101),
  clustering_distance_rows = as.dist(matriz_distancia_clusters),
  clustering_distance_cols = as.dist(matriz_distancia_clusters),
  clustering_method = "complete",
  
  treeheight_col = 0,  # quitamos el dendrograma de las columnas
  treeheight_row = 150,  # ajustamos tamaño del dendrograma de las filas
  # La Barra Lateral  
  annotation_row = anotaciones_df_microGWAS,
  annotation_colors = list("Área Terapéutica" = colores_areas),
  fontsize = 8,
  annotation_names_row = FALSE,
  
  cellheight =20, 
  cellwidth = 10,
  
  show_rownames = TRUE,   # Oculta los 200 nombres de las enfermedades (ilegibles)
  show_colnames = TRUE,
  
  display_numbers = FALSE,
  
  main = paste0("Distancia media entre ", numero_clusters," clusters de enfermedades con variación común (MicroGWAS)"),
  border_color = NA, # Quita bordes para que se vea más limpio
  
))
dev.off()


#--------------------------------------------------
# GRAFO DE REDES--> IGRAPH  #
#--------------------------------------------------

## Vemos distribución de valores
valores <- Matriz_Distancia_MicroGWAS_GWAS[upper.tri(Matriz_Distancia_MicroGWAS_GWAS)]
boxplot(valores,
        boxwex = 0.3,
        main = "Distribución distancia de Manhattan entre enfermedades asociadas a variación común",
        ylab = "Distancia")
summary(valores)

## Hacemos umbral

umbral_p01_microGWAS <- quantile(Matriz_Distancia_MicroGWAS_GWAS, 0.01)

## Hacemos matriz de adyacencia
ºmatriz_adyacencia_microGWAS <- ifelse(Matriz_Distancia_MicroGWAS_GWAS < umbral_p01_microGWAS, Matriz_Distancia_MicroGWAS_GWAS, 0)

## Hacemos la red
red_manhattan_microGWAS <- graph_from_adjacency_matrix(matriz_adyacencia_microGWAS,
                                                  mode = "undirected",
                                                  weighted = TRUE,
                                                  diag = FALSE)

## Añadimos las áreas terapéuticas y el nombre de los traits como atributos / anotaciones
red_manhattan_microGWAS <- anotar_red_igraph(
  red_igraph = red_manhattan_microGWAS, 
  df_anotaciones = traits_MicroGWAS_areas[, c("Rasgo","name", "Nombre_area")],  # los atributos que queremos añadir
  columna_id = "Rasgo" # Le indicamos cómo se llama la columna de unión
)

vertex_attr_names(red_manhattan_microGWAS) # Ver todos los atributos que tiene la nueva red

saveRDS(red_manhattan_microGWAS, file = "./Output/Piloto_Microbiota/Red_Manhattan_microGWAS.rds")
red_manhattan_microGWAS <- readRDS("./Output/Piloto_Microbiota/Red_Manhattan_microGWAS.rds")

## Hacemos subgrupos

componentes_red_manhattan_microGWAS <- components(red_manhattan_microGWAS)
## hacemos lista con los subgrupos   --> en este caso no utilizamos esto
lista_subredes <- decompose(red_manhattan_microGWAS, min.vertices = 4)

# ------------------------------------------------------------------
## Random.walktrap para descomponer la red en módulos/comunidades
# ------------------------------------------------------------------

## Creamos los nodos pie

V(red_manhattan_microGWAS)$shape <- "pie"

# Generamos paleta
areas_brutas <- V(red_manhattan_microGWAS)$Nombre_area
areas_separadas <- unlist(strsplit(areas_brutas, ";\\s*"))
areas_totales <- unique(na.omit(areas_separadas))
paleta_maestra<- generar_paleta_nombrada_grande(areas_totales)

red_manhattan_microGWAS <- preparar_nodos_pie(subred = red_manhattan_microGWAS,
                                              paleta_colores = paleta_maestra,
                                              atributo_areas = "Nombre_area")


# Aplicar Walktrap a toda la red
comunidades_walktrap <- cluster_walktrap(red_manhattan_microGWAS)

# Ver cuántas comunidades se han formado y qué tamaño tienen
tamaño <- sizes(comunidades_walktrap)
comunidades_relevantes <- which(tamaño >= 10)
View(comunidades_relevantes)
cat("Se van a generar", length(comunidades_relevantes), "gráficos individuales.\n")

# Dibujamos cada red para cada comunidad --> Bucle para extraer y dibujar CADA comunidad por separado

ruta_carpeta <- "./Output/Gráficos/MicroGWAS/Red_Interaccion/"

for (id_comunidad in comunidades_relevantes) {
  
  # Extraer nodos de esa comunidad
  nodos_comunidad <- V(red_manhattan_microGWAS)[membership(comunidades_walktrap) == id_comunidad]
  
  # Crear la subred
  subred_individual <- induced_subgraph(red_manhattan_microGWAS, vids = nodos_comunidad)
  
  # Crear nombre del archivo
  nombre_archivo <- paste0(ruta_carpeta, "Comunidad_", id_comunidad, ".pdf")
  
  pdf(file = nombre_archivo, width = 10, height = 7)
  
  dibujar_red_desplazada(
    red_igraph = subred_individual, 
    titulo_grafico = paste("Distancia de Manhattan: Layout_kk; Umbral q1 (Variante comun) Comunidad", id_comunidad, "|", vcount(subred_individual), "rasgos"),
    layout = layout_with_kk(subred_individual),
    atributo_etiqueta = "name"
  ) 
  
  añadir_leyenda_pleiotropia(red_igraph = subred_individual,    ## Añadir leyenda
                             atributo_leyenda = "Nombre_area",
                             paleta_global = paleta_maestra,
                             titulo_leyenda = "Áreas Terapéuticas")
  
  dev.off()
}
cat("¡Proceso terminado! Revisa tu carpeta para ver los PDFs.\n")


# --------------------------------------------------------------
# Ver donde están los rasgos MicroGWAS


rasgos_presentes <- ID_microGWAS[ID_microGWAS %in% V(red_manhattan_microGWAS)$Rasgo]

# Si falta alguno, R te lo dirá:
rasgos_perdidos <- setdiff(ID_microGWAS, rasgos_presentes)
if(length(rasgos_perdidos) > 0) {
  cat("Ojo, estos rasgos no pasaron tu filtro inicial y no están en la red:\n")
  print(rasgos_perdidos)
}

# Buscamos en qué comunidad está cada uno de los que sí sobrevivieron
# Extraemos la membresía global
membresia_global <- membership(comunidades_walktrap)
V(red_manhattan_microGWAS)$Comunidad <- as.numeric(membership(comunidades_walktrap))

# Usamos los nombres para buscar su número de comunidad
posiciones_microGWAS <- match(rasgos_presentes, V(red_manhattan_microGWAS)$Rasgo)
print("Posiciones en la red:")
print(posiciones_microGWAS)

comunidades_con_microGWAS <- membresia_global[posiciones_microGWAS]

Nombre = as.character(vertex_attr(red_manhattan_microGWAS, "name", index = posiciones_microGWAS))
# Lo mostramos en una tabla limpia
resultado_comunidades_microGWAS <- data.frame(
  Rasgo = unlist(rasgos_presentes),
  Nombre = Nombre[posiciones_microGWAS],
  ID_Comunidad = as.numeric(comunidades_con_microGWAS)
)

print("Tus 6 rasgos están en las siguientes comunidades:")
print(resultado_comunidades_microGWAS)

ruta_carpeta <- "./Output/Gráficos/MicroGWAS/Red_Interaccion/"

for (id_comunidad in unique(comunidades_con_microGWAS)) {
  
  # Extraer nodos de esa comunidad
  nodos_comunidad <- V(red_manhattan_microGWAS)[membership(comunidades_walktrap) == id_comunidad]
  
  # Crear la subred
  subred_individual <- induced_subgraph(red_manhattan_microGWAS, vids = nodos_comunidad)
  
  # Crear nombre del archivo
  nombre_archivo <- paste0(ruta_carpeta, "Comunidad_", id_comunidad, ".pdf")
  
  pdf(file = nombre_archivo, width = 10, height = 7)
  
  dibujar_red_desplazada(
    red_igraph = subred_individual, 
    titulo_grafico = paste("Distancia de Manhattan: Layout_kk; Umbral q1 (Variante comun) Comunidad", id_comunidad, "|", vcount(subred_individual), "rasgos"),
    layout = layout_with_kk(subred_individual),
    atributo_etiqueta = "name"
  ) 
  
  añadir_leyenda_pleiotropia(red_igraph = subred_individual,    ## Añadir leyenda
                             atributo_leyenda = "Nombre_area",
                             paleta_global = paleta_maestra,
                             titulo_leyenda = "Áreas Terapéuticas")
  
  dev.off()
}






