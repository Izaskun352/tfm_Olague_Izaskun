

# INTERSECCIONES ENTRE PROTEINAS DEPENDIENTES DE TEJIDO

# LIBRERIAS

library(data.table)
library(tidyverse)
library(ggplot2)
library(hexbin)

# INPUTS

ruta_carpetas <- "./Data/Tissue_Expression_Atlas/association_scores/"
archivos_csv <- list.files(path = ruta_carpetas, pattern = "\\.csv$", full.names = TRUE)

nombres_tejidos <- basename(archivos_csv) %>% 
  str_replace("cohorts_combined_", "") %>% 
  str_replace("_avg_outer_prob.csv", "")

# OUTPUTS

df_conversion <- readRDS("./Data/Tissue_Expression_Atlas/conversion_df_ENSG.rds")
expression_tissues_CV <- readRDS("./Data//Tissue_Expression_Atlas/expression_CV.rds")

carpeta_data_tejidos_CV <-"./Data/Tissue_Expression_Atlas/Data_CV_Tejidos/"
carpeta_plots_tejidos <-"./Data/Tissue_Expression_Atlas/Plots_Tejidos/"
carpeta_data_tejidos_cv_filtrado <- "./Data/Tissue_Expression_Atlas/Data_CV_Tejidos_filtrado/"

# ---- Diccionario genes: simbolo - ENSG ----

conversion <- read.csv2("./Data/Tissue_Expression_Atlas/conversion_df.csv")

conversion_df <- conversion %>%
  separate_wider_delim(
    cols = `from_id.to_id`,   # columna que queremos separar
    delim = ",",                    # separador
    names = c("from_id", "to_id"),   # Nombres de las dos nuevas columnas
    too_few = "align_start"          # Por si acaso alguna fila viene incompleta
  ) %>%
  filter(str_detect(from_id, "^ENSG"))  # nos quedamos solo con los ID de ENSG (De la columna from_id)

# hay simbolos de genes con mas de un ENSG - nos quedamos con el ENSG con un num mas bajo (el primero que se descubrio, el 'estándar')

df_conversion <- conversion_df %>%
  distinct() %>%
  arrange(from_id) %>% # ordenamos por ID de menor a mayor
  
  group_by(to_id) %>%
  slice(1) %>%  # nos quedamos con la primera aparicion de cada uno
  ungroup()

any(duplicated(df_conversion$to_id)) # comprobar si sigue habiendo duplicados

saveRDS(df_conversion, file = "./Data/Tissue_Expression_Atlas/conversion_df_ENSG.rds")

# ---- CALCULO CV y MEDIA GLOBAL (TODOS LOS TEJIDOS) DE CADA INTERACCION ---
# Construimos dataframe base 

cat("Inicializando estructuras base usando el Cerebro...\n")
pos_cerebro <- which(nombres_tejidos == "brain_tumor") 
base_red <- fread(archivos_csv[pos_cerebro], select = c("prot1", "prot2"))  # utilizamos como base el cerebro pq tiene el mayor num de interacciones

acumulador <- data.table(   
  prot1 = base_red$prot1,
  prot2 = base_red$prot2,
  suma_scores = 0,
  suma_cuadrados = 0,
  N = length(archivos_csv) 
)
rm(base_red); gc()

# Buycle para ir sumando tejido a tejido
for(i in 1:length(archivos_csv)) {
  cat("Procesando estadísticas para el tejido:", nombres_tejidos[i], "\n")
  
  tejido_actual <- fread(archivos_csv[i])
  setnames(tejido_actual, old = 3, new = "score_actual")
  acumulador <- merge(acumulador, tejido_actual, by = c("prot1", "prot2"), all.x = TRUE)
  
  # Si en algún tejido esa interacción no viniera (NA), le ponemos un score de 0
  acumulador[is.na(score_actual), score_actual := 0]
  
  # Acumulamos las sumas matemáticas
  acumulador[, suma_scores := suma_scores + score_actual]   # sumamos el score de cada interaccion
  acumulador[, suma_cuadrados := suma_cuadrados + (score_actual^2)]
  
  # Eliminamos la columna temporal 
  acumulador[, score_actual := NULL]
  
  rm(tejido_actual); gc()
}

# Calcular Media, Desviación Estándar y COEFICIENTE DE VARIACIÓN COMUNES (CV)
cat("Calculando Coeficiente de Variación Global...\n")
acumulador[, media_global := suma_scores / N]
acumulador[, desv_global := sqrt((suma_cuadrados - (suma_scores^2 / N)) / (N - 1))]

acumulador[, cv_global := desv_global / (media_global + 0.0001)] # añadimos valor a la media (denominador) para impedir dividir entre 0
acumulador[desv_global == 0, cv_global := 0] # y si la desviación es 0  el CV será 0.

metricas_globales <- acumulador[, .(prot1, prot2, media_global, cv_global)]
rm(acumulador); gc()

saveRDS(metricas_globales, file = "./Data//Tissue_Expression_Atlas/expression_CV.rds")

# GRAFICO DE DENSIDAD HEXAGONL --> media global vs CV

ggplot(expression_tissues_CV, aes(x = media_global, y = cv_global)) +
  geom_hex(bins = 70) +
  scale_fill_viridis_c(option = "plasma", name = "Conteo") + # Paleta científica de alta resolución
  geom_hline(yintercept = 0.4, color = "white", linetype = "dashed", size = 1) + # Umbral de corte de especificidad
  labs(
    title = "Relación entre la Media Global y la Especificidad (CV)",
    x = "Media Global del Score",
    y = "Coeficiente de Variación (CV)"
  ) +
  theme_minimal()

# ---- Calcular media del score de cada tejido ----

carpeta_plots_tejidos <-"./Data/Tissue_Expression_Atlas/Plots_Tejidos/"
dir.create("./Data/Tissue_Expression_Atlas/Plots_Tejidos/", showWarnings = FALSE)

carpeta_data_tejidos_CV <-"./Data/Tissue_Expression_Atlas/Data_CV_Tejidos/"
dir.create("./Data/Tissue_Expression_Atlas/Data_CV_Tejidos/", showWarnings = FALSE)

for(i in 1:length(archivos_csv)) {
  
  cat("Generando análisis y gráfico para el tejido:", nombres_tejidos[i], "\n")
  
  # A. Leer el tejido actual
  tejido_actual <- fread(archivos_csv[i])
  setnames(tejido_actual, old = 3, new = "score_tejido")
  # B. Cruzar con el mapa de CV global que ya tienes en tu entorno (expression_tissues_CV)
  datos_combinados <- merge(tejido_actual, expression_tissues_CV[, .(prot1, prot2, cv_global)], 
                            by = c("prot1", "prot2"), all.x = TRUE)
  datos_combinados[is.na(cv_global), cv_global := 0] # por si acaso
  
  ruta_rds <- paste0(carpeta_data_tejidos_CV, "tabla_CV_", nombres_tejidos[i], ".rds")
  saveRDS(datos_combinados, file = ruta_rds)
  
  # C. Crear el gráfico hexagonal para este tejido concreto
  p <- ggplot(datos_combinados, aes(x = score_tejido, y = cv_global)) +
    geom_hex(bins = 70) +
    scale_fill_gradientn(
      colors = c( "#b4d4e7", "#3b82b6", "#6bb970", "#e65539", "#6b003c"),
      name = "Conteo"
    ) +
    geom_hline(yintercept = 0.4, color = "black", linetype = "dashed", size = 0.8) +
    geom_vline(xintercept = 0.8, color = "black", linetype = "dashed", size = 0.8) +
    labs(
      title = paste("Especificidad (CV Global) vs Score en:", toupper(nombres_tejidos[i])),
      subtitle = "Las interacciones interesantes están en el cuadrante superior derecho",
      x = paste("Score de Asociación en", nombres_tejidos[i]),
      y = "Coeficiente de Variación Global (Especificidad)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 11)
    )
  
  # D. Guardar el gráfico 
  nombre_archivo_plot <- paste0(carpeta_plots_tejidos, nombres_tejidos[i], "_plot_CV.pdf")
  ggsave(nombre_archivo_plot, plot = p, width = 8, height = 6, dpi = 300)
  
  # E. LIMPIEZA RADICAL DE RAM: Borramos todo lo usado en esta vuelta
  rm(tejido_actual, datos_grafico, p); gc()
}

cat("¡Proceso terminado! Tienes los 11 gráficos guardados en la carpeta")

df_conversion <- unique(df_conversion)
duplicados <- df_conversion %>%
  filter(duplicated(to_id) | duplicated(to_id, fromLast = TRUE)) %>%
  arrange(to_id)
head(duplicados, 20)

# ---- FILTRAMOS POR SCORE > 0.8 y CV > 0.4 ----

carpeta_data_tejidos_cv_filtrado <- "./Data/Tissue_Expression_Atlas/Data_CV_Tejidos_filtrado/"
dir.create("./Data/Tissue_Expression_Atlas/Data_CV_Tejidos_filtrado/", showWarnings = FALSE)

archivos_rds <- list.files(path = carpeta_data_tejidos_CV, pattern = "\\.rds$", full.names = TRUE)

for(i in 1:length(archivos_rds)) {
  
  cat("Filtrando interacciones específicas para:", nombres_tejidos[i], "\n")
  
  tabla_actual <- readRDS(archivos_rds[i])
  
  data_filtrada <- tabla_actual %>%
    filter(
      prot1 != prot2,        # Quitar auto-bucles
      score_tejido > 0.8,    
      cv_global > 0.4       
    )
  # añadimos ID ENSG a las proteinas
  data_filtrada <- merge(data_filtrada, df_conversion, 
                         by.x = "prot1", by.y = "to_id", 
                         all.x = TRUE)
  setnames(data_filtrada, "from_id", "ENSG1")
  
  data_filtrada <- merge(data_filtrada, df_conversion, 
                        by.x = "prot2", by.y = "to_id", 
                        all.x = TRUE)
  setnames(data_filtrada, "from_id", "ENSG2")
  
  # reorganizar
  data_filtrada <- data_filtrada[, .(prot1, ENSG1, prot2, ENSG2, score_tejido, cv_global)]
  
  ruta_salida <- paste0(carpeta_data_tejidos_cv_filtrado, nombres_tejidos[i], "_CV_filtrado.rds")
  saveRDS(data_filtrada, ruta_salida)

  rm(tabla_actual, data_filtrada); gc()
}

cat("¡Proceso completado! Todas tus redes filtradas están en './Redes_Filtradas_Finales/'\n")

# ver tamaño de los archivos
archivos_filtrados <- list.files(path = carpeta_data_tejidos_cv_filtrado, pattern = "\\.rds$", full.names = TRUE)
resumen_tamanos <- data.frame(
  Tejido = character(), 
  Num_Interacciones = numeric(), 
  stringsAsFactors = FALSE
)

for(archivo in archivos_filtrados) {
  nombre_tejido <- basename(archivo) %>% 
    str_replace("red_final_ENSG_", "") %>% 
    str_replace("\\.rds$", "")
  
  tabla_actual <- readRDS(archivo)

  num_filas <- nrow(tabla_actual) 
  
  resumen_tamanos <- rbind(resumen_tamanos, data.frame(Tejido = nombre_tejido, Num_Interacciones = num_filas))
}
head(readRDS(archivos_filtrados[1]))
# resultado ordenado de mayor a menor
resumen_tamanos %>% 
  arrange(desc(Num_Interacciones)) %>% 
  print()

