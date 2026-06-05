

## Abrimos script con las librerias

source("scripts/00_setup.R")

#----------------------------------------

## Abrimos archivos necesarios (los inputs)

all.gene.gwas_filter_GP_nonred <- readRDS("./Data./nasertic/input/all.gene.gwas_filter_GP_nonred.rds")
all.gene.gwas_filter_RARE_nonred <- readRDS("./Data./nasertic/input/all.gene.gwas_filter_RARE_nonred.rds")
All_diseases <- readRDS("./Data./nasertic/input/All_diseases.rds")

Diccionario_TherapeuticAreas <- readRDS("./Data./Diccionarios/Diccionario_TherapeuticAreas.rds")

#------------------------------------------------------------
# VARIACIÓN COMÚN 
#------------------------------------------------------------

  ## Número de rasgos
rasgos_comun <- unique(all.gene.gwas_filter_GP_nonred[,3])
length(rasgos_comun)  # 1844 rasgos con variación común
    #summary(rasgos_comun)

  ## Número de genes
genes_comun <-  unique(all.gene.gwas_filter_GP_nonred[,1])
length(genes_comun)  # 10.727 genes asociados a variación común

  ## Genes por rasgo
genes_por_rasgo_GP <- as.data.frame(table (all.gene.gwas_filter_GP_nonred[,3]))
colnames(genes_por_rasgo_GP) <- c("Rasgo", "Numero_Genes")
genes_por_rasgo_GP <- as.matrix(genes_por_rasgo_GP) # Matriz: cuantos genes hay por rasgo

  ## Distribución de área terapéutica

areas_terapeuticas_GP <- left_join(        # Unimos la columna de área terapéutica a cada rasgo (matriz: genes/rasgo)
  as.data.frame(genes_por_rasgo_GP),  # Como usamos left_join convertimos las matrices en dataframes
  as.data.frame(All_diseases),
  by = c("Rasgo" = "ID")                  # aparecen NA!!!!! hay traits que no tienen área terapéutica
)
areas_terapeuticas_GP <- areas_terapeuticas_GP %>%   # Un mismo rasgo - más de un área terapéutica
  separate_rows(therapeuticAreas, sep = ";")              #separamos áreas terapéuticas por ';' y duplicamos filas 

areas_terapeuticas_GP_ordenado <- areas_terapeuticas_GP %>% 
  add_count(therapeuticAreas) %>%                 # ordena el df en función de la variable therapeutic areas
  arrange(desc(n), therapeuticAreas)               # de mas abundante a menos abundante

    ### Vemos que hay valores NA en las áreas- no todos los rasgos están asociados a un área
traits_sin_area <- is.na(areas_terapeuticas_GP_ordenado$name)  # 107 rasgos no asociados a un área

    ### sustituimos NA por 'Not Annotated'
areas_terapeuticas_GP_ordenado[is.na(areas_terapeuticas_GP_ordenado)] <- "Not annotated"

  ## Asignar nombre del área terapéutica

areas_terapeuticas_GP_ordenado <- left_join(areas_terapeuticas_GP_ordenado, Diccionario_TherapeuticAreas, by = c("therapeuticAreas" = "ID_area"))

length(unique(areas_terapeuticas_GP_ordenado$Nombre_area))   # Hay 26 áreas terapéuticas


#------------------------------------------------------------
# VARIACIÓN RARA 
#------------------------------------------------------------

  ## Número de rasgos
  
rasgos_raros <-  unique(all.gene.gwas_filter_RARE_nonred[,3])
length(rasgos_raros)  # 3258 rasgos asociados a variación rara

  ## Número de genes
genes_raros <- unique(all.gene.gwas_filter_RARE_nonred[,1])
length(genes_raros)  # 7091 genes asociados a variación rara

  ## Número de genes por rasgo

genes_por_rasgo_rare <- as.data.frame(table(all.gene.gwas_filter_RARE_nonred[,3]))
colnames(genes_por_rasgo_rare) <- c("Rasgo", "Numero_Genes")
genes_por_rasgo_rare <- as.matrix(genes_por_rasgo_rare) # Matriz: cuantos genes hay por rasgo

  ## Distribución de áreas terapéuticas

areas_terapeuticas_RARE <- left_join(        # Unimos la columna de área terapéutica a cada rasgo
  as.data.frame(genes_por_rasgo_rare),  # Como usamos left_join convertimos las matrices en dataframes
  as.data.frame(All_diseases),
  by = c("Rasgo" = "ID")
)
areas_terapeuticas_RARE <- areas_terapeuticas_RARE %>%   # Un mismo rasgo - más de un área terapéutica
  separate_rows(therapeuticAreas, sep = ";")              #separamos áreas terapéuticas por ; y duplicamos filas 

areas_terapeuticas_RARE_ordenado <- areas_terapeuticas_RARE %>% 
  add_count(therapeuticAreas) %>%                 # ordena el df en función de la variable therapeutic areas
  arrange(desc(n), therapeuticAreas)  

sum(is.na(areas_terapeuticas_RARE_ordenado$therapeuticAreas))
areas_terapeuticas_RARE_ordenado[is.na(areas_terapeuticas_RARE_ordenado)] <- "Not annotated"

areas_terapeuticas_RARE_ordenado <- left_join(
  areas_terapeuticas_RARE_ordenado, Diccionario_TherapeuticAreas, by = c("therapeuticAreas" = "ID_area"))

#------------------------------------------------------------
# TRAITS CON VARIACIÓN COMÚN Y RARA 
#------------------------------------------------------------

# Rasgos con genes con variación común y rara
rasgos_compartidos <- intersect (rasgos_comun, rasgos_raros)
length(rasgos_compartidos)   # 234 traits asociados a variación rara y común


