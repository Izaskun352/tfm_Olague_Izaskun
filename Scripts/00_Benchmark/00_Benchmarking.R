
# LIBRERIAS

if (!require("pacman")) install.packages("pacman")
pacman::p_load(stringr,purrr,dplyr,igraph, pROC,rvest,arrow,ggplot2, tidyverse)

# ====================================
# BENCHMARKING
# ====================================

# -----------------------------------------
# DESCARGAR DEL DATASET DE OPENTARGETS ====

# --- Archivos gold standards ---

#Clinical percedence evidence
url <- "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_clinical_precedence/"
carpeta_entrada<-"./Data/Benchmark/chmbl"

dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))

#Genomics England-PanelApp evidence
url <- "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_genomics_england/"
carpeta_entrada<-"./Data/Benchmark/genomicEngland"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


#ClinGen evidence
url <- "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_clingen/"
carpeta_entrada<-"./Data/Benchmark/clingen"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


#Reactome evidence
url <- "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_reactome/"
carpeta_entrada<-"./Data/Benchmark/reactome"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


#CRISPR screen evidence
url <- "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_crispr_screen/"
carpeta_entrada<-"./Data/Benchmark/crispr"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


#Cancer Gene Census evidence

url<- "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_cancer_gene_census/"
carpeta_entrada<-"./Data/Benchmark/cosmic"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


#Cancer biomarker evidence

url<- "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_cancer_biomarkers/"
carpeta_entrada<-"./Data/Benchmark/biomarkers"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))


#Gene burden evidence

url<-"http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_gene_burden/"
carpeta_entrada<-"./Data/Benchmark/geneburden"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))

#Gene2Phenotype evidence
url<-"http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_gene2phenotype/"
carpeta_entrada<-"./Data/Benchmark/gene2phenotype"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))

#Orphanet evidence
url<- "http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_orphanet/"
carpeta_entrada<-"./Data/Benchmark/orphanet"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))

#ProjectScore evidence
url<-"http://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/evidence_crispr/"
carpeta_entrada<-"./Data/Benchmark/projectScore"
if (dir.exists(carpeta_entrada)) {unlink(carpeta_entrada, recursive = TRUE, force = TRUE)}
dir.create(carpeta_entrada, showWarnings = FALSE, recursive = TRUE)
archivos <- read_html(url) %>% html_elements("a") %>% html_attr("href") %>% .[grep("parquet$", .)]
lapply(archivos, function(f) download.file(paste0(url, f), file.path(carpeta_entrada, f), mode = "wb"))

# -----------------------------------------

# INPUTS

all_gene_GWAS <- readRDS("./Data/nasertic/input/all.gene.gwas_filter_GP_nonred.rds")
carpeta<-"./Data/nasertic/output/RDS_gwas"
seeds <- as.data.frame(readRDS("./Data/nasertic/input/all.gene.gwas_filter_GP_nonred.rds"))## genes semilla
traitstotal<-unique(seeds$disease)

# OUTPUTS

carpeta_destino<-"./Output/ROC"
dir.create(carpeta_destino, showWarnings = FALSE, recursive = TRUE)
gold_standard<- readRDS(paste0(carpeta_destino,"/gold_standard.rds"))

# -----------------------------------------
# --- 1- Filtramos las bases de datos ----

#Clinical percedence evidence
chmbl_ <- arrow::open_dataset("./Data/Benchmark/chmbl") %>%
  filter( clinicalStage %in% c("PHASE_4", "APPROVAL", "PREAPPROVAL", "PHASE_3")) %>%
  select(diseaseId, targetId)%>%
  collect() 

#Genomics England-PanelApp evidence
england_ <- arrow::open_dataset("./Data/Benchmark/genomicEngland") %>%
  filter( score==1) %>%
  select(diseaseId, targetId)%>%
  collect() 

#ClinGen evidence
clingen_ <- arrow::open_dataset("./Data/Benchmark/clingen") %>%
  filter( score==1) %>%
  select(diseaseId, targetId)%>%
  collect() 

#Reactome evidence (No tengo en cuenta el score porque es siempre 1)
reactome_<- arrow::open_dataset("./Data/Benchmark/reactome")%>%
  select(diseaseId, targetId)%>%
  collect() 

#CRISPR screen evidence
crispr_<-arrow::open_dataset("./Data/Benchmark/crispr")%>%
  filter(score>=0.75) %>%
  select(diseaseId, targetId)%>%
  collect()

#Cancer Gene Census evidence
cosmic_<-arrow::open_dataset("./Data/Benchmark/cosmic")%>%
  filter(score>=0.75)%>%
  select(diseaseId, targetId)%>%
  collect()

#Cancer biomarker evidence
biomarker_<-arrow::open_dataset("./Data/Benchmark/biomarkers")%>%
  select(diseaseId, targetId)%>%
  collect()

#Gene burden evidence
geneburden_<-arrow::open_dataset("./Data/Benchmark/geneburden")%>%
  filter(score>=0.75)%>%
  select(diseaseId, targetId)%>%
  collect()

#Gene2Phenotype evidence
gene2phenotype_<-arrow::open_dataset("./Data/Benchmark/gene2phenotype")%>%
  filter(score==1)%>%
  select(diseaseId, targetId)%>%
  collect()

#Orphanet evidence
orphanet_<- arrow::open_dataset("./Data/Benchmark/orphanet")%>%
  filter(score==1)%>%
  select(diseaseId, targetId)%>%
  collect()

#ProjectScore evidence
projectScore_<- arrow::open_dataset("./Data/Benchmark/projectScore")%>%
  filter(score>=0.75)%>%
  select(diseaseId, targetId)%>%
  collect()

datasets_list <- list(
  chmbl = chmbl_,
  england = england_,
  clingen = clingen_,
  reactome = reactome_,
  crispr = crispr_,
  cosmic = cosmic_,
  biomarker = biomarker_,
  geneburden = geneburden_,
  gene2phenotype = gene2phenotype_,
  orphanet = orphanet_,
  projectScore = projectScore_
)


# -----------------------------------------

# -----------------------------------------
# --- 2- Función para quedarnos con genes asociados al rasgo en las db pero NO son genes semilla ----

obtener_targets <- function(dataset, trait_id) {
  seeds_filtrado <- seeds %>%
    filter(disease == trait_id)   # nos quedamos con los genes semilla solo de ese trait
  
  dataset %>%    
    filter(diseaseId == trait_id) %>%  # vemos informaciuón del dataset de ese trait
    select(targetId) %>%   # seleccionamos genes asociados a ese trait en el dataset
    distinct() %>%
    collect() %>%
    anti_join(seeds_filtrado, by = c("targetId" = "gene")) %>%   # se queda solo con los genes en dataset pero no en seeds (elimina los semilla)
    pull(targetId)
}

# -----------------------------------------

# -----------------------------------------
# --- 3- Gold standard ----

#Tabla vacía 
gold_standard <- data.frame(targetId = character(), 
                            diseaseId = character(), 
                            fuente = character(),
                            stringsAsFactors = FALSE) 
# Bucle para cada trait
for( trait in traitstotal){

  tabla<-as.data.frame(readRDS(paste0(carpeta, "/ZSCO.", trait, ".rds")))  # Tabla con resultados propagación
  tabla<-  tabla %>%
    filter(Selected.cluster==1)
  
  gs_list <- map(datasets_list, ~ obtener_targets(.x, trait))  # map recorre cada elemento de la lista (cada database) y aplica la función--> devuelve una lista con los resultados
  # obtenemos los genes que no son seed de cada db para ese trait
  for(fuente in names(gs_list)){   # recorremos cada database
    genes_fuente <- gs_list[[fuente]]   # los genes goldstandars de esa database
    if(length(genes_fuente) >= 10){   # seleccionamos el trait solo si hay minimo 10 genes goldstandard en esa database
      gold_standard <- rbind(gold_standard, data.frame(
        targetId = genes_fuente,
        diseaseId = trait,
        fuente = fuente,
        stringsAsFactors = FALSE
      ))
    }
  }
  
}

saveRDS(gold_standard, file =  paste0(carpeta_destino,"/gold_standard.rds"))

# -----------------------------------------

# -----------------------------------------
# --- 4- ROC ----

fuentes<-c()
evaluar<-gold_standard%>%   # seleccionamos los traits que vamos a evaluar y la fuente que vamos a utilizar
  select(diseaseId, fuente)%>% 
  distinct()

# Crear tabla vacia
tabla_auc <- data.frame(  trait = character(),  
                          fuente = character(),
                          auc = numeric(),
                          n_gs = numeric(),
                          stringsAsFactors = FALSE
)

for (i in 1:nrow(evaluar)){
  trait <- evaluar$diseaseId[i]
  f <- evaluar$fuente[i]
  
  #saco un vector con los gold standarsd del trait.
  gs <- gold_standard %>% 
    filter(diseaseId == trait, fuente==f) %>% 
    pull(targetId)
  
  #Tabla con el nombre del gen su score en pagerank y 1/0 en función si es gs o no
  
  if (file.exists(paste0(carpeta, "/ZSCO.", trait, ".rds"))) {
    propagacion <- as.data.frame(readRDS(paste0(carpeta, "/ZSCO.", trait, ".rds")))
  } 
  propagacion <- propagacion %>%
    mutate(page.rank = as.numeric(page.rank)) %>%
    select(ENSG, page.rank) %>%
    mutate(gs = ifelse(ENSG %in% gs, 1, 0)) # añadimos columna gs(1 = positivos, 0 = negativos)
  roc1 <- roc(propagacion$gs, propagacion$page.rank, quiet = TRUE)  # calculamos roc 
  
  fila <- data.frame(
    trait = trait,
    fuente=f,
    auc = as.numeric(auc(roc1)),  # ponemos valor del área de la curva ROC
    n_gs = length(gs),  # num de genes goldstandard
    stringsAsFactors = FALSE
  )
  tabla_auc <- rbind(tabla_auc, fila) 
}

# -----------------------------------------

# -----------------------------------------
# 5- BOXPLOT ----

# Eliminamos las datasets con < 10 traits

tabla_auc_fltrada <- tabla_auc %>%
  filter(!(fuente %in% c("clingen", "crispr", "gene2phenotype", "orphanet")))


etiqueta <- tabla_auc_fltrada %>%  # Fuente y número de traits que hay en esa fuente
  group_by(fuente) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(fuente, "\n(n:", n, ")"))

tabla_auc_fltrada <- tabla_auc_fltrada %>%     # añadimos la etiqueta al df anterior
  left_join(etiqueta, by = "fuente")

# Gráfico
ggplot(tabla_auc_fltrada, aes(x = auc, y = label)) +
  geom_boxplot(outlier.shape = NA, fill = "burlywood2", color = "tan4", width = 0.6) + 
  geom_jitter(shape = 21, fill = "peru", color = "tan4",    # añadimos los puntos individuales (jitter mueve los puntos horizontalmente para que no se solapen)
              alpha = 0.6, width = 0.2, size = 2) +
  theme_bw() +
  labs(title = "ROC Pagerank",
       x = "AUC",
       y = " ") +
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0.7, linetype = 5, color = "brown4") + 
  xlim(0, 1)

ggsave(paste0(carpeta_destino,"/RocGwas.pdf"), width = 11, height = 8)
ggsave(paste0(carpeta_destino,"/RocGwas.png"), width = 11, height = 8)
