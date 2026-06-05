

source("scripts/00_setup.R")  # Abrimos script con las librerias
source("scripts/03_Funciones_Propias.R")  # Abrimos script con las funciones


#----------------------------------------

## RED DE PROPAGACION 

#--------------------------------------
# INPUTS
#--------------------------------------

Allgene_microGWAS <- as.data.frame(readRDS("./Data/Microbioma/all.gene.gwas_microGWAS.rds"))
nodos_interactoma <- as.data.frame(readRDS("./Data/nasertic/input/all.node.gwas_STRING40_OTAR0924_FILTER.rds"))
string_interactoma <- as.data.frame(readRDS("./Data/nasertic/input/Combined_STRING40_OTAR0924_FILTER.rds"))

edge_string <- string_interactoma  

# read.csv(file = "./Data/Microbioma/microbiome_gwas_master_list.csv")

#--------------------------------------
# OUTPUTS
#--------------------------------------

carpeta_MicroGWAS <- "./Output/MicroGWAS"

# ------------------------------------------------
# FUNCIONES NECESARIAS
# ------------------------------------------------

### Anotamos con padj = 100 los genes GWAS (semilla) de ese trait

anotar_genes_GWAS_trait <- function (lista_nodos,      # El padj de all_genes_GWAS es siempre 100 (no se calcula max padj)
                                     all_genes_GWAS,
                                     target_trait) {
  genes_GWAS <- all_genes_GWAS %>%
    filter(disease %in% target_trait) %>%    # nos quedamos con los genes de la disease indicada
    pull(gene) %>%     # pull() extrae la columna como un vector normal
    unique()           # Nos quedamos con valores únicos por si hay repetidos
  
  lista_nodos <- lista_nodos %>%
    mutate(
      padj = ifelse(ENSG %in% genes_GWAS, 100, 0)
    )
  return(lista_nodos)
}

### Propagación en red

network_propagation <- function (nodos_gwas, # genes gwas
                                 edge_string,  # interactoma
                                 all_nodos) { # si queremos obtener df con todos los nodos
  
  # 1- Propagation / Difusión
  
  ## Creamos red a partir de un df
  net = graph_from_data_frame(d = edge_string,   # d: df con las aristas (interactoma)
                              vertices=nodos_gwas,  # vertices: df con los nodos, no dirigido
                              directed=F) # no dirigido
  
  ## Añadimos atributo weight (score de la interacción en el interactoma)
  E(net)$weight = as.numeric(as.character(edge_string[,"combined_score"])) 
  
  ## Limpiamos la red
  net_clean <- igraph::simplify(net,   
                                remove.loops = T,   # quita conexiones consigo mismo
                                remove.multiple = T ,   # quita si hay varias conexiones exactamente con el mismo par
                                edge.attr.comb = c(weight="max","ignore"))  # nos quedamos con el valor máximo de weight al eliminar conexiones múltiples
  
  ## Calculamos el pagerank personalizado (PPR)
  page_rank <- page_rank(net_clean, 
                         personalized=as.numeric(nodos_gwas[,"padj"]), 
                         weights=E(net_clean)$weight)
  
  nodos_gwas <- cbind(nodos_gwas,page_rank$vector)  # añadimos columna pagerank (añadimos vector con el score para cada nodo)
  colnames(nodos_gwas)[ncol(nodos_gwas)]="page.rank"  # nombramos esta nueva columna como 'page.rank'
  
  degree <- igraph::degree(net_clean)  # calculamos el núm exacto de conexiones directas que tiene cada nodo
  
  # 2- Filtramos los nodos para quedarnos con top25 con mayor score
  
  pagerank_numeric <- as.numeric (nodos_gwas[, "page.rank"])
  threshold <- quantile (pagerank_numeric, probs = 0.75, na.rm = TRUE)
  nodos_filtrado <- nodos_gwas[pagerank_numeric > threshold, ]
  colnames(nodos_filtrado)[1] <- "ENSP"
  
  ## Filtramos tambien las aristas para quedarnos solo con interacciones entre nodos top25
  
  top_nodos <- as.character(nodos_filtrado[, "ENSP"])
  edge_v1 <- as.character(edge_string[, 1])
  edge_v2 <- as.character(edge_string[, 2])
  
  edge_filtrado <- edge_string[edge_v1 %in% top_nodos & edge_v2 %in% top_nodos, ]
  
  ## Quitamos nodos huérfanos
  
  nodos_activos <- unique(c(as.character(edge_filtrado[, 1]), as.character(edge_filtrado[, 2])))
  nodos_filtrado <- nodos_filtrado[as.character(nodos_filtrado[, "ENSP"]) %in% nodos_activos, ]
  
  #3- Creamos la nueva red con los nodos filtrados
  
  ## Cambiamos las aristas/edges
  edge_string <- data.frame(
    source = edge_filtrado[, 1],
    target = edge_filtrado[, 2],
    weight = as.numeric(as.character(edge_filtrado[, "combined_score"]))  #graph_from_data_frame convierte automáticamente la 3ra columna en atributos de la arista
  )
  
  net <- graph_from_data_frame(d = edge_string, 
                               vertices = nodos_filtrado, 
                               directed = FALSE)
  
  ## Limpiamos la nueva red
  net_clean <- igraph::simplify(net, 
                                remove.loops = TRUE,
                                remove.multiple = TRUE,
                                edge.attr.comb = c(weight = "max", "ignore"))
  
  
  # 4- Hacemos clustering (random.walktrap)
  
  cwt <- cluster_walktrap(net_clean,    # hacemos cluistering en la red con los nodos top25
                          weights = E(net_clean)$weight, 
                          steps = 6)   # longitud de los caminos aleatorios (simulará recorridos de 6 saltos consecutivos)
  ## actualizamos df para añadir degree, cluster.walktrap y modularity.walktrap
  
  nodos_filtrado$degree <- igraph::degree(net_clean)
  nodos_filtrado$cluster.walktrap <- as.character(cwt$membership)
  nodos_filtrado$modularity.walktrap <- modularity(cwt)
  ## calculamos frecuencia de los clusters
  cluster_counts <- table(nodos_filtrado$cluster.walktrap)  # name: cluster ; frecuency: frecuencia
  
  ## seleccionamos los nombres de los clusters que tengan > 300 nodos
  
  clusters_to_process <- names(cluster_counts)[cluster_counts >= 300]
  
  ## reclustering hasta que ningún cluster tenga >= 300 nodos ---> bucle
  
  contador_seguridad <- 0
  max_intentos <- 20 # Solo le dejamos intentar romperlo 20 veces
  
  while (length(clusters_to_process) > 0 && contador_seguridad < max_intentos) {
    
    contador_seguridad <- contador_seguridad + 1
    
    current_cluster <- clusters_to_process[1] # Sacamos el primer cluster de la lista
    clusters_to_process <- clusters_to_process[-1] # Lo quitamos de la cola
    
    message("   -> Rompiendo mega-cluster: ", current_cluster, 
            " (Quedan ", length(clusters_to_process), " grandes por procesar...)")
    message("   -> [Intento ", contador_seguridad, "/", max_intentos, "] Rompiendo mega-cluster: ", current_cluster)
    
    nodos_ensp <- nodos_filtrado$ENSP[nodos_filtrado$cluster.walktrap == current_cluster]  # nodos de este cluster
    
    ## Extraemos el subgrafo de solo los nodos de este cluster
    sub_net <- induced_subgraph(net_clean, vids = which(V(net_clean)$name %in% nodos_ensp))  #induced_subgraph va a la red original (net_clean) y recorta solo los nodos de nodos_ensp
    
    # Re-clustering del subgrafo
    cwt_sub <- cluster_walktrap(sub_net, 
                                weights = E(sub_net)$weight, 
                                steps = 6)
    
    ## Renombramos los ID de los clusters (1;2, etc)
    nuevos_sub_clusters <- paste(current_cluster, cwt_sub$membership, sep = ";")
    
    ## buscamos donde (en qué filas) se encuentran los nodos en el df original
    match_idx <- match(V(sub_net)$name, nodos_filtrado$ENSP) 
    
    ## Sobrescribimos con los nuevos datos
    nodos_filtrado$cluster.walktrap[match_idx] <- nuevos_sub_clusters
    mod_value <- modularity(cwt_sub)
    nodos_filtrado$modularity.walktrap[match_idx] <- rep(mod_value, length(match_idx))
    
    ## vemos si sigue habiendo clusters con >= 300 nodos
    
    nuevos_conteos <- table(nuevos_sub_clusters)
    nuevos_grandes <- names(nuevos_conteos)[nuevos_conteos >= 300]  # si el cluster que acabamos de procesar ha creado otro con >=300 nodos
    if (length(nuevos_grandes) > 0) {
      if (contador_seguridad >= max_intentos) {
        message("   -> ⚠️ Límite alcanzado. El cluster ", current_cluster, " es demasiado denso. Se queda como está.")
      } else {
        clusters_to_process <- c(clusters_to_process, nuevos_grandes) # si hay otro con >=300 lo añadimos a los de antes para que se procese
      }
    }
  }
  
  # 5- Devolvemos lo que nos piden en funcion de all_nodos
  if(all_nodos == FALSE){  
    
    return(nodos_filtrado)	# solo le devolvemos los nodos filtrados (el top25 de PPR + clustering al que pertenecen)
    
  }else{
    # Extraemos de la red filtrada todo lo que queremos pegar: ID, Cluster y Grado
    temp <- nodos_filtrado[, c("ENSP", "cluster.walktrap", "degree")] 
    
    # Cruzamos los datos
    nodos_gwas <- merge(nodos_gwas, 
                        temp,
                        by.x = "ENSG",  # la columna de nodos_gwas
                        by.y = "ENSP",  # Nuestra columna renombrada en nodos_filtrados
                        all.x = TRUE)   # fuerza a conservar todas las filas de la matriz original (nodos_gwas)
    
    # los nodos que no estan en el top25 tendrán "NA" en cluster y degree
    
    return(nodos_gwas)   # no devolvemos una matriz --> podemos poner as-matrix()
  }
}

### Estudio de la significancia de clusters (test KS y Fisher)

evaluar_significancia_clusters <- function(df_nodes,      # podriamos añadir tambien si quisieramos cambiar los test estadísticos pero tendríamos que cambiar la función
                                           num_minimo_genes_semilla = 1,
                                           num_minimo_nodes_cluster = 10,
                                           col_cluster = "cluster.walktrap",   # si cambiamos hay que poner .data[[]] + col_cluster, etc
                                           col_padj = "padj",
                                           col_score = "page.rank") {
  df_nodes <- as.data.frame(df_nodes)
  # 1. Selected cluster
  clusters_validos <- df_nodes %>%
    
    group_by(cluster.walktrap) %>%
    summarise(
      all_nodes = n(),  # contamos el total de nodes del cluster
      gwas_nodes = sum(as.numeric(padj) != 0, na.rm = TRUE),   # contamos cuantos nodes son genes semilla
      .groups = "drop"
    ) %>%
    
    filter(all_nodes >= num_minimo_nodes_cluster)  %>%   # que tenga al menos 10 nodes (miembros) el cluster
    filter(gwas_nodes >= num_minimo_genes_semilla) %>%  # que tenga al menos 1 gen semilla
    rename(clust = cluster.walktrap)  # cambiamos nombre para que no se lie con el nombre de la variable del df grande
  
  df_nodes <- df_nodes %>% # Si el cluster del gen está en la lista de los válidos, le ponemos un 1. Si no, un 0.
    mutate( 
      Selected.cluster = ifelse(.data$cluster.walktrap %in% clusters_validos$clust, 1, 0)
    )
  
  # 2. Calculamos significancia de los clusters
  
  
  numero_clusters_validos <- nrow(clusters_validos)
  
  
  # CASO A: 0 Clusters pasaron el filtro previo
  if (numero_clusters_validos == 0) {
    
    message("Ningún cluster superó el filtro inicial. Se devuelven valores a 0.")
    df_nodes <- df_nodes %>%
      mutate(
        Selected.cluster = NA, Selected.fisher = NA, Selected.KS = NA,
        padj.fisher = NA, padj.KS = NA
      )
    return(df_nodes)
    
  } else {    # AL MENOS UN CLUSTER A EVALUAR
    
    message(paste(nrow(clusters_validos), "cluster(s) a evaluar."))
    
    ## Extraemos vectores (evitar warnings 'Unknown column')
    
    PR_vec <- as.numeric(df_nodes[["page.rank"]])
    Clust_vec <- as.character(df_nodes[["cluster.walktrap"]])
    
    ## Variables necesarias
    x_ref_raw <- log10(PR_vec[!is.na(Clust_vec)]) # para calcular significancia KS: log10 del pageRank de todos los nodos
    x_ref_KS <- x_ref_raw[is.finite(x_ref_raw)]  
    
    
    total_gwas <- sum(as.numeric(df_nodes[["padj"]]) != 0, na.rm = TRUE)
    total_nodes <- nrow(df_nodes)  # CUIDADO!!! si hay duplicados no los tiene en cuenta
    
    clusters_evaluados <- clusters_validos %>%
      
      rowwise() %>% # para que se haga fila por fila
      
      mutate(
        # Análisis KS
        
        pval_KS = {
          y_tmp <- PR_vec[Clust_vec == clust & !is.na(Clust_vec)]
          y_log <- log10(y_tmp)
          y_clean <- y_log[is.finite(y_log)] # Limpiamos los infinitos y NAs
          
          # Si quedan datos, hacemos el test en silencio. Si no, devolvemos NA.
          if (length(y_clean) > 0) {
            suppressWarnings(ks.test(x_ref_KS, y_clean, alternative = "greater")$p.value)
          } else {
            NA_real_
          }
        },
        
        # Análisis FIsher
        
        ## Variables para construir la matriz 2x2
        
        gwas_in_cluster = gwas_nodes,  # num de genes GWAS de ESE cluster
        non_gwas_cluster = all_nodes - gwas_nodes,   # genes no GWAS de ese cluster
        gwas_out = total_gwas - gwas_in_cluster,  # genes GWAS FUERA del cluster
        non_gwas_out = (total_nodes - total_gwas) - non_gwas_cluster,  # genes no GWAS FUERA del cluster
        
        matriz_fisher = list(matrix(c(gwas_in_cluster, non_gwas_cluster, gwas_out, non_gwas_out), nrow = 2)),
        pval_fisher = fisher.test(matriz_fisher, alternative = "greater")$p.value
      ) %>%
      
      ungroup() %>%
      
      mutate(
        
        # Ajuste BH
        
        padj.KS = p.adjust(pval_KS, method = "BH"),
        padj.fisher = p.adjust(pval_fisher, method = "BH")
      )
    
    # Pegamos los p-valores calculados
    diccionario <- clusters_evaluados %>% 
      select (cluster.walktrap = clust, padj.KS, padj.fisher)
    
    df_nodes <- df_nodes %>%
      left_join(diccionario, by = "cluster.walktrap")
    
    # CASO B: SOLO TENEMOS UN CLUSTER PARA EVALUAR
    
    if (numero_clusters_validos == 1) {
      
      df_nodes <- df_nodes %>% mutate(
        Selected.cluster = ifelse(!is.na(padj.fisher), "1", "0"),
        Selected.fisher = ifelse(!is.na(padj.fisher), "check", "0"), # ponemos "check" a lo que no sea NA --> pq la significancia estadística es baja al tener pocos valores
        Selected.KS = ifelse(!is.na(padj.KS), "check", "0")
      )
      
    } else {
      
      # CASO C: TENEMOS VARIOS CLUSTERS A EVALUAR
      
      df_nodes <- df_nodes %>% mutate(
        Selected.cluster = ifelse(!is.na(padj.fisher), "1", "0"),  # vemos si el cluster ha llegado al final
        Selected.fisher = ifelse(!is.na(padj.fisher) & padj.fisher <= 0.05, "1", "0"),
        Selected.KS = ifelse(!is.na(padj.KS) & padj.KS <= 0.05, "1", "0")
      )
    }
  }
  
  return(df_nodes)
}

# ------------------------------------------------
# PROCEDIMIENTO
# ------------------------------------------------

# 1- Lo primero es quedarnos solo con los traits que tengan una frecuencia mayor que 1 --> en este caso son todos

table(Allgene_microGWAS[,3])   # todos  

lista_traits_microbioma <- unique(Allgene_microGWAS[,3])

dir.create("./Output/MicroGWAS", recursive = TRUE)  # creamos carpeta
carpeta_destino <- "./Output/MicroGWAS"

# 2- Hacemos iteración para que se haga la propagación por separado para cada trait --> 'elemento para loopear'

resultados_pipeline <- list() # creamos lista vacía para ir guardando los resultados de cada trait

for (trait in lista_traits_microbioma) {
  
  message("Analizando trait: ", trait)
  
  nodos_gwas <- anotar_genes_GWAS_trait(lista_nodos = nodos_interactoma,  ## A partir del df con todos los nodos del interactoma, anotamos los genes que se encuentran aquí y son genes semilla de ese trait en concreto
                          all_genes_GWAS = Allgene_microGWAS,
                          target_trait = trait)
  
  net <- network_propagation(nodos_gwas = nodos_gwas,   ## Corremos la función de propagación para cada trait
                      edge_string = string_interactoma,
                      all_nodos = TRUE)
  
  net_final <- evaluar_significancia_clusters(df_nodes = net)  # Añadir las columnas "Selected.cluster","Selected.fisher","Selected.KS","padj.fisher","padj.KS"
  
  net_final <- net_final %>% ## Se añade la columna con el trait (Siempre el mismo)  --> nombre columna: "Trait"
    mutate(Trait = trait)
  
  #Guardamos cada resultado de cada trait en la carpeta destino
  
  nombre_limpio <- gsub(" ","_", trait)  # limpiamos el nombre (quitamos espacios por _)
  
  nombre_archivo <- paste0 (nombre_limpio, ".rds")  # creamos el nombre del archivo (nombre del trait)
  
  ruta_completa <- file.path(carpeta_destino, nombre_archivo)  # unimos carpeta y archivo
  
  saveRDS(net_final, file = ruta_completa)  # guardamos df en la carpeta
  
  message(trait," guardado con éxito en: ", ruta_completa)
  
  
}


#------------------------------------------------------------------------------
# COMPARACION PARA VER SI LOS RESULTADOS SON SIMILARES
#------------------------------------------------------------------------------

library(mclust)

carpeta_codigo_optimizado <- "./Output/MicroGWAS"

lista_resultados <- list()

carpeta_veridico <- "./Output/RDS_astro"
archivos_codigo_veridico <- list.files(carpeta_veridico, pattern = "\\.rds$")

for (archivo in archivos_codigo_veridico) {
  
  ruta_veridico <- file.path(carpeta_veridico, archivo)
  archivo_nuevo <- gsub("ZSCO.", "", archivo)
  ruta_optimizado <- file.path(carpeta_codigo_optimizado, archivo_nuevo)
  
  df_veridico <- readRDS(ruta_veridico)
  df_veridico <- as.data.frame(df_veridico)
  df_optimizado <- readRDS(ruta_optimizado)
  nombre_trait <- gsub("\\.rds$", "", archivo)
  
  # Asegurarnos de que cruzamos por el identificador del gen (ENSP o ENSG según lo llames)
  col_id <- "ENSG"
  
  # Cruzamos las tablas manteniendo absolutamente TODOS los genes
  comparacion <- inner_join(
    df_veridico %>% select(all_of(col_id), PR_veridico = page.rank, clust_veridico = cluster.walktrap),
    df_optimizado %>% select(all_of(col_id), PR_optimizado = page.rank, clust_optimizado = cluster.walktrap),
    by = col_id
  )
  
  # --- TEST 1: PAGERANK (Spearman) ---
  # Comparamos el ranking de todos los genes de la red
  correlacion_PR <- cor(as.numeric(comparacion$PR_veridico), 
                        as.numeric(comparacion$PR_optimizado), 
                        method = "spearman", use = "complete.obs")
  
  # --- TEST 2: CLUSTERING (ARI) ---
  # Como los genes fuera del top25 pueden tener NA en el cluster, 
  # los cambiamos temporalmente por "Sin_Cluster" para que el ARI los tenga en cuenta
  comparacion <- comparacion %>%
    mutate(
      clust_veridico = ifelse(is.na(clust_veridico), "Sin_Cluster", as.character(clust_veridico)),
      clust_optimizado = ifelse(is.na(clust_optimizado), "Sin_Cluster", as.character(clust_optimizado))
    )
  
  ari_score <- adjustedRandIndex(comparacion$clust_veridico, comparacion$clust_optimizado)
  
  
  # INDICE DE JACCARD   ---> 25%
  
  umbral_veridico <- quantile(as.numeric(comparacion$PR_veridico), probs = 0.75, na.rm = TRUE)
  umbral_optimizado <- quantile(as.numeric(comparacion$PR_optimizado), probs = 0.75, na.rm = TRUE)
  
  # Sacamos las listas de los genes VIP
  genes_top_veridico <- comparacion %>% filter(as.numeric(PR_veridico) > umbral_veridico) %>% pull(col_id)
  genes_top_optimizado <- comparacion %>% filter(as.numeric(PR_optimizado) > umbral_optimizado) %>% pull(col_id)
  
  # Calculamos intersección y unión
  interseccion <- length(intersect(genes_top_veridico, genes_top_optimizado))
  union_genes <- length(union(genes_top_veridico, genes_top_optimizado))
  
  # Porcentaje de Jaccard
  jaccard_elite_pct <- (interseccion / union_genes) * 100
  
  #INDICE DE JACCARD TOP %500
  
  genes_top500_veridico <- comparacion %>% 
    arrange(desc(as.numeric(PR_veridico))) %>% 
    slice_head(n = 30) %>% 
    pull(col_id)
  
  genes_top500_optimizado <- comparacion %>% 
    arrange(desc(as.numeric(PR_optimizado))) %>% 
    slice_head(n = 30) %>% 
    pull(col_id)
  
  interseccion_top500 <- length(intersect(genes_top500_veridico, genes_top500_optimizado))
  union_top500 <- length(union(genes_top500_veridico, genes_top500_optimizado))
  
  jaccard_top500_pct <- (interseccion_top500 / union_top500) * 100
  
  
  # Guardamos los resultados
  resultado_trait <- data.frame(
    Trait = nombre_trait,
    Genes_Totales = nrow(comparacion),
    Correlacion_Global = round(correlacion_PR, 4),
    Jaccard_Top = round(jaccard_elite_pct, 2),
    Jaccard_Top500_Pct = round(jaccard_top500_pct, 2),
    Similitud_Clusters_ARI = round(ari_score, 4)
  )
  
  lista_resultados[[nombre_trait]] <- resultado_trait
}
tabla_validacion_final <- bind_rows(lista_resultados)
print(tabla_validacion_final)









