# tfm_scripts
Este repositorio contiene el código utilizado para el análisis de datos del trabajo de fin de máster.

## 00_Setup.R: 
Principales librerias utilizadas.

## 00_Benchmark 
*00_Benchmarking.R: Benchmark del método de propagación en red. Coge los datos publicados de las bases de datos de genes asociados a enfermedad. Utilizando los resultados de propagación, calcula las curvas ROC y represent gráficamente los resultados de AUC.

## 01:NetworkPropagation: 
* 01_MicroGWAS_NetworkPropagation.R: aplica algoritmo propagación en red a partir del interactoma y los datos de genes #asociados a dfiferentes traits de OTAR. Luego, realiza clustering de genes y calcula los módulos significativos mediante Kolmogorov-Smirnov. 

* Mes1_Datos_GWAS.R: análisis preliminar de los datos de rasgos con variación común y rasgos con variación rara de OTAR. 
 Mes1_Graficos_GWAS: gráficos iniciales sobre los datos de rasgos con variación común y rasgos con variaicón rara de OTAR.

## 02_Manhattan_Jaccard.R
* A_Distancia_Traits_MicroGWAS_GWAS.R: calcula la distancia de Manhattan entre los rasgos asociados a variación común y los rasgos MicroGWAS. Representación gráfica: dendrograma circular y de árbol con las distancias (hclust) + heatmap + red donde cada nodo es un trait.

## 03_Traits_MicroGWAS: 
*Figuras_Traits_MicrogWAS: Representación gráfica datos de los términos EFO de 'microbiome measurement' de GWAS Catalog.

## 04_Clusters_MicroGWAS_GWAS: 
*00_Funciones_Calculo_Intersección: Funciones para calcular el índice de Jaccard entre dos módulos genicos.
*01_Extraccion_Clusters: Funciones + procedimiento de exraer los clusters a partir de los datos de propagación y clustering. Extrae clusters de rasgos de un área terapéutica específica.
*02_Calcular_Intersecciones: Procedimiento para sacar lista con los datos de la intersección de los módulos.
*03_Plot_Red_Pleiotropia.R: Calcula la red con todos ls clusters de rasgos con variación común asociada y los rasgos MicroGWAS. Representa gráficamente las diferentes comunidades formadas.
*05_Plot_Pares_Concretos_Red_Pleiotropia.R: representación gráfica de nodos concretos de la red anterior.
*Expresion_Tisular_Fisher_MicroGWAS: evalúa la expresión de módulos génicos de rasgos MicroGWAs en tejidos específicos utilizando datos de interacción proteína-proteína (PPI), y visualiza los resultados finales mediante heatmaps.
*GOBP_Megaclusters_Micro_GWAS.R: calcula enriquecimiento GOBP de las comunidades de la red con los rasgos de variación común y rasgos MicroGWAS.
*Plot_Components_Red_Pleiotropia.R: Representar una comunidad específica de la red con diferentes filtros.

## 05_Analisis_Clusters: 
*01_Enrich_GOBP: 0_Funciones_GOBP_ funciones para calcular GOBP de una intersección de módulos. A-B: calcula OBP de intersecciones específicas. C: representación gráfca en Heatmap.
*02_Enrich_InterPro: 00: descarga y preparación datos InterPro. 01_ Funcioens para calular enriquecimient Pfam. 2-3: calculo enricquemiento en intersecciones específicas.
*03_Expresión_Tejidos: Cálculo expresion tejidos con datos de Human Proteina ATlas (HPA) (nivel de expresión de cada gen en cada tejido) en los diferentes clusters.
*04_Redes_Genes_Clusters: Plot de las intersecciones entre módulos MicroGWAS y módulos de rasgos neuroinmunes. Colorea los genes semilla y según base de datos de la intersección.
*05_Esencialidad_CRISPR: Cálculo esencialidad genes mediante datos de CRIPSR knockout de DepMap y CRISPRBrain. Calcula esencialidad por intersección del módulo. Representación gráfica en boxplots. 

## 06_Comunidad_cilios: 
*estudio_posicion_variantes.R: estudio de las variantes raras asociadas a los genes presentes en los módulos de la comunidad asociada a cilios. 

## 07_Red_Capas_Finales: 
*01_Funciones_Redes.R: funciones para calcular red, añadir atributos y representarla gráficamente.
*02_TissueExpressionAtlas.R: prepara datos de expresión tisular de interacción proteína-proteína (IPP) para añadir a la red. Agrupa los tejidos en 4 grupos y calcula CV.
*03_CRISPR_KO.R: prepara datod de GRISPR KO para añadir a la red. Agrupa los tejidos en grupos más grandes.
*0_Recursos_Redes.R: prepara datos para la red. Color de las aristas, información de las aristas y datos del dataset de expresión IPP,
A-D: calculo red y representación gráfica de los genes compartidos para cada comunidad seleccionada.
Skin_Haptoglobi: cálculo red y rperesentación gráfica de la comunidad asociada a cilios.

##08_Otros_Analisis: 
*DEG_OTAR.R: analizar qué genes de un módulo están diferencialmente expresados en enfermedades neuropsiquiátricas mediante el dataset de OTAR.
*Estudio_Variantes_codificantes.R: filtrar datset de variantes raras de OTAR para quedarnos solo con las codificantes.
*Expression_atlas_CV.R: analiar y filtrar el dataset de expresion de interaccion-proteina-proteina.
*PyMol.R: código utilizado en PyMol para hacer el estudio y representación de DNAH11.
*Renombrar_Clusters.R: funciones para renombrar los clusters obtenidos en el proceso de clustering tras la propagación.
*Tissue_Expression_Atlas_Red_pleiotropia.R: evalúa la expresión de módulos génicos de las comunidades de la red de pleiotropia en tejidos específicos utilizando datos de interacción proteína-proteína (PPI), y visualiza los resultados finales mediante heatmaps.
