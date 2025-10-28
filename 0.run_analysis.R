# Establecemos el working directory
setwd('/data/')

# Cargamos librerias, creamos carpetas y funciones
source('./1.library_functions.R')

# Pasamos los datos del archivo original que tiene valores por debajo del LOD al propio valor de LOD por columna
source('./2.original_to_lod.R') # se crea el archivo en datos: original_lod.csv

# Realizamos la estadística descriptiva de todo el dataframe
source("./3.descriptive_statistics.R")
run_descriptive_statistics("/data/results/tabla_descriptiva.xlsx")

combinacion02<-combinacion
combinacion02$microcristales02<-ifelse(combinacion02$microcristales=='0' | combinacion02$microcristales=='2',0,1)

# Realizamos la estadística agrupara por columna (deposito, microcristales y sinovitis)
source("./4.inferential_statistics.R")
run_inferential_statistics(combinacion, group_col = "deposito", output_path = "/data/results/Deposito articular/1. Estadística Básica/estadistica_deposito.xlsx")
run_inferential_statistics(combinacion, group_col = "microcristales", output_path = "/data/results/Microcristales/1. Estadística Básica/estadistica_microcristales.xlsx")
run_inferential_statistics(combinacion, group_col = "sinovitis", output_path = "/data/results/Sinovitis ecográfica/1. Estadística Básica/estadistica_sinovitis.xlsx")

# Realizar BoxPlot
source("./5.boxplot.R")
generate_boxplots(combinacion, "deposito", output_dir = "/data/results/Deposito articular/1. Estadística Básica/Boxplot por grupos/")
generate_boxplots(combinacion, "microcristales", output_dir = "/data/results/Microcristales/1. Estadística Básica/Boxplot por grupos/")
generate_boxplots(combinacion, "sinovitis", output_dir = "/data/results/Sinovitis ecográfica/1. Estadística Básica/Boxplot por grupos/")

source("scripts/final/5.1.boxplot_combined.R")

# Realizar VolcanoPlot
source("./6.volcanoplot.R")
generate_volcano_plots(combinacion, group_col = "deposito", output_dir = "/data/results/Deposito articular/2. Volcano Plot/")
generate_volcano_plots(combinacion, group_col = "microcristales", output_dir = "/data/results/Microcristales/2. Volcano Plot/")
generate_volcano_plots(combinacion, group_col = "sinovitis", output_dir = "/data/results/Sinovitis ecográfica/2. Volcano Plot/")
generate_volcano_plots(combinacion02, group_col = "microcristales02", output_dir = "/data/results/")

# Realizar PCA y PCA significativos
source("./7.PCA.R")
run_pca_analysis(group_col = "deposito",  output_dir = "/data/results/Deposito articular/3. PCA/")
run_pca_analysis(group_col = "microcristales",  output_dir = "/data/results/Microcristales/3. PCA/")
run_pca_analysis(group_col = "sinovitis",  output_dir = "/data/results/Sinovitis ecográfica/3. PCA/")

source("./8.corr_heat.R")
# df_correlacion <- read.csv("/data/datos/df_correlacion.csv", sep = ";", header = TRUE)
analizar_synovial_serum(df_correlacion_meta)

