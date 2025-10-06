# Cargamos los grupos
grupos <- read_excel("/data/metadata.xlsx", range = "A1:I61")
colnames(grupos)<-c('cod_paciente','microcristales','suero','leucocito','edad','sexo','ac_urico','deposito','sinovitis')
grupos$cod_paciente<-substr(grupos$cod_paciente,7,8) %>% gsub("^0","",.) %>% paste0('m',.)

# Cargamos las tablas
df <- read.csv("/data/original_lod.csv", sep = ";", header = TRUE)
names(df)[1]<-"cod_paciente"
df$cod_paciente <- gsub(" ", "", df$cod_paciente)

combinacion<-full_join(grupos,df)
write.table(combinacion,file='/data/combinacion.csv',na='',sep=';',row.names = F)

# Cargamos las tablas de suero (hacerlo cuando los cambios en los metadatos se realice)
df_serum <- read.csv("/data/original_serum_lod.csv", sep = ";", header = TRUE)
names(df_serum)[1]<-"cod_paciente"
df_serum$cod_paciente <- gsub(" |Ser", "", df_serum$cod_paciente)

combinacion_serum<-full_join(grupos,df_serum)

combinacion_serum<-filter(combinacion_serum,!is.na(suero))
write.table(combinacion_serum,file='/data/combinacion_serum.csv',na='',sep=';',row.names = F)

# Para adelatar hago la combinación con las propias muestras
df_correlacion<-left_join(df_serum,df)
write.table(df_correlacion,file='/data/df_correlacion.csv',na='',sep=';',row.names = F)

# Seleccionas solo las columnas de interés de la tabla "combinacion"
meta <- dplyr::select(combinacion,
                      cod_paciente,
                      microcristales,
                      suero,
                      deposito,
                      sinovitis)

# Haces join con df_correlacion (que tiene tus biomarcadores)
df_correlacion_meta <- dplyr::left_join(df_correlacion, meta, by = "cod_paciente")
write.table(df_correlacion_meta,file='/data/df_correlacion_meta.csv',na='',sep=';',row.names = F)


list_qual<-c("microcristales",'suero','sexo','deposito','sinovitis')

run_descriptive_statistics <- function(output_path) {
  summary_stats <- calc_singrupo(combinacion,list_qual)
  writexl::write_xlsx(summary_stats, output_path)
  message("Estadística descriptiva guardada en: ", output_path)
}
