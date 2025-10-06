# Leer archivo liquido
original <- read_excel("/data/original.xlsx", skip = 3)

# Eliminar columnas 94 a 97
original <- original[, -c(94:97)]

# Extraer fila con LODs
maxlod <- original[74, ]

maxlod[ , -1] <- lapply(maxlod[ , -1], function(x) ifelse(x < 0, 0, x))

# Crear base de datos sin filas informativas y sin fila de LOD
original2 <- original[-c(1:3, 64:76), ]

# Convertir columnas de biomarcadores a numéricas de forma segura
for (j in 2:93) {
  original2[[j]] <- as.numeric(original2[[j]])
}

# Sustituir valores por debajo del LOD por el LOD
for (i in names(original2)) {
  if (i != 'Assay') {  # Asumiendo que 'Assay' es una columna no biomarcador
    
    # Solo modificar si la columna existe en maxlod
    if (!is.null(maxlod[[i]])) {
      lod_value <- as.numeric(maxlod[[i]])

      # Reemplaza valores por debajo del LOD
      original2[[i]][original2[[i]] < lod_value] <- lod_value
    }
  }
}

original2[2:93]<-round(original2[2:93],3)

write.table(original2,'/data/original_lod.csv',row.names = F,na='',sep=';')

# Leer archivo suero
original_suero <- read_excel("/data/original_serum.xlsx", skip = 3)

# Eliminar columnas 94 a 97
original_suero <- original_suero[, -c(94:95)]

# Extraer fila con LODs
maxlod_suero <- original_suero[43, ]

maxlod_suero[ , -1] <- lapply(maxlod_suero[ , -1], function(x) ifelse(x < 0, 0, x))

# Crear base de datos sin filas informativas y sin fila de LOD
original_suero2 <- original_suero[-c(1:3, 32:46), ]

# Convertir columnas de biomarcadores a numéricas de forma segura
for (j in 2:93) {
  original_suero2[[j]] <- as.numeric(original_suero2[[j]])
}

# Sustituir valores por debajo del LOD por el LOD
for (i in names(original_suero2)) {
  if (i != 'Assay') {  # Asumiendo que 'Assay' es una columna no biomarcador
    
    # Solo modificar si la columna existe en maxlod
    if (!is.null(maxlod_suero[[i]])) {
      lod_value_suero <- as.numeric(maxlod_suero[[i]])
      
      # Reemplaza valores por debajo del LOD
      original_suero2[[i]][original_suero2[[i]] < lod_value_suero] <- lod_value_suero
      colnames(original_suero2)[colnames(original_suero2) == i] <- paste0(i, '_serum')
      
    }
  }
}

original_suero2[2:93]<-round(original_suero2[2:93],3)

write.table(original_suero2,'/data/original_serum_lod.csv',row.names = F,na='',sep=';')
