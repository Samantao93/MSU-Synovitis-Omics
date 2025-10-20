path_sign='/boxplot/significativas/'

# estadistica <- read_excel("results/estadistica.xlsx")
# estadistica_filter<-filter(estadistica,test<0.001) # Filtramos 0.01 y 0.001.

# files_filter<-paste0(path_sign,estadistica_filter$col,'_boxplot.svg')

files_filter <- list.files(path_sign, pattern = "\\.svg$", full.names = TRUE)
files_filter <- files_filter[!grepl("leucocito", files_filter)]
n <- length(files_filter)
cols <- ceiling(sqrt(n))
rows <- ceiling(n / cols)
tile_str <- paste0(cols, "x", rows)

imgs <- lapply(files_filter, function(f) {
  image_read_svg(f)
})

combined <- image_montage(
  image_join(imgs),
  tile = tile_str,
  geometry = "+12+12",     
  bg = "white"             
)

# Guardar salida (elige uno)
# out_png <- file.path(path_sign, "combinado_auto.png")
image_write(image_convert(combined, "jpeg"),
            path   = file.path(path_sign, "combinado_auto.jpg"), # 01 y 001
            format = "jpeg",
            quality = 100,      
            density = "300x300")

