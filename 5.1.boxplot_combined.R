for (i in c("Sinovitis ecográfica","Microcristales")){
  path_sign=paste0('/datos/sam_user/reuma/results/',i,'/1. Estadística Básica/Boxplot por grupos/significativas/')
  if (i == 'Sinovitis ecográfica'){
    estadistica <- read_excel(paste0("results/",i,"/1. Estadística Básica/estadistica_sinovitis.xlsx"))
  } else {
    estadistica <- read_excel(paste0("results/",i,"/1. Estadística Básica/estadistica_microcristales.xlsx"))
  }
  estadistica<-filter(estadistica,col != 'leucocito' & col != 'edad')
  
  
  if (i == 'Sinovitis ecográfica'){
    estadistica_filter_01<-filter(estadistica,test<=0.01)
    estadistica_filter_05<-filter(estadistica,test>0.01 & test<=0.05)
    bands <- list(p_05_01 = estadistica_filter_05,
                  p_01 = estadistica_filter_01)
  } else {
    estadistica_filter_001<-filter(estadistica,test<=0.001)
    estadistica_filter_05_001<-filter(estadistica,test>0.001 & test<=0.05)
    bands <- list(p_05_001 = estadistica_filter_05_001,
                  p_001 = estadistica_filter_001)
  }
  
  for(band in names(bands)){
    x <- bands[[band]]
    n <- nrow(x)
    cols <- ceiling(sqrt(n))
    rows <- ceiling(n / cols)
    tile_str <- paste0(cols, "x", rows)

    files_filter<-paste0(path_sign,x$col,'_boxplot.svg')

    imgs <- lapply(files_filter, function(f) {
      image_read_svg(f)
    })

    combined <- image_montage(
      image_join(imgs),
      tile = tile_str,
      geometry = "+12+12",
      bg = "white"
    )

    image_write(image_convert(combined, "jpeg"),
                path   = file.path(path_sign, paste0("combinado_auto_",band,".jpg")),
                format = "jpeg",
                quality = 100,
                density = "300x300")
  }
}
