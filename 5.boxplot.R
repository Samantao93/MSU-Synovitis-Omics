generate_boxplots <- function(data, group_col, output_dir) {
  path_out <- file.path(output_dir)
  path_sign <- file.path(path_out, "significativas")
  
  # Cargar resultados inferenciales
  final_result <- get(paste0("final_result_", group_col), envir = .GlobalEnv)
  
  # Definir columnas a eliminar según el grupo
  all_group_vars <- c("deposito", "microcristales", "sinovitis", "suero", "sexo")
  cols_to_remove <- setdiff(all_group_vars, group_col)
  
  # Filtrar datos: quitar NA en la variable de agrupación
  data <- data[!is.na(data[[group_col]]), ]
  data[[group_col]] <- factor(data[[group_col]])
  
  if (group_col == "deposito") {
    labels <- c("0" = "No", "1" = "Yes")
    colors <- c("0" = "#999999", "1" = "#E69F00")
    legend_title <- "Deposits"
  } else if (group_col == "sinovitis") {
    labels <- c("0" = "No", "1" = "Yes")
    colors <- c("0" = "#999999", "1" = "#E69F00")
    legend_title <- "Synovitis"
  } else if (group_col == "microcristales") {
    labels <- c("0" = "NO", "1" = "MSU", "2" = "CPP")
    colors <- c("0" = "#999999", "1" = "#E69F00", "2" = "#56B4E9")
    legend_title <- "Microcrystals"
  } else {
    labels <- waiver()
    colors <- NULL
    legend_title <- waiver()
  }
  
  # Variables numéricas para graficar
  numeric_vars <- colnames(data)[
    sapply(data, function(x) is.numeric(x) || is.integer(x)) &
      !(colnames(data) %in% c(group_col, "cod_paciente", "Assay", cols_to_remove))
  ]
  
  for (i in numeric_vars) {
    # Obtener p-valor si existe
    if (i %in% final_result$col) {
      row_info <- final_result %>% filter(col == i)
      p_value <- row_info$test
    } else {
      p_value <- NA
    }
    
    # Etiqueta de significancia
    signif_label <- ifelse(is.na(p_value), "ns",
                           ifelse(p_value < 0.001, "***",
                                  ifelse(p_value < 0.01, "**",
                                         ifelse(p_value < 0.05, "*", "ns"))))
    
    # Ruta de guardado según significancia
    output_file <- if (!is.na(p_value) && p_value < 0.05) {
      file.path(path_sign, paste0(i, "_boxplot.png"))
    } else {
      file.path(path_out, paste0(i, "_boxplot.png"))
    }
    
    # Coordenada y para anotación
    # Coordenada y para anotación y ajuste del eje Y
    max_val <- max(data[[i]], na.rm = TRUE)
    min_val <- min(data[[i]], na.rm = TRUE)
    annot_y <- max_val + 0.05 * abs(max_val)
    y_lower_limit <- min_val - 5
    
    
    # Crear gráfico base
    p <- ggplot(data, aes_string(x = group_col, y = i, colour = group_col)) +
      geom_boxplot() +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      scale_color_manual(values = colors, labels = labels) +
      scale_y_continuous(
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.15))
      ) +
      labs(
        title = paste0(i),
        x = group_col,
        y = paste0("NPX"),
        colour = legend_title
      )
    
    levels_group <- levels(data[[group_col]])
    
    if (group_col == "microcristales") {
      # Comparaciones posthoc si existen y significativas
      if (!is.na(p_value) && p_value < 0.05) {
        posthoc_str <- row_info$posthoc
        if (!is.na(posthoc_str)) {
          comparisons <- strsplit(posthoc_str, ";")[[1]]
          comparisons <- trimws(comparisons)
          
          parsed_comparisons <- list()
          annotations <- c()
          
          for (comp in comparisons) {
            parts <- strsplit(comp, ":")[[1]]
            groups <- trimws(strsplit(parts[1], "vs")[[1]])
            pval <- as.numeric(trimws(parts[2]))
            
            if (!is.na(pval) && pval < 0.05) {
              signif_label <- ifelse(pval < 0.001, "***",
                                     ifelse(pval < 0.01, "**",
                                            ifelse(pval < 0.05, "*", "ns")))
              
              parsed_comparisons <- append(parsed_comparisons, list(groups))
              annotations <- c(annotations, signif_label)
            }
          }
          
          # Añadir solo si hay comparaciones significativas
          if (length(parsed_comparisons) > 0) {
            p <- p + ggsignif::geom_signif(
              comparisons = parsed_comparisons,
              annotations = annotations,
              step_increase = 0.1,
              tip_length = 0.01,
              textsize = 5
            )
          }
        }
      }
    } else {
      # Comparación general entre extremos
      if (length(levels_group) >= 2) {
        comparison <- list(c(levels_group[1], levels_group[length(levels_group)]))
        p <- p + ggsignif::geom_signif(
          comparisons = comparison,
          annotations = signif_label,
          tip_length = 0.01,
          y_position = annot_y,
          textsize = 5
        )
      }
    }
    
    # Guardar gráfico
    ggsave(filename = output_file, plot = p, width = 7, height = 5, dpi = 300)
  }
  
  message("✅ Boxplots guardados en: ", path_out, " y ", path_sign)
}