generate_volcano_plots <- function(data, group_col, output_dir) {
  # Eliminar únicamente los NA de la columna de agrupación
  data <- data[!is.na(data[[group_col]]), ]
  
  # Lista de columnas que quieres eliminar, excepto si coinciden con group_col
  all_group_vars <- c("cod_paciente", "deposito", "microcristales", "sinovitis", "suero", "sexo", "leucocito", "edad")
  cols_to_remove <- setdiff(all_group_vars, group_col)
  
  # Eliminar esas columnas
  data2 <- dplyr::select(data, -any_of(cols_to_remove))
  
  # Convertir a factor y generar nombres válidos
  data2[[group_col]] <- factor(data2[[group_col]])
  group_levels <- levels(data2[[group_col]])
  valid_names <- make.names(group_levels)
  
  # Crear diseño del modelo
  design <- model.matrix(~ 0 + data2[[group_col]])
  colnames(design) <- valid_names
  
  levels_list <- colnames(design)
  contrast_list <- list()
  contrast_names <- c()
  
  if (length(valid_names) == 2) {
    contrast_list[[1]] <- paste0(valid_names[2], " - ", valid_names[1])
    contrast_names <- paste0("Group", group_levels[2], "_vs_", group_levels[1])
  } else {
    for (i in 2:length(valid_names)) {
      for (j in 1:(i - 1)) {
        contrast_list <- c(contrast_list, paste0(valid_names[j], " - ", valid_names[i]))
        contrast_names <- c(contrast_names, paste0("Group", group_levels[j], "_vs_", group_levels[i]))
      }
    }
  }
  
  contrast.matrix <- makeContrasts(contrasts = unlist(contrast_list), levels = design)
  colnames(contrast.matrix) <- contrast_names
  
  expr_data <- t(data2[, !(names(data2) %in% group_col)])
  
  fit <- lmFit(expr_data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Diccionarios de etiquetas legibles para valores de agrupación
  label_map <- list(
    sinovitis = c("X0" = "No", "X1" = "Yes"),
    deposito = c("X0" = "No", "X1" = "Yes"),
    microcristales = c("X0" = "NO", "X1" = "MSU", "X2" = "CPP"),
    microcristales02 = c("X0" = "NO+CPP", "X1" = "MSU")
  )
  
  # Función para traducir nombres del contraste a etiquetas legibles
  format_group_label <- function(group_var, level) {
    if (group_var %in% names(label_map)) {
      translated <- label_map[[group_var]][as.character(level)]
      if (!is.null(translated)) return(translated)
    }
    return(level)  # Devuelve sin traducir si no está en el diccionario
  }
  
  for (k in seq_along(contrast_names)) {
    res <- topTable(fit2, coef = contrast_names[k], number = Inf)
    res$minus_log10_pvalue <- -log10(res$P.Value)
    res <- rownames_to_column(res, "Gene")
    
    # Exportar CSV de resultados brutos
    out_csv_raw <- file.path(output_dir, paste0(contrast_names[k], "_raw.csv"))
    write.table(res, out_csv_raw, row.names = FALSE, sep = ";", dec = ",", na = "")
    
    res_sign <- filter(res, adj.P.Val < 0.05)
    
    # Exportar CSV de resultados significativos
    out_csv <- file.path(output_dir, paste0(contrast_names[k], "_significant.csv"))
    write.table(res_sign, out_csv, row.names = FALSE, sep = ";", dec = ",", na = "")
    
    # 🔹 Función interna para generar el volcano plot
    plot_volcano <- function(res, limits_x = NULL, suffix = "") {
      # Extraer niveles originales del contraste (ej. "1 - 0")
      raw_contrast <- contrast_list[[k]]
      levels_split <- unlist(strsplit(raw_contrast, " - "))
      label1 <- format_group_label(group_col, levels_split[1])
      label2 <- format_group_label(group_col, levels_split[2])
      plot_title <- paste("Volcano Plot:", label1, "vs", label2)
      
      p <- ggplot(res, aes(x = logFC, y = minus_log10_pvalue)) +
        geom_point(aes(color = abs(logFC) > 1 & adj.P.Val < 0.05,
                       shape = adj.P.Val < 0.05),
                   size = 2, alpha = 0.7) +
        scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
        geom_text_repel(aes(label = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, as.character(Gene), "")),
                        size = 3, box.padding = 0.5, point.padding = 0.5,
                        segment.color = 'grey50', max.overlaps = Inf) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right"
        ) +
        labs(
          title = plot_title,
          x = "Log2 Fold Change",
          y = "-Log10(p-value)",
          color = "Significant",
          shape = "Significant"
        )
      
      if (!is.null(limits_x)) {
        p <- p + scale_x_continuous(limits = limits_x)
      }
      
      ggsave(
        filename = file.path(output_dir, paste0(contrast_names[k], "_volcano_plot", suffix, ".png")),
        plot = p, width = 10, height = 8, dpi = 300
      )
    }
    
    # 🔹 Guardar volcano plots
    plot_volcano(res, limits_x = NULL, suffix = "")
    plot_volcano(res, limits_x = c(-10, 10), suffix = "_limited")
  }
  
  message("✅ Volcano plots generados en ", output_dir)
}