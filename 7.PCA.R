run_pca_analysis <- function(group_col, output_dir) {
  message(glue("✅ Iniciando análisis PCA para: {group_col}"))
  
  # Diccionarios de etiquetas legibles
  label_map <- list(
    sinovitis = c("0" = "No", "1" = "Yes"),
    deposito = c("0" = "No", "1" = "Yes"),
    microcristales = c("0" = "NO", "1" = "MSU", "2" = "CPP")
  )
  
  # color_map <- list(
  #   sinovitis = c("0" = "#999999", "1" = "#E69F00"),
  #   deposito = c("0" = "#999999", "1" = "#E69F00"),
  #   microcristales = c("0" = "#999999", "1" = "#E69F00", "2" = "#56B4E9")
  # )
  
  legend_titles <- list(
    sinovitis = "Synovitis",
    deposito = "Deposits",
    microcristales = "Microcrystals"
  )
  
  format_label <- function(group_var, values) {
    if (group_var %in% names(label_map)) {
      out <- label_map[[group_var]][as.character(values)]
      out[is.na(out)] <- as.character(values[is.na(out)])
      return(out)
    }
    return(as.character(values))
  }
  
  group_labels <- function(values) {
    format_label(group_col, values)
  }
  
  # group_colors <- if (group_col %in% names(color_map)) color_map[[group_col]] else NULL
  legend_title <- if (group_col %in% names(legend_titles)) legend_titles[[group_col]] else waiver()
  
  # --- Configuración ---
  all_group_vars <- c("deposito", "microcristales", "sinovitis", "suero", "sexo", "leucocito", "edad")
  results_path <- file.path("resultados", group_col, glue("estadistica_{group_col}.xlsx"))
  
  # --- Datos totales ---
  data_filtrada <- combinacion %>% filter(!is.na(.data[[group_col]]))
  cols_to_remove <- setdiff(all_group_vars, group_col)
  
  mat <- data_filtrada %>%
    dplyr::select(-any_of(cols_to_remove)) %>%
    relocate(cod_paciente) %>%
    column_to_rownames("cod_paciente") %>%
    scale() %>%
    as.data.frame() %>%
    dplyr::select(where(~ !all(is.na(.))))
  
  metadatos <- data_filtrada %>%
    dplyr::select(cod_paciente, all_of(group_col)) %>%
    column_to_rownames("cod_paciente") %>%
    .[rownames(mat), , drop = FALSE]
  
  group_vals <- metadatos[[group_col]]
  group_labels_vals <- format_label(group_col, group_vals)
  
  res.pca <- PCA(mat, graph = FALSE)
  
  contrib_total <- as.data.frame(res.pca$var$contrib) %>%
    round(2) %>%
    tibble::rownames_to_column("variable")
  
  write.table(contrib_total,
              file = file.path(output_dir, '1. Con todos', glue("contribuciones_variables_total_{group_col}.csv")),
              na = "", row.names = FALSE, sep = ";")
  
  fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
  ggsave(file.path(output_dir, "1. Con todos", glue("fviz_screeplot_total_{group_col}.png")),
         width = 1920, height = 1920, units = "px")
  
  fviz_pca_ind(res.pca,
               col.ind = group_labels_vals,
               # palette = group_colors,
               addEllipses = TRUE,
               legend.title = legend_title)
  ggsave(file.path(output_dir, "1. Con todos", glue("fviz_pca_ind_total_{group_col}.png")),
         width = 1920, height = 1920, units = "px")
  
  message("✅ PCA total completado.")
  
  # --- PCA con variables significativas ---
  if (!file.exists(results_path)) {
    warning(glue("❌ Archivo de resultados no encontrado: {results_path}"))
    return(invisible(NULL))
  }
  
  message("✅ Iniciando PCA con variables significativas (p < 0.05)")
  
  results_sign <- read_excel(results_path) %>% filter(test < 0.05)
  col_sig <- results_sign$col
  col_sig <- col_sig[col_sig %in% colnames(data_filtrada)]
  col_sig <- col_sig[sapply(data_filtrada[, col_sig], is.numeric)]
  
  if (length(col_sig) > 2) {
    mat_filtrado <- data_filtrada %>%
      dplyr::select(cod_paciente, all_of(col_sig)) %>%
      dplyr::select(-any_of(cols_to_remove)) %>%
      relocate(cod_paciente) %>%
      column_to_rownames("cod_paciente") %>%
      scale() %>%
      as.data.frame() %>%
      dplyr::select(where(~ !all(is.na(.))))
    
    metadatos_filt <- data_filtrada %>%
      dplyr::select(cod_paciente, all_of(group_col)) %>%
      column_to_rownames("cod_paciente") %>%
      .[rownames(mat_filtrado), , drop = FALSE]
    
    group_vals_filt <- metadatos_filt[[group_col]]
    group_labels_filt <- format_label(group_col, group_vals_filt)
    
    res.pca.filtrado <- PCA(mat_filtrado, graph = FALSE)
    
    contrib_filtrado <- as.data.frame(res.pca.filtrado$var$contrib) %>%
      round(2) %>%
      tibble::rownames_to_column("variable")
    
    write.table(contrib_filtrado,
                file = file.path(output_dir, "2. Significativos por estadística", glue("contribuciones_variables_filtrado_{group_col}.csv")),
                na = "", row.names = FALSE, sep = ";")
    
    fviz_screeplot(res.pca.filtrado, addlabels = TRUE, ylim = c(0, 70))
    ggsave(file.path(output_dir, "2. Significativos por estadística", glue("fviz_screeplot_filtrado_{group_col}.png")),
           width = 1920, height = 1920, units = "px")
    
    fviz_pca_ind(res.pca.filtrado,
                 col.ind = group_labels_filt,
                 # palette = group_colors,
                 addEllipses = TRUE,
                 legend.title = legend_title)
    ggsave(file.path(output_dir, "2. Significativos por estadística", glue("fviz_pca_ind_filtrado_{group_col}.png")),
           width = 1920, height = 1920, units = "px")
    
    message("✅ PCA filtrado completado.")
  } else {
    warning("❌ No hay suficientes variables numéricas significativas para PCA filtrado.")
  }
  
  message(glue("🎯 Análisis PCA finalizado para: {group_col}"))
}
