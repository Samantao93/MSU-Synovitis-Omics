# R - Estadística por grupo

run_inferential_statistics <- function(df, group_col, output_path) {

  # Eliminar otras columnas de agrupación
  group_vars <- c("deposito", "microcristales", "sinovitis")
  vars_to_remove <- setdiff(group_vars, group_col)
  df <- df[, !(colnames(df) %in% c('cod_paciente',vars_to_remove))]

  columns <- colnames(df)[!colnames(df) %in% c(group_col, "cod_paciente")]

  resultados <- lapply(columns, function(col) {
    if (names(df[col]) %in% list_qual) {
      tryCatch({
        qualitative_data(col, group_col)
      }, error = function(e) {
        message("Error en cualitativa ", col, ": ", e$message)
        return(NULL)
      })
    } else {
      tryCatch({
        quantitative_data(col, group_col)
      }, error = function(e) {
        message("Error en cuantitativa ", col, ": ", e$message)
        return(NULL)
      })
    }
  })

  resultados_validos <- Filter(function(x) !is.null(x) && nrow(x) > 0, resultados)
  final_result <- do.call(plyr::rbind.fill, resultados_validos)

  writexl::write_xlsx(final_result, output_path)
  message("Estadística inferencial guardada en: ", output_path)
  assign(paste0("final_result_", group_col), final_result, envir = .GlobalEnv)
}
