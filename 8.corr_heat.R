analizar_synovial_serum <- function(
    df,
    id_col = "cod_paciente",
    serum_suffix = "_serum$",
    out_dir = "results",
    thr = 0.6,                  # umbral Spearman
    top_n = NA_integer_,        # top-N por Spearman
    icc_thr = 0.5,              # umbral ICC
    icc_top_n = NA_integer_     # top-N por ICC
) {
  
  es_valida <- function(x) {
    if (!is.numeric(x)) {
      x_s <- suppressWarnings(as.numeric(as.character(x)))
      if (!all(is.na(x_s))) x <- x_s
    }
    sum(!is.na(x)) >= 3 && stats::sd(x, na.rm = TRUE) > 0
  }
  
  # Heatmap genérico
  heat_pairs <- function(df_pairs, value_col, titulo, archivo_png, limits=c(0,1), midpoint=0.5, text_size = 4) {
    if (nrow(df_pairs) == 0) return(invisible(NULL))
    stopifnot(value_col %in% names(df_pairs))
    plt_df <- data.frame(
      yvar = df_pairs$biomarcador_synovial,
      xvar = df_pairs$biomarcador_serum,
      value = df_pairs[[value_col]],
      label = ifelse(is.finite(df_pairs[[value_col]]), sprintf("%.2f", df_pairs[[value_col]]), "NA"),
      stringsAsFactors = FALSE
    )
    n_x <- length(unique(plt_df$xvar))
    n_y <- length(unique(plt_df$yvar))
    w <- max(8, min(24, n_x * 0.6))
    h <- max(6, min(24, n_y * 0.6))
    p <- ggplot(plt_df, aes(x = xvar, y = yvar, fill = value)) +
      geom_tile(color = "white") +
      geom_text(aes(label = label), size = text_size) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "pink",
                           midpoint = midpoint, limit = limits, na.value = "grey90") +
      theme_minimal(base_size = 12) +
      labs(title = titulo, x = "", y = "", fill = value_col) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))
      )
    ggsave(file.path(out_dir, archivo_png), plot = p, width = w, height = h, dpi = 300)
  }
  
  # Trocea un data.frame en K partes (por filas)
  split_in_k <- function(df_in, k = 5) {
    if (nrow(df_in) == 0) return(list())
    k <- max(1, as.integer(k))
    n <- nrow(df_in)
    chunk <- ceiling(n / k)
    idx_starts <- seq(1, n, by = chunk)
    lapply(seq_along(idx_starts), function(i) {
      i1 <- idx_starts[i]
      i2 <- min(n, i1 + chunk - 1)
      df_in[i1:i2, , drop = FALSE]
    })
  }
  
  # ---------- Selección columnas ----------
  syn_cols   <- df[, !grepl(serum_suffix, names(df)) & !grepl(paste0("^", id_col, "$"), names(df)), drop = FALSE]
  serum_cols <- df[,  grepl(serum_suffix, names(df)), drop = FALSE]
  if (ncol(syn_cols) > 0)   syn_cols   <- syn_cols[,   sapply(syn_cols,   es_valida),  drop = FALSE]
  if (ncol(serum_cols) > 0) serum_cols <- serum_cols[, sapply(serum_cols, es_valida),  drop = FALSE]
  
  mismos <- intersect(names(syn_cols), gsub(serum_suffix, "", names(serum_cols)))
  if (length(mismos) == 0) {
    message("No se encontraron pares m y m_serum válidos.")
    return(invisible(FALSE))
  }
  
  # ---------- Spearman + ICC ----------
  df_pairs <- do.call(rbind, lapply(mismos, function(m) {
    x <- df[[m]]
    y <- df[[paste0(m, "_serum")]]
    if (!is.numeric(x)) x <- suppressWarnings(as.numeric(as.character(x)))
    if (!is.numeric(y)) y <- suppressWarnings(as.numeric(as.character(y)))
    keep <- stats::complete.cases(x, y)
    x2 <- x[keep]; y2 <- y[keep]
    
    ct <- tryCatch(suppressWarnings(cor.test(x2, y2, method = "spearman")), error = function(e) NULL)
    rho <- if (!is.null(ct)) unname(ct$estimate) else NA_real_
    pvl <- if (!is.null(ct)) ct$p.value else NA_real_
    
    icc_val <- NA_real_; icc_lo <- NA_real_; icc_hi <- NA_real_
    if (sum(keep) >= 2) {
      dat <- data.frame(synovial = x2, serum = y2)
      icc_fit <- tryCatch(irr::icc(dat, model = "twoway", type = "agreement", unit = "single"),
                          error = function(e) NULL)
      if (!is.null(icc_fit)) {
        icc_val <- icc_fit$value
        icc_lo  <- icc_fit$lbound
        icc_hi  <- icc_fit$ubound
      }
    }
    
    data.frame(
      biomarcador_synovial = paste0(m, "_synovial"),
      biomarcador_serum    = paste0(m, "_serum"),
      n_usable             = sum(keep),
      correlation          = round(rho, 2),
      p_value              = round(signif(pvl, 3), 3),
      ICC                  = round(icc_val, 2),
      ICC_L95              = round(icc_lo, 2),
      ICC_U95              = round(icc_hi, 2),
      stringsAsFactors     = FALSE
    )
  }))
  
  # --- Guardar y graficar (Spearman) ---
  thr_str <- gsub("\\.","", as.character(thr))
  write.csv(df_pairs, file.path(out_dir, "pairs_spearman_icc_all.csv"), row.names = FALSE)
  
  # --- Guardar y graficar (Spearman) ---
  # si filtras SOLO por p-valor:
  df_pairs_thr <- subset(df_pairs, !is.na(correlation) & !is.na(p_value) & p_value <= 0.05)
  
  # (OPCIONAL) si prefieres p-valor + tamaño de efecto:
  # df_pairs_thr <- subset(df_pairs, !is.na(correlation) & !is.na(p_value) &
  #                        p_value <= 0.05 & abs(correlation) >= thr)
  
  # Orden y top-N
  df_pairs_thr <- df_pairs_thr[order(-df_pairs_thr$correlation), , drop = FALSE]
  if (!is.na(top_n) && is.finite(top_n) && top_n > 0) df_pairs_thr <- head(df_pairs_thr, top_n)
  
  # Guardar CSV con nombre coherente
  write.csv(df_pairs, file.path(out_dir, "pairs_spearman_icc_all.csv"), row.names = FALSE)
  write.csv(df_pairs_thr, file.path(out_dir, "pairs_spearman_p_le_005.csv"), row.names = FALSE)
  
  # Heatmaps GLOBAL en 5 partes (r puede ser negativo)
  parts_cor <- split_in_k(df_pairs, k = 5)
  for (i in seq_along(parts_cor)) {
    heat_pairs(parts_cor[[i]], "correlation",
               titulo = paste0("Synovial vs Serum – Spearman (all) – Parte ", i, "/", length(parts_cor)),
               archivo_png = paste0("heat_cor_mismo_biomarcador_all_part", i, ".png"),
               limits = c(-1, 1), midpoint = 0, text_size = 4.2)
  }
  
  # Heatmap filtrado por p ≤ 0.05 (uno solo)
  if (nrow(df_pairs_thr) > 0) {
    heat_pairs(df_pairs_thr, "correlation",
               titulo = if (!is.na(top_n) && top_n > 0)
                 paste0("Spearman p ≤ 0.05 (Top ", top_n, ")")
               else
                 "Spearman p ≤ 0.05",
               archivo_png = if (!is.na(top_n) && top_n > 0)
                 "heat_cor_mismo_biomarcador_p_le_005_top.png"
               else
                 "heat_cor_mismo_biomarcador_p_le_005.png",
               limits = c(-1, 1), midpoint = 0, text_size = 5)
  }
  
  
  # --- Guardar y graficar (ICC) ---
  icc_str <- gsub("\\.","", as.character(icc_thr))
  df_pairs_icc_thr <- subset(df_pairs, !is.na(ICC) & ICC >= icc_thr)
  df_pairs_icc_thr <- df_pairs_icc_thr[order(-df_pairs_icc_thr$ICC), , drop = FALSE]
  if (!is.na(icc_top_n) && is.finite(icc_top_n) && icc_top_n > 0) df_pairs_icc_thr <- head(df_pairs_icc_thr, icc_top_n)
  write.csv(df_pairs_icc_thr, file.path(out_dir, paste0("pairs_icc_ge_", icc_str, ".csv")), row.names = FALSE)
  
  # Heatmaps GLOBAL de ICC en 5 partes
  parts_icc <- split_in_k(df_pairs, k = 5)
  for (i in seq_along(parts_icc)) {
    heat_pairs(parts_icc[[i]], "ICC",
               titulo = paste0("Synovial vs Serum – ICC (agreement, single) – Parte ", i, "/", length(parts_icc)),
               archivo_png = paste0("heat_icc_mismo_biomarcador_all_part", i, ".png"),
               limits = c(0, 1), midpoint = max(0.75, icc_thr), text_size = 4.2)
  }
  # Heatmap ICC filtrado por umbral (uno solo)
  if (nrow(df_pairs_icc_thr) > 0) {
    heat_pairs(df_pairs_icc_thr, "ICC",
               titulo = if (!is.na(icc_top_n) && icc_top_n > 0)
                 paste0("ICC ≥ ", icc_thr, " (Top ", icc_top_n, ")")
               else paste0("ICC ≥ ", icc_thr),
               archivo_png = if (!is.na(icc_top_n) && icc_top_n > 0)
                 paste0("heat_icc_ge_", icc_str, "_top", icc_top_n, ".png")
               else
                 paste0("heat_icc_ge_", icc_str, ".png"),
               limits = c(max(icc_thr, min(df_pairs_icc_thr$ICC, na.rm = TRUE)), 1),
               midpoint = (icc_thr + 1) / 2,
               text_size = 5)
  } else {
    message(paste0("No hay pares con ICC ≥ ", icc_thr, " para el heatmap filtrado."))
  }
  
  # ---------- ICC estratificado por microcristales / deposito / sinovitis ----------
  strata_vars <- c("microcristales", "deposito", "sinovitis")
  strata_vars <- strata_vars[strata_vars %in% names(df)]
  
  # Comparación robusta por texto
  as_key <- function(x) trimws(as.character(x))
  
  if (length(strata_vars) > 0) {
    icc_strata_list <- list()
    
    for (sv in strata_vars) {
      # Fuerza niveles 0/1/2 en microcristales; en resto, toma niveles presentes (como texto)
      if (sv == "microcristales") {
        sv_levels <- c("0","1","2")
      } else {
        sv_levels <- stats::na.omit(unique(as_key(df[[sv]])))
      }
      
      # (Opcional) diagnóstico de tamaños por nivel
      msg_tab <- table(as_key(df[[sv]]), useNA = "ifany")
      message(sprintf("[Estrato %s] %s",
                      sv,
                      paste(sprintf("%s:%d", names(msg_tab), as.integer(msg_tab)), collapse=" | ")))
      
      for (lvl in sv_levels) {
        df_sub <- df[as_key(df[[sv]]) == lvl, , drop = FALSE]
        message(sprintf("  - %s == %s -> n filas: %d", sv, lvl, nrow(df_sub)))
        if (nrow(df_sub) < 2) next  # ICC necesita >=2 sujetos
        
        sub_pairs <- do.call(rbind, lapply(mismos, function(m) {
          x <- df_sub[[m]]
          y <- df_sub[[paste0(m, "_serum")]]
          if (!is.numeric(x)) x <- suppressWarnings(as.numeric(as.character(x)))
          if (!is.numeric(y)) y <- suppressWarnings(as.numeric(as.character(y)))
          keep <- stats::complete.cases(x, y)
          x2 <- x[keep]; y2 <- y[keep]
          n_use <- sum(keep)
          
          icc_val <- icc_lo <- icc_hi <- NA_real_
          if (n_use >= 2) {
            dat <- data.frame(synovial = x2, serum = y2)
            icc_fit <- tryCatch(irr::icc(dat, model = "twoway", type = "agreement", unit = "single"),
                                error = function(e) NULL)
            if (!is.null(icc_fit)) {
              icc_val <- icc_fit$value
              icc_lo  <- icc_fit$lbound
              icc_hi  <- icc_fit$ubound
            }
          }
          data.frame(
            estrato_var = sv,
            estrato_val = lvl,                       # guardamos como texto
            biomarcador_synovial = paste0(m, "_synovial"),
            biomarcador_serum    = paste0(m, "_serum"),
            n_usable             = n_use,
            ICC                  = round(icc_val, 2),
            ICC_L95              = round(icc_lo, 2),
            ICC_U95              = round(icc_hi, 2),
            stringsAsFactors     = FALSE
          )
        }))
        
        if (!is.null(sub_pairs) && nrow(sub_pairs) > 0) {
          icc_strata_list[[paste0(sv, "==", lvl)]] <- sub_pairs
        }
      }
    }
    
    if (length(icc_strata_list) > 0) {
      df_icc_strata <- dplyr::bind_rows(icc_strata_list)
      
      # Guarda TODO (sin filtrar) para ver qué hay en 0/1/2 aunque no pasen el umbral
      write.csv(df_icc_strata,
                file.path(out_dir, "icc_estratificado_all.csv"),
                row.names = FALSE)
      
      # Filtrado por umbral ICC
      icc_str <- gsub("\\.","", as.character(icc_thr))
      df_icc_strata_thr <- subset(df_icc_strata, !is.na(ICC) & ICC >= icc_thr)
      write.csv(df_icc_strata_thr,
                file.path(out_dir, paste0("icc_estratificado_ge_", icc_str, ".csv")),
                row.names = FALSE)
      
      # (Opcional) Heatmaps por estrato SIN filtrar (para que veas 0 y 2 aunque no superen el umbral)
      split_all <- split(df_icc_strata,
                         interaction(df_icc_strata$estrato_var,
                                     df_icc_strata$estrato_val, drop=TRUE))
      for (nm in names(split_all)) {
        df_s <- split_all[[nm]]
        if (nrow(df_s) == 0) next
        safe_nm <- gsub("[^A-Za-z0-9_\\-]+","_", nm)
        heat_pairs(df_s, "ICC",
                   titulo = paste0("ICC por estrato (ALL) – ", nm),
                   archivo_png = paste0("heat_icc_", safe_nm, "_ALL.png"),
                   limits = c(0, 1), midpoint = max(0.75, icc_thr))
      }
      
      # Heatmaps por estrato FILTRADOS por ICC ≥ icc_thr
      split_thr <- split(df_icc_strata_thr,
                         interaction(df_icc_strata_thr$estrato_var,
                                     df_icc_strata_thr$estrato_val, drop=TRUE))
      for (nm in names(split_thr)) {
        df_s <- split_thr[[nm]]
        if (nrow(df_s) == 0) next
        safe_nm <- gsub("[^A-Za-z0-9_\\-]+","_", nm)
        heat_pairs(df_s, "ICC",
                   titulo = paste0("ICC ≥ ", icc_thr, " por estrato – ", nm),
                   archivo_png = paste0("heat_icc_", safe_nm, "_ge_", icc_str, ".png"),
                   limits = c(0, 1), midpoint = max(0.75, icc_thr))
      }
    } else {
      message("No se pudo calcular ICC estratificado (sin subgrupos válidos).")
    }
  }
  
  
  invisible(TRUE)
}

