library(readxl)
library(dplyr)
library(nortest)
library(tibble)
library(GGally)
library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(FactoMineR)
library(factoextra)
library(plyr)
library(rstatix)
library(writexl)
library(MASS)
library(caret)
library(stringr)
library(ggrepel)
library(ggsignif)
library(glue)
library(reshape2)
library(irr)
library(magick)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# Lista de variables de agrupación
group_vars <- c("Deposito articular", "Microcristales", "Sinovitis ecográfica")

# Estructura de carpetas por cada variable
subfolders <- list(
  "1. Estadística Básica/Boxplot por grupos/significativas",
  "2. Volcano Plot",
  "3. PCA/1. Con todos",
  "3. PCA/2. Significativos por estadística"
)

# Carpeta raíz donde guardar todo
root_dir <- "./results"

# Crear carpetas
for (group in group_vars) {
  for (sub in subfolders) {
    dir_path <- file.path(root_dir, group, sub)
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

message("✅ Estructura de carpetas creada correctamente en: ", normalizePath(root_dir))

############## Funciones ##############
calc_singrupo<-function(x,list_qual) {
  # Dataframe with result
  empty_df <- data.frame(
    column = character(),
    total = numeric(),
    stat = character(),
    medida = numeric(),
    sd_iq = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 2:ncol(x)){
    # print(i)
    if(names(x[i]) %in% list_qual){
      column_quali<-paste0(names(x[i]),'*')
      total_quali<-sum(!is.na(x[[i]]))
      stat_quali<-'Frec. Absolute/Relative'
      medida_quali<-table(x[[i]])
      sd_quali<-round((prop.table(table(x[[i]])))*100,2)
      
      res <- data.frame(
        column=column_quali,
        total=total_quali,
        stat=stat_quali,
        freqAbs=medida_quali,
        freqRel=sd_quali,
        stringsAsFactors = TRUE
      )
      empty_df<-rbind.fill(empty_df,res)
      
      
    } else {
      if (length(unique(x[[i]])) <= 1) {
        stat <- "Constante"
        medida <- round(unique(x[[i]])[1], 2)
        sd_iq <- 0
      } else if (lillie.test(x[[i]])$p.value < 0.05) {
        stat <- "Mediana (Q1-Q3)"
        medida <- round(median(x[[i]], na.rm = TRUE), 2)
        quantiles <- quantile(x[[i]], probs = c(0.25, 0.75), na.rm = TRUE)
        q25 <- round(quantiles[1], 2)
        q75 <- round(quantiles[2], 2)
        sd_iq <- paste0(q25, " - ", q75)
      } else {
        stat <- "Media (sd)"
        medida <- round(mean(x[[i]], na.rm = TRUE), 2)
        sd_iq <- round(sd(x[[i]], na.rm = TRUE), 2)
      }
      column<-names(x[i])
      total<-sum(!is.na(x[[i]]))
      
      res2 <- data.frame(
        column,
        total,
        stat,
        medida,
        sd_iq,
        stringsAsFactors = FALSE
      )
      empty_df<-rbind.fill(empty_df,res2)
    }
  }
  return(empty_df)
}

### Agrupando por cristales
calculate_stat <- function(x) {
  n_valid <- sum(!is.na(x))
  
  if (n_valid <= 4) {
    return(c(total=n_valid,stat = "Mediana",medida=round(median(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2)))
  } else {
    test_result <- lillie.test(x)
    
    if (test_result$p.value < 0.05) {
      return(c(total=n_valid,stat = "Mediana",medida= round(median(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2)))
    } else {
      return(c(total=n_valid,stat = "Media",medida= round(mean(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2)))
    }
  }
}

calculate_factor <- function(p) {
  freq<-table(p)
  absFreq<-as.integer(freq)
  relFreq<-round(100*prop.table(freq),2)
  return(list(c(AbsFreq=absFreq,RelFreq=relFreq)))
}

filtrar_columnas_menos4 <- function(df, threshold = 5) {
  filtered_df <- df[, sapply(df, function(col) sum(!is.na(col)) < threshold)]
  return(filtered_df)
}

filtrar_columnas_mas4 <- function(df, threshold = 5) {
  filtered_df <- df[, sapply(df, function(col) sum(!is.na(col)) > threshold)]
  return(filtered_df)
}


qualitative_data <- function(col, group, data = combinacion) {
  print(paste("Procesando columna:", col))
  
  contingency_table <- table(data[[group]], data[[col]])
  expected <- suppressWarnings(chisq.test(contingency_table)$expected)
  
  if (any(expected < 5)) {
    fisher_result <- fisher.test(contingency_table)
    return(data.frame(
      col = col,
      test_type = "Fisher's Exact",
      test = round(fisher_result$p.value,3),
      posthoc = NA,
      parametric = NA,
      message = NA,
      stringsAsFactors = FALSE
    ))
  } else {
    chi_result <- chisq.test(contingency_table)
    group_levels <- length(unique(data[[group]]))
    
    if (chi_result$p.value < 0.05 && group_levels > 2) {
      comparisons <- rownames(pairwise_result$p.value)
      posthoc_df <- data.frame(
        comparison = rep(comparisons, times = ncol(pairwise_result$p.value)),
        group = rep(colnames(pairwise_result$p.value), each = length(comparisons)),
        p.value = as.vector(round(pairwise_result$p.value,3))
      )
      posthoc_df <- na.omit(posthoc_df)
      posthoc_df <- posthoc_df[posthoc_df$comparison != posthoc_df$group, ]
      posthoc_df$comparison <- paste(posthoc_df$comparison, "vs", posthoc_df$group)
      posthoc_df <- unique(posthoc_df[, c("comparison", "p.value")])
      posthoc_text <- paste(apply(posthoc_df, 1, paste, collapse = ": "), collapse = "; ")
      
    } else {
      posthoc_text <- NA
    }
    
    return(data.frame(
      col = col,
      test_type = "Chi-squared",
      test = round(chi_result$p.value,3),
      posthoc = posthoc_text,
      parametric = NA,
      message = NA,
      stringsAsFactors = FALSE
    ))
  }
}

quantitative_data <- function(col, group, data = combinacion, posthoc_method = "bonferroni") {
  print(paste("Procesando columna:", col))
  
  data_clean <- na.omit(data[, c(col, group)])
  group_levels <- unique(data_clean[[group]])
  num_groups <- length(group_levels)
  
  is_normal <- sapply(group_levels, function(g) {
    vals <- data_clean[data_clean[[group]] == g, col]
    if (length(vals) >= 3) {
      shapiro.test(vals)$p.value > 0.05
    } else FALSE
  })
  
  parametric <- all(is_normal)
  has_ties <- sum(duplicated(data_clean[[col]])) > 0
  
  if (num_groups == 2) {
    if (parametric) {
      test <- t.test(reformulate(group, col), data = data_clean)
      method <- "t-test"
    } else {
      test <- wilcox.test(reformulate(group, col), data = data_clean, exact = !has_ties)
      method <- "Mann-Whitney U"
    }
    return(data.frame(
      col = col,
      test_type = method,
      test = round(test$p.value,3),
      posthoc = NA,
      parametric = parametric,
      message = NA,
      stringsAsFactors = FALSE
    ))
  } else if (num_groups > 2) {
    if (parametric) {
      model <- aov(reformulate(group, col), data = data_clean)
      summary_result <- summary(model)
      p_val <- summary_result[[1]][["Pr(>F)"]][1]
      posthoc <- if (!is.na(p_val) && p_val < 0.05) {
        ph <- TukeyHSD(model)
        paste(capture.output(print(ph)), collapse = "\n")
      } else NA
      method <- "ANOVA"
      return(data.frame(
        col = col,
        test_type = method,
        test = round(p_val,3),
        posthoc = posthoc,
        parametric = parametric,
        message = NA,
        stringsAsFactors = FALSE
      ))
    } else {
      test <- kruskal.test(reformulate(group, col), data = data_clean)
      posthoc <- if (test$p.value < 0.05) {
        ph <- pairwise.wilcox.test(data_clean[[col]], data_clean[[group]], p.adjust.method = posthoc_method)
        ph_matrix <- round(ph$p.value,3)
        comparisons <- rownames(ph_matrix)
        posthoc_df <- data.frame(
          comparison = rep(comparisons, times = ncol(ph_matrix)),
          group = rep(colnames(ph_matrix), each = length(comparisons)),
          p.value = as.vector(ph_matrix)
        )
        posthoc_df <- na.omit(posthoc_df)
        posthoc_df <- posthoc_df[posthoc_df$comparison != posthoc_df$group, ]
        posthoc_df$comparison <- paste(posthoc_df$comparison, "vs", posthoc_df$group)
        posthoc_df <- unique(posthoc_df[, c("comparison", "p.value")])
        
        posthoc_text <- paste(apply(posthoc_df, 1, paste, collapse = ": "), collapse = "; ")
        
      } else NA
      method <- "Kruskal-Wallis"
      return(data.frame(
        col = col,
        test_type = method,
        test = round(test$p.value,3),
        posthoc = posthoc,
        parametric = parametric,
        message = NA,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    return(data.frame(
      col = col,
      test_type = NA,
      test = NA,
      posthoc = NA,
      parametric = parametric,
      message = "Número de grupos insuficiente",
      stringsAsFactors = FALSE
    ))
  }

}
