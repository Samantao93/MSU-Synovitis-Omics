# MSU-Synovitis-Omics

# Repository structure

```
MSU-Synovitis-Omics/
├── 0.run_analysis.R            # Orchestrates the full pipeline
├── 1.library_functions.R       # Helper/library functions used by the pipeline
├── 2.original_to_lod.R         # Data preprocessing (e.g., handling LOD)
├── 3.descriptive_statistics.R  # Summary stats
├── 4.inferential_statistics.R  # Hypothesis testing / modeling
├── 5.boxplot.R                 # Boxplot generation
├── 6.volcanoplot.R             # Volcano plot generation
├── 7.PCA.R                     # Principal Component Analysis
├── 8.corr_heat.R               # Correlation heatmap
├── Definitive.ipynb            # Interactive analysis notebook
├── README.md                   # You are here
└── environment.yml             # (Optional) Conda env for Python side (for reproducibility)
```

# Requirements

## R Packages (R 4.5.0)

| Package        | Version |
|----------------|----------|
| dplyr          | 1.1.4    |
| EnhancedVolcano| 1.26.0   |
| FactoMineR     | 2.11     |
| factoextra     | 1.0.7    |
| ggplot2        | 3.5.2    |
| ggpubr         | 0.6.1    |
| gridExtra      | 2.3      |
| limma          | 3.64.1   |
| readr          | 2.1.5    |
| readxl         | 1.4.5    |
| tidyr          | 1.3.1    |

---

## Python Packages (Python 3.12)

| Package       | Version |
|----------------|----------|
| catboost      | 1.2.8    |
| lightgbm      | 4.6.0    |
| matplotlib    | 3.10.0   |
| numpy         | 2.0.1    |
| pandas        | 2.2.3    |
| scikit-learn  | 1.6.1    |
| seaborn       | 0.13.2   |
| xgboost       | 3.0.0    |

# Steps of the pipeline (R)

- Data preprocessing (2.original_to_lod.R)
Converts original measurements to a format suitable for analysis, including handling values near/below the Limit of Detection (LOD).

- Descriptive statistics (3.descriptive_statistics.R)
Computes summary metrics to characterize cohorts/conditions.

- Inferential statistics (4.inferential_statistics.R)
Performs group comparisons (e.g., differential expression/abundance) and produces tables for downstream plotting.

- Visualizations

  - Boxplots (5.boxplot.R): distribution comparisons by group.
  - Volcano plots (6.volcanoplot.R): effect size vs significance for rapid hit selection (uses limma/EnhancedVolcano).
  - PCA (7.PCA.R): global structure and separation (via FactoMineR/factoextra).
  - Correlation heatmaps (8.corr_heat.R): feature–feature or sample–sample relationships.
 
# Steps of the pipeline (Python)
- Data preprocessing
  - Removes leakage columns.
  - Fills missing values with the column median.
  - Drops rows with missing targets and the first (ID) column.

- Models evaluated
Logistic Regression, Random Forest, AdaBoost, SVM, Gaussian Naive Bayes, KNN, XGBoost, LightGBM, and CatBoost.

- Feature selection methods
all, f_classif, chi2, RFE, lasso, and ridge — each tested with 5, 10, 15, and 20 selected features.

- Metrics reported
Accuracy, F1 (macro/weighted), Precision/VPP (macro), Recall (macro/weighted), AUC (macro/weighted; OVO), Cohen’s Kappa, and VPN.
