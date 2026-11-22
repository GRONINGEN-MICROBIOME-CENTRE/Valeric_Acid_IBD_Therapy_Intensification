#!/usr/bin/env Rscript

required_packages <- c("dplyr", "purrr", "broom", "readr", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(broom)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(FALSE)
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', args, value=TRUE)), fixed=TRUE))
root <- dirname(script)
input <- Sys.getenv('VALERIC_DATA_DIR', file.path(root, 'data', 'private'))
input <- normalizePath(input, mustWork=TRUE)
input_dir <- input
results_dir <- file.path(root, 'outputs', 'microbiome')
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
message("Input directory: ", input_dir)
message("Output directory: ", results_dir)

read_matrix <- function(filename) {
  x <- read_csv(file.path(input_dir, filename), show_col_types = FALSE)
  ids <- x$ResearchID
  stopifnot(!anyNA(ids), !anyDuplicated(ids))
  x <- as.data.frame(x[, setdiff(names(x), "ResearchID"), drop = FALSE])
  rownames(x) <- ids
  x
}

df_scfa_log2 <- read_matrix("SCFA_log2.csv")
df_species_clr <- read_matrix("species_CLR.csv")
df_pathways_clr <- read_matrix("pathways_CLR.csv")
metadata <- as.data.frame(
  read_csv(file.path(input_dir, "metadata.csv"), show_col_types = FALSE)
)
cox_results <- read_csv(
  file.path(root, "outputs", "tables", "Cox_results_metagenomic_features.csv"),
  show_col_types = FALSE
)
pathway_labels <- read_csv(
  file.path(input_dir, "pathway_labels.csv"),
  show_col_types = FALSE
)

stopifnot(!anyDuplicated(metadata$ResearchID))
stopifnot(all(metadata$ResearchID %in% rownames(df_scfa_log2)))
stopifnot(all(metadata$ResearchID %in% rownames(df_species_clr)))
stopifnot(all(metadata$ResearchID %in% rownames(df_pathways_clr)))

df_scfa_log2 <- df_scfa_log2[metadata$ResearchID, , drop = FALSE]
df_species_clr <- df_species_clr[metadata$ResearchID, , drop = FALSE]
df_pathways_clr <- df_pathways_clr[metadata$ResearchID, , drop = FALSE]
rownames(metadata) <- metadata$ResearchID

confounders <- c(
  "ActiveDisease_at_sampling",
  "diagnosis_at_sampling",
  "age_at_sampling",
  "resection_prior_sampling",
  "clinical_BowelMovementADayDef"
)
stopifnot(all(confounders %in% names(metadata)))

metadata <- metadata %>%
  mutate(
    ActiveDisease_at_sampling = factor(ActiveDisease_at_sampling),
    diagnosis_at_sampling = factor(diagnosis_at_sampling),
    resection_prior_sampling = factor(resection_prior_sampling)
  )

qc_summary <- tibble(
  item = c(
    "participants_metadata",
    "participants_SCFA",
    "participants_species",
    "participants_pathways",
    "SCFA_count",
    "species_count",
    "pathway_count",
    "Cox_result_rows"
  ),
  value = c(
    nrow(metadata),
    nrow(df_scfa_log2),
    nrow(df_species_clr),
    nrow(df_pathways_clr),
    ncol(df_scfa_log2),
    ncol(df_species_clr),
    ncol(df_pathways_clr),
    nrow(cox_results)
  )
)
write_csv(qc_summary, file.path(results_dir, "00_qc_summary.csv"))

run_manova_layer <- function(feature_data, data_layer) {
  map_dfr(names(feature_data), function(feature_name) {
    model_data <- data.frame(
      df_scfa_log2,
      feature = feature_data[[feature_name]],
      metadata[, confounders, drop = FALSE],
      check.names = FALSE
    )
    model_data <- model_data[complete.cases(model_data), , drop = FALSE]

    if (nrow(model_data) < 10 || dplyr::n_distinct(model_data$feature) < 2) {
      return(tibble(
        data_layer = data_layer,
        feature = feature_name,
        n = nrow(model_data),
        Pillai = NA_real_,
        F_value = NA_real_,
        df1 = NA_real_,
        df2 = NA_real_,
        p_value = NA_real_
      ))
    }

    response <- paste0(
      "cbind(", paste(names(df_scfa_log2), collapse = ", "), ")"
    )
    model_formula <- as.formula(
      paste(response, "~ feature +", paste(confounders, collapse = " + "))
    )
    fit <- manova(model_formula, data = model_data)
    test <- summary(fit, test = "Pillai")$stats["feature", ]

    tibble(
      data_layer = data_layer,
      feature = feature_name,
      n = nrow(model_data),
      Pillai = unname(test["Pillai"]),
      F_value = unname(test["approx F"]),
      df1 = unname(test["num Df"]),
      df2 = unname(test["den Df"]),
      p_value = unname(test["Pr(>F)"])
    )
  }) %>%
    mutate(p_adjusted_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

message("Running MANOVA for 283 species...")
manova_species <- run_manova_layer(df_species_clr, "Species")
message("Running MANOVA for 310 pathways...")
manova_pathways <- run_manova_layer(df_pathways_clr, "Pathways")

manova_all <- bind_rows(manova_species, manova_pathways) %>%
  arrange(data_layer, p_value)
manova_selected <- manova_all %>%
  filter(p_adjusted_BH < 0.05, Pillai > 0.20)

write_csv(manova_all, file.path(results_dir, "01_MANOVA_all_features.csv"))
write_csv(
  manova_selected,
  file.path(results_dir, "02_MANOVA_selected_FDR005_Pillai020.csv")
)

get_scfa_coefficients <- function(selected_features, feature_data, data_layer) {
  map_dfr(selected_features, function(feature_name) {
    model_data <- data.frame(
      df_scfa_log2,
      feature = feature_data[[feature_name]],
      metadata[, confounders, drop = FALSE],
      check.names = FALSE
    )
    model_data <- model_data[complete.cases(model_data), , drop = FALSE]
    response <- paste0(
      "cbind(", paste(names(df_scfa_log2), collapse = ", "), ")"
    )
    model_formula <- as.formula(
      paste(response, "~ feature +", paste(confounders, collapse = " + "))
    )
    fit <- lm(model_formula, data = model_data)

    broom::tidy(fit) %>%
      filter(term == "feature") %>%
      transmute(
        data_layer = data_layer,
        feature = feature_name,
        SCFA = response,
        n = nrow(model_data),
        beta = estimate,
        standard_error = std.error,
        t_value = statistic,
        p_value = p.value,
        CI_lower = estimate - 1.96 * std.error,
        CI_upper = estimate + 1.96 * std.error
      )
  }) %>%
    group_by(data_layer, feature) %>%
    mutate(FDR_within_feature_8SCFA = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

branch_definitions <- list(nominal_p_lt_0.05 = cox_results %>% filter(p_value < 0.05))


for (branch_name in names(branch_definitions)) {
  message("Running branch: ", branch_name)
  branch_dir <- file.path(results_dir, branch_name)
  dir.create(branch_dir, recursive = TRUE, showWarnings = FALSE)

  cox_pass <- branch_definitions[[branch_name]]
  cox_feature_summary <- cox_pass %>%
    group_by(dataset, feature) %>%
    summarise(
      outcomes_passing = paste(sort(unique(outcome)), collapse = ";"),
      min_cox_p = min(p_value, na.rm = TRUE),
      min_cox_FDR = min(p_adj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(data_layer = if_else(dataset == "species", "Species", "Pathways"))

  selected_features <- manova_selected %>%
    inner_join(
      cox_feature_summary,
      by = c("data_layer", "feature")
    ) %>%
    arrange(data_layer, p_adjusted_BH, desc(Pillai))

  species_features <- selected_features %>%
    filter(data_layer == "Species") %>%
    pull(feature)
  pathway_features <- selected_features %>%
    filter(data_layer == "Pathways") %>%
    pull(feature)

  coefficients <- bind_rows(
    get_scfa_coefficients(species_features, df_species_clr, "Species"),
    get_scfa_coefficients(pathway_features, df_pathways_clr, "Pathways")
  ) %>%
    left_join(
      selected_features %>%
        select(
          data_layer, feature, Pillai, MANOVA_p = p_value,
          MANOVA_FDR = p_adjusted_BH, outcomes_passing, min_cox_p, min_cox_FDR
        ),
      by = c("data_layer", "feature")
    )

  write_csv(selected_features, file.path(branch_dir, "02_MANOVA_Cox_intersection.csv"))
  write_csv(coefficients, file.path(branch_dir, "03_SCFA_specific_coefficients_all.csv"))
}
