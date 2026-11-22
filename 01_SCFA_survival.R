#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
})

if (!requireNamespace("writexl", quietly = TRUE)) {
  stop("The writexl package is required.")
}

args <- commandArgs(FALSE)
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', args, value=TRUE)), fixed=TRUE))
root <- dirname(script)
input <- Sys.getenv('VALERIC_DATA_DIR', file.path(root, 'data', 'private'))
input <- normalizePath(input, mustWork=TRUE)
input_dir <- input
output_dir <- file.path(root, 'outputs', 'tables')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
metadata <- read_csv(
  file.path(input_dir, "metadata.csv"),
  show_col_types = FALSE
)
scfa <- read_csv(
  file.path(input_dir, "SCFA_log2.csv"),
  show_col_types = FALSE
)
analysis_data <- inner_join(metadata, scfa, by = "ResearchID")

outcomes <- tibble::tribble(
  ~outcome, ~outcome_label, ~time_variable, ~event_variable,
  "therapy_intensification", "Overall therapy intensification",
  "time_to_therapy_intensification", "event_therapy_intensification",
  "medication_intensification", "Medication intensification",
  "time_to_medication_intensification", "event_medication_intensification",
  "intensive_escalation", "Intensive escalation",
  "time_to_intensive_escalation_m", "event_intensive_escalation"
)

scfa_lookup <- tibble(
  feature = setdiff(names(scfa), "ResearchID")
) %>%
  mutate(
    SCFA = case_when(
      feature == "X2.Methylbutyric_acid" ~ "2-Methylbutyric acid",
      TRUE ~ gsub("_", " ", feature)
    )
  )

adjustment_variables <- c(
  "diagnosis_at_sampling",
  "age_at_sampling",
  "Sex",
  "BMI_at_sampling",
  "resection_prior_sampling",
  "ActiveDisease_at_sampling",
  "clinical_BowelMovementADayDef"
)

median_results <- list()
continuous_results <- list()
row_index <- 1L

for (outcome_index in seq_len(nrow(outcomes))) {
  outcome_row <- outcomes[outcome_index, ]

  for (feature_name in scfa_lookup$feature) {
    # Median-stratified Cox and log-rank analyses underlying Figure 2 and
    # Supplementary Figure S1. The HR compares High with Low concentration.
    median_data <- analysis_data %>%
      transmute(
        time = .data[[outcome_row$time_variable]],
        event = as.numeric(.data[[outcome_row$event_variable]]),
        SCFA_value = .data[[feature_name]]
      ) %>%
      filter(complete.cases(.))

    median_cutpoint <- median(median_data$SCFA_value)
    median_data <- median_data %>%
      mutate(
        group = factor(
          if_else(SCFA_value >= median_cutpoint, "High", "Low"),
          levels = c("Low", "High")
        )
      )

    median_cox <- coxph(
      Surv(time, event) ~ group,
      data = median_data,
      ties = "efron"
    )
    median_summary <- summary(median_cox)
    logrank <- survdiff(Surv(time, event) ~ group, data = median_data)
    group_summary <- median_data %>%
      group_by(group) %>%
      summarise(n = n(), events = sum(event), .groups = "drop")

    median_results[[row_index]] <- tibble(
      outcome = outcome_row$outcome,
      outcome_label = outcome_row$outcome_label,
      feature = feature_name,
      median_log2_concentration = median_cutpoint,
      n_total = nrow(median_data),
      events_total = sum(median_data$event),
      n_low = group_summary$n[group_summary$group == "Low"],
      events_low = group_summary$events[group_summary$group == "Low"],
      n_high = group_summary$n[group_summary$group == "High"],
      events_high = group_summary$events[group_summary$group == "High"],
      comparison = "High versus Low",
      HR = median_summary$conf.int[1, "exp(coef)"],
      CI_low = median_summary$conf.int[1, "lower .95"],
      CI_high = median_summary$conf.int[1, "upper .95"],
      Cox_p_value = median_summary$coefficients[1, "Pr(>|z|)"],
      logrank_p_value = pchisq(logrank$chisq, df = 1, lower.tail = FALSE)
    )

    # Fully adjusted Cox analysis with SCFA modelled continuously per one-unit
    # increase in log2 concentration.
    continuous_data <- analysis_data %>%
      transmute(
        time = .data[[outcome_row$time_variable]],
        event = as.numeric(.data[[outcome_row$event_variable]]),
        SCFA_value = .data[[feature_name]],
        across(all_of(adjustment_variables))
      ) %>%
      filter(complete.cases(.))

    continuous_formula <- as.formula(
      paste0(
        "Surv(time, event) ~ SCFA_value + ",
        paste(adjustment_variables, collapse = " + ")
      )
    )
    continuous_cox <- coxph(
      continuous_formula,
      data = continuous_data,
      ties = "efron"
    )
    continuous_summary <- summary(continuous_cox)

    continuous_results[[row_index]] <- tibble(
      outcome = outcome_row$outcome,
      outcome_label = outcome_row$outcome_label,
      feature = feature_name,
      n = nrow(continuous_data),
      events = sum(continuous_data$event),
      exposure_scale = "Per 1-unit higher log2 SCFA concentration",
      adjusted_HR = continuous_summary$conf.int[1, "exp(coef)"],
      CI_low = continuous_summary$conf.int[1, "lower .95"],
      CI_high = continuous_summary$conf.int[1, "upper .95"],
      p_value = continuous_summary$coefficients[1, "Pr(>|z|)"],
      adjustment_set = paste(adjustment_variables, collapse = "; ")
    )

    row_index <- row_index + 1L
  }
}

s3a <- bind_rows(median_results) %>%
  left_join(scfa_lookup, by = "feature") %>%
  group_by(outcome) %>%
  mutate(
    Cox_FDR_within_outcome_8SCFAs = p.adjust(Cox_p_value, method = "BH"),
    logrank_FDR_within_outcome_8SCFAs =
      p.adjust(logrank_p_value, method = "BH")
  ) %>%
  ungroup() %>%
  select(
    outcome, outcome_label, feature, SCFA, median_log2_concentration,
    n_total, events_total, n_low, events_low, n_high, events_high,
    comparison, HR, CI_low, CI_high, Cox_p_value,
    Cox_FDR_within_outcome_8SCFAs, logrank_p_value,
    logrank_FDR_within_outcome_8SCFAs
  ) %>%
  arrange(outcome, SCFA)

s3b <- bind_rows(continuous_results) %>%
  left_join(scfa_lookup, by = "feature") %>%
  group_by(outcome) %>%
  mutate(FDR_within_outcome_8SCFAs = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  select(
    outcome, outcome_label, feature, SCFA, n, events, exposure_scale,
    adjusted_HR, CI_low, CI_high, p_value, FDR_within_outcome_8SCFAs,
    adjustment_set
  ) %>%
  arrange(outcome, SCFA)

workbook_path <- file.path(
  output_dir,
  "Supplementary_Table_S3_SCFA_survival_analysis.xlsx"
)
writexl::write_xlsx(
  list(
    S3A_median_stratified = s3a,
    S3B_continuous_adjusted = s3b
  ),
  workbook_path
)

write_csv(
  s3a,
  file.path(output_dir, "S3A_SCFA_median_stratified_survival.csv")
)
write_csv(
  s3b,
  file.path(output_dir, "S3B_SCFA_continuous_adjusted_Cox.csv")
)

message("Supplementary Table S3 written to: ", workbook_path)

write_csv(s3a, file.path(output_dir, 'Supplementary_Table_S2_KM_median_stratified_SCFA.csv'))
write_csv(s3b, file.path(output_dir, 'Supplementary_Table_S3_continuous_SCFA_Cox.csv'))
capture.output(sessionInfo(), file=file.path(output_dir, 'sessionInfo.txt'))
