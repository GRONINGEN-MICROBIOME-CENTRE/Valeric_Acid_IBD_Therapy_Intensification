#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(FALSE)
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', args, value=TRUE)), fixed=TRUE))
root <- dirname(script)
input <- Sys.getenv('VALERIC_DATA_DIR', file.path(root, 'data', 'private'))
input <- normalizePath(input, mustWork=TRUE)
analysis_dir <- file.path(root, 'outputs')
results_dir <- file.path(root, 'outputs', 'microbiome')
input_dir <- input
output_dir <- file.path(results_dir, 'figure3')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
heatmap_palette <- c(low='#2166AC', mid='#FFFFFF', high='#B2182B')
output_stem <- 'Figure3_SCFA_microbiome_outcome_heatmap'

coefficients <- read_csv(
  file.path(
    results_dir, "nominal_p_lt_0.05",
    "03_SCFA_specific_coefficients_all.csv"
  ),
  show_col_types = FALSE
) %>%
  rename(FDR = FDR_within_feature_8SCFA)

cox_results <- read_csv(
  file.path(root, "outputs", "tables", "Cox_results_metagenomic_features.csv"),
  show_col_types = FALSE
)
pathway_labels <- read_csv(
  file.path(input_dir, "pathway_labels.csv"),
  show_col_types = FALSE
)

# Apply the historical display rule to the FDR <0.05 rerun coefficients.
heat_data <- coefficients %>%
  group_by(data_layer, feature) %>%
  filter(any(!is.na(beta) & !is.na(FDR) & abs(beta) > 0.10 & FDR < 1)) %>%
  ungroup() %>%
  left_join(
    pathway_labels %>% select(pwy_id, label_publication),
    by = c("feature" = "pwy_id")
  ) %>%
  mutate(
    feature_label = case_when(
      data_layer == "Pathways" & !is.na(label_publication) ~ label_publication,
      data_layer == "Species" ~ gsub("_", " ", gsub("^s__", "", feature)),
      TRUE ~ gsub("_", " ", feature)
    ),
    feature_label = if_else(
      data_layer == "Species",
      paste0(
        toupper(substr(feature_label, 1, 1)),
        substr(feature_label, 2, nchar(feature_label))
      ),
      feature_label
    ),
    # Avoid device-specific font warnings while preserving meaning.
    feature_label = gsub("β", "beta", feature_label, fixed = TRUE),
    SCFA_label = case_when(
      SCFA == "X2.Methylbutyric_acid" ~ "2-Methylbutyric acid",
      TRUE ~ gsub("_", " ", SCFA)
    ),
    stars = case_when(
      FDR < 0.01 ~ "***",
      FDR < 0.05 ~ "**",
      FDR < 0.10 ~ "*",
      TRUE ~ ""
    ),
    data_layer = factor(data_layer, levels = c("Species", "Pathways"))
  )

row_order <- heat_data %>%
  group_by(data_layer, feature, feature_label) %>%
  summarise(max_abs_beta = max(abs(beta), na.rm = TRUE), .groups = "drop") %>%
  arrange(data_layer, desc(max_abs_beta), feature_label)

species_order <- row_order %>%
  filter(data_layer == "Species") %>%
  pull(feature_label)
pathway_order <- row_order %>%
  filter(data_layer == "Pathways") %>%
  pull(feature_label)
row_levels <- rev(c(pathway_order, species_order))

feature_lookup <- heat_data %>%
  distinct(feature, feature_label, data_layer)

scfa_order <- c(
  "Acetic acid",
  "Propionic acid",
  "Butyric acid",
  "Isobutyric acid",
  "Valeric acid",
  "Isovaleric acid",
  "2-Methylbutyric acid",
  "Hexanoic acid"
)
outcome_labels <- c(
  therapy_intensification = "Overall\nintensification",
  medication_intensification = "Medication\nintensification",
  intensive_escalation = "Intensive\nescalation"
)
x_order <- c(scfa_order, " ", unname(outcome_labels))

heat_data <- heat_data %>%
  mutate(
    feature_label = factor(feature_label, levels = row_levels),
    x_label = factor(SCFA_label, levels = x_order)
  )

outcome_data <- cox_results %>%
  filter(feature %in% feature_lookup$feature) %>%
  inner_join(feature_lookup, by = "feature") %>%
  mutate(
    outcome_label = unname(outcome_labels[outcome]),
    x_label = factor(outcome_label, levels = x_order),
    feature_label = factor(feature_label, levels = row_levels),
    nominally_significant = p_value < 0.05,
    log_HR = log(HR),
    absolute_log_HR = abs(log_HR),
    direction = case_when(
      !nominally_significant ~ "Not nominally significant",
      HR > 1 ~ "Higher intensification risk",
      HR < 1 ~ "Lower intensification risk",
      TRUE ~ "Not nominally significant"
    ),
    direction = factor(
      direction,
      levels = c(
        "Higher intensification risk",
        "Lower intensification risk",
        "Not nominally significant"
      )
    )
  )

beta_limit <- max(abs(heat_data$beta), na.rm = TRUE)

p <- ggplot() +
  # Light background cells visually group the three Cox outcome columns.
  geom_tile(
    data = outcome_data,
    aes(x = x_label, y = feature_label),
    fill = "#F4F4F2",
    color = "white",
    linewidth = 0.35
  ) +
  geom_tile(
    data = heat_data,
    aes(x = x_label, y = feature_label, fill = beta),
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    data = heat_data,
    aes(x = x_label, y = feature_label, label = stars),
    color = "#202020",
    size = 3.0,
    fontface = "bold"
  ) +
  # Non-significant Cox associations are small neutral open circles.
  geom_point(
    data = outcome_data %>% filter(!nominally_significant),
    aes(x = x_label, y = feature_label),
    shape = 1,
    color = "#B8B8B8",
    size = 1.5,
    stroke = 0.55
  ) +
  # Significant Cox associations use direction-specific triangles.
  geom_point(
    data = outcome_data %>% filter(nominally_significant),
    aes(
      x = x_label,
      y = feature_label,
      shape = direction,
      color = direction,
      size = absolute_log_HR
    ),
    fill = "white",
    stroke = 1.15
  ) +
  scale_fill_gradient2(
    low = unname(heatmap_palette["low"]),
    mid = unname(heatmap_palette["mid"]),
    high = unname(heatmap_palette["high"]),
    midpoint = 0,
    limits = c(-beta_limit, beta_limit),
    name = expression("Adjusted SCFA " * beta)
  ) +
  scale_shape_manual(
    values = c(
      "Higher intensification risk" = 24,
      "Lower intensification risk" = 25
    ),
    name = "Cox association",
    guide = guide_legend(
      override.aes = list(
        size = 4,
        color = c("#C45A52", "#4B7F52")
      )
    )
  ) +
  scale_color_manual(
    values = c(
      "Higher intensification risk" = "#C45A52",
      "Lower intensification risk" = "#4B7F52"
    ),
    name = "Cox association",
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(2.3, 5.0),
    breaks = c(0.1, 0.3, 0.5),
    name = expression("|log(HR)|")
  ) +
  scale_x_discrete(drop = FALSE) +
  facet_grid(
    data_layer ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  geom_vline(xintercept = 8.5, color = "#666666", linewidth = 0.5) +
  geom_vline(xintercept = 9.5, color = "#666666", linewidth = 0.5) +
  labs(
    x = NULL,
    y = NULL,
    title = "Microbial associations with faecal SCFAs and therapy intensification",
    subtitle = paste0(
      "SCFA tiles show adjusted beta coefficients; stars show within-feature 8-SCFA FDR. ",
      "Cox symbols are shown for nominal p < 0.05."
    ),
    caption = paste0(
      "Upward triangle: HR > 1; downward triangle: HR < 1; ",
      "symbol size: |log(HR)|; grey open circle: p >= 0.05."
    )
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(32, "mm")),
    size = guide_legend(
      order = 3,
      override.aes = list(
        shape = 24,
        color = "#555555",
        fill = "white",
        stroke = 1.0
      )
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(1.5, "mm"),
    panel.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "#2B2B2B"
    ),
    axis.text.y = element_text(face = "bold", size = 7.8, color = "#2B2B2B"),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      color = "#333333",
      margin = margin(r = 7)
    ),
    strip.background = element_rect(fill = "#ECEAE4", color = NA),
    plot.title = element_text(face = "bold", size = 15, color = "#222222"),
    plot.subtitle = element_text(size = 9.5, color = "#4A4A4A"),
    plot.caption = element_text(size = 8.5, color = "#555555", hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 12, 8, 8)
  )

write_csv(
  heat_data,
  file.path(output_dir, "SCFA_heatmap_data.csv")
)
write_csv(
  outcome_data %>%
    select(
      data_layer, feature, feature_label, outcome, outcome_label,
      HR, CI_low, CI_high, p_value, p_adj, nominally_significant,
      log_HR, absolute_log_HR, direction
    ),
  file.path(output_dir, "Cox_outcome_annotation_data.csv")
)

ggsave(
  file.path(output_dir, paste0(output_stem, ".pdf")),
  p,
  width = 12,
  height = 15,
  limitsize = FALSE
)
ggsave(
  file.path(output_dir, paste0(output_stem, ".png")),
  p,
  width = 12,
  height = 15,
  dpi = 320,
  limitsize = FALSE
)

message("Publication-style heatmap written to: ", output_dir)
