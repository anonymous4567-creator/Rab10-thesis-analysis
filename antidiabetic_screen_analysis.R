# ============================================================
# FINAL ANTIDIABETIC SCREEN ANALYSIS: ALL 35 CONDITIONS — RUN THIS VERSION
#
# Input file:
#   Selected through a file-selection window when the script runs
#
# Statistical unit:
#   well (image fields are averaged within each well)
#
# Included conditions:
#   all 30 compounds and all 5 controls
#
# Dunnett reference:
#   Starved 1 h
# ============================================================


# 1. PACKAGES -------------------------------------------------

packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tidyr",
  "purrr",
  "ggplot2",
  "emmeans",
  "forcats",
  "janitor",
  "writexl"
)

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))


# 2. INPUT AND OUTPUT PATHS -----------------------------------

message("Select the complete image-level data1 CSV file.")
input_file <- file.choose()

output_folder <- file.path(
  dirname(input_file),
  "final_full_35_conditions"
)
figure_folder <- file.path(output_folder, "figures")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)


# 3. IMPORT IMAGE-LEVEL DATA ----------------------------------

image_data <- read_csv(
  input_file,
  locale = locale(encoding = "Windows-1252"),
  show_col_types = FALSE
) |>
  clean_names()

required_columns <- c(
  "plate",
  "well",
  "field",
  "screen_id",
  "treatment",
  "type",
  "p_rab10_int_den",
  "total_rab10_int_den",
  "map2_int_den"
)

missing_columns <- setdiff(required_columns, names(image_data))

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}


# 4. CLEAN VALUES AND CALCULATE IMAGE-LEVEL RATIOS ------------

image_data <- image_data |>
  mutate(
    plate = str_to_upper(str_squish(as.character(plate))),
    well = str_to_upper(str_squish(as.character(well))),
    screen_id = str_squish(as.character(screen_id)),
    treatment = str_squish(
      str_replace_all(as.character(treatment), "\u00A0", " ")
    ),
    type = str_to_lower(str_squish(as.character(type))),
    field = as.integer(field),
    across(
      c(p_rab10_int_den, total_rab10_int_den, map2_int_den),
      as.numeric
    ),
    p_rab10_total_ratio = if_else(
      total_rab10_int_den > 0,
      p_rab10_int_den / total_rab10_int_den,
      NA_real_
    )
  )


# 5. AVERAGE IMAGE FIELDS WITHIN EACH WELL --------------------

well_data <- image_data |>
  group_by(plate, well, screen_id, treatment, type) |>
  summarise(
    n_fields = n_distinct(field),
    raw_p_rab10 = mean(p_rab10_int_den, na.rm = TRUE),
    raw_total_rab10 = mean(total_rab10_int_den, na.rm = TRUE),
    raw_map2 = mean(map2_int_den, na.rm = TRUE),
    raw_p_rab10_total = mean(p_rab10_total_ratio, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  well_data,
  file.path(output_folder, "well_level_data_all_35_conditions.csv")
)


# 6. VERIFY THE INCLUDED CONDITIONS ---------------------------

condition_check <- well_data |>
  group_by(treatment, type) |>
  summarise(
    n_wells = n(),
    plates = paste(sort(unique(plate)), collapse = ", "),
    .groups = "drop"
  )

write_csv(
  condition_check,
  file.path(output_folder, "included_conditions.csv")
)

if (n_distinct(well_data$treatment) != 35) {
  warning(
    "The input contains ",
    n_distinct(well_data$treatment),
    " conditions rather than the expected 35. Check included_conditions.csv."
  )
}

if (!"Starved 1 h" %in% well_data$treatment) {
  stop('The Dunnett reference "Starved 1 h" was not found.')
}


# 7. DESCRIPTIVE STATISTICS -----------------------------------

outcomes <- c(
  "raw_p_rab10",
  "raw_total_rab10",
  "raw_map2",
  "raw_p_rab10_total"
)

descriptive_statistics <- well_data |>
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "outcome",
    values_to = "value"
  ) |>
  group_by(treatment, type, outcome) |>
  summarise(
    n_wells = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    sem = sd / sqrt(n_wells),
    minimum = min(value, na.rm = TRUE),
    maximum = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  descriptive_statistics,
  file.path(output_folder, "descriptive_statistics_all_35_conditions.csv")
)


# 8. DUNNETT TEST AGAINST STARVED 1 h -------------------------

run_dunnett <- function(data, outcome_name) {

  model_data <- data |>
    select(treatment, type, all_of(outcome_name)) |>
    rename(value = all_of(outcome_name)) |>
    filter(!is.na(value))

  level_order <- c(
    "Starved 1 h",
    sort(setdiff(unique(model_data$treatment), "Starved 1 h"))
  )

  model_data <- model_data |>
    mutate(treatment = factor(treatment, levels = level_order))

  model <- lm(value ~ treatment, data = model_data)
  estimated_means <- emmeans(model, ~ treatment)

  results <- contrast(
    estimated_means,
    method = "trt.vs.ctrl",
    ref = 1,
    adjust = "dunnett"
  ) |>
    summary(infer = c(TRUE, TRUE)) |>
    as.data.frame() |>
    as_tibble()

  comparison_levels <- level_order[-1]

  if (nrow(results) != length(comparison_levels)) {
    stop("The Dunnett result rows do not match the treatment levels.")
  }

  results |>
    mutate(
      outcome = outcome_name,
      treatment = comparison_levels,
      significance = case_when(
        p.value < 0.0001 ~ "****",
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    select(
      outcome,
      treatment,
      estimate,
      SE,
      df,
      lower.CL,
      upper.CL,
      t.ratio,
      p.value,
      significance
    )
}

dunnett_results <- run_dunnett(
  well_data,
  "raw_p_rab10_total"
) |>
  left_join(
    condition_check |>
      select(treatment, type, n_wells),
    by = "treatment"
  ) |>
  arrange(p.value)

write_csv(
  dunnett_results,
  file.path(output_folder, "dunnett_pRab10_total_all_35_conditions.csv")
)


# 9. ORDER AND COLOURS FOR FIGURES ----------------------------

control_order <- c(
  "Naïve (no starvation)",
  "No starvation + Metformin 2 µM",
  "Starved 1 h",
  "Starved 1 h + MLi-2 100 nM",
  "Starved 1 h + BDNF"
)

control_order <- control_order[control_order %in% well_data$treatment]

compound_order <- well_data |>
  filter(type == "compound") |>
  distinct(treatment) |>
  arrange(treatment) |>
  pull(treatment) |>
  as.character()

treatment_order <- c(control_order, compound_order)

well_data <- well_data |>
  mutate(treatment = factor(treatment, levels = treatment_order))

compound_colours <- setNames(
  grDevices::hcl.colors(
    length(compound_order),
    palette = "Dark 3"
  ),
  compound_order
)

colour_values <- c(
  setNames(rep("black", length(control_order)), control_order),
  compound_colours
)


# 10. FIGURE FUNCTION -----------------------------------------

make_plot <- function(outcome_name, y_label, file_name) {

  plot_data <- well_data |>
    select(treatment, type, all_of(outcome_name)) |>
    rename(value = all_of(outcome_name)) |>
    filter(!is.na(value))

  plot_object <- ggplot(
    plot_data,
    aes(
      x = treatment,
      y = value,
      colour = treatment,
      shape = type
    )
  ) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0),
      size = 2.4,
      alpha = 0.8
    ) +
    stat_summary(
      aes(colour = treatment),
      fun = mean,
      geom = "crossbar",
      width = 0.45,
      linewidth = 0.45
    ) +
    scale_colour_manual(values = colour_values, guide = "none") +
    labs(
      x = NULL,
      y = y_label,
      subtitle = paste(
        "Each point represents one well; crossbar = mean.",
        "Compound measurements are from Plates A and B; controls are from Plate A."
      ),
      shape = "Type"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      legend.position = "top"
    )

  ggsave(
    filename = file.path(figure_folder, file_name),
    plot = plot_object,
    width = 16,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  plot_object
}


# 11. CREATE FINAL RAW AND RATIO FIGURES ----------------------

make_plot(
  "raw_p_rab10",
  "Raw pRab10 fluorescence intensity (a.u.)",
  "01_raw_pRab10_all_35_conditions.png"
)

make_plot(
  "raw_total_rab10",
  "Raw total Rab10 fluorescence intensity (a.u.)",
  "02_raw_total_Rab10_all_35_conditions.png"
)

make_plot(
  "raw_map2",
  "Raw MAP2 fluorescence intensity (a.u.)",
  "03_raw_MAP2_all_35_conditions.png"
)

make_plot(
  "raw_p_rab10_total",
  "pRab10/total Rab10 ratio",
  "04_pRab10_total_Rab10_ratio_all_35_conditions.png"
)


# 12. COMBINED EXCEL OUTPUT -----------------------------------

write_xlsx(
  list(
    included_conditions = condition_check,
    well_level_data = well_data,
    descriptive_statistics = descriptive_statistics,
    dunnett_results = dunnett_results
  ),
  file.path(output_folder, "analysis_all_35_conditions.xlsx")
)


# 13. FINISH --------------------------------------------------

cat(
  "\nAnalysis completed.\n",
  "Conditions included:", n_distinct(well_data$treatment), "\n",
  "Well-level observations:", nrow(well_data), "\n",
  "Dunnett reference: Starved 1 h\n",
  "Outputs:", output_folder, "\n"
)
