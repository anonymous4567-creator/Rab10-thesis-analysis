# ============================================================
# MAKE A TWO-PANEL PORTRAIT pRab10/TOTAL RAB10 RATIO PLOT
#
# When prompted, select:
#   well_level_data_all_35_conditions.csv
#
# This changes only the figure layout. It does not alter the
# data or rerun the statistical analysis.
# ============================================================


# 1. PACKAGES -------------------------------------------------

packages <- c("readr", "dplyr", "stringr", "ggplot2")

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))


# 2. SELECT THE COMPLETED WELL-LEVEL FILE ---------------------

message(
  "Select well_level_data_all_35_conditions.csv from the final_full_35_conditions folder."
)

input_file <- file.choose()

well_data <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "treatment",
  "type",
  "raw_p_rab10_total"
)

missing_columns <- setdiff(required_columns, names(well_data))

if (length(missing_columns) > 0) {
  stop(
    "This is not the required well-level file. Select well_level_data_all_35_conditions.csv."
  )
}

well_data <- well_data |>
  mutate(
    treatment = str_squish(as.character(treatment)),
    type = str_to_lower(str_squish(as.character(type))),
    raw_p_rab10_total = as.numeric(raw_p_rab10_total)
  )


# 3. DIVIDE THE CONDITIONS BETWEEN TWO PANELS ----------------

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

split_point <- ceiling(length(compound_order) / 2)
first_compounds <- compound_order[seq_len(split_point)]
second_compounds <- compound_order[(split_point + 1):length(compound_order)]

panel_a_order <- c(control_order, first_compounds)
panel_b_order <- c("Starved 1 h", second_compounds)

panel_a <- well_data |>
  filter(treatment %in% panel_a_order) |>
  mutate(
    panel = "A. Controls and compounds A–L",
    plot_id = paste0("A___", treatment)
  )

panel_b <- well_data |>
  filter(treatment %in% panel_b_order) |>
  mutate(
    panel = "B. Starved reference and compounds M–V",
    plot_id = paste0("B___", treatment)
  )

plot_data <- bind_rows(panel_a, panel_b)

plot_levels <- c(
  paste0("A___", panel_a_order),
  paste0("B___", panel_b_order)
)

plot_data <- plot_data |>
  mutate(
    panel = factor(
      panel,
      levels = c(
        "A. Controls and compounds A–L",
        "B. Starved reference and compounds M–V"
      )
    ),
    plot_id = factor(plot_id, levels = plot_levels)
  )


# 4. COLOUR PALETTE -------------------------------------------

compound_colours <- setNames(
  grDevices::hcl.colors(
    length(compound_order),
    palette = "Dark 3"
  ),
  compound_order
)

treatment_colours <- c(
  setNames(rep("#252525", length(control_order)), control_order),
  compound_colours
)

# Map the same treatment colour into both panels.
plot_colour_values <- c(
  setNames(
    treatment_colours[panel_a_order],
    paste0("A___", panel_a_order)
  ),
  setNames(
    treatment_colours[panel_b_order],
    paste0("B___", panel_b_order)
  )
)


# 5. CREATE THE TWO-PANEL FIGURE ------------------------------

two_panel_plot <- ggplot(
  plot_data,
  aes(
    x = plot_id,
    y = raw_p_rab10_total,
    colour = plot_id,
    shape = type
  )
) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0),
    size = 3.0,
    alpha = 0.85
  ) +
  stat_summary(
    aes(colour = plot_id),
    fun = mean,
    geom = "crossbar",
    width = 0.48,
    linewidth = 0.65
  ) +
  facet_wrap(
    ~ panel,
    ncol = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(labels) {
      clean_labels <- str_remove(labels, "^[AB]___")
      recode(
        clean_labels,
        "Naïve (no starvation)" = "Naïve",
        "No starvation + Metformin 2 µM" = "Metformin control",
        "Starved 1 h + MLi-2 100 nM" = "Starved + MLi-2",
        "Starved 1 h + BDNF" = "Starved + BDNF",
        .default = clean_labels
      )
    },
    expand = expansion(add = c(0.6, 0.4))
  ) +
  scale_colour_manual(
    values = plot_colour_values,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(compound = 16, control = 17),
    labels = c(compound = "Compound", control = "Control")
  ) +
  scale_y_continuous(
    breaks = seq(0.6, 1.4, 0.2),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  coord_cartesian(ylim = c(0.55, 1.45)) +
  labs(
    x = NULL,
    y = "pRab10/total Rab10 ratio",
    shape = "Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8.4,
      colour = "black"
    ),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.title.y = element_text(
      size = 10.5,
      margin = margin(r = 12)
    ),
    legend.position = "top",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10.5,
      hjust = 0
    ),
    panel.spacing = grid::unit(1.0, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  )


# 6. SAVE BESIDE THE WELL-LEVEL CSV ---------------------------

output_file <- file.path(
  dirname(input_file),
  "07_FINAL_FIXED_two_panel_pRab10_total_Rab10_ratio.png"
)

ggsave(
  filename = output_file,
  plot = two_panel_plot,
  width = 8.27,
  height = 11.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

print(two_panel_plot)

cat(
  "\nTwo-panel figure saved to:\n",
  output_file,
  "\n"
)
