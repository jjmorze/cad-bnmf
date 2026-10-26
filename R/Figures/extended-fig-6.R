# ==============================================================================
#  OUTCOME-SPECIFIC PRS LOLLIPOP / RANGE PLOT
# ==============================================================================
#
# Eligibility:
#   - Include outcome if >=1 process-specific PRS is Bonferroni significant
#   - Bonferroni threshold = 0.05 / 32
#
# For each included outcome:
#   - Total CAD PRS
#   - pPRS with highest OR
#   - pPRS with lowest OR
#
# Visual encoding:
#   - Diseases on y-axis
#   - Odds ratios on log-scaled x-axis
#   - Grey60 line between lowest and highest pPRS estimates
#   - Total CAD PRS = black triangle
#   - pPRS = colored circles
#   - Filled marker = P_Bonferroni < 0.05
#   - Hollow marker = P_Bonferroni >= 0.05
#
# Ordering:
#   - Top to bottom by signed discrepancy:
#
#       log(highest pPRS OR) - log(Total CAD PRS OR)
#
# Panels:
#   - Alzheimer's disease + vascular dementia panel on top
#   - Main panel below
#   - No facet/panel headers
#   - Relative heights approximately 1:10
#
# Legends:
#   - Polygenic Risk Scores:
#       3 columns x 3 rows
#   - Significance:
#       P_Bonferroni >= 0.05
#       P_Bonferroni < 0.05
#
#   Both legend blocks are manually drawn as ordinary ggplots.
#   No cowplot legend extraction is used.
#
# ==============================================================================


# ==============================================================================
#  PACKAGES
# ==============================================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(grid)


# ==============================================================================
#  READ DATA
# ==============================================================================

results <- read.csv(
  "figures/figure-5/logistic_results_ukbb_100626.csv",
  stringsAsFactors = FALSE
)


# ==============================================================================
#  USER SETTINGS
# ==============================================================================


# ------------------------------------------------------------------------------
# Bonferroni significance
# ------------------------------------------------------------------------------

bonferroni_denominator <- 32

bonferroni_threshold <- 0.05 / bonferroni_denominator


message(
  "Bonferroni raw-P threshold: ",
  signif(
    bonferroni_threshold,
    6
  )
)


# ------------------------------------------------------------------------------
# Relative panel heights
# ------------------------------------------------------------------------------

neurological_panel_height <- 1

main_panel_height <- 10


# ------------------------------------------------------------------------------
# Typography
# ------------------------------------------------------------------------------

base_text_size_pt <- 11

legend_text_size_pt <- 10.5

legend_title_size_pt <- 11


# ------------------------------------------------------------------------------
# Point and line geometry
# ------------------------------------------------------------------------------

point_size <- 3.2

point_stroke <- 0.80

range_linewidth <- 2.2

reference_linewidth <- 0.45

gridline_linewidth <- 0.25

panel_border_linewidth <- 0.4


# ------------------------------------------------------------------------------
# Colors
# ------------------------------------------------------------------------------

range_line_colour <- "grey60"

gridline_colour <- "grey70"


# ------------------------------------------------------------------------------
# Figure dimensions
# ------------------------------------------------------------------------------

figure_width <- 9

figure_height <- 10

figure_dpi <- 1200


# ==============================================================================
#  PRS COLORS
# ==============================================================================

cluster_cols <- c(
  "ChEL" = "#B22D4C",
  "TRL" = "#E07A8E",
  "Lp(a)" = "#7A1F3A",
  "MUNW" = "#4A7C4E",
  "Obesity" = "#95C075",
  "Endothelial" = "#2D5F8B",
  "Coagulation" = "#7AAED1",
  "Cardiac remodeling" = "#5D4787"
)


all_score_cols <- c(
  "Total CAD PRS" = "black",
  cluster_cols
)


# ==============================================================================
#  PRS SHAPES
# ==============================================================================

# Shape 24 = filled triangle
# Shape 21 = filled circle

all_score_shapes <- c(
  "Total CAD PRS" = 24,
  "ChEL" = 21,
  "TRL" = 21,
  "Lp(a)" = 21,
  "MUNW" = 21,
  "Obesity" = 21,
  "Endothelial" = 21,
  "Coagulation" = 21,
  "Cardiac remodeling" = 21
)


# ==============================================================================
#  MAP INTERNAL PRS VARIABLES TO DISPLAY LABELS
# ==============================================================================

prs_label_map <- c(
  total_196_prs_std   = "Total CAD PRS",
  chol_prs_std        = "ChEL",
  trl_prs_std         = "TRL",
  lpa_prs_std         = "Lp(a)",
  mdwo_prs_std        = "MUNW",
  obesity_prs_std     = "Obesity",
  endothelial_prs_std = "Endothelial",
  coagulation_prs_std = "Coagulation",
  cardiac_prs_std     = "Cardiac remodeling"
)


total_prs_variable <- "total_196_prs_std"


pprs_variables <- setdiff(
  names(prs_label_map),
  total_prs_variable
)


# ==============================================================================
#  CHECK REQUIRED COLUMNS
# ==============================================================================

required_columns <- c(
  "outcome_short",
  "outcome_long",
  "predictor_short",
  "or",
  "lb",
  "ub",
  "p_value"
)


missing_columns <- setdiff(
  required_columns,
  names(results)
)


if (length(missing_columns) > 0) {
  
  stop(
    "The following required columns are missing: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ==============================================================================
#  PREPARE DATA
# ==============================================================================

results <- results %>%
  filter(
    predictor_short %in% names(prs_label_map)
  ) %>%
  mutate(
    
    prs = unname(
      prs_label_map[predictor_short]
    ),
    
    
    # Raw-P equivalent of Bonferroni-adjusted P < 0.05
    
    significant = (
      !is.na(p_value) &
        p_value < bonferroni_threshold
    ),
    
    
    # Explicit Bonferroni-adjusted P value
    
    p_bonferroni_plot = pmin(
      p_value * bonferroni_denominator,
      1
    )
  )


# ==============================================================================
#  IDENTIFY ELIGIBLE OUTCOMES
# ==============================================================================

# Outcome is included only if >=1 process-specific PRS satisfies:
#
#   P < 0.05 / 32
#
# Total CAD PRS significance alone does not qualify the outcome.

significant_outcomes <- results %>%
  filter(
    predictor_short %in% pprs_variables
  ) %>%
  group_by(
    outcome_short,
    outcome_long
  ) %>%
  summarise(
    
    any_significant_pprs = any(
      significant,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  filter(
    any_significant_pprs
  )


message(
  "Number of outcomes with >=1 Bonferroni-significant pPRS: ",
  nrow(significant_outcomes)
)


# ==============================================================================
#  PROCESS-SPECIFIC PRS RESULTS
# ==============================================================================

pprs_results <- results %>%
  semi_join(
    significant_outcomes,
    by = c(
      "outcome_short",
      "outcome_long"
    )
  ) %>%
  filter(
    predictor_short %in% pprs_variables
  )


# ==============================================================================
#  HIGHEST pPRS
# ==============================================================================

highest_pprs <- pprs_results %>%
  group_by(
    outcome_short,
    outcome_long
  ) %>%
  arrange(
    desc(or),
    predictor_short,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  mutate(
    selection = "Highest pPRS"
  )


# ==============================================================================
#  LOWEST pPRS
# ==============================================================================

lowest_pprs <- pprs_results %>%
  group_by(
    outcome_short,
    outcome_long
  ) %>%
  arrange(
    or,
    predictor_short,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  mutate(
    selection = "Lowest pPRS"
  )


# ==============================================================================
#  TOTAL CAD PRS
# ==============================================================================

total_prs_results <- results %>%
  semi_join(
    significant_outcomes,
    by = c(
      "outcome_short",
      "outcome_long"
    )
  ) %>%
  filter(
    predictor_short == total_prs_variable
  ) %>%
  mutate(
    selection = "Total CAD PRS"
  )


# ==============================================================================
#  COMBINE DISPLAYED ESTIMATES
# ==============================================================================

plot_data <- bind_rows(
  total_prs_results,
  highest_pprs,
  lowest_pprs
)


# ==============================================================================
#  CHECK THREE ESTIMATES PER OUTCOME
# ==============================================================================

invalid_outcomes <- plot_data %>%
  count(
    outcome_short,
    outcome_long
  ) %>%
  filter(
    n != 3
  )


if (nrow(invalid_outcomes) > 0) {
  
  warning(
    "Some outcomes do not contain exactly three selected estimates."
  )
  
  print(
    invalid_outcomes
  )
}


# ==============================================================================
#  GREY pPRS RANGE
# ==============================================================================

# Grey line spans:
#
#   lowest pPRS OR -------- highest pPRS OR
#
# It does NOT represent a confidence interval.

range_data <- highest_pprs %>%
  select(
    outcome_short,
    outcome_long,
    high_prs = prs,
    high_or = or
  ) %>%
  left_join(
    
    lowest_pprs %>%
      select(
        outcome_short,
        outcome_long,
        low_prs = prs,
        low_or = or
      ),
    
    by = c(
      "outcome_short",
      "outcome_long"
    )
  )


# ==============================================================================
#  CALCULATE SIGNED DISCREPANCY FOR ORDERING
# ==============================================================================

# Signed discrepancy:
#
#   log(highest pPRS OR) - log(Total CAD PRS OR)
#
# Largest positive value appears at the top.

outcome_order_data <- highest_pprs %>%
  select(
    outcome_short,
    outcome_long,
    highest_pprs_or = or
  ) %>%
  left_join(
    
    total_prs_results %>%
      select(
        outcome_short,
        outcome_long,
        total_prs_or = or
      ),
    
    by = c(
      "outcome_short",
      "outcome_long"
    )
  ) %>%
  mutate(
    
    signed_discrepancy =
      log(highest_pprs_or) -
      log(total_prs_or)
  )


# ==============================================================================
#  DEFINE NEUROLOGICAL PANEL
# ==============================================================================

special_outcomes <- c(
  "ad_ukb_ind",
  "vd_ukb_ind"
)


outcome_order_data <- outcome_order_data %>%
  mutate(
    
    panel = if_else(
      outcome_short %in% special_outcomes,
      "neurological",
      "main"
    )
  )


# ==============================================================================
#  ORDER OUTCOMES WITHIN PANELS
# ==============================================================================

# ggplot discrete y levels are drawn bottom-to-top.
#
# Therefore ascending factor levels result in the largest discrepancy
# being plotted at the top.

main_order <- outcome_order_data %>%
  filter(
    panel == "main"
  ) %>%
  arrange(
    signed_discrepancy
  ) %>%
  pull(
    outcome_long
  )


neurological_order <- outcome_order_data %>%
  filter(
    panel == "neurological"
  ) %>%
  arrange(
    signed_discrepancy
  ) %>%
  pull(
    outcome_long
  )


# ==============================================================================
#  PREPARE MAIN PANEL DATA
# ==============================================================================

main_plot_data <- plot_data %>%
  filter(
    !outcome_short %in% special_outcomes
  ) %>%
  left_join(
    
    outcome_order_data %>%
      select(
        outcome_short,
        signed_discrepancy
      ),
    
    by = "outcome_short"
  ) %>%
  mutate(
    
    outcome_long = factor(
      outcome_long,
      levels = main_order
    ),
    
    
    point_fill = if_else(
      significant,
      unname(
        all_score_cols[prs]
      ),
      "white"
    )
  )


main_range_data <- range_data %>%
  filter(
    !outcome_short %in% special_outcomes
  ) %>%
  mutate(
    
    outcome_long = factor(
      outcome_long,
      levels = main_order
    )
  )


# ==============================================================================
#  PREPARE NEUROLOGICAL PANEL DATA
# ==============================================================================

neurological_plot_data <- plot_data %>%
  filter(
    outcome_short %in% special_outcomes
  ) %>%
  left_join(
    
    outcome_order_data %>%
      select(
        outcome_short,
        signed_discrepancy
      ),
    
    by = "outcome_short"
  ) %>%
  mutate(
    
    outcome_long = factor(
      outcome_long,
      levels = neurological_order
    ),
    
    
    point_fill = if_else(
      significant,
      unname(
        all_score_cols[prs]
      ),
      "white"
    )
  )


neurological_range_data <- range_data %>%
  filter(
    outcome_short %in% special_outcomes
  ) %>%
  mutate(
    
    outcome_long = factor(
      outcome_long,
      levels = neurological_order
    )
  )


# ==============================================================================
#  COMMON DATA-PANEL FUNCTION
# ==============================================================================

make_lollipop_panel <- function(
    point_data,
    range_data,
    show_x_title = TRUE
) {
  
  ggplot() +
    
    
    # --------------------------------------------------------------------------
  # Grey pPRS range
  # --------------------------------------------------------------------------
  
  geom_segment(
    data = range_data,
    
    aes(
      x = low_or,
      xend = high_or,
      y = outcome_long,
      yend = outcome_long
    ),
    
    colour = range_line_colour,
    linewidth = range_linewidth,
    lineend = "round"
  ) +
    
    
    # --------------------------------------------------------------------------
  # OR = 1 reference
  # --------------------------------------------------------------------------
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = reference_linewidth,
    colour = "grey35"
  ) +
    
    
    # --------------------------------------------------------------------------
  # Selected PRS estimates
  # --------------------------------------------------------------------------
  
  geom_point(
    data = point_data,
    
    aes(
      x = or,
      y = outcome_long,
      colour = prs,
      shape = prs
    ),
    
    fill = point_data$point_fill,
    
    size = point_size,
    stroke = point_stroke,
    
    show.legend = FALSE
  ) +
    
    
    # --------------------------------------------------------------------------
  # Colors
  # --------------------------------------------------------------------------
  
  scale_colour_manual(
    values = all_score_cols,
    limits = names(all_score_cols),
    drop = FALSE
  ) +
    
    
    # --------------------------------------------------------------------------
  # Shapes
  # --------------------------------------------------------------------------
  
  scale_shape_manual(
    values = all_score_shapes,
    limits = names(all_score_shapes),
    drop = FALSE
  ) +
    
    
    # --------------------------------------------------------------------------
  # Log OR axis
  # --------------------------------------------------------------------------
  
  scale_x_log10(
    
    breaks = c(
      0.5,
      0.75,
      1,
      1.5,
      2,
      3,
      5,
      10
    ),
    
    labels = c(
      "0.5",
      "0.75",
      "1.0",
      "1.5",
      "2.0",
      "3.0",
      "5.0",
      "10.0"
    ),
    
    expand = expansion(
      mult = c(
        0.04,
        0.06
      )
    )
  ) +
    
    
    labs(
      
      x = if (
        show_x_title
      ) {
        "Odds ratio for upper vs. bottom decile"
      } else {
        NULL
      },
      
      y = NULL
    ) +
    
    
    theme_linedraw(
      base_family = "Arial",
      base_size = base_text_size_pt
    ) +
    
    
    theme(
      
      text = element_text(
        family = "Arial",
        size = base_text_size_pt,
        colour = "black"
      ),
      
      
      axis.title.x = element_text(
        family = "Arial",
        size = base_text_size_pt,
        colour = "black",
        
        margin = margin(
          t = 7
        )
      ),
      
      
      axis.text = element_text(
        family = "Arial",
        size = base_text_size_pt,
        colour = "black"
      ),
      
      
      axis.text.y = element_text(
        family = "Arial",
        size = base_text_size_pt,
        colour = "black",
        
        margin = margin(
          r = 5
        )
      ),
      
      
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.3
      ),
      
      
      axis.ticks.length = grid::unit(
        2,
        "pt"
      ),
      
      
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = panel_border_linewidth
      ),
      
      
      panel.grid.major.y = element_blank(),
      
      panel.grid.minor = element_blank(),
      
      
      panel.grid.major.x = element_line(
        colour = gridline_colour,
        linewidth = gridline_linewidth
      ),
      
      
      legend.position = "none",
      
      
      plot.margin = margin(
        t = 3,
        r = 5,
        b = 3,
        l = 4
      )
    )
}


# ==============================================================================
#  BUILD NEUROLOGICAL PANEL
# ==============================================================================

neurological_plot <- make_lollipop_panel(
  point_data = neurological_plot_data,
  range_data = neurological_range_data,
  show_x_title = FALSE
)


# ==============================================================================
#  BUILD MAIN PANEL
# ==============================================================================

main_plot <- make_lollipop_panel(
  point_data = main_plot_data,
  range_data = main_range_data,
  show_x_title = TRUE
)


# ==============================================================================
#  CUSTOM POLYGENIC RISK SCORE LEGEND DATA
# ==============================================================================

# Manual 3 x 3 grid:
#
# Row 1:
#   Total CAD PRS | ChEL | TRL
#
# Row 2:
#   Lp(a) | MUNW | Obesity
#
# Row 3:
#   Endothelial | Coagulation | Cardiac remodeling

prs_legend_data <- tibble(
  
  prs = c(
    "Total CAD PRS",
    "ChEL",
    "TRL",
    "Lp(a)",
    "MUNW",
    "Obesity",
    "Endothelial",
    "Coagulation",
    "Cardiac remodeling"
  ),
  
  col = c(
    1, 2, 3,
    1, 2, 3,
    1, 2, 3
  ),
  
  row = c(
    3, 3, 3,
    2, 2, 2,
    1, 1, 1
  )
) %>%
  mutate(
    
    colour = unname(
      all_score_cols[prs]
    ),
    
    shape = unname(
      all_score_shapes[prs]
    )
  )


# ==============================================================================
#  CUSTOM POLYGENIC RISK SCORE LEGEND
# ==============================================================================

prs_legend_plot <- ggplot(
  prs_legend_data
) +
  
  
  geom_point(
    aes(
      x = col,
      y = row,
      colour = prs,
      fill = prs,
      shape = prs
    ),
    
    size = point_size,
    stroke = point_stroke
  ) +
  
  
  geom_text(
    aes(
      x = col + 0.12,
      y = row,
      label = prs
    ),
    
    family = "Arial",
    size = legend_text_size_pt / ggplot2::.pt,
    colour = "black",
    hjust = 0,
    vjust = 0.5
  ) +
  
  
  scale_colour_manual(
    values = all_score_cols,
    guide = "none"
  ) +
  
  
  scale_fill_manual(
    values = all_score_cols,
    guide = "none"
  ) +
  
  
  scale_shape_manual(
    values = all_score_shapes,
    guide = "none"
  ) +
  
  
  annotate(
    "text",
    
    x = 2,
    y = 4.05,
    
    label = "Polygenic Risk Scores",
    
    family = "Arial",
    fontface = "bold",
    size = legend_title_size_pt / ggplot2::.pt,
    
    hjust = 0.5
  ) +
  
  
  scale_x_continuous(
    limits = c(
      0.8,
      4.0
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  
  scale_y_continuous(
    limits = c(
      0.6,
      4.3
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  
  coord_cartesian(
    clip = "off"
  ) +
  
  
  theme_void(
    base_family = "Arial"
  ) +
  
  
  theme(
    
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )


# ==============================================================================
#  CUSTOM SIGNIFICANCE LEGEND DATA
# ==============================================================================

significance_legend_data <- tibble(
  
  x = c(
    1,
    1
  ),
  
  y = c(
    2,
    1
  ),
  
  fill_value = c(
    "white",
    "black"
  ),
  
  label = c(
    "P[Bonferroni] >= 0.05",
    "P[Bonferroni] < 0.05"
  )
)


# ==============================================================================
#  CUSTOM SIGNIFICANCE LEGEND
# ==============================================================================

significance_legend_plot <- ggplot(
  significance_legend_data
) +
  
  
  # Hollow / filled circles
  
  geom_point(
    aes(
      x = x,
      y = y
    ),
    
    shape = 21,
    
    fill = significance_legend_data$fill_value,
    
    colour = "black",
    
    size = point_size,
    
    stroke = point_stroke
  ) +
  
  
  # Mathematical labels:
  #
  # P_Bonferroni >= 0.05
  # P_Bonferroni < 0.05
  
  geom_text(
    aes(
      x = x + 0.16,
      y = y,
      label = label
    ),
    
    parse = TRUE,
    
    family = "Arial",
    
    size = legend_text_size_pt / ggplot2::.pt,
    
    colour = "black",
    
    hjust = 0,
    
    vjust = 0.5
  ) +
  
  
  annotate(
    "text",
    
    x = 1.65,
    y = 3.0,
    
    label = "Significance",
    
    family = "Arial",
    
    fontface = "bold",
    
    size = legend_title_size_pt / ggplot2::.pt,
    
    hjust = 0.5
  ) +
  
  
  scale_x_continuous(
    limits = c(
      0.85,
      2.45
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  
  scale_y_continuous(
    limits = c(
      0.55,
      3.25
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  
  coord_cartesian(
    clip = "off"
  ) +
  
  
  theme_void(
    base_family = "Arial"
  ) +
  
  
  theme(
    
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )


# ==============================================================================
#  COMBINE CUSTOM LEGENDS
# ==============================================================================

# Both are ordinary ggplot objects.
#
# No get_legend(), no cowplot, and therefore no internal PDF/null device.

legend_core <- patchwork::wrap_plots(
  
  prs_legend_plot,
  significance_legend_plot,
  
  nrow = 1,
  
  widths = c(
    3.6,
    1.5
  )
)


# ==============================================================================
#  CENTER LEGEND BAND
# ==============================================================================

# Add empty spacers at left and right.

legend_band <- patchwork::wrap_plots(
  
  patchwork::plot_spacer(),
  legend_core,
  patchwork::plot_spacer(),
  
  nrow = 1,
  
  widths = c(
    0.25,
    5,
    0.25
  )
)


# ==============================================================================
#  COMBINE DATA PANELS
# ==============================================================================

data_panels <- patchwork::wrap_plots(
  
  neurological_plot,
  main_plot,
  
  ncol = 1,
  
  heights = c(
    neurological_panel_height,
    main_panel_height
  )
)


# ==============================================================================
#  LOCK DATA PANELS
# ==============================================================================

data_panels_locked <- patchwork::wrap_elements(
  full = data_panels
)


# ==============================================================================
#  LOCK LEGEND BAND
# ==============================================================================

legend_band_locked <- patchwork::wrap_elements(
  full = legend_band
)


# ==============================================================================
#  FINAL FIGURE
# ==============================================================================

combined_plot <- patchwork::wrap_plots(
  
  data_panels_locked,
  legend_band_locked,
  
  ncol = 1,
  
  heights = c(
    1,
    0.13
  )
)


# ==============================================================================
#  DISPLAY
# ==============================================================================

print(
  combined_plot
)


# ==============================================================================
#  SAVE PNG
# ==============================================================================

ggsave(
  filename = "figures/ext-data-figure-6/extended-fig-6.png",
  
  plot = combined_plot,
  
  width = figure_width,
  
  height = figure_height,
  
  units = "in",
  
  dpi = figure_dpi,
  
  bg = "white"
)



# ==============================================================================
#  END
# ==============================================================================