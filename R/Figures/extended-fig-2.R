# ==============================================================================
#  EXTENDED DATA FIGURE 3
# ==============================================================================
#
# Panel a:
#   Sex-stratified PRS associations among White participants
#   Facets: Male, Female
#
# Panel b:
#   Ethnicity-stratified PRS associations
#   Facets: White, Asian, Black, Mixed/Other
#
# Layout:
#   - Two complete panels in one row
#   - Adjustable relative widths for panels a and b
#   - Adjustable spacing between subgroup facets
#   - Adjustable spacing between forest plots and P-interaction columns
#   - Fixed-width P-interaction columns
#
# Estimates:
#   Odds ratios comparing approximately the 95th versus 5th percentile
#   of each standardized PRS
#
# ==============================================================================

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Desktop – Jakub’s Mac mini/Projects/UK Biobank/PPRS/bnmf-cad-clusters/bnmf-cad-clusters/figures/ext-data-figure-3")

# ==============================================================================
#  PACKAGES
# ==============================================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(grid)


# ==============================================================================
#  USER SETTINGS
# ==============================================================================

# ------------------------------------------------------------------------------
# Typography
# ------------------------------------------------------------------------------

base_text_size_pt <- 11
facet_text_size_pt <- 12
panel_label_size_pt <- 14
p_value_text_size_pt <- 10.5


# ------------------------------------------------------------------------------
# Alternating row fills in forest plots only
# ------------------------------------------------------------------------------

stripe_grey <- "grey94"
stripe_white <- "white"


# ------------------------------------------------------------------------------
# Forest-plot geometry
# ------------------------------------------------------------------------------

point_size <- 2.7
ci_linewidth <- 0.65
reference_linewidth <- 0.45
panel_border_linewidth <- 0.4


# ------------------------------------------------------------------------------
# Panel a x-axis: sex
# ------------------------------------------------------------------------------

sex_x_limits <- c(
  0.90,
  2.00
)

sex_x_breaks <- c(
  1.00,
  1.20,
  1.50,
  2.00
)

sex_x_labels <- c(
  "1.0",
  "1.2",
  "1.5",
  "2.0"
)


# ------------------------------------------------------------------------------
# Panel b x-axis: ethnicity
# ------------------------------------------------------------------------------

race_x_limits <- c(
  0.60,
  3.50
)

race_x_breaks <- c(
  1.0,
  2.0,
  3.0
)

race_x_labels <- c(
  "1.0",
  "2.0",
  "3.0"
)


# ------------------------------------------------------------------------------
# Width and spacing controls
# ------------------------------------------------------------------------------

# Relative widths of the two complete panels
sex_panel_width <- 1.00
race_panel_width <- 1.75

# Space between subgroup facets inside each forest plot
sex_subgroup_gap_cm <- 0.30
race_subgroup_gap_cm <- 0.30

# Space between the forest facets and P-interaction column
forest_p_gap_cm <- 0.12

# Fixed physical width of each P-interaction column
p_interaction_width_cm <- 1.75

# Figure dimensions
figure_width <- 18
figure_height <- 5.8
figure_dpi <- 1200


# ==============================================================================
#  PRS ORDER AND COLORS
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

prs_order <- names(cluster_cols)

# Internal order from bottom to top
prs_y_order <- rev(prs_order)

n_prs <- length(prs_order)


# ==============================================================================
#  MAP PREDICTOR VARIABLES TO DISPLAY LABELS
# ==============================================================================

prs_label_map <- c(
  chol_prs_std        = "ChEL",
  trl_prs_std         = "TRL",
  lpa_prs_std         = "Lp(a)",
  mdwo_prs_std        = "MUNW",
  obesity_prs_std     = "Obesity",
  endothelial_prs_std = "Endothelial",
  coagulation_prs_std = "Coagulation",
  cardiac_prs_std     = "Cardiac remodeling"
)


# ==============================================================================
#  CHECK REQUIRED INPUT OBJECTS
# ==============================================================================

required_objects <- c(
  "sex_interaction_white",
  "race_interaction_whole_cohort"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "The following result objects are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}


required_result_columns <- c(
  "predictor_short",
  "subgroup_level",
  "or",
  "lb",
  "ub",
  "p_interaction"
)

missing_sex_columns <- setdiff(
  required_result_columns,
  names(sex_interaction_white)
)

missing_race_columns <- setdiff(
  required_result_columns,
  names(race_interaction_whole_cohort)
)

if (length(missing_sex_columns) > 0) {
  stop(
    "sex_interaction_white is missing columns: ",
    paste(missing_sex_columns, collapse = ", ")
  )
}

if (length(missing_race_columns) > 0) {
  stop(
    "race_interaction_whole_cohort is missing columns: ",
    paste(missing_race_columns, collapse = ", ")
  )
}


# ==============================================================================
#  HELPER: FORMAT P-VALUES IN SCIENTIFIC NOTATION
# ==============================================================================

format_p_scientific <- function(p) {
  ifelse(
    is.na(p),
    "",
    formatC(
      p,
      format = "e",
      digits = 2
    )
  )
}


# ==============================================================================
#  HELPER: ADD NUMERIC PRS ROW POSITIONS
# ==============================================================================

add_prs_positions <- function(data) {
  data %>%
    mutate(
      prs = unname(
        prs_label_map[predictor_short]
      ),
      
      # y = 1 is the lowest row and y = n_prs is the highest row
      y = match(
        prs,
        prs_y_order
      )
    )
}


# ==============================================================================
#  ALTERNATING STRIPE DATA
# ==============================================================================

stripe_data <- tibble(
  y = seq_len(n_prs),
  
  ymin = y - 0.5,
  ymax = y + 0.5,
  
  stripe_fill = if_else(
    y %% 2L == 1L,
    stripe_grey,
    stripe_white
  )
)


# ==============================================================================
#  SHARED FOREST-PLOT THEME
# ==============================================================================

theme_extended_fig3 <- function(
    subgroup_gap_cm = 0.20,
    plot_margin = margin(
      t = 4,
      r = 0,
      b = 2,
      l = 4
    )
) {
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
      
      axis.title = element_text(
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
      
      strip.background = element_rect(
        fill = "black",
        colour = "black",
        linewidth = panel_border_linewidth
      ),
      
      strip.text = element_text(
        family = "Arial",
        face = "bold",
        size = facet_text_size_pt,
        colour = "white",
        margin = margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
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
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.spacing.x = grid::unit(
        subgroup_gap_cm,
        "cm"
      ),
      
      legend.position = "none",
      
      plot.title = element_blank(),
      
      plot.tag = element_text(
        family = "Arial",
        face = "bold",
        size = panel_label_size_pt,
        colour = "black"
      ),
      
      plot.tag.position = c(
        0,
        1
      ),
      
      plot.margin = plot_margin
    )
}


# ==============================================================================
#  PREPARE SEX-INTERACTION RESULTS
# ==============================================================================

sex_plot_data <- sex_interaction_white %>%
  filter(
    predictor_short %in% names(prs_label_map),
    subgroup_level %in% c(
      "Male",
      "Female"
    )
  ) %>%
  add_prs_positions() %>%
  mutate(
    subgroup = factor(
      subgroup_level,
      levels = c(
        "Male",
        "Female"
      )
    )
  )


# ==============================================================================
#  PREPARE ETHNICITY-INTERACTION RESULTS
# ==============================================================================

race_plot_data <- race_interaction_whole_cohort %>%
  filter(
    predictor_short %in% names(prs_label_map),
    subgroup_level %in% c(
      "White",
      "Asian",
      "Black",
      "Mixed/Other"
    )
  ) %>%
  add_prs_positions() %>%
  mutate(
    subgroup = factor(
      subgroup_level,
      levels = c(
        "White",
        "Asian",
        "Black",
        "Mixed/Other"
      )
    )
  )


# ==============================================================================
#  CHECK FOR INVALID OR DUPLICATE RESULTS
# ==============================================================================

if (any(is.na(sex_plot_data$y))) {
  stop(
    "At least one predictor in sex_interaction_white could not be mapped ",
    "to a PRS display label."
  )
}

if (any(is.na(race_plot_data$y))) {
  stop(
    "At least one predictor in race_interaction_whole_cohort could not be ",
    "mapped to a PRS display label."
  )
}


sex_duplicates <- sex_plot_data %>%
  count(
    predictor_short,
    subgroup_level
  ) %>%
  filter(
    n > 1
  )

race_duplicates <- race_plot_data %>%
  count(
    predictor_short,
    subgroup_level
  ) %>%
  filter(
    n > 1
  )

if (nrow(sex_duplicates) > 0) {
  stop(
    "Duplicate sex-interaction rows were detected. ",
    "Filter sex_interaction_white to one model and one exclusion set."
  )
}

if (nrow(race_duplicates) > 0) {
  stop(
    "Duplicate ethnicity-interaction rows were detected. ",
    "Filter race_interaction_whole_cohort to one model and one exclusion set."
  )
}


# ==============================================================================
#  CHECK FOR MISSING PRS RESULTS
# ==============================================================================

missing_sex_prs <- setdiff(
  prs_order,
  unique(sex_plot_data$prs)
)

missing_race_prs <- setdiff(
  prs_order,
  unique(race_plot_data$prs)
)

if (length(missing_sex_prs) > 0) {
  warning(
    "Missing sex-interaction estimates for: ",
    paste(missing_sex_prs, collapse = ", ")
  )
}

if (length(missing_race_prs) > 0) {
  warning(
    "Missing ethnicity-interaction estimates for: ",
    paste(missing_race_prs, collapse = ", ")
  )
}


# ==============================================================================
#  PREPARE P-INTERACTION DATA
# ==============================================================================

sex_p_data <- sex_plot_data %>%
  distinct(
    predictor_short,
    prs,
    y,
    p_interaction
  ) %>%
  mutate(
    p_label = format_p_scientific(
      p_interaction
    ),
    
    p_strip = "italic(P)[interaction]"
  )


race_p_data <- race_plot_data %>%
  distinct(
    predictor_short,
    prs,
    y,
    p_interaction
  ) %>%
  mutate(
    p_label = format_p_scientific(
      p_interaction
    ),
    
    p_strip = "italic(P)[interaction]"
  )


# ==============================================================================
#  FOREST-PLOT FUNCTION
# ==============================================================================

make_forest_plot <- function(
    plot_data,
    x_limits,
    x_breaks,
    x_labels,
    panel_tag,
    subgroup_gap_cm = 0.20,
    show_y_labels = TRUE,
    left_margin_pt = 4,
    right_margin_pt = 0
) {
  p <- ggplot(
    plot_data,
    aes(
      x = or,
      y = y,
      colour = prs
    )
  ) +
    
    geom_rect(
      data = stripe_data,
      aes(
        xmin = x_limits[1],
        xmax = x_limits[2],
        ymin = ymin,
        ymax = ymax,
        fill = stripe_fill
      ),
      inherit.aes = FALSE,
      colour = NA
    ) +
    
    scale_fill_identity() +
    
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = reference_linewidth,
      colour = "grey35"
    ) +
    
    geom_errorbarh(
      aes(
        xmin = lb,
        xmax = ub
      ),
      height = 0,
      linewidth = ci_linewidth,
      na.rm = TRUE
    ) +
    
    geom_point(
      size = point_size,
      na.rm = TRUE
    ) +
    
    facet_grid(
      rows = vars(),
      cols = vars(subgroup),
      scales = "fixed",
      space = "fixed",
      drop = FALSE
    ) +
    
    scale_colour_manual(
      values = cluster_cols,
      limits = prs_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      breaks = seq_len(n_prs),
      labels = prs_y_order,
      limits = c(
        0.5,
        n_prs + 0.5
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    scale_x_log10(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    coord_cartesian(
      xlim = x_limits,
      ylim = c(
        0.5,
        n_prs + 0.5
      ),
      clip = "on"
    ) +
    
    labs(
      x = "Odds ratio for CAD",
      y = NULL,
      tag = panel_tag
    ) +
    
    theme_extended_fig3(
      subgroup_gap_cm = subgroup_gap_cm,
      
      plot_margin = margin(
        t = 4,
        r = right_margin_pt,
        b = 2,
        l = left_margin_pt
      )
    )
  
  if (!show_y_labels) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  p
}


# ==============================================================================
#  P-INTERACTION COLUMN FUNCTION
# ==============================================================================

make_p_column <- function(
    p_data,
    left_margin_pt = 0,
    right_margin_pt = 0
) {
  ggplot(
    p_data,
    aes(
      x = 1,
      y = y,
      label = p_label
    )
  ) +
    
    geom_text(
      family = "Arial",
      size = p_value_text_size_pt / ggplot2::.pt,
      colour = "black",
      hjust = 0.5,
      vjust = 0.5,
      na.rm = TRUE
    ) +
    
    facet_grid(
      rows = vars(),
      cols = vars(p_strip),
      labeller = label_parsed,
      drop = FALSE
    ) +
    
    scale_x_continuous(
      limits = c(
        0.5,
        1.5
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    scale_y_continuous(
      breaks = seq_len(n_prs),
      limits = c(
        0.5,
        n_prs + 0.5
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    coord_cartesian(
      ylim = c(
        0.5,
        n_prs + 0.5
      ),
      clip = "on"
    ) +
    
    labs(
      x = NULL,
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
      
      strip.background = element_rect(
        fill = NA,
        colour = NA
      ),
      
      strip.text = element_text(
        family = "Arial",
        face = "bold",
        size = facet_text_size_pt,
        colour = "black",
        margin = margin(
          t = 5,
          r = 5,
          b = 5,
          l = 5
        )
      ),
      
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.border = element_blank(),
      
      legend.position = "none",
      
      plot.margin = margin(
        t = 4,
        r = right_margin_pt,
        b = 2,
        l = left_margin_pt
      )
    )
}


# ==============================================================================
#  BUILD PANEL A: SEX
# ==============================================================================

panel_a_forest <- make_forest_plot(
  plot_data = sex_plot_data,
  
  x_limits = sex_x_limits,
  x_breaks = sex_x_breaks,
  x_labels = sex_x_labels,
  
  panel_tag = "a",
  
  subgroup_gap_cm = sex_subgroup_gap_cm,
  
  show_y_labels = TRUE,
  
  left_margin_pt = 4,
  right_margin_pt = 0
)


panel_a_p <- make_p_column(
  p_data = sex_p_data,
  
  left_margin_pt = 0,
  right_margin_pt = 0
)


panel_a_inner <- patchwork::wrap_plots(
  panel_a_forest,
  patchwork::plot_spacer(),
  panel_a_p,
  
  nrow = 1,
  
  widths = grid::unit.c(
    grid::unit(
      1,
      "null"
    ),
    
    grid::unit(
      forest_p_gap_cm,
      "cm"
    ),
    
    grid::unit(
      p_interaction_width_cm,
      "cm"
    )
  )
)


# Keep the nested panel as one item in the outer layout
panel_a <- patchwork::wrap_elements(
  full = panel_a_inner
)


# ==============================================================================
#  BUILD PANEL B: ETHNICITY
# ==============================================================================

panel_b_forest <- make_forest_plot(
  plot_data = race_plot_data,
  
  x_limits = race_x_limits,
  x_breaks = race_x_breaks,
  x_labels = race_x_labels,
  
  panel_tag = "b",
  
  subgroup_gap_cm = race_subgroup_gap_cm,
  
  show_y_labels = TRUE,
  
  left_margin_pt = 4,
  right_margin_pt = 0
)


panel_b_p <- make_p_column(
  p_data = race_p_data,
  
  left_margin_pt = 0,
  right_margin_pt = 0
)


panel_b_inner <- patchwork::wrap_plots(
  panel_b_forest,
  patchwork::plot_spacer(),
  panel_b_p,
  
  nrow = 1,
  
  widths = grid::unit.c(
    grid::unit(
      1,
      "null"
    ),
    
    grid::unit(
      forest_p_gap_cm,
      "cm"
    ),
    
    grid::unit(
      p_interaction_width_cm,
      "cm"
    )
  )
)


# Keep the nested panel as one item in the outer layout
panel_b <- patchwork::wrap_elements(
  full = panel_b_inner
)


# ==============================================================================
#  COMBINE PANELS IN ONE ROW
# ==============================================================================

extended_fig_3 <- patchwork::wrap_plots(
  panel_a,
  panel_b,
  
  nrow = 1,
  
  widths = c(
    sex_panel_width,
    race_panel_width
  )
)


# ==============================================================================
#  DISPLAY FIGURE
# ==============================================================================

print(
  extended_fig_3
)


# ==============================================================================
#  SAVE FIGURE
# ==============================================================================

ggsave(
  filename = "extended-fig-3.png",
  plot = extended_fig_3,
  width = 13,
  height = 4,
  units = "in",
  dpi = 1200,
  bg = "white"
)



# ==============================================================================
#  END OF SCRIPT
# ==============================================================================