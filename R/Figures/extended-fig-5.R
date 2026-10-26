# ==============================================================================
# Faceted volcano plots:
# Associations of eight process-specific PRSs with Olink proteins
#
# Input:
#   ukbb_olink_ols_all_markers
#
# Required columns:
#   outcome
#   outcome_category
#   predictor
#   beta
#   p_value
#
# Analysis:
#   - Retains outcome_category == "Proteomics"
#   - Retains the eight process-specific PRSs
#   - Uses outcome as the protein label
#   - Removes the leading "olink_" prefix from protein labels
#   - Applies a Bonferroni threshold of 0.05 / 2923 to each PRS panel
#   - Prints the number of Bonferroni-significant proteins in each facet strip
#   - Removes major and minor panel gridlines
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(tibble)
library(scales)
library(grid)
library(readr)

# ==============================================================================
# User-adjustable settings
# ==============================================================================

# ------------------------------------------------------------------------------
# Bonferroni correction
# ------------------------------------------------------------------------------

bonferroni_alpha <- 0.05
bonferroni_n_proteins <- 2923

bonferroni_threshold <- bonferroni_alpha / bonferroni_n_proteins
bonferroni_log10_threshold <- -log10(bonferroni_threshold)

# ------------------------------------------------------------------------------
# Protein labels
# ------------------------------------------------------------------------------

# Maximum number of proteins labelled in each panel
n_labels_per_panel <- 10

# Label only proteins significant at the Bonferroni threshold
label_significant_only <- TRUE

# Method used to select labelled proteins:
#   "p_value"  = lowest p-values
#   "abs_beta" = largest absolute standardized beta coefficients
label_ranking <- "p_value"

protein_label_wrap <- 18
protein_label_size <- 10 / ggplot2::.pt

label_max_overlaps <- Inf
label_box_padding <- 0.35
label_point_padding <- 0.15
label_segment_colour <- "grey45"
label_segment_width <- 0.30
label_seed <- 1

# ------------------------------------------------------------------------------
# Points
# ------------------------------------------------------------------------------

point_size <- 1.25
point_alpha <- 0.70

colour_non_significant <- "grey75"
colour_negative <- "#2166AC"
colour_positive <- "#B2182B"

# ------------------------------------------------------------------------------
# Axes
# ------------------------------------------------------------------------------

# TRUE gives all panels identical x-axis limits
use_common_x_limits <- TRUE

# Additional proportional padding around the largest absolute coefficient
x_axis_padding <- 0.08

# Optional manually specified symmetric x-axis limit
# Example:
# manual_x_limit <- 0.20
#
# Set to NULL for automatic calculation.
manual_x_limit <- NULL

# Additional upper-axis space for labels
y_axis_upper_expansion <- 0.10

# ------------------------------------------------------------------------------
# Facet layout and text
# ------------------------------------------------------------------------------

facet_ncol <- 4
facet_nrow <- 2

# Facet-box text size
facet_label_size <- 11

# All other plot text
plot_text_size <- 10

# Border thickness for panel and black facet boxes
box_linewidth <- 0.25

# ------------------------------------------------------------------------------
# Reference lines
# ------------------------------------------------------------------------------

threshold_line_colour <- "black"
threshold_line_type <- "dashed"
threshold_line_width <- 0.45

zero_line_colour <- "grey35"
zero_line_type <- "solid"
zero_line_width <- 0.35

# ------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------

plot_width <- 16
plot_height <- 9

output_png <- "olink_protein_volcano_plots_eight_process_prs.png"
output_pdf <- "olink_protein_volcano_plots_eight_process_prs.pdf"

output_significant_csv <-
  "significant_olink_proteins_by_process_prs.csv"

# ==============================================================================
# Validate settings
# ==============================================================================

if (!label_ranking %in% c("p_value", "abs_beta")) {
  stop(
    'label_ranking must be either "p_value" or "abs_beta".'
  )
}

if (
  !is.numeric(n_labels_per_panel) ||
  length(n_labels_per_panel) != 1 ||
  is.na(n_labels_per_panel) ||
  n_labels_per_panel < 0
) {
  stop(
    "n_labels_per_panel must be one non-negative numeric value."
  )
}

if (
  !is.numeric(bonferroni_n_proteins) ||
  length(bonferroni_n_proteins) != 1 ||
  is.na(bonferroni_n_proteins) ||
  bonferroni_n_proteins <= 0
) {
  stop(
    "bonferroni_n_proteins must be one positive numeric value."
  )
}

# ==============================================================================
# Process-specific PRS definitions
# ==============================================================================

cluster_map <- tibble::tribble(
  ~predictor,                   ~cluster_label,          ~cluster_order,
  "prs_ldl",                    "ChEL",                               1,
  "prs_trl",                    "TRL",                                2,
  "prs_lpa",                    "Lp(a)",                              3,
  "prs_mdwo",                   "MUNW",                               4,
  "prs_obesity",                "Obesity",                            5,
  "prs_endothelial",            "Endothelial",                        6,
  "prs_coagulation",            "Coagulation",                        7,
  "prs_cardiac_remodeling",     "Cardiac remodeling",                 8
)

cluster_predictors <- cluster_map$predictor

# ==============================================================================
# Check required input columns
# ==============================================================================

required_columns <- c(
  "outcome",
  "outcome_category",
  "predictor",
  "beta",
  "p_value"
)

missing_columns <- setdiff(
  required_columns,
  names(ukbb_olink_ols_all_markers)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "ukbb_olink_ols_all_markers is missing the following columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ==============================================================================
# Prepare Olink protein data
# ==============================================================================

olink_volcano_data <- ukbb_olink_ols_all_markers %>%
  mutate(
    outcome = trimws(as.character(outcome)),
    outcome_category = trimws(as.character(outcome_category)),
    predictor = trimws(as.character(predictor)),
    beta = suppressWarnings(as.numeric(beta)),
    p_value = suppressWarnings(as.numeric(p_value))
  ) %>%
  filter(
    outcome_category == "Proteomics",
    predictor %in% cluster_predictors,
    !is.na(outcome),
    outcome != "",
    !is.na(beta),
    is.finite(beta),
    !is.na(p_value),
    is.finite(p_value),
    p_value >= 0,
    p_value <= 1
  ) %>%
  mutate(
    # Use outcome and remove only a leading olink_ prefix
    protein_label = sub(
      pattern = "^olink_",
      replacement = "",
      x = outcome
    )
  )

if (nrow(olink_volcano_data) == 0) {
  stop(
    paste0(
      "No observations remained after filtering to ",
      'outcome_category == "Proteomics" and the eight process-specific PRSs.'
    )
  )
}

# ==============================================================================
# Check predictor availability
# ==============================================================================

missing_cluster_predictors <- cluster_map %>%
  anti_join(
    olink_volcano_data %>%
      distinct(predictor),
    by = "predictor"
  )

if (nrow(missing_cluster_predictors) > 0) {
  cat(
    "\nThe following process-specific PRSs were not found in the Proteomics data:\n"
  )
  
  print(
    missing_cluster_predictors,
    n = Inf,
    width = Inf
  )
}

# ==============================================================================
# Check duplicate protein-PRS combinations
# ==============================================================================

duplicate_combinations <- olink_volcano_data %>%
  count(
    predictor,
    outcome,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )

if (nrow(duplicate_combinations) > 0) {
  cat(
    "\nDuplicate protein-PRS combinations were detected:\n"
  )
  
  print(
    duplicate_combinations,
    n = Inf,
    width = Inf
  )
  
  stop(
    paste0(
      "The data contain more than one row for at least one protein-PRS ",
      "combination. Filter the input to one model, transformation, exclusion ",
      "set, or adjustment set before creating the volcano plots."
    )
  )
}

# ==============================================================================
# Calculate volcano-plot variables
# ==============================================================================

# p-values equal to zero cannot be transformed directly.
smallest_positive_p <- .Machine$double.xmin

olink_volcano_data <- olink_volcano_data %>%
  left_join(
    cluster_map,
    by = "predictor"
  ) %>%
  mutate(
    p_value_plot = pmax(
      p_value,
      smallest_positive_p
    ),
    
    minus_log10_p = -log10(p_value_plot),
    
    significant = p_value < bonferroni_threshold,
    
    association_group = case_when(
      significant & beta < 0 ~
        "Bonferroni P significant, Negative",
      
      significant & beta > 0 ~
        "Bonferroni P significant, Positive",
      
      significant & beta == 0 ~
        "Bonferroni P significant, Positive",
      
      TRUE ~
        "Not Bonferroni P significant"
    ),
    
    cluster_label = factor(
      cluster_label,
      levels = cluster_map$cluster_label
    )
  )

# ==============================================================================
# Count tested and significant proteins in each cluster
# ==============================================================================

cluster_counts <- olink_volcano_data %>%
  group_by(
    predictor,
    cluster_label,
    cluster_order
  ) %>%
  summarise(
    n_proteins_tested = n_distinct(outcome),
    n_significant = n_distinct(
      outcome[significant]
    ),
    .groups = "drop"
  ) %>%
  arrange(
    cluster_order
  ) %>%
  mutate(
    facet_label = paste0(
      as.character(cluster_label),
      "\nN significant proteins = ",
      format(
        n_significant,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      )
    )
  )

cat(
  "\nBonferroni threshold:",
  format(
    bonferroni_threshold,
    scientific = TRUE,
    digits = 4
  ),
  "\n"
)

cat(
  "Corresponding -log10(P) threshold:",
  round(
    bonferroni_log10_threshold,
    3
  ),
  "\n\n"
)

cat(
  "Protein counts by process-specific PRS:\n\n"
)

print(
  cluster_counts,
  n = Inf,
  width = Inf
)

# ==============================================================================
# Check whether every cluster contains 2,923 proteins
# ==============================================================================

unexpected_protein_counts <- cluster_counts %>%
  filter(
    n_proteins_tested != bonferroni_n_proteins
  )

if (nrow(unexpected_protein_counts) > 0) {
  cat(
    paste0(
      "\nWarning: some panels do not contain exactly ",
      bonferroni_n_proteins,
      " unique proteins:\n"
    )
  )
  
  print(
    unexpected_protein_counts,
    n = Inf,
    width = Inf
  )
}

# ==============================================================================
# Add two-line facet labels to the plotting data
# ==============================================================================

facet_label_map <- cluster_counts %>%
  select(
    predictor,
    facet_label,
    cluster_order
  )

olink_volcano_data <- olink_volcano_data %>%
  left_join(
    facet_label_map,
    by = c(
      "predictor",
      "cluster_order"
    )
  ) %>%
  mutate(
    facet_label = factor(
      facet_label,
      levels = cluster_counts$facet_label
    )
  )

# ==============================================================================
# Select proteins to label
# ==============================================================================

if (label_significant_only) {
  protein_label_candidates <- olink_volcano_data %>%
    filter(
      significant
    )
} else {
  protein_label_candidates <- olink_volcano_data
}

if (label_ranking == "p_value") {
  volcano_label_data <- protein_label_candidates %>%
    group_by(
      predictor
    ) %>%
    arrange(
      p_value,
      desc(abs(beta)),
      .by_group = TRUE
    ) %>%
    slice_head(
      n = n_labels_per_panel
    ) %>%
    ungroup()
} else {
  volcano_label_data <- protein_label_candidates %>%
    group_by(
      predictor
    ) %>%
    arrange(
      desc(abs(beta)),
      p_value,
      .by_group = TRUE
    ) %>%
    slice_head(
      n = n_labels_per_panel
    ) %>%
    ungroup()
}

volcano_label_data <- volcano_label_data %>%
  mutate(
    protein_label_plot = stringr::str_wrap(
      protein_label,
      width = protein_label_wrap
    )
  )

# ==============================================================================
# Define a common symmetric x-axis range
# ==============================================================================

if (!is.null(manual_x_limit)) {
  if (
    !is.numeric(manual_x_limit) ||
    length(manual_x_limit) != 1 ||
    is.na(manual_x_limit) ||
    manual_x_limit <= 0
  ) {
    stop(
      "manual_x_limit must be NULL or one positive numeric value."
    )
  }
  
  common_x_limit <- manual_x_limit
  
} else {
  
  common_x_limit <- max(
    abs(olink_volcano_data$beta),
    na.rm = TRUE
  )
  
  common_x_limit <- common_x_limit * (
    1 + x_axis_padding
  )
}

# ==============================================================================
# Legend colours and order
# ==============================================================================

association_levels <- c(
  "Bonferroni P significant, Negative",
  "Not Bonferroni P significant",
  "Bonferroni P significant, Positive"
)

association_colours <- c(
  "Bonferroni P significant, Negative" = colour_negative,
  "Not Bonferroni P significant" = colour_non_significant,
  "Bonferroni P significant, Positive" = colour_positive
)

olink_volcano_data <- olink_volcano_data %>%
  mutate(
    association_group = factor(
      association_group,
      levels = association_levels
    )
  )

# ==============================================================================
# Create faceted volcano plot
# ==============================================================================

olink_cluster_volcano_plot <- ggplot(
  olink_volcano_data,
  aes(
    x = beta,
    y = minus_log10_p
  )
) +
  
  # Vertical zero-effect reference line
  geom_vline(
    xintercept = 0,
    colour = zero_line_colour,
    linetype = zero_line_type,
    linewidth = zero_line_width
  ) +
  
  # Bonferroni threshold
  geom_hline(
    yintercept = bonferroni_log10_threshold,
    colour = threshold_line_colour,
    linetype = threshold_line_type,
    linewidth = threshold_line_width
  ) +
  
  # Protein points
  geom_point(
    aes(
      colour = association_group
    ),
    size = point_size,
    alpha = point_alpha,
    na.rm = TRUE
  ) +
  
  # Labels for selected proteins
  ggrepel::geom_text_repel(
    data = volcano_label_data,
    aes(
      label = protein_label_plot
    ),
    family = "Arial",
    size = protein_label_size,
    colour = "black",
    box.padding = label_box_padding,
    point.padding = label_point_padding,
    min.segment.length = 0,
    segment.colour = label_segment_colour,
    segment.linewidth = label_segment_width,
    max.overlaps = label_max_overlaps,
    seed = label_seed,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    vars(facet_label),
    ncol = facet_ncol,
    nrow = facet_nrow,
    scales = if (use_common_x_limits) {
      "fixed"
    } else {
      "free_x"
    }
  ) +
  
  scale_colour_manual(
    values = association_colours,
    breaks = association_levels,
    labels = association_levels,
    drop = FALSE,
    name = "Association"
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.02,
        y_axis_upper_expansion
      )
    )
  ) +
  
  labs(
    x = "Standardized beta",
    y = expression(-log[10](italic(P))),
    title = NULL,
    subtitle = NULL,
    caption = NULL
  ) +
  
  coord_cartesian(
    xlim = if (use_common_x_limits) {
      c(
        -common_x_limit,
        common_x_limit
      )
    } else {
      NULL
    },
    clip = "off"
  ) +
  
  theme_linedraw(
    base_size = plot_text_size,
    base_family = "Arial"
  ) +
  
  theme(
    text = element_text(
      family = "Arial",
      size = plot_text_size,
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = plot_text_size,
      colour = "black"
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = plot_text_size,
      colour = "black"
    ),
    
    # Remove all panel gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = box_linewidth
    ),
    
    strip.background = element_rect(
      fill = "black",
      colour = "black",
      linewidth = box_linewidth
    ),
    
    strip.text = element_text(
      family = "Arial",
      face = "bold",
      colour = "white",
      size = facet_label_size,
      lineheight = 1.10,
      margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = plot_text_size
    ),
    
    legend.text = element_text(
      family = "Arial",
      size = plot_text_size
    ),
    
    legend.key.width = grid::unit(
      1.5,
      "lines"
    ),
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 20,
      b = 10,
      l = 10
    )
  )

# ==============================================================================
# Print plot
# ==============================================================================

print(
  olink_cluster_volcano_plot
)

# ==============================================================================
# Print significant proteins
# ==============================================================================

significant_olink_proteins <- olink_volcano_data %>%
  filter(
    significant
  ) %>%
  arrange(
    cluster_order,
    p_value
  ) %>%
  select(
    cluster_label,
    predictor,
    outcome,
    protein_label,
    beta,
    p_value,
    minus_log10_p
  )

cat(
  "\nBonferroni-significant Olink protein associations:\n\n"
)

print(
  significant_olink_proteins,
  n = Inf,
  width = Inf
)

# ==============================================================================
# Save significant-protein table
# ==============================================================================

readr::write_csv(
  significant_olink_proteins,
  output_significant_csv
)

# ==============================================================================
# Save figure
# ==============================================================================

ggsave(
  filename = "extended-fig-5.png",
  plot = olink_cluster_volcano_plot,
  width = 12,
  height = 7,
  units = "in",
  dpi = 1200,
  bg = "white",
  limitsize = FALSE
)

