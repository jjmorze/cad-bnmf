# ==============================================================================
# Two-row horizontal heatmap of remaining selected downstream biomarkers
#
# Row 1:
#   Abdominal MRI, Anthropometry, DXA body composition
#
# Row 2:
#   Biochemistry, Blood count, Cardiac MRI, ECG, Urine biomarkers
#
# Input:
#   ukbb_selected_biomarker_associations
#
# Required columns:
#   outcome
#   predictor
#   beta
#   p_value
#   p_fdr
#
# Notes:
#   - Pulse rate is excluded.
#   - Biomarkers are uniquely identified using outcome.
#   - Existing p_fdr values are used; FDR is not recalculated.
#   - Tile colour uses sign(beta) * sqrt(abs(beta)).
#   - Legend labels show the original beta scale.
#   - Original beta values are printed only for significant associations.
#   - Coefficient labels are white when abs(beta) > 0.20.
#   - Non-significant tiles are grey80 without coefficient labels.
#   - Category boxes are positioned above the heatmaps.
#   - Adjustable whitespace separates pathway PRSs from Total PRS.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)
library(tibble)
library(stringr)
library(patchwork)

# ==============================================================================
# User-adjustable settings
# ==============================================================================

significance_method <- "Bonferroni"
# Options:
#   "FDR"
#   "Bonferroni"

fdr_alpha <- 0.05

bonferroni_alpha <- 0.05
bonferroni_n_prs <- 1
bonferroni_n_biomarkers <- 83

# ------------------------------------------------------------------------------
# Tile dimensions
# ------------------------------------------------------------------------------

tile_width <- 0.92
tile_height <- 0.92

# ------------------------------------------------------------------------------
# Text sizes
# ------------------------------------------------------------------------------

tile_text_size <- 3.5
biomarker_label_size <- 10
prs_label_size <- 10
category_label_size <- 11
legend_title_size <- 10
legend_text_size <- 10

# Use white coefficient labels when absolute beta exceeds this value
white_label_threshold <- 0.20

# Wrap long biomarker labels
biomarker_label_wrap <- 28

# Biomarker-axis label angle
biomarker_label_angle <- 45

# ------------------------------------------------------------------------------
# Spacing controls
# ------------------------------------------------------------------------------

# Horizontal space between category panels
category_spacing <- 0.65

# Vertical whitespace between the two complete heatmap rows, in points
heatmap_row_spacing_pt <- 4

# Additional whitespace above the upper heatmap, in points
# Increase this to move the upper black category boxes downward.
top_row_top_whitespace_pt <- 20

# Padding within the black category boxes
category_box_top_padding_pt <- 10
category_box_bottom_padding_pt <- 10
category_box_side_padding_pt <- 10

# Adjustable vertical gap between pathway PRSs and Total PRS.
# 0 gives ordinary row spacing.
# Suggested values: 0.25 to 1.00.
total_prs_row_gap <- 0.25

# Relative heights of the two heatmap rows
row_relative_heights <- c(1, 1)

# ------------------------------------------------------------------------------
# Colour scale
# ------------------------------------------------------------------------------

# Maximum absolute original beta represented by the colour scale
beta_colour_limit <- 0.30

# ------------------------------------------------------------------------------
# Output dimensions
# ------------------------------------------------------------------------------

plot_width <- 26
plot_height <- 15

# ==============================================================================
# Validate significance method
# ==============================================================================

if (!significance_method %in% c("FDR", "Bonferroni")) {
  stop(
    'significance_method must be either "FDR" or "Bonferroni".'
  )
}

bonferroni_threshold <- bonferroni_alpha / (
  bonferroni_n_prs * bonferroni_n_biomarkers
)

# ==============================================================================
# PRSs to retain, order, labels, and numeric y positions
# ==============================================================================

prs_keep <- c(
  "prs_ldl",
  "prs_trl",
  "prs_lpa",
  "prs_mdwo",
  "prs_obesity",
  "prs_endothelial",
  "prs_coagulation",
  "prs_cardiac_remodeling",
  "prs_total_cluster"
)

prs_label_map <- c(
  prs_ldl = "ChEL",
  prs_trl = "TRL",
  prs_lpa = "Lp(a)",
  prs_mdwo = "MUNW",
  prs_obesity = "Obesity",
  prs_endothelial = "Endothelial dysfunction",
  prs_coagulation = "Coagulation",
  prs_cardiac_remodeling = "Cardiac remodeling",
  prs_total_cluster = "Total CAD PRS"
)

# Numeric y positions create an adjustable gap before Total PRS.
#
# The first pathway PRS is shown at the top because the y scale is reversed.
prs_position_map <- c(
  prs_ldl = 1,
  prs_trl = 2,
  prs_lpa = 3,
  prs_mdwo = 4,
  prs_obesity = 5,
  prs_endothelial = 6,
  prs_coagulation = 7,
  prs_cardiac_remodeling = 8,
  prs_total_cluster = 9 + total_prs_row_gap
)

prs_breaks <- unname(
  prs_position_map[prs_keep]
)

prs_labels <- unname(
  prs_label_map[prs_keep]
)

# ==============================================================================
# Manual biomarker labels, categories, and order
#
# bp_mean_pulse_r is deliberately omitted.
# ==============================================================================

remaining_biomarker_map <- tibble::tribble(
  ~outcome,                      ~biomarker_label,                            ~outcome_category,      ~category_order, ~biomarker_order, ~plot_row,
  
  # ---------------------------------------------------------------------------
  # Row 1: Abdominal MRI
  # ---------------------------------------------------------------------------
  
  "amri_asat",                   "ASAT volume",      "Abdominal MRI",          1,  1, 1,
  "amri_asat_bmiadj",            "ASAT (BMIadj)",  "Abdominal MRI",          1,  2, 1,
  "amri_gfat",                   "GFAT volume",                          "Abdominal MRI",          1,  3, 1,
  "amri_gfat_bmiadj",            "GFAT volume (BMIadj)",           "Abdominal MRI",          1,  4, 1,
  "amri_liver_pdff_fr",          "Liver PDFF",                        "Abdominal MRI",          1,  5, 1,
  "amri_liver_pdff_fr_bmiadj",   "Liver PDFF (BMIadj)",         "Abdominal MRI",          1,  6, 1,
  "amri_tat",                    "Total adipose tissue volume",                "Abdominal MRI",          1,  7, 1,
  "amri_vat",                    "VAT volume",             "Abdominal MRI",          1,  8, 1,
  "amri_vat_bmiadj",             "VAT (BMIadj)",     "Abdominal MRI",          1,  9, 1,
  
  # ---------------------------------------------------------------------------
  # Row 1: Anthropometry
  # ---------------------------------------------------------------------------
  
  "anthro_hc",                   "Hip circumference",                          "Anthropometry",          2,  1, 1,
  "anthro_hc_bmiadj",            "Hip circumference (BMIadj)",           "Anthropometry",          2,  2, 1,
  "anthro_height",               "Height",                            "Anthropometry",          2,  3, 1,
  "anthro_wc",                   "Waist circumference",                        "Anthropometry",          2,  4, 1,
  "anthro_wc_bmiadj",            "Waist circumference (BMIadj)",         "Anthropometry",          2,  5, 1,
  "anthro_whr",                  "Waist-to-hip ratio",                         "Anthropometry",          2,  6, 1,
  "anthro_whr_bmiadj",           "Waist-to-hip ratio (BMIadj)",          "Anthropometry",          2,  7, 1,
  "anthro_wt",                   "Weight",                                     "Anthropometry",          2,  8, 1,
  
  # ---------------------------------------------------------------------------
  # Row 1: DXA body composition
  # ---------------------------------------------------------------------------
  
  "dxa_ag_ratio",                "Android-to-gynoid fat ratio",                "DXA body composition",   3,  1, 1,
  "dxa_android_fat",             "Android fat mass",                           "DXA body composition",   3,  2, 1,
  "dxa_android_fat_frac",        "Android fat fraction",                       "DXA body composition",   3,  3, 1,
  "dxa_gynoid_fat",              "Gynoid fat mass",                            "DXA body composition",   3,  4, 1,
  "dxa_gynoid_fat_frac",         "Gynoid fat fraction",                        "DXA body composition",   3,  5, 1,
  "dxa_legs_fat",                "Leg fat mass",                               "DXA body composition",   3,  6, 1,
  "dxa_legs_fat_frac",           "Leg fat fraction",                           "DXA body composition",   3,  7, 1,
  "dxa_total_bone",              "Total bone mass",                            "DXA body composition",   3,  8, 1,
  "dxa_total_fat",               "Total fat mass",                             "DXA body composition",   3,  9, 1,
  "dxa_total_fat_frac",          "Total fat fraction",                         "DXA body composition",   3, 10, 1,
  "dxa_total_lean",              "Total lean mass",                            "DXA body composition",   3, 11, 1,
  "dxa_trunk_fat",               "Trunk fat mass",                             "DXA body composition",   3, 12, 1,
  "dxa_trunk_fat_frac",          "Trunk fat fraction",                         "DXA body composition",   3, 13, 1,
  "dxa_trunk_leg_fat_ratio",     "Trunk-to-leg fat ratio",                     "DXA body composition",   3, 14, 1,
  "dxa_vat_mass",                "VAT mass",               "DXA body composition",   3, 15, 1,
  "dxa_vat_vol",                 "VAT volume",             "DXA body composition",   3, 16, 1,
  
  # ---------------------------------------------------------------------------
  # Row 2: Biochemistry
  # ---------------------------------------------------------------------------
  
  "conven_albumin",              "Albumin",                                    "Biochemistry",           4,  1, 2,
  "conven_alt",                  "Alanine aminotransferase",                    "Biochemistry",           4,  2, 2,
  "conven_apoa_therapy_corr",    "Apolipoprotein A",                           "Biochemistry",           4,  3, 2,
  "conven_ca",                   "Calcium",                                    "Biochemistry",           4,  4, 2,
  "conven_igf1",                 "Insulin growth factor-1",                                      "Biochemistry",           4,  5, 2,
  "conven_nonhdlc_therapy_corr", "Non-HDL cholesterol",                        "Biochemistry",           4,  6, 2,
  "conven_phos",                 "Phosphate",                                  "Biochemistry",           4,  7, 2,
  "conven_rchol_therapy_corr",   "Remnant cholesterol",                        "Biochemistry",           4,  8, 2,
  "conven_tchol_therapy_corr",   "Total cholesterol",                          "Biochemistry",           4,  9, 2,
  "conven_urea",                 "Urea",                                       "Biochemistry",           4, 10, 2,
  
  # ---------------------------------------------------------------------------
  # Row 2: Blood count
  # ---------------------------------------------------------------------------
  
  "blood_hgb",                   "Hemoglobin",                                "Blood count",            5,  1, 2,
  "blood_lymph",                 "Lymphocyte count",                           "Blood count",            5,  2, 2,
  "blood_neut",                  "Neutrophil count",                           "Blood count",            5,  3, 2,
  "blood_plt",                   "Platelet count",                             "Blood count",            5,  4, 2,
  "blood_rbc",                   "Red blood cell count",                       "Blood count",            5,  5, 2,
  "blood_rdw",                   "Red blood cell distribution width",                "Blood count",            5,  6, 2,
  
  # ---------------------------------------------------------------------------
  # Row 2: Cardiac MRI
  # ---------------------------------------------------------------------------
  
  "cmr_lvedv_bsaadj",            "LVEDV index",              "Cardiac MRI",            6,  1, 2,
  "cmr_lvef",                    "LV ejection fraction",                       "Cardiac MRI",            6,  2, 2,
  "cmr_lvesv_bsaadj",            "LVESV index",               "Cardiac MRI",            6,  3, 2,
  "cmr_lvmass_bsaadj",           "LV mass index",                              "Cardiac MRI",            6,  4, 2,
  "cmr_lvmvr",                   "LV mass-to-volume ratio",                     "Cardiac MRI",            6,  5, 2,
  "cmr_lvsv_bsaadj",             "LV stroke volume index",                     "Cardiac MRI",            6,  6, 2,
  "cmr_pericardial_fat",         "Pericardial fat",                            "Cardiac MRI",            6,  7, 2,
  "cmr_pericardial_fat_bmiadj",  "Pericardial fat (BMIadj)",             "Cardiac MRI",            6,  8, 2,
  
  # ---------------------------------------------------------------------------
  # Row 2: ECG
  # ---------------------------------------------------------------------------
  
  "ecg_p_dur",                   "P-wave duration",                            "ECG",                    7,  1, 2,
  "ecg_qrs_dur",                 "QRS duration",                               "ECG",                    7,  2, 2,
  "ecg_qtc",                     "QTc interval",                               "ECG",                    7,  3, 2,
  
  # ---------------------------------------------------------------------------
  # Row 2: Urine biomarkers
  # ---------------------------------------------------------------------------
  
  "urine_k",                     "Urine potassium",                            "Urinary",       8,  1, 2,
  "urine_na",                    "Urine sodium",                               "Urinary",       8,  2, 2,
  "urine_uacr",                  "UACR",           "Urinary",       8,  3, 2
)

# ==============================================================================
# Check required input columns
# ==============================================================================

required_columns <- c(
  "outcome",
  "predictor",
  "beta",
  "p_value",
  "p_fdr"
)

missing_columns <- setdiff(
  required_columns,
  names(ukbb_selected_biomarker_associations)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "ukbb_selected_biomarker_associations is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ==============================================================================
# Prepare source data
# ==============================================================================

source_data <- ukbb_selected_biomarker_associations %>%
  mutate(
    outcome = trimws(as.character(outcome)),
    predictor = trimws(as.character(predictor)),
    beta = suppressWarnings(as.numeric(beta)),
    p_value = suppressWarnings(as.numeric(p_value)),
    p_fdr = suppressWarnings(as.numeric(p_fdr))
  ) %>%
  filter(
    outcome %in% remaining_biomarker_map$outcome,
    predictor %in% prs_keep
  ) %>%
  select(
    -any_of(
      c(
        "outcome_long",
        "outcome_category"
      )
    )
  ) %>%
  left_join(
    remaining_biomarker_map,
    by = "outcome"
  )

# ==============================================================================
# Check for duplicate outcome-PRS combinations
# ==============================================================================

duplicate_combinations <- source_data %>%
  count(
    outcome,
    predictor,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )

if (nrow(duplicate_combinations) > 0) {
  print(
    duplicate_combinations,
    n = Inf,
    width = Inf
  )
  
  stop(
    paste0(
      "Duplicate outcome-PRS combinations were found. ",
      "Filter the input data to one model and one exclusion set before plotting."
    )
  )
}

# ==============================================================================
# Check mapped biomarkers not found in results
# ==============================================================================

missing_mapped_biomarkers <- remaining_biomarker_map %>%
  anti_join(
    source_data %>%
      distinct(outcome),
    by = "outcome"
  )

if (nrow(missing_mapped_biomarkers) > 0) {
  cat(
    "\nMapped biomarkers not found in the results:\n"
  )
  
  print(
    missing_mapped_biomarkers,
    n = Inf,
    width = Inf
  )
}

# ==============================================================================
# Biomarker metadata
# ==============================================================================

biomarker_metadata <- remaining_biomarker_map %>%
  semi_join(
    source_data %>%
      distinct(outcome),
    by = "outcome"
  ) %>%
  arrange(
    plot_row,
    category_order,
    biomarker_order
  )

# ==============================================================================
# Create complete biomarker-PRS grid
# ==============================================================================

heatmap_data <- tidyr::crossing(
  outcome = biomarker_metadata$outcome,
  predictor = prs_keep
) %>%
  left_join(
    source_data,
    by = c(
      "outcome",
      "predictor"
    )
  ) %>%
  select(
    -any_of(
      c(
        "biomarker_label",
        "outcome_category",
        "category_order",
        "biomarker_order",
        "plot_row"
      )
    )
  ) %>%
  left_join(
    biomarker_metadata,
    by = "outcome"
  )

# ==============================================================================
# Define significance, transformed beta, and coefficient labels
# ==============================================================================

heatmap_data <- heatmap_data %>%
  mutate(
    significant = case_when(
      significance_method == "FDR" ~
        !is.na(p_fdr) &
        p_fdr < fdr_alpha,
      
      significance_method == "Bonferroni" ~
        !is.na(p_value) &
        p_value < bonferroni_threshold,
      
      TRUE ~ FALSE
    ),
    
    # Signed square-root transformation used for colour
    beta_sqrt = sign(beta) * sqrt(abs(beta)),
    
    # Non-significant associations are grey
    plotted_beta = if_else(
      significant %in% TRUE,
      beta_sqrt,
      NA_real_
    ),
    
    # Original beta printed only for significant associations
    beta_label = case_when(
      !significant %in% TRUE ~ "",
      is.na(beta) ~ "",
      TRUE ~ sub(
        "^(-?)0\\.",
        "\\1.",
        sprintf("%.2f", beta)
      )
    ),
    
    beta_label_colour = case_when(
      !significant %in% TRUE ~ NA_character_,
      is.na(beta) ~ NA_character_,
      abs(beta) > white_label_threshold ~ "white",
      TRUE ~ "black"
    ),
    
    # Numeric PRS position allows an adjustable gap before Total PRS
    prs_position = unname(
      prs_position_map[predictor]
    )
  )

# ==============================================================================
# Print biomarkers included in each row
# ==============================================================================

biomarkers_in_plot <- biomarker_metadata %>%
  select(
    plot_row,
    outcome_category,
    outcome,
    biomarker_label
  )

cat(
  "\nBiomarkers included in the two-row heatmap:\n\n"
)

print(
  biomarkers_in_plot,
  n = Inf,
  width = Inf
)

# ==============================================================================
# Check missing biomarker-PRS combinations
# ==============================================================================

missing_combinations <- heatmap_data %>%
  filter(
    is.na(beta)
  ) %>%
  select(
    plot_row,
    outcome_category,
    outcome,
    biomarker_label,
    predictor
  )

cat(
  "\nNumber of missing biomarker-PRS combinations:",
  nrow(missing_combinations),
  "\n"
)

if (nrow(missing_combinations) > 0) {
  print(
    missing_combinations,
    n = Inf,
    width = Inf
  )
}

# ==============================================================================
# Signed square-root colour scale
#
# Break positions are transformed.
# Labels show original beta values.
# ==============================================================================

legend_beta_breaks <- c(
  -0.30,
  -0.10,
  -0.01,
  0.01,
  0.10,
  0.30
)

legend_transformed_breaks <- sign(
  legend_beta_breaks
) * sqrt(
  abs(legend_beta_breaks)
)

legend_beta_labels <- c(
  "-.30",
  "-.10",
  "-.01",
  ".01",
  ".10",
  ".30"
)

transformed_colour_limit <- sqrt(
  beta_colour_limit
)

# ==============================================================================
# Plot-building function
# ==============================================================================

make_heatmap_row <- function(
    data,
    row_number,
    show_legend = TRUE,
    top_margin_pt = 0,
    bottom_margin_pt = 0
) {
  
  row_metadata <- biomarker_metadata %>%
    filter(
      plot_row == row_number
    ) %>%
    arrange(
      category_order,
      biomarker_order
    )
  
  row_categories <- row_metadata %>%
    distinct(
      outcome_category,
      category_order
    ) %>%
    arrange(
      category_order
    ) %>%
    pull(
      outcome_category
    )
  
  row_biomarker_levels <- row_metadata %>%
    pull(
      outcome
    )
  
  row_label_lookup <- setNames(
    row_metadata$biomarker_label,
    row_metadata$outcome
  )
  
  plot_data <- data %>%
    filter(
      plot_row == row_number
    ) %>%
    mutate(
      outcome_category = factor(
        outcome_category,
        levels = row_categories
      ),
      
      biomarker_outcome = factor(
        outcome,
        levels = row_biomarker_levels
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = biomarker_outcome,
      y = prs_position,
      fill = plotted_beta
    )
  ) +
    
    geom_tile(
      width = tile_width,
      height = tile_height,
      colour = "white",
      linewidth = 0.35
    ) +
    
    geom_text(
      aes(
        label = beta_label,
        colour = beta_label_colour
      ),
      family = "Arial",
      size = tile_text_size,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    
    scale_colour_identity() +
    
    # Category boxes above the heatmap
    facet_grid(
      cols = vars(outcome_category),
      scales = "free_x",
      space = "free_x"
    ) +
    
    scale_x_discrete(
      labels = function(x) {
        stringr::str_wrap(
          unname(
            row_label_lookup[x]
          ),
          width = biomarker_label_wrap
        )
      },
      expand = expansion(
        add = 0
      )
    ) +
    
    # Reverse the numeric scale so ChEL is at the top.
    # Numeric positions create a gap before Total PRS.
    scale_y_reverse(
      breaks = prs_breaks,
      labels = prs_labels,
      limits = c(
        max(prs_breaks) + tile_height / 2,
        min(prs_breaks) - tile_height / 2
      ),
      expand = expansion(
        mult = 0,
        add = 0
      )
    ) +
    
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      
      limits = c(
        -transformed_colour_limit,
        transformed_colour_limit
      ),
      
      breaks = legend_transformed_breaks,
      labels = legend_beta_labels,
      
      oob = scales::squish,
      na.value = "grey80",
      name = "Standardized\nbeta"
    ) +
    
    labs(
      x = NULL,
      y = NULL,
      title = NULL,
      subtitle = NULL,
      caption = NULL
    ) +
    
    guides(
      fill = if (show_legend) {
        guide_colourbar(
          title.position = "top",
          title.hjust = 0.5
        )
      } else {
        "none"
      }
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_linedraw(
      base_size = 11,
      base_family = "Arial"
    ) +
    
    theme(
      text = element_text(
        family = "Arial"
      ),
      
      axis.title = element_blank(),
      
      axis.text.x = element_text(
        family = "Arial",
        size = biomarker_label_size,
        colour = "black",
        angle = biomarker_label_angle,
        hjust = 1,
        vjust = 1
      ),
      
      axis.text.y = element_text(
        family = "Arial",
        size = prs_label_size,
        colour = "black"
      ),
      
      axis.ticks.y = element_blank(),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.35
      ),
      
      # Category boxes above the heatmap
      strip.background.x = element_rect(
        fill = "black",
        colour = "black",
        linewidth = 0.35
      ),
      
      strip.text.x = element_text(
        family = "Arial",
        face = "bold",
        colour = "white",
        size = category_label_size,
        margin = margin(
          t = category_box_top_padding_pt,
          r = category_box_side_padding_pt,
          b = category_box_bottom_padding_pt,
          l = category_box_side_padding_pt,
          unit = "pt"
        )
      ),
      
      strip.placement = "inside",
      
      panel.spacing.x = grid::unit(
        category_spacing,
        "lines"
      ),
      
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      
      # Adjustable outer spacing for each complete heatmap row
      plot.margin = margin(
        t = top_margin_pt,
        r = 5,
        b = bottom_margin_pt,
        l = 5,
        unit = "pt"
      ),
      
      legend.title = element_text(
        family = "Arial",
        face = "bold",
        size = legend_title_size
      ),
      
      legend.text = element_text(
        family = "Arial",
        size = legend_text_size
      ),
      
      legend.position = if (show_legend) {
        "right"
      } else {
        "none"
      }
    )
}

# ==============================================================================
# Build the two heatmap rows
# ==============================================================================

heatmap_row_1 <- make_heatmap_row(
  data = heatmap_data,
  row_number = 1,
  show_legend = TRUE,
  top_margin_pt = top_row_top_whitespace_pt,
  bottom_margin_pt = heatmap_row_spacing_pt / 2
)

heatmap_row_2 <- make_heatmap_row(
  data = heatmap_data,
  row_number = 2,
  show_legend = TRUE,
  top_margin_pt = heatmap_row_spacing_pt / 2,
  bottom_margin_pt = 0
)

# ==============================================================================
# Combine into a two-row figure
# ==============================================================================

remaining_biomarker_heatmap <- (
  heatmap_row_1 /
    heatmap_row_2
) +
  plot_layout(
    heights = row_relative_heights,
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

# ==============================================================================
# Print
# ==============================================================================

print(
  remaining_biomarker_heatmap
)

# ==============================================================================
# Save
# ==============================================================================

ggsave(
  filename = "extended-fig-4-100926.png",
  plot = remaining_biomarker_heatmap,
  width = 13,
  height = 10,
  units = "in",
  dpi = 1200,
  bg = "white",
  limitsize = FALSE
)

