# ==============================================================================
# Shiny App: Four-panel biomarker / process-specific PRS figure builder
#
# Panel A: selected downstream-biomarker heatmap
# Panel B: MUNW-versus-obesity effect-size scatter plot
# Forest plot (displayed as Panel D) with selectable PRSs
# Protein heatmap (displayed as Panel C)
# ==============================================================================

options(shiny.maxRequestSize = 1024 * 1024^2)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(DT)
library(patchwork)
library(grid)
library(jsonlite)
library(tibble)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

stop_user <- function(message) stop(message, call. = FALSE)

# ==============================================================================
# Predictor aliases, labels, and colours
# ==============================================================================

PREDICTOR_ALIASES <- c(
  chol_prs_std = "prs_ldl",
  trl_prs_std = "prs_trl",
  lpa_prs_std = "prs_lpa",
  mdwo_prs_std = "prs_mdwo",
  obesity_prs_std = "prs_obesity",
  obesity_prs = "prs_obesity",
  endothelial_prs_std = "prs_endothelial",
  cardiac_prs_std = "prs_cardiac_remodeling",
  coagulation_prs_std = "prs_coagulation",
  total_196_prs_std = "prs_total",
  total_prs_std = "prs_total",
  prs_total_cluster = "prs_total"
)

DEFAULT_PREDICTOR_LABELS <- c(
  prs_ldl = "ChEL",
  prs_trl = "TRL",
  prs_lpa = "Lp(a)",
  prs_mdwo = "MUNW",
  prs_obesity = "Obesity",
  prs_endothelial = "Endothelial",
  prs_coagulation = "Coagulation",
  prs_cardiac_remodeling = "Cardiac remodeling",
  prs_total = "Total PRS"
)

DEFAULT_PREDICTOR_COLORS <- c(
  prs_ldl = "#B22D4C",
  prs_trl = "#E07A8E",
  prs_lpa = "#7A1F3A",
  prs_mdwo = "#4A7C4E",
  prs_obesity = "#95C075",
  prs_endothelial = "#2D5F8B",
  prs_coagulation = "#7AAED1",
  prs_cardiac_remodeling = "#5D4787",
  prs_total = "grey20"
)

ALL_PRS <- names(DEFAULT_PREDICTOR_LABELS)

resolve_predictors <- function(x) {
  x <- trimws(as.character(x))
  mapped <- unname(PREDICTOR_ALIASES[x])
  out <- ifelse(!is.na(mapped), mapped, x)
  unique(out[nzchar(out)])
}

resolve_predictor_column <- function(x) {
  x <- trimws(as.character(x))
  mapped <- unname(PREDICTOR_ALIASES[x])
  ifelse(!is.na(mapped), mapped, x)
}

parse_lines <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character())
  ans <- trimws(unlist(strsplit(x, "\n", fixed = TRUE)))
  ans[nzchar(ans)]
}

parse_named_mapping <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(c())
  lines <- trimws(unlist(strsplit(text, "\n", fixed = TRUE)))
  lines <- lines[nzchar(lines)]
  parsed <- lapply(lines, function(line) {
    pieces <- unlist(strsplit(line, "=", fixed = TRUE))
    if (length(pieces) < 2) return(NULL)
    c(key = trimws(pieces[1]), value = trimws(paste(pieces[-1], collapse = "=")))
  })
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (!length(parsed)) return(c())
  setNames(
    vapply(parsed, function(x) x[["value"]], character(1)),
    vapply(parsed, function(x) x[["key"]], character(1))
  )
}

parse_numeric_vector <- function(x) {
  vals <- suppressWarnings(as.numeric(trimws(unlist(strsplit(x %||% "", ",", fixed = TRUE)))))
  vals[is.finite(vals)]
}

parse_xy_offset_mapping <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(tibble::tibble(outcome = character(), label_dx = numeric(), label_dy = numeric()))
  }

  lines <- trimws(unlist(strsplit(text, "\n", fixed = TRUE)))
  lines <- lines[nzchar(lines)]

  parsed <- lapply(lines, function(line) {
    pieces <- unlist(strsplit(line, "=", fixed = TRUE))
    if (length(pieces) < 2) return(NULL)

    outcome <- trimws(pieces[1])
    offsets <- trimws(unlist(strsplit(paste(pieces[-1], collapse = "="), ",", fixed = TRUE)))
    if (!nzchar(outcome) || length(offsets) < 1) return(NULL)

    dx <- suppressWarnings(as.numeric(offsets[1]))
    dy <- if (length(offsets) >= 2) suppressWarnings(as.numeric(offsets[2])) else NA_real_

    tibble::tibble(
      outcome = outcome,
      label_dx = if (is.finite(dx)) dx else NA_real_,
      label_dy = if (is.finite(dy)) dy else NA_real_
    )
  })

  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (!length(parsed)) {
    return(tibble::tibble(outcome = character(), label_dx = numeric(), label_dy = numeric()))
  }

  dplyr::bind_rows(parsed) %>%
    dplyr::distinct(outcome, .keep_all = TRUE)
}

format_beta_ticks <- function(x) {
  formatC(x, format = "f", digits = 2)
}

merge_map <- function(default, custom) {
  output <- default
  if (length(custom)) output[names(custom)] <- unname(custom)
  output
}

label_lookup <- function(map, keys) {
  out <- unname(map[keys])
  bad <- is.na(out) | !nzchar(out)
  out[bad] <- keys[bad]
  out
}

color_lookup <- function(map, keys) {
  out <- unname(map[keys])
  bad <- is.na(out) | !nzchar(out)
  out[bad] <- "grey40"
  out
}

safe_read_csv <- function(path) {
  if (tolower(tools::file_ext(path)) == "gz") {
    read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# ==============================================================================
# Default panel selections
# ==============================================================================

PRS_KEEP_A <- c(
  "prs_ldl", "prs_trl", "prs_lpa", "prs_mdwo", "prs_obesity",
  "prs_endothelial", "prs_coagulation", "prs_cardiac_remodeling", "prs_total"
)

BIOMARKER_MAP_A <- tibble::tribble(
  ~outcome,                   ~biomarker_label,                 ~block,                 ~biomarker_order,
  "conven_apob_therapy_corr", "Apolipoprotein B",               "Lipid metabolism",      1,
  "conven_ldl_therapy_corr",  "LDL cholesterol",                "Lipid metabolism",      2,
  "conven_tg_therapy_corr",   "Triglycerides",                  "Lipid metabolism",      3,
  "conven_lpa",               "Lipoprotein(a)",                 "Lipid metabolism",      4,
  "conven_hdlc_therapy_corr", "HDL cholesterol",                "Lipid metabolism",      5,
  "anthro_bmi",               "Body mass index",                "Adiposity",              6,
  "conven_hba1c",             "Glycated haemoglobin A1c",       "Glycaemic status",       7,
  "conven_glu",               "Blood glucose concentration",    "Glycaemic status",       8,
  "conven_alp",               "Alkaline phosphatase",           "Liver and metabolism",   9,
  "conven_aspat",             "Aspartate aminotransferase",     "Liver and metabolism",  10,
  "conven_ggt",               "Gamma-glutamyl transferase",     "Liver and metabolism",  11,
  "conven_tbili",             "Total bilirubin",                "Liver and metabolism",  12,
  "conven_urate",             "Uric acid",                      "Liver and metabolism",  13,
  "conven_egfr_creat",        "eGFR (creatinine)",              "Kidney function",       14,
  "conven_egfr_cysc",         "eGFR (cystatin C)",              "Kidney function",       15,
  "bp_mean_sbp",              "Systolic BP",                    "Blood pressure",        16,
  "bp_mean_dbp",              "Diastolic BP",                   "Blood pressure",        17,
  "conven_crp",               "C-reactive protein",             "Inflammation",          18,
  "blood_wbc",                "White blood cell count",         "Inflammation",          19,
  "blood_nlr",                "Neutrophil-to-lymphocyte ratio", "Inflammation",          20
)

BLOCK_ORDER_A <- c(
  "Lipid metabolism", "Adiposity", "Glycaemic status",
  "Liver and metabolism", "Kidney function", "Blood pressure", "Inflammation"
)

DEFAULT_B_OUTCOMES <- c(
  "anthro_wt",
  "dxa_total_fat",
  "anthro_hc",
  "anthro_wc",
  "anthro_whr",
  "anthro_whr_bmiadj",
  "anthro_whtr_bmiadj",
  "anthro_hc_bmiadj",
  "dxa_trunk_leg_fat_ratio",
  "dxa_ag_ratio",
  "dxa_gynoid_fat",
  "dxa_legs_fat",
  "dxa_legs_fat_frac",
  "amri_gfat_bmiadj",
  "amri_vat",
  "amri_vat_bmiadj",
  "amri_asat_bmiadj",
  "amri_liver_pdff_fr",
  "amri_liver_pdff_fr_bmiadj"
)

DEFAULT_B_LABELS <- c(
  anthro_wt = "Weight",
  dxa_total_fat = "Total fat mass",
  anthro_hc = "Hip circumference",
  anthro_wc = "Waist circumference",
  anthro_whr = "Waist-to-hip ratio",
  anthro_whr_bmiadj = "Waist-to-hip ratio adjusted for BMI",
  anthro_whtr_bmiadj = "Waist-to-height ratio adjusted for BMI",
  anthro_hc_bmiadj = "Hip circumference adjusted for BMI",
  dxa_trunk_leg_fat_ratio = "Trunk-to-leg fat ratio",
  dxa_ag_ratio = "Android-to-gynoid fat ratio",
  dxa_gynoid_fat = "Gynoid fat mass",
  dxa_legs_fat = "Legs fat mass",
  dxa_legs_fat_frac = "Legs fat fraction",
  amri_gfat_bmiadj = "Gluteofemoral fat adjusted for BMI",
  amri_vat = "Visceral adipose tissue volume",
  amri_vat_bmiadj = "Visceral adipose tissue adjusted for BMI",
  amri_asat_bmiadj = "Abdominal subcutaneous adipose tissue adjusted for BMI",
  amri_liver_pdff_fr = "Liver PDFF (FR mean)",
  amri_liver_pdff_fr_bmiadj = "Liver PDFF adjusted for BMI"
)

DEFAULT_C_OUTCOMES <- c(
  "cmr_lvmass",
  "cmr_lvmvr",
  "cmr_lvedv",
  "cmr_rvedv",
  "cmr_lvef",
  "olink_NTproBNP"
)

DEFAULT_C_LABELS <- c(
  cmr_lvmass = "Left ventricular mass",
  cmr_lvmvr = "Left ventricular mass-to-volume ratio",
  cmr_lvedv = "Left ventricular end-diastolic volume",
  cmr_rvedv = "Right ventricular end-diastolic volume",
  cmr_lvef = "Left ventricular ejection fraction",
  olink_NTproBNP = "NT-proBNP"
)

DEFAULT_D_GROUPS <- c(
  olink_ABO = "ABO-VWF axis",
  olink_VWF = "ABO-VWF axis",
  olink_F2 = "Coagulation cascade",
  olink_F7 = "Coagulation cascade",
  olink_F11 = "Coagulation cascade",
  olink_VCAM1 = "Endothelial activation",
  olink_ICAM1 = "Endothelial activation",
  olink_SELE = "Endothelial activation",
  olink_NOS3 = "Endothelial activation",
  olink_NCR1 = "Lymphocyte immunity",
  olink_FASLG = "Lymphocyte immunity",
  olink_CD244 = "Lymphocyte immunity",
  olink_LTA = "Lymphocyte immunity"
)

DEFAULT_D_LABELS <- c(
  olink_ABO = "Histo-blood group ABO system transferase",
  olink_VWF = "von Willebrand factor",
  olink_F2 = "Prothrombin",
  olink_F7 = "Coagulation factor VII",
  olink_F11 = "Coagulation factor XI",
  olink_VCAM1 = "Vascular cell adhesion protein 1",
  olink_ICAM1 = "Intercellular adhesion molecule 1",
  olink_SELE = "E-selectin",
  olink_NOS3 = "Endothelial nitric oxide synthase",
  olink_NCR1 = "Natural cytotoxicity triggering receptor 1",
  olink_FASLG = "Fas ligand",
  olink_CD244 = "Natural killer cell receptor 2B4",
  olink_LTA = "Lymphotoxin-alpha"
)

# ==============================================================================
# Input preparation
# ==============================================================================

coerce_association_data <- function(dat) {
  required <- c("outcome", "predictor", "beta", "p_value")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop_user(paste("Input data is missing:", paste(missing, collapse = ", ")))

  optional_defaults <- list(
    study = "UKBB", outcome_long = NA_character_, outcome_category = NA_character_,
    predictor_long = NA_character_, model = "", exclusion = "", n = NA_real_,
    se = NA_real_, lb = NA_real_, ub = NA_real_, p_fdr = NA_real_
  )
  for (nm in names(optional_defaults)) {
    if (!nm %in% names(dat)) dat[[nm]] <- optional_defaults[[nm]]
  }

  dat %>%
    mutate(
      across(where(is.character), trimws),
      predictor_original = predictor,
      predictor = resolve_predictor_column(predictor),
      across(c(n, beta, se, lb, ub, p_value, p_fdr), ~ suppressWarnings(as.numeric(.x))),
      lb = if_else(is.na(lb) & is.finite(beta) & is.finite(se), beta - 1.96 * se, lb),
      ub = if_else(is.na(ub) & is.finite(beta) & is.finite(se), beta + 1.96 * se, ub),
      lb_plot = pmin(lb, ub, beta, na.rm = TRUE),
      ub_plot = pmax(lb, ub, beta, na.rm = TRUE),
      lb_plot = if_else(is.infinite(lb_plot), NA_real_, lb_plot),
      ub_plot = if_else(is.infinite(ub_plot), NA_real_, ub_plot)
    )
}

filter_model_rows <- function(dat, study_filter = "", model_filter = "", exclusion_filter = "") {
  out <- dat
  if (nzchar(trimws(study_filter))) out <- out %>% filter(study == trimws(study_filter))
  if (nzchar(trimws(model_filter))) out <- out %>% filter(model == trimws(model_filter))
  if (nzchar(trimws(exclusion_filter))) out <- out %>% filter(exclusion == trimws(exclusion_filter))
  out
}

# Prefer a single original variable when several aliases map to the same canonical PRS.
# Total PRS aliases are collapsed to the canonical prs_total label.
# prs_total_cluster is the preferred/default source variable for Total PRS.
PREDICTOR_SOURCE_PRIORITY <- c(
  prs_total_cluster = 1,
  total_196_prs_std = 2,
  total_prs_std = 3,
  prs_total = 4,
  cardiac_prs_std = 1,
  prs_cardiac_remodeling = 2,
  endothelial_prs_std = 1,
  prs_endothelial = 2
)

resolve_alias_duplicates <- function(dat, panel_name) {
  if (!nrow(dat)) return(dat)

  split_rows <- split(dat, interaction(dat$outcome, dat$predictor, drop = TRUE, lex.order = TRUE))

  resolved <- lapply(split_rows, function(x) {
    # Remove fully identical duplicated rows first.
    x <- distinct(x)
    if (nrow(x) <= 1) return(x)

    original_names <- unique(x$predictor_original)

    # Only auto-resolve when duplication is caused by different source-variable
    # aliases that all map to the same canonical predictor.
    if (length(original_names) > 1 && all(table(x$predictor_original) == 1)) {
      priority <- unname(PREDICTOR_SOURCE_PRIORITY[x$predictor_original])
      priority[is.na(priority)] <- 1000
      best <- which(priority == min(priority))

      if (length(best) == 1) return(x[best, , drop = FALSE])
    }

    x
  })

  out <- bind_rows(resolved)
  assert_unique_pairs(out, panel_name)
  out
}

assert_unique_pairs <- function(dat, panel_name) {
  duplicates <- dat %>% count(outcome, predictor, name = "n_rows") %>% filter(n_rows > 1)
  if (nrow(duplicates)) {
    example_pairs <- duplicates %>%
      slice_head(n = 5) %>%
      transmute(pair = paste0(outcome, " + ", predictor, " (", n_rows, " rows)")) %>%
      pull(pair)

    stop_user(
      paste0(
        panel_name,
        " still has duplicate outcome-PRS combinations after filtering. Examples: ",
        paste(example_pairs, collapse = "; "),
        ". Check whether the input contains repeated rows within the selected study/model/exclusion set."
      )
    )
  }
  invisible(dat)
}

make_global_outcome_labels <- function(dat) {
  dat %>%
    distinct(outcome, outcome_long) %>%
    group_by(outcome) %>%
    summarise(label = {
      x <- outcome_long[!is.na(outcome_long) & nzchar(outcome_long)]
      if (length(x)) first(x) else first(outcome)
    }, .groups = "drop") %>%
    {setNames(.$label, .$outcome)}
}

# ==============================================================================
# Panel A: heatmap
# ==============================================================================

plot_panel_a <- function(dat, significance_method = "Bonferroni", fdr_alpha = 0.05,
                         bonferroni_alpha = 0.05, bonferroni_n_prs = 9,
                         bonferroni_n_biomarkers = 84, block_spacing = 0.45,
                         total_prs_gap = 0.55, tile_width = 0.92, tile_height = 0.92,
                         tile_text_size = 3.1, white_label_threshold = 0.20,
                         axis_label_size = 10, legend_title_size = 10, legend_text_size = 10,
                         beta_colour_limit = 0.30, legend_position = "right", colour_ticks = c(-0.30, -0.10, -0.01, 0, 0.01, 0.10, 0.30),
                         panel_border_colour = "grey45", panel_border_width = 0.35, margin_left = 5.5,
                         margin_right = 5.5, margin_top = 5.5, margin_bottom = 5.5) {
  if (!significance_method %in% c("FDR", "Bonferroni")) stop_user("Panel A significance method must be FDR or Bonferroni.")

  source <- dat %>%
    filter(outcome %in% BIOMARKER_MAP_A$outcome, predictor %in% PRS_KEEP_A) %>%
    # For Total PRS, use prs_total_cluster as the sole/default source variable.
    # Other source aliases that canonicalize to prs_total are excluded here to
    # prevent duplicate outcome-PRS combinations in Panel A.
    filter(predictor != "prs_total" | predictor_original == "prs_total_cluster")
  assert_unique_pairs(source, "Panel A")

  threshold <- bonferroni_alpha / (bonferroni_n_prs * bonferroni_n_biomarkers)

  # Numeric y positions preserve an adjustable gap before the Total PRS row.
  prs_positions <- c(
    setNames(rev(seq_len(8)), PRS_KEEP_A[1:8]),
    prs_total = 0 - total_prs_gap
  )

  heatmap <- crossing(outcome = BIOMARKER_MAP_A$outcome, predictor = PRS_KEEP_A) %>%
    left_join(source %>% select(outcome, predictor, beta, p_value, p_fdr), by = c("outcome", "predictor")) %>%
    left_join(BIOMARKER_MAP_A, by = "outcome") %>%
    mutate(
      significant = if (significance_method == "FDR") {
        !is.na(p_fdr) & p_fdr < fdr_alpha
      } else {
        !is.na(p_value) & p_value < threshold
      },
      beta_sqrt = sign(beta) * sqrt(abs(beta)),
      plotted_beta = if_else(significant, beta_sqrt, NA_real_),
      beta_label = case_when(
        !significant | is.na(beta) ~ "",
        TRUE ~ sub("^(-?)0\\.", "\\1.", sprintf("%.2f", beta))
      ),
      beta_label_colour = case_when(
        !significant | is.na(beta) ~ NA_character_,
        abs(beta) > white_label_threshold ~ "white",
        TRUE ~ "black"
      ),
      prs_position = unname(prs_positions[predictor]),
      block = factor(block, levels = BLOCK_ORDER_A),
      biomarker_label = factor(
        biomarker_label,
        levels = BIOMARKER_MAP_A %>% arrange(biomarker_order) %>% pull(biomarker_label)
      )
    )

  legend_beta_breaks <- colour_ticks[colour_ticks >= -beta_colour_limit & colour_ticks <= beta_colour_limit]
  if (!length(legend_beta_breaks)) legend_beta_breaks <- waiver()
  legend_breaks <- if (inherits(legend_beta_breaks, "waiver")) waiver() else sign(legend_beta_breaks) * sqrt(abs(legend_beta_breaks))
  legend_labels <- if (inherits(legend_beta_breaks, "waiver")) waiver() else format_beta_ticks(legend_beta_breaks)
  transformed_limit <- sqrt(beta_colour_limit)

  ggplot(heatmap, aes(x = biomarker_label, y = prs_position, fill = plotted_beta)) +
    geom_tile(width = tile_width, height = tile_height, colour = "white", linewidth = 0.35) +
    geom_text(
      aes(label = beta_label, colour = beta_label_colour),
      family = "Arial", size = tile_text_size, na.rm = TRUE, show.legend = FALSE
    ) +
    scale_colour_identity() +
    facet_grid(cols = vars(block), scales = "free_x", space = "free_x") +
    scale_y_continuous(
      breaks = unname(prs_positions[PRS_KEEP_A]),
      labels = unname(DEFAULT_PREDICTOR_LABELS[PRS_KEEP_A]),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-transformed_limit, transformed_limit), breaks = legend_breaks,
      labels = legend_labels, oob = scales::squish, na.value = "grey80",
      name = "Standardized\nbeta"
    ) +
    labs(x = NULL, y = NULL) +
    theme_linedraw(base_size = 11, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      axis.title = element_blank(),
      axis.text.x = element_text(size = axis_label_size, colour = "black", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = axis_label_size, colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = panel_border_colour, fill = NA, linewidth = panel_border_width),
      strip.background = element_blank(), strip.text = element_blank(),
      panel.spacing.x = unit(block_spacing, "lines"),
      legend.title = element_text(face = "bold", size = legend_title_size, hjust = 0.5),
      legend.text = element_text(size = legend_text_size),
      legend.position = legend_position,
      plot.margin = margin(margin_top, margin_right, margin_bottom, margin_left)
    )
}

# ==============================================================================
# Panel B: MUNW-versus-obesity scatter plot
# ==============================================================================

plot_panel_b <- function(dat, outcomes, outcome_labels, point_color = "#699760",
                         x_axis_label = "Beta for MUNW PRS",
                         y_axis_label = "Beta for obesity PRS",
                         x_limits = c(NA_real_, NA_real_), y_limits = c(NA_real_, NA_real_),
                         x_breaks = numeric(0), y_breaks = numeric(0),
                         legend_x_offset = 0, legend_y_offset = 0, legend_ncol = 1,
                         legend_title_hjust = 0.5,
                         point_size = 3.5, error_linewidth = 0.35,
                         label_size = 3.2, label_jitter_x = 0.005,
                         label_jitter_y = 0.005, label_offsets = NULL,
                         show_error_bars = TRUE,
                         show_labels = TRUE, gridlines = TRUE,
                         axis_text_size = 10, axis_title_size = 11,
                         legend_title_size = 9, legend_text_size = 9,
                         panel_border_colour = "grey45", panel_border_width = 0.35,
                         margin_left = 5.5, margin_right = 5.5,
                         margin_top = 5.5, margin_bottom = 5.5) {
  if (!length(outcomes)) stop_user("Panel B outcome list is empty.")

  source <- dat %>%
    filter(outcome %in% outcomes, predictor %in% c("prs_mdwo", "prs_obesity"))
  assert_unique_pairs(source, "Panel B")

  wide <- source %>%
    select(outcome, predictor, beta, lb_plot, ub_plot) %>%
    pivot_wider(names_from = predictor, values_from = c(beta, lb_plot, ub_plot), names_sep = "__")

  required <- c("beta__prs_mdwo", "beta__prs_obesity")
  missing <- setdiff(required, names(wide))
  if (length(missing)) stop_user("Panel B requires both MUNW and obesity PRS estimates for the selected outcomes.")

  if (is.null(label_offsets)) {
    label_offsets <- tibble::tibble(
      outcome = character(),
      label_dx = numeric(),
      label_dy = numeric()
    )
  }

  wide <- wide %>%
    left_join(label_offsets, by = "outcome") %>%
    mutate(
      label = label_lookup(outcome_labels, outcome),
      label_dx = coalesce(label_dx, label_jitter_x),
      label_dy = coalesce(label_dy, label_jitter_y),
      label_x = beta__prs_mdwo + label_dx,
      label_y = beta__prs_obesity + label_dy,
      modality = case_when(
        startsWith(outcome, "anthro_") ~ "Anthropometry",
        startsWith(outcome, "dxa_") ~ "DXA",
        startsWith(outcome, "amri_") ~ "Abdominal MRI",
        TRUE ~ "Other"
      ),
      modality = factor(modality, levels = c("Anthropometry", "DXA", "Abdominal MRI", "Other"))
    )

  modality_shapes <- c("Anthropometry" = 21, "DXA" = 22, "Abdominal MRI" = 24, "Other" = 23)

  p <- ggplot(wide, aes(x = beta__prs_mdwo, y = beta__prs_obesity))
  if (isTRUE(show_error_bars) && all(c("lb_plot__prs_mdwo", "ub_plot__prs_mdwo") %in% names(wide))) {
    p <- p + geom_segment(
      aes(x = lb_plot__prs_mdwo, xend = ub_plot__prs_mdwo,
          y = beta__prs_obesity, yend = beta__prs_obesity),
      linewidth = error_linewidth, color = point_color, na.rm = TRUE
    )
  }
  if (isTRUE(show_error_bars) && all(c("lb_plot__prs_obesity", "ub_plot__prs_obesity") %in% names(wide))) {
    p <- p + geom_segment(
      aes(x = beta__prs_mdwo, xend = beta__prs_mdwo,
          y = lb_plot__prs_obesity, yend = ub_plot__prs_obesity),
      linewidth = error_linewidth, color = point_color, na.rm = TRUE
    )
  }

  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey45") +
    geom_point(
      aes(shape = modality), fill = point_color, color = "black",
      size = point_size, stroke = 0.35
    )

  if (isTRUE(show_labels)) {
    p <- p + geom_text(aes(x = label_x, y = label_y, label = label), size = label_size, check_overlap = FALSE)
  }

  p +
    labs(x = x_axis_label, y = y_axis_label, shape = "Measurement modality") +
    scale_shape_manual(
      values = modality_shapes,
      drop = TRUE,
      guide = guide_legend(
        ncol = max(1, as.integer(legend_ncol)),
        byrow = TRUE,
        title.position = "top",
        title.hjust = legend_title_hjust
      )
    ) +
    scale_x_continuous(breaks = if (length(x_breaks)) x_breaks else waiver(), expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = if (length(y_breaks)) y_breaks else waiver(), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    theme_linedraw(base_size = 10, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(size = axis_text_size, color = "black"),
      axis.title.x = element_text(size = axis_title_size, color = "black"),
      axis.title.y = element_text(size = axis_title_size, color = "black", margin = margin(r = 3)),
      panel.grid.major = if (isTRUE(gridlines)) element_line(color = "grey92", linewidth = 0.25) else element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = panel_border_colour, linewidth = panel_border_width),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = element_text(size = legend_title_size, face = "bold", hjust = legend_title_hjust),
      legend.text = element_text(size = legend_text_size),
      legend.key.width = unit(1.2, "lines"),
      legend.margin = margin(t = 2 - legend_y_offset, r = -legend_x_offset, b = legend_y_offset, l = legend_x_offset),
      plot.margin = margin(margin_top, margin_right, margin_bottom, margin_left)
    )
}

# ==============================================================================
# Forest plot (displayed as Panel D)
# ==============================================================================

plot_panel_c <- function(dat, outcomes, predictors, outcome_labels,
                         predictor_labels, predictor_colors,
                         x_axis_label = "Beta estimate", x_limits = c(NA_real_, NA_real_),
                         x_breaks = numeric(0), legend_x_offset = 0, legend_y_offset = 0, legend_ncol = 1,
                         legend_title_hjust = 0.5,
                         point_size = 3.4, whisker_linewidth = 0.55,
                         zero_linewidth = 0.35, significance_threshold = 0.05,
                         outcome_text_size = 9, outcome_wrap_width = 32,
                         axis_text_size = 9,
                         axis_title_size = 10, legend_title_size = 9,
                         legend_text_size = 9, dodge_width = 0.65,
                         direct_label_size = 3.0, direct_label_hjust = -0.15,
                         alternating_fill = "grey96",
                         row_background_height = 1.00,
                         y_bottom_padding = 0.10,
                         y_top_padding = 0.10,
                         panel_border_colour = "grey45", panel_border_width = 0.35,
                         margin_left = 5.5, margin_right = 5.5,
                         margin_top = 5.5, margin_bottom = 5.5) {
  if (!length(outcomes)) stop_user("Panel C outcome list is empty.")
  if (!length(predictors)) stop_user("Panel C PRS list is empty.")

  source <- dat %>%
    filter(outcome %in% outcomes, predictor %in% predictors) %>%
    resolve_alias_duplicates("Panel C")
  if (!nrow(source)) stop_user("No data found for Panel C selections.")

  wrapped_outcome_labels <- stringr::str_wrap(
    label_lookup(outcome_labels, outcomes),
    width = max(1, as.integer(outcome_wrap_width))
  )

  source <- source %>%
    mutate(
      outcome_f = factor(outcome, levels = rev(outcomes), labels = rev(wrapped_outcome_labels)),
      predictor_f = factor(predictor, levels = predictors, labels = label_lookup(predictor_labels, predictors)),
      predictor_label = label_lookup(predictor_labels, predictor),
      color = color_lookup(predictor_colors, predictor),
      significant = !is.na(p_value) & p_value < significance_threshold,
      fill_color = if_else(significant, color, "white")
    )

  background <- tibble(
    outcome_f = factor(rev(wrapped_outcome_labels), levels = levels(source$outcome_f)),
    row_id = seq_along(outcomes),
    y_numeric = seq_along(outcomes)
  ) %>% filter(row_id %% 2 == 1)

  dodge <- position_dodge(width = dodge_width)


  predictor_colour_values <- setNames(
    color_lookup(predictor_colors, predictors),
    label_lookup(predictor_labels, predictors)
  )

  ggplot(source, aes(y = outcome_f, x = beta, group = predictor_f)) +
    geom_rect(
      data = background,
      aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = y_numeric - row_background_height / 2,
        ymax = y_numeric + row_background_height / 2
      ),
      inherit.aes = FALSE, fill = alternating_fill, color = NA
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = zero_linewidth) +
    geom_errorbarh(
      aes(xmin = lb_plot, xmax = ub_plot, color = predictor_f),
      height = 0, linewidth = whisker_linewidth, position = dodge, na.rm = TRUE
    ) +
    geom_point(
      aes(fill = fill_color, color = predictor_f), shape = 21, size = point_size,
      stroke = 0.3, position = dodge, na.rm = TRUE
    ) +
    scale_color_manual(
      values = predictor_colour_values,
      name = "Polygenic risk score",
      breaks = label_lookup(predictor_labels, predictors),
      drop = FALSE,
      guide = guide_legend(
        ncol = max(1, as.integer(legend_ncol)),
        byrow = TRUE,
        title.position = "top",
        title.hjust = legend_title_hjust,
        override.aes = list(
          shape = 21,
          size = point_size,
          fill = unname(predictor_colour_values),
          colour = unname(predictor_colour_values),
          stroke = 0.3,
          linewidth = whisker_linewidth
        )
      )
    ) +
    scale_fill_identity() +
    scale_x_continuous(breaks = if (length(x_breaks)) x_breaks else waiver(), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(
      expand = expansion(add = c(y_bottom_padding, y_top_padding))
    ) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(x = x_axis_label, y = NULL) +
    theme_linedraw(base_size = 10, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      axis.text.y = element_text(size = outcome_text_size, color = "black"),
      axis.text.x = element_text(size = axis_text_size, color = "black"),
      axis.title.x = element_text(size = axis_title_size),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_rect(fill = NA, color = panel_border_colour, linewidth = panel_border_width),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = element_text(size = legend_title_size, face = "bold", hjust = legend_title_hjust),
      legend.text = element_text(size = legend_text_size),
      legend.key.width = unit(1.2, "lines"),
      legend.margin = margin(t = 2 - legend_y_offset, r = -legend_x_offset, b = legend_y_offset, l = legend_x_offset),
      plot.margin = margin(margin_top, margin_right, margin_bottom, margin_left)
    )
}

# ============================================================================== 
# Protein heatmap (displayed as Panel C)
# ============================================================================== 

plot_panel_d <- function(dat, group_map, predictors, outcome_labels,
                         predictor_labels, predictor_colors,
                         significance_method = "Bonferroni", fdr_alpha = 0.05,
                         bonferroni_alpha = 0.05, bonferroni_n_tests = 2923,
                         group_spacing = 0.45, tile_width = 0.92,
                         tile_height = 0.92, tile_text_size = 3.0,
                         white_label_threshold = 0.20,
                         beta_colour_limit = 0.30,
                         protein_text_size = 9, prs_text_size = 9,
                         legend_title_size = 9, legend_text_size = 9,
                         colourbar_width_cm = 4.5,
                         label_wrap_width = 34,
                         legend_position = "bottom",
                         colour_ticks = c(-0.30, -0.10, -0.01, 0, 0.01, 0.10, 0.30),
                         legend_x_offset = 0, legend_y_offset = 0, legend_title_hjust = 0.5,
                         panel_border_colour = "grey45", panel_border_width = 0.35,
                         margin_left = 5.5, margin_right = 5.5,
                         margin_top = 5.5, margin_bottom = 5.5) {
  if (!length(group_map)) stop_user("Panel D protein/group list is empty.")
  if (!length(predictors)) stop_user("Panel D PRS list is empty.")
  if (!significance_method %in% c("FDR", "Bonferroni")) stop_user("Panel D significance method must be FDR or Bonferroni.")

  protein_map <- tibble(
    outcome = names(group_map),
    group = unname(group_map),
    protein_order = seq_along(group_map)
  )
  group_order <- unique(protein_map$group)

  source <- dat %>%
    filter(outcome %in% protein_map$outcome, predictor %in% predictors) %>%
    filter(predictor != "prs_total" | predictor_original == "prs_total_cluster") %>%
    mutate(.source_row = row_number()) %>%
    group_by(outcome, predictor) %>%
    arrange(desc(!is.na(beta)), desc(!is.na(p_value)), desc(!is.na(p_fdr)), .source_row, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-.source_row)

  assert_unique_pairs(source, "Panel D")
  if (!nrow(source)) stop_user("No data found for Panel D selections.")

  threshold <- bonferroni_alpha / bonferroni_n_tests

  heatmap <- crossing(outcome = protein_map$outcome, predictor = predictors) %>%
    left_join(source %>% select(outcome, predictor, beta, p_value, p_fdr), by = c("outcome", "predictor")) %>%
    left_join(protein_map, by = "outcome") %>%
    mutate(
      significant = if (significance_method == "FDR") {
        !is.na(p_fdr) & p_fdr < fdr_alpha
      } else {
        !is.na(p_value) & p_value < threshold
      },
      beta_sqrt = sign(beta) * sqrt(abs(beta)),
      plotted_beta = if_else(significant, beta_sqrt, NA_real_),
      beta_label = case_when(
        !significant | is.na(beta) ~ "",
        TRUE ~ sub("^(-?)0\\.", "\\1.", sprintf("%.2f", beta))
      ),
      beta_label_colour = case_when(
        !significant | is.na(beta) ~ NA_character_,
        abs(beta) > white_label_threshold ~ "white",
        TRUE ~ "black"
      ),
      predictor_f = factor(predictor, levels = predictors, labels = label_lookup(predictor_labels, predictors)),
      group = factor(group, levels = group_order),
      protein_key = factor(outcome, levels = rev(protein_map$outcome))
    )

  protein_axis_labels <- setNames(
    str_wrap(label_lookup(outcome_labels, protein_map$outcome), width = label_wrap_width),
    protein_map$outcome
  )

  transformed_colour_limit <- sqrt(beta_colour_limit)
  legend_beta_breaks <- colour_ticks[colour_ticks >= -beta_colour_limit & colour_ticks <= beta_colour_limit]
  if (!length(legend_beta_breaks)) legend_beta_breaks <- waiver()
  legend_transformed_breaks <- if (inherits(legend_beta_breaks, "waiver")) waiver() else sign(legend_beta_breaks) * sqrt(abs(legend_beta_breaks))
  legend_beta_labels <- if (inherits(legend_beta_breaks, "waiver")) waiver() else format_beta_ticks(legend_beta_breaks)

  ggplot(heatmap, aes(x = predictor_f, y = protein_key, fill = plotted_beta)) +
    geom_tile(width = tile_width, height = tile_height, colour = "white", linewidth = 0.35) +
    geom_text(
      data = filter(heatmap, beta_label_colour == "black"),
      aes(label = beta_label), family = "Arial", size = tile_text_size,
      colour = "black", na.rm = TRUE, show.legend = FALSE
    ) +
    geom_text(
      data = filter(heatmap, beta_label_colour == "white"),
      aes(label = beta_label), family = "Arial", size = tile_text_size,
      colour = "white", na.rm = TRUE, show.legend = FALSE
    ) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free_y") +
    # Match Panel A: retain a visible white perimeter around the outer tiles.
    # The tile width/height controls continue to determine the gutters between cells.
    scale_x_discrete(expand = expansion(add = c(0.60, 0.60))) +
    scale_y_discrete(
      labels = protein_axis_labels,
      expand = expansion(add = c(0.60, 0.60))
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-transformed_colour_limit, transformed_colour_limit),
      breaks = legend_transformed_breaks,
      labels = legend_beta_labels,
      oob = scales::squish,
      na.value = "grey80",
      name = "Standardized beta",
      guide = guide_colourbar(
        title.position = "top",
        title.hjust = legend_title_hjust,
        barwidth = unit(colourbar_width_cm, "cm"),
        barheight = unit(0.35, "cm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_linedraw(base_size = 10, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      axis.title = element_blank(),
      axis.text.x = element_text(size = prs_text_size, colour = "black", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = protein_text_size, colour = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = panel_border_colour, linewidth = panel_border_width),
      strip.background = element_blank(), strip.text = element_blank(),
      panel.spacing.y = unit(group_spacing, "lines"),
      legend.position = legend_position,
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.title = element_text(size = legend_title_size, face = "bold", hjust = legend_title_hjust),
      legend.text = element_text(size = legend_text_size),
      legend.key.width = unit(1.2, "lines"),
      legend.margin = margin(t = 2 - legend_y_offset, r = -legend_x_offset, b = legend_y_offset, l = legend_x_offset),
      plot.margin = margin(margin_top, margin_right, margin_bottom, margin_left)
    )
}

# ==============================================================================
# Shared panel wrappers and legends
# ==============================================================================

make_prs_legend <- function(predictors, predictor_labels, predictor_colors,
                            title = "PRS", point_size = 3, text_size = 9,
                            title_size = 9, nrow = 1, x_spacing = 1.35) {
  if (!length(predictors)) return(ggplot() + theme_void())
  legend_data <- tibble(
    predictor = predictors,
    label = label_lookup(predictor_labels, predictors),
    color = color_lookup(predictor_colors, predictors),
    item = seq_along(predictors)
  ) %>%
    mutate(row = ((item - 1) %% nrow) + 1, col = floor((item - 1) / nrow) + 1,
           x = col * x_spacing, y = -row)

  ggplot(legend_data, aes(x, y)) +
    annotate("text", x = 0.25, y = -1, label = title, hjust = 0,
             fontface = "bold", size = title_size / ggplot2::.pt) +
    geom_point(aes(fill = color, color = color), shape = 21, size = point_size, stroke = 0.3) +
    geom_text(aes(label = label), nudge_x = 0.18, hjust = 0, size = text_size / ggplot2::.pt) +
    scale_fill_identity() + scale_color_identity() +
    coord_cartesian(xlim = c(0.15, max(legend_data$x) + 1.1), ylim = c(-nrow - 0.5, 0), clip = "off") +
    theme_void()
}

extract_legend_grob <- function(plot) {
  gt <- ggplotGrob(plot)
  grob_names <- vapply(gt$grobs, function(x) x$name %||% "", character(1))
  idx <- which(grepl("^guide-box", grob_names))
  if (!length(idx)) return(nullGrob())
  gt$grobs[[idx[[1]]]]
}

wrap_panel <- function(plot, tag, tag_size = 16, tag_x_offset = 0, tag_y_offset = 0,
                       legend = NULL, legend_height = 0.10) {
  plot <- plot + labs(tag = tag) +
    theme(
      plot.tag = element_text(face = "bold", size = tag_size, family = "Arial"),
      plot.tag.position = c(0 + tag_x_offset, 1 + tag_y_offset)
    )
  if (is.null(legend)) return(plot)
  plot / legend + plot_layout(heights = c(1, legend_height))
}

free_left_side <- function(plot) {
  if ("free" %in% getNamespaceExports("patchwork")) {
    patchwork::free(plot, side = "l")
  } else {
    plot
  }
}

# ==============================================================================
# Example data
# ==============================================================================

make_example_data <- function() {
  set.seed(1804)
  outcomes <- unique(c(BIOMARKER_MAP_A$outcome, DEFAULT_B_OUTCOMES, DEFAULT_C_OUTCOMES, names(DEFAULT_D_GROUPS)))
  expand_grid(study = "UKBB", outcome = outcomes, predictor = ALL_PRS) %>%
    mutate(
      outcome_long = outcome,
      outcome_long = if_else(outcome %in% names(DEFAULT_C_LABELS), unname(DEFAULT_C_LABELS[outcome]), outcome_long),
      outcome_long = if_else(outcome %in% names(DEFAULT_D_LABELS), unname(DEFAULT_D_LABELS[outcome]), outcome_long),
      outcome_category = "Example",
      predictor_long = unname(DEFAULT_PREDICTOR_LABELS[predictor]),
      model = "Main model", exclusion = "Main analysis", n = 29515,
      signal = case_when(
        predictor == "prs_mdwo" & outcome %in% DEFAULT_B_OUTCOMES ~ runif(n(), -0.04, 0.14),
        predictor == "prs_obesity" & outcome %in% DEFAULT_B_OUTCOMES ~ runif(n(), -0.02, 0.18),
        predictor == "prs_cardiac_remodeling" & outcome %in% DEFAULT_C_OUTCOMES ~ runif(n(), 0.05, 0.18),
        predictor == "prs_endothelial" & outcome %in% DEFAULT_C_OUTCOMES ~ runif(n(), -0.03, 0.10),
        predictor == "prs_total" & outcome %in% DEFAULT_C_OUTCOMES ~ runif(n(), 0.03, 0.13),
        outcome %in% names(DEFAULT_D_GROUPS) & predictor %in% c("prs_coagulation", "prs_endothelial", "prs_total") ~ runif(n(), -0.12, 0.30),
        outcome %in% BIOMARKER_MAP_A$outcome ~ runif(n(), -0.25, 0.25),
        TRUE ~ rnorm(n(), 0, 0.03)
      ),
      beta = signal + rnorm(n(), 0, 0.015),
      se = runif(n(), 0.008, 0.025),
      lb = beta - 1.96 * se, ub = beta + 1.96 * se,
      p_value = 2 * pnorm(-abs(beta / se)),
      p_fdr = p.adjust(p_value, method = "fdr")
    ) %>%
    select(study, outcome, outcome_long, outcome_category, predictor, predictor_long,
           model, exclusion, n, beta, se, lb, ub, p_value, p_fdr)
}

# ==============================================================================
# Settings import/export
# ==============================================================================

SETTINGS_INPUT_TYPES <- c(
  study_filter = "text", model_filter = "text", exclusion_filter = "text",
  predictor_label_map = "text_area", predictor_color_map = "text_area",
  a_significance_method = "select", a_fdr_alpha = "numeric", a_bonf_alpha = "numeric",
  a_bonf_n_prs = "numeric", a_bonf_n_biomarkers = "numeric", a_block_spacing = "numeric",
  a_total_gap = "numeric", a_tile_width = "numeric", a_tile_height = "numeric",
  a_tile_text_size = "numeric", a_white_threshold = "numeric", a_axis_text_size = "numeric",
  a_legend_title_size = "numeric", a_legend_text_size = "numeric", a_beta_limit = "numeric", a_colour_ticks = "text", a_legend_position = "select",
  b_outcomes = "text_area", b_label_map = "text_area", b_x_label = "text", b_y_label = "text",
  b_xmin = "numeric", b_xmax = "numeric", b_ymin = "numeric", b_ymax = "numeric",
  b_x_ticks = "text", b_y_ticks = "text", b_legend_x = "numeric", b_legend_y = "numeric", b_legend_ncol = "numeric", b_legend_title_hjust = "numeric",
  b_point_color = "text", b_point_size = "numeric", b_error_width = "numeric",
  b_label_size = "numeric", b_label_jitter_x = "numeric", b_label_jitter_y = "numeric",
  b_show_error = "checkbox", b_show_labels = "checkbox", b_gridlines = "checkbox",
  b_axis_text_size = "numeric", b_axis_title_size = "numeric",
  b_legend_title_size = "numeric", b_legend_text_size = "numeric",
  c_outcomes = "text_area", c_predictors = "text_area", c_label_map = "text_area",
  c_x_label = "text", c_xmin = "numeric", c_xmax = "numeric", c_x_ticks = "text",
  c_legend_x = "numeric", c_legend_y = "numeric", c_legend_ncol = "numeric", c_legend_title_hjust = "numeric", c_point_size = "numeric",
  c_whisker_width = "numeric", c_p_threshold = "numeric", c_dodge_width = "numeric",
  c_outcome_text_size = "numeric", c_label_wrap = "numeric", c_axis_text_size = "numeric", c_axis_title_size = "numeric",
  c_legend_title_size = "numeric", c_legend_text_size = "numeric",
  d_groups = "text_area", d_predictors = "text_area", d_label_map = "text_area",
  d_significance_method = "select", d_fdr_alpha = "numeric", d_bonf_alpha = "numeric",
  d_bonf_n_tests = "numeric", d_group_gap = "numeric", d_label_wrap = "numeric",
  d_tile_width = "numeric", d_tile_height = "numeric", d_tile_text_size = "numeric",
  d_white_threshold = "numeric", d_beta_limit = "numeric", d_colour_ticks = "text", d_protein_text_size = "numeric",
  d_prs_text_size = "numeric", d_legend_title_size = "numeric", d_legend_text_size = "numeric",
  d_colourbar_width = "numeric", d_legend_x = "numeric", d_legend_y = "numeric", d_legend_title_hjust = "numeric", d_legend_position = "select",
  panel_border_colour = "text", panel_border_width = "numeric",
  panel_tag_size = "numeric",
  panel_a_tag_x = "numeric", panel_a_tag_y = "numeric",
  panel_b_tag_x = "numeric", panel_b_tag_y = "numeric",
  panel_c_tag_x = "numeric", panel_c_tag_y = "numeric",
  panel_d_tag_x = "numeric", panel_d_tag_y = "numeric",
  col_b_width = "numeric", col_c_width = "numeric", col_d_width = "numeric",
  c_legend_height = "numeric", row_top = "numeric", row_bottom = "numeric",
  margin_a_left = "numeric", margin_b_left = "numeric", margin_c_left = "numeric",
  margin_d_left = "numeric", margin_right = "numeric", margin_d_right = "numeric", margin_top = "numeric", margin_bottom = "numeric",
  png_width = "numeric", png_height = "numeric", png_dpi = "numeric"
)

collect_settings <- function(input) {
  vals <- setNames(vector("list", length(SETTINGS_INPUT_TYPES)), names(SETTINGS_INPUT_TYPES))
  for (id in names(SETTINGS_INPUT_TYPES)) vals[id] <- list(input[[id]])
  list(app = "four_panel_biomarker_prs_app", settings_version = 6,
       exported_at = as.character(Sys.time()), settings = vals)
}

apply_settings <- function(settings, session) {
  for (id in intersect(names(settings), names(SETTINGS_INPUT_TYPES))) {
    type <- SETTINGS_INPUT_TYPES[[id]]
    value <- settings[[id]]
    if (type == "numeric") {
      updateNumericInput(session, id, value = if (is.null(value)) NA_real_ else as.numeric(value[[1]]))
    } else if (type == "checkbox") {
      updateCheckboxInput(session, id, value = isTRUE(value[[1]]))
    } else if (type == "select") {
      updateSelectInput(session, id, selected = as.character(value[[1]]))
    } else if (type == "text_area") {
      updateTextAreaInput(session, id, value = if (is.null(value)) "" else paste(unlist(value), collapse = "\n"))
    } else {
      updateTextInput(session, id, value = if (is.null(value)) "" else paste(unlist(value), collapse = "\n"))
    }
  }
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  titlePanel("Four-panel biomarker / process-specific PRS figure builder"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tabsetPanel(
        tabPanel(
          "Data and labels",
          fileInput("assoc_file", "Upload association CSV or CSV.GZ", accept = c(".csv", ".gz")),
          checkboxInput("use_example", "Use example data", TRUE),
          textInput("study_filter", "Study filter; blank = all", "UKBB"),
          textInput("model_filter", "Model filter; blank = all", ""),
          textInput("exclusion_filter", "Exclusion filter; blank = all", ""),
          tags$hr(),
          fileInput("settings_file", "Import settings JSON", accept = ".json"),
          actionButton("load_settings", "Apply settings"),
          downloadButton("download_settings", "Export settings"),
          tags$hr(),
          textAreaInput(
            "predictor_label_map", "PRS labels: predictor = label",
            paste(names(DEFAULT_PREDICTOR_LABELS), DEFAULT_PREDICTOR_LABELS, sep = " = ", collapse = "\n"), rows = 9
          ),
          textAreaInput(
            "predictor_color_map", "PRS colours: predictor = colour",
            paste(names(DEFAULT_PREDICTOR_COLORS), DEFAULT_PREDICTOR_COLORS, sep = " = ", collapse = "\n"), rows = 9
          )
        ),
        tabPanel(
          "Panel A",
          helpText(
            "Biomarker rows are selected using the outcome column. Choose one study, model, and exclusion set so that each outcome-PRS combination is unique."
          ),
          selectInput("a_study_filter", "Study", choices = character(0)),
          selectInput("a_model_filter", "Model", choices = character(0)),
          selectInput("a_exclusion_filter", "Exclusion set", choices = character(0)),
          selectInput("a_significance_method", "Significance method", c("Bonferroni", "FDR"), "Bonferroni"),
          fluidRow(column(6, numericInput("a_fdr_alpha", "FDR alpha", 0.05, min = 0, max = 1)),
                   column(6, numericInput("a_bonf_alpha", "Bonferroni alpha", 0.05, min = 0, max = 1))),
          fluidRow(column(6, numericInput("a_bonf_n_prs", "Number of PRSs", 9, min = 1)),
                   column(6, numericInput("a_bonf_n_biomarkers", "Number of biomarkers", 84, min = 1))),
          fluidRow(column(6, numericInput("a_block_spacing", "Block spacing", 0.45, step = 0.05)),
                   column(6, numericInput("a_total_gap", "Gap before Total PRS", 0.55, step = 0.05))),
          fluidRow(column(6, numericInput("a_tile_width", "Tile width", 0.92, step = 0.02)),
                   column(6, numericInput("a_tile_height", "Tile height", 0.92, step = 0.02))),
          h4("Typography"),
          fluidRow(column(6, numericInput("a_tile_text_size", "Coefficient text size", 3.1, step = 0.1)),
                   column(6, numericInput("a_axis_text_size", "Axis text size", 10, step = 0.5))),
          fluidRow(column(6, numericInput("a_legend_title_size", "Legend title size", 10, step = 0.5)),
                   column(6, numericInput("a_legend_text_size", "Legend tick text size", 10, step = 0.5))),
          fluidRow(column(6, numericInput("a_white_threshold", "White-label |beta| threshold", 0.20, step = 0.01)),
                   column(6, numericInput("a_beta_limit", "Colour-scale |beta| limit", 0.30, step = 0.01))),
          textInput("a_colour_ticks", "Colour-bar ticks, comma separated", "-.30, -.10, -.01, 0, .01, .10, .30"),
          selectInput("a_legend_position", "Legend position", c("right", "bottom", "left", "top"), "right")
        ),
        tabPanel(
          "Panel B",
          textAreaInput("b_outcomes", "Outcomes, in plotting order", paste(DEFAULT_B_OUTCOMES, collapse = "\n"), rows = 12),
          textAreaInput(
            "b_label_map", "Outcome labels: outcome = label",
            paste(names(DEFAULT_B_LABELS), DEFAULT_B_LABELS, sep = " = ", collapse = "\n"),
            rows = 12
          ),
          textInput("b_x_label", "X-axis label", "Beta for MUNW PRS"),
          textInput("b_y_label", "Y-axis label", "Beta for obesity PRS"),
          fluidRow(column(3, numericInput("b_xmin", "X min", NA)), column(3, numericInput("b_xmax", "X max", NA)),
                   column(3, numericInput("b_ymin", "Y min", NA)), column(3, numericInput("b_ymax", "Y max", NA))),
          fluidRow(column(6, textInput("b_x_ticks", "X-axis ticks", "")),
                   column(6, textInput("b_y_ticks", "Y-axis ticks", ""))),
          fluidRow(column(4, numericInput("b_legend_x", "Legend x offset (mm)", 0, step = 1)),
                   column(4, numericInput("b_legend_y", "Legend y offset (mm)", 0, step = 1)),
                   column(4, numericInput("b_legend_ncol", "Legend columns", 1, min = 1, max = 4, step = 1))),
          numericInput("b_legend_title_hjust", "Legend title horizontal alignment (0-1)", 0.5, min = 0, max = 1, step = 0.05),
          fluidRow(column(4, textInput("b_point_color", "Point colour", "#699760")),
                   column(4, numericInput("b_point_size", "Point size", 3.5, step = 0.1)),
                   column(4, numericInput("b_error_width", "CI line width", 0.35, step = 0.05))),
          fluidRow(column(4, numericInput("b_label_size", "Label size", 3.2, step = 0.1)),
                   column(4, numericInput("b_label_jitter_x", "Default label x offset", 0.005, step = 0.005)),
                   column(4, numericInput("b_label_jitter_y", "Default label y offset", 0.005, step = 0.005))),
          textAreaInput(
            "b_label_offsets",
            "Per-outcome label offsets: outcome = x offset, y offset",
            value = "",
            rows = 7,
            placeholder = paste(
              "anthro_wt = 0.010, 0.005",
              "dxa_total_fat = -0.015, 0.010",
              sep = "\n"
            )
          ),
          helpText("Specific offsets override the defaults. Leave an offset blank or non-numeric to use the default for that direction."),
          checkboxInput("b_show_error", "Show horizontal and vertical 95% CIs", TRUE),
          checkboxInput("b_show_labels", "Show outcome labels", TRUE),
          checkboxInput("b_gridlines", "Show major gridlines", TRUE),
          h4("Typography"),
          fluidRow(column(6, numericInput("b_axis_text_size", "Axis tick size", 10, step = 0.5)),
                   column(6, numericInput("b_axis_title_size", "Axis title size", 11, step = 0.5))),
          fluidRow(column(6, numericInput("b_legend_title_size", "Legend title size", 9, step = 0.5)),
                   column(6, numericInput("b_legend_text_size", "Legend text size", 9, step = 0.5)))
        ),
        tabPanel(
          "Panel D",
          textAreaInput("c_outcomes", "Cardiac traits, in plotting order", paste(DEFAULT_C_OUTCOMES, collapse = "\n"), rows = 8),
          textAreaInput("c_predictors", "PRSs", paste(c("prs_cardiac_remodeling", "prs_endothelial", "prs_total"), collapse = "\n"), rows = 4),
          textAreaInput("c_label_map", "Outcome labels: outcome = label", paste(names(DEFAULT_C_LABELS), DEFAULT_C_LABELS, sep = " = ", collapse = "\n"), rows = 8),
          textInput("c_x_label", "X-axis label", "Beta estimate"),
          fluidRow(column(6, numericInput("c_xmin", "X min", NA)), column(6, numericInput("c_xmax", "X max", NA))),
          textInput("c_x_ticks", "X-axis ticks, comma separated", ""),
          fluidRow(column(4, numericInput("c_legend_x", "Legend x offset (mm)", 0, step = 1)),
                   column(4, numericInput("c_legend_y", "Legend y offset (mm)", 0, step = 1)),
                   column(4, numericInput("c_legend_ncol", "Legend columns", 1, min = 1, max = 4, step = 1))),
          numericInput("c_legend_title_hjust", "Legend title horizontal alignment (0-1)", 0.5, min = 0, max = 1, step = 0.05),
          fluidRow(column(4, numericInput("c_point_size", "Point size", 3.4, step = 0.1)),
                   column(4, numericInput("c_whisker_width", "CI line width", 0.55, step = 0.05)),
                   column(4, numericInput("c_dodge_width", "PRS separation", 0.65, step = 0.05))),
          fluidRow(
            column(4, numericInput("c_row_background_height", "Alternating box height", 1.00, min = 0.1, max = 1.5, step = 0.02)),
            column(4, numericInput("c_y_bottom_padding", "Space below last outcome", 0.10, min = 0, step = 0.05)),
            column(4, numericInput("c_y_top_padding", "Space above first outcome", 0.10, min = 0, step = 0.05))
          ),
          numericInput("c_p_threshold", "Filled-head p-value threshold", 0.05, min = 0, max = 1, step = 0.001),
          h4("Typography"),
          fluidRow(column(4, numericInput("c_outcome_text_size", "Outcome text size", 9, step = 0.5)),
                   column(4, numericInput("c_label_wrap", "Outcome wrap width", 32, min = 5, step = 1)),
                   column(4, numericInput("c_axis_text_size", "X-axis tick size", 9, step = 0.5))),
          fluidRow(column(6, numericInput("c_axis_title_size", "X-axis title size", 10, step = 0.5)),
                   column(6, numericInput("c_legend_title_size", "Legend title size", 9, step = 0.5))),
          numericInput("c_legend_text_size", "Legend text size", 9, step = 0.5)
        ),
        tabPanel(
          "Panel C",
          textAreaInput("d_groups", "Protein panels: outcome = panel", paste(names(DEFAULT_D_GROUPS), DEFAULT_D_GROUPS, sep = " = ", collapse = "\n"), rows = 14),
          textAreaInput("d_predictors", "PRSs", paste(c("prs_coagulation", "prs_endothelial", "prs_total"), collapse = "\n"), rows = 4),
          textAreaInput("d_label_map", "Protein labels: outcome = full protein name", paste(names(DEFAULT_D_LABELS), DEFAULT_D_LABELS, sep = " = ", collapse = "\n"), rows = 14),
          selectInput("d_significance_method", "Significance method", c("Bonferroni", "FDR"), "Bonferroni"),
          fluidRow(column(4, numericInput("d_fdr_alpha", "FDR alpha", 0.05)),
                   column(4, numericInput("d_bonf_alpha", "Bonferroni alpha", 0.05)),
                   column(4, numericInput("d_bonf_n_tests", "Number of tests", 2923, min = 1))),
          fluidRow(column(4, numericInput("d_tile_width", "Tile width", 0.92, step = 0.02)),
                   column(4, numericInput("d_tile_height", "Tile height", 0.92, step = 0.02)),
                   column(4, numericInput("d_tile_text_size", "Coefficient text size", 3.0, step = 0.1))),
          fluidRow(column(4, numericInput("d_group_gap", "Panel gap", 0.45, step = 0.05)),
                   column(4, numericInput("d_label_wrap", "Protein wrap width", 34, min = 5)),
                   column(4, numericInput("d_protein_text_size", "Protein text size", 9, step = 0.5))),
          h4("Typography"),
          fluidRow(column(4, numericInput("d_prs_text_size", "PRS text size", 9, step = 0.5)),
                   column(4, numericInput("d_legend_title_size", "Legend title size", 9, step = 0.5)),
                   column(4, numericInput("d_legend_text_size", "Legend tick text size", 9, step = 0.5))),
          numericInput("d_colourbar_width", "Colour-bar width (cm)", 4.5, min = 0.5, step = 0.25),
          fluidRow(column(6, numericInput("d_white_threshold", "White-label |beta| threshold", 0.20, step = 0.01)),
                   column(6, numericInput("d_beta_limit", "Colour-scale |beta| limit", 0.30, step = 0.01))),
          textInput("d_colour_ticks", "Colour-bar ticks, comma separated", "-.30, -.10, -.01, 0, .01, .10, .30"),
          fluidRow(column(6, numericInput("d_legend_x", "Legend x offset (mm)", 0, step = 1)),
                   column(6, numericInput("d_legend_y", "Legend y offset (mm)", 0, step = 1))),
          numericInput("d_legend_title_hjust", "Legend title horizontal alignment (0-1)", 0.5, min = 0, max = 1, step = 0.05),
          selectInput("d_legend_position", "Legend position", c("bottom", "right", "left", "top"), "bottom")
        ),
        tabPanel(
          "Layout and export",
          fluidRow(column(3, numericInput("panel_tag_size", "Panel-letter size", 16, step = 0.5)),
                   column(3, numericInput("col_b_width", "Panel B relative width", 1, min = 0.1, step = 0.05)),
                   column(3, numericInput("col_c_width", "Panel C relative width", 1, min = 0.1, step = 0.05)),
                   column(3, numericInput("col_d_width", "Panel D relative width", 1, min = 0.1, step = 0.05))),
          h4("Shared plotting-area border"),
          fluidRow(column(6, textInput("panel_border_colour", "Border colour", "grey45")),
                   column(6, numericInput("panel_border_width", "Border width", 0.35, min = 0, step = 0.05))),
          h4("Panel-label offsets"),
          helpText("Offsets are relative to the default upper-left position. Positive x moves right; positive y moves up."),
          fluidRow(column(3, numericInput("panel_a_tag_x", "A x offset", 0, step = 0.01)),
                   column(3, numericInput("panel_a_tag_y", "A y offset", 0, step = 0.01)),
                   column(3, numericInput("panel_b_tag_x", "B x offset", 0, step = 0.01)),
                   column(3, numericInput("panel_b_tag_y", "B y offset", 0, step = 0.01))),
          fluidRow(column(3, numericInput("panel_c_tag_x", "C x offset", 0, step = 0.01)),
                   column(3, numericInput("panel_c_tag_y", "C y offset", 0, step = 0.01)),
                   column(3, numericInput("panel_d_tag_x", "D x offset", 0, step = 0.01)),
                   column(3, numericInput("panel_d_tag_y", "D y offset", 0, step = 0.01))),
          fluidRow(column(4, numericInput("row_top", "Top-row height", 1, step = 0.05)),
                   column(4, numericInput("row_bottom", "Bottom-row height", 1, step = 0.05)),
                   column(4, numericInput("c_legend_height", "Shared lower legend height", 0.22, min = 0.08, max = 0.50, step = 0.01))),
          h4("Panel left margins"),
          fluidRow(column(3, numericInput("margin_a_left", "A", 5.5, step = 1)),
                   column(3, numericInput("margin_b_left", "B", 5.5, step = 1)),
                   column(3, numericInput("margin_c_left", "C", 5.5, step = 1)),
                   column(3, numericInput("margin_d_left", "D", 5.5, step = 1))),
          fluidRow(column(3, numericInput("margin_right", "Shared right margin", 5.5, step = 1)),
                   column(3, numericInput("margin_d_right", "Panel D right margin", 0, step = 1)),
                   column(3, numericInput("margin_top", "Top margin", 5.5, step = 1)),
                   column(3, numericInput("margin_bottom", "Bottom margin", 5.5, step = 1))),
          tags$hr(),
          fluidRow(column(4, numericInput("png_width", "PNG width (in)", 16, step = 0.5)),
                   column(4, numericInput("png_height", "PNG height (in)", 14, step = 0.5)),
                   column(4, numericInput("png_dpi", "PNG DPI", 600, min = 72))),
          downloadButton("download_png", "Download PNG"),
          downloadButton("download_pdf", "Download PDF")
        )
      )
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Plot", plotOutput("combined_plot", height = "1250px")),
        tabPanel("Data preview", DTOutput("data_table")),
        tabPanel("Available names", verbatimTextOutput("available_names"))
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {
  raw_data <- reactive({
    if (isTRUE(input$use_example)) return(coerce_association_data(make_example_data()))
    if (is.null(input$assoc_file)) stop_user("Upload an association file or enable example data.")
    coerce_association_data(safe_read_csv(input$assoc_file$datapath))
  })

  assoc_data <- reactive({
    filter_model_rows(raw_data(), input$study_filter, input$model_filter, input$exclusion_filter)
  })

  # Panel A uses its own exact analysis-stratum selectors. Biomarkers themselves
  # are always identified by the outcome column inside plot_panel_a().
  observeEvent(raw_data(), {
    dat <- raw_data()

    study_choices <- sort(unique(dat$study[!is.na(dat$study)]))
    model_choices <- sort(unique(dat$model[!is.na(dat$model)]))
    exclusion_choices <- sort(unique(dat$exclusion[!is.na(dat$exclusion)]))

    labelled_choices <- function(x) {
      labels <- ifelse(nzchar(x), x, "(blank)")
      stats::setNames(x, labels)
    }

    current_study <- isolate(input$a_study_filter)
    current_model <- isolate(input$a_model_filter)
    current_exclusion <- isolate(input$a_exclusion_filter)

    default_study <- if ("UKBB" %in% study_choices) "UKBB" else if (length(study_choices)) study_choices[[1]] else ""
    default_model <- if (length(model_choices)) model_choices[[1]] else ""
    default_exclusion <- if (length(exclusion_choices)) exclusion_choices[[1]] else ""

    updateSelectInput(
      session, "a_study_filter",
      choices = labelled_choices(study_choices),
      selected = if (!is.null(current_study) && current_study %in% study_choices) current_study else default_study
    )
    updateSelectInput(
      session, "a_model_filter",
      choices = labelled_choices(model_choices),
      selected = if (!is.null(current_model) && current_model %in% model_choices) current_model else default_model
    )
    updateSelectInput(
      session, "a_exclusion_filter",
      choices = labelled_choices(exclusion_choices),
      selected = if (!is.null(current_exclusion) && current_exclusion %in% exclusion_choices) current_exclusion else default_exclusion
    )
  }, ignoreInit = FALSE)

  panel_a_data <- reactive({
    dat <- raw_data()
    req(
      !is.null(input$a_study_filter),
      !is.null(input$a_model_filter),
      !is.null(input$a_exclusion_filter)
    )

    dat %>%
      filter(
        study == input$a_study_filter,
        model == input$a_model_filter,
        exclusion == input$a_exclusion_filter
      )
  })

  global_labels <- reactive(make_global_outcome_labels(assoc_data()))
  predictor_labels <- reactive(merge_map(DEFAULT_PREDICTOR_LABELS, parse_named_mapping(input$predictor_label_map)))
  predictor_colors <- reactive(merge_map(DEFAULT_PREDICTOR_COLORS, parse_named_mapping(input$predictor_color_map)))

  b_outcomes <- reactive(parse_lines(input$b_outcomes))
  b_labels <- reactive(merge_map(global_labels(), parse_named_mapping(input$b_label_map)))
  b_label_offsets <- reactive(parse_xy_offset_mapping(input$b_label_offsets))
  c_outcomes <- reactive(parse_lines(input$c_outcomes))
  c_predictors <- reactive(resolve_predictors(parse_lines(input$c_predictors)))
  c_labels <- reactive(merge_map(merge_map(global_labels(), DEFAULT_C_LABELS), parse_named_mapping(input$c_label_map)))
  d_groups <- reactive(parse_named_mapping(input$d_groups))
  d_predictors <- reactive(resolve_predictors(parse_lines(input$d_predictors)))
  d_labels <- reactive(merge_map(merge_map(global_labels(), DEFAULT_D_LABELS), parse_named_mapping(input$d_label_map)))

  panel_a <- reactive({
    wrap_panel(
      plot_panel_a(
        panel_a_data(), input$a_significance_method, input$a_fdr_alpha,
        input$a_bonf_alpha, input$a_bonf_n_prs, input$a_bonf_n_biomarkers,
        input$a_block_spacing, input$a_total_gap, input$a_tile_width,
        input$a_tile_height, input$a_tile_text_size, input$a_white_threshold,
        input$a_axis_text_size, input$a_legend_title_size, input$a_legend_text_size,
        input$a_beta_limit, input$a_legend_position,
        parse_numeric_vector(input$a_colour_ticks), input$panel_border_colour, input$panel_border_width, input$margin_a_left, input$margin_right, input$margin_top, input$margin_bottom
      ), "a", input$panel_tag_size, input$panel_a_tag_x, input$panel_a_tag_y
    )
  })

  panel_b_raw <- reactive({
    plot_panel_b(
      assoc_data(), b_outcomes(), b_labels(), input$b_point_color,
      input$b_x_label, input$b_y_label,
      c(input$b_xmin, input$b_xmax), c(input$b_ymin, input$b_ymax),
      parse_numeric_vector(input$b_x_ticks), parse_numeric_vector(input$b_y_ticks),
      input$b_legend_x, input$b_legend_y, input$b_legend_ncol, input$b_legend_title_hjust, input$b_point_size, input$b_error_width, input$b_label_size,
      input$b_label_jitter_x, input$b_label_jitter_y, b_label_offsets(),
      input$b_show_error, input$b_show_labels, input$b_gridlines,
      axis_text_size = input$b_axis_text_size, axis_title_size = input$b_axis_title_size,
      legend_title_size = input$b_legend_title_size, legend_text_size = input$b_legend_text_size,
      panel_border_colour = input$panel_border_colour, panel_border_width = input$panel_border_width,
      margin_left = input$margin_b_left, margin_right = input$margin_right,
      margin_top = input$margin_top, margin_bottom = input$margin_bottom
    )
  })

  panel_b <- reactive({
    wrap_panel(
      panel_b_raw() + theme(legend.position = "none"),
      "b", input$panel_tag_size, input$panel_b_tag_x, input$panel_b_tag_y
    )
  })

  panel_b_legend <- reactive({
    wrap_elements(full = extract_legend_grob(panel_b_raw()))
  })

  panel_c_raw <- reactive({
    plot_panel_c(
      assoc_data(), c_outcomes(), c_predictors(), c_labels(),
      predictor_labels(), predictor_colors(), input$c_x_label,
      c(input$c_xmin, input$c_xmax), parse_numeric_vector(input$c_x_ticks),
      input$c_legend_x, input$c_legend_y, input$c_legend_ncol, input$c_legend_title_hjust, input$c_point_size,
      input$c_whisker_width, significance_threshold = input$c_p_threshold,
      dodge_width = input$c_dodge_width,
      outcome_text_size = input$c_outcome_text_size, outcome_wrap_width = input$c_label_wrap,
      axis_text_size = input$c_axis_text_size,
      axis_title_size = input$c_axis_title_size, legend_title_size = input$c_legend_title_size,
      legend_text_size = input$c_legend_text_size,
      row_background_height = input$c_row_background_height,
      y_bottom_padding = input$c_y_bottom_padding,
      y_top_padding = input$c_y_top_padding,
      panel_border_colour = input$panel_border_colour, panel_border_width = input$panel_border_width,
      margin_left = input$margin_d_left, margin_right = input$margin_d_right,
      margin_top = input$margin_top, margin_bottom = input$margin_bottom
    )
  })

  panel_c <- reactive({
    wrap_panel(
      panel_c_raw() + theme(legend.position = "none"),
      "d", input$panel_tag_size, input$panel_d_tag_x, input$panel_d_tag_y
    )
  })

  panel_c_legend <- reactive({
    wrap_elements(full = extract_legend_grob(panel_c_raw()))
  })

  panel_d_raw <- reactive({
    plot_panel_d(
      assoc_data(), d_groups(), d_predictors(), d_labels(),
      predictor_labels(), predictor_colors(),
      significance_method = input$d_significance_method,
      fdr_alpha = input$d_fdr_alpha, bonferroni_alpha = input$d_bonf_alpha,
      bonferroni_n_tests = input$d_bonf_n_tests,
      group_spacing = input$d_group_gap, tile_width = input$d_tile_width,
      tile_height = input$d_tile_height, tile_text_size = input$d_tile_text_size,
      white_label_threshold = input$d_white_threshold,
      beta_colour_limit = input$d_beta_limit,
      protein_text_size = input$d_protein_text_size, prs_text_size = input$d_prs_text_size,
      legend_title_size = input$d_legend_title_size, legend_text_size = input$d_legend_text_size,
      colourbar_width_cm = input$d_colourbar_width,
      label_wrap_width = input$d_label_wrap,
      legend_position = "bottom",
      colour_ticks = parse_numeric_vector(input$d_colour_ticks),
      legend_x_offset = input$d_legend_x, legend_y_offset = input$d_legend_y, legend_title_hjust = input$d_legend_title_hjust,
      panel_border_colour = input$panel_border_colour, panel_border_width = input$panel_border_width,
      margin_left = input$margin_c_left, margin_right = input$margin_right,
      margin_top = input$margin_top, margin_bottom = input$margin_bottom
    )
  })

  panel_d <- reactive({
    wrap_panel(
      panel_d_raw() + theme(legend.position = "none"),
      "c", input$panel_tag_size, input$panel_c_tag_x, input$panel_c_tag_y
    )
  })

  panel_d_legend <- reactive({
    wrap_elements(full = extract_legend_grob(panel_d_raw()))
  })

  combined_plot <- reactive({
    footer_height <- max(0.03, input$c_legend_height)
    main_height <- 1

    panel_b_stack <- free_left_side(panel_b()) / panel_b_legend() +
      plot_layout(heights = c(main_height, footer_height))

    panel_c_stack <- panel_c() / panel_c_legend() +
      plot_layout(heights = c(main_height, footer_height))

    panel_d_stack <- panel_d() / panel_d_legend() +
      plot_layout(heights = c(main_height, footer_height))

    bottom_row <- (panel_b_stack | panel_d_stack | panel_c_stack) +
      plot_layout(widths = c(input$col_b_width, input$col_c_width, input$col_d_width))

    panel_a() / bottom_row +
      plot_layout(heights = c(input$row_top, input$row_bottom))
  })

  output$combined_plot <- renderPlot(print(combined_plot()), bg = "white")

  output$data_table <- renderDT({
    datatable(head(assoc_data(), 500), options = list(scrollX = TRUE, pageLength = 12))
  })

  output$available_names <- renderPrint({
    dat <- assoc_data()
    cat("Predictors after alias resolution:\n")
    print(sort(unique(dat$predictor)))
    cat("\nOriginal predictor names:\n")
    print(sort(unique(dat$predictor_original)))
    cat("\nOutcomes:\n")
    print(sort(unique(dat$outcome)))
    cat("\nStudies:\n")
    print(sort(unique(dat$study)))
    cat("\nModels:\n")
    print(sort(unique(dat$model)))
    cat("\nExclusions:\n")
    print(sort(unique(dat$exclusion)))
  })

  observeEvent(input$load_settings, {
    if (is.null(input$settings_file)) {
      showNotification("Choose a JSON file first.", type = "warning")
      return()
    }
    tryCatch({
      object <- jsonlite::fromJSON(input$settings_file$datapath, simplifyVector = FALSE)
      apply_settings(object$settings %||% object, session)
      showNotification("Settings imported.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 8))
  })

  output$download_settings <- downloadHandler(
    filename = function() paste0("four_panel_prs_settings_", Sys.Date(), ".json"),
    content = function(file) jsonlite::write_json(collect_settings(input), file, pretty = TRUE, auto_unbox = TRUE, null = "null")
  )

  output$download_png <- downloadHandler(
    filename = function() paste0("four_panel_biomarker_prs_", Sys.Date(), ".png"),
    content = function(file) ggsave(file, combined_plot(), width = input$png_width,
                                    height = input$png_height, dpi = input$png_dpi, bg = "white")
  )

  output$download_pdf <- downloadHandler(
    filename = function() paste0("four_panel_biomarker_prs_", Sys.Date(), ".pdf"),
    content = function(file) ggsave(file, combined_plot(), width = input$png_width,
                                    height = input$png_height, device = cairo_pdf, bg = "white")
  )
}

shinyApp(ui, server)
