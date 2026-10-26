# ==============================================================================
# Shiny App: Three-panel pPRS Forest Plot Builder
# Panel B combined forest plot + Panel C selectable coefficient dot matrix
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
library(jsonlite)

# ==============================================================================
# Defaults
# ==============================================================================

KEEP_PREDICTORS <- c(
  "coagulation_prs_std",
  "mdwo_prs_std",
  "trl_prs_std",
  "endothelial_prs_std",
  "obesity_prs_std",
  "chol_prs_std",
  "lpa_prs_std",
  "cardiac_prs_std"
)

DEFAULT_PREDICTOR_ORDER <- c(
  "chol_prs_std",
  "trl_prs_std",
  "lpa_prs_std",
  "mdwo_prs_std",
  "obesity_prs_std",
  "endothelial_prs_std",
  "coagulation_prs_std",
  "cardiac_prs_std"
)

DEFAULT_PREDICTOR_LABELS <- c(
  chol_prs_std = "ChEL",
  trl_prs_std = "TRL",
  lpa_prs_std = "Lp(a)",
  mdwo_prs_std = "MUNW",
  obesity_prs_std = "Obesity",
  endothelial_prs_std = "Endothelial dysfunction",
  cardiac_prs_std = "Cardiac remodeling",
  coagulation_prs_std = "Coagulation"
)

DEFAULT_PREDICTOR_COLORS <- c(
  chol_prs_std = "#B22D4C",
  trl_prs_std = "#E07A8E",
  lpa_prs_std = "#7A1F3A",
  mdwo_prs_std = "#4A7C4E",
  obesity_prs_std = "#95C075",
  endothelial_prs_std = "#2D5F8B",
  coagulation_prs_std = "#7AAED1",
  cardiac_prs_std = "#5D4787"
)

DEFAULT_STUDY_ORDER <- c("UKBB", "SCAPIS", "MGBB", "Pooled")

DEFAULT_DOT_LINEAR_OUTCOMES <- c("ct_cac_tot_log", "ct_sis")
PANEL_C_NON_CAC_ADJUSTED_OUTCOME <- "ct_cac_tot_log"

DEFAULT_DOT_LOGISTIC_OUTCOMES <- c(
  "ct_50_stenosis_any",
  "ct_athero_noncalc",
  "ct_av_calcification"
)

DEFAULT_DOT_OUTCOME_LABELS <- c(
  ct_cac_tot = "CACS",
  ct_cac_tot_log = "CACS",
  ct_sis = "SIS",
  ct_50_stenosis_any = "Stenosis \u226550%",
  ct_athero_noncalc = "Non-calcified plaque",
  ct_av_calcification = "AV calcification"
)

shape_choices <- c(
  "circle filled" = 21,
  "square filled" = 22,
  "diamond filled" = 23,
  "triangle filled" = 24,
  "circle solid" = 16,
  "square solid" = 15,
  "diamond solid" = 18
)

linetype_choices <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

# ==============================================================================
# Example data
# ==============================================================================

generate_example_meta_data <- function() {
  set.seed(101)
  
  tidyr::expand_grid(
    meta_name = c("CAD broad", "Example secondary meta-analysis"),
    outcome_short = "CAD broad",
    predictor_short = KEEP_PREDICTORS,
    study = DEFAULT_STUDY_ORDER
  ) %>%
    mutate(
      row_type = if_else(study == "Pooled", "Pooled", "Study"),
      outcome_long = "Coronary artery disease",
      outcome_category = "Disease Prevalence",
      predictor_long = paste0(unname(DEFAULT_PREDICTOR_LABELS[predictor_short]), " PRS"),
      model = "Model 1: Age + Sex + PCs",
      exclusion = "Main analysis",
      meta_model = if_else(row_type == "Pooled", "Random effects", NA_character_),
      method = if_else(row_type == "Pooled", "REML", NA_character_),
      measure_scale = "or",
      scale_z = 3.509966639,
      test_tail = "two-sided",
      ci_tail = "two-sided",
      k = if_else(row_type == "Pooled", 3L, NA_integer_),
      n = if_else(row_type == "Study", sample(20000:460000, n(), replace = TRUE), NA_integer_),
      n_events = if_else(row_type == "Study", sample(500:30000, n(), replace = TRUE), NA_integer_),
      n_total = if_else(row_type == "Pooled", 486990L, NA_integer_),
      n_events_total = if_else(row_type == "Pooled", 76231L, NA_integer_),
      log_or = case_when(
        predictor_short == "chol_prs_std" ~ rnorm(n(), log(1.42), 0.07),
        predictor_short == "trl_prs_std" ~ rnorm(n(), log(1.20), 0.06),
        predictor_short == "lpa_prs_std" ~ rnorm(n(), log(1.25), 0.06),
        predictor_short == "mdwo_prs_std" ~ rnorm(n(), log(1.15), 0.06),
        predictor_short == "obesity_prs_std" ~ rnorm(n(), log(1.10), 0.05),
        predictor_short == "endothelial_prs_std" ~ rnorm(n(), log(1.12), 0.05),
        predictor_short == "cardiac_prs_std" ~ rnorm(n(), log(1.55), 0.08),
        predictor_short == "coagulation_prs_std" ~ rnorm(n(), log(1.08), 0.05),
        TRUE ~ rnorm(n(), log(1.10), 0.05)
      ),
      se = if_else(row_type == "Pooled", runif(n(), 0.025, 0.040), runif(n(), 0.035, 0.075)),
      or = exp(log_or),
      estimate = or,
      lb = exp(log_or - 1.96 * se),
      ub = exp(log_or + 1.96 * se),
      p_value = 2 * pnorm(-abs(log_or / se)),
      i2 = if_else(row_type == "Pooled", runif(n(), 0, 65), NA_real_),
      tau2 = if_else(row_type == "Pooled", runif(n(), 0, 0.01), NA_real_),
      qe = if_else(row_type == "Pooled", runif(n(), 0, 5), NA_real_),
      qe_p_value = if_else(row_type == "Pooled", runif(n(), 0.05, 1), NA_real_),
      studies = if_else(row_type == "Pooled", "UKBB; SCAPIS; MGBB", NA_character_),
      outcomes = outcome_short,
      outcome_longs = outcome_long
    ) %>%
    select(-log_or)
}

generate_example_imaging_data <- function() {
  set.seed(202)
  
  tidyr::expand_grid(
    study = "SCAPIS",
    outcome_short = c("ct_athero_any", "us_carotid_plaque_presence"),
    predictor_short = c(
      KEEP_PREDICTORS,
      "total_526_prs_std",
      "total_196_prs_std",
      "total_sum_cluster_prs_std"
    )
  ) %>%
    mutate(
      outcome_long = if_else(
        outcome_short == "ct_athero_any",
        "Coronary plaque presence",
        "Carotid plaque"
      ),
      outcome_category = "Imaging",
      predictor_long = paste0(predictor_short, " PRS"),
      model = "Model 1: Age + Sex + PCs",
      exclusion = "Main analysis (European ancestry)",
      n = 24708L,
      n_events = if_else(outcome_short == "ct_athero_any", 11000L, 13652L),
      log_or = case_when(
        predictor_short == "chol_prs_std" ~ log(1.12),
        predictor_short == "trl_prs_std" ~ log(1.08),
        predictor_short == "lpa_prs_std" ~ log(1.07),
        predictor_short == "mdwo_prs_std" ~ log(1.05),
        predictor_short == "obesity_prs_std" ~ log(1.02),
        predictor_short == "endothelial_prs_std" ~ log(1.04),
        predictor_short == "cardiac_prs_std" ~ log(1.06),
        predictor_short == "coagulation_prs_std" ~ log(1.02),
        TRUE ~ log(1.20)
      ) + rnorm(n(), 0, 0.015),
      se = runif(n(), 0.012, 0.018),
      or = exp(log_or),
      estimate = or,
      lb = exp(log_or - 1.96 * se),
      ub = exp(log_or + 1.96 * se),
      p_value = 2 * pnorm(-abs(log_or / se))
    ) %>%
    select(-log_or)
}

generate_example_dot_linear_data <- function() {
  set.seed(303)
  
  tidyr::expand_grid(
    study = "SCAPIS",
    outcome_short = c("ct_cac_tot", "ct_sis"),
    predictor_short = KEEP_PREDICTORS
  ) %>%
    mutate(
      beta = case_when(
        predictor_short == "chol_prs_std" ~ 0.08,
        predictor_short == "trl_prs_std" ~ 0.04,
        predictor_short == "lpa_prs_std" ~ 0.05,
        predictor_short == "mdwo_prs_std" ~ 0.03,
        predictor_short == "obesity_prs_std" ~ 0.02,
        predictor_short == "endothelial_prs_std" ~ 0.025,
        predictor_short == "cardiac_prs_std" ~ 0.06,
        predictor_short == "coagulation_prs_std" ~ 0.01,
        TRUE ~ 0.02
      ) + rnorm(n(), 0, 0.015),
      se = runif(n(), 0.012, 0.025),
      estimate = beta,
      lb = beta - 1.96 * se,
      ub = beta + 1.96 * se,
      p_value = 2 * pnorm(-abs(beta / se))
    )
}


generate_example_dot_linear_non_cac_adjusted_data <- function() {
  set.seed(304)
  
  tidyr::expand_grid(
    study = "SCAPIS",
    outcome_short = PANEL_C_NON_CAC_ADJUSTED_OUTCOME,
    predictor_short = KEEP_PREDICTORS
  ) %>%
    mutate(
      beta = case_when(
        predictor_short == "chol_prs_std" ~ 0.10,
        predictor_short == "trl_prs_std" ~ 0.05,
        predictor_short == "lpa_prs_std" ~ 0.06,
        predictor_short == "mdwo_prs_std" ~ 0.035,
        predictor_short == "obesity_prs_std" ~ 0.025,
        predictor_short == "endothelial_prs_std" ~ 0.03,
        predictor_short == "cardiac_prs_std" ~ 0.07,
        predictor_short == "coagulation_prs_std" ~ 0.015,
        TRUE ~ 0.02
      ) + rnorm(n(), 0, 0.015),
      se = runif(n(), 0.012, 0.025),
      estimate = beta,
      lb = beta - 1.96 * se,
      ub = beta + 1.96 * se,
      p_value = 2 * pnorm(-abs(beta / se))
    )
}

generate_example_dot_logistic_data <- function() {
  set.seed(404)
  
  tidyr::expand_grid(
    study = "SCAPIS",
    outcome_short = c("ct_50_stenosis_any", "ct_athero_noncalc", "ct_av_calcification"),
    predictor_short = KEEP_PREDICTORS
  ) %>%
    mutate(
      log_or = case_when(
        predictor_short == "chol_prs_std" ~ log(1.12),
        predictor_short == "trl_prs_std" ~ log(1.07),
        predictor_short == "lpa_prs_std" ~ log(1.10),
        predictor_short == "mdwo_prs_std" ~ log(1.04),
        predictor_short == "obesity_prs_std" ~ log(1.02),
        predictor_short == "endothelial_prs_std" ~ log(1.05),
        predictor_short == "cardiac_prs_std" ~ log(1.09),
        predictor_short == "coagulation_prs_std" ~ log(1.01),
        TRUE ~ log(1.04)
      ) + rnorm(n(), 0, 0.015),
      se = runif(n(), 0.012, 0.025),
      or = exp(log_or),
      estimate = or,
      lb = exp(log_or - 1.96 * se),
      ub = exp(log_or + 1.96 * se),
      p_value = 2 * pnorm(-abs(log_or / se))
    )
}

# ==============================================================================
# Helpers
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

safe_read_csv <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "gz") {
    read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

parse_lines <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character())
  out <- trimws(unlist(strsplit(x, "\n")))
  out[nzchar(out)]
}

parse_named_mapping <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(c())
  
  lines <- trimws(unlist(strsplit(text, "\n")))
  lines <- lines[nzchar(lines)]
  
  vals <- lapply(lines, function(ln) {
    parts <- unlist(strsplit(ln, "=", fixed = TRUE))
    if (length(parts) < 2) return(NULL)
    
    c(
      key = trimws(parts[1]),
      value = trimws(paste(parts[-1], collapse = "="))
    )
  })
  
  vals <- vals[!vapply(vals, is.null, logical(1))]
  if (!length(vals)) return(c())
  
  setNames(
    vapply(vals, function(x) x[["value"]], character(1)),
    vapply(vals, function(x) x[["key"]], character(1))
  )
}

parse_breaks <- function(x) {
  vals <- trimws(unlist(strsplit(x %||% "", ",")))
  vals <- vals[nzchar(vals)]
  suppressWarnings(as.numeric(vals))
}

clean_character_columns <- function(dat) {
  dat %>%
    mutate(across(where(is.character), ~ trimws(.x)))
}

wrap_if_needed <- function(x, width) {
  if (is.null(width) || is.na(width) || width <= 0) return(x)
  stringr::str_wrap(x, width = width)
}

first_existing_numeric <- function(dat, candidates) {
  existing <- candidates[candidates %in% names(dat)]
  
  if (!length(existing)) {
    return(rep(NA_real_, nrow(dat)))
  }
  
  suppressWarnings(as.numeric(dat[[existing[1]]]))
}

coerce_forest_columns <- function(dat) {
  dat <- clean_character_columns(dat)
  p_raw <- first_existing_numeric(dat, c("p_value", "p", "pval", "p.value"))
  i2_raw <- first_existing_numeric(dat, c("i2", "I2", "i_squared", "i.squared"))
  
  dat %>%
    mutate(
      or = if ("or" %in% names(.)) {
        suppressWarnings(as.numeric(.data$or))
      } else if ("estimate" %in% names(.)) {
        suppressWarnings(as.numeric(.data$estimate))
      } else {
        NA_real_
      },
      lb = if ("lb" %in% names(.)) suppressWarnings(as.numeric(.data$lb)) else NA_real_,
      ub = if ("ub" %in% names(.)) suppressWarnings(as.numeric(.data$ub)) else NA_real_,
      p_value = p_raw,
      i2 = i2_raw,
      or_plot = or,
      lb_plot = pmin(lb, ub, or, na.rm = TRUE),
      ub_plot = pmax(lb, ub, or, na.rm = TRUE),
      lb_plot = if_else(is.infinite(lb_plot), NA_real_, lb_plot),
      ub_plot = if_else(is.infinite(ub_plot), NA_real_, ub_plot)
    )
}

coerce_dot_columns <- function(dat, regression_type = c("linear", "logistic")) {
  regression_type <- match.arg(regression_type)
  
  dat <- clean_character_columns(dat)
  
  shiny::validate(
    shiny::need("outcome_short" %in% names(dat), "Dot-matrix data is missing outcome_short."),
    shiny::need("predictor_short" %in% names(dat), "Dot-matrix data is missing predictor_short.")
  )
  
  p_raw <- first_existing_numeric(dat, c("p_value", "p", "pval", "p.value"))
  
  if (regression_type == "linear") {
    effect <- first_existing_numeric(dat, c("beta", "estimate", "coef", "coefficient"))
    effect_lb <- first_existing_numeric(dat, c("lb", "ci_low", "lower", "lower_ci", "conf.low"))
    effect_ub <- first_existing_numeric(dat, c("ub", "ci_high", "upper", "upper_ci", "conf.high"))
  } else {
    if ("log_or" %in% names(dat)) {
      effect <- first_existing_numeric(dat, c("log_or"))
    } else {
      or_raw <- first_existing_numeric(dat, c("or", "estimate", "OR"))
      effect <- log(or_raw)
    }
    
    lb_raw <- first_existing_numeric(dat, c("lb", "ci_low", "lower", "lower_ci", "conf.low"))
    ub_raw <- first_existing_numeric(dat, c("ub", "ci_high", "upper", "upper_ci", "conf.high"))
    
    effect_lb <- ifelse(is.finite(lb_raw) & lb_raw > 0, log(lb_raw), NA_real_)
    effect_ub <- ifelse(is.finite(ub_raw) & ub_raw > 0, log(ub_raw), NA_real_)
  }
  
  dat %>%
    mutate(
      regression_type = regression_type,
      effect = as.numeric(effect),
      effect_lb = as.numeric(effect_lb),
      effect_ub = as.numeric(effect_ub),
      p_value = p_raw,
      sig = case_when(
        is.finite(effect_lb) & is.finite(effect_ub) ~ effect_lb > 0 | effect_ub < 0,
        is.finite(p_value) ~ p_value < 0.05,
        TRUE ~ FALSE
      ),
      abs_effect = abs(effect)
    )
}

validate_meta_data <- function(dat) {
  required <- c(
    "meta_name",
    "row_type",
    "study",
    "predictor_short",
    "meta_model",
    "or",
    "lb",
    "ub"
  )
  
  missing <- setdiff(required, names(dat))
  
  shiny::validate(
    shiny::need(
      length(missing) == 0,
      paste("Meta-analysis data is missing columns:", paste(missing, collapse = ", "))
    )
  )
  
  dat
}

validate_imaging_data <- function(dat) {
  required <- c(
    "study",
    "outcome_short",
    "predictor_short",
    "or",
    "lb",
    "ub"
  )
  
  missing <- setdiff(required, names(dat))
  
  shiny::validate(
    shiny::need(
      length(missing) == 0,
      paste("Imaging data is missing columns:", paste(missing, collapse = ", "))
    )
  )
  
  dat
}

make_x_mapper <- function(axis_min, axis_max, plot_x0, plot_x1) {
  force(axis_min)
  force(axis_max)
  force(plot_x0)
  force(plot_x1)
  
  function(x) {
    x <- as.numeric(x)
    x_clip <- pmin(pmax(x, axis_min), axis_max)
    
    scales::rescale(
      log10(x_clip),
      to = c(plot_x0, plot_x1),
      from = log10(c(axis_min, axis_max))
    )
  }
}

make_axis_df <- function(axis_breaks, axis_min, axis_max, plot_x0, plot_x1, y_axis) {
  mapper <- make_x_mapper(axis_min, axis_max, plot_x0, plot_x1)
  
  tibble(
    value = axis_breaks,
    x = mapper(axis_breaks),
    y = y_axis,
    label = as.character(axis_breaks)
  )
}

make_dot_outcome_labels <- function(outcomes, custom_label_map = c()) {
  default <- unname(DEFAULT_DOT_OUTCOME_LABELS[outcomes])
  default <- ifelse(is.na(default), outcomes, default)
  names(default) <- outcomes
  
  if (length(custom_label_map)) {
    matched <- intersect(names(custom_label_map), outcomes)
    default[matched] <- unname(custom_label_map[matched])
  }
  
  default
}

# ==============================================================================
# Data preparation
# ==============================================================================

prepare_meta_plot_data <- function(dat,
                                   meta_name_selected,
                                   meta_model_selected,
                                   predictor_order,
                                   predictor_labels,
                                   predictor_colors,
                                   study_order = DEFAULT_STUDY_ORDER,
                                   group_gap = 0.65,
                                   pathway_wrap_width = 18) {
  dat <- coerce_forest_columns(dat) %>%
    validate_meta_data() %>%
    mutate(
      meta_name = trimws(as.character(meta_name)),
      row_type = trimws(as.character(row_type)),
      study = trimws(as.character(study)),
      predictor_short = trimws(as.character(predictor_short)),
      meta_model = if_else(is.na(meta_model), "", trimws(as.character(meta_model)))
    )
  
  predictor_order <- intersect(predictor_order, KEEP_PREDICTORS)
  
  shiny::validate(
    shiny::need(
      length(predictor_order) > 0,
      "Panel a PRS order contains none of the eight allowed predictor_short values."
    )
  )
  
  filtered <- dat %>%
    filter(
      .data$meta_name == meta_name_selected,
      .data$predictor_short %in% predictor_order
    ) %>%
    filter(
      .data$row_type == "Study" |
        (.data$row_type == "Pooled" & .data$meta_model == meta_model_selected)
    ) %>%
    group_by(predictor_short, study, row_type) %>%
    slice(1) %>%
    ungroup()
  
  shiny::validate(
    shiny::need(
      nrow(filtered) > 0,
      "No rows remain after filtering meta-analysis data. Check meta_name, meta_model, study, row_type, and predictor_short."
    )
  )
  
  template <- tidyr::expand_grid(
    predictor_short = predictor_order,
    study = study_order
  ) %>%
    mutate(
      row_type = if_else(study == "Pooled", "Pooled", "Study"),
      predictor_id = match(predictor_short, predictor_order),
      study_id = match(study, study_order),
      y = rev(seq_len(n())) +
        rep(
          rev(seq(0, by = group_gap, length.out = length(predictor_order))),
          each = length(study_order)
        )
    )
  
  plot_df <- template %>%
    left_join(
      filtered,
      by = c("predictor_short", "study", "row_type"),
      suffix = c("", ".data")
    ) %>%
    mutate(
      predictor_label = unname(predictor_labels[predictor_short]),
      predictor_label = if_else(is.na(predictor_label), predictor_short, predictor_label),
      predictor_label = wrap_if_needed(predictor_label, pathway_wrap_width),
      predictor_color = unname(predictor_colors[predictor_short]),
      predictor_color = if_else(is.na(predictor_color), "grey40", predictor_color),
      study_label = study,
      pooled = study == "Pooled",
      sig = !is.na(lb_plot) & !is.na(ub_plot) & (lb_plot > 1 | ub_plot < 1),
      ci_color = if_else(sig, predictor_color, "grey55"),
      point_color = predictor_color
    )
  
  label_df <- plot_df %>%
    group_by(predictor_short, predictor_label) %>%
    summarise(y = mean(range(y, na.rm = TRUE)), .groups = "drop")
  
  list(rows = plot_df, labels = label_df)
}

prepare_combined_imaging_forest_data <- function(dat,
                                                 outcome_b,
                                                 outcome_c,
                                                 outcome_b_label,
                                                 outcome_c_label,
                                                 predictor_order,
                                                 predictor_labels,
                                                 predictor_colors,
                                                 group_gap = 0.26,
                                                 pathway_wrap_width = 18) {
  dat <- coerce_forest_columns(dat) %>%
    validate_imaging_data() %>%
    mutate(
      study = trimws(as.character(study)),
      outcome_short = trimws(as.character(outcome_short)),
      predictor_short = trimws(as.character(predictor_short))
    )
  
  predictor_order <- intersect(predictor_order, KEEP_PREDICTORS)
  
  shiny::validate(
    shiny::need(
      length(predictor_order) > 0,
      "Panel b PRS order contains none of the eight allowed predictor_short values."
    )
  )
  
  selected_outcomes <- c(outcome_b, outcome_c)
  
  filtered <- dat %>%
    filter(
      .data$outcome_short %in% selected_outcomes,
      .data$predictor_short %in% predictor_order
    ) %>%
    group_by(outcome_short, predictor_short) %>%
    slice(1) %>%
    ungroup()
  
  shiny::validate(
    shiny::need(
      nrow(filtered) > 0,
      "No rows remain for combined imaging forest plot. Check outcome_short and predictor_short."
    )
  )
  
  outcome_label_map <- c(
    setNames(outcome_b_label, outcome_b),
    setNames(outcome_c_label, outcome_c)
  )
  
  template <- tidyr::expand_grid(
    predictor_short = predictor_order,
    outcome_short = selected_outcomes
  ) %>%
    mutate(
      predictor_id = match(predictor_short, predictor_order),
      outcome_id = match(outcome_short, selected_outcomes),
      y_base = rev(seq_along(predictor_order))[predictor_id] +
        rev(seq(0, by = group_gap, length.out = length(predictor_order)))[predictor_id],
      y = y_base + case_when(
        outcome_id == 1L ~ 0.11,
        outcome_id == 2L ~ -0.11,
        TRUE ~ 0
      )
    )
  
  template %>%
    left_join(filtered, by = c("predictor_short", "outcome_short")) %>%
    mutate(
      outcome_label = unname(outcome_label_map[outcome_short]),
      outcome_label = if_else(is.na(outcome_label), outcome_short, outcome_label),
      predictor_label = unname(predictor_labels[predictor_short]),
      predictor_label = if_else(is.na(predictor_label), predictor_short, predictor_label),
      predictor_label = wrap_if_needed(predictor_label, pathway_wrap_width),
      predictor_color = unname(predictor_colors[predictor_short]),
      predictor_color = if_else(is.na(predictor_color), "grey40", predictor_color),
      sig = !is.na(lb_plot) & !is.na(ub_plot) & (lb_plot > 1 | ub_plot < 1),
      ci_color = if_else(sig, predictor_color, "grey55"),
      point_color = predictor_color
    )
}

prepare_dot_matrix_data <- function(linear_dat,
                                    logistic_dat,
                                    linear_outcomes,
                                    logistic_outcomes,
                                    predictor_order,
                                    predictor_labels,
                                    predictor_colors,
                                    dot_outcome_labels = c(),
                                    pathway_wrap_width = 18,
                                    outcome_wrap_width = 18) {
  predictor_order <- intersect(predictor_order, KEEP_PREDICTORS)
  
  shiny::validate(
    shiny::need(
      length(predictor_order) > 0,
      "Panel c PRS order contains none of the eight allowed predictor_short values."
    )
  )
  
  linear_outcomes <- linear_outcomes %||% character()
  logistic_outcomes <- logistic_outcomes %||% character()
  selected_outcomes <- c(linear_outcomes, logistic_outcomes)
  
  shiny::validate(
    shiny::need(
      length(selected_outcomes) > 0,
      "Select at least one linear or logistic outcome for panel c."
    )
  )
  
  linear_clean <- coerce_dot_columns(linear_dat, regression_type = "linear") %>%
    filter(outcome_short %in% linear_outcomes)
  
  logistic_clean <- coerce_dot_columns(logistic_dat, regression_type = "logistic") %>%
    filter(outcome_short %in% logistic_outcomes)
  
  dat <- bind_rows(linear_clean, logistic_clean) %>%
    filter(
      predictor_short %in% predictor_order,
      outcome_short %in% selected_outcomes
    ) %>%
    group_by(predictor_short, outcome_short) %>%
    slice(1) %>%
    ungroup()
  
  shiny::validate(
    shiny::need(
      nrow(dat) > 0,
      "No rows remain for panel c dot matrix. Check selected outcomes and predictor_short values."
    )
  )
  
  outcome_labels <- make_dot_outcome_labels(
    outcomes = selected_outcomes,
    custom_label_map = dot_outcome_labels
  )
  
  template <- tidyr::expand_grid(
    predictor_short = predictor_order,
    outcome_short = selected_outcomes
  ) %>%
    mutate(
      predictor_id = match(predictor_short, predictor_order),
      outcome_id = match(outcome_short, selected_outcomes),
      y = rev(seq_along(predictor_order))[predictor_id],
      x = outcome_id
    )
  
  template %>%
    left_join(dat, by = c("predictor_short", "outcome_short")) %>%
    mutate(
      predictor_label = unname(predictor_labels[predictor_short]),
      predictor_label = if_else(is.na(predictor_label), predictor_short, predictor_label),
      predictor_label = wrap_if_needed(predictor_label, pathway_wrap_width),
      predictor_color = unname(predictor_colors[predictor_short]),
      predictor_color = if_else(is.na(predictor_color), "grey40", predictor_color),
      outcome_label = unname(outcome_labels[outcome_short]),
      outcome_label = if_else(is.na(outcome_label), outcome_short, outcome_label),
      outcome_label = wrap_if_needed(outcome_label, outcome_wrap_width),
      abs_effect = if_else(is.finite(abs_effect), abs_effect, NA_real_),
      sig = if_else(is.na(sig), FALSE, sig),
      dot_color = if_else(sig, predictor_color, "grey40"),
      dot_size_value = if_else(sig & is.finite(abs_effect), abs_effect, NA_real_)
    )
}

# ==============================================================================
# Plot functions
# ==============================================================================

format_panel_a_ci <- function(lb, ub) {
  ifelse(
    is.finite(lb) & is.finite(ub),
    sprintf("%.2f–%.2f", lb, ub),
    ""
  )
}

plot_meta_forest <- function(plot_data,
                             title = "CAD meta-analysis",
                             panel_label = "a",
                             y_axis_title = "Process-specific PRS",
                             axis_min = 0.80,
                             axis_max = 1.80,
                             axis_breaks = c(0.80, 1.00, 1.25, 1.50, 1.80),
                             plot_x0 = 0,
                             plot_x1 = 1,
                             pathway_x = -0.68,
                             study_x = -0.28,
                             y_axis_title_x = -0.92,
                             y_axis_title_y = NA_real_,
                             xlim_left = -1.05,
                             xlim_right = 2.28,
                             study_point_size = 2.2,
                             pooled_point_size = 2.8,
                             study_shape = 21,
                             pooled_shape = 23,
                             ci_linewidth = 0.45,
                             ref_linewidth = 0.35,
                             pooled_bg = "grey92",
                             pathway_label_size = 3.2,
                             study_label_size = 2.9,
                             axis_text_size = 3,
                             axis_title_size = 3.2,
                             y_axis_title_size = 3.4,
                             title_size = 3.6,
                             panel_label_size = 5,
                             xlab = "Odds ratio",
                             pathway_label_col = "black",
                             study_label_col = "grey20",
                             title_col = "black",
                             axis_text_col = "black",
                             axis_title_col = "black",
                             y_axis_title_col = "black",
                             panel_label_col = "black",
                             ref_line_col = "grey40",
                             plot_box_col = "grey40",
                             plot_box_linewidth = 0.45,
                             tick_len = 0.16,
                             tick_text_pad = 0.35,
                             axis_title_pad = 0.85,
                             pooled_bar_height = 0.88,
                             plot_y_pad = 0.55,
                             title_y_pad = 1.35,
                             pathway_y_nudge = 0,
                             study_y_nudge = 0,
                             title_x = 0,
                             title_y_nudge = 0,
                             title_hjust = 0,
                             title_vjust = 0.5,
                             panel_label_x = -0.78,
                             panel_label_y = NA_real_,
                             panel_label_y_nudge = 0,
                             panel_label_hjust = 0,
                             panel_label_vjust = 0.5,
                             show_tick_guides = TRUE,
                             tick_guide_col = "grey88",
                             tick_guide_linewidth = 0.25,
                             tick_guide_linetype = "solid",
                             stats_or_x = 1.22,
                             stats_p_x = 1.62,
                             stats_i2_x = 2.02,
                             stats_header_size = 2.9,
                             stats_text_size = 2.75,
                             stats_text_col = "black") {
  rows <- plot_data$rows
  labels <- plot_data$labels
  
  map_x <- make_x_mapper(axis_min, axis_max, plot_x0, plot_x1)
  
  rows2 <- rows %>%
    mutate(
      x = map_x(or_plot),
      x_lb = map_x(lb_plot),
      x_ub = map_x(ub_plot),
      or_text = if_else(is.finite(or_plot), sprintf("%.2f", or_plot), ""),
      ci_text = format_panel_a_ci(lb_plot, ub_plot),
      i2_text = if_else(
        pooled & is.finite(i2),
        sprintf("%.0f%%", i2),
        ""
      )
    )
  
  labels2 <- labels %>%
    mutate(
      x = pathway_x,
      y = y + pathway_y_nudge
    )
  
  y_plot_min <- min(rows2$y, na.rm = TRUE) - plot_y_pad
  y_plot_max <- max(rows2$y, na.rm = TRUE) + plot_y_pad
  y_title <- y_plot_max + title_y_pad
  y_axis <- y_plot_min
  y_tick_end <- y_axis - tick_len
  y_tick_text <- y_tick_end - tick_text_pad
  y_x_axis_title_pos <- y_tick_text - axis_title_pad
  
  if (is.na(y_axis_title_y)) {
    y_axis_title_y <- mean(c(y_plot_min, y_plot_max))
  }
  
  panel_label_y_final <- if (is.na(panel_label_y)) {
    y_title + panel_label_y_nudge
  } else {
    panel_label_y
  }
  
  axis_df <- make_axis_df(axis_breaks, axis_min, axis_max, plot_x0, plot_x1, y_axis)
  
  guide_df <- axis_df %>%
    filter(
      abs(value - 1) > .Machine$double.eps^0.5,
      value < axis_max - .Machine$double.eps^0.5
    )
  
  pooled_bg_df <- rows2 %>%
    filter(pooled) %>%
    transmute(
      xmin = plot_x0,
      xmax = plot_x1,
      ymin = y - pooled_bar_height / 2,
      ymax = y + pooled_bar_height / 2
    )
  
  border_df <- tibble(
    x = c(plot_x0, plot_x1, plot_x1, plot_x0, plot_x0),
    y = c(y_plot_min, y_plot_min, y_plot_max, y_plot_max, y_plot_min)
  )
  
  p <- ggplot() +
    geom_rect(
      data = pooled_bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = pooled_bg,
      color = NA
    ) +
    geom_path(
      data = border_df,
      aes(x = x, y = y),
      color = plot_box_col,
      linewidth = plot_box_linewidth,
      lineend = "butt"
    )
  
  if (isTRUE(show_tick_guides)) {
    p <- p +
      geom_segment(
        data = guide_df,
        aes(x = x, xend = x, y = y_plot_min, yend = y_plot_max),
        color = tick_guide_col,
        linewidth = tick_guide_linewidth,
        linetype = tick_guide_linetype
      )
  }
  
  p +
    geom_segment(
      aes(x = map_x(1), xend = map_x(1), y = y_plot_min, yend = y_plot_max),
      linewidth = ref_linewidth,
      linetype = "dashed",
      color = ref_line_col
    ) +
    geom_segment(
      data = rows2,
      aes(x = x_lb, xend = x_ub, y = y, yend = y, color = ci_color),
      linewidth = ci_linewidth,
      lineend = "round",
      na.rm = TRUE
    ) +
    geom_point(
      data = rows2 %>% filter(!pooled),
      aes(x = x, y = y, fill = point_color, color = point_color),
      shape = study_shape,
      size = study_point_size,
      stroke = 0.25,
      na.rm = TRUE
    ) +
    geom_point(
      data = rows2 %>% filter(pooled),
      aes(x = x, y = y, fill = point_color, color = point_color),
      shape = pooled_shape,
      size = pooled_point_size,
      stroke = 0.30,
      na.rm = TRUE
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    geom_text(
      data = labels2,
      aes(x = x, y = y, label = predictor_label),
      hjust = 0,
      size = pathway_label_size,
      fontface = "bold",
      color = pathway_label_col,
      lineheight = 0.9
    ) +
    geom_text(
      data = rows2,
      aes(
        x = study_x,
        y = y + study_y_nudge,
        label = study_label,
        fontface = if_else(pooled, "bold", "plain")
      ),
      hjust = 0,
      size = study_label_size,
      color = study_label_col
    ) +
    annotate("text", x = stats_or_x, y = y_plot_max + 0.20, label = "OR",
             size = stats_header_size, fontface = "bold", color = stats_text_col, hjust = 0.5) +
    annotate("text", x = stats_p_x, y = y_plot_max + 0.20, label = "95%CI",
             size = stats_header_size, fontface = "bold", color = stats_text_col, hjust = 0.5) +
    annotate("text", x = stats_i2_x, y = y_plot_max + 0.20, label = expression(bold(I^2)),
             size = stats_header_size, fontface = "bold", color = stats_text_col, hjust = 0.5) +
    geom_text(data = rows2, aes(x = stats_or_x, y = y, label = or_text,
                                  fontface = if_else(pooled, "bold", "plain")),
              size = stats_text_size, color = stats_text_col, hjust = 0.5, na.rm = TRUE) +
    geom_text(data = rows2,
              aes(x = stats_p_x, y = y, label = ci_text,
                  fontface = if_else(pooled, "bold", "plain")),
              family = "Arial", size = stats_text_size, color = stats_text_col,
              hjust = 0.5, parse = FALSE, na.rm = TRUE) +
    geom_text(data = rows2, aes(x = stats_i2_x, y = y, label = i2_text,
                                  fontface = if_else(pooled, "bold", "plain")),
              size = stats_text_size, color = stats_text_col, hjust = 0.5, na.rm = TRUE) +
    geom_segment(
      data = axis_df,
      aes(x = x, xend = x, y = y_axis, yend = y_tick_end),
      color = axis_text_col,
      linewidth = 0.35
    ) +
    geom_text(
      data = axis_df,
      aes(x = x, y = y_tick_text, label = label),
      size = axis_text_size,
      color = axis_text_col,
      hjust = 0.5,
      vjust = 1
    ) +
    annotate(
      "text",
      x = mean(c(plot_x0, plot_x1)),
      y = y_x_axis_title_pos,
      label = xlab,
      size = axis_title_size,
      color = axis_title_col,
      hjust = 0.5,
      vjust = 1
    ) +
    annotate(
      "text",
      x = y_axis_title_x,
      y = y_axis_title_y,
      label = y_axis_title,
      size = y_axis_title_size,
      color = y_axis_title_col,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ) +
    annotate(
      "text",
      x = panel_label_x,
      y = panel_label_y_final,
      label = panel_label,
      size = panel_label_size,
      fontface = "bold",
      color = panel_label_col,
      hjust = panel_label_hjust,
      vjust = panel_label_vjust
    ) +
    annotate(
      "text",
      x = title_x,
      y = y_title + title_y_nudge,
      label = title,
      size = title_size,
      fontface = "bold",
      color = title_col,
      hjust = title_hjust,
      vjust = title_vjust
    ) +
    coord_cartesian(
      xlim = c(xlim_left, xlim_right),
      ylim = c(y_x_axis_title_pos - 0.25, y_title + 0.35),
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(8, 8, 8, 8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

plot_combined_imaging_forest <- function(rows,
                                         title = "Imaging outcomes",
                                         panel_label = "b",
                                         y_axis_title = "Process-specific PRS",
                                         axis_min = 0.90,
                                         axis_max = 1.30,
                                         axis_breaks = c(0.90, 1.00, 1.10, 1.20, 1.30),
                                         plot_x0 = 0,
                                         plot_x1 = 1,
                                         pathway_x = -0.52,
                                         y_axis_title_x = -0.76,
                                         y_axis_title_y = NA_real_,
                                         xlim_left = -0.85,
                                         xlim_right = 1.04,
                                         point_size = 2.4,
                                         outcome_shape_1 = 21,
                                         outcome_shape_2 = 24,
                                         ci_linewidth = 0.5,
                                         ref_linewidth = 0.35,
                                         pathway_label_size = 3.0,
                                         axis_text_size = 3,
                                         axis_title_size = 3.2,
                                         y_axis_title_size = 3.4,
                                         title_size = 3.4,
                                         panel_label_size = 5,
                                         xlab = "Odds ratio",
                                         pathway_label_col = "black",
                                         title_col = "black",
                                         axis_text_col = "black",
                                         axis_title_col = "black",
                                         y_axis_title_col = "black",
                                         panel_label_col = "black",
                                         ref_line_col = "grey40",
                                         plot_box_col = "grey40",
                                         plot_box_linewidth = 0.45,
                                         tick_len = 0.12,
                                         tick_text_pad = 0.25,
                                         axis_title_pad = 0.58,
                                         plot_y_pad = 0.45,
                                         title_y_pad = 1.35,
                                         pathway_y_nudge = 0,
                                         title_x = 0,
                                         title_y_nudge = 0,
                                         title_hjust = 0,
                                         title_vjust = 0.5,
                                         panel_label_x = -0.60,
                                         panel_label_y = NA_real_,
                                         panel_label_y_nudge = 0,
                                         panel_label_hjust = 0,
                                         panel_label_vjust = 0.5,
                                         legend_x = 0.02,
                                         legend_y = 1.02,
                                         legend_just_x = 0,
                                         legend_just_y = 0,
                                         legend_text_size = 8,
                                         legend_symbol_col = "grey60",
                                         show_tick_guides = TRUE,
                                         tick_guide_col = "grey88",
                                         tick_guide_linewidth = 0.25,
                                         tick_guide_linetype = "solid",
                                         alternating_fill = "grey95") {
  map_x <- make_x_mapper(axis_min, axis_max, plot_x0, plot_x1)
  
  rows2 <- rows %>%
    mutate(
      x = map_x(or_plot),
      x_lb = map_x(lb_plot),
      x_ub = map_x(ub_plot)
    )
  
  label_df <- rows2 %>%
    group_by(predictor_short, predictor_label) %>%
    summarise(y = mean(range(y, na.rm = TRUE)), .groups = "drop") %>%
    mutate(
      x = pathway_x,
      y = y + pathway_y_nudge
    )
  
  y_plot_min <- min(rows2$y, na.rm = TRUE) - plot_y_pad
  y_plot_max <- max(rows2$y, na.rm = TRUE) + plot_y_pad
  y_title <- y_plot_max + title_y_pad
  y_axis <- y_plot_min
  y_tick_end <- y_axis - tick_len
  y_tick_text <- y_tick_end - tick_text_pad
  y_x_axis_title_pos <- y_tick_text - axis_title_pad
  
  if (is.na(y_axis_title_y)) {
    y_axis_title_y <- mean(c(y_plot_min, y_plot_max))
  }
  
  panel_label_y_final <- if (is.na(panel_label_y)) {
    y_title + panel_label_y_nudge
  } else {
    panel_label_y
  }
  
  axis_df <- make_axis_df(axis_breaks, axis_min, axis_max, plot_x0, plot_x1, y_axis)
  
  guide_df <- axis_df %>%
    filter(
      abs(value - 1) > .Machine$double.eps^0.5,
      value < axis_max - .Machine$double.eps^0.5
    )
  
  border_df <- tibble(
    x = c(plot_x0, plot_x1, plot_x1, plot_x0, plot_x0),
    y = c(y_plot_min, y_plot_min, y_plot_max, y_plot_max, y_plot_min)
  )

  stripe_df <- rows2 %>%
    group_by(predictor_short) %>%
    summarise(ymin = min(y) - 0.48, ymax = max(y) + 0.48, .groups = "drop") %>%
    arrange(desc(ymax)) %>%
    mutate(show = row_number() %% 2 == 0) %>%
    filter(show) %>%
    transmute(xmin = plot_x0, xmax = plot_x1, ymin, ymax)
  
  shape_values <- setNames(
    c(outcome_shape_1, outcome_shape_2)[seq_along(unique(rows2$outcome_label))],
    unique(rows2$outcome_label)
  )
  
  p <- ggplot() +
    geom_rect(
      data = stripe_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = alternating_fill, color = NA
    ) +
    geom_path(
      data = border_df,
      aes(x = x, y = y),
      color = plot_box_col,
      linewidth = plot_box_linewidth,
      lineend = "butt"
    )
  
  if (isTRUE(show_tick_guides)) {
    p <- p +
      geom_segment(
        data = guide_df,
        aes(x = x, xend = x, y = y_plot_min, yend = y_plot_max),
        color = tick_guide_col,
        linewidth = tick_guide_linewidth,
        linetype = tick_guide_linetype
      )
  }
  
  p +
    geom_segment(
      aes(x = map_x(1), xend = map_x(1), y = y_plot_min, yend = y_plot_max),
      linewidth = ref_linewidth,
      linetype = "dashed",
      color = ref_line_col
    ) +
    geom_segment(
      data = rows2,
      aes(x = x_lb, xend = x_ub, y = y, yend = y, color = ci_color),
      linewidth = ci_linewidth,
      lineend = "round",
      na.rm = TRUE
    ) +
    geom_point(
      data = rows2,
      aes(x = x, y = y, fill = point_color, color = point_color, shape = outcome_label),
      size = point_size,
      stroke = 0.25,
      na.rm = TRUE
    ) +
    scale_shape_manual(values = shape_values) +
    scale_color_identity() +
    scale_fill_identity() +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = predictor_label),
      hjust = 0,
      size = pathway_label_size,
      fontface = "bold",
      color = pathway_label_col,
      lineheight = 0.9
    ) +
    geom_segment(
      data = axis_df,
      aes(x = x, xend = x, y = y_axis, yend = y_tick_end),
      color = axis_text_col,
      linewidth = 0.35
    ) +
    geom_text(
      data = axis_df,
      aes(x = x, y = y_tick_text, label = label),
      size = axis_text_size,
      color = axis_text_col,
      hjust = 0.5,
      vjust = 1
    ) +
    annotate(
      "text",
      x = mean(c(plot_x0, plot_x1)),
      y = y_x_axis_title_pos,
      label = xlab,
      size = axis_title_size,
      color = axis_title_col,
      hjust = 0.5,
      vjust = 1
    ) +
    annotate(
      "text",
      x = y_axis_title_x,
      y = y_axis_title_y,
      label = y_axis_title,
      size = y_axis_title_size,
      color = y_axis_title_col,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ) +
    annotate(
      "text",
      x = panel_label_x,
      y = panel_label_y_final,
      label = panel_label,
      size = panel_label_size,
      fontface = "bold",
      color = panel_label_col,
      hjust = panel_label_hjust,
      vjust = panel_label_vjust
    ) +
    annotate(
      "text",
      x = title_x,
      y = y_title + title_y_nudge,
      label = title,
      size = title_size,
      fontface = "bold",
      color = title_col,
      hjust = title_hjust,
      vjust = title_vjust
    ) +
    coord_cartesian(
      xlim = c(xlim_left, xlim_right),
      ylim = c(y_x_axis_title_pos - 0.18, y_title + 0.35),
      clip = "off"
    ) +
    guides(
      shape = guide_legend(
        title = "Plaque",
        title.position = "left",
        title.hjust = 0,
        override.aes = list(fill = legend_symbol_col, color = legend_symbol_col)
      )
    ) +
    theme_void() +
    theme(
      legend.position = c(legend_x, legend_y),
      legend.justification = c(legend_just_x, legend_just_y),
      legend.direction = "horizontal",
      legend.title = element_text(size = legend_text_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(8, 8, 8, 8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

plot_dot_matrix <- function(rows,
                            title = "CAC-adjusted imaging associations",
                            panel_label = "c",
                            y_axis_title = "Process-specific PRS",
                            pathway_x = 0.28,
                            y_axis_title_x = -0.10,
                            xlim_left = -0.20,
                            xlim_right = NULL,
                            plot_y_pad = 0.50,
                            title_y_pad = 1.15,
                            x_label_y_pad = -0.55,
                            column_label_y = NA_real_,
                            pathway_label_size = 3.0,
                            column_label_size = 2.8,
                            column_label_angle = 0,
                            column_label_hjust = 0.5,
                            column_label_vjust = 1,
                            y_axis_title_size = 3.4,
                            title_size = 3.4,
                            panel_label_size = 5,
                            panel_label_x = -0.05,
                            panel_label_y = NA_real_,
                            title_x = 1,
                            title_hjust = 0,
                            title_vjust = 0.5,
                            panel_label_hjust = 0,
                            panel_label_vjust = 0.5,
                            dot_min_size = 1.2,
                            dot_max_size = 6.0,
                            dot_scale_multiplier = 1.0,
                            dot_legend_breaks_n = 3,
                            dot_legend_title = "Regression coefficient",
                            nonsig_dot_size = 1.4,
                            nonsig_dot_col = "grey40",
                            show_grid = TRUE,
                            grid_col = "grey88",
                            grid_linewidth = 0.35,
                            border_col = "grey80",
                            border_linewidth = 0.45,
                            legend_x = 0.50,
                            legend_y = 1.04,
                            legend_just_x = 0.5,
                            legend_just_y = 0,
                            legend_title_size = 8,
                            legend_text_size = 8,
                            legend_symbol_col = "grey60",
                            pathway_label_col = "black",
                            y_axis_title_col = "black",
                            title_col = "black",
                            panel_label_col = "black",
                            column_label_col = "black") {
  rows2 <- rows %>%
    mutate(
      x = as.numeric(x),
      y = as.numeric(y),
      dot_size_plot = if_else(
        sig & is.finite(dot_size_value),
        dot_size_value * dot_scale_multiplier,
        NA_real_
      )
    )
  
  n_cols <- length(unique(rows2$outcome_short))
  
  plot_x_min <- 0.5
  plot_x_max <- n_cols + 0.5
  
  if (is.null(xlim_right) || !is.finite(xlim_right)) {
    xlim_right <- plot_x_max + 0.10
  }
  
  plot_y_min <- min(rows2$y, na.rm = TRUE) - plot_y_pad
  plot_y_max <- max(rows2$y, na.rm = TRUE) + plot_y_pad
  y_title <- plot_y_max + title_y_pad
  
  panel_label_y_final <- if (is.na(panel_label_y)) {
    y_title
  } else {
    panel_label_y
  }
  
  y_column_label <- if (is.na(column_label_y)) {
    plot_y_min + x_label_y_pad
  } else {
    column_label_y
  }
  
  y_bottom <- min(y_column_label - 0.35, plot_y_min - 0.10)
  
  label_df <- rows2 %>%
    distinct(predictor_short, predictor_label, y) %>%
    mutate(x = pathway_x)
  
  column_df <- rows2 %>%
    distinct(outcome_short, outcome_label, x) %>%
    arrange(x)
  
  grid_v <- tibble(
    x = seq_len(n_cols),
    xend = seq_len(n_cols),
    y = plot_y_min,
    yend = plot_y_max
  )
  
  grid_h <- tibble(
    y = sort(unique(rows2$y)),
    yend = sort(unique(rows2$y)),
    x = plot_x_min,
    xend = plot_x_max
  )
  
  border_df <- tibble(
    x = c(plot_x_min, plot_x_max, plot_x_max, plot_x_min, plot_x_min),
    y = c(plot_y_min, plot_y_min, plot_y_max, plot_y_max, plot_y_min)
  )
  
  max_abs <- suppressWarnings(max(rows2$dot_size_plot, na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 0.6
  legend_max <- max(max_abs, 0.6)
  
  p <- ggplot()
  
  if (isTRUE(show_grid)) {
    p <- p +
      geom_segment(
        data = grid_v,
        aes(x = x, xend = xend, y = y, yend = yend),
        color = grid_col,
        linewidth = grid_linewidth
      ) +
      geom_segment(
        data = grid_h,
        aes(x = x, xend = xend, y = y, yend = yend),
        color = grid_col,
        linewidth = grid_linewidth
      )
  }
  
  p +
    geom_path(
      data = border_df,
      aes(x = x, y = y),
      color = border_col,
      linewidth = border_linewidth,
      lineend = "butt"
    ) +
    geom_point(
      data = rows2 %>% filter(!sig),
      aes(x = x, y = y),
      shape = 16,
      size = nonsig_dot_size,
      color = nonsig_dot_col,
      na.rm = TRUE
    ) +
    geom_point(
      data = rows2 %>% filter(sig),
      aes(x = x, y = y, size = dot_size_plot, fill = dot_color, color = dot_color),
      shape = 21,
      stroke = 0.25,
      na.rm = TRUE
    ) +
    scale_size_continuous(
      name = dot_legend_title,
      range = c(dot_min_size, dot_max_size),
      limits = c(0, legend_max),
      breaks = c(0.2, 0.4, 0.6)
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = predictor_label),
      hjust = 0,
      size = pathway_label_size,
      fontface = "bold",
      color = pathway_label_col,
      lineheight = 0.9
    ) +
    geom_text(
      data = column_df,
      aes(x = x, y = y_column_label, label = outcome_label),
      hjust = column_label_hjust,
      vjust = column_label_vjust,
      angle = column_label_angle,
      size = column_label_size,
      color = column_label_col,
      lineheight = 0.9
    ) +
    annotate(
      "text",
      x = y_axis_title_x,
      y = mean(c(plot_y_min, plot_y_max)),
      label = y_axis_title,
      size = y_axis_title_size,
      color = y_axis_title_col,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ) +
    annotate(
      "text",
      x = panel_label_x,
      y = panel_label_y_final,
      label = panel_label,
      size = panel_label_size,
      fontface = "bold",
      color = panel_label_col,
      hjust = panel_label_hjust,
      vjust = panel_label_vjust
    ) +
    annotate(
      "text",
      x = title_x,
      y = y_title,
      label = title,
      size = title_size,
      fontface = "bold",
      color = title_col,
      hjust = title_hjust,
      vjust = title_vjust
    ) +
    coord_cartesian(
      xlim = c(xlim_left, xlim_right),
      ylim = c(y_bottom, y_title + 0.35),
      clip = "off"
    ) +
    guides(
      size = guide_legend(
        title.position = "left",
        title.hjust = 0,
        override.aes = list(fill = legend_symbol_col, color = legend_symbol_col)
      )
    ) +
    theme_void() +
    theme(
      legend.position = c(legend_x, legend_y),
      legend.justification = c(legend_just_x, legend_just_y),
      legend.direction = "horizontal",
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(8, 8, 8, 8),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

combine_three_panels <- function(plot_a,
                                 plot_b,
                                 plot_c,
                                 column_width_a = 1.25,
                                 column_width_bc = 1.00,
                                 row_height_b = 1,
                                 row_height_c = 1,
                                 gap_width_ab = 0.05,
                                 gap_height_bc = 0.05) {
  right_col <-
    plot_b /
    patchwork::plot_spacer() /
    plot_c +
    plot_layout(heights = c(row_height_b, gap_height_bc, row_height_c))
  
  plot_a |
    patchwork::plot_spacer() |
    right_col +
    plot_layout(widths = c(column_width_a, gap_width_ab, column_width_bc))
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  titlePanel("Three-panel pPRS Forest Plot Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      tabsetPanel(
        id = "control_tabs",
        
        tabPanel(
          "Data",
          checkboxInput("use_example", "Use example data", TRUE),
          fluidRow(
            column(6, downloadButton("download_settings", "Export settings (JSON)", width = "100%")),
            column(6, fileInput("settings_file", "Import settings (JSON)", accept = ".json", buttonLabel = "Import", placeholder = "No JSON selected"))
          ),
          fileInput("meta_file", "Upload meta-analysis CSV", accept = c(".csv", ".gz")),
          fileInput("imaging_file", "Upload SCAPIS imaging logistic CSV", accept = c(".csv", ".gz")),
          fileInput("dot_linear_file", "Upload CAC-adjusted linear CSV", accept = c(".csv", ".gz")),
          fileInput(
            "dot_linear_non_cac_file",
            "Upload non-CAC-adjusted linear CSV (source for ct_cac_tot_log)",
            accept = c(".csv", ".gz")
          ),
          fileInput("dot_logistic_file", "Upload CAC-adjusted logistic CSV", accept = c(".csv", ".gz")),
          tags$hr(),
          uiOutput("meta_name_ui"),
          selectInput(
            "meta_model",
            "Meta-analysis model for pooled row",
            choices = c("Fixed effect", "Random effects"),
            selected = "Fixed effect"
          ),
          tags$hr(),
          uiOutput("coronary_outcome_ui"),
          uiOutput("carotid_outcome_ui"),
          textInput("coronary_outcome_label", "Panel b symbol 1 label", "Coronary plaque"),
          textInput("carotid_outcome_label", "Panel b symbol 2 label", "Carotid plaque"),
          tags$hr(),
          uiOutput("dot_linear_outcomes_ui"),
          uiOutput("dot_logistic_outcomes_ui"),
          textAreaInput(
            "dot_outcome_label_map",
            "Panel c column labels: outcome_short = label",
            value = paste0(
              names(DEFAULT_DOT_OUTCOME_LABELS),
              " = ",
              unname(DEFAULT_DOT_OUTCOME_LABELS),
              collapse = "\n"
            ),
            rows = 5
          ),
          helpText(
            "Panel b combines two selected imaging outcomes. ",
            "Panel c uses selected linear outcomes as coefficients and selected logistic outcomes as log(OR)."
          )
        ),
        
        tabPanel(
          "PRS order/labels/colors",
          textAreaInput(
            "order_a",
            "Panel a PRS order",
            value = paste(DEFAULT_PREDICTOR_ORDER, collapse = "\n"),
            rows = 8
          ),
          textAreaInput(
            "order_b",
            "Panel b PRS order",
            value = paste(DEFAULT_PREDICTOR_ORDER, collapse = "\n"),
            rows = 8
          ),
          textAreaInput(
            "order_c",
            "Panel c PRS order",
            value = paste(DEFAULT_PREDICTOR_ORDER, collapse = "\n"),
            rows = 8
          ),
          tags$hr(),
          textAreaInput(
            "label_map",
            "PRS labels: predictor_short = label",
            value = paste0(
              names(DEFAULT_PREDICTOR_LABELS),
              " = ",
              unname(DEFAULT_PREDICTOR_LABELS),
              collapse = "\n"
            ),
            rows = 8
          ),
          textAreaInput(
            "color_map",
            "PRS colors: predictor_short = hex/name",
            value = paste0(
              names(DEFAULT_PREDICTOR_COLORS),
              " = ",
              unname(DEFAULT_PREDICTOR_COLORS),
              collapse = "\n"
            ),
            rows = 8
          ),
          numericInput("pathway_wrap_width", "PRS label wrap after characters", 18, min = 1, step = 1),
          numericInput("outcome_wrap_width", "Panel c outcome label wrap after characters", 18, min = 1, step = 1)
        ),
        
        tabPanel(
          "Axes",
          h4("Panel a"),
          fluidRow(
            column(6, numericInput("axis_min_a", "x min", 0.80, step = 0.01)),
            column(6, numericInput("axis_max_a", "x max", 1.80, step = 0.01))
          ),
          textInput("axis_breaks_a", "x breaks", "0.80,1.00,1.25,1.50,1.80"),
          textInput("xlab_a", "x-axis label", "Odds ratio"),
          textInput("y_axis_title_a", "vertical y-axis title", "Process-specific PRS"),
          
          h4("Panel b"),
          fluidRow(
            column(6, numericInput("axis_min_b", "x min", 0.90, step = 0.01)),
            column(6, numericInput("axis_max_b", "x max", 1.30, step = 0.01))
          ),
          textInput("axis_breaks_b", "x breaks", "0.90,1.00,1.10,1.20,1.30"),
          textInput("xlab_b", "x-axis label", "Odds ratio"),
          textInput("y_axis_title_b", "vertical y-axis title", "Process-specific PRS"),
          
          h4("Panel c"),
          textInput("y_axis_title_c", "vertical y-axis title", "Process-specific PRS")
        ),
        
        tabPanel(
          "Manual positions",
          h4("Panel a horizontal positions"),
          fluidRow(
            column(4, numericInput("pathway_x_a", "Pathway x", -0.68, step = 0.02)),
            column(4, numericInput("study_x_a", "Study x", -0.28, step = 0.02)),
            column(4, numericInput("panel_label_x_a", "Panel letter x", -0.78, step = 0.02))
          ),
          fluidRow(
            column(6, numericInput("panel_label_y_a", "Panel letter y; NA = auto", NA, step = 0.05)),
            column(6, numericInput("panel_label_y_nudge_a", "Panel letter y nudge", 0, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("xlim_left_a", "xlim left", -1.05, step = 0.02)),
            column(4, numericInput("xlim_right_a", "xlim right", 2.28, step = 0.02)),
            column(4, numericInput("title_x_a", "Title x", 0, step = 0.02))
          ),
          h4("Panel a statistics table positions"),
          fluidRow(
            column(4, numericInput("stats_or_x_a", "OR x", 1.22, step = 0.02)),
            column(4, numericInput("stats_p_x_a", "95% CI x", 1.62, step = 0.02)),
            column(4, numericInput("stats_i2_x_a", "I2 x", 2.02, step = 0.02))
          ),
          fluidRow(
            column(6, numericInput("stats_header_size_a", "Statistics header size", 2.9, step = 0.1)),
            column(6, numericInput("stats_text_size_a", "Statistics text size", 2.75, step = 0.1))
          ),
          fluidRow(
            column(4, numericInput("y_axis_title_x_a", "Y-title x", -0.92, step = 0.02)),
            column(4, numericInput("y_axis_title_y_a", "Y-title y; NA = center", NA, step = 0.1)),
            column(4, numericInput("title_y_nudge_a", "Title y nudge", 0, step = 0.05))
          ),
          
          h4("Panel b horizontal positions"),
          fluidRow(
            column(4, numericInput("pathway_x_b", "Pathway x", -0.52, step = 0.02)),
            column(4, numericInput("panel_label_x_b", "Panel letter x", -0.60, step = 0.02)),
            column(4, numericInput("title_x_b", "Title x", 0, step = 0.02))
          ),
          fluidRow(
            column(6, numericInput("panel_label_y_b", "Panel letter y; NA = auto", NA, step = 0.05)),
            column(6, numericInput("panel_label_y_nudge_b", "Panel letter y nudge", 0, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("xlim_left_b", "xlim left", -0.85, step = 0.02)),
            column(4, numericInput("xlim_right_b", "xlim right", 1.04, step = 0.02)),
            column(4, numericInput("title_y_nudge_b", "Title y nudge", 0, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("y_axis_title_x_b", "Y-title x", -0.76, step = 0.02)),
            column(6, numericInput("y_axis_title_y_b", "Y-title y; NA = center", NA, step = 0.1))
          ),
          
          h4("Panel c positions"),
          fluidRow(
            column(4, numericInput("pathway_x_c", "Pathway x", 0.28, step = 0.02)),
            column(4, numericInput("panel_label_x_c", "Panel letter x", -0.05, step = 0.02)),
            column(4, numericInput("title_x_c", "Title x", 1.00, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("panel_label_y_c", "Panel letter y; NA = auto", NA, step = 0.05)),
            column(6, helpText("Panel C letter y is absolute. Use NA for automatic placement at title height."))
          ),
          fluidRow(
            column(4, numericInput("xlim_left_c", "xlim left", -0.20, step = 0.05)),
            column(4, numericInput("xlim_right_c", "xlim right; NA = auto", NA, step = 0.05)),
            column(4, numericInput("y_axis_title_x_c", "Y-title x", -0.10, step = 0.02))
          ),
          
          h4("Vertical positions"),
          fluidRow(
            column(4, numericInput("meta_group_gap", "Panel a group gap", 0.65, step = 0.05)),
            column(4, numericInput("imaging_group_gap", "Panel b row gap", 0.26, step = 0.05)),
            column(4, numericInput("pooled_bar_height", "Pooled bar height", 0.88, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("plot_y_pad_a", "Panel a plot y pad", 0.55, step = 0.05)),
            column(4, numericInput("title_y_pad_a", "Panel a title y pad", 1.35, step = 0.05)),
            column(4, numericInput("pathway_y_nudge_a", "Panel a pathway y nudge", 0, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("plot_y_pad_b", "Panel b plot y pad", 0.45, step = 0.05)),
            column(4, numericInput("title_y_pad_b", "Panel b title y pad", 1.35, step = 0.05)),
            column(4, numericInput("pathway_y_nudge_b", "Panel b pathway y nudge", 0, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("plot_y_pad_c", "Panel c plot y pad", 0.50, step = 0.05)),
            column(4, numericInput("title_y_pad_c", "Panel c title y pad", 1.15, step = 0.05)),
            column(4, numericInput("column_label_y_pad_c", "Panel c column label y pad; used if label y is NA", -0.55, step = 0.05))
          ),
          numericInput("column_label_y_c", "Panel c column label y; NA = use y pad", NA, step = 0.05)
        ),
        
        tabPanel(
          "Style",
          h4("Points and lines"),
          fluidRow(
            column(4, numericInput("study_point_size", "Study point size", 2.2, step = 0.1)),
            column(4, numericInput("pooled_point_size", "Pooled point size", 2.8, step = 0.1)),
            column(4, numericInput("imaging_point_size", "Panel b point size", 2.4, step = 0.1))
          ),
          fluidRow(
            column(4, selectInput("study_shape", "Study shape", choices = shape_choices, selected = 21)),
            column(4, selectInput("pooled_shape", "Pooled shape", choices = shape_choices, selected = 23)),
            column(4, selectInput("panel_b_shape_1", "Panel b symbol 1", choices = shape_choices, selected = 22))
          ),
          fluidRow(
            column(4, selectInput("panel_b_shape_2", "Panel b symbol 2", choices = shape_choices, selected = 24)),
            column(4, numericInput("ci_linewidth", "CI linewidth", 0.5, step = 0.05)),
            column(4, numericInput("ref_linewidth", "Reference linewidth", 0.35, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("plot_box_linewidth", "Box linewidth", 0.45, step = 0.05)),
            column(4, numericInput("dot_min_size", "Panel c significant min dot size", 1.2, step = 0.1)),
            column(4, numericInput("dot_max_size", "Panel c significant max dot size", 6.0, step = 0.1))
          ),
          fluidRow(
            column(4, numericInput("nonsig_dot_size", "Panel c non-significant dot size", 1.4, step = 0.1)),
            column(4, textInput("nonsig_dot_col", "Panel c non-significant dot color", "grey40")),
            column(4, numericInput("dot_scale_multiplier", "Panel c significant dot scale multiplier", 1.0, step = 0.1))
          ),
          fluidRow(
            column(6, numericInput("dot_legend_breaks_n", "Panel c dot legend breaks", 3, min = 2, step = 1)),
            column(6, checkboxInput("show_dot_grid", "Show Panel c gridlines", TRUE))
          ),
          
          h4("Panel b legend"),
          fluidRow(
            column(3, numericInput("legend_b_x", "x", 0.02, step = 0.02)),
            column(3, numericInput("legend_b_y", "y", 1.02, step = 0.02)),
            column(3, numericInput("legend_b_just_x", "just x", 0, step = 0.1)),
            column(3, numericInput("legend_b_just_y", "just y", 0, step = 0.1))
          ),
          numericInput("legend_b_text_size", "Panel b legend text size", 8, step = 0.5),
          textInput("legend_b_symbol_col", "Panel b legend symbol color", "grey60"),
          textInput("panel_b_alternating_fill", "Panel b alternating row color", "grey95"),
          
          h4("Panel c legend"),
          textInput("dot_legend_title", "Panel c legend title", "Regression coefficient"),
          fluidRow(
            column(3, numericInput("legend_c_x", "x", 0.50, step = 0.02)),
            column(3, numericInput("legend_c_y", "y", 1.04, step = 0.02)),
            column(3, numericInput("legend_c_just_x", "just x", 0.5, step = 0.1)),
            column(3, numericInput("legend_c_just_y", "just y", 0, step = 0.1))
          ),
          fluidRow(
            column(6, numericInput("legend_c_title_size", "Panel c legend title size", 8, step = 0.5)),
            column(6, numericInput("legend_c_text_size", "Panel c legend text size", 8, step = 0.5))
          ),
          textInput("legend_c_symbol_col", "Panel c legend symbol color", "grey60"),
          
          h4("Vertical tick guide lines"),
          checkboxInput("show_tick_guides_a", "Show vertical tick guide lines in panel a", TRUE),
          checkboxInput("show_tick_guides_b", "Show vertical tick guide lines in panel b", TRUE),
          fluidRow(
            column(4, textInput("tick_guide_col_a", "Panel a tick guide color", "grey88")),
            column(4, numericInput("tick_guide_linewidth_a", "Panel a tick guide linewidth", 0.25, step = 0.05)),
            column(4, selectInput("tick_guide_linetype_a", "Panel a tick guide linetype", choices = linetype_choices, selected = "solid"))
          ),
          fluidRow(
            column(4, textInput("tick_guide_col_b", "Panel b tick guide color", "grey88")),
            column(4, numericInput("tick_guide_linewidth_b", "Panel b tick guide linewidth", 0.25, step = 0.05)),
            column(4, selectInput("tick_guide_linetype_b", "Panel b tick guide linetype", choices = linetype_choices, selected = "solid"))
          ),
          
          h4("Colors"),
          fluidRow(
            column(4, textInput("pooled_bg", "Pooled row background", "grey92")),
            column(4, textInput("ref_line_col", "Reference line color", "grey40")),
            column(4, textInput("plot_box_col", "Box color", "grey40"))
          ),
          fluidRow(
            column(4, textInput("pathway_label_col", "Pathway label color", "black")),
            column(4, textInput("study_label_col", "Study label color", "grey20")),
            column(4, textInput("axis_text_col", "Axis text color", "black"))
          ),
          fluidRow(
            column(4, textInput("axis_title_col", "Axis title color", "black")),
            column(4, textInput("y_axis_title_col", "Y-axis title color", "black")),
            column(4, textInput("title_col", "Title color", "black"))
          ),
          fluidRow(
            column(3, textInput("panel_label_col", "Panel letter color", "black")),
            column(3, textInput("dot_grid_col", "Panel c grid color", "grey88")),
            column(3, numericInput("dot_grid_linewidth", "Panel c grid linewidth", 0.35, step = 0.05)),
            column(3, textInput("dot_border_col", "Panel c border color", "grey80"))
          ),
          
          h4("Font sizes"),
          fluidRow(
            column(4, numericInput("pathway_label_size_a", "Panel a pathway size", 3.2, step = 0.1)),
            column(4, numericInput("pathway_label_size_b", "Panel b pathway size", 3.0, step = 0.1)),
            column(4, numericInput("pathway_label_size_c", "Panel c pathway size", 3.0, step = 0.1))
          ),
          fluidRow(
            column(4, numericInput("study_label_size", "Study label size", 2.9, step = 0.1)),
            column(4, numericInput("axis_text_size", "Axis tick size", 3.0, step = 0.1)),
            column(4, numericInput("axis_title_size", "Axis title size", 3.2, step = 0.1))
          ),
          fluidRow(
            column(4, numericInput("y_axis_title_size", "Y-axis title size", 3.4, step = 0.1)),
            column(4, numericInput("title_size", "Panel title size", 3.6, step = 0.1)),
            column(4, numericInput("panel_label_size", "Panel letter size", 5, step = 0.5))
          ),
          
          h4("Panel c column labels"),
          fluidRow(
            column(3, numericInput("column_label_size_c", "Size", 2.8, step = 0.1)),
            column(3, numericInput("column_label_angle_c", "Angle", 0, step = 5)),
            column(3, numericInput("column_label_hjust_c", "hjust", 0.5, step = 0.1)),
            column(3, numericInput("column_label_vjust_c", "vjust", 1.0, step = 0.1))
          ),
          
          h4("Panel title text and justification"),
          textInput("title_a", "Panel a title", "CAD meta-analysis"),
          textInput("title_b", "Panel b title", "Imaging outcomes"),
          textInput("title_c", "Panel c title", "CAC-adjusted imaging associations"),
          fluidRow(
            column(6, numericInput("title_hjust", "Title hjust", 0, step = 0.1)),
            column(6, numericInput("title_vjust", "Title vjust", 0.5, step = 0.1))
          ),
          fluidRow(
            column(6, numericInput("panel_label_hjust", "Panel letter hjust", 0, step = 0.1)),
            column(6, numericInput("panel_label_vjust", "Panel letter vjust", 0.5, step = 0.1))
          )
        ),
        
        tabPanel(
          "Layout",
          fluidRow(
            column(6, numericInput("column_width_a", "Column 1 width", 1.25, step = 0.05)),
            column(6, numericInput("column_width_bc", "Column 2 width", 1.00, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("row_height_b", "Panel b height", 1.00, step = 0.05)),
            column(6, numericInput("row_height_c", "Panel c height", 1.00, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("gap_width_ab", "Whitespace between left and right columns", 0.05, step = 0.01)),
            column(6, numericInput("gap_height_bc", "Whitespace between panels b and c", 0.05, step = 0.01))
          )
        ),
        
        tabPanel(
          "Export",
          numericInput("png_width", "PNG width (in)", 13, step = 0.5),
          numericInput("png_height", "PNG height (in)", 9, step = 0.5),
          numericInput("png_dpi", "PNG DPI", 300, min = 72, step = 1),
          textInput("png_bg", "PNG background", "white"),
          downloadButton("download_png", "Download PNG")
        )
      )
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Plot", plotOutput("combined_plot", height = "950px")),
        tabPanel("Meta-analysis preview", DTOutput("meta_table")),
        tabPanel("Imaging preview", DTOutput("imaging_table")),
        tabPanel("Panel C linear preview", DTOutput("dot_linear_table")),
        tabPanel("Panel C logistic preview", DTOutput("dot_logistic_table")),
        tabPanel("Panel C prepared data", DTOutput("dot_matrix_table"))
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {
  SETTINGS_TEXTAREA_IDS <- c("dot_outcome_label_map", "order_a", "order_b", "order_c", "label_map", "color_map")
  SETTINGS_TEXT_IDS <- c("coronary_outcome_label", "carotid_outcome_label", "axis_breaks_a", "xlab_a", "y_axis_title_a", "axis_breaks_b", "xlab_b", "y_axis_title_b", "y_axis_title_c", "nonsig_dot_col", "panel_b_alternating_fill", "legend_b_symbol_col", "legend_c_symbol_col", "dot_legend_title", "tick_guide_col_a", "tick_guide_col_b", "pooled_bg", "ref_line_col", "plot_box_col", "pathway_label_col", "study_label_col", "axis_text_col", "axis_title_col", "y_axis_title_col", "title_col", "panel_label_col", "dot_grid_col", "dot_border_col", "title_a", "title_b", "title_c", "png_bg", "dot_outcome_label_map", "order_a", "order_b", "order_c", "label_map", "color_map")
  SETTINGS_NUMERIC_IDS <- c("pathway_wrap_width", "outcome_wrap_width", "axis_min_a", "axis_max_a", "axis_min_b", "axis_max_b", "pathway_x_a", "study_x_a", "panel_label_x_a", "panel_label_y_a", "panel_label_y_nudge_a", "xlim_left_a", "xlim_right_a", "title_x_a", "stats_or_x_a", "stats_p_x_a", "stats_i2_x_a", "stats_header_size_a", "stats_text_size_a", "y_axis_title_x_a", "y_axis_title_y_a", "title_y_nudge_a", "pathway_x_b", "panel_label_x_b", "title_x_b", "panel_label_y_b", "panel_label_y_nudge_b", "xlim_left_b", "xlim_right_b", "title_y_nudge_b", "y_axis_title_x_b", "y_axis_title_y_b", "pathway_x_c", "panel_label_x_c", "title_x_c", "panel_label_y_c", "xlim_left_c", "xlim_right_c", "y_axis_title_x_c", "meta_group_gap", "imaging_group_gap", "pooled_bar_height", "plot_y_pad_a", "title_y_pad_a", "pathway_y_nudge_a", "plot_y_pad_b", "title_y_pad_b", "pathway_y_nudge_b", "plot_y_pad_c", "title_y_pad_c", "column_label_y_pad_c", "column_label_y_c", "study_point_size", "pooled_point_size", "imaging_point_size", "ci_linewidth", "ref_linewidth", "plot_box_linewidth", "dot_min_size", "dot_max_size", "nonsig_dot_size", "dot_scale_multiplier", "dot_legend_breaks_n", "legend_b_x", "legend_b_y", "legend_b_just_x", "legend_b_just_y", "legend_b_text_size", "legend_c_x", "legend_c_y", "legend_c_just_x", "legend_c_just_y", "legend_c_title_size", "legend_c_text_size", "tick_guide_linewidth_a", "tick_guide_linewidth_b", "dot_grid_linewidth", "pathway_label_size_a", "pathway_label_size_b", "pathway_label_size_c", "study_label_size", "axis_text_size", "axis_title_size", "y_axis_title_size", "title_size", "panel_label_size", "column_label_size_c", "column_label_angle_c", "column_label_hjust_c", "column_label_vjust_c", "title_hjust", "title_vjust", "panel_label_hjust", "panel_label_vjust", "column_width_a", "column_width_bc", "row_height_b", "row_height_c", "gap_width_ab", "gap_height_bc", "png_width", "png_height", "png_dpi")
  SETTINGS_CHECKBOX_IDS <- c("use_example", "show_dot_grid", "show_tick_guides_a", "show_tick_guides_b")
  SETTINGS_SELECT_IDS <- c("study_shape", "pooled_shape", "panel_b_shape_1", "panel_b_shape_2", "tick_guide_linetype_a", "tick_guide_linetype_b", "meta_name", "meta_model", "coronary_outcome", "carotid_outcome", "dot_linear_outcomes", "dot_logistic_outcomes")
  SETTINGS_TEXT_IDS <- setdiff(SETTINGS_TEXT_IDS, SETTINGS_TEXTAREA_IDS)
  SETTINGS_IDS <- unique(c(SETTINGS_TEXT_IDS, SETTINGS_TEXTAREA_IDS, SETTINGS_NUMERIC_IDS, SETTINGS_CHECKBOX_IDS, SETTINGS_SELECT_IDS))

  # Embedded defaults from figure_3_shiny_settings_2026-08-03(1).json.
  # The legacy stats_p_x_a setting now controls the 95% CI column x-position.
  DEFAULT_SETTINGS <- list(
    coronary_outcome_label = "Coronary plaque",
    carotid_outcome_label = "Carotid plaque",
    axis_breaks_a = "1.0, 1.2, 1.5, 2.0, 3.0",
    xlab_a = "Odds ratio for CAD",
    y_axis_title_a = "Process-specific PRS",
    axis_breaks_b = "1.0, 1.2, 1.5, 2.0",
    xlab_b = "Odds ratio for plaque presence",
    y_axis_title_b = "Process-specific PRS",
    y_axis_title_c = "Process-specific PRS",
    nonsig_dot_col = "grey40",
    panel_b_alternating_fill = "grey95",
    legend_b_symbol_col = "grey60",
    legend_c_symbol_col = "grey60",
    dot_legend_title = "Regression coefficient",
    tick_guide_col_a = "grey88",
    tick_guide_col_b = "grey88",
    pooled_bg = "grey92",
    ref_line_col = "grey40",
    plot_box_col = "grey40",
    pathway_label_col = "black",
    study_label_col = "grey20",
    axis_text_col = "black",
    axis_title_col = "black",
    y_axis_title_col = "black",
    title_col = "black",
    panel_label_col = "black",
    dot_grid_col = "grey88",
    dot_border_col = "grey40",
    title_a = "CAD meta-analysis",
    title_b = "Imaging outcomes",
    title_c = "CAC-adjusted imaging associations",
    png_bg = "white",
    dot_outcome_label_map = "ct_cac_tot = Coronary calcium score\nct_cac_tot_log = Coronary calcium score\nct_sis = Segment involvement score\nct_50_stenosis_any = Stenosis ≥50%\nct_athero_noncalc = Non-calcified plaque\nct_av_calcification = Aortic valve calcification",
    order_a = "chol_prs_std\ntrl_prs_std\nlpa_prs_std\nmdwo_prs_std\nobesity_prs_std\nendothelial_prs_std\ncoagulation_prs_std\ncardiac_prs_std",
    order_b = "chol_prs_std\ntrl_prs_std\nlpa_prs_std\nmdwo_prs_std\nobesity_prs_std\nendothelial_prs_std\ncoagulation_prs_std\ncardiac_prs_std",
    order_c = "chol_prs_std\ntrl_prs_std\nlpa_prs_std\nmdwo_prs_std\nobesity_prs_std\nendothelial_prs_std\ncoagulation_prs_std\ncardiac_prs_std",
    label_map = "chol_prs_std = ChEL\ntrl_prs_std = TRL\nlpa_prs_std = Lp(a)\nmdwo_prs_std = MUNW\nobesity_prs_std = Obesity\nendothelial_prs_std = Endothelial dysfunction\ncardiac_prs_std = Cardiac remodeling\ncoagulation_prs_std = Coagulation",
    color_map = "chol_prs_std = #B22D4C\ntrl_prs_std = #E07A8E\nlpa_prs_std = #7A1F3A\nmdwo_prs_std = #4A7C4E\nobesity_prs_std = #95C075\nendothelial_prs_std = #2D5F8B\ncoagulation_prs_std = #7AAED1\ncardiac_prs_std = #5D4787",
    pathway_wrap_width = 16,
    outcome_wrap_width = 18,
    axis_min_a = 0.9,
    axis_max_a = 3.2,
    axis_min_b = 0.94,
    axis_max_b = 2.1,
    pathway_x_a = -0.89,
    study_x_a = -0.34,
    panel_label_x_a = -1.1,
    panel_label_y_a = NA,
    panel_label_y_nudge_a = 0,
    xlim_left_a = -1,
    xlim_right_a = 0.8,
    title_x_a = 0,
    stats_or_x_a = 1.18,
    stats_p_x_a = 1.58,
    stats_i2_x_a = 1.98,
    stats_header_size_a = 2.9,
    stats_text_size_a = 2.75,
    y_axis_title_x_a = -1.1,
    y_axis_title_y_a = NA,
    title_y_nudge_a = 0,
    pathway_x_b = -0.52,
    panel_label_x_b = -0.82,
    title_x_b = 0,
    panel_label_y_b = 11.3,
    panel_label_y_nudge_b = 0,
    xlim_left_b = -0.3,
    xlim_right_b = 1.04,
    title_y_nudge_b = 0,
    y_axis_title_x_b = -0.76,
    y_axis_title_y_b = NA,
    pathway_x_c = -1.32,
    panel_label_x_c = -2.8,
    title_x_c = 1,
    panel_label_y_c = NA,
    xlim_left_c = -0.4,
    xlim_right_c = NA,
    y_axis_title_x_c = -2.5,
    meta_group_gap = 0.65,
    imaging_group_gap = 0.26,
    pooled_bar_height = 0.88,
    plot_y_pad_a = 0.55,
    title_y_pad_a = 1.35,
    pathway_y_nudge_a = 0,
    plot_y_pad_b = 0.45,
    title_y_pad_b = 1.35,
    pathway_y_nudge_b = 0,
    plot_y_pad_c = 0.5,
    title_y_pad_c = 1.15,
    column_label_y_pad_c = -0.55,
    column_label_y_c = NA,
    study_point_size = 2.2,
    pooled_point_size = 2.8,
    imaging_point_size = 2.4,
    ci_linewidth = 0.5,
    ref_linewidth = 0.35,
    plot_box_linewidth = 0.45,
    dot_min_size = 1.2,
    dot_max_size = 8,
    nonsig_dot_size = 1.4,
    dot_scale_multiplier = 1,
    dot_legend_breaks_n = 3,
    legend_b_x = 0.05,
    legend_b_y = 0.87,
    legend_b_just_x = 0,
    legend_b_just_y = 0,
    legend_b_text_size = 8,
    legend_c_x = 0.45,
    legend_c_y = 0.85,
    legend_c_just_x = 0.5,
    legend_c_just_y = 0,
    legend_c_title_size = 8,
    legend_c_text_size = 8,
    tick_guide_linewidth_a = 0.25,
    tick_guide_linewidth_b = 0.25,
    dot_grid_linewidth = 0.35,
    pathway_label_size_a = 3.2,
    pathway_label_size_b = 3,
    pathway_label_size_c = 3,
    study_label_size = 2.9,
    axis_text_size = 3,
    axis_title_size = 3.2,
    y_axis_title_size = 3.4,
    title_size = 3.6,
    panel_label_size = 4,
    column_label_size_c = 2.8,
    column_label_angle_c = 40,
    column_label_hjust_c = 0.85,
    column_label_vjust_c = 0,
    title_hjust = 0,
    title_vjust = 0.5,
    panel_label_hjust = 0,
    panel_label_vjust = 0.5,
    column_width_a = 1.25,
    column_width_bc = 1,
    row_height_b = 1,
    row_height_c = 1,
    gap_width_ab = 0.05,
    gap_height_bc = -0.1,
    png_width = 10,
    png_height = 9,
    png_dpi = 1200,
    use_example = FALSE,
    show_dot_grid = TRUE,
    show_tick_guides_a = TRUE,
    show_tick_guides_b = TRUE,
    study_shape = "21",
    pooled_shape = "23",
    panel_b_shape_1 = "22",
    panel_b_shape_2 = "24",
    tick_guide_linetype_a = "solid",
    tick_guide_linetype_b = "solid",
    meta_name = "CAD broad",
    meta_model = "Fixed effect",
    coronary_outcome = "ct_athero_any",
    carotid_outcome = "us_carotid_plaque_presence",
    dot_linear_outcomes = c("ct_cac_tot_log", "ct_sis"),
    dot_logistic_outcomes = c("ct_50_stenosis_any", "ct_athero_noncalc", "ct_av_calcification")
  )

  apply_settings <- function(settings, notify = FALSE) {
    for (id in intersect(names(settings), SETTINGS_TEXT_IDS))
      updateTextInput(session, id, value = as.character(settings[[id]] %||% ""))
    for (id in intersect(names(settings), SETTINGS_TEXTAREA_IDS))
      updateTextAreaInput(session, id, value = as.character(settings[[id]] %||% ""))
    for (id in intersect(names(settings), SETTINGS_NUMERIC_IDS)) {
      value <- suppressWarnings(as.numeric(settings[[id]]))
      updateNumericInput(session, id, value = if (length(value)) value[1] else NA_real_)
    }
    for (id in intersect(names(settings), SETTINGS_CHECKBOX_IDS))
      updateCheckboxInput(session, id, value = isTRUE(settings[[id]]))
    select_settings <- settings[intersect(names(settings), SETTINGS_SELECT_IDS)]
    session$onFlushed(function() {
      for (id in names(select_settings)) {
        value <- select_settings[[id]]
        if (!is.null(value)) updateSelectInput(session, id, selected = value)
      }
    }, once = TRUE)
    if (isTRUE(notify)) showNotification("Settings imported from JSON.", type = "message")
  }

  session$onFlushed(function() apply_settings(DEFAULT_SETTINGS, notify = FALSE), once = TRUE)

  output$download_settings <- downloadHandler(
    filename = function() paste0("figure_3_shiny_settings_", Sys.Date(), ".json"),
    content = function(file) {
      values <- reactiveValuesToList(input)
      values <- values[intersect(SETTINGS_IDS, names(values))]
      jsonlite::write_json(
        list(version = 1, exported_at = format(Sys.time(), tz = "UTC"), settings = values),
        path = file, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"
      )
    }
  )

  observeEvent(input$settings_file, {
    req(input$settings_file$datapath)
    imported <- tryCatch(
      jsonlite::read_json(input$settings_file$datapath, simplifyVector = TRUE),
      error = function(e) NULL
    )
    if (is.null(imported)) {
      showNotification("The selected file is not valid JSON.", type = "error")
      return()
    }
    settings <- imported$settings %||% imported
    if (!is.list(settings)) {
      showNotification("The JSON file does not contain a settings object.", type = "error")
      return()
    }
    apply_settings(settings, notify = TRUE)
  }, ignoreInit = TRUE)

  meta_data <- reactive({
    if (isTRUE(input$use_example)) {
      return(generate_example_meta_data())
    }
    
    shiny::validate(
      shiny::need(!is.null(input$meta_file), "Please upload a meta-analysis CSV or use example data.")
    )
    
    safe_read_csv(input$meta_file$datapath)
  })
  
  imaging_data <- reactive({
    if (isTRUE(input$use_example)) {
      return(generate_example_imaging_data())
    }
    
    shiny::validate(
      shiny::need(!is.null(input$imaging_file), "Please upload a SCAPIS imaging CSV or use example data.")
    )
    
    safe_read_csv(input$imaging_file$datapath)
  })
  
  dot_linear_data <- reactive({
    if (isTRUE(input$use_example)) {
      return(generate_example_dot_linear_data())
    }
    
    shiny::validate(
      shiny::need(
        !is.null(input$dot_linear_file),
        "Please upload results_imaging_linear_cac_adjusted_scapis.csv or use example data."
      )
    )
    
    safe_read_csv(input$dot_linear_file$datapath)
  })
  
  dot_linear_non_cac_data <- reactive({
    if (isTRUE(input$use_example)) {
      return(generate_example_dot_linear_non_cac_adjusted_data())
    }
    
    shiny::validate(
      shiny::need(
        !is.null(input$dot_linear_non_cac_file),
        "Please upload results_imaging_linear.csv (non-CAC-adjusted) for ct_cac_tot_log or use example data."
      )
    )
    
    safe_read_csv(input$dot_linear_non_cac_file$datapath)
  })
  
  dot_logistic_data <- reactive({
    if (isTRUE(input$use_example)) {
      return(generate_example_dot_logistic_data())
    }
    
    shiny::validate(
      shiny::need(
        !is.null(input$dot_logistic_file),
        "Please upload results_imaging_logistic_cac_adjusted_scapis.csv or use example data."
      )
    )
    
    safe_read_csv(input$dot_logistic_file$datapath)
  })
  
  predictor_labels <- reactive({
    x <- parse_named_mapping(input$label_map)
    if (!length(x)) DEFAULT_PREDICTOR_LABELS else x
  })
  
  predictor_colors <- reactive({
    x <- parse_named_mapping(input$color_map)
    if (!length(x)) DEFAULT_PREDICTOR_COLORS else x
  })
  
  order_a <- reactive({
    x <- parse_lines(input$order_a)
    shiny::validate(shiny::need(length(x) > 0, "Panel a PRS order is empty."))
    intersect(x, KEEP_PREDICTORS)
  })
  
  order_b <- reactive({
    x <- parse_lines(input$order_b)
    shiny::validate(shiny::need(length(x) > 0, "Panel b PRS order is empty."))
    intersect(x, KEEP_PREDICTORS)
  })
  
  order_c <- reactive({
    x <- parse_lines(input$order_c)
    shiny::validate(shiny::need(length(x) > 0, "Panel c PRS order is empty."))
    intersect(x, KEEP_PREDICTORS)
  })
  
  output$meta_name_ui <- renderUI({
    dat <- meta_data()
    
    shiny::validate(
      shiny::need("meta_name" %in% names(dat), "meta_name column not found.")
    )
    
    choices <- sort(unique(trimws(as.character(dat$meta_name))))
    
    selectInput(
      "meta_name",
      "Meta-analysis name",
      choices = choices,
      selected = choices[1]
    )
  })
  
  observe({
    dat <- meta_data()
    
    shiny::validate(
      shiny::need("meta_model" %in% names(dat), "meta_model column not found.")
    )
    
    choices <- trimws(as.character(dat$meta_model))
    choices <- choices[!is.na(choices) & nzchar(choices)]
    choices <- sort(unique(choices))
    
    if (!length(choices)) {
      choices <- c("Fixed effect", "Random effects")
    }
    
    current <- isolate(input$meta_model)
    preferred <- if (!is.null(current) && current %in% choices) {
      current
    } else if ("Fixed effect" %in% choices) {
      "Fixed effect"
    } else {
      choices[1]
    }
    
    updateSelectInput(
      session,
      "meta_model",
      choices = choices,
      selected = preferred
    )
  })
  
  output$coronary_outcome_ui <- renderUI({
    dat <- imaging_data()

    if (!"outcome_short" %in% names(dat)) {
      return(tags$div(
        class = "text-danger",
        "outcome_short column not found in imaging data."
      ))
    }

    outcomes <- trimws(as.character(dat$outcome_short))
    outcomes <- sort(unique(outcomes[!is.na(outcomes) & nzchar(outcomes)]))

    if (!length(outcomes)) {
      return(tags$div(
        class = "text-danger",
        "No non-missing outcome_short values were found in imaging data."
      ))
    }

    selected_outcome <- if ("ct_athero_any" %in% outcomes) {
      "ct_athero_any"
    } else {
      outcomes[[1L]]
    }

    selectInput(
      "coronary_outcome",
      "Panel b outcome, symbol 1",
      choices = outcomes,
      selected = selected_outcome
    )
  })

  output$carotid_outcome_ui <- renderUI({
    dat <- imaging_data()

    if (!"outcome_short" %in% names(dat)) {
      return(tags$div(
        class = "text-danger",
        "outcome_short column not found in imaging data."
      ))
    }

    outcomes <- trimws(as.character(dat$outcome_short))
    outcomes <- sort(unique(outcomes[!is.na(outcomes) & nzchar(outcomes)]))

    if (!length(outcomes)) {
      return(tags$div(
        class = "text-danger",
        "No non-missing outcome_short values were found in imaging data."
      ))
    }

    selected_outcome <- if ("us_carotid_plaque_presence" %in% outcomes) {
      "us_carotid_plaque_presence"
    } else if (length(outcomes) >= 2L) {
      outcomes[[2L]]
    } else {
      outcomes[[1L]]
    }

    selectInput(
      "carotid_outcome",
      "Panel b outcome, symbol 2",
      choices = outcomes,
      selected = selected_outcome
    )
  })
  
  output$dot_linear_outcomes_ui <- renderUI({
    adjusted_dat <- dot_linear_data()
    non_cac_dat <- dot_linear_non_cac_data()
    
    shiny::validate(
      shiny::need("outcome_short" %in% names(adjusted_dat), "outcome_short column not found in CAC-adjusted linear panel C data."),
      shiny::need("outcome_short" %in% names(non_cac_dat), "outcome_short column not found in non-CAC-adjusted linear panel C data.")
    )
    
    adjusted_outcomes <- trimws(as.character(adjusted_dat$outcome_short))
    adjusted_outcomes <- adjusted_outcomes[!is.na(adjusted_outcomes) & nzchar(adjusted_outcomes)]
    adjusted_outcomes <- setdiff(adjusted_outcomes, PANEL_C_NON_CAC_ADJUSTED_OUTCOME)
    
    non_cac_outcomes <- trimws(as.character(non_cac_dat$outcome_short))
    non_cac_outcomes <- non_cac_outcomes[!is.na(non_cac_outcomes) & nzchar(non_cac_outcomes)]
    non_cac_outcomes <- intersect(non_cac_outcomes, PANEL_C_NON_CAC_ADJUSTED_OUTCOME)
    
    outcomes <- unique(c(non_cac_outcomes, sort(unique(adjusted_outcomes))))
    shiny::validate(
      shiny::need(length(outcomes) > 0, "No valid linear outcomes were found for panel C."),
      shiny::need(
        PANEL_C_NON_CAC_ADJUSTED_OUTCOME %in% outcomes,
        paste0(PANEL_C_NON_CAC_ADJUSTED_OUTCOME, " was not found in the non-CAC-adjusted linear file.")
      )
    )
    
    default <- intersect(DEFAULT_DOT_LINEAR_OUTCOMES, outcomes)
    if (!length(default)) default <- outcomes
    
    selectInput(
      "dot_linear_outcomes",
      "Panel c linear columns",
      choices = outcomes,
      selected = default,
      multiple = TRUE
    )
  })
  
  output$dot_logistic_outcomes_ui <- renderUI({
    dat <- dot_logistic_data()
    
    shiny::validate(
      shiny::need("outcome_short" %in% names(dat), "outcome_short column not found in logistic panel C data.")
    )
    
    outcomes <- sort(unique(trimws(as.character(dat$outcome_short))))
    default <- intersect(DEFAULT_DOT_LOGISTIC_OUTCOMES, outcomes)
    if (!length(default)) default <- outcomes
    
    selectInput(
      "dot_logistic_outcomes",
      "Panel c logistic columns",
      choices = outcomes,
      selected = default,
      multiple = TRUE
    )
  })
  
  axis_breaks_a <- reactive({
    vals <- parse_breaks(input$axis_breaks_a)
    vals <- vals[is.finite(vals) & vals > 0]
    shiny::validate(shiny::need(length(vals) > 0, "Panel a needs at least one positive x-axis break."))
    vals
  })
  
  axis_breaks_b <- reactive({
    vals <- parse_breaks(input$axis_breaks_b)
    vals <- vals[is.finite(vals) & vals > 0]
    shiny::validate(shiny::need(length(vals) > 0, "Panel b needs at least one positive x-axis break."))
    vals
  })
  
  filtered_meta_plot_data <- reactive({
    req(input$meta_name, input$meta_model)
    
    prepare_meta_plot_data(
      dat = meta_data(),
      meta_name_selected = input$meta_name,
      meta_model_selected = input$meta_model,
      predictor_order = order_a(),
      predictor_labels = predictor_labels(),
      predictor_colors = predictor_colors(),
      study_order = DEFAULT_STUDY_ORDER,
      group_gap = input$meta_group_gap,
      pathway_wrap_width = input$pathway_wrap_width
    )
  })
  
  filtered_combined_imaging_data <- reactive({
    req(input$coronary_outcome, input$carotid_outcome)
    
    prepare_combined_imaging_forest_data(
      dat = imaging_data(),
      outcome_b = input$coronary_outcome,
      outcome_c = input$carotid_outcome,
      outcome_b_label = input$coronary_outcome_label,
      outcome_c_label = input$carotid_outcome_label,
      predictor_order = order_b(),
      predictor_labels = predictor_labels(),
      predictor_colors = predictor_colors(),
      group_gap = input$imaging_group_gap,
      pathway_wrap_width = input$pathway_wrap_width
    )
  })
  
  filtered_dot_matrix_data <- reactive({
    req(input$dot_linear_outcomes, input$dot_logistic_outcomes)
    
    adjusted_linear <- dot_linear_data()
    non_cac_linear <- dot_linear_non_cac_data()
    
    shiny::validate(
      shiny::need("outcome_short" %in% names(adjusted_linear), "outcome_short column not found in CAC-adjusted linear data."),
      shiny::need("outcome_short" %in% names(non_cac_linear), "outcome_short column not found in non-CAC-adjusted linear data.")
    )
    
    # ct_cac_tot_log must come exclusively from the non-CAC-adjusted source.
    adjusted_linear <- adjusted_linear %>%
      filter(trimws(as.character(outcome_short)) != PANEL_C_NON_CAC_ADJUSTED_OUTCOME)
    non_cac_linear <- non_cac_linear %>%
      filter(trimws(as.character(outcome_short)) == PANEL_C_NON_CAC_ADJUSTED_OUTCOME)
    
    shiny::validate(
      shiny::need(
        nrow(non_cac_linear) > 0,
        paste0(PANEL_C_NON_CAC_ADJUSTED_OUTCOME, " was not found in the non-CAC-adjusted linear file.")
      )
    )
    
    combined_linear <- bind_rows(adjusted_linear, non_cac_linear)
    
    prepare_dot_matrix_data(
      linear_dat = combined_linear,
      logistic_dat = dot_logistic_data(),
      linear_outcomes = input$dot_linear_outcomes,
      logistic_outcomes = input$dot_logistic_outcomes,
      predictor_order = order_c(),
      predictor_labels = predictor_labels(),
      predictor_colors = predictor_colors(),
      dot_outcome_labels = parse_named_mapping(input$dot_outcome_label_map),
      pathway_wrap_width = input$pathway_wrap_width,
      outcome_wrap_width = input$outcome_wrap_width
    )
  })
  
  panel_a <- reactive({
    shiny::validate(
      shiny::need(input$axis_min_a > 0 && input$axis_max_a > input$axis_min_a, "Panel a axis limits are invalid.")
    )
    
    plot_meta_forest(
      plot_data = filtered_meta_plot_data(),
      title = "",
      panel_label = "a",
      y_axis_title = input$y_axis_title_a,
      axis_min = input$axis_min_a,
      axis_max = input$axis_max_a,
      axis_breaks = axis_breaks_a(),
      pathway_x = input$pathway_x_a,
      study_x = input$study_x_a,
      y_axis_title_x = input$y_axis_title_x_a,
      y_axis_title_y = input$y_axis_title_y_a,
      xlim_left = input$xlim_left_a,
      xlim_right = input$xlim_right_a,
      study_point_size = input$study_point_size,
      pooled_point_size = input$pooled_point_size,
      study_shape = as.numeric(input$study_shape),
      pooled_shape = as.numeric(input$pooled_shape),
      ci_linewidth = input$ci_linewidth,
      ref_linewidth = input$ref_linewidth,
      pooled_bg = input$pooled_bg,
      pathway_label_size = input$pathway_label_size_a,
      study_label_size = input$study_label_size,
      axis_text_size = input$axis_text_size,
      axis_title_size = input$axis_title_size,
      y_axis_title_size = input$y_axis_title_size,
      title_size = input$title_size,
      panel_label_size = input$panel_label_size,
      xlab = input$xlab_a,
      pathway_label_col = input$pathway_label_col,
      study_label_col = input$study_label_col,
      title_col = input$title_col,
      axis_text_col = input$axis_text_col,
      axis_title_col = input$axis_title_col,
      y_axis_title_col = input$y_axis_title_col,
      panel_label_col = input$panel_label_col,
      ref_line_col = input$ref_line_col,
      plot_box_col = input$plot_box_col,
      plot_box_linewidth = input$plot_box_linewidth,
      pooled_bar_height = input$pooled_bar_height,
      plot_y_pad = input$plot_y_pad_a,
      title_y_pad = input$title_y_pad_a,
      pathway_y_nudge = input$pathway_y_nudge_a,
      title_x = input$title_x_a,
      title_y_nudge = input$title_y_nudge_a,
      title_hjust = input$title_hjust,
      title_vjust = input$title_vjust,
      panel_label_x = input$panel_label_x_a,
      panel_label_y = input$panel_label_y_a,
      panel_label_y_nudge = input$panel_label_y_nudge_a,
      panel_label_hjust = input$panel_label_hjust,
      panel_label_vjust = input$panel_label_vjust,
      show_tick_guides = input$show_tick_guides_a,
      tick_guide_col = input$tick_guide_col_a,
      tick_guide_linewidth = input$tick_guide_linewidth_a,
      tick_guide_linetype = input$tick_guide_linetype_a,
      stats_or_x = input$stats_or_x_a,
      stats_p_x = input$stats_p_x_a,
      stats_i2_x = input$stats_i2_x_a,
      stats_header_size = input$stats_header_size_a,
      stats_text_size = input$stats_text_size_a
    )
  })
  
  panel_b <- reactive({
    shiny::validate(
      shiny::need(input$axis_min_b > 0 && input$axis_max_b > input$axis_min_b, "Panel b axis limits are invalid.")
    )
    
    plot_combined_imaging_forest(
      rows = filtered_combined_imaging_data(),
      title = "",
      panel_label = "b",
      y_axis_title = input$y_axis_title_b,
      axis_min = input$axis_min_b,
      axis_max = input$axis_max_b,
      axis_breaks = axis_breaks_b(),
      pathway_x = input$pathway_x_b,
      y_axis_title_x = input$y_axis_title_x_b,
      y_axis_title_y = input$y_axis_title_y_b,
      xlim_left = input$xlim_left_b,
      xlim_right = input$xlim_right_b,
      point_size = input$imaging_point_size,
      outcome_shape_1 = as.numeric(input$panel_b_shape_1),
      outcome_shape_2 = as.numeric(input$panel_b_shape_2),
      ci_linewidth = input$ci_linewidth,
      ref_linewidth = input$ref_linewidth,
      pathway_label_size = input$pathway_label_size_b,
      axis_text_size = input$axis_text_size,
      axis_title_size = input$axis_title_size,
      y_axis_title_size = input$y_axis_title_size,
      title_size = input$title_size,
      panel_label_size = input$panel_label_size,
      xlab = input$xlab_b,
      pathway_label_col = input$pathway_label_col,
      title_col = input$title_col,
      axis_text_col = input$axis_text_col,
      axis_title_col = input$axis_title_col,
      y_axis_title_col = input$y_axis_title_col,
      panel_label_col = input$panel_label_col,
      ref_line_col = input$ref_line_col,
      plot_box_col = input$plot_box_col,
      plot_box_linewidth = input$plot_box_linewidth,
      plot_y_pad = input$plot_y_pad_b,
      title_y_pad = input$title_y_pad_b,
      pathway_y_nudge = input$pathway_y_nudge_b,
      title_x = input$title_x_b,
      title_y_nudge = input$title_y_nudge_b,
      title_hjust = input$title_hjust,
      title_vjust = input$title_vjust,
      panel_label_x = input$panel_label_x_b,
      panel_label_y = input$panel_label_y_b,
      panel_label_y_nudge = input$panel_label_y_nudge_b,
      panel_label_hjust = input$panel_label_hjust,
      panel_label_vjust = input$panel_label_vjust,
      legend_x = input$legend_b_x,
      legend_y = input$legend_b_y,
      legend_just_x = input$legend_b_just_x,
      legend_just_y = input$legend_b_just_y,
      legend_text_size = input$legend_b_text_size,
      legend_symbol_col = input$legend_b_symbol_col,
      show_tick_guides = input$show_tick_guides_b,
      tick_guide_col = input$tick_guide_col_b,
      tick_guide_linewidth = input$tick_guide_linewidth_b,
      tick_guide_linetype = input$tick_guide_linetype_b,
      alternating_fill = input$panel_b_alternating_fill
    )
  })
  
  panel_c <- reactive({
    plot_dot_matrix(
      rows = filtered_dot_matrix_data(),
      title = "",
      panel_label = "c",
      y_axis_title = input$y_axis_title_c,
      pathway_x = input$pathway_x_c,
      y_axis_title_x = input$y_axis_title_x_c,
      xlim_left = input$xlim_left_c,
      xlim_right = input$xlim_right_c,
      plot_y_pad = input$plot_y_pad_c,
      title_y_pad = input$title_y_pad_c,
      x_label_y_pad = input$column_label_y_pad_c,
      column_label_y = input$column_label_y_c,
      pathway_label_size = input$pathway_label_size_c,
      column_label_size = input$column_label_size_c,
      column_label_angle = input$column_label_angle_c,
      column_label_hjust = input$column_label_hjust_c,
      column_label_vjust = input$column_label_vjust_c,
      y_axis_title_size = input$y_axis_title_size,
      title_size = input$title_size,
      panel_label_size = input$panel_label_size,
      panel_label_x = input$panel_label_x_c,
      panel_label_y = input$panel_label_y_c,
      title_x = input$title_x_c,
      title_hjust = input$title_hjust,
      title_vjust = input$title_vjust,
      panel_label_hjust = input$panel_label_hjust,
      panel_label_vjust = input$panel_label_vjust,
      dot_min_size = input$dot_min_size,
      dot_max_size = input$dot_max_size,
      dot_scale_multiplier = input$dot_scale_multiplier,
      dot_legend_breaks_n = input$dot_legend_breaks_n,
      dot_legend_title = input$dot_legend_title,
      nonsig_dot_size = input$nonsig_dot_size,
      nonsig_dot_col = input$nonsig_dot_col,
      show_grid = input$show_dot_grid,
      grid_col = input$dot_grid_col,
      grid_linewidth = input$dot_grid_linewidth,
      border_col = input$dot_border_col,
      border_linewidth = input$plot_box_linewidth,
      legend_x = input$legend_c_x,
      legend_y = input$legend_c_y,
      legend_just_x = input$legend_c_just_x,
      legend_just_y = input$legend_c_just_y,
      legend_title_size = input$legend_c_title_size,
      legend_text_size = input$legend_c_text_size,
      legend_symbol_col = input$legend_c_symbol_col,
      pathway_label_col = input$pathway_label_col,
      y_axis_title_col = input$y_axis_title_col,
      title_col = input$title_col,
      panel_label_col = input$panel_label_col,
      column_label_col = input$axis_text_col
    )
  })
  
  combined_plot <- reactive({
    combine_three_panels(
      plot_a = panel_a(),
      plot_b = panel_b(),
      plot_c = panel_c(),
      column_width_a = input$column_width_a,
      column_width_bc = input$column_width_bc,
      row_height_b = input$row_height_b,
      row_height_c = input$row_height_c,
      gap_width_ab = input$gap_width_ab,
      gap_height_bc = input$gap_height_bc
    )
  })
  
  output$combined_plot <- renderPlot({
    print(combined_plot())
  }, bg = "white")
  
  output$meta_table <- renderDT({
    datatable(head(meta_data(), 300), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$imaging_table <- renderDT({
    datatable(head(imaging_data(), 300), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$dot_linear_table <- renderDT({
    datatable(head(dot_linear_data(), 300), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$dot_logistic_table <- renderDT({
    datatable(head(dot_logistic_data(), 300), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$dot_matrix_table <- renderDT({
    datatable(head(filtered_dot_matrix_data(), 300), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$download_png <- downloadHandler(
    filename = function() {
      paste0("three_panel_prs_forest_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = combined_plot(),
        width = input$png_width,
        height = input$png_height,
        dpi = input$png_dpi,
        bg = input$png_bg
      )
    }
  )
}

shinyApp(ui = ui, server = server)