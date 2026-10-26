# ==============================================================================
# Shiny App: Circular PRS Forest Plot Builder
# ==============================================================================

options(shiny.maxRequestSize = 1024 * 1024^2)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(DT)
library(ggnewscale)

# ==============================================================================
# Example data
# ==============================================================================

generate_example_data <- function() {
  outcome_labels <- c(
    as_ukb_ind = "Aortic stenosis",
    vhd_ukb_ind = "Valvular heart disease",
    vte_ukb_ind = "Venous thromboembolism",
    arrhythmia_non_af_ukb_ind = "Non-AF arrhythmia",
    af_ukb_ind = "Atrial fibrillation",
    hf_ukb_ind = "Heart failure",
    aaa_ukb_ind = "Abdominal aortic aneurysm",
    pad_ukb_ind = "Peripheral artery disease",
    is_ukb_ind = "Ischemic stroke",
    cvdmort_ukb_death = "Cardiovascular mortality",
    acmort_ukb_death = "All-cause mortality",
    dm2_ukb_ind = "Type 2 diabetes",
    nafld_ukb_ind = "Non-alcoholic fatty liver disease",
    biliary_pancreatic_ukb_ind = "Biliary and pancreatic disease",
    ckd_ukb_ind = "Chronic kidney disease",
    anemia_ukb_ind = "Anemia",
    thyroid_disease_ukb_ind = "Thyroid disease",
    copd_ukb_ind = "Chronic obstructive pulmonary disease",
    ad_ukb_ind = "Alzheimer's disease"
  )
  
  predictor_labels <- c(
    chol_prs_std = "ChEL",
    trl_prs_std = "TRL",
    lpa_prs_std = "Lp(a)",
    mdwo_prs_std = "MUNW",
    obesity_prs_std = "Obesity",
    endothelial_prs_std = "Endothelial dysfunction",
    coagulation_prs_std = "Coagulation",
    cardiac_prs_std = "Cardiac remodeling"
  )
  
  predictor_colors <- c(
    chol_prs_std = "#B22D4C",
    trl_prs_std = "#E07A8E",
    lpa_prs_std = "#7A1F3A",
    mdwo_prs_std = "#4A7C4E",
    obesity_prs_std = "#95C075",
    endothelial_prs_std = "#2D5F8B",
    coagulation_prs_std = "#7AAED1",
    cardiac_prs_std = "#5D4787"
  )
  
  north <- c(
    "as_ukb_ind", "vhd_ukb_ind", "vte_ukb_ind",
    "arrhythmia_non_af_ukb_ind", "af_ukb_ind", "hf_ukb_ind",
    "aaa_ukb_ind", "pad_ukb_ind", "is_ukb_ind", "cvdmort_ukb_death"
  )
  
  south <- c(
    "acmort_ukb_death", "dm2_ukb_ind", "nafld_ukb_ind",
    "biliary_pancreatic_ukb_ind", "ckd_ukb_ind", "anemia_ukb_ind",
    "thyroid_disease_ukb_ind", "copd_ukb_ind", "ad_ukb_ind"
  )
  
  predictors <- names(predictor_labels)
  
  set.seed(123)
  dat <- tidyr::expand_grid(
    study = "UKB",
    outcome_short = c(north, south),
    predictor_short = predictors
  ) %>%
    mutate(
      outcome_long = unname(outcome_labels[outcome_short]),
      predictor_long = unname(predictor_labels[predictor_short]),
      model = "Model 1: Age + Sex + PCs",
      exclusion = "Main analysis (only Whites)",
      n = 458905,
      n_events = sample(2500:40000, n(), replace = TRUE),
      log_or = rnorm(n(), mean = 0.08, sd = 0.07),
      se = runif(n(), min = 0.02, max = 0.06),
      or = exp(log_or),
      lb = exp(log_or - 1.96 * se),
      ub = exp(log_or + 1.96 * se),
      p_value = 2 * pnorm(-abs(log_or / se)),
      p_bonferroni = pmin(p_value * n(), 1)
    )
  
  list(
    data = dat,
    outcome_labels = outcome_labels,
    predictor_labels = predictor_labels,
    predictor_colors = predictor_colors,
    north = north,
    south = south,
    predictors = predictors
  )
}

# ==============================================================================
# Helpers
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

deg2rad <- function(x) x * pi / 180

safe_read_csv <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "gz") {
    read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

parse_lines <- function(x) {
  x <- trimws(x)
  x <- x[nzchar(x)]
  x
}

parse_named_mapping <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(c())
  lines <- unlist(strsplit(text, "\n"))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  vals <- lapply(lines, function(ln) {
    parts <- unlist(strsplit(ln, "=", fixed = TRUE))
    if (length(parts) < 2) return(NULL)
    c(key = trimws(parts[1]), value = trimws(paste(parts[-1], collapse = "=")))
  })
  vals <- vals[!vapply(vals, is.null, logical(1))]
  if (length(vals) == 0) return(c())
  setNames(vapply(vals, function(x) x[["value"]], character(1)),
           vapply(vals, function(x) x[["key"]], character(1)))
}

make_default_palette <- function(n) {
  base <- c(
    "#C4510B", "#E09000", "#D8AB1B", "#F7E798",
    "#96D2EA", "#3D9AC4", "#2070B0", "#2177B5",
    "#6A3D9A", "#B15928", "#1B9E77", "#D95F02"
  )
  if (n <= length(base)) base[seq_len(n)] else colorRampPalette(base)(n)
}

make_predictor_color_inputs <- function(ids, cols) {
  tagList(
    lapply(seq_along(ids), function(i) {
      fluidRow(
        column(7, div(style = "padding-top:8px;", ids[i])),
        column(5, textInput(paste0("pred_col_", i), NULL, value = cols[i]))
      )
    })
  )
}

make_arc_polygon <- function(start_deg, end_deg, r_inner, r_outer, n = 180) {
  ang <- seq(start_deg, end_deg, length.out = n)
  tibble(
    x = c(r_inner * cos(deg2rad(ang)), rev(r_outer * cos(deg2rad(ang)))),
    y = c(r_inner * sin(deg2rad(ang)), rev(r_outer * sin(deg2rad(ang))))
  )
}

make_arc_line <- function(start_deg, end_deg, r, n = 240) {
  ang <- seq(start_deg, end_deg, length.out = n)
  tibble(
    x = r * cos(deg2rad(ang)),
    y = r * sin(deg2rad(ang))
  )
}

make_radial_line <- function(angle_deg, r_start, r_end) {
  tibble(
    x = c(r_start, r_end) * cos(deg2rad(angle_deg)),
    y = c(r_start, r_end) * sin(deg2rad(angle_deg))
  )
}

build_half_layout <- function(outcomes,
                              predictors,
                              start_deg,
                              end_deg,
                              sector_gap_deg = 2,
                              predictor_gap_deg = 0.3) {
  n_out <- length(outcomes)
  n_pred <- length(predictors)
  total_span <- end_deg - start_deg
  sector_span <- (total_span - sector_gap_deg * (n_out - 1)) / n_out
  slot_span <- (sector_span - predictor_gap_deg * (n_pred - 1)) / n_pred
  
  sector_tbl <- bind_rows(lapply(seq_along(outcomes), function(i) {
    sec_start <- start_deg + (i - 1) * (sector_span + sector_gap_deg)
    sec_end <- sec_start + sector_span
    tibble(
      outcome_short = outcomes[i],
      sector_id = i,
      sector_start = sec_start,
      sector_end = sec_end,
      sector_mid = (sec_start + sec_end) / 2
    )
  }))
  
  slot_tbl <- sector_tbl %>%
    slice(rep(1:n(), each = n_pred)) %>%
    mutate(
      predictor_short = rep(predictors, times = n_out),
      predictor_id = rep(seq_len(n_pred), times = n_out),
      slot_start = sector_start + (predictor_id - 1) * (slot_span + predictor_gap_deg),
      slot_end = slot_start + slot_span,
      angle_deg = (slot_start + slot_end) / 2
    )
  
  list(sectors = sector_tbl, slots = slot_tbl)
}

get_transform_fun <- function(transform_name, factor = 1) {
  factor <- factor %||% 1
  root_power <- if (isTRUE(is.finite(factor)) && factor > 0) factor else 2
  switch(
    transform_name,
    identity = function(x) x,
    log = function(x) x,
    log10 = function(x) sign(x) * log10(abs(x) + 1),
    sqrt = function(x) sign(x) * abs(x)^(1 / root_power),
    root10 = function(x) sign(x) * abs(x)^(1 / 10),
    exponential = function(x) sign(x) * (exp(factor * abs(x)) - 1),
    function(x) x
  )
}

# ==============================================================================
# Plot function
# ==============================================================================

plot_circular_prs_forest <- function(data,
                                     north_outcomes,
                                     south_outcomes,
                                     predictor_order,
                                     outcome_labels,
                                     predictor_labels,
                                     predictor_palette,
                                     bonferroni_denominator = 152,
                                     north_start = 10,
                                     north_end = 170,
                                     south_start = 190,
                                     south_end = 350,
                                     sector_gap_deg = 2.1,
                                     predictor_gap_deg = 0.6,
                                     density_scale = 12,
                                     or_limits = c(0.74, 2.08),
                                     or_breaks = c(0.9, 1.0, 1.1, 1.3, 1.6, 2),
                                     axis_transform = "log",
                                     axis_transform_factor = 2,
                                     r_inner = 2.2,
                                     r_outer = 6.8,
                                     r_label = 7.45,
                                     label_width = 14,
                                     label_angle_offset = 90,
                                     label_radius_offset = 0.45,
                                     label_angle_north_extra = 180,
                                     label_angle_south_extra = 0,
                                     point_size = 3.3,
                                     point_stroke = 0.35,
                                     ci_linewidth = 0.65,
                                     grid_linewidth = 0.35,
                                     tick_len = 0,
                                     tick_pad = 0,
                                     tick_angle_offset = 0,
                                     tick_text_angle = 0,
                                     tick_text_size = 4.5,
                                     tick_text_radial_extra = -0.4,
                                     tick_vjust_north = -0.6,
                                     tick_vjust_south = 1.6,
                                     outcome_text_size = 4.5,
                                     sector_fill_a = "grey92",
                                     sector_fill_b = "white",
                                     grey_sector_width_factor = 1.03,
                                     nonsig_col = "grey65",
                                     notch_half_height = 0.10,
                                     legend_title = "Process-specific PRS",
                                     legend_text_size = 11,
                                     legend_title_size = 13.5,
                                     legend_direction = "vertical",
                                     legend_point_size = 4,
                                     legend_x = -0.8,
                                     legend_y_top = 0.8,
                                     legend_spacing = 0.25,
                                     legend_title_x_offset = 0.7,
                                     background_fill = "white",
                                     piece_border_col = "grey60",
                                     piece_border_linewidth = 0.6,
                                     show_piece_border = TRUE,
                                     outer_ring_enabled = TRUE,
                                     outer_ring_threshold = 2,
                                     outer_ring_min = 6,
                                     outer_ring_max = 8,
                                     outer_ring_breaks = c(6, 8),
                                     outer_ring_outer = 7.7,
                                     outer_ring_gap = 0.22,
                                     outer_ring_border_col = "grey60",
                                     outer_ring_border_linewidth = 0.6,
                                     outer_ring_label_bump = 0.5,
                                     outer_tick_len = 0,
                                     outer_tick_text_size = 4.5,
                                     outer_tick_text_radial_extra = -0.4,
                                     north_axis_title = "OR for upper vs. bottom decile",
                                     south_axis_title = "OR for upper vs. bottom decile",
                                     north_axis_title_x = 5.7,
                                     north_axis_title_y = 0,
                                     south_axis_title_x = -5.7,
                                     south_axis_title_y = 0,
                                     north_axis_title_angle = 0,
                                     south_axis_title_angle = 0,
                                     axis_title_size = 5) {
  
  stopifnot(all(c("outcome_short", "predictor_short", "or", "lb", "ub", "p_value") %in% colnames(data)))
  stopifnot(length(bonferroni_denominator) == 1, is.finite(bonferroni_denominator), bonferroni_denominator > 0)
  
  transform_fun <- get_transform_fun(axis_transform, axis_transform_factor)
  
  tangential_rotation <- function(angle_deg, offset = 90) {
    angle_plot <- angle_deg + offset
    angle_plot <- ifelse(angle_plot > 180, angle_plot - 360, angle_plot)
    angle_plot <- ifelse(angle_plot < -180, angle_plot + 360, angle_plot)
    tibble(angle_plot = angle_plot)
  }
  
  log_or_to_r <- function(x, limits = or_limits, r0 = r_inner, r1 = r_outer) {
    lx <- log(x)
    ll <- log(limits)
    tx <- transform_fun(lx)
    tl <- transform_fun(ll)
    scales::rescale(tx, to = c(r0, r1), from = range(tl))
  }
  
  north_layout <- build_half_layout(
    outcomes = north_outcomes,
    predictors = predictor_order,
    start_deg = north_start,
    end_deg = north_end,
    sector_gap_deg = sector_gap_deg,
    predictor_gap_deg = predictor_gap_deg / density_scale
  )
  
  south_layout <- build_half_layout(
    outcomes = south_outcomes,
    predictors = predictor_order,
    start_deg = south_start,
    end_deg = south_end,
    sector_gap_deg = sector_gap_deg,
    predictor_gap_deg = predictor_gap_deg / density_scale
  )
  
  sector_layout <- bind_rows(
    north_layout$sectors %>% mutate(half = "north"),
    south_layout$sectors %>% mutate(half = "south")
  )
  
  slot_layout <- bind_rows(
    north_layout$slots %>% mutate(half = "north"),
    south_layout$slots %>% mutate(half = "south")
  )
  
  plot_df <- data %>%
    filter(
      outcome_short %in% c(north_outcomes, south_outcomes),
      predictor_short %in% predictor_order
    ) %>%
    inner_join(slot_layout, by = c("outcome_short", "predictor_short")) %>%
    mutate(
      outcome_label = dplyr::recode(outcome_short, !!!outcome_labels),
      outcome_label = stringr::str_wrap(outcome_label, width = label_width),
      predictor_label = dplyr::recode(predictor_short, !!!predictor_labels),
      p_bonferroni = pmin(p_value * bonferroni_denominator, 1),
      bonferroni_sig = !is.na(p_value) & p_value < (0.05 / bonferroni_denominator),
      nominal_sig = !is.na(p_value) & p_value < 0.05,
      r_est = log_or_to_r(or),
      r_lb = pmax(log_or_to_r(lb), r_inner),
      r_ub = pmin(log_or_to_r(ub), r_outer),
      point_col = unname(predictor_palette[predictor_short]),
      point_draw_col = ifelse(nominal_sig, point_col, nonsig_col),
      ci_draw_col = ifelse(bonferroni_sig, point_col, nonsig_col)
    )
  
  ci_df <- plot_df %>%
    mutate(
      is_extreme = outer_ring_enabled & (or >= outer_ring_threshold),
      x = r_est * cos(deg2rad(angle_deg)),
      y = r_est * sin(deg2rad(angle_deg)),
      x_lb = r_lb * cos(deg2rad(angle_deg)),
      y_lb = r_lb * sin(deg2rad(angle_deg)),
      x_ub = r_ub * cos(deg2rad(angle_deg)),
      y_ub = r_ub * sin(deg2rad(angle_deg))
    )
  
  outer_ring_inner <- r_outer + outer_ring_gap
  outer_limits <- c(outer_ring_min, outer_ring_max)
  if (!is.finite(outer_limits[1]) || !is.finite(outer_limits[2]) || outer_limits[2] <= outer_limits[1]) {
    outer_limits <- c(outer_ring_threshold, max(outer_ring_threshold * 1.1, outer_ring_threshold + 0.1))
  }
  outer_or_to_r <- function(x) scales::rescale(x, to = c(outer_ring_inner, outer_ring_outer), from = outer_limits)
  
  outer_ring_breaks <- outer_ring_breaks %||% pretty(outer_limits, n = 3)
  outer_ring_breaks <- outer_ring_breaks[outer_ring_breaks >= outer_limits[1] & outer_ring_breaks <= outer_limits[2]]
  if (!length(outer_ring_breaks)) {
    outer_ring_breaks <- outer_limits
  }
  outer_ring_breaks <- sort(unique(outer_ring_breaks))
  
  ci_outer_df <- ci_df %>%
    filter(is_extreme) %>%
    mutate(
      r_est_outer = outer_or_to_r(or),
      r_lb_outer = outer_or_to_r(pmin(pmax(lb, outer_ring_min), outer_ring_max)),
      r_ub_outer = outer_or_to_r(pmin(pmax(ub, outer_ring_min), outer_ring_max)),
      x_outer = r_est_outer * cos(deg2rad(angle_deg)),
      y_outer = r_est_outer * sin(deg2rad(angle_deg)),
      x_lb_outer = r_lb_outer * cos(deg2rad(angle_deg)),
      y_lb_outer = r_lb_outer * sin(deg2rad(angle_deg)),
      x_ub_outer = r_ub_outer * cos(deg2rad(angle_deg)),
      y_ub_outer = r_ub_outer * sin(deg2rad(angle_deg))
    )
  
  ci_inner_df <- ci_df %>%
    filter(!is_extreme)
  
  sector_poly <- sector_layout %>%
    arrange(half, sector_id) %>%
    mutate(
      fill_group = rep(c(sector_fill_a, sector_fill_b), length.out = n()),
      fill_width_factor = ifelse(fill_group == sector_fill_a, grey_sector_width_factor, 1),
      fill_half_width = (sector_end - sector_start) / 2 * fill_width_factor,
      half_start = ifelse(half == "north", north_start, south_start),
      half_end = ifelse(half == "north", north_end, south_end),
      fill_start = pmax(sector_mid - fill_half_width, half_start),
      fill_end = pmin(sector_mid + fill_half_width, half_end)
    ) %>%
    rowwise() %>%
    mutate(poly = list(make_arc_polygon(fill_start, fill_end, r_inner, r_outer))) %>%
    ungroup() %>%
    unnest(poly)
  
  grid_df <- bind_rows(
    bind_rows(lapply(or_breaks, function(v) {
      make_arc_line(north_start, north_end, log_or_to_r(v)) %>% mutate(or_break = v, half = "north")
    })),
    bind_rows(lapply(or_breaks, function(v) {
      make_arc_line(south_start, south_end, log_or_to_r(v)) %>% mutate(or_break = v, half = "south")
    }))
  )
  
  label_df <- sector_layout %>%
    mutate(
      outcome_label = dplyr::recode(outcome_short, !!!outcome_labels),
      outcome_label = stringr::str_wrap(outcome_label, width = label_width),
      r_lab = r_label + label_radius_offset + ifelse(outer_ring_enabled, outer_ring_label_bump, 0),
      angle_add = ifelse(half == "north", label_angle_north_extra, label_angle_south_extra)
    ) %>%
    bind_cols(tangential_rotation(.$sector_mid, offset = label_angle_offset + .$angle_add)) %>%
    mutate(
      x = r_lab * cos(deg2rad(sector_mid)),
      y = r_lab * sin(deg2rad(sector_mid))
    )
  
  tick_angles <- c(north_start, north_end, south_start, south_end) + tick_angle_offset
  tick_df <- expand.grid(or_break = or_breaks, end_angle = tick_angles) %>%
    as_tibble() %>%
    mutate(
      r = log_or_to_r(or_break),
      x = (r + 0.01) * cos(deg2rad(end_angle)),
      y = (r + 0.01) * sin(deg2rad(end_angle)),
      xend = (r + 0.01 + tick_len) * cos(deg2rad(end_angle)),
      yend = (r + 0.01 + tick_len) * sin(deg2rad(end_angle))
    )
  
  tick_lab_df <- tick_df %>%
    mutate(
      lx = (r + tick_pad) * cos(deg2rad(end_angle)),
      ly = (r + tick_pad) * sin(deg2rad(end_angle)),
      lab = as.character(or_break),
      angle = tick_text_angle,
      vjust = ifelse(sin(deg2rad(end_angle)) >= 0, tick_vjust_north, tick_vjust_south),
      hjust = 0.5,
      ly = ly + ifelse(sin(deg2rad(end_angle)) >= 0, tick_text_radial_extra, -tick_text_radial_extra)
    )
  
  outer_grid_df <- if (outer_ring_enabled) {
    bind_rows(
      bind_rows(lapply(outer_ring_breaks, function(v) {
        make_arc_line(north_start, north_end, outer_or_to_r(v)) %>% mutate(or_break = v, half = "north")
      })),
      bind_rows(lapply(outer_ring_breaks, function(v) {
        make_arc_line(south_start, south_end, outer_or_to_r(v)) %>% mutate(or_break = v, half = "south")
      }))
    )
  } else {
    tibble()
  }
  
  outer_tick_df <- if (outer_ring_enabled) {
    expand.grid(or_break = outer_ring_breaks, end_angle = tick_angles) %>%
      as_tibble() %>%
      mutate(
        r = outer_or_to_r(or_break),
        x = (r + 0.01) * cos(deg2rad(end_angle)),
        y = (r + 0.01) * sin(deg2rad(end_angle)),
        xend = (r + 0.01 + outer_tick_len) * cos(deg2rad(end_angle)),
        yend = (r + 0.01 + outer_tick_len) * sin(deg2rad(end_angle))
      )
  } else {
    tibble()
  }
  
  outer_tick_lab_df <- if (outer_ring_enabled) {
    outer_tick_df %>%
      mutate(
        lx = (r + outer_tick_len + tick_pad) * cos(deg2rad(end_angle)),
        ly = (r + outer_tick_len + tick_pad) * sin(deg2rad(end_angle)),
        lab = as.character(or_break),
        angle = tick_text_angle,
        vjust = ifelse(sin(deg2rad(end_angle)) >= 0, tick_vjust_north, tick_vjust_south),
        hjust = 0.5,
        ly = ly + ifelse(sin(deg2rad(end_angle)) >= 0, outer_tick_text_radial_extra, -outer_tick_text_radial_extra)
      )
  } else {
    tibble()
  }
  
  axis_title_df <- tibble(
    x = c(north_axis_title_x, south_axis_title_x),
    y = c(north_axis_title_y, south_axis_title_y),
    label = c(north_axis_title, south_axis_title),
    angle = c(north_axis_title_angle, south_axis_title_angle)
  )
  
  notch_df <- tibble(
    xmin = -outer_ring_outer - 2,
    xmax = outer_ring_outer + 2,
    ymin = -notch_half_height,
    ymax = notch_half_height
  )
  
  piece_border_df <- bind_rows(
    make_arc_line(north_start, north_end, r_inner) %>% mutate(border_id = "north_inner"),
    make_arc_line(north_start, north_end, r_outer) %>% mutate(border_id = "north_outer"),
    make_arc_line(south_start, south_end, r_inner) %>% mutate(border_id = "south_inner"),
    make_arc_line(south_start, south_end, r_outer) %>% mutate(border_id = "south_outer"),
    make_radial_line(north_start, r_inner, r_outer) %>% mutate(border_id = "north_start"),
    make_radial_line(north_end, r_inner, r_outer) %>% mutate(border_id = "north_end"),
    make_radial_line(south_start, r_inner, r_outer) %>% mutate(border_id = "south_start"),
    make_radial_line(south_end, r_inner, r_outer) %>% mutate(border_id = "south_end")
  )
  
  outer_ring_sector_poly <- if (outer_ring_enabled) {
    sector_layout %>%
      arrange(half, sector_id) %>%
      mutate(
        fill_group = rep(c(sector_fill_a, sector_fill_b), length.out = n()),
        fill_width_factor = ifelse(fill_group == sector_fill_a, grey_sector_width_factor, 1),
        fill_half_width = (sector_end - sector_start) / 2 * fill_width_factor,
      half_start = ifelse(half == "north", north_start, south_start),
      half_end = ifelse(half == "north", north_end, south_end),
      fill_start = pmax(sector_mid - fill_half_width, half_start),
      fill_end = pmin(sector_mid + fill_half_width, half_end)
      ) %>%
      rowwise() %>%
      mutate(poly = list(make_arc_polygon(fill_start, fill_end, outer_ring_inner, outer_ring_outer))) %>%
      ungroup() %>%
      unnest(poly)
  } else {
    tibble(x = numeric(), y = numeric(), outcome_short = character(), half = character(), fill_group = character())
  }
  
  outer_ring_arc_df <- if (outer_ring_enabled) bind_rows(
    make_arc_line(north_start, north_end, outer_ring_inner) %>% mutate(border_id = "north_oring_inner"),
    make_arc_line(north_start, north_end, outer_ring_outer) %>% mutate(border_id = "north_oring_outer"),
    make_arc_line(south_start, south_end, outer_ring_inner) %>% mutate(border_id = "south_oring_inner"),
    make_arc_line(south_start, south_end, outer_ring_outer) %>% mutate(border_id = "south_oring_outer")
  ) else tibble(x = numeric(), y = numeric(), border_id = character())
  
  outer_ring_side_df <- if (outer_ring_enabled) bind_rows(
    make_radial_line(north_start, outer_ring_inner, outer_ring_outer) %>% mutate(border_id = "north_oring_start"),
    make_radial_line(north_end, outer_ring_inner, outer_ring_outer) %>% mutate(border_id = "north_oring_end"),
    make_radial_line(south_start, outer_ring_inner, outer_ring_outer) %>% mutate(border_id = "south_oring_start"),
    make_radial_line(south_end, outer_ring_inner, outer_ring_outer) %>% mutate(border_id = "south_oring_end")
  ) else tibble(x = numeric(), y = numeric(), border_id = character())
  
  legend_df <- tibble(
    predictor_short = predictor_order,
    predictor_label = dplyr::recode(predictor_order, !!!predictor_labels),
    col = predictor_palette[predictor_order],
    x = 0,
    y = seq(0.9, by = -0.25, length.out = length(predictor_order))
  )
  
  if (legend_direction == "vertical") {
    legend_df <- legend_df %>%
      mutate(x = legend_x, y = legend_y_top - (row_number() - 1) * legend_spacing)
  } else {
    legend_df <- legend_df %>%
      mutate(x = legend_x + (row_number() - 1) * legend_spacing, y = legend_y_top)
  }
  
  p <- ggplot() +
    geom_polygon(
      data = sector_poly,
      aes(x = x, y = y, group = interaction(outcome_short, half), fill = fill_group),
      color = NA
    ) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    geom_path(
      data = grid_df,
      aes(x = x, y = y, group = interaction(or_break, half)),
      color = "grey75",
      linewidth = grid_linewidth,
      lineend = "butt"
    ) +
    geom_rect(
      data = notch_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = background_fill,
      color = NA
    )
  
  if (show_piece_border) {
    p <- p + geom_path(
      data = piece_border_df,
      aes(x = x, y = y, group = border_id),
      inherit.aes = FALSE,
      color = piece_border_col,
      linewidth = piece_border_linewidth,
      lineend = "butt"
    )
  }
  
  if (outer_ring_enabled && nrow(outer_ring_sector_poly) > 0) {
    p <- p + geom_polygon(
      data = outer_ring_sector_poly,
      aes(x = x, y = y, group = interaction(outcome_short, half), fill = fill_group),
      inherit.aes = FALSE,
      color = NA
    )
  }
  
  if (outer_ring_enabled && nrow(outer_grid_df) > 0) {
    p <- p + geom_path(
      data = outer_grid_df,
      aes(x = x, y = y, group = interaction(or_break, half)),
      color = "grey75",
      linewidth = grid_linewidth,
      lineend = "butt"
    )
  }
  
  if (outer_ring_enabled && nrow(outer_ring_arc_df) > 0) {
    p <- p + geom_path(
      data = outer_ring_arc_df,
      aes(x = x, y = y, group = border_id),
      inherit.aes = FALSE,
      color = outer_ring_border_col,
      linewidth = outer_ring_border_linewidth,
      lineend = "butt"
    )
  }
  
  if (outer_ring_enabled && nrow(outer_ring_side_df) > 0) {
    p <- p + geom_path(
      data = outer_ring_side_df,
      aes(x = x, y = y, group = border_id),
      inherit.aes = FALSE,
      color = outer_ring_border_col,
      linewidth = outer_ring_border_linewidth,
      lineend = "butt"
    )
  }
  
  p <- p +
    geom_segment(
      data = ci_inner_df,
      aes(x = x_lb, y = y_lb, xend = x_ub, yend = y_ub, color = ci_draw_col),
      linewidth = ci_linewidth,
      lineend = "round",
      show.legend = FALSE
    ) +
    scale_color_identity() +
    geom_point(
      data = ci_inner_df,
      aes(x = x, y = y, fill = point_draw_col, color = point_draw_col),
      shape = 21,
      size = point_size,
      stroke = point_stroke,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_color_identity()
  
  if (outer_ring_enabled && nrow(ci_outer_df) > 0) {
    p <- p +
      geom_segment(
        data = ci_outer_df,
        aes(x = x_lb_outer, y = y_lb_outer, xend = x_ub_outer, yend = y_ub_outer, color = ci_draw_col),
        linewidth = ci_linewidth,
        lineend = "round",
        show.legend = FALSE
      ) +
      geom_point(
        data = ci_outer_df,
        aes(x = x_outer, y = y_outer, fill = point_draw_col, color = point_draw_col),
        shape = 21,
        size = point_size,
        stroke = point_stroke,
        show.legend = FALSE
      )
  }
  
  p <- p +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = outcome_label, angle = angle_plot),
      size = outcome_text_size,
      lineheight = 0.9,
      hjust = 0.5,
      vjust = 0.5
    ) +
    geom_segment(
      data = tick_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = grid_linewidth,
      color = "black"
    ) +
    geom_text(
      data = tick_lab_df,
      aes(x = lx, y = ly, label = lab, hjust = hjust, vjust = vjust, angle = angle),
      size = tick_text_size
    ) +
    coord_equal(clip = "off") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(10, 20, 10, 20),
      plot.background = element_rect(fill = background_fill, color = NA),
      panel.background = element_rect(fill = background_fill, color = NA)
    )
  
  if (outer_ring_enabled && nrow(outer_tick_df) > 0) {
    p <- p +
      geom_segment(
        data = outer_tick_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        linewidth = grid_linewidth,
        color = "black"
      ) +
      geom_text(
        data = outer_tick_lab_df,
        aes(x = lx, y = ly, label = lab, hjust = hjust, vjust = vjust, angle = angle),
        size = outer_tick_text_size
      ) +
      geom_text(
        data = axis_title_df,
        aes(x = x, y = y, label = label, angle = angle),
        inherit.aes = FALSE,
        size = axis_title_size,
        hjust = 0.5,
        vjust = 0.5
      )
  }
  
  if (legend_direction == "vertical") {
    text_x <- legend_df$x + 0.18
    text_y <- legend_df$y
    text_hjust <- 0
    text_vjust <- 0.5
    title_x <- legend_x + legend_title_x_offset
    title_y <- max(legend_df$y) + 0.35
  } else {
    text_x <- legend_df$x
    text_y <- legend_df$y - 0.2
    text_hjust <- 0.5
    text_vjust <- -1.1
    title_x <- legend_title_x_offset
    title_y <- 0.35
  }
  
  p +
    geom_point(
      data = legend_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = legend_point_size,
      stroke = point_stroke,
      fill = legend_df$col,
      color = legend_df$col
    ) +
    geom_text(
      data = transform(legend_df, text_x = text_x, text_y = text_y),
      aes(x = text_x, y = text_y, label = predictor_label),
      inherit.aes = FALSE,
      hjust = text_hjust,
      vjust = text_vjust,
      size = legend_text_size / 3.2
    ) +
    annotate(
      "text",
      x = title_x,
      y = title_y,
      label = legend_title,
      size = legend_title_size / 3.2,
      fontface = "bold"
    )
}

# ==============================================================================
# UI helpers
# ==============================================================================

make_prs_color_controls <- function(ids, cols) {
  tagList(
    h4("Predictor colors"),
    make_predictor_color_inputs(ids, cols)
  )
}

# ==============================================================================
# UI
# ==============================================================================

example <- generate_example_data()

ui <- fluidPage(
  titlePanel("Circular PRS Forest Plot Shiny App"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tabsetPanel(
        id = "control_tabs",
        tabPanel(
          "Data",
          checkboxInput("use_example", "Use example data", TRUE),
          fileInput("file", "Upload CSV"),
          tags$hr(),
          textAreaInput("north_outcomes", "North outcomes (one per line)",
                        value = paste(example$north, collapse = "\n"), rows = 8),
          textAreaInput("south_outcomes", "South outcomes (one per line)",
                        value = paste(example$south, collapse = "\n"), rows = 8),
          textAreaInput("predictor_order", "Predictor order (one per line)",
                        value = paste(example$predictors, collapse = "\n"), rows = 8),
          tags$hr(),
          textAreaInput("outcome_map", "Outcome labels (one per line: code = label)", rows = 10,
                        value = paste0(names(example$outcome_labels), " = ", unname(example$outcome_labels), collapse = "\n")),
          textAreaInput("predictor_map", "Predictor labels (one per line: code = label)", rows = 8,
                        value = paste0(names(example$predictor_labels), " = ", unname(example$predictor_labels), collapse = "\n"))
        ),
        tabPanel(
          "Geometry",
          fluidRow(
            column(6, numericInput("north_start", "North start", 10, step = 1)),
            column(6, numericInput("north_end", "North end", 170, step = 1))
          ),
          fluidRow(
            column(6, numericInput("south_start", "South start", 190, step = 1)),
            column(6, numericInput("south_end", "South end", 350, step = 1))
          ),
          fluidRow(
            column(4, numericInput("sector_gap_deg", "Sector gap (deg)", 2.1, step = 0.1)),
            column(4, numericInput("predictor_gap_deg", "Predictor gap (deg)", 0.6, step = 0.05)),
            column(4, numericInput("density_scale", "Density scale", 12, step = 0.1))
          ),
          fluidRow(
            column(4, numericInput("r_inner", "Inner radius", 2.2, step = 0.1)),
            column(4, numericInput("r_outer", "Outer radius", 6.8, step = 0.1)),
            column(4, numericInput("r_label", "Label radius", 7.45, step = 0.1))
          ),
          numericInput("notch_half_height", "Notch half-height", 0.10, step = 0.01),
          checkboxInput("outer_ring_enabled", "Enable outer ring for extreme estimates", TRUE),
          numericInput("outer_ring_threshold", "Outer ring threshold (OR)", 2, step = 0.1),
          numericInput("outer_ring_min", "Outer ring axis min", 6, step = 0.1),
          numericInput("outer_ring_max", "Outer ring axis max", 8, step = 0.1),
          textInput("outer_ring_breaks", "Outer ring breaks (comma-separated, blank = auto)", "6, 8"),
          numericInput("outer_ring_outer", "Outer ring outer radius", 7.7, step = 0.1),
          numericInput("outer_ring_gap", "Outer ring gap", 0.22, step = 0.01),
          numericInput("outer_ring_label_bump", "Outer label bump", 0.5, step = 0.05)
        ),
        tabPanel(
          "Axes",
          selectInput("axis_transform", "Axis transformation",
                      choices = c("log", "log10", "sqrt", "root10", "identity", "exponential"),
                      selected = "log"),
          numericInput("axis_transform_factor", "Transform factor / root basis", 2, step = 0.1),
          helpText("For exponential: factor controls curvature. For sqrt: factor = 2 gives square root, 3 gives cube-root-like compression, 4 gives fourth-root-like compression."),
          fluidRow(
            column(6, numericInput("axis_min", "Axis min", 0.74, step = 0.01)),
            column(6, numericInput("axis_max", "Axis max", 2.08, step = 0.01))
          ),
          textInput("axis_breaks", "Axis breaks (comma-separated)", "0.9, 1.0, 1.1, 1.3, 1.6, 2"),
          fluidRow(
            column(4, numericInput("tick_len", "Tick length", 0, step = 0.01)),
            column(4, numericInput("tick_pad", "Tick pad", 0, step = 0.01)),
            column(4, numericInput("tick_angle_offset", "Tick angle offset", 0, step = 0.5))
          ),
          fluidRow(
            column(4, numericInput("tick_text_size", "Tick text size", 4.5, step = 0.1)),
            column(4, numericInput("tick_text_angle", "Tick text angle", 0, step = 1)),
            column(4, numericInput("tick_text_radial_extra", "Tick text radial offset", -0.4, step = 0.01))
          ),
          fluidRow(
            column(4, numericInput("outer_tick_len", "Outer tick length", 0, step = 0.01)),
            column(4, numericInput("outer_tick_text_size", "Outer tick text size", 4.5, step = 0.1)),
            column(4, numericInput("outer_tick_text_radial_extra", "Outer text radial offset", -0.4, step = 0.01))
          ),
          fluidRow(
            column(4, numericInput("axis_title_size", "Axis title size", 5, step = 0.1)),
            column(4, textInput("north_axis_title", "North axis title", "OR for upper vs. bottom decile")),
            column(4, textInput("south_axis_title", "South axis title", "OR for upper vs. bottom decile"))
          ),
          fluidRow(
            column(3, numericInput("north_axis_title_x", "North title x", 5.7, step = 0.05)),
            column(3, numericInput("north_axis_title_y", "North title y", 0, step = 0.05)),
            column(3, numericInput("north_axis_title_angle", "North title angle", 0, step = 1)),
            column(3, tags$div())
          ),
          fluidRow(
            column(3, numericInput("south_axis_title_x", "South title x", -5.7, step = 0.05)),
            column(3, numericInput("south_axis_title_y", "South title y", 0, step = 0.05)),
            column(3, numericInput("south_axis_title_angle", "South title angle", 0, step = 1)),
            column(3, tags$div())
          ),
          fluidRow(
            column(6, numericInput("tick_vjust_north", "Tick vjust north", -0.6, step = 0.1)),
            column(6, numericInput("tick_vjust_south", "Tick vjust south", 1.6, step = 0.1))
          )
        ),
        tabPanel(
          "Labels",
          fluidRow(
            column(4, numericInput("label_width", "Wrap width", 14, step = 1)),
            column(4, numericInput("label_angle_offset", "Base label angle", 90, step = 1)),
            column(4, numericInput("label_radius_offset", "Label radius offset", 0.45, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("label_angle_north_extra", "North label angle extra", 180, step = 1)),
            column(6, numericInput("label_angle_south_extra", "South label angle extra", 0, step = 1))
          ),
          numericInput("outcome_text_size", "Outcome text size", 4.5, step = 0.1)
        ),
        tabPanel(
          "Markers & Lines",
          numericInput("bonferroni_denominator", "Bonferroni denominator", 152, min = 1, step = 1),
          helpText("Bonferroni significance threshold = 0.05 / denominator."),
          fluidRow(
            column(4, numericInput("point_size", "Point size", 3.3, step = 0.1)),
            column(4, numericInput("point_stroke", "Point stroke", 0.35, step = 0.05)),
            column(4, numericInput("ci_linewidth", "CI linewidth", 0.65, step = 0.05))
          ),
          fluidRow(
            column(6, numericInput("grid_linewidth", "Grid linewidth", 0.35, step = 0.05)),
            column(6, textInput("nonsig_col", "Non-significant color", "grey65"))
          )
        ),
        tabPanel(
          "Colors & Legend",
          uiOutput("predictor_colors_ui"),
          tags$hr(),
          fluidRow(
            column(4, textInput("sector_fill_a", "Sector fill A", "grey92")),
            column(4, textInput("sector_fill_b", "Sector fill B", "white")),
            column(4, numericInput("grey_sector_width_factor", "Grey-sector width factor", 1.03, min = 0.90, max = 1.10, step = 0.005))
          ),
          textInput("background_fill", "Background fill", "white"),
          checkboxInput("show_piece_border", "Show north/south piece border", TRUE),
          fluidRow(
            column(6, textInput("piece_border_col", "Piece border color", "grey60")),
            column(6, numericInput("piece_border_linewidth", "Piece border linewidth", 0.6, step = 0.05))
          ),
          fluidRow(
            column(6, textInput("outer_ring_border_col", "Outer ring border color", "grey60")),
            column(6, numericInput("outer_ring_border_linewidth", "Outer ring border linewidth", 0.6, step = 0.05))
          ),
          tags$hr(),
          textInput("legend_title", "Legend title", "Process-specific PRS"),
          selectInput("legend_direction", "Legend direction", choices = c("vertical", "horizontal"), selected = "vertical"),
          fluidRow(
            column(3, numericInput("legend_x", "Legend x", -0.8, step = 0.05)),
            column(3, numericInput("legend_y_top", "Legend top y", 0.8, step = 0.05)),
            column(3, numericInput("legend_spacing", "Legend spacing", 0.25, step = 0.01)),
            column(3, numericInput("legend_title_x_offset", "Legend title x offset", 0.7, step = 0.05))
          ),
          fluidRow(
            column(4, numericInput("legend_text_size", "Legend text size", 11, step = 0.5)),
            column(4, numericInput("legend_title_size", "Legend title size", 13.5, step = 0.5)),
            column(4, numericInput("legend_point_size", "Legend point size", 4, step = 0.1))
          )
        ),
        tabPanel(
          "Export",
          numericInput("png_width", "PNG width (in)", 10, step = 0.5),
          numericInput("png_height", "PNG height (in)", 10, step = 0.5),
          numericInput("png_dpi", "PNG DPI", 300, min = 72, step = 1),
          textInput("png_bg", "PNG background", "white"),
          downloadButton("download_png", "Download PNG"),
          tags$hr(),
          h4("Plot settings (JSON)"),
          fileInput("settings_json", "Import plot settings", accept = c(".json", "application/json")),
          actionButton("import_settings", "Apply imported settings"),
          downloadButton("download_settings", "Export plot settings")
        )
      )
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot", height = "950px")),
        tabPanel("Data preview", DTOutput("data_table"))
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {
  # Inputs included in JSON plot-setting import/export. Data files and PNG export
  # dimensions are intentionally excluded.
  numeric_setting_ids <- c(
    "north_start", "north_end", "south_start", "south_end",
    "sector_gap_deg", "predictor_gap_deg", "density_scale",
    "r_inner", "r_outer", "r_label", "notch_half_height",
    "outer_ring_threshold", "outer_ring_min", "outer_ring_max",
    "outer_ring_outer", "outer_ring_gap", "outer_ring_label_bump",
    "axis_transform_factor", "axis_min", "axis_max", "tick_len", "tick_pad",
    "tick_angle_offset", "tick_text_size", "tick_text_angle",
    "tick_text_radial_extra", "outer_tick_len", "outer_tick_text_size",
    "outer_tick_text_radial_extra", "axis_title_size", "north_axis_title_x",
    "north_axis_title_y", "north_axis_title_angle", "south_axis_title_x",
    "south_axis_title_y", "south_axis_title_angle", "tick_vjust_north",
    "tick_vjust_south", "label_width", "label_angle_offset",
    "label_radius_offset", "label_angle_north_extra", "label_angle_south_extra",
    "outcome_text_size", "point_size", "point_stroke", "ci_linewidth",
    "bonferroni_denominator", "grid_linewidth", "grey_sector_width_factor", "piece_border_linewidth", "outer_ring_border_linewidth",
    "legend_x", "legend_y_top", "legend_spacing", "legend_title_x_offset",
    "legend_text_size", "legend_title_size", "legend_point_size"
  )
  textarea_setting_ids <- c(
    "north_outcomes", "south_outcomes", "predictor_order", "outcome_map", "predictor_map"
  )
  text_setting_ids <- c(
    "outer_ring_breaks", "axis_breaks", "north_axis_title", "south_axis_title",
    "nonsig_col", "sector_fill_a", "sector_fill_b", "background_fill",
    "piece_border_col", "outer_ring_border_col", "legend_title"
  )
  checkbox_setting_ids <- c("outer_ring_enabled", "show_piece_border")
  select_setting_ids <- c("axis_transform", "legend_direction")
  color_override <- reactiveVal(NULL)
  current_data <- reactive({
    if (isTRUE(input$use_example)) {
      return(generate_example_data()$data)
    }
    shiny::validate(shiny::need(!is.null(input$file), "Please upload a CSV file or use example data."))
    safe_read_csv(input$file$datapath)
  })
  
  north_outcomes <- reactive(parse_lines(unlist(strsplit(input$north_outcomes, "\n"))))
  south_outcomes <- reactive(parse_lines(unlist(strsplit(input$south_outcomes, "\n"))))
  predictor_order <- reactive(parse_lines(unlist(strsplit(input$predictor_order, "\n"))))
  parsed_outcome_labels <- reactive(parse_named_mapping(input$outcome_map))
  parsed_predictor_labels <- reactive(parse_named_mapping(input$predictor_map))
  
  output$predictor_colors_ui <- renderUI({
    ids <- predictor_order()
    fallback <- setNames(make_default_palette(length(ids)), ids)
    imported <- color_override()
    cols <- unname(example$predictor_colors[ids])
    cols[is.na(cols)] <- fallback[ids][is.na(cols)]
    if (!is.null(imported)) {
      hit <- ids %in% names(imported)
      cols[hit] <- unname(imported[ids[hit]])
    }
    make_prs_color_controls(ids, cols)
  })
  
  predictor_palette <- reactive({
    ids <- predictor_order()
    fallback <- setNames(make_default_palette(length(ids)), ids)
    defaults <- unname(example$predictor_colors[ids])
    defaults[is.na(defaults)] <- fallback[ids][is.na(defaults)]
    imported <- color_override()
    if (!is.null(imported)) {
      hit <- ids %in% names(imported)
      defaults[hit] <- unname(imported[ids[hit]])
    }
    setNames(
      vapply(seq_along(ids), function(i) input[[paste0("pred_col_", i)]] %||% defaults[i], character(1)),
      ids
    )
  })
  
  axis_breaks <- reactive({
    as.numeric(trimws(unlist(strsplit(input$axis_breaks, ","))))
  })
  
  outer_axis_breaks <- reactive({
    vals <- trimws(unlist(strsplit(input$outer_ring_breaks %||% "", ",")))
    vals <- vals[nzchar(vals)]
    if (!length(vals)) return(NULL)
    as.numeric(vals)
  })
  
  plot_settings <- reactive({
    values <- c(
      setNames(lapply(numeric_setting_ids, function(id) input[[id]]), numeric_setting_ids),
      setNames(lapply(textarea_setting_ids, function(id) input[[id]]), textarea_setting_ids),
      setNames(lapply(text_setting_ids, function(id) input[[id]]), text_setting_ids),
      setNames(lapply(checkbox_setting_ids, function(id) input[[id]]), checkbox_setting_ids),
      setNames(lapply(select_setting_ids, function(id) input[[id]]), select_setting_ids)
    )
    list(
      format = "figure5-circular-prs-plot-settings",
      version = 1,
      settings = values,
      predictor_colors = as.list(predictor_palette())
    )
  })

  output$download_settings <- downloadHandler(
    filename = function() paste0("figure5_plot_settings_", Sys.Date(), ".json"),
    content = function(file) {
      jsonlite::write_json(plot_settings(), file, pretty = TRUE, auto_unbox = TRUE, null = "null")
    }
  )

  observeEvent(input$import_settings, {
    req(input$settings_json)
    cfg <- tryCatch(
      jsonlite::read_json(input$settings_json$datapath, simplifyVector = TRUE),
      error = function(e) {
        showNotification(paste("Could not read JSON:", e$message), type = "error")
        NULL
      }
    )
    req(cfg)
    vals <- cfg$settings %||% cfg
    for (id in intersect(numeric_setting_ids, names(vals))) {
      updateNumericInput(session, id, value = as.numeric(vals[[id]]))
    }
    for (id in intersect(textarea_setting_ids, names(vals))) {
      updateTextAreaInput(session, id, value = as.character(vals[[id]]))
    }
    for (id in intersect(text_setting_ids, names(vals))) {
      updateTextInput(session, id, value = as.character(vals[[id]]))
    }
    for (id in intersect(checkbox_setting_ids, names(vals))) {
      updateCheckboxInput(session, id, value = isTRUE(vals[[id]]))
    }
    for (id in intersect(select_setting_ids, names(vals))) {
      updateSelectInput(session, id, selected = as.character(vals[[id]]))
    }
    if (!is.null(cfg$predictor_colors)) {
      imported_cols <- unlist(cfg$predictor_colors, use.names = TRUE)
      color_override(imported_cols)
      session$onFlushed(function() {
        ids <- predictor_order()
        for (i in seq_along(ids)) {
          if (ids[i] %in% names(imported_cols)) {
            updateTextInput(session, paste0("pred_col_", i), value = imported_cols[[ids[i]]])
          }
        }
      }, once = TRUE)
    }
    showNotification("Plot settings imported.", type = "message")
  })

  current_plot <- reactive({
    plot_circular_prs_forest(
      data = current_data(),
      north_outcomes = north_outcomes(),
      south_outcomes = south_outcomes(),
      predictor_order = predictor_order(),
      outcome_labels = parsed_outcome_labels(),
      predictor_labels = parsed_predictor_labels(),
      predictor_palette = predictor_palette(),
      bonferroni_denominator = input$bonferroni_denominator,
      north_start = input$north_start,
      north_end = input$north_end,
      south_start = input$south_start,
      south_end = input$south_end,
      sector_gap_deg = input$sector_gap_deg,
      predictor_gap_deg = input$predictor_gap_deg,
      density_scale = input$density_scale,
      or_limits = c(input$axis_min, input$axis_max),
      or_breaks = axis_breaks(),
      axis_transform = input$axis_transform,
      axis_transform_factor = input$axis_transform_factor,
      r_inner = input$r_inner,
      r_outer = input$r_outer,
      r_label = input$r_label,
      label_width = input$label_width,
      label_angle_offset = input$label_angle_offset,
      label_radius_offset = input$label_radius_offset,
      label_angle_north_extra = input$label_angle_north_extra,
      label_angle_south_extra = input$label_angle_south_extra,
      point_size = input$point_size,
      point_stroke = input$point_stroke,
      ci_linewidth = input$ci_linewidth,
      grid_linewidth = input$grid_linewidth,
      tick_len = input$tick_len,
      tick_pad = input$tick_pad,
      tick_angle_offset = input$tick_angle_offset,
      tick_text_angle = input$tick_text_angle,
      tick_text_size = input$tick_text_size,
      tick_text_radial_extra = input$tick_text_radial_extra,
      tick_vjust_north = input$tick_vjust_north,
      tick_vjust_south = input$tick_vjust_south,
      outcome_text_size = input$outcome_text_size,
      sector_fill_a = input$sector_fill_a,
      sector_fill_b = input$sector_fill_b,
      grey_sector_width_factor = input$grey_sector_width_factor,
      nonsig_col = input$nonsig_col,
      notch_half_height = input$notch_half_height,
      legend_title = input$legend_title,
      legend_text_size = input$legend_text_size,
      legend_title_size = input$legend_title_size,
      legend_direction = input$legend_direction,
      legend_point_size = input$legend_point_size,
      legend_x = input$legend_x,
      legend_y_top = input$legend_y_top,
      legend_spacing = input$legend_spacing,
      legend_title_x_offset = input$legend_title_x_offset,
      background_fill = input$background_fill,
      piece_border_col = input$piece_border_col,
      piece_border_linewidth = input$piece_border_linewidth,
      show_piece_border = input$show_piece_border,
      outer_ring_enabled = input$outer_ring_enabled,
      outer_ring_threshold = input$outer_ring_threshold,
      outer_ring_min = input$outer_ring_min,
      outer_ring_max = input$outer_ring_max,
      outer_ring_breaks = outer_axis_breaks(),
      outer_ring_outer = input$outer_ring_outer,
      outer_ring_gap = input$outer_ring_gap,
      outer_ring_border_col = input$outer_ring_border_col,
      outer_ring_border_linewidth = input$outer_ring_border_linewidth,
      outer_ring_label_bump = input$outer_ring_label_bump,
      outer_tick_len = input$outer_tick_len,
      outer_tick_text_size = input$outer_tick_text_size,
      outer_tick_text_radial_extra = input$outer_tick_text_radial_extra,
      north_axis_title = input$north_axis_title,
      south_axis_title = input$south_axis_title,
      north_axis_title_x = input$north_axis_title_x,
      north_axis_title_y = input$north_axis_title_y,
      south_axis_title_x = input$south_axis_title_x,
      south_axis_title_y = input$south_axis_title_y,
      north_axis_title_angle = input$north_axis_title_angle,
      south_axis_title_angle = input$south_axis_title_angle,
      axis_title_size = input$axis_title_size
    )
  })
  
  output$plot <- renderPlot({
    print(current_plot())
  }, bg = "white")
  
  output$data_table <- renderDT({
    datatable(head(current_data(), 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$download_png <- downloadHandler(
    filename = function() paste0("circular_prs_forest_", Sys.Date(), ".png"),
    content = function(file) {
      ggsave(
        filename = file,
        plot = current_plot(),
        width = input$png_width,
        height = input$png_height,
        dpi = input$png_dpi,
        bg = input$png_bg
      )
    }
  )
}

shinyApp(ui = ui, server = server)
