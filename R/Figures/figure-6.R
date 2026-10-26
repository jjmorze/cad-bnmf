# ==============================================================================
#  FOUR-PANEL FIGURE
#  a) UKBB: MDWO/MUNW pPRS x waist circumference interaction curve for MI
#  b) UKBB: Total PRS x waist circumference interaction curve for MI
#  c) UKBB: Sex-stratified categorical MDWO/MUNW pPRS x waist circumference HR bars
#  d) SCAPIS: Linear interaction plot for ct_sis ~ mdwo_prs_std x wc
#
#  Main layout fix:
#  - Final assembly uses cowplot.
#  - Panel C and D are aligned at the top only, not at the bottom.
#  - This prevents angled x-axis labels in Panel C from pushing down
#    the x-axis title of Panel D.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(forcats)
library(survival)
library(broom)
library(cowplot)

# ------------------------------------------------------------------------------
#  1. Input data objects
# ------------------------------------------------------------------------------

ukbb_data <- dataset_ukbb

# Change this if your SCAPIS object has a different name
scapis_data <- dataset_scapis

# ------------------------------------------------------------------------------
#  2. Global settings
# ------------------------------------------------------------------------------

pc_covars <- paste0("pc", 1:10)

# Adjustment set for Panels A and B
ukbb_full_covars <- c(
  "age",
  "sex",
  pc_covars,
  "bmi",
  "smoking_status",
  "conven_ldl_therapy_corr",
  "log_conven_crp",
  "log_conven_hba1c",
  "sbp_mean_therapy_corr",
  "conven_rchol_therapy_corr",
  "egfr_ckd_cystatin"
)

# Simple model covariates for Panel C
panel_c_covars <- c(
  "age",
  pc_covars,
  "bmi"
)

# Sex coding
male_value <- "Male"
female_value <- "Female"

# Waist cut-offs for Panel C
wc_cut_men <- 94
wc_cut_women <- 80

# One-sided 95% CI multiplier for Panel C
panel_c_z_ci <- qnorm(0.95)

# Position for Panel A/B Pinteraction labels
panel_ab_p_label_y <- 1.14

# Fixed Panel A/B visible y-axis limits and breaks
panel_ab_y_limits <- c(0.90, 1.15)
panel_ab_y_breaks <- seq(0.90, 1.15, by = 0.05)

# Fixed Panel C y-axis limits
panel_c_y_min <- 0
panel_c_y_max <- 2.4

# Panel C left/right white space
# First value = space left of first column
# Second value = space right of fourth column
panel_c_x_expand <- c(0.75, 0.75)

# Adjustable Panel C x-axis label angle
# Use 0, 30, 45, or 60 depending on label fit
panel_c_x_text_angle <- 45
panel_c_x_text_hjust <- ifelse(panel_c_x_text_angle == 0, 0.5, 1)
panel_c_x_text_vjust <- ifelse(panel_c_x_text_angle == 0, 0.5, 1)

# Colors for Panels A/B/D
mdwo_col <- "#4A7C4E"
mdwo_fill <- "#C8DDBE"

total_col <- "#444444"
total_fill <- "#D0D0D0"

scapis_low_col <- "#8DBA7D"
scapis_high_col <- "#4A7C4E"

# Four ordered shades for Panel C
bar_fill_low_prs_low_wc <- "#C8DDBE"
bar_fill_low_prs_high_wc <- "#4A7C4E"
bar_fill_high_prs_low_wc <- "#C8DDBE"
bar_fill_high_prs_high_wc <- "#4A7C4E"

# Panel-label settings for final cowplot canvas
panel_label_size <- 12
panel_label_x_left <- 0.006
panel_label_x_right <- 0.506
panel_label_y_top <- 0.994
panel_label_y_bottom <- 0.498

# Unicode superscript labels
label_5th <- "5ᵗʰ"
label_95th <- "95ᵗʰ"

# ------------------------------------------------------------------------------
#  3. Prepare UKBB analytic dataset
# ------------------------------------------------------------------------------

ukbb_mi_raw <- ukbb_data %>%
  mutate(
    log_conven_crp = ifelse(
      !is.na(conven_crp) & conven_crp > 0,
      log(conven_crp),
      NA_real_
    ),
    log_conven_hba1c = ifelse(
      !is.na(conven_hba1c) & conven_hba1c > 0,
      log(conven_hba1c),
      NA_real_
    )
  ) %>%
  filter(
    race == "White",
    ascvd_ukb_prev != 1
  )

mdwo_mu <- mean(ukbb_mi_raw$mdwo_prs_std, na.rm = TRUE)
mdwo_sd <- sd(ukbb_mi_raw$mdwo_prs_std, na.rm = TRUE)

total_mu <- mean(ukbb_mi_raw$total_196_prs_std, na.rm = TRUE)
total_sd <- sd(ukbb_mi_raw$total_196_prs_std, na.rm = TRUE)

wc_mu <- mean(ukbb_mi_raw$wc, na.rm = TRUE)
wc_sd <- sd(ukbb_mi_raw$wc, na.rm = TRUE)
wc_median <- median(ukbb_mi_raw$wc, na.rm = TRUE)

ukbb_mi <- ukbb_mi_raw %>%
  mutate(
    mdwo_prs_z = (mdwo_prs_std - mdwo_mu) / mdwo_sd,
    total_prs_z = (total_196_prs_std - total_mu) / total_sd,
    wc_z = (wc - wc_mu) / wc_sd
  )

# ------------------------------------------------------------------------------
#  4. Helper functions
# ------------------------------------------------------------------------------

extract_p_int <- function(fit, interaction_term) {
  
  coef_summ <- summary(fit)$coefficients
  
  if (!interaction_term %in% rownames(coef_summ)) {
    return(NA_real_)
  }
  
  coef_summ[interaction_term, "Pr(>|z|)"]
}

format_p_int_plotmath <- function(p, digits = 3) {
  
  if (is.na(p)) {
    return("italic(P)[interaction] == NA")
  }
  
  if (p < 0.001 && digits <= 3) {
    return("italic(P)[interaction] < 0.001")
  }
  
  paste0(
    "italic(P)[interaction] == ",
    formatC(p, format = "f", digits = digits)
  )
}

get_conditional_prs_hr_from_z_model <- function(
    fit,
    interaction_term,
    biomarker_original_seq,
    biomarker_mu,
    biomarker_sd,
    biomarker_median,
    biomarker_label,
    prs_label
) {
  
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)
  
  if (!interaction_term %in% names(beta)) {
    stop(
      "Interaction term not found in model coefficients: ",
      interaction_term
    )
  }
  
  biomarker_z_seq <- (
    biomarker_original_seq - biomarker_mu
  ) / biomarker_sd
  
  biomarker_median_z <- (
    biomarker_median - biomarker_mu
  ) / biomarker_sd
  
  out <- lapply(
    seq_along(biomarker_original_seq),
    function(i) {
      
      z <- biomarker_z_seq[i]
      
      contrast <- rep(0, length(beta))
      names(contrast) <- names(beta)
      
      contrast[interaction_term] <- z - biomarker_median_z
      
      log_hr <- sum(contrast * beta)
      
      se_log_hr <- sqrt(
        as.numeric(
          t(contrast) %*% vc %*% contrast
        )
      )
      
      data.frame(
        biomarker_original = biomarker_original_seq[i],
        biomarker_z = z,
        biomarker_median = biomarker_median,
        log_hr = log_hr,
        se = se_log_hr,
        hr = exp(log_hr),
        lb = exp(log_hr - 1.96 * se_log_hr),
        ub = exp(log_hr + 1.96 * se_log_hr),
        biomarker = biomarker_label,
        prs = prs_label,
        stringsAsFactors = FALSE
      )
    }
  )
  
  bind_rows(out)
}

# ------------------------------------------------------------------------------
#  5. Panel A model: MDWO/MUNW pPRS x WC
# ------------------------------------------------------------------------------

fit_mdwo_wc_ukbb <- coxph(
  as.formula(
    paste0(
      "Surv(mi_ukb_survtime, mi_ukb_ind) ~ ",
      "mdwo_prs_z * wc_z + ",
      paste(ukbb_full_covars, collapse = " + ")
    )
  ),
  data = ukbb_mi,
  ties = "efron",
  na.action = na.omit
)

p_mdwo_wc <- extract_p_int(
  fit = fit_mdwo_wc_ukbb,
  interaction_term = "mdwo_prs_z:wc_z"
)

p_mdwo_wc_label <- format_p_int_plotmath(
  p = p_mdwo_wc,
  digits = 4
)

# ------------------------------------------------------------------------------
#  6. Panel B model: Total PRS x WC
# ------------------------------------------------------------------------------

fit_total_wc_ukbb <- coxph(
  as.formula(
    paste0(
      "Surv(mi_ukb_survtime, mi_ukb_ind) ~ ",
      "total_prs_z * wc_z + ",
      paste(ukbb_full_covars, collapse = " + ")
    )
  ),
  data = ukbb_mi,
  ties = "efron",
  na.action = na.omit
)

p_total_wc <- extract_p_int(
  fit = fit_total_wc_ukbb,
  interaction_term = "total_prs_z:wc_z"
)

p_total_wc_label <- format_p_int_plotmath(
  p = p_total_wc,
  digits = 3
)

# ------------------------------------------------------------------------------
#  7. Generate Panel A/B curves
# ------------------------------------------------------------------------------

wc_grid <- seq(
  quantile(
    ukbb_mi$wc,
    0.01,
    na.rm = TRUE
  ),
  quantile(
    ukbb_mi$wc,
    0.99,
    na.rm = TRUE
  ),
  length.out = 100
)

mdwo_wc_curve <- get_conditional_prs_hr_from_z_model(
  fit = fit_mdwo_wc_ukbb,
  interaction_term = "mdwo_prs_z:wc_z",
  biomarker_original_seq = wc_grid,
  biomarker_mu = wc_mu,
  biomarker_sd = wc_sd,
  biomarker_median = wc_median,
  biomarker_label = "Waist circumference",
  prs_label = "MUNW pPRS"
)

total_wc_curve <- get_conditional_prs_hr_from_z_model(
  fit = fit_total_wc_ukbb,
  interaction_term = "total_prs_z:wc_z",
  biomarker_original_seq = wc_grid,
  biomarker_mu = wc_mu,
  biomarker_sd = wc_sd,
  biomarker_median = wc_median,
  biomarker_label = "Waist circumference",
  prs_label = "Total PRS"
)

# ------------------------------------------------------------------------------
#  8. Shared theme for Panel A/B line plots
# ------------------------------------------------------------------------------

line_panel_theme <- theme_bw(
  base_size = 11
) +
  theme(
    plot.title = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.25
    ),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.35
    )
  )

# ------------------------------------------------------------------------------
#  9. Panel A plot
# ------------------------------------------------------------------------------

panel_a <- ggplot(
  mdwo_wc_curve,
  aes(
    x = biomarker_original,
    y = hr
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_vline(
    xintercept = wc_median,
    linetype = "dotted",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_ribbon(
    aes(
      ymin = lb,
      ymax = ub
    ),
    fill = mdwo_fill,
    alpha = 0.45
  ) +
  geom_line(
    color = mdwo_col,
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = panel_ab_p_label_y,
    label = p_mdwo_wc_label,
    hjust = 1.05,
    vjust = 0.5,
    size = 3.4,
    parse = TRUE
  ) +
  scale_y_log10(
    breaks = panel_ab_y_breaks,
    labels = sprintf(
      "%.2f",
      panel_ab_y_breaks
    )
  ) +
  coord_cartesian(
    ylim = panel_ab_y_limits
  ) +
  labs(
    x = "Waist circumference, cm",
    y = "HR for 1 s.d. higher MUNW pPRS"
  ) +
  line_panel_theme

# ------------------------------------------------------------------------------
#  10. Panel B plot
# ------------------------------------------------------------------------------

panel_b <- ggplot(
  total_wc_curve,
  aes(
    x = biomarker_original,
    y = hr
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_vline(
    xintercept = wc_median,
    linetype = "dotted",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_ribbon(
    aes(
      ymin = lb,
      ymax = ub
    ),
    fill = total_fill,
    alpha = 0.45
  ) +
  geom_line(
    color = total_col,
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = panel_ab_p_label_y,
    label = p_total_wc_label,
    hjust = 1.05,
    vjust = 0.5,
    size = 3.4,
    parse = TRUE
  ) +
  scale_y_log10(
    breaks = panel_ab_y_breaks,
    labels = sprintf(
      "%.2f",
      panel_ab_y_breaks
    )
  ) +
  coord_cartesian(
    ylim = panel_ab_y_limits
  ) +
  labs(
    x = "Waist circumference, cm",
    y = "HR for 1 s.d. higher total PRS"
  ) +
  line_panel_theme

# ------------------------------------------------------------------------------
#  11. Panel C data/model helper
# ------------------------------------------------------------------------------

make_panel_c_data_by_sex <- function(
    data,
    sex_value,
    sex_label,
    wc_cutoff,
    time_var = "mi_ukb_survtime",
    event_var = "mi_ukb_ind"
) {
  
  dat_sex_raw <- data %>%
    filter(
      sex == sex_value
    )
  
  prs_p05 <- quantile(
    dat_sex_raw$mdwo_prs_std,
    0.05,
    na.rm = TRUE
  )
  
  prs_p95 <- quantile(
    dat_sex_raw$mdwo_prs_std,
    0.95,
    na.rm = TRUE
  )
  
  dat_sex <- dat_sex_raw %>%
    mutate(
      mdwo_wc_group = case_when(
        mdwo_prs_std < prs_p05 &
          wc < wc_cutoff ~
          "Low pPRS / low WC",
        
        mdwo_prs_std < prs_p05 &
          wc >= wc_cutoff ~
          "Low pPRS / high WC",
        
        mdwo_prs_std > prs_p95 &
          wc < wc_cutoff ~
          "High pPRS / low WC",
        
        mdwo_prs_std > prs_p95 &
          wc >= wc_cutoff ~
          "High pPRS / high WC",
        
        TRUE ~ NA_character_
      ),
      mdwo_wc_group = factor(
        mdwo_wc_group,
        levels = c(
          "Low pPRS / low WC",
          "Low pPRS / high WC",
          "High pPRS / low WC",
          "High pPRS / high WC"
        )
      )
    ) %>%
    filter(
      !is.na(mdwo_wc_group)
    ) %>%
    select(
      all_of(
        c(
          time_var,
          event_var,
          "mdwo_wc_group",
          panel_c_covars
        )
      )
    ) %>%
    na.omit()
  
  model_formula <- as.formula(
    paste0(
      "Surv(",
      time_var,
      ", ",
      event_var,
      ") ~ ",
      "mdwo_wc_group + ",
      paste(
        panel_c_covars,
        collapse = " + "
      )
    )
  )
  
  fit <- coxph(
    formula = model_formula,
    data = dat_sex,
    ties = "efron"
  )
  
  coef_summ <- summary(fit)$coefficients
  
  ref_row <- data.frame(
    sex = sex_label,
    wc_cutoff = wc_cutoff,
    prs_p05 = prs_p05,
    prs_p95 = prs_p95,
    group = "Low pPRS / low WC",
    beta = 0,
    se = NA_real_,
    hr = 1,
    lb = NA_real_,
    ub = NA_real_,
    p_value = NA_real_,
    n = sum(
      dat_sex$mdwo_wc_group ==
        "Low pPRS / low WC",
      na.rm = TRUE
    ),
    events = sum(
      dat_sex[[event_var]] == 1 &
        dat_sex$mdwo_wc_group ==
        "Low pPRS / low WC",
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  nonref_terms <- rownames(coef_summ)[
    grepl(
      "^mdwo_wc_group",
      rownames(coef_summ)
    )
  ]
  
  nonref_rows <- lapply(
    nonref_terms,
    function(term_i) {
      
      group_i <- gsub(
        "^mdwo_wc_group",
        "",
        term_i
      )
      
      beta_i <- coef_summ[
        term_i,
        "coef"
      ]
      
      se_i <- coef_summ[
        term_i,
        "se(coef)"
      ]
      
      data.frame(
        sex = sex_label,
        wc_cutoff = wc_cutoff,
        prs_p05 = prs_p05,
        prs_p95 = prs_p95,
        group = group_i,
        beta = beta_i,
        se = se_i,
        hr = exp(beta_i),
        lb = exp(
          beta_i -
            panel_c_z_ci * se_i
        ),
        ub = exp(
          beta_i +
            panel_c_z_ci * se_i
        ),
        p_value = coef_summ[
          term_i,
          "Pr(>|z|)"
        ],
        n = sum(
          dat_sex$mdwo_wc_group ==
            group_i,
          na.rm = TRUE
        ),
        events = sum(
          dat_sex[[event_var]] == 1 &
            dat_sex$mdwo_wc_group ==
            group_i,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  ) %>%
    bind_rows()
  
  out <- bind_rows(
    ref_row,
    nonref_rows
  ) %>%
    mutate(
      sex = factor(
        sex,
        levels = c(
          "Men",
          "Women"
        )
      ),
      group = factor(
        group,
        levels = c(
          "Low pPRS / low WC",
          "Low pPRS / high WC",
          "High pPRS / low WC",
          "High pPRS / high WC"
        )
      ),
      group_label = recode(
        as.character(group),
        "Low pPRS / low WC" =
          paste0(
            "<",
            label_5th,
            " pPRS\nlow WC"
          ),
        "Low pPRS / high WC" =
          paste0(
            "<",
            label_5th,
            " pPRS\nhigh WC"
          ),
        "High pPRS / low WC" =
          paste0(
            ">",
            label_95th,
            " pPRS\nlow WC"
          ),
        "High pPRS / high WC" =
          paste0(
            ">",
            label_95th,
            " pPRS\nhigh WC"
          )
      ),
      group_label = factor(
        group_label,
        levels = c(
          paste0(
            "<",
            label_5th,
            " pPRS\nlow WC"
          ),
          paste0(
            "<",
            label_5th,
            " pPRS\nhigh WC"
          ),
          paste0(
            ">",
            label_95th,
            " pPRS\nlow WC"
          ),
          paste0(
            ">",
            label_95th,
            " pPRS\nhigh WC"
          )
        )
      ),
      fill_group = case_when(
        group == "Low pPRS / low WC" ~
          "<5th pPRS / low WC",
        group == "Low pPRS / high WC" ~
          "<5th pPRS / high WC",
        group == "High pPRS / low WC" ~
          ">95th pPRS / low WC",
        group == "High pPRS / high WC" ~
          ">95th pPRS / high WC"
      ),
      hr_label = ifelse(
        group == "Low pPRS / low WC",
        "Ref.",
        paste0(
          sprintf("%.2f", hr),
          "\n(",
          sprintf("%.2f", lb),
          "-",
          sprintf("%.2f", ub),
          ")"
        )
      )
    )
  
  list(
    fit = fit,
    data = dat_sex,
    results = out,
    prs_p05 = prs_p05,
    prs_p95 = prs_p95
  )
}

panel_c_men <- make_panel_c_data_by_sex(
  data = ukbb_mi,
  sex_value = male_value,
  sex_label = "Men",
  wc_cutoff = wc_cut_men
)

panel_c_women <- make_panel_c_data_by_sex(
  data = ukbb_mi,
  sex_value = female_value,
  sex_label = "Women",
  wc_cutoff = wc_cut_women
)

panel_c_dat <- bind_rows(
  panel_c_men$results,
  panel_c_women$results
)

panel_c_label_offset <- 0.06 * panel_c_y_max

panel_c_dat <- panel_c_dat %>%
  mutate(
    label_y = case_when(
      group == "Low pPRS / low WC" ~
        hr + panel_c_label_offset,
      
      !is.na(ub) ~
        ub + panel_c_label_offset,
      
      TRUE ~ NA_real_
    ),
    label_y = pmin(
      label_y,
      panel_c_y_max * 0.94
    )
  )

panel_c_breaks <- pretty(
  c(
    panel_c_y_min,
    panel_c_y_max
  ),
  n = 5
)

panel_c_breaks <- panel_c_breaks[
  panel_c_breaks > 0 &
    panel_c_breaks <= panel_c_y_max
]

# ------------------------------------------------------------------------------
#  12. Panel C plot
# ------------------------------------------------------------------------------

panel_c <- ggplot(
  panel_c_dat,
  aes(
    x = group_label,
    y = hr,
    fill = fill_group
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_col(
    width = 0.58,
    color = "black",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(
      ymin = lb,
      ymax = ub
    ),
    width = 0.10,
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  geom_text(
    aes(
      y = label_y,
      label = hr_label
    ),
    size = 2.55,
    lineheight = 0.88,
    vjust = 0,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "<5th pPRS / low WC" =
        bar_fill_low_prs_low_wc,
      
      "<5th pPRS / high WC" =
        bar_fill_low_prs_high_wc,
      
      ">95th pPRS / low WC" =
        bar_fill_high_prs_low_wc,
      
      ">95th pPRS / high WC" =
        bar_fill_high_prs_high_wc
    )
  ) +
  scale_x_discrete(
    expand = expansion(
      add = panel_c_x_expand
    )
  ) +
  scale_y_continuous(
    limits = c(
      panel_c_y_min,
      panel_c_y_max
    ),
    breaks = panel_c_breaks,
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  facet_grid(
    . ~ sex
  ) +
  labs(
    x = NULL,
    y = "HR for MI"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.25
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.35
    ),
    strip.background = element_rect(
      fill = "black",
      color = "black",
      linewidth = 0.35
    ),
    strip.text = element_text(
      face = "bold",
      color = "white"
    ),
    axis.text.x = element_text(
      color = "black",
      size = 8,
      angle = panel_c_x_text_angle,
      hjust = panel_c_x_text_hjust,
      vjust = panel_c_x_text_vjust
    ),
    plot.margin = margin(
      t = 5.5,
      r = 12,
      b = 5.5,
      l = 5.5
    )
  )

# ------------------------------------------------------------------------------
#  13. Prepare SCAPIS data for Panel D
# ------------------------------------------------------------------------------

scapis_lm_dat <- scapis_data %>%
  select(
    ct_sis,
    age,
    sex,
    mdwo_prs_std,
    wc,
    all_of(pc_covars),
    bmi
  ) %>%
  na.omit()

# ------------------------------------------------------------------------------
#  14. Fit SCAPIS linear interaction model
# ------------------------------------------------------------------------------

fit_scapis_mdwo_wc <- lm(
  ct_sis ~
    age +
    sex +
    mdwo_prs_std * wc +
    pc1 +
    pc2 +
    pc3 +
    pc4 +
    pc5 +
    pc6 +
    pc7 +
    pc8 +
    pc9 +
    pc10 +
    bmi,
  data = scapis_lm_dat
)

p_scapis_int <- summary(
  fit_scapis_mdwo_wc
)$coefficients[
  "mdwo_prs_std:wc",
  "Pr(>|t|)"
]

p_scapis_int_label <- format_p_int_plotmath(
  p = p_scapis_int,
  digits = 3
)

# ------------------------------------------------------------------------------
#  15. Generate SCAPIS prediction grid
# ------------------------------------------------------------------------------

scapis_wc_grid <- seq(
  quantile(
    scapis_lm_dat$wc,
    0.01,
    na.rm = TRUE
  ),
  quantile(
    scapis_lm_dat$wc,
    0.99,
    na.rm = TRUE
  ),
  length.out = 100
)

scapis_mdwo_p05 <- quantile(
  scapis_lm_dat$mdwo_prs_std,
  0.05,
  na.rm = TRUE
)

scapis_mdwo_p95 <- quantile(
  scapis_lm_dat$mdwo_prs_std,
  0.95,
  na.rm = TRUE
)

get_reference_value <- function(x) {
  
  if (is.factor(x)) {
    return(
      levels(x)[1]
    )
  }
  
  if (is.character(x)) {
    return(
      names(
        sort(
          table(x),
          decreasing = TRUE
        )
      )[1]
    )
  }
  
  mean(
    x,
    na.rm = TRUE
  )
}

scapis_ref <- scapis_lm_dat %>%
  summarise(
    age = mean(
      age,
      na.rm = TRUE
    ),
    bmi = mean(
      bmi,
      na.rm = TRUE
    ),
    sex = get_reference_value(
      sex
    ),
    pc1 = mean(
      pc1,
      na.rm = TRUE
    ),
    pc2 = mean(
      pc2,
      na.rm = TRUE
    ),
    pc3 = mean(
      pc3,
      na.rm = TRUE
    ),
    pc4 = mean(
      pc4,
      na.rm = TRUE
    ),
    pc5 = mean(
      pc5,
      na.rm = TRUE
    ),
    pc6 = mean(
      pc6,
      na.rm = TRUE
    ),
    pc7 = mean(
      pc7,
      na.rm = TRUE
    ),
    pc8 = mean(
      pc8,
      na.rm = TRUE
    ),
    pc9 = mean(
      pc9,
      na.rm = TRUE
    ),
    pc10 = mean(
      pc10,
      na.rm = TRUE
    )
  )

scapis_pred_grid <- expand.grid(
  wc = scapis_wc_grid,
  mdwo_prs_std = c(
    scapis_mdwo_p05,
    scapis_mdwo_p95
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    mdwo_group = ifelse(
      mdwo_prs_std ==
        scapis_mdwo_p05,
      
      paste0(
        "<",
        label_5th,
        " centile MUNW pPRS"
      ),
      
      paste0(
        ">",
        label_95th,
        " centile MUNW pPRS"
      )
    )
  ) %>%
  bind_cols(
    scapis_ref[
      rep(
        1,
        nrow(.)
      ),
    ]
  )

if (
  is.factor(
    scapis_lm_dat$sex
  )
) {
  
  scapis_pred_grid$sex <- factor(
    scapis_pred_grid$sex,
    levels = levels(
      scapis_lm_dat$sex
    )
  )
}

scapis_pred <- predict(
  fit_scapis_mdwo_wc,
  newdata = scapis_pred_grid,
  se.fit = TRUE
)

scapis_pred_grid <- scapis_pred_grid %>%
  mutate(
    fit = as.numeric(
      scapis_pred$fit
    ),
    se = as.numeric(
      scapis_pred$se.fit
    ),
    lb = fit - 1.96 * se,
    ub = fit + 1.96 * se,
    mdwo_group = factor(
      mdwo_group,
      levels = c(
        paste0(
          "<",
          label_5th,
          " centile MUNW pPRS"
        ),
        paste0(
          ">",
          label_95th,
          " centile MUNW pPRS"
        )
      )
    )
  )

# ------------------------------------------------------------------------------
#  16. Panel D plot
# ------------------------------------------------------------------------------

panel_d <- ggplot(
  scapis_pred_grid,
  aes(
    x = wc,
    y = fit,
    color = mdwo_group,
    fill = mdwo_group
  )
) +
  geom_ribbon(
    aes(
      ymin = lb,
      ymax = ub
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = p_scapis_int_label,
    hjust = 1.05,
    vjust = 1.5,
    size = 3.4,
    parse = TRUE
  ) +
  scale_color_manual(
    values = setNames(
      c(
        scapis_low_col,
        scapis_high_col
      ),
      c(
        paste0(
          "<",
          label_5th,
          " centile MUNW pPRS"
        ),
        paste0(
          ">",
          label_95th,
          " centile MUNW pPRS"
        )
      )
    )
  ) +
  scale_fill_manual(
    values = setNames(
      c(
        scapis_low_col,
        scapis_high_col
      ),
      c(
        paste0(
          "<",
          label_5th,
          " centile MUNW pPRS"
        ),
        paste0(
          ">",
          label_95th,
          " centile MUNW pPRS"
        )
      )
    )
  ) +
  labs(
    x = "Waist circumference, cm",
    y = "Predicted segment involvement score",
    color = NULL,
    fill = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.25
    ),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.35
    ),
    plot.margin = margin(
      t = 5.5,
      r = 5.5,
      b = 5.5,
      l = 5.5
    )
  )

# ------------------------------------------------------------------------------
#  17. Combine with cowplot
#
#  Critical layout point:
#  - Panel C and D are aligned horizontally by the top axis only: axis = "t".
#  - This avoids forcing Panel D to adopt Panel C's larger bottom space
#    caused by angled x-axis labels.
# ------------------------------------------------------------------------------

top_row <- cowplot::plot_grid(
  panel_a,
  panel_b,
  nrow = 1,
  rel_widths = c(
    1,
    1
  ),
  align = "hv",
  axis = "tblr"
)

bottom_row <- cowplot::plot_grid(
  panel_c,
  panel_d,
  nrow = 1,
  rel_widths = c(
    1,
    1
  ),
  align = "h",
  axis = "t"
)

fig_body <- cowplot::plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(
    1,
    1.08
  ),
  align = "v",
  axis = "lr"
)

fig_mdwo_wc_four_panel <- cowplot::ggdraw() +
  cowplot::draw_plot(
    fig_body,
    x = 0,
    y = 0,
    width = 1,
    height = 1
  ) +
  cowplot::draw_label(
    "a",
    x = panel_label_x_left,
    y = panel_label_y_top,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = panel_label_size
  ) +
  cowplot::draw_label(
    "b",
    x = panel_label_x_right,
    y = panel_label_y_top,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = panel_label_size
  ) +
  cowplot::draw_label(
    "c",
    x = panel_label_x_left,
    y = panel_label_y_bottom,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = panel_label_size
  ) +
  cowplot::draw_label(
    "d",
    x = panel_label_x_right,
    y = panel_label_y_bottom,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = panel_label_size
  )

fig_mdwo_wc_four_panel

# ------------------------------------------------------------------------------
#  18. Save figure
# ------------------------------------------------------------------------------

ggsave(
  filename = "figure-6.png",
  plot = fig_mdwo_wc_four_panel,
  width = 8,
  height = 7,
  bg = "white",
  dpi = 1200
)

# ------------------------------------------------------------------------------
#  19. Optional: print model summaries and key results
# ------------------------------------------------------------------------------

summary(
  fit_mdwo_wc_ukbb
)

summary(
  fit_total_wc_ukbb
)

summary(
  panel_c_men$fit
)

summary(
  panel_c_women$fit
)

summary(
  fit_scapis_mdwo_wc
)

panel_c_dat

data.frame(
  dataset = "SCAPIS",
  group = c(
    "<5th centile MUNW pPRS",
    ">95th centile MUNW pPRS"
  ),
  mdwo_prs_std_value = c(
    scapis_mdwo_p05,
    scapis_mdwo_p95
  )
)