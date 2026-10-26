# ==============================================================================
#  ONE-CANVAS MULTIPANEL FIGURE
#  a) Stage 1 forest plot
#  b) Fully adjusted conditional HR curves
#
#  Current specifications:
#  - Patchwork one-canvas figure
#  - Panel A: Stage 1 PRS × biomarker interaction forest plot
#  - Panel B: LDL panel left, TRL/remnant cholesterol panel right
#  - No panel c
#  - Manual bold panel labels a and b, font size 18
#  - P[interaction] labels with subscript, placed at HR = 1.06
#  - Native LDL/ChEL y-axis title outside plot panel
#  - Panel B1 left side freed from Panel A alignment
#  - Panel A vertical gridlines black and drawn above grey bands
#  - Non-significant estimates grey60
# ==============================================================================

library(dplyr)
library(ggplot2)
library(forcats)
library(survival)
library(patchwork)

# ------------------------------------------------------------------------------
#  1. Inputs
# ------------------------------------------------------------------------------

stage1_results <- ukbb_cad_interactions_stage1

prs_order <- c(
  "ChEL pPRS",
  "TRL pPRS",
  "Lp(a) pPRS",
  "MUNW pPRS",
  "Obesity pPRS",
  "Endothelial dysfunction pPRS",
  "Cardiac remodeling pPRS",
  "Coagulation pPRS"
)

prs_vars_no_total <- c(
  "chol_prs_std",
  "trl_prs_std",
  "lpa_prs_std",
  "mdwo_prs_std",
  "obesity_prs_std",
  "endothelial_prs_std",
  "cardiac_prs_std",
  "coagulation_prs_std"
)

prs_labels_no_total <- c(
  chol_prs_std        = "ChEL pPRS",
  trl_prs_std         = "TRL pPRS",
  lpa_prs_std         = "Lp(a) pPRS",
  mdwo_prs_std        = "MUNW pPRS",
  obesity_prs_std     = "Obesity pPRS",
  endothelial_prs_std = "Endothelial dysfunction pPRS",
  cardiac_prs_std     = "Cardiac remodeling pPRS",
  coagulation_prs_std = "Coagulation pPRS"
)

prs_colors <- c(
  chol_prs_std        = "#B22D4C",
  trl_prs_std         = "#E07A8E",
  lpa_prs_std         = "#7A1F3A",
  mdwo_prs_std        = "#4A7C4E",
  obesity_prs_std     = "#95C075",
  endothelial_prs_std = "#2D5F8B",
  cardiac_prs_std     = "#7AAED1",
  coagulation_prs_std = "#5D4787"
)

panel_a_interactions <- c(
  "conven_ldl_therapy_corr",
  "log_conven_crp",
  "log_conven_hba1c",
  "sbp_mean_therapy_corr",
  "wc",
  "conven_rchol_therapy_corr"
)

panel_a_interaction_labels <- c(
  conven_ldl_therapy_corr   = "LDL-C",
  log_conven_crp            = "log(CRP)",
  log_conven_hba1c          = "log(HbA1c)",
  sbp_mean_therapy_corr     = "SBP",
  wc                        = "WC",
  conven_rchol_therapy_corr = "RChol"
)

interaction_p_threshold <- 0.05 / 8

# Panel A shared border width for strips and plot boxes
panel_a_box_lwd <- 0.35

# ------------------------------------------------------------------------------
#  2. Panel A data
# ------------------------------------------------------------------------------

panel_a_dat <- stage1_results %>%
  filter(
    predictor_short %in% prs_vars_no_total,
    interaction_short %in% panel_a_interactions
  ) %>%
  mutate(
    predictor_long = recode(
      predictor_short,
      !!!prs_labels_no_total
    ),
    interaction_long = recode(
      interaction_short,
      !!!panel_a_interaction_labels
    ),
    predictor_long = factor(
      predictor_long,
      levels = rev(prs_order)
    ),
    interaction_long = factor(
      interaction_long,
      levels = unname(panel_a_interaction_labels)
    ),
    is_significant = !is.na(p_value) &
      p_value < interaction_p_threshold,
    plot_color = ifelse(
      is_significant,
      prs_colors[predictor_short],
      "grey60"
    )
  )

panel_a_bands <- expand.grid(
  predictor_long = levels(panel_a_dat$predictor_long),
  interaction_long = levels(panel_a_dat$interaction_long),
  stringsAsFactors = FALSE
) %>%
  mutate(
    predictor_long = factor(
      predictor_long,
      levels = levels(panel_a_dat$predictor_long)
    ),
    interaction_long = factor(
      interaction_long,
      levels = levels(panel_a_dat$interaction_long)
    ),
    y_id = as.numeric(predictor_long),
    band = y_id %% 2 == 1
  ) %>%
  filter(
    band
  )

# ------------------------------------------------------------------------------
#  3. Panel A forest plot
# ------------------------------------------------------------------------------

panel_a <- ggplot(
  panel_a_dat,
  aes(
    x = hr,
    y = predictor_long
  )
) +
  geom_rect(
    data = panel_a_bands,
    aes(
      ymin = y_id - 0.5,
      ymax = y_id + 0.5,
      xmin = -Inf,
      xmax = Inf
    ),
    inherit.aes = FALSE,
    fill = "grey95",
    color = NA
  ) +
  geom_vline(
    xintercept = c(
      0.96,
      1.04
    ),
    color = "black",
    linewidth = 0.25
  ) +
  geom_vline(
    xintercept = 1.00,
    color = "black",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_errorbarh(
    aes(
      xmin = lb,
      xmax = ub,
      color = plot_color
    ),
    height = 0,
    linewidth = 0.45
  ) +
  geom_point(
    aes(
      fill = plot_color,
      color = plot_color
    ),
    shape = 22,
    size = 2.4,
    stroke = 0.25
  ) +
  scale_x_continuous(
    limits = c(
      0.955,
      1.045
    ),
    breaks = c(
      0.96,
      1.00,
      1.04
    ),
    labels = c(
      "0.96",
      "1.00",
      "1.04"
    )
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  facet_grid(
    . ~ interaction_long,
    scales = "fixed",
    space = "fixed"
  ) +
  labs(
    x = "Interaction HR per 1 s.d. pPRS × 1 s.d. biomarker",
    y = NULL,
    tag = "a"
  ) +
  coord_cartesian(
    xlim = c(
      0.955,
      1.045
    ),
    clip = "off"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    plot.title = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(
      1.1,
      "lines"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = panel_a_box_lwd
    ),
    strip.background = element_rect(
      fill = "black",
      color = "black",
      linewidth = panel_a_box_lwd
    ),
    strip.text = element_text(
      face = "bold",
      color = "white"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    plot.tag = element_text(
      face = "bold",
      size = 18,
      color = "black"
    ),
    plot.tag.position = c(
      0,
      1
    ),
    plot.margin = margin(
      5.5,
      5.5,
      5.5,
      5.5
    )
  )

# ------------------------------------------------------------------------------
#  4. Prepare MI dataset
# ------------------------------------------------------------------------------

cox_dat_mi_raw <- dataset_ukbb %>%
  mutate(
    log_conven_crp = ifelse(
      !is.na(conven_crp) &
        conven_crp > 0,
      log(conven_crp),
      NA_real_
    ),
    log_conven_hba1c = ifelse(
      !is.na(conven_hba1c) &
        conven_hba1c > 0,
      log(conven_hba1c),
      NA_real_
    )
  ) %>%
  filter(
    race == "White",
    ascvd_ukb_prev != 1
  )

# Scaling parameters from the analytic subset

trl_mu <- mean(
  cox_dat_mi_raw$trl_prs_std,
  na.rm = TRUE
)

trl_sd <- sd(
  cox_dat_mi_raw$trl_prs_std,
  na.rm = TRUE
)

chol_mu <- mean(
  cox_dat_mi_raw$chol_prs_std,
  na.rm = TRUE
)

chol_sd <- sd(
  cox_dat_mi_raw$chol_prs_std,
  na.rm = TRUE
)

rchol_mu <- mean(
  cox_dat_mi_raw$conven_rchol_therapy_corr,
  na.rm = TRUE
)

rchol_sd <- sd(
  cox_dat_mi_raw$conven_rchol_therapy_corr,
  na.rm = TRUE
)

rchol_median <- median(
  cox_dat_mi_raw$conven_rchol_therapy_corr,
  na.rm = TRUE
)

ldl_mu <- mean(
  cox_dat_mi_raw$conven_ldl_therapy_corr,
  na.rm = TRUE
)

ldl_sd <- sd(
  cox_dat_mi_raw$conven_ldl_therapy_corr,
  na.rm = TRUE
)

ldl_median <- median(
  cox_dat_mi_raw$conven_ldl_therapy_corr,
  na.rm = TRUE
)

cox_dat_mi <- cox_dat_mi_raw %>%
  mutate(
    trl_prs_z = (
      trl_prs_std - trl_mu
    ) / trl_sd,
    
    chol_prs_z = (
      chol_prs_std - chol_mu
    ) / chol_sd,
    
    conven_rchol_z = (
      conven_rchol_therapy_corr -
        rchol_mu
    ) / rchol_sd,
    
    conven_ldl_z = (
      conven_ldl_therapy_corr -
        ldl_mu
    ) / ldl_sd
  )

# ------------------------------------------------------------------------------
#  5. Fully adjusted Cox model:
#     TRL pPRS × therapy-corrected remnant cholesterol
# ------------------------------------------------------------------------------

fit_trl_rchol_full <- coxph(
  Surv(
    mi_ukb_survtime,
    mi_ukb_ind
  ) ~
    trl_prs_z * conven_rchol_z +
    age +
    sex +
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
    bmi +
    smoking_status +
    conven_ldl_therapy_corr +
    log_conven_crp +
    log_conven_hba1c +
    sbp_mean_therapy_corr +
    egfr_ckd_cystatin +
    wc,
  data = cox_dat_mi,
  na.action = na.omit
)

# ------------------------------------------------------------------------------
#  6. Fully adjusted Cox model:
#     ChEL pPRS × therapy-corrected LDL-C
# ------------------------------------------------------------------------------

fit_chel_ldl_full <- coxph(
  Surv(
    mi_ukb_survtime,
    mi_ukb_ind
  ) ~
    chol_prs_z * conven_ldl_z +
    age +
    sex +
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
    bmi +
    smoking_status +
    conven_rchol_therapy_corr +
    log_conven_crp +
    log_conven_hba1c +
    egfr_ckd_cystatin +
    sbp_mean_therapy_corr +
    wc,
  data = cox_dat_mi,
  na.action = na.omit
)

# ------------------------------------------------------------------------------
#  7. Extract interaction P-values
# ------------------------------------------------------------------------------

extract_p_int <- function(
    fit,
    interaction_term
) {
  
  coef_summ <- summary(
    fit
  )$coefficients
  
  if (
    !interaction_term %in%
    rownames(coef_summ)
  ) {
    return(
      NA_real_
    )
  }
  
  coef_summ[
    interaction_term,
    "Pr(>|z|)"
  ]
}

format_p_int_plotmath <- function(p) {
  
  if (
    is.na(p)
  ) {
    return(
      "italic(P)[interaction] == NA"
    )
  }
  
  if (
    p < 0.001
  ) {
    return(
      "italic(P)[interaction] < 0.001"
    )
  }
  
  paste0(
    "italic(P)[interaction] == ",
    formatC(
      p,
      format = "f",
      digits = 2
    )
  )
}

p_trl_rchol <- extract_p_int(
  fit = fit_trl_rchol_full,
  interaction_term =
    "trl_prs_z:conven_rchol_z"
)

p_chel_ldl <- extract_p_int(
  fit = fit_chel_ldl_full,
  interaction_term =
    "chol_prs_z:conven_ldl_z"
)

p_trl_rchol_label <-
  format_p_int_plotmath(
    p_trl_rchol
  )

p_chel_ldl_label <-
  format_p_int_plotmath(
    p_chel_ldl
  )

# ------------------------------------------------------------------------------
#  8. Helper:
#     Conditional HR for 1-SD PRS across biomarker levels.
#
#     Model scale:
#       z(PRS) * z(biomarker)
#
#     Plot scale:
#       x-axis uses original biomarker values.
#
#     Reference:
#       HR = 1 at median biomarker value.
#
#     Therefore:
#       log HR relative to median =
#       beta_interaction * (z_x - z_median)
# ------------------------------------------------------------------------------

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
  
  beta <- stats::coef(
    fit
  )
  
  vc <- stats::vcov(
    fit
  )
  
  if (
    !interaction_term %in%
    names(beta)
  ) {
    stop(
      "Interaction term not found in model coefficients: ",
      interaction_term
    )
  }
  
  biomarker_z_seq <- (
    biomarker_original_seq -
      biomarker_mu
  ) / biomarker_sd
  
  biomarker_median_z <- (
    biomarker_median -
      biomarker_mu
  ) / biomarker_sd
  
  out <- lapply(
    seq_along(
      biomarker_original_seq
    ),
    function(i) {
      
      z <- biomarker_z_seq[i]
      
      contrast <- rep(
        0,
        length(beta)
      )
      
      names(
        contrast
      ) <- names(
        beta
      )
      
      contrast[
        interaction_term
      ] <- z - biomarker_median_z
      
      log_hr <- sum(
        contrast * beta
      )
      
      se_log_hr <- sqrt(
        as.numeric(
          t(contrast) %*%
            vc %*%
            contrast
        )
      )
      
      data.frame(
        biomarker_original =
          biomarker_original_seq[i],
        
        biomarker_z =
          z,
        
        biomarker_median =
          biomarker_median,
        
        log_hr =
          log_hr,
        
        se =
          se_log_hr,
        
        hr =
          exp(log_hr),
        
        lb =
          exp(
            log_hr -
              1.96 * se_log_hr
          ),
        
        ub =
          exp(
            log_hr +
              1.96 * se_log_hr
          ),
        
        biomarker =
          biomarker_label,
        
        prs =
          prs_label,
        
        stringsAsFactors =
          FALSE
      )
    }
  )
  
  bind_rows(
    out
  )
}

# ------------------------------------------------------------------------------
#  9. Generate Panel B curves
# ------------------------------------------------------------------------------

rchol_grid <- seq(
  quantile(
    cox_dat_mi$conven_rchol_therapy_corr,
    0.01,
    na.rm = TRUE
  ),
  quantile(
    cox_dat_mi$conven_rchol_therapy_corr,
    0.99,
    na.rm = TRUE
  ),
  length.out = 100
)

ldl_grid <- seq(
  quantile(
    cox_dat_mi$conven_ldl_therapy_corr,
    0.01,
    na.rm = TRUE
  ),
  quantile(
    cox_dat_mi$conven_ldl_therapy_corr,
    0.99,
    na.rm = TRUE
  ),
  length.out = 100
)

trl_rchol_curve <-
  get_conditional_prs_hr_from_z_model(
    fit =
      fit_trl_rchol_full,
    
    interaction_term =
      "trl_prs_z:conven_rchol_z",
    
    biomarker_original_seq =
      rchol_grid,
    
    biomarker_mu =
      rchol_mu,
    
    biomarker_sd =
      rchol_sd,
    
    biomarker_median =
      rchol_median,
    
    biomarker_label =
      "Remnant cholesterol",
    
    prs_label =
      "TRL pPRS"
  )

chel_ldl_curve <-
  get_conditional_prs_hr_from_z_model(
    fit =
      fit_chel_ldl_full,
    
    interaction_term =
      "chol_prs_z:conven_ldl_z",
    
    biomarker_original_seq =
      ldl_grid,
    
    biomarker_mu =
      ldl_mu,
    
    biomarker_sd =
      ldl_sd,
    
    biomarker_median =
      ldl_median,
    
    biomarker_label =
      "LDL-C",
    
    prs_label =
      "ChEL pPRS"
  )

# ------------------------------------------------------------------------------
#  10. Panel B aesthetics
# ------------------------------------------------------------------------------

rchol_col <- "#E07A8E"
rchol_fill <- "#F4C5CD"

ldl_col <- "#B22D4C"
ldl_fill <- "#E6A6B5"

panel_b_y_limits <- range(
  c(
    trl_rchol_curve$lb,
    trl_rchol_curve$ub,
    chel_ldl_curve$lb,
    chel_ldl_curve$ub,
    1.06
  ),
  na.rm = TRUE
)

panel_b_y_limits <- c(
  floor(
    panel_b_y_limits[1] *
      100
  ) / 100,
  
  ceiling(
    panel_b_y_limits[2] *
      100
  ) / 100
)

panel_b_theme <- theme_bw(
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
    ),
    plot.margin = margin(
      5.5,
      5.5,
      5.5,
      5.5
    )
  )

# ------------------------------------------------------------------------------
#  11. Panel B1: ChEL pPRS by LDL-C
#      LEFT panel
#
#      Native y-axis title is used.
#      patchwork::free() will prevent Panel A's long y-axis labels from
#      forcing this plot's y-axis title far to the left.
# ------------------------------------------------------------------------------

panel_b1 <- ggplot(
  chel_ldl_curve,
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
    xintercept = ldl_median,
    linetype = "dotted",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_ribbon(
    aes(
      ymin = lb,
      ymax = ub
    ),
    fill = ldl_fill,
    alpha = 0.45
  ) +
  geom_line(
    color = ldl_col,
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1.06,
    label = p_chel_ldl_label,
    hjust = 1.05,
    vjust = 0.5,
    size = 3.4,
    parse = TRUE
  ) +
  scale_y_log10(
    limits = panel_b_y_limits
  ) +
  labs(
    x =
      "Therapy-corrected LDL-C, mmol/L",
    
    y =
      "HR for 1 s.d. higher ChEL pPRS",
    
    tag =
      "b"
  ) +
  panel_b_theme +
  theme(
    axis.title.y = element_text(
      margin = margin(
        r = 4
      )
    ),
    plot.tag = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),
    plot.tag.position = c(
      0,
      1
    ),
    plot.margin = margin(
      5.5,
      5.5,
      5.5,
      5.5
    )
  )

# ------------------------------------------------------------------------------
#  12. Panel B2: TRL pPRS by remnant cholesterol
#      RIGHT panel
# ------------------------------------------------------------------------------

panel_b2 <- ggplot(
  trl_rchol_curve,
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
    xintercept = rchol_median,
    linetype = "dotted",
    linewidth = 0.35,
    color = "black"
  ) +
  geom_ribbon(
    aes(
      ymin = lb,
      ymax = ub
    ),
    fill = rchol_fill,
    alpha = 0.45
  ) +
  geom_line(
    color = rchol_col,
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = Inf,
    y = 1.06,
    label = p_trl_rchol_label,
    hjust = 1.05,
    vjust = 0.5,
    size = 3.4,
    parse = TRUE
  ) +
  scale_y_log10(
    limits = panel_b_y_limits
  ) +
  labs(
    x =
      "Therapy-corrected remnant cholesterol, mmol/L",
    
    y =
      "HR for 1 s.d. higher TRL pPRS"
  ) +
  panel_b_theme +
  theme(
    axis.title.y = element_text(
      margin = margin(
        r = 2
      )
    )
  )

# ------------------------------------------------------------------------------
#  13. Combine Panel B
#
#      Important:
#      free(..., side = "l", type = "space") means that Panel B1 does not
#      inherit the large left-side spacing required by Panel A's long labels.
# ------------------------------------------------------------------------------

panel_b <- patchwork::free(
  panel_b1,
  side = "l",
  type = "space"
) +
  panel_b2 +
  plot_layout(
    ncol = 2
  )

# ------------------------------------------------------------------------------
#  14. Combine into one canvas
#
#      Manual panel tags only:
#      panel_a has tag "a"
#      panel_b1 has tag "b"
#      panel_b2 has no tag
#      Therefore, no panel c is generated.
# ------------------------------------------------------------------------------

fig_prs_interactions <-
  panel_a /
  panel_b +
  plot_layout(
    heights = c(
      1.25,
      1
    )
  )

fig_prs_interactions

# ------------------------------------------------------------------------------
#  15. Save one-canvas figure
# ------------------------------------------------------------------------------

ggsave(
  filename = "ext-data-fig-7.png",
  plot = fig_prs_interactions,
  width = 10,
  height = 7,
  dpi = 1200
)