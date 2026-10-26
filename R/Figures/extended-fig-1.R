library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# -----------------------------
# Set the work directory
# -----------------------------

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Desktop – Jakub’s Mac mini/Projects/UK Biobank/PPRS/bnmf-cad-clusters/bnmf-cad-clusters/figures/ext-data-figure-1")

# -----------------------------
# Load data
# -----------------------------

weights_df <- read_excel("weights_df.xlsx")

weight_cols <- paste0("X", 1:8)

global_cutoff <- 0.941758

cluster_cutoffs <- c(
  X1 = 0.4495459,
  X2 = 0.6491895,
  X3 = 0.3195119,
  X4 = 0.7870852,
  X5 = 0.6185486,
  X6 = 0.3837583,
  X7 = 0.2300184,
  X8 = 0.6937431
)

cluster_names <- c(
  X1 = "Coagulation",
  X2 = "MUNW",
  X3 = "TRL",
  X4 = "Endothelial",
  X5 = "Obesity",
  X6 = "ChEL",
  X7 = "Lp(a)",
  X8 = "Cardiac remodeling"
)

# Cardiac remodeling and Coagulation have swapped colours.
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

# Facet order follows the order of cluster_cols.
# Cardiac remodeling is therefore displayed last.
cluster_order <- names(cluster_cols)

# -----------------------------
# Long-format data
# -----------------------------

weights_long <- weights_df %>%
  pivot_longer(
    cols = all_of(weight_cols),
    names_to = "cluster_id",
    values_to = "weight"
  ) %>%
  mutate(
    cluster = cluster_names[cluster_id],
    cluster = factor(cluster, levels = cluster_order),
    cluster_cutoff = cluster_cutoffs[cluster_id]
  )

# -----------------------------
# Panel 1: global scree data
# -----------------------------

global_min <- min(weights_long$weight, na.rm = TRUE)
global_max <- max(weights_long$weight, na.rm = TRUE)

global_scree_df <- weights_long %>%
  mutate(
    weight_norm = rescale(
      weight,
      to = c(0, 1),
      from = c(global_min, global_max)
    )
  ) %>%
  arrange(desc(weight_norm)) %>%
  mutate(
    order_raw = row_number(),
    order_scaled = if (n() == 1) {
      0
    } else {
      (order_raw - 1) / (n() - 1)
    },
    point_group = if_else(
      weight >= global_cutoff,
      "selected",
      "not_selected"
    ),
    panel = factor(
      "Global cut-off",
      levels = "Global cut-off"
    )
  )

global_cutoff_x <- global_scree_df %>%
  filter(weight >= global_cutoff) %>%
  summarise(
    xintercept = if (n() == 0) {
      NA_real_
    } else {
      max(order_scaled, na.rm = TRUE)
    }
  ) %>%
  pull(xintercept)

global_annotation_df <- tibble(
  panel = factor(
    "Global cut-off",
    levels = "Global cut-off"
  ),
  x = 0.98,
  y = 1,
  label = paste0(
    "bold('Cut-off:')~",
    global_cutoff
  )
)

# -----------------------------
# Shared theme
# -----------------------------

scree_theme <- theme_linedraw(
  base_size = 11,
  base_family = "Arial"
) +
  theme(
    strip.background = element_rect(
      fill = "black",
      color = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(
      family = "Arial",
      face = "plain",
      color = "white",
      size = 11
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      family = "Arial"
    ),
    axis.text = element_text(
      family = "Arial"
    ),
    legend.position = "none",
    plot.title = element_blank()
  )

# -----------------------------
# Panel 1: global scree plot
# -----------------------------

p_global <- ggplot(
  global_scree_df,
  aes(
    x = order_scaled,
    y = weight_norm
  )
) +
  geom_point(
    aes(color = point_group),
    size = 1.4,
    alpha = 0.9
  ) +
  geom_vline(
    xintercept = global_cutoff_x,
    linetype = "dotted",
    color = "grey20",
    linewidth = 0.45
  ) +
  geom_text(
    data = global_annotation_df,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.15,
    size = 3.45,
    color = "grey20",
    family = "Arial",
    parse = TRUE
  ) +
  facet_grid(. ~ panel) +
  scale_color_manual(
    values = c(
      selected = "grey20",
      not_selected = "grey80"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 1),
    trans = "sqrt",
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 1),
    expand = expansion(
      mult = c(0.02, 0.12)
    )
  ) +
  labs(
    x = "Rank",
    y = "Weight"
  ) +
  scree_theme

# -----------------------------
# Panels 2-9: cluster-specific scree data
# -----------------------------

cluster_scree_df <- weights_long %>%
  group_by(cluster_id, cluster) %>%
  arrange(desc(weight), .by_group = TRUE) %>%
  mutate(
    order_raw = row_number(),
    
    order_scaled = if (n() == 1) {
      0
    } else {
      (order_raw - 1) / (n() - 1)
    },
    
    weight_scaled = if (
      max(weight, na.rm = TRUE) ==
      min(weight, na.rm = TRUE)
    ) {
      0
    } else {
      rescale(
        weight,
        to = c(0, 1),
        from = range(weight, na.rm = TRUE)
      )
    },
    
    cluster_cutoff = first(
      cluster_cutoffs[cluster_id]
    ),
    
    cluster_cutoff_scaled = if (
      max(weight, na.rm = TRUE) ==
      min(weight, na.rm = TRUE)
    ) {
      0
    } else {
      rescale(
        first(cluster_cutoffs[cluster_id]),
        to = c(0, 1),
        from = range(weight, na.rm = TRUE)
      )
    },
    
    global_cutoff_scaled = if (
      max(weight, na.rm = TRUE) ==
      min(weight, na.rm = TRUE)
    ) {
      0
    } else {
      rescale(
        global_cutoff,
        to = c(0, 1),
        from = range(weight, na.rm = TRUE)
      )
    },
    
    point_group = if_else(
      weight >= cluster_cutoff,
      as.character(cluster),
      "not_selected"
    )
  ) %>%
  ungroup() %>%
  mutate(
    cluster = factor(
      cluster,
      levels = cluster_order
    ),
    point_group = factor(
      point_group,
      levels = c(
        cluster_order,
        "not_selected"
      )
    )
  )

cutoff_df <- cluster_scree_df %>%
  group_by(cluster_id, cluster) %>%
  summarise(
    cluster_cutoff = first(cluster_cutoff),
    global_cutoff = global_cutoff,
    
    global_cutoff_x = if (
      sum(
        weight >= global_cutoff,
        na.rm = TRUE
      ) == 0
    ) {
      NA_real_
    } else {
      max(
        order_scaled[weight >= global_cutoff],
        na.rm = TRUE
      )
    },
    
    cluster_cutoff_x = if (
      sum(
        weight >= first(cluster_cutoff),
        na.rm = TRUE
      ) == 0
    ) {
      NA_real_
    } else {
      max(
        order_scaled[
          weight >= first(cluster_cutoff)
        ],
        na.rm = TRUE
      )
    },
    
    .groups = "drop"
  )

# Variant-count annotations:
# X = number of variants >= global cutoff within a given cluster
# Y = number of variants >= cluster-specific cutoff within a given cluster

annotation_df <- cluster_scree_df %>%
  group_by(cluster_id, cluster) %>%
  summarise(
    n_global = sum(
      weight >= global_cutoff,
      na.rm = TRUE
    ),
    
    n_cluster = sum(
      weight >= first(cluster_cutoff),
      na.rm = TRUE
    ),
    
    cluster_cutoff = first(cluster_cutoff),
    
    x = 0.98,
    y_variants = 1,
    y_cutoff = 1,
    
    label_variants = paste0(
      n_global,
      " %->% ",
      n_cluster,
      "~'variants'"
    ),
    
    label_cutoff = paste0(
      "bold('Cut-off:')~",
      signif(cluster_cutoff, 4)
    ),
    
    .groups = "drop"
  )

# -----------------------------
# Panels 2-9: cluster-specific scree plot
# -----------------------------

p_clusters <- ggplot(
  cluster_scree_df,
  aes(
    x = order_scaled,
    y = weight_scaled
  )
) +
  geom_point(
    aes(color = point_group),
    size = 1.4,
    alpha = 0.9
  ) +
  geom_vline(
    data = cutoff_df,
    aes(xintercept = global_cutoff_x),
    inherit.aes = FALSE,
    linetype = "dotted",
    color = "grey20",
    linewidth = 0.45
  ) +
  geom_vline(
    data = cutoff_df,
    aes(
      xintercept = cluster_cutoff_x,
      color = cluster
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.45
  ) +
  geom_text(
    data = annotation_df,
    aes(
      x = x,
      y = y_variants,
      label = label_variants
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1.15,
    size = 3.25,
    color = "grey20",
    family = "Arial",
    parse = TRUE
  ) +
  geom_text(
    data = annotation_df,
    aes(
      x = x,
      y = y_cutoff,
      label = label_cutoff
    ),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 2.65,
    size = 3.25,
    color = "grey20",
    family = "Arial",
    parse = TRUE
  ) +
  facet_wrap(
    ~ cluster,
    ncol = 4,
    scales = "free_y",
    drop = FALSE
  ) +
  scale_color_manual(
    values = c(
      cluster_cols,
      not_selected = "grey80"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 1),
    trans = "sqrt",
    expand = expansion(
      mult = c(0.04, 0.01)
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 1),
    expand = expansion(
      mult = c(0.03, 0.18)
    )
  ) +
  labs(
    x = "Rank",
    y = "Weight"
  ) +
  scree_theme

# -----------------------------
# Combine:
# Panel 1 left 1/3,
# Panels 2-9 right 2/3
# -----------------------------

final_plot <- p_global + p_clusters +
  plot_layout(
    widths = c(1, 2)
  )

final_plot

# -----------------------------
# Save outputs
# -----------------------------

ggsave(
  filename = "extended-data-fig-1.png",
  plot = final_plot,
  width = 12,
  height = 4,
  units = "in",
  dpi = 1200
)

# ============================================================
# Variant counts and median counts
# ============================================================

library(dplyr)
library(tidyr)
library(readxl)

weights_df <- read_excel("weights_df.xlsx")

weight_cols <- paste0("X", 1:8)

global_cutoff <- 0.941758

cluster_cutoffs <- c(
  X1 = 0.4495459,
  X2 = 0.6491895,
  X3 = 0.3195119,
  X4 = 0.7870852,
  X5 = 0.6185486,
  X6 = 0.3837583,
  X7 = 0.2300184,
  X8 = 0.6937431
)

cluster_names <- c(
  X1 = "Coagulation",
  X2 = "MUNW",
  X3 = "TRL",
  X4 = "Endothelial",
  X5 = "Obesity",
  X6 = "ChEL",
  X7 = "Lp(a)",
  X8 = "Cardiac remodeling"
)

cluster_order <- c(
  "ChEL",
  "TRL",
  "Lp(a)",
  "MUNW",
  "Obesity",
  "Endothelial",
  "Coagulation",
  "Cardiac remodeling"
)

variant_counts <- weights_df %>%
  pivot_longer(
    cols = all_of(weight_cols),
    names_to = "cluster_id",
    values_to = "weight"
  ) %>%
  mutate(
    cluster = cluster_names[cluster_id],
    cluster = factor(
      cluster,
      levels = cluster_order
    ),
    cluster_cutoff = cluster_cutoffs[cluster_id]
  ) %>%
  group_by(cluster_id, cluster) %>%
  summarise(
    n_global_cutoff = sum(
      weight >= global_cutoff,
      na.rm = TRUE
    ),
    n_cluster_specific_cutoff = sum(
      weight >= first(cluster_cutoff),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(cluster)

median_counts <- variant_counts %>%
  summarise(
    median_n_global_cutoff = median(
      n_global_cutoff,
      na.rm = TRUE
    ),
    median_n_cluster_specific_cutoff = median(
      n_cluster_specific_cutoff,
      na.rm = TRUE
    )
  )

variant_counts
median_counts