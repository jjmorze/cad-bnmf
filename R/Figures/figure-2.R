# ==============================================================================
# Integrated Shiny App
# Panel A: radial bNMF cluster plots
# Panel B: cell-type enrichment dot matrix
# ==============================================================================

options(shiny.maxRequestSize = 1024 * 1024^2)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(grid)
library(gridExtra)
library(DT)

# ==============================================================================
# CONSTANTS
# ==============================================================================

radial_clusters <- paste0("X", 1:8)

default_radial_cluster_order <- c(
  "X6", "X7", "X3", "X2",
  "X1", "X4", "X8", "X5"
)

default_dot_atlas_order <- c(
  "Coronary plaque",
  "Heart atlas",
  "Whole-body (CATLAS)"
)

default_dot_cluster_order <- c(
  "ChEL",
  "TRL",
  "Lp(a)",
  "MUNW",
  "Obesity",
  "Endothelial",
  "Coagulation",
  "Cardiac"
)

default_coronary_cells <- c(
  "Endothelial",
  "SMC",
  "Pericyte",
  "Fibroblast",
  "Macrophage",
  "Mast",
  "T",
  "Plasma"
)

default_heart_cells <- c(
  "aCM",
  "vCM",
  "FB",
  "EC",
  "SM",
  "MAC",
  "AD",
  "LC",
  "NR"
)

default_catlas_cells <- c(
  "hepatocyte",
  "melanocyte",
  "vascular_smooth_muscle_2",
  "vascular_smooth_muscle_1",
  "pericyte_general_3",
  "zona_fasciculata_cortical_cell",
  "endothelial_cell_myocardial",
  "airway_goblet_cell",
  "endothelial_cell_general_1",
  "endothelial_cell_general_2",
  "transitional_zone_cortical_cell",
  "zona_glomerulosa_cortical_cell",
  "luteal_cell_ovarian",
  "endocardial_cell",
  "cardiac_fibroblasts"
)

default_cell_label_overrides <- paste(
  "hepatocyte = Hepatocyte",
  "melanocyte = Melanocyte",
  "vascular_smooth_muscle_2 = VSMC 2 (artery tibial)",
  "vascular_smooth_muscle_1 = VSMC 1 (artery aorta)",
  "pericyte_general_3 = Pericyte 3",
  "zona_fasciculata_cortical_cell = Zona fasciculata",
  "endothelial_cell_myocardial = Endothelial (myocardial)",
  "airway_goblet_cell = Goblet cell (esophagus)",
  "endothelial_cell_general_1 = Endothelial 1",
  "endothelial_cell_general_2 = Endothelial 2",
  "transitional_zone_cortical_cell = Transitional zone cell",
  "zona_glomerulosa_cortical_cell = Zona glomerulosa",
  "luteal_cell_ovarian = Luteal cell",
  "endocardial_cell = Endocardial cell",
  "cardiac_fibroblasts = Cardiac fibroblast",
  sep = "\n"
)

# Saved figure settings embedded as the application default.
default_app_config <- list(
  "format_version" = 3,
  "generated" = "2026-07-26 12:31:18.025522",
  "general" = list(
    "variant_n_text" = "X1 = 7\nX2 = 41\nX3 = 33\nX4 = 13\nX5 = 32\nX6 = 34\nX7 = 26\nX8 = 50\n",
    "gene_labels_text" = "",
    "trait_labels_text" = "bmi = BMI\ngfatadjbmi = aGFAT\nhipcadjbmifemale = aHCf\nhipcadjbmimale = aHCm\nasatadjbmi = aASAT\nvatadjbmi = aVAT\nwhradjbmi = aWHR\nwbffm = WBFFM\nhdl = HDL-C\nldl = LDL-C\nlpa = Lp(a)\nlogtg = TG\nace = ACE\nace2 = ACE2\nagt = AGT\ndbp = DBP\nren = Ren\nsbp = SBP\nbnp = BNP\nepat = EPAT\ngal3 = Gal-3\ngdf15 = GDF-15\nhr = HR\nfabp3 = FABP3\ntnt = TnT\nlvmbsa = LVMBSA\nlvef = LVEF\nlvedvbsa = aLVEDV\nlvesvbsa = aLVESV\nntprobnp = NTproBNP\nprinterval = PR\nqrsduration = QRSd\nqtcinterval = QTc\nil1r1 = sST2\naptt = aPTT\nddimer = D-dimer\nf7 = FVII\nf8 = FVIII\nmpv = MPV\npai1 = PAI-1\nplt = PLT\npdw = PDW\npf4 = PF4\nprothrombintime = PT\ntpa = tPA\nvwf = vWF\narterstiff = AS\nsele = SELE\nicam1 = ICAM-1\nselp = SELP\nsvbsa = aSV\nvcam1 = VCAM-1\n2hrgluadjbmi = a2hrGlu\nadipoq = AdipoQ\nfiadjbmi = aFI\nhba1cadjbmi = aHbA1c\nheight = Height\nifcadjbmi = aIFC\nisiadjbmi = aISI\nleptin = Leptin\nretn = Retn\ngpx1 = GPx1\nhscrp = hsCRP\nil1b = IL-1b\nil6 = IL-6\nlympho = Lympho\nmcp1 = MCP-1\nmpo = MPO\nneutro = Neutro\nrbdw = RBDW\ntnf = TNFa\ncreatinine = Creat\ncystatinc = Cystatin C\nmicroalburi = Microalb\nkuri = Urine K\nnauri = Urine Na\nsod1 = SOD-1\nsod3 = SOD-3\nurea = Urea\nuacr = UACR\nalt = ALT\nalbumin = Albumin\nalp = ALP\nast = AST\nggt = GGT\nliverfat = PDFF\ntotalbilirubin = Bilirubin",
    "top_panel_height" = 2.1,
    "bottom_panel_height" = 1.2,
    "panel_a_label" = "a",
    "panel_a_label_x" = 0.01,
    "panel_a_label_y" = 0.94,
    "panel_a_label_size" = 24,
    "panel_a_label_colour" = "black",
    "panel_b_label" = "b",
    "panel_b_label_x" = 0.01,
    "panel_b_label_y" = 1.075,
    "panel_b_label_size" = 24,
    "panel_b_label_colour" = "black",
    "preview_height_px" = 1100,
    "cluster_order" = "X6,X3,X7,X2,X5,X4,X1,X8",
    "panel_ncol" = 4,
    "panel_nrow" = 2,
    "panel_h_spacing" = -0.1,
    "panel_v_spacing" = -0.15,
    "plot_scale" = 1.3,
    "panel_background" = "transparent",
    "show_titles" = TRUE,
    "title_size" = 15.1,
    "title_wrap_lines" = 1,
    "title_fontface" = "bold",
    "title_color" = "black",
    "title_mapping_text" = "X1 = Coagulation\nX2 = Metabolically unhealthy normal weight\nX3 = Triglyceride-rich lipoproteins\nX4 = Endothelial dysfunction\nX5 = Obesity\nX6 = Cholesterol-enriched lipoproteins\nX7 = Lipoprotein(a)\nX8 = Cardiac remodeling",
    "title_x_default" = 0.5,
    "title_y_default" = 0.83,
    "title_x_mapping_text" = "",
    "title_y_mapping_text" = "",
    "propagate_source_cluster" = "X1",
    "propagate_setting_group" = "All visual and filtering settings",
    "dot_atlas_order" = "Coronary plaque,Heart atlas,Whole-body (CATLAS)",
    "dot_atlas_labels" = "Coronary plaque = Coronary Plaque (Turner 2022)\nHeart atlas = Heart Atlas (CATlas)\nWhole-body (CATLAS) = Human Multi-tissues (CATlas)",
    "dot_atlas_gap" = 1.2,
    "dot_cluster_order" = "ChEL,TRL,Lp(a),MUNW,Obesity,Endothelial,Coagulation,Cardiac",
    "dot_cluster_labels" = "ChEL = ChEL\nTRL = TRL\nLp(a) = Lp(a)\nMUNW = MUNW\nObesity = Obesity\nEndothelial = Endothelial dysfunction\nCoagulation = Coagulation\nCardiac = Cardiac remodeling",
    "dot_coronary_cells" = paste(default_coronary_cells, collapse = "\n"),
    "dot_heart_cells" = paste(default_heart_cells, collapse = "\n"),
    "dot_catlas_cells" = paste(default_catlas_cells, collapse = "\n"),
    "dot_cell_labels" = default_cell_label_overrides,
    "dot_use_cell_atlas_name" = TRUE,
    "dot_nominal_alpha" = 0.05,
    "dot_bonferroni_n_tests" = 128,
    "dot_size_nonsignificant" = 2.8,
    "dot_size_nominal" = 5,
    "dot_size_bonferroni" = 7.5,
    "dot_outline_nonsignificant" = 0.5,
    "dot_outline_significant" = 1,
    "dot_outline_nonsignificant_colour" = "grey25",
    "dot_outline_significant_colour" = "grey25",
    "dot_significance_legend_colour" = "grey60",
    "dot_z_min" = -5,
    "dot_z_max" = 5,
    "dot_z_breaks" = "-5,-2.5,0,2.5,5",
    "dot_colour_low" = "#4682B4",
    "dot_colour_mid" = "white",
    "dot_colour_high" = "#CB4154",
    "dot_z_legend_title" = "Z-score",
    "dot_significance_legend_title" = "Significance",
    "dot_tile_fill" = "grey93",
    "dot_tile_border" = "white",
    "dot_tile_border_width" = 0.65,
    "dot_tile_width" = 0.92,
    "dot_tile_height" = 0.92,
    "dot_x_text_angle" = 45,
    "dot_x_text_size" = 12,
    "dot_x_text_face" = "plain",
    "dot_x_text_hjust" = 1,
    "dot_x_text_vjust" = 1,
    "dot_x_text_colour" = "grey20",
    "dot_y_axis_title" = "",
    "dot_y_axis_title_size" = 11,
    "dot_y_axis_title_face" = "bold",
    "dot_y_axis_title_colour" = "grey20",
    "dot_y_axis_title_margin" = 10,
    "dot_y_text_size" = 12,
    "dot_y_text_face" = "plain",
    "dot_y_text_colour" = "grey20",
    "dot_show_atlas_headers" = TRUE,
    "dot_show_atlas_brackets" = TRUE,
    "dot_atlas_header_size" = 15,
    "dot_atlas_header_face" = "bold",
    "dot_atlas_header_colour" = "grey20",
    "dot_atlas_bracket_offset" = 0.75,
    "dot_atlas_label_offset" = 1.2,
    "dot_atlas_bracket_colour" = "grey45",
    "dot_atlas_bracket_width" = 0.5,
    "dot_show_atlas_separator" = TRUE,
    "dot_atlas_gap_line_width" = 4,
    "dot_atlas_separator_width" = 0.45,
    "dot_atlas_separator_colour" = "grey65",
    "dot_atlas_separator_linetype" = "dashed",
    "dot_base_text_size" = 11,
    "dot_legend_position" = "right",
    "dot_legend_title_size" = 12,
    "dot_legend_text_size" = 11,
    "dot_legend_title_face" = "bold",
    "dot_legend_key_height" = 0.7,
    "dot_colorbar_width" = 1,
    "dot_colorbar_height" = 4,
    "dot_margin_top" = -60,
    "dot_margin_right" = 15,
    "dot_margin_bottom" = 45,
    "dot_margin_left" = 15,
    "dot_x_expand_left" = 0.01,
    "dot_x_expand_right" = 0.01,
    "dot_y_expand_bottom" = 0.04,
    "dot_y_expand_top" = 0.16,
    "export_width" = 17,
    "export_height" = 15,
    "export_dpi" = 1200
  ),
  "radial_clusters" = list(
    "X1" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "manual",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "f8 = 0.16",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X2" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "global",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "hscrp = 0.18\nldl = 0.18\nhdl = 0.18\nbmi = 0.18\nsbp = 0.24\ndbp = 0.24",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X3" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "manual",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "logtg = 0.14",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X4" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "manual",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "vcam1 = 0.3\ncystatinc = 0.33\nlvmbsa = 0.3",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X5" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "manual",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "cystatinc = 0.33\nlvmbsa = 0.3",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X6" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "global",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "hscrp = 0.18\nldl = 0.18\nhdl = 0.18\nbmi = 0.18\nsbp = 0.24\ndbp = 0.24",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X7" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "manual",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "lpa = 0.18\nlvmbsa = 0.28",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    ),
    "X8" = list(
      "n_top_genes" = 10,
      "gene_weight_threshold" = 0,
      "n_top_traits" = 10,
      "trait_weight_threshold" = 0,
      "gene_radius" = 0.5,
      "inner_radius_offset" = 0,
      "bar_scale" = 0.8,
      "bar_width" = 0.23,
      "pos_range" = list(16, 170),
      "neg_range" = list(190, 346),
      "axis_intervals" = 3,
      "show_axis" = TRUE,
      "axis_label_size" = 5.5,
      "outer_radius_extend" = 0,
      "show_outer_circle" = TRUE,
      "outer_circle_color" = "grey80",
      "outer_circle_linewidth" = 0.25,
      "variant_circle_scale" = 1,
      "variant_circle_alpha" = 1,
      "variant_circle_transform" = "sqrt",
      "variant_label_size" = 7,
      "variant_label_x" = 0,
      "variant_label_y" = 0,
      "variant_label_wrap_width" = 12,
      "variant_label_color" = "black",
      "variant_label_fontface" = "bold",
      "variant_label_two_lines" = TRUE,
      "trait_label_size" = 6,
      "trait_label_spacing_mode" = "manual",
      "trait_label_offset_short" = 0.18,
      "trait_label_offset_medium" = 0.24,
      "trait_label_offset_long" = 0.28,
      "trait_label_manual_offsets" = "ntprobnp = 0.32\nlvmbsa = 0.3\nlvesvbsa = 0.3 ",
      "trait_label_angle_adjust" = 0,
      "trait_labels_radial" = TRUE,
      "color_gene" = "#3EB489",
      "color_positive" = "#CB4154",
      "color_negative" = "#4682B4",
      "color_axis" = "gray70",
      "color_trait_label" = "gray20",
      "x_offset" = 0,
      "y_offset" = 0
    )
  )
)


# ==============================================================================
# RADIAL SETTING GROUPS
# ==============================================================================

radial_setting_groups <- list(
  "All visual and filtering settings" = c(
    "n_top_genes",
    "gene_weight_threshold",
    "n_top_traits",
    "trait_weight_threshold",
    "gene_radius",
    "inner_radius_offset",
    "bar_scale",
    "bar_width",
    "pos_range",
    "neg_range",
    "axis_intervals",
    "show_axis",
    "axis_label_size",
    "outer_radius_extend",
    "show_outer_circle",
    "outer_circle_color",
    "outer_circle_linewidth",
    "variant_circle_scale",
    "variant_circle_alpha",
    "variant_circle_transform",
    "variant_label_size",
    "variant_label_x",
    "variant_label_y",
    "variant_label_wrap_width",
    "variant_label_color",
    "variant_label_fontface",
    "variant_label_two_lines",
    "trait_label_size",
    "trait_label_spacing_mode",
    "trait_label_offset_short",
    "trait_label_offset_medium",
    "trait_label_offset_long",
    "trait_label_manual_offsets",
    "trait_label_angle_adjust",
    "trait_labels_radial",
    "color_gene",
    "color_positive",
    "color_negative",
    "color_axis",
    "color_trait_label",
    "x_offset",
    "y_offset"
  ),
  
  "Thresholds and top-N" = c(
    "n_top_genes",
    "gene_weight_threshold",
    "n_top_traits",
    "trait_weight_threshold"
  ),
  
  "Geometry and axes" = c(
    "gene_radius",
    "inner_radius_offset",
    "bar_scale",
    "bar_width",
    "pos_range",
    "neg_range",
    "axis_intervals",
    "show_axis",
    "axis_label_size",
    "outer_radius_extend",
    "show_outer_circle",
    "outer_circle_color",
    "outer_circle_linewidth"
  ),
  
  "Central circle and label" = c(
    "variant_circle_scale",
    "variant_circle_alpha",
    "variant_circle_transform",
    "variant_label_size",
    "variant_label_x",
    "variant_label_y",
    "variant_label_wrap_width",
    "variant_label_color",
    "variant_label_fontface",
    "variant_label_two_lines"
  ),
  
  "Trait labels" = c(
    "trait_label_size",
    "trait_label_spacing_mode",
    "trait_label_offset_short",
    "trait_label_offset_medium",
    "trait_label_offset_long",
    "trait_label_manual_offsets",
    "trait_label_angle_adjust",
    "trait_labels_radial"
  ),
  
  "Colours" = c(
    "outer_circle_color",
    "variant_label_color",
    "color_gene",
    "color_positive",
    "color_negative",
    "color_axis",
    "color_trait_label"
  ),
  
  "Plot position" = c(
    "x_offset",
    "y_offset"
  )
)

radial_numeric_suffixes <- c(
  "n_top_genes",
  "gene_weight_threshold",
  "n_top_traits",
  "trait_weight_threshold",
  "gene_radius",
  "inner_radius_offset",
  "bar_scale",
  "bar_width",
  "axis_intervals",
  "axis_label_size",
  "outer_radius_extend",
  "outer_circle_linewidth",
  "variant_circle_scale",
  "variant_circle_alpha",
  "variant_label_size",
  "variant_label_x",
  "variant_label_y",
  "variant_label_wrap_width",
  "trait_label_size",
  "trait_label_offset_short",
  "trait_label_offset_medium",
  "trait_label_offset_long",
  "trait_label_angle_adjust",
  "x_offset",
  "y_offset"
)

radial_checkbox_suffixes <- c(
  "show_axis",
  "show_outer_circle",
  "variant_label_two_lines",
  "trait_labels_radial"
)

radial_select_suffixes <- c(
  "variant_circle_transform",
  "variant_label_fontface",
  "trait_label_spacing_mode"
)

radial_text_suffixes <- c(
  "outer_circle_color",
  "variant_label_color",
  "trait_label_manual_offsets",
  "color_gene",
  "color_positive",
  "color_negative",
  "color_axis",
  "color_trait_label"
)

radial_slider_suffixes <- c(
  "pos_range",
  "neg_range"
)

all_radial_setting_suffixes <- unique(
  unlist(
    radial_setting_groups,
    use.names = FALSE
  )
)

# ==============================================================================
# EXAMPLE DATA
# ==============================================================================

generate_example_circular_data <- function() {
  set.seed(123)
  
  n_genes <- 15
  
  gene_weight_df <- data.frame(
    VAR_ID = paste0("var_", seq_len(n_genes)),
    rsID = paste0("rs", sample(1000000:9999999, n_genes)),
    gene = paste0("Gene_", LETTERS[seq_len(n_genes)]),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(8)) {
    gene_weight_df[[paste0("X", i)]] <- runif(n_genes, 0.3, 1)
  }
  
  base_traits <- c(
    "hscrp",
    "ldl",
    "hdl",
    "bmi",
    "sbp",
    "dbp",
    "glucose",
    "hba1c",
    "trig",
    "chol"
  )
  
  trait_weight_df <- data.frame(
    trait = c(
      paste0(base_traits, "_pos"),
      paste0(base_traits, "_neg")
    ),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(8)) {
    pos_weights <- runif(length(base_traits), 0.3, 1)
    neg_weights <- runif(length(base_traits), -1, -0.3)
    trait_weight_df[[paste0("X", i)]] <- c(pos_weights, neg_weights)
  }
  
  list(
    gene = gene_weight_df,
    trait = trait_weight_df
  )
}

generate_example_dot_data <- function() {
  set.seed(321)
  
  atlas_cells <- list(
    "Coronary plaque" = default_coronary_cells,
    "Heart atlas" = default_heart_cells,
    "Whole-body (CATLAS)" = default_catlas_cells
  )
  
  example_df <- bind_rows(
    lapply(names(atlas_cells), function(atlas_name) {
      expand.grid(
        atlas = atlas_name,
        cell = atlas_cells[[atlas_name]],
        cluster = default_dot_cluster_order,
        stringsAsFactors = FALSE
      )
    })
  )
  
  example_df %>%
    mutate(
      group = atlas,
      cell_atlas_name = cell,
      z = rnorm(n(), mean = 0.6, sd = 1.8),
      p_nominal = pmin(1, 2 * pnorm(-abs(z)))
    )
}

# ==============================================================================
# HELPERS
# ==============================================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) {
    y
  } else {
    x
  }
}

safe_read_csv <- function(path) {
  extension <- tolower(tools::file_ext(path))
  
  if (extension == "gz") {
    read.csv(
      gzfile(path),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}

parse_named_mapping <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(list())
  }
  
  lines <- trimws(unlist(strsplit(text, "\n")))
  lines <- lines[nzchar(lines)]
  
  output <- list()
  
  for (line in lines) {
    parts <- unlist(strsplit(line, "=", fixed = TRUE))
    
    if (length(parts) >= 2) {
      key <- trimws(parts[1])
      value <- trimws(paste(parts[-1], collapse = "="))
      
      if (nzchar(key)) {
        output[[key]] <- value
      }
    }
  }
  
  output
}

parse_numeric_mapping <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(list())
  }
  
  lines <- trimws(unlist(strsplit(text, "\n")))
  lines <- lines[nzchar(lines)]
  
  output <- list()
  
  for (line in lines) {
    parts <- unlist(strsplit(line, "=", fixed = TRUE))
    
    if (length(parts) >= 2) {
      key <- trimws(parts[1])
      
      value <- suppressWarnings(
        as.numeric(
          trimws(
            paste(
              parts[-1],
              collapse = "="
            )
          )
        )
      )
      
      if (nzchar(key) && is.finite(value)) {
        output[[key]] <- value
      }
    }
  }
  
  output
}

parse_comma_vector <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(character(0))
  }
  
  values <- trimws(unlist(strsplit(text, ",", fixed = TRUE)))
  values[nzchar(values)]
}

parse_line_vector <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(character(0))
  }
  
  values <- trimws(unlist(strsplit(text, "\n")))
  values[nzchar(values)]
}

safe_colour <- function(value, fallback = "black") {
  if (
    is.null(value) ||
    length(value) == 0 ||
    is.na(value) ||
    !nzchar(value)
  ) {
    return(fallback)
  }
  
  test <- try(
    grDevices::col2rgb(value),
    silent = TRUE
  )
  
  if (inherits(test, "try-error")) {
    fallback
  } else {
    value
  }
}

wrap_text_width <- function(text, width) {
  width <- suppressWarnings(as.numeric(width))
  
  if (
    is.null(text) ||
    length(text) == 0 ||
    is.na(text) ||
    !nzchar(text)
  ) {
    return(text)
  }
  
  if (!is.finite(width) || width <= 0) {
    return(text)
  }
  
  stringr::str_wrap(
    text,
    width = width
  )
}

# ==============================================================================
# RADIAL DEFAULTS
# ==============================================================================

cluster_defaults <- function() {
  list(
    n_top_genes = 10,
    gene_weight_threshold = 0,
    
    n_top_traits = 10,
    trait_weight_threshold = 0,
    
    gene_radius = 0.6,
    inner_radius_offset = 0.1,
    bar_scale = 0.8,
    bar_width = 0.18,
    
    pos_range = c(10, 170),
    neg_range = c(190, 350),
    
    axis_intervals = 3,
    show_axis = TRUE,
    axis_label_size = 2.5,
    outer_radius_extend = 0,
    
    show_outer_circle = TRUE,
    outer_circle_color = "grey60",
    outer_circle_linewidth = 0.5,
    
    variant_circle_scale = 1,
    variant_circle_alpha = 1,
    variant_circle_transform = "linear",
    
    variant_label_size = 3.5,
    variant_label_x = 0,
    variant_label_y = 0,
    variant_label_wrap_width = 12,
    variant_label_color = "white",
    variant_label_fontface = "bold",
    variant_label_two_lines = TRUE,
    
    trait_label_size = 4,
    trait_label_spacing_mode = "global",
    trait_label_offset_short = 0.18,
    trait_label_offset_medium = 0.24,
    trait_label_offset_long = 0.28,
    trait_label_manual_offsets = "",
    trait_label_angle_adjust = 0,
    trait_labels_radial = TRUE,
    
    color_gene = "#3EB489",
    color_positive = "#CB4154",
    color_negative = "#4682B4",
    color_axis = "gray70",
    color_trait_label = "gray20",
    
    x_offset = 0,
    y_offset = 0
  )
}

get_cluster_settings_from_input <- function(input, cluster) {
  defaults <- cluster_defaults()
  
  list(
    n_top_genes = input[[paste0(cluster, "_n_top_genes")]] %||% defaults$n_top_genes,
    gene_weight_threshold = input[[paste0(cluster, "_gene_weight_threshold")]] %||% defaults$gene_weight_threshold,
    n_top_traits = input[[paste0(cluster, "_n_top_traits")]] %||% defaults$n_top_traits,
    trait_weight_threshold = input[[paste0(cluster, "_trait_weight_threshold")]] %||% defaults$trait_weight_threshold,
    
    gene_radius = input[[paste0(cluster, "_gene_radius")]] %||% defaults$gene_radius,
    inner_radius_offset = input[[paste0(cluster, "_inner_radius_offset")]] %||% defaults$inner_radius_offset,
    bar_scale = input[[paste0(cluster, "_bar_scale")]] %||% defaults$bar_scale,
    bar_width = input[[paste0(cluster, "_bar_width")]] %||% defaults$bar_width,
    
    pos_range = input[[paste0(cluster, "_pos_range")]] %||% defaults$pos_range,
    neg_range = input[[paste0(cluster, "_neg_range")]] %||% defaults$neg_range,
    
    axis_intervals = input[[paste0(cluster, "_axis_intervals")]] %||% defaults$axis_intervals,
    show_axis = input[[paste0(cluster, "_show_axis")]] %||% defaults$show_axis,
    axis_label_size = input[[paste0(cluster, "_axis_label_size")]] %||% defaults$axis_label_size,
    outer_radius_extend = input[[paste0(cluster, "_outer_radius_extend")]] %||% defaults$outer_radius_extend,
    
    show_outer_circle = input[[paste0(cluster, "_show_outer_circle")]] %||% defaults$show_outer_circle,
    outer_circle_color = input[[paste0(cluster, "_outer_circle_color")]] %||% defaults$outer_circle_color,
    outer_circle_linewidth = input[[paste0(cluster, "_outer_circle_linewidth")]] %||% defaults$outer_circle_linewidth,
    
    variant_circle_scale = input[[paste0(cluster, "_variant_circle_scale")]] %||% defaults$variant_circle_scale,
    variant_circle_alpha = input[[paste0(cluster, "_variant_circle_alpha")]] %||% defaults$variant_circle_alpha,
    variant_circle_transform = input[[paste0(cluster, "_variant_circle_transform")]] %||% defaults$variant_circle_transform,
    
    variant_label_size = input[[paste0(cluster, "_variant_label_size")]] %||% defaults$variant_label_size,
    variant_label_x = input[[paste0(cluster, "_variant_label_x")]] %||% defaults$variant_label_x,
    variant_label_y = input[[paste0(cluster, "_variant_label_y")]] %||% defaults$variant_label_y,
    variant_label_wrap_width = input[[paste0(cluster, "_variant_label_wrap_width")]] %||% defaults$variant_label_wrap_width,
    variant_label_color = input[[paste0(cluster, "_variant_label_color")]] %||% defaults$variant_label_color,
    variant_label_fontface = input[[paste0(cluster, "_variant_label_fontface")]] %||% defaults$variant_label_fontface,
    variant_label_two_lines = input[[paste0(cluster, "_variant_label_two_lines")]] %||% defaults$variant_label_two_lines,
    
    trait_label_size = input[[paste0(cluster, "_trait_label_size")]] %||% defaults$trait_label_size,
    trait_label_spacing_mode = input[[paste0(cluster, "_trait_label_spacing_mode")]] %||% defaults$trait_label_spacing_mode,
    trait_label_offset_short = input[[paste0(cluster, "_trait_label_offset_short")]] %||% defaults$trait_label_offset_short,
    trait_label_offset_medium = input[[paste0(cluster, "_trait_label_offset_medium")]] %||% defaults$trait_label_offset_medium,
    trait_label_offset_long = input[[paste0(cluster, "_trait_label_offset_long")]] %||% defaults$trait_label_offset_long,
    trait_label_manual_offsets = input[[paste0(cluster, "_trait_label_manual_offsets")]] %||% defaults$trait_label_manual_offsets,
    trait_label_angle_adjust = input[[paste0(cluster, "_trait_label_angle_adjust")]] %||% defaults$trait_label_angle_adjust,
    trait_labels_radial = input[[paste0(cluster, "_trait_labels_radial")]] %||% defaults$trait_labels_radial,
    
    color_gene = input[[paste0(cluster, "_color_gene")]] %||% defaults$color_gene,
    color_positive = input[[paste0(cluster, "_color_positive")]] %||% defaults$color_positive,
    color_negative = input[[paste0(cluster, "_color_negative")]] %||% defaults$color_negative,
    color_axis = input[[paste0(cluster, "_color_axis")]] %||% defaults$color_axis,
    color_trait_label = input[[paste0(cluster, "_color_trait_label")]] %||% defaults$color_trait_label,
    
    x_offset = input[[paste0(cluster, "_x_offset")]] %||% defaults$x_offset,
    y_offset = input[[paste0(cluster, "_y_offset")]] %||% defaults$y_offset
  )
}

# ==============================================================================
# RADIAL PLOT
# ==============================================================================

create_cluster_plot <- function(
    cluster,
    gene_df,
    trait_df,
    settings,
    trait_labels,
    gene_labels,
    variant_n_map = list(),
    for_panel = TRUE
) {
  if (is.null(gene_df) || is.null(trait_df)) {
    return(NULL)
  }
  
  if (!(cluster %in% names(gene_df)) || !(cluster %in% names(trait_df))) {
    return(NULL)
  }
  
  if (!("gene" %in% names(gene_df)) || !("trait" %in% names(trait_df))) {
    return(NULL)
  }
  
  s <- settings
  
  genes <- gene_df %>%
    select(
      gene,
      weight = all_of(cluster)
    ) %>%
    mutate(
      weight = as.numeric(weight)
    ) %>%
    filter(!is.na(weight)) %>%
    filter(abs(weight) >= s$gene_weight_threshold) %>%
    arrange(desc(abs(weight))) %>%
    head(s$n_top_genes)
  
  if (nrow(genes) == 0) {
    return(NULL)
  }
  
  n_variants_observed <- nrow(genes)
  n_variants_for_circle <- variant_n_map[[cluster]] %||% n_variants_observed
  
  manual_offset_map <- parse_numeric_mapping(
    s$trait_label_manual_offsets
  )
  
  traits_all <- trait_df %>%
    select(
      trait,
      weight = all_of(cluster)
    ) %>%
    mutate(
      weight = as.numeric(weight)
    ) %>%
    filter(!is.na(weight)) %>%
    mutate(
      type = ifelse(
        grepl("_pos$", trait),
        "positive",
        "negative"
      ),
      
      base_trait = gsub(
        "_pos$|_neg$",
        "",
        trait
      ),
      
      display_label = vapply(
        base_trait,
        function(x) {
          trait_labels[[x]] %||% x
        },
        character(1)
      ),
      
      original_weight = abs(weight)
    )
  
  pos_traits <- traits_all %>%
    filter(type == "positive") %>%
    filter(abs(weight) >= s$trait_weight_threshold) %>%
    arrange(desc(abs(weight))) %>%
    head(s$n_top_traits)
  
  neg_traits <- traits_all %>%
    filter(type == "negative") %>%
    filter(abs(weight) >= s$trait_weight_threshold) %>%
    arrange(desc(abs(weight))) %>%
    head(s$n_top_traits)
  
  if (nrow(pos_traits) == 0 && nrow(neg_traits) == 0) {
    return(NULL)
  }
  
  traits_combined <- bind_rows(
    pos_traits,
    neg_traits
  )
  
  original_min <- min(
    traits_combined$original_weight
  )
  
  original_max <- max(
    traits_combined$original_weight
  )
  
  if (!is.finite(original_min) || !is.finite(original_max)) {
    return(NULL)
  }
  
  if (original_max == original_min) {
    original_max <- original_min + 1
  }
  
  gene_radius <- s$gene_radius
  
  inner_radius <- gene_radius +
    s$inner_radius_offset
  
  effective_bar_scale <- s$bar_scale +
    s$outer_radius_extend
  
  if (nrow(pos_traits) > 0) {
    pos_traits$angle <- seq(
      s$pos_range[1],
      s$pos_range[2],
      length.out = nrow(pos_traits)
    )
    
    pos_traits$angle_rad <- pos_traits$angle * pi / 180
    
    pos_traits$bar_weight <- (
      pos_traits$original_weight -
        original_min
    ) / (
      original_max -
        original_min
    )
  }
  
  if (nrow(neg_traits) > 0) {
    neg_traits$angle <- seq(
      s$neg_range[1],
      s$neg_range[2],
      length.out = nrow(neg_traits)
    )
    
    neg_traits$angle_rad <- neg_traits$angle * pi / 180
    
    neg_traits$bar_weight <- (
      neg_traits$original_weight -
        original_min
    ) / (
      original_max -
        original_min
    )
  }
  
  traits <- bind_rows(
    pos_traits,
    neg_traits
  )
  
  get_global_trait_offset <- function(label) {
    if (nchar(label) <= 3) {
      s$trait_label_offset_short
    } else if (nchar(label) <= 6) {
      s$trait_label_offset_medium
    } else {
      s$trait_label_offset_long
    }
  }
  
  get_manual_trait_offset <- function(
    trait_name,
    base_trait,
    display_label
  ) {
    possible_keys <- c(
      trait_name,
      base_trait,
      display_label
    )
    
    matched_keys <- possible_keys[
      possible_keys %in% names(manual_offset_map)
    ]
    
    if (length(matched_keys) > 0) {
      as.numeric(
        manual_offset_map[[matched_keys[1]]]
      )
    } else {
      get_global_trait_offset(
        display_label
      )
    }
  }
  
  traits$label_offset <- vapply(
    seq_len(nrow(traits)),
    function(i) {
      if (identical(s$trait_label_spacing_mode, "manual")) {
        get_manual_trait_offset(
          trait_name = traits$trait[i],
          base_trait = traits$base_trait[i],
          display_label = traits$display_label[i]
        )
      } else {
        get_global_trait_offset(
          traits$display_label[i]
        )
      }
    },
    numeric(1)
  )
  
  traits <- traits %>%
    mutate(
      r_start = inner_radius,
      
      r_end = inner_radius +
        bar_weight *
        effective_bar_scale,
      
      label_radius = r_end +
        label_offset,
      
      x_label = cos(angle_rad) *
        label_radius,
      
      y_label = sin(angle_rad) *
        label_radius,
      
      label_angle = if (isTRUE(s$trait_labels_radial)) {
        ifelse(
          angle > 90 &
            angle < 270,
          angle + 180,
          angle
        ) +
          s$trait_label_angle_adjust
      } else {
        s$trait_label_angle_adjust
      }
    )
  
  create_bar_polygon <- function(
    angle_rad,
    r_start,
    r_end,
    width
  ) {
    half_width <- width / 2
    
    angles <- c(
      angle_rad - half_width,
      angle_rad + half_width,
      angle_rad + half_width,
      angle_rad - half_width
    )
    
    radii <- c(
      r_start,
      r_start,
      r_end,
      r_end
    )
    
    data.frame(
      x = cos(angles) * radii,
      y = sin(angles) * radii
    )
  }
  
  make_circle_df <- function(
    radius,
    n = 300
  ) {
    angles <- seq(
      0,
      2 * pi,
      length.out = n
    )
    
    data.frame(
      x = cos(angles) * radius,
      y = sin(angles) * radius
    )
  }
  
  bar_polygons <- do.call(
    rbind,
    lapply(seq_len(nrow(traits)), function(i) {
      polygon <- create_bar_polygon(
        angle_rad = traits$angle_rad[i],
        r_start = traits$r_start[i],
        r_end = traits$r_end[i],
        width = s$bar_width
      )
      
      polygon$id <- i
      polygon$type <- traits$type[i]
      
      polygon
    })
  )
  
  outer_circle_radius <- max(
    traits$r_end
  )
  
  if (isTRUE(s$show_axis)) {
    axis_intervals <- max(
      2,
      round(s$axis_intervals)
    )
    
    scaled_values <- seq(
      0,
      1,
      length.out = axis_intervals
    )
    
    original_values <- original_min +
      scaled_values *
      (
        original_max -
          original_min
      )
    
    axis_ticks <- data.frame(
      scaled_value = scaled_values,
      original_value = original_values
    ) %>%
      mutate(
        radius = inner_radius +
          scaled_value *
          effective_bar_scale
      )
    
    axis_circles <- do.call(
      rbind,
      lapply(axis_ticks$radius, function(radius_value) {
        circle <- make_circle_df(
          radius_value,
          n = 100
        )
        
        circle$radius <- radius_value
        circle
      })
    )
    
    axis_labels <- axis_ticks %>%
      mutate(
        x = radius + 0.05,
        y = 0,
        label = sprintf(
          "%.2f",
          original_value
        )
      )
  }
  
  maximum_label_offset <- max(
    traits$label_offset,
    na.rm = TRUE
  )
  
  if (!is.finite(maximum_label_offset)) {
    maximum_label_offset <- s$trait_label_offset_long
  }
  
  max_radius <- inner_radius +
    effective_bar_scale +
    maximum_label_offset +
    0.4
  
  size_factor <- if (for_panel) {
    0.7
  } else {
    1
  }
  
  make_variant_label <- function(
    n,
    two_lines = TRUE,
    wrap_width = 12
  ) {
    variant_word <- ifelse(
      n == 1,
      "variant",
      "variants"
    )
    
    if (isTRUE(two_lines)) {
      paste0(
        n,
        "\n",
        variant_word
      )
    } else {
      str_wrap(
        paste0(
          n,
          " ",
          variant_word
        ),
        width = wrap_width
      )
    }
  }
  
  variant_circle_df <- data.frame(
    x = s$variant_label_x,
    y = s$variant_label_y,
    
    label = make_variant_label(
      n = n_variants_for_circle,
      two_lines = s$variant_label_two_lines,
      wrap_width = s$variant_label_wrap_width
    ),
    
    stringsAsFactors = FALSE
  )
  
  variant_n_values <- suppressWarnings(
    as.numeric(
      unlist(
        variant_n_map,
        use.names = FALSE
      )
    )
  )
  
  variant_n_values <- variant_n_values[
    is.finite(variant_n_values) &
      variant_n_values > 0
  ]
  
  max_variants_for_scaling <- if (
    length(variant_n_values) > 0
  ) {
    max(variant_n_values)
  } else {
    n_variants_for_circle
  }
  
  transform_variant_n <- function(
    x,
    transform = "linear"
  ) {
    x <- pmax(x, 0)
    
    if (identical(transform, "sqrt")) {
      sqrt(x)
    } else if (identical(transform, "log")) {
      log1p(x)
    } else {
      x
    }
  }
  
  transformed_n <- transform_variant_n(
    n_variants_for_circle,
    s$variant_circle_transform
  )
  
  transformed_max_n <- transform_variant_n(
    max_variants_for_scaling,
    s$variant_circle_transform
  )
  
  if (!is.finite(transformed_max_n) || transformed_max_n <= 0) {
    transformed_max_n <- 1
  }
  
  max_variant_circle_radius <- gene_radius * 0.80
  
  variant_circle_radius <- max_variant_circle_radius *
    (
      transformed_n /
        transformed_max_n
    ) *
    s$variant_circle_scale
  
  variant_circle_radius <- min(
    max_variant_circle_radius,
    max(
      0,
      variant_circle_radius
    )
  )
  
  variant_circle_poly <- make_circle_df(
    variant_circle_radius
  )
  
  plot_object <- ggplot()
  
  if (isTRUE(s$show_axis)) {
    plot_object <- plot_object +
      geom_path(
        data = axis_circles,
        
        aes(
          x = x,
          y = y,
          group = radius
        ),
        
        color = safe_colour(
          s$color_axis,
          "gray70"
        ),
        
        linewidth = 0.2,
        linetype = "dotted"
      )
  }
  
  plot_object <- plot_object +
    geom_path(
      data = make_circle_df(
        inner_radius,
        n = 100
      ),
      
      aes(
        x = x,
        y = y
      ),
      
      color = "gray50",
      linewidth = 0.4
    ) +
    
    geom_polygon(
      data = bar_polygons,
      
      aes(
        x = x,
        y = y,
        group = id,
        fill = type
      ),
      
      color = "white",
      linewidth = 0.15
    )
  
  if (isTRUE(s$show_outer_circle)) {
    plot_object <- plot_object +
      geom_path(
        data = make_circle_df(
          outer_circle_radius
        ),
        
        aes(
          x = x,
          y = y
        ),
        
        color = safe_colour(
          s$outer_circle_color,
          "grey60"
        ),
        
        linewidth = s$outer_circle_linewidth
      )
  }
  
  plot_object <- plot_object +
    geom_polygon(
      data = variant_circle_poly,
      
      aes(
        x = x,
        y = y
      ),
      
      color = safe_colour(
        s$color_gene,
        "#3EB489"
      ),
      
      fill = safe_colour(
        s$color_gene,
        "#3EB489"
      ),
      
      alpha = s$variant_circle_alpha,
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = variant_circle_df,
      
      aes(
        x = x,
        y = y,
        label = label
      ),
      
      size = s$variant_label_size *
        size_factor,
      
      color = safe_colour(
        s$variant_label_color,
        "white"
      ),
      
      fontface = s$variant_label_fontface,
      lineheight = 0.9,
      hjust = 0.5,
      vjust = 0.5
    ) +
    
    geom_text(
      data = traits,
      
      aes(
        x = x_label,
        y = y_label,
        label = display_label,
        angle = label_angle
      ),
      
      size = s$trait_label_size *
        size_factor,
      
      color = safe_colour(
        s$color_trait_label,
        "gray20"
      )
    )
  
  if (isTRUE(s$show_axis)) {
    plot_object <- plot_object +
      geom_text(
        data = axis_labels,
        
        aes(
          x = x,
          y = y,
          label = label
        ),
        
        size = s$axis_label_size *
          size_factor,
        
        hjust = 0,
        color = "gray30"
      )
  }
  
  plot_object +
    scale_fill_manual(
      values = c(
        positive = safe_colour(
          s$color_positive,
          "#CB4154"
        ),
        
        negative = safe_colour(
          s$color_negative,
          "#4682B4"
        )
      ),
      
      guide = "none"
    ) +
    
    coord_fixed(
      xlim = c(
        -max_radius,
        max_radius
      ),
      
      ylim = c(
        -max_radius,
        max_radius
      )
    ) +
    
    theme_void() +
    
    theme(
      panel.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      
      plot.background = element_rect(
        fill = "transparent",
        color = NA
      ),
      
      plot.margin = margin(
        5,
        5,
        5,
        5
      ),
      
      legend.position = "none"
    )
}

# ==============================================================================
# RADIAL PANEL ASSEMBLY
# ==============================================================================

create_panel_with_titles <- function(
    plots,
    cluster_order,
    cluster_settings,
    ncol,
    nrow,
    h_spacing,
    v_spacing,
    plot_scale,
    title_positions,
    show_titles,
    title_size,
    title_fontface,
    title_color,
    bg_color = "transparent"
) {
  if (length(plots) == 0) {
    return(NULL)
  }
  
  background_grob <- if (
    identical(bg_color, "transparent") ||
    is.na(bg_color)
  ) {
    rectGrob(
      gp = gpar(
        fill = NA,
        col = NA
      )
    )
  } else {
    rectGrob(
      gp = gpar(
        fill = bg_color,
        col = NA
      )
    )
  }
  
  ncol <- max(
    1,
    round(ncol)
  )
  
  nrow <- max(
    1,
    round(nrow)
  )
  
  cell_width <- 1 / ncol
  cell_height <- 1 / nrow
  
  scale_width <- cell_width *
    plot_scale
  
  scale_height <- cell_height *
    plot_scale
  
  x_step <- cell_width *
    (
      1 +
        h_spacing
    )
  
  y_step <- cell_height *
    (
      1 +
        v_spacing
    )
  
  total_width <- x_step *
    (
      ncol -
        1
    )
  
  total_height <- y_step *
    (
      nrow -
        1
    )
  
  left_x <- 0.5 -
    total_width / 2
  
  top_y <- 0.5 +
    total_height / 2
  
  plot_grobs <- list()
  title_grobs <- list()
  
  for (i in seq_along(cluster_order)) {
    cluster_name <- cluster_order[i]
    
    if (i > length(plots)) {
      break
    }
    
    plot_column <- ((i - 1) %% ncol) + 1
    plot_row <- ((i - 1) %/% ncol) + 1
    
    base_x <- left_x +
      (
        plot_column -
          1
      ) *
      x_step
    
    base_y <- top_y -
      (
        plot_row -
          1
      ) *
      y_step
    
    x_offset <- cluster_settings[[cluster_name]]$x_offset %||% 0
    y_offset <- cluster_settings[[cluster_name]]$y_offset %||% 0
    
    final_x <- base_x +
      x_offset
    
    final_y <- base_y +
      y_offset
    
    plot_viewport <- viewport(
      x = unit(
        final_x,
        "npc"
      ),
      
      y = unit(
        final_y,
        "npc"
      ),
      
      width = unit(
        scale_width,
        "npc"
      ),
      
      height = unit(
        scale_height,
        "npc"
      ),
      
      just = c(
        "center",
        "center"
      )
    )
    
    plot_grobs[[length(plot_grobs) + 1]] <- gTree(
      children = gList(
        plots[[i]]
      ),
      
      vp = plot_viewport
    )
    
    if (isTRUE(show_titles)) {
      title_row <- title_positions[
        title_positions$Cluster == cluster_name,
        ,
        drop = FALSE
      ]
      
      if (nrow(title_row) > 0) {
        title_x <- final_x +
          (
            title_row$X[1] -
              0.5
          ) *
          scale_width
        
        title_y <- final_y +
          (
            title_row$Y[1] -
              0.5
          ) *
          scale_height
        
        title_grobs[[length(title_grobs) + 1]] <- textGrob(
          label = title_row$Title[1],
          
          x = unit(
            title_x,
            "npc"
          ),
          
          y = unit(
            title_y,
            "npc"
          ),
          
          just = "center",
          
          gp = gpar(
            fontsize = title_size,
            fontface = title_fontface,
            
            col = safe_colour(
              title_color,
              "black"
            ),
            
            lineheight = 0.9
          )
        )
      }
    }
  }
  
  gTree(
    children = do.call(
      gList,
      
      c(
        list(background_grob),
        plot_grobs,
        title_grobs
      )
    )
  )
}

# ==============================================================================
# DOT-PLOT INPUT NORMALISATION
# ==============================================================================

normalize_dot_input <- function(source_df) {
  source_names <- names(source_df)

  # New gchromVAR/PPA long-format input.
  new_schema <- c(
    "Cluster", "atlas", "cell_type", "description",
    "z", "p"
  )

  if (all(new_schema %in% source_names)) {
    return(
      source_df %>%
        transmute(
          atlas_raw = trimws(as.character(atlas)),
          atlas = dplyr::recode(
            atlas_raw,
            "Coronary Plaque (Turner 2022)" = "Coronary plaque",
            "Heart Atlas (CATLAS)" = "Heart atlas",
            "Human Multi-tissue (CATLAS, adult)" = "Whole-body (CATLAS)",
            .default = atlas_raw
          ),
          cell = trimws(as.character(cell_type)),
          cell_atlas_name = trimws(as.character(description)),
          cluster_raw = trimws(as.character(Cluster)),
          cluster = dplyr::recode(
            cluster_raw,
            "Endothelial dysfunction" = "Endothelial",
            "Cardiac remodeling" = "Cardiac",
            .default = cluster_raw
          ),
          z = as.numeric(z),
          p_nominal = as.numeric(p),
          source_bonferroni_threshold = if (
            "bonferroni_threshold" %in% source_names
          ) {
            as.numeric(source_df$bonferroni_threshold)
          } else {
            NA_real_
          }
        )
    )
  }

  # Backwards-compatible support for the previous app schema.
  old_schema <- c(
    "atlas", "cell", "cell_atlas_name", "cluster",
    "z", "p_nominal"
  )

  if (all(old_schema %in% source_names)) {
    return(
      source_df %>%
        mutate(
          atlas = trimws(as.character(atlas)),
          cell = trimws(as.character(cell)),
          cell_atlas_name = trimws(as.character(cell_atlas_name)),
          cluster = trimws(as.character(cluster)),
          z = as.numeric(z),
          p_nominal = as.numeric(p_nominal),
          source_bonferroni_threshold = NA_real_
        )
    )
  }

  stop(
    paste(
      "Unsupported dot-matrix schema.",
      "For the revised data use columns:",
      "Cluster, atlas, cell_type, description, z, p",
      "(extra columns such as q_fdr_cluster are allowed but ignored)."
    )
  )
}

# ==============================================================================
# DOT-PLOT DATA PREPARATION
# ==============================================================================

prepare_dot_plot_data <- function(
    source_df,
    settings
) {
  source_df <- normalize_dot_input(source_df)

  required_columns <- c(
    "atlas", "cell", "cell_atlas_name", "cluster",
    "z", "p_nominal"
  )

  missing_columns <- setdiff(required_columns, names(source_df))
  if (length(missing_columns) > 0) {
    stop(
      "The normalized dot-matrix data are missing these columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  source_df <- source_df %>%
    filter(
      atlas %in% settings$atlas_order,
      cluster %in% settings$cluster_order
    )

  if (nrow(source_df) == 0) {
    stop("No dot-matrix rows remain after atlas and cluster filtering.")
  }

  duplicated_plot_rows <- source_df %>%
    count(atlas, cell, cluster, name = "n") %>%
    filter(n > 1)

  if (nrow(duplicated_plot_rows) > 0) {
    stop(
      "More than one row exists for at least one atlas-cell-cluster combination."
    )
  }

  # Plaque and heart: retain every available cell.
  always_show_atlases <- c("Coronary plaque", "Heart atlas")

  retained_standard_cells <- source_df %>%
    filter(atlas %in% always_show_atlases) %>%
    distinct(atlas, cell)

  # Whole-body CATLAS: retain cells significant for >=1 process according to
  # the user-selected rule. The full eight-process profile is then plotted.
  catlas_df <- source_df %>%
    filter(atlas == "Whole-body (CATLAS)")

  bonferroni_cutoff <- settings$nominal_alpha / settings$bonferroni_n_tests

  if (identical(settings$catlas_filter_mode, "bonferroni")) {
    retained_catlas_cells <- catlas_df %>%
      group_by(atlas, cell) %>%
      summarise(
        keep = any(!is.na(p_nominal) & p_nominal < bonferroni_cutoff),
        .groups = "drop"
      ) %>%
      filter(keep) %>%
      select(atlas, cell)
  } else {
    retained_catlas_cells <- catlas_df %>%
      group_by(atlas, cell) %>%
      summarise(
        keep = any(!is.na(p_nominal) & p_nominal < settings$nominal_alpha),
        .groups = "drop"
      ) %>%
      filter(keep) %>%
      select(atlas, cell)
  }

  retained_cells <- bind_rows(
    retained_standard_cells,
    retained_catlas_cells
  ) %>%
    distinct(atlas, cell)

  plot_df <- source_df %>%
    inner_join(retained_cells, by = c("atlas", "cell"))

  if (nrow(plot_df) == 0) {
    stop(
      "No cells satisfy the current dot-matrix display criteria. ",
      "Try nominal P or decrease the Bonferroni denominator / increase alpha."
    )
  }

  # Build cell order. Text-box entries are optional ordering preferences, not
  # filters: every retained cell is appended if it is not explicitly listed.
  available_cells <- plot_df %>%
    distinct(atlas, cell, cell_atlas_name)

  build_order <- function(atlas_name) {
    available <- available_cells %>%
      filter(atlas == atlas_name)

    if (nrow(available) == 0) {
      return(data.frame())
    }

    requested <- settings$cell_order[[atlas_name]] %||% character(0)
    requested <- unique(trimws(as.character(requested)))
    requested <- requested[nzchar(requested)]
    requested <- requested[requested %in% available$cell]

    remaining <- available$cell[!available$cell %in% requested]
    final_order <- c(requested, remaining)

    data.frame(
      atlas = atlas_name,
      cell = final_order,
      cell_order = seq_along(final_order),
      stringsAsFactors = FALSE
    )
  }

  selected_atlas_cells <- bind_rows(
    lapply(settings$atlas_order, build_order)
  )

  active_atlas_order <- settings$atlas_order[
    settings$atlas_order %in% selected_atlas_cells$atlas
  ]

  position_tables <- vector("list", length(active_atlas_order))
  current_x <- 0

  for (i in seq_along(active_atlas_order)) {
    atlas_name <- active_atlas_order[i]
    atlas_table <- selected_atlas_cells %>%
      filter(atlas == atlas_name) %>%
      arrange(cell_order)

    if (i > 1) current_x <- current_x + settings$atlas_gap

    atlas_table$x_pos <- current_x + seq_len(nrow(atlas_table))
    current_x <- max(atlas_table$x_pos)
    position_tables[[i]] <- atlas_table
  }

  cell_position_table <- bind_rows(position_tables)

  cell_label_source <- available_cells %>%
    distinct(atlas, cell, .keep_all = TRUE)

  cell_position_table <- cell_position_table %>%
    left_join(cell_label_source, by = c("atlas", "cell")) %>%
    mutate(
      # Revised dataset descriptions are the preferred human-readable labels.
      default_label = ifelse(
        !is.na(cell_atlas_name) & nzchar(cell_atlas_name),
        cell_atlas_name,
        cell
      ),
      display_label = default_label
    )

  # Custom labels support either:
  #   cell_type = My label
  # or atlas-specific:
  #   Whole-body (CATLAS)|cell_type = My label
  if (length(settings$cell_labels) > 0) {
    global_key <- cell_position_table$cell
    atlas_key <- paste(cell_position_table$atlas, cell_position_table$cell, sep = "|")

    global_idx <- match(global_key, names(settings$cell_labels))
    atlas_idx <- match(atlas_key, names(settings$cell_labels))

    has_global <- !is.na(global_idx)
    has_atlas <- !is.na(atlas_idx)

    cell_position_table$display_label[has_global] <- unname(
      settings$cell_labels[global_idx[has_global]]
    )
    cell_position_table$display_label[has_atlas] <- unname(
      settings$cell_labels[atlas_idx[has_atlas]]
    )
  }

  plot_df <- plot_df %>%
    inner_join(
      cell_position_table %>%
        select(atlas, cell, x_pos, display_label),
      by = c("atlas", "cell")
    )

  # Significance is based only on the raw P value. Bonferroni is applied
  # once using the user-specified denominator: alpha / N tests.
  plot_df <- plot_df %>%
    mutate(
      bonferroni_cutoff = settings$nominal_alpha / settings$bonferroni_n_tests,
      significance_plot = case_when(
        !is.na(p_nominal) & p_nominal < bonferroni_cutoff ~ "Bonferroni P < 0.05",
        !is.na(p_nominal) & p_nominal < settings$nominal_alpha ~ "Nominal P < 0.05",
        TRUE ~ "P \u2265 0.05"
      ),
      significance_plot = factor(
        significance_plot,
        levels = c("P \u2265 0.05", "Nominal P < 0.05", "Bonferroni P < 0.05")
      )
    )

  cluster_position_table <- data.frame(
    cluster = settings$cluster_order,
    y_pos = rev(seq_along(settings$cluster_order)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      display_label_raw = vapply(
        cluster,
        function(cluster_name) {
          settings$cluster_labels[[cluster_name]] %||% cluster_name
        },
        character(1)
      ),
      display_label = if (isTRUE(settings$wrap_y_labels)) {
        stringr::str_wrap(display_label_raw, width = settings$y_label_wrap_width)
      } else {
        display_label_raw
      }
    )

  plot_df <- plot_df %>%
    inner_join(cluster_position_table, by = "cluster")

  nonsignificant_df <- plot_df %>%
    filter(significance_plot == "P \u2265 0.05")

  significant_df <- plot_df %>%
    filter(significance_plot != "P \u2265 0.05")

  grid_df <- expand.grid(
    x_pos = cell_position_table$x_pos,
    y_pos = cluster_position_table$y_pos,
    stringsAsFactors = FALSE
  )

  atlas_ranges <- cell_position_table %>%
    group_by(atlas) %>%
    summarise(
      xmin = min(x_pos) - 0.45,
      xmax = max(x_pos) + 0.45,
      xmid = mean(range(x_pos)),
      .groups = "drop"
    ) %>%
    mutate(
      atlas = factor(atlas, levels = active_atlas_order),
      atlas_label = vapply(
        as.character(atlas),
        function(atlas_name) {
          settings$atlas_labels[[atlas_name]] %||% atlas_name
        },
        character(1)
      )
    ) %>%
    arrange(atlas)

  atlas_separators <- numeric(0)
  if (nrow(atlas_ranges) > 1) {
    atlas_separators <- vapply(
      seq_len(nrow(atlas_ranges) - 1),
      function(i) {
        mean(c(atlas_ranges$xmax[i], atlas_ranges$xmin[i + 1]))
      },
      numeric(1)
    )
  }

  list(
    plot_df = plot_df,
    nonsignificant_df = nonsignificant_df,
    significant_df = significant_df,
    grid_df = grid_df,
    cell_position_table = cell_position_table,
    cluster_position_table = cluster_position_table,
    atlas_ranges = atlas_ranges,
    atlas_separators = atlas_separators,
    bonferroni_cutoff = bonferroni_cutoff,
    n_catlas_cells_retained = nrow(retained_catlas_cells)
  )
}

# ==============================================================================
# DOT-PLOT FUNCTION
# ==============================================================================

create_coloured_dot_heatmap <- function(
    prepared,
    settings
) {
  nonsignificant_df <- prepared$nonsignificant_df
  significant_df <- prepared$significant_df
  grid_df <- prepared$grid_df
  cell_position_table <- prepared$cell_position_table
  cluster_position_table <- prepared$cluster_position_table
  atlas_ranges <- prepared$atlas_ranges
  atlas_separators <- prepared$atlas_separators
  
  x_breaks <- cell_position_table$x_pos
  x_labels <- cell_position_table$display_label
  
  y_breaks <- cluster_position_table$y_pos
  y_labels <- cluster_position_table$display_label
  
  y_top <- max(
    cluster_position_table$y_pos
  )
  
  atlas_bracket_y <- y_top +
    settings$atlas_bracket_offset
  
  atlas_label_y <- y_top +
    settings$atlas_label_offset
  
  plot_object <- ggplot() +
    geom_tile(
      data = grid_df,
      
      aes(
        x = x_pos,
        y = y_pos
      ),
      
      fill = safe_colour(
        settings$tile_fill,
        "grey93"
      ),
      
      colour = safe_colour(
        settings$tile_border_colour,
        "white"
      ),
      
      linewidth = settings$tile_border_width,
      width = settings$tile_width,
      height = settings$tile_height
    ) +
    
    geom_point(
      data = nonsignificant_df,
      
      aes(
        x = x_pos,
        y = y_pos,
        size = significance_plot,
        fill = z
      ),
      
      shape = 21,
      
      colour = safe_colour(
        settings$nonsignificant_outline_colour,
        "grey25"
      ),
      
      stroke = settings$nonsignificant_outline_width
    ) +
    
    geom_point(
      data = significant_df,
      
      aes(
        x = x_pos,
        y = y_pos,
        size = significance_plot,
        fill = z
      ),
      
      shape = 21,
      
      colour = safe_colour(
        settings$significant_outline_colour,
        "grey25"
      ),
      
      stroke = settings$significant_outline_width
    )
  
  if (length(atlas_separators) > 0) {
    plot_object <- plot_object +
      geom_vline(
        xintercept = atlas_separators,
        colour = "white",
        linewidth = settings$atlas_gap_line_width
      )
    
    if (isTRUE(settings$show_atlas_separator)) {
      plot_object <- plot_object +
        geom_vline(
          xintercept = atlas_separators,
          
          colour = safe_colour(
            settings$atlas_separator_colour,
            "grey65"
          ),
          
          linewidth = settings$atlas_separator_width,
          linetype = settings$atlas_separator_linetype
        )
    }
  }
  
  if (isTRUE(settings$show_atlas_brackets)) {
    plot_object <- plot_object +
      annotate(
        geom = "segment",
        x = atlas_ranges$xmin,
        xend = atlas_ranges$xmax,
        y = atlas_bracket_y,
        yend = atlas_bracket_y,
        
        colour = safe_colour(
          settings$atlas_bracket_colour,
          "grey45"
        ),
        
        linewidth = settings$atlas_bracket_width
      )
  }
  
  if (isTRUE(settings$show_atlas_headers)) {
    plot_object <- plot_object +
      annotate(
        geom = "text",
        x = atlas_ranges$xmid,
        y = atlas_label_y,
        label = atlas_ranges$atlas_label,
        
        size = settings$atlas_header_size /
          ggplot2::.pt,
        
        fontface = settings$atlas_header_fontface,
        
        colour = safe_colour(
          settings$atlas_header_colour,
          "grey20"
        )
      )
  }
  
  plot_object +
    scale_fill_gradient2(
      low = safe_colour(
        settings$divergent_colour_low,
        "steelblue4"
      ),
      
      mid = safe_colour(
        settings$divergent_colour_mid,
        "white"
      ),
      
      high = safe_colour(
        settings$divergent_colour_high,
        "firebrick3"
      ),
      
      midpoint = 0,
      
      limits = c(
        settings$z_min,
        settings$z_max
      ),
      
      oob = scales::squish,
      breaks = settings$z_breaks,
      labels = settings$z_breaks,
      name = settings$z_legend_title
    ) +
    
    scale_size_manual(
      values = settings$dot_sizes[
        c(
          "P \u2265 0.05",
          "Nominal P < 0.05",
          "Bonferroni P < 0.05"
        )
      ],
      
      breaks = c(
        "P \u2265 0.05",
        "Nominal P < 0.05",
        "Bonferroni P < 0.05"
      ),
      
      name = settings$significance_legend_title,
      drop = FALSE
    ) +
    
    labs(
      x = NULL,
      y = settings$y_axis_title
    ) +
    
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      
      expand = expansion(
        mult = c(
          settings$x_expand_left,
          settings$x_expand_right
        )
      )
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      
      expand = expansion(
        mult = c(
          settings$y_expand_bottom,
          settings$y_expand_top
        )
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_minimal(
      base_size = settings$base_text_size
    ) +
    
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      axis.title.x = element_blank(),
      
      axis.title.y = element_text(
        size = settings$y_axis_title_size,
        face = settings$y_axis_title_face,
        
        colour = safe_colour(
          settings$y_axis_title_colour,
          "grey20"
        ),
        
        angle = settings$y_axis_title_angle,
        hjust = settings$y_axis_title_hjust,
        vjust = settings$y_axis_title_vjust,
        
        margin = margin(
          r = settings$y_axis_title_margin
        )
      ),
      
      axis.ticks = element_blank(),
      
      axis.text.x = element_text(
        angle = settings$x_text_angle,
        hjust = settings$x_text_hjust,
        vjust = settings$x_text_vjust,
        size = settings$x_text_size,
        
        colour = safe_colour(
          settings$x_text_colour,
          "grey20"
        ),
        
        face = settings$x_text_face
      ),
      
      axis.text.y = element_text(
        size = settings$y_text_size,
        
        colour = safe_colour(
          settings$y_text_colour,
          "grey20"
        ),
        
        hjust = settings$y_text_hjust,
        vjust = settings$y_text_vjust,
        face = settings$y_text_face,
        lineheight = settings$y_text_lineheight
      ),
      
      legend.position = settings$legend_position,
      
      legend.title = element_text(
        size = settings$legend_title_size,
        face = settings$legend_title_face
      ),
      
      legend.text = element_text(
        size = settings$legend_text_size
      ),
      
      legend.key.height = unit(
        settings$legend_key_height,
        "cm"
      ),
      
      plot.margin = margin(
        t = settings$plot_margin_top,
        r = settings$plot_margin_right,
        b = settings$plot_margin_bottom,
        l = settings$plot_margin_left,
        unit = "pt"
      )
    ) +
    
    guides(
      fill = guide_colorbar(
        order = 1,
        title.position = "top",
        barwidth = settings$colorbar_width,
        barheight = settings$colorbar_height
      ),
      
      size = guide_legend(
        order = 2,
        title.position = "top",
        
        override.aes = list(
          shape = 16,
          
          colour = safe_colour(
            settings$significance_legend_colour,
            "grey60"
          ),
          
          fill = safe_colour(
            settings$significance_legend_colour,
            "grey60"
          ),
          
          stroke = 0
        )
      )
    )
}

# ==============================================================================
# COMBINED FIGURE
# ==============================================================================

create_combined_figure <- function(
    radial_grob,
    dot_grob,
    top_height,
    bottom_height,
    panel_a_label,
    panel_a_x,
    panel_a_y,
    panel_a_size,
    panel_a_colour,
    panel_b_label,
    panel_b_x,
    panel_b_y,
    panel_b_size,
    panel_b_colour
) {
  top_height <- max(
    0.01,
    top_height
  )
  
  bottom_height <- max(
    0.01,
    bottom_height
  )
  
  total_height <- top_height +
    bottom_height
  
  top_fraction <- top_height /
    total_height
  
  bottom_fraction <- bottom_height /
    total_height
  
  combined_base <- arrangeGrob(
    radial_grob,
    dot_grob,
    
    ncol = 1,
    
    heights = c(
      top_height,
      bottom_height
    )
  )
  
  label_grobs <- list()
  
  if (
    !is.null(panel_a_label) &&
    nzchar(panel_a_label)
  ) {
    panel_a_global_y <- bottom_fraction +
      top_fraction *
      panel_a_y
    
    label_grobs[[length(label_grobs) + 1]] <- textGrob(
      label = panel_a_label,
      
      x = unit(
        panel_a_x,
        "npc"
      ),
      
      y = unit(
        panel_a_global_y,
        "npc"
      ),
      
      just = c(
        "left",
        "top"
      ),
      
      gp = gpar(
        fontsize = panel_a_size,
        fontface = "bold",
        
        col = safe_colour(
          panel_a_colour,
          "black"
        )
      )
    )
  }
  
  if (
    !is.null(panel_b_label) &&
    nzchar(panel_b_label)
  ) {
    panel_b_global_y <- bottom_fraction *
      panel_b_y
    
    label_grobs[[length(label_grobs) + 1]] <- textGrob(
      label = panel_b_label,
      
      x = unit(
        panel_b_x,
        "npc"
      ),
      
      y = unit(
        panel_b_global_y,
        "npc"
      ),
      
      just = c(
        "left",
        "top"
      ),
      
      gp = gpar(
        fontsize = panel_b_size,
        fontface = "bold",
        
        col = safe_colour(
          panel_b_colour,
          "black"
        )
      )
    )
  }
  
  gTree(
    children = do.call(
      gList,
      
      c(
        list(combined_base),
        label_grobs
      )
    )
  )
}

# ==============================================================================
# RADIAL CLUSTER UI
# ==============================================================================

make_cluster_controls <- function(cluster) {
  defaults <- cluster_defaults()
  
  tagList(
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_n_top_genes"),
          "Top variants/SNPs",
          defaults$n_top_genes,
          min = 1,
          step = 1
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_gene_weight_threshold"),
          "Variant threshold",
          defaults$gene_weight_threshold,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_n_top_traits"),
          "Top traits per direction",
          defaults$n_top_traits,
          min = 1,
          step = 1
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_trait_weight_threshold"),
          "Trait threshold",
          defaults$trait_weight_threshold,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_gene_radius"),
          "Inner guide radius",
          defaults$gene_radius,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_inner_radius_offset"),
          "Inner radius offset",
          defaults$inner_radius_offset,
          step = 0.01
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_bar_scale"),
          "Bar scale",
          defaults$bar_scale,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_bar_width"),
          "Bar width radians",
          defaults$bar_width,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_outer_radius_extend"),
          "Outer-radius extension",
          defaults$outer_radius_extend,
          step = 0.01
        )
      )
    ),
    
    fluidRow(
      column(
        6,
        
        sliderInput(
          paste0(cluster, "_pos_range"),
          "Positive angle range",
          min = 0,
          max = 360,
          value = defaults$pos_range
        )
      ),
      
      column(
        6,
        
        sliderInput(
          paste0(cluster, "_neg_range"),
          "Negative angle range",
          min = 0,
          max = 360,
          value = defaults$neg_range
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        checkboxInput(
          paste0(cluster, "_show_axis"),
          "Show axis",
          defaults$show_axis
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_axis_intervals"),
          "Axis intervals",
          defaults$axis_intervals,
          min = 2,
          step = 1
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_axis_label_size"),
          "Axis-label size",
          defaults$axis_label_size,
          step = 0.1
        )
      )
    ),
    
    tags$hr(),
    
    h5("Outermost circle"),
    
    fluidRow(
      column(
        4,
        
        checkboxInput(
          paste0(cluster, "_show_outer_circle"),
          "Show outer circle",
          defaults$show_outer_circle
        )
      ),
      
      column(
        4,
        
        textInput(
          paste0(cluster, "_outer_circle_color"),
          "Outer-circle colour",
          defaults$outer_circle_color
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_outer_circle_linewidth"),
          "Outer-circle width",
          defaults$outer_circle_linewidth,
          min = 0,
          step = 0.05
        )
      )
    ),
    
    tags$hr(),
    
    h5("Central variant circle"),
    
    fluidRow(
      column(
        4,
        
        selectInput(
          paste0(cluster, "_variant_circle_transform"),
          "Circle-size transform",
          
          choices = c(
            "linear",
            "sqrt",
            "log"
          ),
          
          selected = defaults$variant_circle_transform
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_variant_circle_scale"),
          "Circle-radius scale",
          defaults$variant_circle_scale,
          min = 0.01,
          step = 0.05
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_variant_circle_alpha"),
          "Circle alpha",
          defaults$variant_circle_alpha,
          min = 0,
          max = 1,
          step = 0.05
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        checkboxInput(
          paste0(cluster, "_variant_label_two_lines"),
          "Variant on second line",
          defaults$variant_label_two_lines
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_variant_label_size"),
          "Central-label size",
          defaults$variant_label_size,
          min = 0.1,
          step = 0.1
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_variant_label_wrap_width"),
          "Label wrap width",
          defaults$variant_label_wrap_width,
          min = 1,
          step = 1
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_variant_label_x"),
          "Central-label X",
          defaults$variant_label_x,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_variant_label_y"),
          "Central-label Y",
          defaults$variant_label_y,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        textInput(
          paste0(cluster, "_variant_label_color"),
          "Central-label colour",
          defaults$variant_label_color
        )
      )
    ),
    
    selectInput(
      paste0(cluster, "_variant_label_fontface"),
      "Central-label face",
      
      choices = c(
        "plain",
        "bold",
        "italic",
        "bold.italic"
      ),
      
      selected = defaults$variant_label_fontface
    ),
    
    tags$hr(),
    
    h5("Trait labels"),
    
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_trait_label_size"),
          "Trait-label size",
          defaults$trait_label_size,
          step = 0.1
        )
      ),
      
      column(
        4,
        
        selectInput(
          paste0(cluster, "_trait_label_spacing_mode"),
          "Label-spacing mode",
          
          choices = c(
            "Global short/medium/long" = "global",
            "Manual per trait" = "manual"
          ),
          
          selected = defaults$trait_label_spacing_mode
        )
      ),
      
      column(
        4,
        
        checkboxInput(
          paste0(cluster, "_trait_labels_radial"),
          "Radial labels",
          defaults$trait_labels_radial
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_trait_label_offset_short"),
          "Short-label offset",
          defaults$trait_label_offset_short,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_trait_label_offset_medium"),
          "Medium-label offset",
          defaults$trait_label_offset_medium,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_trait_label_offset_long"),
          "Long-label offset",
          defaults$trait_label_offset_long,
          step = 0.01
        )
      )
    ),
    
    textAreaInput(
      paste0(cluster, "_trait_label_manual_offsets"),
      "Manual trait offsets: trait = offset",
      rows = 8,
      
      value = paste(
        "hscrp = 0.18",
        "ldl = 0.18",
        "hdl = 0.18",
        "bmi = 0.18",
        "sbp = 0.24",
        "dbp = 0.24",
        sep = "\n"
      )
    ),
    
    helpText(
      paste(
        "Manual mode accepts the original trait name,",
        "the name without _pos/_neg, or the displayed label.",
        "Traits not listed use the global offset."
      )
    ),
    
    numericInput(
      paste0(cluster, "_trait_label_angle_adjust"),
      "Trait-label angle adjustment",
      defaults$trait_label_angle_adjust,
      step = 1
    ),
    
    tags$hr(),
    
    h5("Position and colours"),
    
    fluidRow(
      column(
        4,
        
        numericInput(
          paste0(cluster, "_x_offset"),
          "X offset",
          defaults$x_offset,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        numericInput(
          paste0(cluster, "_y_offset"),
          "Y offset",
          defaults$y_offset,
          step = 0.01
        )
      ),
      
      column(
        4,
        
        textInput(
          paste0(cluster, "_color_gene"),
          "Variant-circle colour",
          defaults$color_gene
        )
      )
    ),
    
    fluidRow(
      column(
        4,
        
        textInput(
          paste0(cluster, "_color_positive"),
          "Positive colour",
          defaults$color_positive
        )
      ),
      
      column(
        4,
        
        textInput(
          paste0(cluster, "_color_negative"),
          "Negative colour",
          defaults$color_negative
        )
      ),
      
      column(
        4,
        
        textInput(
          paste0(cluster, "_color_axis"),
          "Axis colour",
          defaults$color_axis
        )
      )
    ),
    
    textInput(
      paste0(cluster, "_color_trait_label"),
      "Trait-label colour",
      defaults$color_trait_label
    )
  )
}

# ==============================================================================
# UI
# ==============================================================================

ui <- fluidPage(
  titlePanel(
    "Integrated radial and cell-enrichment figure"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      tabsetPanel(
        id = "control_tabs",
        
        tabPanel(
          "Data",
          
          h4("Panel A: radial data"),
          
          checkboxInput(
            "use_example_circular",
            "Use example radial data",
            TRUE
          ),
          
          fileInput(
            "gene_file",
            "Variant/gene weights CSV"
          ),
          
          fileInput(
            "trait_file",
            "Trait weights CSV"
          ),
          
          textAreaInput(
            "variant_n_text",
            "Variants per cluster: Cluster = N",
            rows = 8,
            
            value = paste(
              paste0(
                "X",
                1:8,
                " = 10"
              ),
              
              collapse = "\n"
            )
          ),
          
          helpText(
            paste(
              "These values determine the number displayed",
              "inside each central radial circle and its relative size."
            )
          ),
          
          tags$hr(),
          
          h4("Panel B: dot-matrix data"),
          
          checkboxInput(
            "use_example_dot",
            "Use example dot-matrix data",
            FALSE
          ),
          
          fileInput(
            "dot_file",
            "Dot-matrix CSV"
          ),
          
          helpText(
            paste(
              "Revised input columns:",
              "Cluster, atlas, cell_type, description, z, p.",
              "Extra columns (including q_fdr_cluster) are allowed but ignored.",
              "The previous app schema is still accepted."
            )
          ),
          
          tags$hr(),
          
          h4("Radial label remapping"),
          
          textAreaInput(
            "gene_labels_text",
            "Variant/gene labels: old = new",
            rows = 3,
            value = ""
          ),
          
          textAreaInput(
            "trait_labels_text",
            "Trait labels: old = new",
            rows = 10,
            
            value = paste(
              "hscrp = hs-CRP",
              "ldl = LDL-C",
              "hdl = HDL-C",
              "bmi = BMI",
              "sbp = Systolic BP",
              "dbp = Diastolic BP",
              "glucose = Glucose",
              "hba1c = HbA1c",
              "trig = Triglycerides",
              "chol = Total Cholesterol",
              sep = "\n"
            )
          )
        ),
        
        tabPanel(
          "Combined layout",
          
          numericInput(
            "top_panel_height",
            "Panel A relative height",
            2.1,
            min = 0.05,
            step = 0.05
          ),
          
          numericInput(
            "bottom_panel_height",
            "Panel B relative height",
            1.2,
            min = 0.05,
            step = 0.05
          ),
          
          tags$hr(),
          
          h4("Panel A label"),
          
          textInput(
            "panel_a_label",
            "Label",
            "A"
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "panel_a_label_x",
                "X",
                0.01,
                min = 0,
                max = 1,
                step = 0.005
              )
            ),
            
            column(
              4,
              
              numericInput(
                "panel_a_label_y",
                "Y within Panel A",
                0.99,
                min = 0,
                max = 1,
                step = 0.005
              )
            ),
            
            column(
              4,
              
              numericInput(
                "panel_a_label_size",
                "Size",
                18,
                min = 1,
                step = 0.5
              )
            )
          ),
          
          textInput(
            "panel_a_label_colour",
            "Colour",
            "black"
          ),
          
          tags$hr(),
          
          h4("Panel B label"),
          
          textInput(
            "panel_b_label",
            "Label",
            "B"
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "panel_b_label_x",
                "X",
                0.01,
                min = 0,
                max = 1,
                step = 0.005
              )
            ),
            
            column(
              4,
              
              numericInput(
                "panel_b_label_y",
                "Y within Panel B",
                0.99,
                min = 0,
                max = 1,
                step = 0.005
              )
            ),
            
            column(
              4,
              
              numericInput(
                "panel_b_label_size",
                "Size",
                18,
                min = 1,
                step = 0.5
              )
            )
          ),
          
          textInput(
            "panel_b_label_colour",
            "Colour",
            "black"
          ),
          
          tags$hr(),
          
          numericInput(
            "preview_height_px",
            "Preview height in pixels",
            1100,
            min = 300,
            step = 50
          )
        ),
        
        tabPanel(
          "Radial layout",
          
          textInput(
            "cluster_order",
            "Cluster order, comma-separated",
            
            paste(
              default_radial_cluster_order,
              collapse = ","
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "panel_ncol",
                "Columns",
                4,
                min = 1,
                step = 1
              )
            ),
            
            column(
              6,
              
              numericInput(
                "panel_nrow",
                "Rows",
                2,
                min = 1,
                step = 1
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "panel_h_spacing",
                "Horizontal spacing",
                0,
                min = -0.8,
                max = 1,
                step = 0.01
              )
            ),
            
            column(
              6,
              
              numericInput(
                "panel_v_spacing",
                "Vertical spacing",
                -0.15,
                min = -0.8,
                max = 1,
                step = 0.01
              )
            )
          ),
          
          numericInput(
            "plot_scale",
            "Radial-plot scale",
            1,
            min = 0.05,
            step = 0.01
          ),
          
          selectInput(
            "panel_background",
            "Radial-panel background",
            
            choices = c(
              "transparent",
              "white",
              "gray95",
              "gray90"
            ),
            
            selected = "transparent"
          ),
          
          checkboxInput(
            "show_titles",
            "Show radial titles",
            TRUE
          ),
          
          numericInput(
            "title_size",
            "Radial-title size",
            12,
            min = 1,
            step = 0.1
          ),
          
          selectInput(
            "title_fontface",
            "Title face",
            
            choices = c(
              "plain",
              "bold",
              "italic",
              "bold.italic"
            ),
            
            selected = "bold"
          ),
          
          textInput(
            "title_color",
            "Title colour",
            "black"
          ),
          
          textAreaInput(
            "title_mapping_text",
            "Radial titles: Cluster = Title",
            rows = 8,
            
            value = paste(
              "X1 = Obesity",
              "X2 = MDWO",
              "X3 = Endothelial",
              "X4 = Cardiac remodeling",
              "X5 = Coagulation",
              "X6 = LDL",
              "X7 = TRL",
              "X8 = Lp(a)",
              sep = "\n"
            )
          ),
          
          numericInput(
            "title_wrap_width_default",
            "Default title wrap width",
            16,
            min = 1,
            step = 1
          ),
          
          textAreaInput(
            "title_wrap_width_mapping_text",
            "Individual title wrap widths: Cluster = width",
            rows = 8,
            
            value = paste(
              "X1 = 16",
              "X2 = 16",
              "X3 = 16",
              "X4 = 14",
              "X5 = 16",
              "X6 = 16",
              "X7 = 16",
              "X8 = 16",
              sep = "\n"
            )
          ),
          
          numericInput(
            "title_x_default",
            "Default title X",
            0.5,
            step = 0.01
          ),
          
          numericInput(
            "title_y_default",
            "Default title Y",
            0.86,
            step = 0.01
          ),
          
          textAreaInput(
            "title_x_mapping_text",
            "Optional title X overrides: Cluster = X",
            rows = 4,
            value = ""
          ),
          
          textAreaInput(
            "title_y_mapping_text",
            "Optional title Y overrides: Cluster = Y",
            rows = 4,
            value = ""
          )
        ),
        
        tabPanel(
          "Radial clusters",
          
          wellPanel(
            h4("Propagate radial settings"),
            
            fluidRow(
              column(
                4,
                
                selectInput(
                  "propagate_source_cluster",
                  "Source cluster",
                  choices = radial_clusters,
                  selected = "X1"
                )
              ),
              
              column(
                5,
                
                selectInput(
                  "propagate_setting_group",
                  "Settings to copy",
                  choices = names(radial_setting_groups),
                  selected = "All visual and filtering settings"
                )
              ),
              
              column(
                3,
                
                tags$div(
                  style = "margin-top: 25px;",
                  
                  actionButton(
                    "propagate_cluster_settings",
                    "Propagate to other 7 clusters",
                    class = "btn-primary"
                  )
                )
              )
            ),
            
            helpText(
              paste(
                "Copies the selected settings from the source cluster",
                "to all other radial clusters.",
                "The values in Variants per cluster are not changed."
              )
            ),
            
            verbatimTextOutput(
              "propagation_status"
            )
          ),
          
          tabsetPanel(
            id = "cluster_tabs",
            
            tabPanel(
              "X1",
              make_cluster_controls("X1")
            ),
            
            tabPanel(
              "X2",
              make_cluster_controls("X2")
            ),
            
            tabPanel(
              "X3",
              make_cluster_controls("X3")
            ),
            
            tabPanel(
              "X4",
              make_cluster_controls("X4")
            ),
            
            tabPanel(
              "X5",
              make_cluster_controls("X5")
            ),
            
            tabPanel(
              "X6",
              make_cluster_controls("X6")
            ),
            
            tabPanel(
              "X7",
              make_cluster_controls("X7")
            ),
            
            tabPanel(
              "X8",
              make_cluster_controls("X8")
            )
          )
        ),
        
        tabPanel(
          "Dot data layout",
          
          textInput(
            "dot_atlas_order",
            "Atlas order, comma-separated",
            
            paste(
              default_dot_atlas_order,
              collapse = ","
            )
          ),
          
          textAreaInput(
            "dot_atlas_labels",
            "Atlas headers: stored value = displayed label",
            rows = 5,
            
            value = paste(
              "Coronary plaque = Coronary Plaque (Turner et al. 2022)",
              "Heart atlas = Heart Atlas (CATlas)",
              "Whole-body (CATLAS) = Human Multi-tissues (CATlas)",
              sep = "\n"
            )
          ),
          
          numericInput(
            "dot_atlas_gap",
            "Gap between atlas groups",
            1.2,
            min = 0,
            step = 0.1
          ),
          
          tags$hr(),
          
          textInput(
            "dot_cluster_order",
            "Cluster order, comma-separated",
            
            paste(
              default_dot_cluster_order,
              collapse = ","
            )
          ),
          
          textAreaInput(
            "dot_cluster_labels",
            "Cluster labels: stored value = displayed label",
            rows = 8,
            
            value = paste(
              "ChEL = ChEL",
              "TRL = TRL",
              "Lp(a) = Lp(a)",
              "MUNW = MUNW",
              "Obesity = Obesity",
              "Endothelial = Endothelial dysfunction",
              "Coagulation = Coagulation",
              "Cardiac = Cardiac remodeling",
              sep = "\n"
            )
          ),
          
          tags$hr(),
          
          textAreaInput(
            "dot_coronary_cells",
            "Coronary-plaque cell order (optional; all cells are always shown)",
            rows = 8,
            
            value = paste(
              default_coronary_cells,
              collapse = "\n"
            )
          ),
          
          textAreaInput(
            "dot_heart_cells",
            "Heart-atlas cell order (optional; all cells are always shown)",
            rows = 12,
            
            value = paste(
              default_heart_cells,
              collapse = "\n"
            )
          ),
          
          textAreaInput(
            "dot_catlas_cells",
            "Whole-body CATLAS cell order (optional; qualifying cells not listed are appended)",
            rows = 18,
            
            value = paste(
              default_catlas_cells,
              collapse = "\n"
            )
          ),
          
          textAreaInput(
            "dot_cell_labels",
            "Cell-label overrides: cell_type = label (or atlas|cell_type = label)",
            rows = 5,
            
            value = default_cell_label_overrides
          ),
          
          checkboxInput(
            "dot_use_cell_atlas_name",
            "Use dataset description as default cell label",
            TRUE
          ),
          helpText(
            "For the revised dataset, description is used as the default display label. Custom labels override it."
          )
        ),
        
        tabPanel(
          "Dot appearance",

          h4("Whole-body CATLAS filtering"),

          selectInput(
            "dot_catlas_filter_mode",
            "Cells to show from whole-body CATLAS",
            choices = c(
              "Nominal P < alpha in at least one process" = "nominal",
              "Bonferroni P < alpha / N in at least one process" = "bonferroni"
            ),
            selected = "bonferroni"
          ),

          numericInput(
            "dot_bonferroni_n_tests",
            "Bonferroni denominator (N tests)",
            value = 128,
            min = 1,
            step = 1
          ),

          helpText(
            "Plaque and heart atlas always show every available cell. The whole-body filter retains a cell if at least one process meets the selected threshold; all process dots are then displayed for that cell."
          ),

          textOutput("dot_bonferroni_status"),

          tags$hr(),
          h4("Point significance and appearance"),

          numericInput(
            "dot_nominal_alpha",
            "Nominal alpha",
            0.05,
            min = 0,
            max = 1,
            step = 0.001
          ),
          
          tags$hr(),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_size_nonsignificant",
                "P \u2265 0.05 size",
                1.8,
                min = 0,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_size_nominal",
                "Nominal size",
                4,
                min = 0,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_size_bonferroni",
                "Bonferroni size",
                6.5,
                min = 0,
                step = 0.1
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_outline_nonsignificant",
                "Non-significant stroke",
                0.5,
                min = 0,
                step = 0.05
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_outline_significant",
                "Significant stroke",
                1,
                min = 0,
                step = 0.05
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              textInput(
                "dot_outline_nonsignificant_colour",
                "Non-significant outline",
                "grey25"
              )
            ),
            
            column(
              6,
              
              textInput(
                "dot_outline_significant_colour",
                "Significant outline",
                "grey25"
              )
            )
          ),
          
          textInput(
            "dot_significance_legend_colour",
            "Significance-legend colour",
            "grey60"
          ),
          
          tags$hr(),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_z_min",
                "Z-score minimum",
                -5,
                step = 0.5
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_z_max",
                "Z-score maximum",
                5,
                step = 0.5
              )
            )
          ),
          
          textInput(
            "dot_z_breaks",
            "Legend breaks, comma-separated",
            "-5,-2.5,0,2.5,5"
          ),
          
          fluidRow(
            column(
              4,
              
              textInput(
                "dot_colour_low",
                "Low colour",
                "steelblue4"
              )
            ),
            
            column(
              4,
              
              textInput(
                "dot_colour_mid",
                "Mid colour",
                "white"
              )
            ),
            
            column(
              4,
              
              textInput(
                "dot_colour_high",
                "High colour",
                "firebrick3"
              )
            )
          ),
          
          textInput(
            "dot_z_legend_title",
            "Z-score legend title",
            "Z-score"
          ),
          
          textInput(
            "dot_significance_legend_title",
            "Significance legend title",
            "Significance"
          ),
          
          tags$hr(),
          
          fluidRow(
            column(
              4,
              
              textInput(
                "dot_tile_fill",
                "Tile fill",
                "grey93"
              )
            ),
            
            column(
              4,
              
              textInput(
                "dot_tile_border",
                "Tile border",
                "white"
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_tile_border_width",
                "Border width",
                0.65,
                min = 0,
                step = 0.05
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_tile_width",
                "Tile width",
                0.92,
                min = 0,
                step = 0.01
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_tile_height",
                "Tile height",
                0.92,
                min = 0,
                step = 0.01
              )
            )
          )
        ),
        
        tabPanel(
          "Dot labels",
          
          h4("X-axis cell labels"),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_x_text_angle",
                "Angle",
                45,
                step = 1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_x_text_size",
                "Size",
                9,
                min = 1,
                step = 0.5
              )
            ),
            
            column(
              4,
              
              selectInput(
                "dot_x_text_face",
                "Face",
                
                choices = c(
                  "plain",
                  "bold",
                  "italic",
                  "bold.italic"
                ),
                
                selected = "plain"
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_x_text_hjust",
                "Horizontal justification",
                1,
                step = 0.1
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_x_text_vjust",
                "Vertical justification",
                1,
                step = 0.1
              )
            )
          ),
          
          textInput(
            "dot_x_text_colour",
            "X-axis text colour",
            "grey20"
          ),
          
          tags$hr(),
          
          h4("Y-axis title"),
          
          textInput(
            "dot_y_axis_title",
            "Title",
            "Cell-type enrichment cluster"
          ),
          
          fluidRow(
            column(
              3,
              
              numericInput(
                "dot_y_axis_title_size",
                "Size",
                11,
                min = 1,
                step = 0.5
              )
            ),
            
            column(
              3,
              
              selectInput(
                "dot_y_axis_title_face",
                "Face",
                
                choices = c(
                  "plain",
                  "bold",
                  "italic",
                  "bold.italic"
                ),
                
                selected = "bold"
              )
            ),
            
            column(
              3,
              
              textInput(
                "dot_y_axis_title_colour",
                "Colour",
                "grey20"
              )
            ),
            
            column(
              3,
              
              numericInput(
                "dot_y_axis_title_angle",
                "Angle",
                90,
                step = 1
              )
            )
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_y_axis_title_hjust",
                "Vertical position",
                0.5,
                step = 0.05
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_y_axis_title_vjust",
                "Horizontal position",
                1,
                step = 0.05
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_y_axis_title_margin",
                "Right margin",
                10,
                min = 0,
                step = 1
              )
            )
          ),
          
          helpText(
            paste(
              "For a 90-degree title, horizontal justification",
              "moves the title vertically along the y-axis."
            )
          ),
          
          tags$hr(),
          
          h4("Y-axis cluster labels"),
          
          checkboxInput(
            "dot_wrap_y_labels",
            "Wrap y-axis labels",
            TRUE
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_y_label_wrap_width",
                "Wrap width",
                16,
                min = 1,
                step = 1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_y_text_size",
                "Size",
                10,
                min = 1,
                step = 0.5
              )
            ),
            
            column(
              4,
              
              selectInput(
                "dot_y_text_face",
                "Face",
                
                choices = c(
                  "plain",
                  "bold",
                  "italic",
                  "bold.italic"
                ),
                
                selected = "plain"
              )
            )
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_y_text_hjust",
                "Horizontal justification",
                1,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_y_text_vjust",
                "Vertical justification",
                0.5,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_y_text_lineheight",
                "Wrapped-line height",
                0.9,
                min = 0.1,
                step = 0.05
              )
            )
          ),
          
          textInput(
            "dot_y_text_colour",
            "Y-axis text colour",
            "grey20"
          ),
          
          tags$hr(),
          
          h4("Atlas headers"),
          
          checkboxInput(
            "dot_show_atlas_headers",
            "Show atlas headers",
            TRUE
          ),
          
          checkboxInput(
            "dot_show_atlas_brackets",
            "Show atlas brackets",
            TRUE
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_atlas_header_size",
                "Header size",
                11,
                min = 1,
                step = 0.5
              )
            ),
            
            column(
              4,
              
              selectInput(
                "dot_atlas_header_face",
                "Header face",
                
                choices = c(
                  "plain",
                  "bold",
                  "italic",
                  "bold.italic"
                ),
                
                selected = "bold"
              )
            ),
            
            column(
              4,
              
              textInput(
                "dot_atlas_header_colour",
                "Header colour",
                "grey20"
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_atlas_bracket_offset",
                "Bracket Y offset",
                0.75,
                step = 0.05
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_atlas_label_offset",
                "Header Y offset",
                1.2,
                step = 0.05
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              textInput(
                "dot_atlas_bracket_colour",
                "Bracket colour",
                "grey45"
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_atlas_bracket_width",
                "Bracket width",
                0.5,
                min = 0,
                step = 0.05
              )
            )
          ),
          
          tags$hr(),
          
          h4("Atlas separators"),
          
          checkboxInput(
            "dot_show_atlas_separator",
            "Show dashed separator",
            TRUE
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_atlas_gap_line_width",
                "White-gap width",
                4,
                min = 0,
                step = 0.25
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_atlas_separator_width",
                "Separator width",
                0.45,
                min = 0,
                step = 0.05
              )
            ),
            
            column(
              4,
              
              textInput(
                "dot_atlas_separator_colour",
                "Separator colour",
                "grey65"
              )
            )
          ),
          
          selectInput(
            "dot_atlas_separator_linetype",
            "Separator linetype",
            
            choices = c(
              "solid",
              "dashed",
              "dotted",
              "dotdash",
              "longdash",
              "twodash"
            ),
            
            selected = "dashed"
          )
        ),
        
        tabPanel(
          "Dot plot layout",
          
          numericInput(
            "dot_base_text_size",
            "Base text size",
            11,
            min = 1,
            step = 0.5
          ),
          
          selectInput(
            "dot_legend_position",
            "Legend position",
            
            choices = c(
              "right",
              "left",
              "top",
              "bottom",
              "none"
            ),
            
            selected = "right"
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_legend_title_size",
                "Legend-title size",
                10,
                min = 1,
                step = 0.5
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_legend_text_size",
                "Legend-text size",
                9,
                min = 1,
                step = 0.5
              )
            )
          ),
          
          selectInput(
            "dot_legend_title_face",
            "Legend-title face",
            
            choices = c(
              "plain",
              "bold",
              "italic",
              "bold.italic"
            ),
            
            selected = "bold"
          ),
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "dot_legend_key_height",
                "Legend-key height cm",
                0.7,
                min = 0,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_colorbar_width",
                "Colour-bar width",
                1,
                min = 0,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "dot_colorbar_height",
                "Colour-bar height",
                4,
                min = 0,
                step = 0.25
              )
            )
          ),
          
          tags$hr(),
          
          fluidRow(
            column(
              3,
              
              numericInput(
                "dot_margin_top",
                "Top margin",
                65,
                min = 0,
                step = 1
              )
            ),
            
            column(
              3,
              
              numericInput(
                "dot_margin_right",
                "Right margin",
                15,
                min = 0,
                step = 1
              )
            ),
            
            column(
              3,
              
              numericInput(
                "dot_margin_bottom",
                "Bottom margin",
                15,
                min = 0,
                step = 1
              )
            ),
            
            column(
              3,
              
              numericInput(
                "dot_margin_left",
                "Left margin",
                15,
                min = 0,
                step = 1
              )
            )
          ),
          
          tags$hr(),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_x_expand_left",
                "X expansion left",
                0.01,
                min = 0,
                step = 0.01
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_x_expand_right",
                "X expansion right",
                0.01,
                min = 0,
                step = 0.01
              )
            )
          ),
          
          fluidRow(
            column(
              6,
              
              numericInput(
                "dot_y_expand_bottom",
                "Y expansion bottom",
                0.04,
                min = 0,
                step = 0.01
              )
            ),
            
            column(
              6,
              
              numericInput(
                "dot_y_expand_top",
                "Y expansion top",
                0.16,
                min = 0,
                step = 0.01
              )
            )
          )
        ),
        
        tabPanel(
          "Export and JSON",
          
          fluidRow(
            column(
              4,
              
              numericInput(
                "export_width",
                "Width in",
                18,
                min = 1,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "export_height",
                "Height in",
                14,
                min = 1,
                step = 0.1
              )
            ),
            
            column(
              4,
              
              numericInput(
                "export_dpi",
                "DPI",
                600,
                min = 72,
                step = 1
              )
            )
          ),
          
          downloadButton(
            "download_png",
            "Download combined PNG"
          ),
          
          tags$br(),
          tags$br(),
          
          downloadButton(
            "download_pdf",
            "Download combined PDF"
          ),
          
          tags$hr(),
          
          downloadButton(
            "download_json",
            "Export settings JSON"
          ),
          
          tags$br(),
          tags$br(),
          
          fileInput(
            "json_file",
            "Import settings JSON",
            
            accept = c(
              ".json",
              "application/json"
            )
          ),
          
          actionButton(
            "import_json_btn",
            "Apply imported JSON settings"
          ),
          
          tags$br(),
          tags$br(),
          
          verbatimTextOutput(
            "json_status"
          )
        )
      )
    ),
    
    mainPanel(
      width = 8,
      
      tabsetPanel(
        tabPanel(
          "Combined figure",
          uiOutput("combined_plot_ui")
        ),
        
        tabPanel(
          "Panel A only",
          
          plotOutput(
            "radial_plot",
            height = "850px"
          )
        ),
        
        tabPanel(
          "Panel B only",
          
          plotOutput(
            "dot_plot",
            height = "750px"
          )
        ),
        
        tabPanel(
          "Data preview",
          
          h4("Variant/gene weights"),
          DTOutput("gene_table"),
          
          h4("Trait weights"),
          DTOutput("trait_table"),
          
          h4("Dot-matrix data"),
          DTOutput("dot_table")
        ),
        
        tabPanel(
          "Processed dot results",
          DTOutput("dot_results_table")
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(
    input,
    output,
    session
) {
  propagation_message <- reactiveVal(
    "No radial settings propagated yet."
  )
  
  json_message <- reactiveVal(
    "No JSON settings imported."
  )
  
  output$propagation_status <- renderText({
    propagation_message()
  })
  
  output$json_status <- renderText({
    json_message()
  })
  
  # ---------------------------------------------------------------------------
  # SAFE UPDATE FUNCTIONS
  # ---------------------------------------------------------------------------
  
  update_numeric_input_safe <- function(
    input_id,
    value
  ) {
    numeric_value <- suppressWarnings(
      as.numeric(
        unlist(value)
      )
    )
    
    if (
      length(numeric_value) > 0 &&
      is.finite(numeric_value[1])
    ) {
      updateNumericInput(
        session,
        input_id,
        value = numeric_value[1]
      )
    }
  }
  
  update_checkbox_input_safe <- function(
    input_id,
    value
  ) {
    updateCheckboxInput(
      session,
      input_id,
      value = isTRUE(value)
    )
  }
  
  update_select_input_safe <- function(
    input_id,
    value
  ) {
    updateSelectInput(
      session,
      input_id,
      selected = as.character(value)
    )
  }
  
  update_text_input_safe <- function(
    input_id,
    value
  ) {
    updateTextInput(
      session,
      input_id,
      value = as.character(value %||% "")
    )
  }
  
  update_text_area_safe <- function(
    input_id,
    value
  ) {
    updateTextAreaInput(
      session,
      input_id,
      value = as.character(value %||% "")
    )
  }
  
  update_slider_input_safe <- function(
    input_id,
    value
  ) {
    numeric_value <- suppressWarnings(
      as.numeric(
        unlist(value)
      )
    )
    
    if (
      length(numeric_value) >= 2 &&
      all(
        is.finite(
          numeric_value[1:2]
        )
      )
    ) {
      updateSliderInput(
        session,
        input_id,
        value = numeric_value[1:2]
      )
    }
  }
  
  update_radial_setting <- function(
    cluster,
    suffix,
    value
  ) {
    input_id <- paste0(
      cluster,
      "_",
      suffix
    )
    
    if (suffix %in% radial_numeric_suffixes) {
      update_numeric_input_safe(
        input_id,
        value
      )
    } else if (suffix %in% radial_checkbox_suffixes) {
      update_checkbox_input_safe(
        input_id,
        value
      )
    } else if (suffix %in% radial_select_suffixes) {
      update_select_input_safe(
        input_id,
        value
      )
    } else if (suffix == "trait_label_manual_offsets") {
      update_text_area_safe(
        input_id,
        value
      )
    } else if (suffix %in% radial_text_suffixes) {
      update_text_input_safe(
        input_id,
        value
      )
    } else if (suffix %in% radial_slider_suffixes) {
      update_slider_input_safe(
        input_id,
        value
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # PROPAGATE RADIAL SETTINGS
  # ---------------------------------------------------------------------------
  
  observeEvent(
    input$propagate_cluster_settings,
    {
      source_cluster <- input$propagate_source_cluster
      selected_group <- input$propagate_setting_group
      
      suffixes_to_copy <- radial_setting_groups[[selected_group]]
      
      shiny::validate(
        shiny::need(
          source_cluster %in% radial_clusters,
          "Select a valid source cluster."
        ),
        
        shiny::need(
          length(suffixes_to_copy) > 0,
          "Select at least one settings group."
        )
      )
      
      target_clusters <- setdiff(
        radial_clusters,
        source_cluster
      )
      
      for (target_cluster in target_clusters) {
        for (suffix in suffixes_to_copy) {
          source_input_id <- paste0(
            source_cluster,
            "_",
            suffix
          )
          
          source_value <- input[[source_input_id]]
          
          update_radial_setting(
            cluster = target_cluster,
            suffix = suffix,
            value = source_value
          )
        }
      }
      
      propagation_message(
        paste0(
          "Copied “",
          selected_group,
          "” from ",
          source_cluster,
          " to ",
          paste(
            target_clusters,
            collapse = ", "
          ),
          "."
        )
      )
    }
  )
  
  # ---------------------------------------------------------------------------
  # DATA
  # ---------------------------------------------------------------------------
  
  circular_data <- reactive({
    if (isTRUE(input$use_example_circular)) {
      return(
        generate_example_circular_data()
      )
    }
    
    shiny::validate(
      shiny::need(
        !is.null(input$gene_file),
        "Upload the variant/gene weights CSV."
      ),
      
      shiny::need(
        !is.null(input$trait_file),
        "Upload the trait weights CSV."
      )
    )
    
    list(
      gene = safe_read_csv(
        input$gene_file$datapath
      ),
      
      trait = safe_read_csv(
        input$trait_file$datapath
      )
    )
  })
  
  dot_data <- reactive({
    if (isTRUE(input$use_example_dot)) {
      return(generate_example_dot_data())
    }

    if (!is.null(input$dot_file)) {
      return(safe_read_csv(input$dot_file$datapath))
    }

    local_dot_file <- "gchromvar_ppa_matrix_long.csv"
    if (file.exists(local_dot_file)) {
      return(safe_read_csv(local_dot_file))
    }

    shiny::validate(
      shiny::need(
        FALSE,
        paste(
          "Upload the revised dot-matrix CSV, or place",
          "gchromvar_ppa_matrix_long.csv in the same working directory as the app."
        )
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # RADIAL REACTIVES
  # ---------------------------------------------------------------------------
  
  cluster_order_reactive <- reactive({
    value <- parse_comma_vector(
      input$cluster_order
    )
    
    if (length(value) == 0) {
      default_radial_cluster_order
    } else {
      value
    }
  })
  
  gene_labels_reactive <- reactive({
    parse_named_mapping(
      input$gene_labels_text
    )
  })
  
  trait_labels_reactive <- reactive({
    parse_named_mapping(
      input$trait_labels_text
    )
  })
  
  variant_n_reactive <- reactive({
    parse_numeric_mapping(
      input$variant_n_text
    )
  })
  
  title_positions_reactive <- reactive({
    title_mapping <- parse_named_mapping(
      input$title_mapping_text
    )
    
    title_wrap_mapping <- parse_numeric_mapping(
      input$title_wrap_width_mapping_text
    )
    
    x_mapping <- parse_numeric_mapping(
      input$title_x_mapping_text
    )
    
    y_mapping <- parse_numeric_mapping(
      input$title_y_mapping_text
    )
    
    raw_titles <- vapply(
      radial_clusters,
      
      function(cluster_name) {
        title_mapping[[cluster_name]] %||%
          cluster_name
      },
      
      character(1)
    )
    
    wrap_widths <- vapply(
      radial_clusters,
      
      function(cluster_name) {
        title_wrap_mapping[[cluster_name]] %||%
          input$title_wrap_width_default
      },
      
      numeric(1)
    )
    
    wrapped_titles <- vapply(
      seq_along(radial_clusters),
      
      function(i) {
        wrap_text_width(
          raw_titles[i],
          wrap_widths[i]
        )
      },
      
      character(1)
    )
    
    x_positions <- vapply(
      radial_clusters,
      
      function(cluster_name) {
        x_mapping[[cluster_name]] %||%
          input$title_x_default
      },
      
      numeric(1)
    )
    
    y_positions <- vapply(
      radial_clusters,
      
      function(cluster_name) {
        y_mapping[[cluster_name]] %||%
          input$title_y_default
      },
      
      numeric(1)
    )
    
    data.frame(
      Cluster = radial_clusters,
      Title = wrapped_titles,
      WrapWidth = wrap_widths,
      X = x_positions,
      Y = y_positions,
      stringsAsFactors = FALSE
    )
  })
  
  cluster_settings_reactive <- reactive({
    setNames(
      lapply(radial_clusters, function(cluster_name) {
        get_cluster_settings_from_input(
          input,
          cluster_name
        )
      }),
      
      radial_clusters
    )
  })
  
  radial_panel_grob <- reactive({
    data_object <- circular_data()
    cluster_order <- cluster_order_reactive()
    settings_list <- cluster_settings_reactive()
    
    shiny::validate(
      shiny::need(
        "gene" %in% names(data_object$gene),
        "The radial gene-weight data need a column named gene."
      ),
      
      shiny::need(
        "trait" %in% names(data_object$trait),
        "The radial trait-weight data need a column named trait."
      )
    )
    
    plot_items <- lapply(cluster_order, function(cluster_name) {
      if (!(cluster_name %in% names(settings_list))) {
        return(NULL)
      }
      
      plot_object <- create_cluster_plot(
        cluster = cluster_name,
        gene_df = data_object$gene,
        trait_df = data_object$trait,
        settings = settings_list[[cluster_name]],
        trait_labels = trait_labels_reactive(),
        gene_labels = gene_labels_reactive(),
        variant_n_map = variant_n_reactive(),
        for_panel = TRUE
      )
      
      if (is.null(plot_object)) {
        return(NULL)
      }
      
      list(
        cluster = cluster_name,
        grob = ggplotGrob(plot_object)
      )
    })
    
    plot_items <- Filter(
      Negate(is.null),
      plot_items
    )
    
    shiny::validate(
      shiny::need(
        length(plot_items) > 0,
        
        paste(
          "No radial plots could be created.",
          "Check the radial cluster columns and filtering thresholds."
        )
      )
    )
    
    valid_clusters <- vapply(
      plot_items,
      function(item) {
        item$cluster
      },
      character(1)
    )
    
    valid_grobs <- lapply(
      plot_items,
      function(item) {
        item$grob
      }
    )
    
    create_panel_with_titles(
      plots = valid_grobs,
      cluster_order = valid_clusters,
      cluster_settings = settings_list,
      ncol = input$panel_ncol,
      nrow = input$panel_nrow,
      h_spacing = input$panel_h_spacing,
      v_spacing = input$panel_v_spacing,
      plot_scale = input$plot_scale,
      title_positions = title_positions_reactive(),
      show_titles = input$show_titles,
      title_size = input$title_size,
      title_fontface = input$title_fontface,
      title_color = input$title_color,
      bg_color = input$panel_background
    )
  })
  
  # ---------------------------------------------------------------------------
  # DOT SETTINGS
  # ---------------------------------------------------------------------------
  
  dot_settings_reactive <- reactive({
    z_breaks <- suppressWarnings(
      as.numeric(
        parse_comma_vector(
          input$dot_z_breaks
        )
      )
    )
    
    z_breaks <- z_breaks[
      is.finite(z_breaks)
    ]
    
    if (length(z_breaks) == 0) {
      z_breaks <- c(
        -5,
        -2.5,
        0,
        2.5,
        5
      )
    }
    
    atlas_order <- parse_comma_vector(
      input$dot_atlas_order
    )
    
    if (length(atlas_order) == 0) {
      atlas_order <- default_dot_atlas_order
    }
    
    cluster_order <- parse_comma_vector(
      input$dot_cluster_order
    )
    
    if (length(cluster_order) == 0) {
      cluster_order <- default_dot_cluster_order
    }
    
    z_min <- input$dot_z_min
    z_max <- input$dot_z_max
    
    if (
      !is.finite(z_min) ||
      !is.finite(z_max) ||
      z_min >= z_max
    ) {
      z_min <- -5
      z_max <- 5
    }
    
    list(
      nominal_alpha = input$dot_nominal_alpha,
      catlas_filter_mode = input$dot_catlas_filter_mode %||% "bonferroni",
      bonferroni_n_tests = max(1, as.numeric(input$dot_bonferroni_n_tests %||% 128)),
      
      atlas_order = atlas_order,
      
      atlas_labels = parse_named_mapping(
        input$dot_atlas_labels
      ),
      
      atlas_gap = input$dot_atlas_gap,
      cluster_order = cluster_order,
      
      cluster_labels = parse_named_mapping(
        input$dot_cluster_labels
      ),
      
      cell_order = list(
        "Coronary plaque" = parse_line_vector(input$dot_coronary_cells),
        "Heart atlas" = parse_line_vector(input$dot_heart_cells),
        "Whole-body (CATLAS)" = parse_line_vector(input$dot_catlas_cells)
      ),
      
      cell_labels = parse_named_mapping(
        input$dot_cell_labels
      ),
      
      use_cell_atlas_name_as_label = input$dot_use_cell_atlas_name,
      
      dot_sizes = c(
        "P \u2265 0.05" = input$dot_size_nonsignificant,
        "Nominal P < 0.05" = input$dot_size_nominal,
        "Bonferroni P < 0.05" = input$dot_size_bonferroni
      ),
      
      nonsignificant_outline_colour = input$dot_outline_nonsignificant_colour,
      nonsignificant_outline_width = input$dot_outline_nonsignificant,
      significant_outline_colour = input$dot_outline_significant_colour,
      significant_outline_width = input$dot_outline_significant,
      significance_legend_colour = input$dot_significance_legend_colour,
      
      z_min = z_min,
      z_max = z_max,
      z_breaks = z_breaks,
      
      divergent_colour_low = input$dot_colour_low,
      divergent_colour_mid = input$dot_colour_mid,
      divergent_colour_high = input$dot_colour_high,
      
      z_legend_title = input$dot_z_legend_title,
      significance_legend_title = input$dot_significance_legend_title,
      
      tile_fill = input$dot_tile_fill,
      tile_border_colour = input$dot_tile_border,
      tile_border_width = input$dot_tile_border_width,
      tile_width = input$dot_tile_width,
      tile_height = input$dot_tile_height,
      
      x_text_angle = input$dot_x_text_angle,
      x_text_size = input$dot_x_text_size,
      x_text_face = input$dot_x_text_face,
      x_text_hjust = input$dot_x_text_hjust,
      x_text_vjust = input$dot_x_text_vjust,
      x_text_colour = input$dot_x_text_colour,
      
      y_axis_title = input$dot_y_axis_title,
      y_axis_title_size = input$dot_y_axis_title_size,
      y_axis_title_face = input$dot_y_axis_title_face,
      y_axis_title_colour = input$dot_y_axis_title_colour,
      y_axis_title_angle = input$dot_y_axis_title_angle,
      y_axis_title_hjust = input$dot_y_axis_title_hjust,
      y_axis_title_vjust = input$dot_y_axis_title_vjust,
      y_axis_title_margin = input$dot_y_axis_title_margin,
      
      wrap_y_labels = input$dot_wrap_y_labels,
      y_label_wrap_width = input$dot_y_label_wrap_width,
      y_text_size = input$dot_y_text_size,
      y_text_face = input$dot_y_text_face,
      y_text_colour = input$dot_y_text_colour,
      y_text_hjust = input$dot_y_text_hjust,
      y_text_vjust = input$dot_y_text_vjust,
      y_text_lineheight = input$dot_y_text_lineheight,
      
      show_atlas_headers = input$dot_show_atlas_headers,
      show_atlas_brackets = input$dot_show_atlas_brackets,
      atlas_header_size = input$dot_atlas_header_size,
      atlas_header_fontface = input$dot_atlas_header_face,
      atlas_header_colour = input$dot_atlas_header_colour,
      atlas_bracket_offset = input$dot_atlas_bracket_offset,
      atlas_label_offset = input$dot_atlas_label_offset,
      atlas_bracket_colour = input$dot_atlas_bracket_colour,
      atlas_bracket_width = input$dot_atlas_bracket_width,
      
      show_atlas_separator = input$dot_show_atlas_separator,
      atlas_gap_line_width = input$dot_atlas_gap_line_width,
      atlas_separator_width = input$dot_atlas_separator_width,
      atlas_separator_colour = input$dot_atlas_separator_colour,
      atlas_separator_linetype = input$dot_atlas_separator_linetype,
      
      base_text_size = input$dot_base_text_size,
      legend_position = input$dot_legend_position,
      legend_title_size = input$dot_legend_title_size,
      legend_title_face = input$dot_legend_title_face,
      legend_text_size = input$dot_legend_text_size,
      legend_key_height = input$dot_legend_key_height,
      colorbar_width = input$dot_colorbar_width,
      colorbar_height = input$dot_colorbar_height,
      
      plot_margin_top = input$dot_margin_top,
      plot_margin_right = input$dot_margin_right,
      plot_margin_bottom = input$dot_margin_bottom,
      plot_margin_left = input$dot_margin_left,
      
      x_expand_left = input$dot_x_expand_left,
      x_expand_right = input$dot_x_expand_right,
      y_expand_bottom = input$dot_y_expand_bottom,
      y_expand_top = input$dot_y_expand_top
    )
  })
  
  prepared_dot_data <- reactive({
    prepare_dot_plot_data(
      source_df = dot_data(),
      settings = dot_settings_reactive()
    )
  })

  output$dot_bonferroni_status <- renderText({
    settings <- dot_settings_reactive()
    cutoff <- settings$nominal_alpha / settings$bonferroni_n_tests
    paste0(
      "Current Bonferroni raw-P cutoff: ",
      format(cutoff, digits = 4, scientific = TRUE),
      " (", settings$nominal_alpha, " / ", settings$bonferroni_n_tests, ")."
    )
  })
  
  dot_plot_reactive <- reactive({
    create_coloured_dot_heatmap(
      prepared = prepared_dot_data(),
      settings = dot_settings_reactive()
    )
  })
  
  dot_grob_reactive <- reactive({
    ggplotGrob(
      dot_plot_reactive()
    )
  })
  
  # ---------------------------------------------------------------------------
  # COMBINED FIGURE
  # ---------------------------------------------------------------------------
  
  combined_grob_reactive <- reactive({
    create_combined_figure(
      radial_grob = radial_panel_grob(),
      dot_grob = dot_grob_reactive(),
      top_height = input$top_panel_height,
      bottom_height = input$bottom_panel_height,
      panel_a_label = input$panel_a_label,
      panel_a_x = input$panel_a_label_x,
      panel_a_y = input$panel_a_label_y,
      panel_a_size = input$panel_a_label_size,
      panel_a_colour = input$panel_a_label_colour,
      panel_b_label = input$panel_b_label,
      panel_b_x = input$panel_b_label_x,
      panel_b_y = input$panel_b_label_y,
      panel_b_size = input$panel_b_label_size,
      panel_b_colour = input$panel_b_label_colour
    )
  })
  
  # ---------------------------------------------------------------------------
  # PREVIEWS
  # ---------------------------------------------------------------------------
  
  output$combined_plot_ui <- renderUI({
    plotOutput(
      "combined_plot",
      
      height = paste0(
        input$preview_height_px,
        "px"
      )
    )
  })
  
  output$combined_plot <- renderPlot({
    grid.newpage()
    
    grid.draw(
      combined_grob_reactive()
    )
  }, bg = "white")
  
  output$radial_plot <- renderPlot({
    grid.newpage()
    
    grid.draw(
      radial_panel_grob()
    )
  }, bg = "white")
  
  output$dot_plot <- renderPlot({
    print(
      dot_plot_reactive()
    )
  }, bg = "white")
  
  # ---------------------------------------------------------------------------
  # TABLES
  # ---------------------------------------------------------------------------
  
  output$gene_table <- renderDT({
    datatable(
      head(
        circular_data()$gene,
        100
      ),
      
      options = list(
        scrollX = TRUE,
        pageLength = 5
      )
    )
  })
  
  output$trait_table <- renderDT({
    datatable(
      head(
        circular_data()$trait,
        100
      ),
      
      options = list(
        scrollX = TRUE,
        pageLength = 5
      )
    )
  })
  
  output$dot_table <- renderDT({
    datatable(
      head(
        dot_data(),
        200
      ),
      
      options = list(
        scrollX = TRUE,
        pageLength = 8
      )
    )
  })
  
  output$dot_results_table <- renderDT({
    results <- prepared_dot_data()$plot_df %>%
      select(
        atlas,
        cell,
        cell_atlas_name,
        cluster,
        z,
        p_nominal,
        q_plot,
        n_selected_tests_per_cluster,
        significance_plot,
        x_pos,
        y_pos
      ) %>%
      arrange(
        x_pos,
        desc(y_pos)
      )
    
    datatable(
      results,
      
      options = list(
        scrollX = TRUE,
        pageLength = 15
      )
    )
  })
  
  # ---------------------------------------------------------------------------
  # FIGURE EXPORT
  # ---------------------------------------------------------------------------
  
  output$download_png <- downloadHandler(
    filename = function() {
      paste0(
        "integrated_radial_dot_figure_",
        Sys.Date(),
        ".png"
      )
    },
    
    content = function(file) {
      png(
        filename = file,
        width = input$export_width,
        height = input$export_height,
        units = "in",
        res = input$export_dpi,
        bg = "white"
      )
      
      on.exit(
        dev.off(),
        add = TRUE
      )
      
      grid.newpage()
      
      grid.draw(
        combined_grob_reactive()
      )
    }
  )
  
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0(
        "integrated_radial_dot_figure_",
        Sys.Date(),
        ".pdf"
      )
    },
    
    content = function(file) {
      pdf(
        file = file,
        width = input$export_width,
        height = input$export_height,
        bg = "white",
        useDingbats = FALSE
      )
      
      on.exit(
        dev.off(),
        add = TRUE
      )
      
      grid.newpage()
      
      grid.draw(
        combined_grob_reactive()
      )
    }
  )
  
  # ---------------------------------------------------------------------------
  # JSON SETTINGS
  # ---------------------------------------------------------------------------
  
  general_config_ids <- c(
    "variant_n_text",
    "gene_labels_text",
    "trait_labels_text",
    
    "top_panel_height",
    "bottom_panel_height",
    
    "panel_a_label",
    "panel_a_label_x",
    "panel_a_label_y",
    "panel_a_label_size",
    "panel_a_label_colour",
    
    "panel_b_label",
    "panel_b_label_x",
    "panel_b_label_y",
    "panel_b_label_size",
    "panel_b_label_colour",
    
    "preview_height_px",
    
    "cluster_order",
    "panel_ncol",
    "panel_nrow",
    "panel_h_spacing",
    "panel_v_spacing",
    "plot_scale",
    "panel_background",
    
    "show_titles",
    "title_size",
    "title_fontface",
    "title_color",
    "title_mapping_text",
    "title_wrap_width_default",
    "title_wrap_width_mapping_text",
    "title_x_default",
    "title_y_default",
    "title_x_mapping_text",
    "title_y_mapping_text",
    
    "propagate_source_cluster",
    "propagate_setting_group",
    
    "dot_atlas_order",
    "dot_atlas_labels",
    "dot_atlas_gap",
    
    "dot_cluster_order",
    "dot_cluster_labels",
    
    "dot_coronary_cells",
    "dot_heart_cells",
    "dot_catlas_cells",
    "dot_cell_labels",
    "dot_use_cell_atlas_name",
    "dot_catlas_filter_mode",
    "dot_bonferroni_n_tests",
    
    "dot_nominal_alpha",
    
    "dot_size_nonsignificant",
    "dot_size_nominal",
    "dot_size_bonferroni",
    
    "dot_outline_nonsignificant",
    "dot_outline_significant",
    "dot_outline_nonsignificant_colour",
    "dot_outline_significant_colour",
    "dot_significance_legend_colour",
    
    "dot_z_min",
    "dot_z_max",
    "dot_z_breaks",
    "dot_colour_low",
    "dot_colour_mid",
    "dot_colour_high",
    "dot_z_legend_title",
    "dot_significance_legend_title",
    
    "dot_tile_fill",
    "dot_tile_border",
    "dot_tile_border_width",
    "dot_tile_width",
    "dot_tile_height",
    
    "dot_x_text_angle",
    "dot_x_text_size",
    "dot_x_text_face",
    "dot_x_text_hjust",
    "dot_x_text_vjust",
    "dot_x_text_colour",
    
    "dot_y_axis_title",
    "dot_y_axis_title_size",
    "dot_y_axis_title_face",
    "dot_y_axis_title_colour",
    "dot_y_axis_title_angle",
    "dot_y_axis_title_hjust",
    "dot_y_axis_title_vjust",
    "dot_y_axis_title_margin",
    
    "dot_wrap_y_labels",
    "dot_y_label_wrap_width",
    "dot_y_text_size",
    "dot_y_text_face",
    "dot_y_text_colour",
    "dot_y_text_hjust",
    "dot_y_text_vjust",
    "dot_y_text_lineheight",
    
    "dot_show_atlas_headers",
    "dot_show_atlas_brackets",
    "dot_atlas_header_size",
    "dot_atlas_header_face",
    "dot_atlas_header_colour",
    "dot_atlas_bracket_offset",
    "dot_atlas_label_offset",
    "dot_atlas_bracket_colour",
    "dot_atlas_bracket_width",
    
    "dot_show_atlas_separator",
    "dot_atlas_gap_line_width",
    "dot_atlas_separator_width",
    "dot_atlas_separator_colour",
    "dot_atlas_separator_linetype",
    
    "dot_base_text_size",
    "dot_legend_position",
    "dot_legend_title_size",
    "dot_legend_text_size",
    "dot_legend_title_face",
    "dot_legend_key_height",
    "dot_colorbar_width",
    "dot_colorbar_height",
    
    "dot_margin_top",
    "dot_margin_right",
    "dot_margin_bottom",
    "dot_margin_left",
    
    "dot_x_expand_left",
    "dot_x_expand_right",
    "dot_y_expand_bottom",
    "dot_y_expand_top",
    
    "export_width",
    "export_height",
    "export_dpi"
  )
  
  collect_config <- reactive({
    general_settings <- setNames(
      lapply(general_config_ids, function(input_id) {
        input[[input_id]]
      }),
      
      general_config_ids
    )
    
    cluster_settings <- setNames(
      lapply(radial_clusters, function(cluster_name) {
        setNames(
          lapply(all_radial_setting_suffixes, function(suffix) {
            input[[paste0(cluster_name, "_", suffix)]]
          }),
          
          all_radial_setting_suffixes
        )
      }),
      
      radial_clusters
    )
    
    list(
      format_version = 6,
      generated = as.character(Sys.time()),
      general = general_settings,
      radial_clusters = cluster_settings
    )
  })
  
  output$download_json <- downloadHandler(
    filename = function() {
      paste0(
        "integrated_figure_settings_",
        Sys.Date(),
        ".json"
      )
    },
    
    content = function(file) {
      jsonlite::write_json(
        collect_config(),
        path = file,
        pretty = TRUE,
        auto_unbox = TRUE,
        null = "null",
        na = "null"
      )
    }
  )
  
  # Apply the supplied settings JSON automatically when a new session starts.
  # New controls introduced in format v6 keep their UI defaults if absent from
  # the saved v3 configuration.
  session$onFlushed(function() {
    general_settings <- default_app_config$general %||% list()

    for (input_id in names(general_settings)) {
      value <- general_settings[[input_id]]

      if (input_id %in% text_area_ids) {
        update_text_area_safe(input_id, value)
      } else if (input_id %in% checkbox_ids) {
        update_checkbox_input_safe(input_id, value)
      } else if (input_id %in% select_ids) {
        update_select_input_safe(input_id, value)
      } else if (input_id %in% text_ids) {
        update_text_input_safe(input_id, value)
      } else {
        update_numeric_input_safe(input_id, value)
      }
    }

    cluster_config <- default_app_config$radial_clusters %||% list()
    valid_clusters <- intersect(names(cluster_config), radial_clusters)

    for (cluster_name in valid_clusters) {
      cluster_values <- cluster_config[[cluster_name]]
      for (suffix in names(cluster_values)) {
        update_radial_setting(
          cluster = cluster_name,
          suffix = suffix,
          value = cluster_values[[suffix]]
        )
      }
    }

    json_message("Loaded embedded default settings from settings_072626.json.")
  }, once = TRUE)

  imported_config <- eventReactive(
    input$import_json_btn,
    {
      shiny::validate(
        shiny::need(
          !is.null(input$json_file),
          "Select a JSON file first."
        )
      )
      
      jsonlite::read_json(
        input$json_file$datapath,
        simplifyVector = FALSE
      )
    },
    
    ignoreInit = TRUE
  )
  
  text_area_ids <- c(
    "variant_n_text",
    "gene_labels_text",
    "trait_labels_text",
    "title_mapping_text",
    "title_wrap_width_mapping_text",
    "title_x_mapping_text",
    "title_y_mapping_text",
    "dot_atlas_labels",
    "dot_cluster_labels",
    "dot_coronary_cells",
    "dot_heart_cells",
    "dot_catlas_cells",
    "dot_cell_labels"
  )
  
  checkbox_ids <- c(
    "show_titles",
    "dot_use_cell_atlas_name",
    "dot_wrap_y_labels",
    "dot_show_atlas_headers",
    "dot_show_atlas_brackets",
    "dot_show_atlas_separator"
  )
  
  select_ids <- c(
    "panel_background",
    "title_fontface",
    "propagate_source_cluster",
    "propagate_setting_group",
    "dot_catlas_filter_mode",
    "dot_x_text_face",
    "dot_y_axis_title_face",
    "dot_y_text_face",
    "dot_atlas_header_face",
    "dot_atlas_separator_linetype",
    "dot_legend_position",
    "dot_legend_title_face"
  )
  
  text_ids <- c(
    "panel_a_label",
    "panel_a_label_colour",
    "panel_b_label",
    "panel_b_label_colour",
    "cluster_order",
    "title_color",
    "dot_atlas_order",
    "dot_cluster_order",
    "dot_outline_nonsignificant_colour",
    "dot_outline_significant_colour",
    "dot_significance_legend_colour",
    "dot_z_breaks",
    "dot_colour_low",
    "dot_colour_mid",
    "dot_colour_high",
    "dot_z_legend_title",
    "dot_significance_legend_title",
    "dot_tile_fill",
    "dot_tile_border",
    "dot_x_text_colour",
    "dot_y_axis_title",
    "dot_y_axis_title_colour",
    "dot_y_text_colour",
    "dot_atlas_header_colour",
    "dot_atlas_bracket_colour",
    "dot_atlas_separator_colour"
  )
  
  observeEvent(
    imported_config(),
    {
      config <- imported_config()
      
      tryCatch(
        {
          general_settings <- config$general %||% list()
          
          for (input_id in names(general_settings)) {
            value <- general_settings[[input_id]]
            
            if (input_id %in% text_area_ids) {
              update_text_area_safe(
                input_id,
                value
              )
            } else if (input_id %in% checkbox_ids) {
              update_checkbox_input_safe(
                input_id,
                value
              )
            } else if (input_id %in% select_ids) {
              update_select_input_safe(
                input_id,
                value
              )
            } else if (input_id %in% text_ids) {
              update_text_input_safe(
                input_id,
                value
              )
            } else {
              update_numeric_input_safe(
                input_id,
                value
              )
            }
          }
          
          cluster_config <- config$radial_clusters %||% list()
          
          valid_clusters <- intersect(
            names(cluster_config),
            radial_clusters
          )
          
          for (cluster_name in valid_clusters) {
            cluster_values <- cluster_config[[cluster_name]]
            
            for (suffix in names(cluster_values)) {
              update_radial_setting(
                cluster = cluster_name,
                suffix = suffix,
                value = cluster_values[[suffix]]
              )
            }
          }
          
          json_message(
            paste0(
              "Imported settings from: ",
              input$json_file$name
            )
          )
        },
        
        error = function(error_object) {
          json_message(
            paste0(
              "JSON import failed: ",
              error_object$message
            )
          )
        }
      )
    }
  )
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(
  ui = ui,
  server = server
)
