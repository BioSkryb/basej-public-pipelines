suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
# ggplot2 before ggtree so ggplot2 utilities (waiver, themes) resolve for ggtree/ggfun layers.
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(ggtree))
suppressMessages(library(ape))
suppressMessages(library(paletteer))
suppressMessages(library(reshape2))
suppressMessages(library(tidytree))
suppressMessages(library(egg))

# read.table() fails on placement TSVs when a row has extra/missing tabs (e.g. indel REF/ALT
# or ragged rbind of unfiltered + phylo blocks). fread(fill=TRUE) tolerates variable field counts.
read_placement_tsv <- function(path) {
  dt <- data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    na.strings = c("", "NA"),
    data.table = FALSE,
    stringsAsFactors = FALSE
  )
  req <- c("Chr", "Pos", "Ref", "Alt", "Branch")
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop("placement TSV missing columns: ", paste(miss, collapse = ", "))
  }
  ok <- stats::complete.cases(dt[, req, drop = FALSE])
  dt <- dt[ok, , drop = FALSE]
  dt$Chr <- as.character(dt$Chr)
  dt$Pos <- as.character(dt$Pos)
  dt$Ref <- as.character(dt$Ref)
  dt$Alt <- as.character(dt$Alt)
  dt
}

merge_vep_relevance_into_vaf <- function(df_vaf, file_vep) {
  if (missing(file_vep)) return(df_vaf)
  fv <- as.character(file_vep)
  if (!length(fv) || is.na(fv[1]) || !nzchar(fv[1])) return(df_vaf)
  path <- fv[1]
  if (!file.exists(path)) return(df_vaf)
  sz <- file.info(path)$size
  if (is.na(sz) || sz < 2L) return(df_vaf)
  dt <- tryCatch(
    data.table::fread(
      path, sep = "\t", header = TRUE, quote = "",
      na.strings = c("", "NA"), data.table = FALSE, stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(dt) || !nrow(dt) || !"VARIANT_ID" %in% names(dt)) return(df_vaf)
  want <- intersect(c(
    "Verdict", "Genes_affected", "NumberGenesAffected", "Existing_variation_combined",
    "Pathogenicity_evidence", "Consequence_class", "CLIN_SIG_summary", "IMPACT", "Consequence"
  ), names(dt))
  if (!length(want)) return(df_vaf)
  sub <- dt[, c("VARIANT_ID", want), drop = FALSE]
  names(sub)[1] <- "VariantId"
  merge(df_vaf, sub, by = "VariantId", all.x = TRUE)
}

read_vep_verdict_table <- function(file_vep) {
  fv <- as.character(file_vep)
  if (!length(fv) || is.na(fv[1]) || !nzchar(fv[1])) {
    return(data.frame(VARIANT_ID = character(0), Verdict = character(0)))
  }
  path <- fv[1]
  if (!file.exists(path) || file.info(path)$size < 2L) {
    return(data.frame(VARIANT_ID = character(0), Verdict = character(0)))
  }
  dt <- tryCatch(
    data.table::fread(
      path, sep = "\t", header = TRUE, quote = "",
      na.strings = c("", "NA"), data.table = FALSE, stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(dt) || !nrow(dt) || !"VARIANT_ID" %in% names(dt)) {
    return(data.frame(VARIANT_ID = character(0), Verdict = character(0)))
  }
  if (!"Verdict" %in% names(dt)) {
    return(data.frame(VARIANT_ID = character(0), Verdict = character(0)))
  }
  unique(dt[, c("VARIANT_ID", "Verdict"), drop = FALSE])
}

normalize_placed_variants <- function(df_placed) {
  df_placed$VariantId <- paste(df_placed$Chr, df_placed$Pos, df_placed$Ref, df_placed$Alt, sep = "_")
  if ("provenance" %in% names(df_placed)) {
    prov_ord <- c("phylogeny_filtered_variant_placement", "unfiltered_variant_placement")
    rk <- match(df_placed$provenance, prov_ord)
    rk[is.na(rk)] <- length(prov_ord) + 1L
    df_placed <- df_placed[order(rk), , drop = FALSE]
  }
  df_placed <- df_placed[!duplicated(df_placed$VariantId), , drop = FALSE]
  df_placed
}

# Map ape node id (Branch) to all tip.label values in the clade below that node; sorted, "; "-joined.
append_descendant_tip_labels <- function(df_placed, phy) {
  n <- nrow(df_placed)
  if (!n || !"Branch" %in% names(df_placed)) {
    df_placed$DescendantTips <- rep(NA_character_, n)
    return(df_placed)
  }
  br <- suppressWarnings(as.integer(df_placed$Branch))
  u <- unique(br[!is.na(br)])
  if (!length(u)) {
    df_placed$DescendantTips <- rep(NA_character_, n)
    return(df_placed)
  }
  max_node <- max(phy$edge)
  ntip <- ape::Ntip(phy)
  lab_map <- vector("character", length(u))
  names(lab_map) <- as.character(u)
  for (node in u) {
    key <- as.character(node)
    if (is.na(node) || node < 1L || node > max_node) {
      lab_map[key] <- NA_character_
      next
    }
    # Tips are numbered 1..ntip; extract.clade() only accepts internal nodes (node > ntip).
    if (node <= ntip) {
      tl <- phy$tip.label[node]
      lab_map[key] <- if (is.na(tl) || !nzchar(tl)) NA_character_ else as.character(tl)
      next
    }
    sub <- tryCatch(ape::extract.clade(phy, node), error = function(e) NULL)
    if (is.null(sub) || !length(sub$tip.label)) {
      lab_map[key] <- NA_character_
    } else {
      lab_map[key] <- paste(sort(sub$tip.label), collapse = "; ")
    }
  }
  df_placed$DescendantTips <- unname(lab_map[as.character(br)])
  df_placed
}

# log10(p_else_where); p==0 or non-finite -> -20; raster scale squish keeps fill within [-20, 0]
log10_p_else_capped <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  y <- log10(p)
  y[p <= 0 | !is.finite(p) | is.na(p)] <- -20
  y[y < -20] <- -20
  y[y > 0] <- 0
  y
}

theme_ohchibi <- function(size_axis_text.x = 12, size_axis_text.y = 12, size_axis_title.x = 13,
                         size_axis_title.y = 13, angle_text.x = 0, angle_text.y = 0,
                         legend_proportion_size = 1, size_title_text = 13, size_legend_text = 12,
                         size_lines_panel = 0.3, size_panel_border = 0.75, x_hjust = 0.5,
                         y_vjust = 0.5, font_family = "Helvetica", font_face = "plain",
                         size_ticks = 0.5) {
  theme(
    strip.background.x = element_blank(), strip.background.y = element_blank(),
    strip.text.x = element_text(size = size_axis_title.x),
    strip.text.y = element_text(size = size_axis_title.y),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "grey89", size = 0.25),
    panel.grid.major.y = element_line(color = "grey89", size = 0.25),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = size_panel_border),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", size = size_ticks),
    axis.text.x = element_text(
      family = font_family, face = font_face, size = size_axis_text.x,
      colour = "black", hjust = x_hjust, angle = angle_text.x
    ),
    axis.text.y = element_text(
      family = font_family, face = font_face, size = size_axis_text.y,
      colour = "black", vjust = y_vjust, angle = angle_text.y
    ),
    axis.title.x = element_text(family = font_family, face = font_face, size = size_axis_title.x, colour = "black"),
    axis.title.y = element_text(family = font_family, face = font_face, size = size_axis_title.y, colour = "black"),
    legend.background = element_blank(),
    legend.key.size = unit(legend_proportion_size, "line"),
    legend.title = element_text(size = size_title_text, family = font_family, face = font_face, colour = "black"),
    legend.key = element_blank(),
    legend.text = element_text(size = size_legend_text, family = font_family, face = font_face, colour = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = size_axis_title.x)
  )
}

oh.save.pdf <- function(p = NULL, outname = NULL, outdir = "figures/", width = 20,
                       height = 15, family = "Arial", fallback_resolution = 1200,
                       antialias = "default", pointsize = 12) {
  dir.create(outdir, showWarnings = FALSE)
  myfilename <- paste(outdir, outname, sep = "/")
  cairo_pdf(
    filename = myfilename, onefile = FALSE, fallback_resolution = fallback_resolution,
    width = width, height = height, family = family, antialias = antialias,
    pointsize = pointsize
  )
  print(p)
  dev.off()
}

extract_tree_data <- function(tree_disp, displayorder = TRUE) {
  td_out <- tree_disp$data
  if (displayorder) {
    td_out <- dplyr::arrange(td_out, y)
  }
  td_out
}

paleta_gt <- c("#1082A2", "#A0CC2C", "#12284C", "#DD14D3")
names(paleta_gt) <- c("0/0", "0/1", "1/1", "No coverage")

# Verdict strip: fixed colors and legend order (canonical names)
VERDICT_FILL <- c(
  Pathogenic = "#E41A1CFF",
  Deleterious = "#984EA3FF",
  Benign = "#377EB8FF",
  Other = "#A65628FF",
  Uncertain = "#999999FF"
)
VERDICT_LEVELS <- names(VERDICT_FILL)

canon_verdict_factor <- function(v) {
  v <- trimws(as.character(v))
  vl <- tolower(v)
  is_empty <- is.na(v) | vl == ""
  out <- rep("Other", length(v))
  out[vl == "pathogenic"] <- "Pathogenic"
  out[vl == "deleterious"] <- "Deleterious"
  out[vl == "benign"] <- "Benign"
  out[is_empty | vl == "unknown" | vl == "uncertain"] <- "Uncertain"
  factor(out, levels = VERDICT_LEVELS)
}

# Shared facet theme: free x, zero horizontal strip spacing (align with heatmap)
theme_facet_heatmap_x <- function() {
  theme(
    panel.spacing.x = unit(0, "line"),
    panel.spacing.y = unit(0, "line")
  )
}

# Tight vertical stacking under heatmap (egg row gaps + ggplot margins)
theme_track_strip <- function() {
  theme(
    plot.margin = margin(0, 3, 0, 3, "pt"),
    legend.margin = margin(0, 0, 0, 0, "pt"),
    legend.box.margin = margin(0, 0, 0, 0, "pt")
  )
}

# One row per (VariantId, LabelFacet) with provenance, Verdict, log10 p_else; VariantId factor order global.
# heatmap_facet_levels: left-to-right facet order (internal panels then Tips). Tips x-order uses tree tip order.
build_annotation_track_df <- function(
    df_vaf_int, df_placed_norm, file_vep,
    heatmap_facet_levels,
    tips_facet_label = NA_character_,
    tree_tip_labels = NULL,
    order_labels_displayed = NULL) {
  meta <- df_placed_norm
  meta <- meta[!duplicated(meta$VariantId), , drop = FALSE]
  meta$PosNum <- suppressWarnings(as.numeric(meta$Pos))

  vep_v <- read_vep_verdict_table(file_vep)
  if (nrow(vep_v)) {
    vep_v <- unique(vep_v[, c("VARIANT_ID", "Verdict"), drop = FALSE])
    vep_v <- vep_v[!duplicated(vep_v$VARIANT_ID), , drop = FALSE]
  }

  meta_cols <- intersect(
    names(meta),
    c("VariantId", "provenance", "p_else_where", "Chr", "Pos", "Branch", "PosNum", "DescendantTips")
  )
  d <- df_vaf_int %>%
    dplyr::distinct(.data$VariantId, .data$LabelFacet, .data$PlacedNodeId) %>%
    dplyr::left_join(meta[, meta_cols, drop = FALSE], by = "VariantId")
  if (nrow(vep_v)) {
    d <- dplyr::left_join(d, vep_v, by = c("VariantId" = "VARIANT_ID"))
  } else {
    d$Verdict <- NA_character_
  }

  d$provenance <- ifelse(is.na(d$provenance), "unknown", as.character(d$provenance))
  d$Verdict <- ifelse(is.na(d$Verdict) | !nzchar(as.character(d$Verdict)), "Unknown", as.character(d$Verdict))
  d$Verdict <- canon_verdict_factor(d$Verdict)
  prov_u <- unique(d$provenance)
  prov_lv <- c(
    intersect(c("phylogeny_filtered_variant_placement", "unfiltered_variant_placement", "unknown"), prov_u),
    setdiff(prov_u, c("phylogeny_filtered_variant_placement", "unfiltered_variant_placement", "unknown"))
  )
  d$provenance <- factor(d$provenance, levels = prov_lv)
  d$log10_p_else <- log10_p_else_capped(d$p_else_where)

  d$LabelFacet <- factor(as.character(d$LabelFacet), levels = heatmap_facet_levels)
  facet_present <- unique(as.character(d$LabelFacet))
  facet_lv <- heatmap_facet_levels[heatmap_facet_levels %in% facet_present]
  stray <- setdiff(facet_present, facet_lv)
  if (length(stray)) {
    facet_lv <- c(facet_lv, stray)
  }

  use_tips_order <- !is.na(tips_facet_label) && nzchar(tips_facet_label) &&
    !is.null(tree_tip_labels) && !is.null(order_labels_displayed)

  ord_ids <- c()
  for (lf in facet_lv) {
    sub <- d[as.character(d$LabelFacet) == lf, , drop = FALSE]
    if (use_tips_order && lf == tips_facet_label) {
      pni <- suppressWarnings(as.integer(sub$PlacedNodeId))
      ntip <- length(tree_tip_labels)
      tip_rank <- rep(NA_integer_, nrow(sub))
      ok <- !is.na(pni) & pni >= 1L & pni <= ntip
      tip_rank[ok] <- match(tree_tip_labels[pni[ok]], order_labels_displayed)
      tip_rank[is.na(tip_rank)] <- .Machine$integer.max
      sub <- sub[order(tip_rank, sub$Chr, sub$PosNum, sub$VariantId), , drop = FALSE]
    } else {
      sub <- sub[order(sub$PlacedNodeId, sub$Chr, sub$PosNum, sub$VariantId), , drop = FALSE]
    }
    ord_ids <- c(ord_ids, unique(as.character(sub$VariantId)))
  }
  d$VariantId <- factor(d$VariantId, levels = ord_ids)
  d
}

apply_variant_factor_to_heatmap_df <- function(df_int, variant_levels) {
  df_int$VariantId <- factor(df_int$VariantId, levels = variant_levels)
  df_int
}

provenance_fill_scale <- function() {
  manual <- c(
    phylogeny_filtered_variant_placement = "#d9a500",
    unfiltered_variant_placement = "#013b75",
    unknown = "#B8B8B8"
  )
  scale_fill_manual(values = manual, name = "Provenance", na.value = "#DDDDDD", drop = FALSE)
}

verdict_fill_scale <- function() {
  ggplot2::scale_fill_manual(
    values = VERDICT_FILL,
    breaks = VERDICT_LEVELS,
    limits = VERDICT_LEVELS,
    drop = FALSE,
    name = "Verdict",
    na.value = VERDICT_FILL[["Uncertain"]]
  )
}

# Three faceted panels: same LabelFacet, same VariantId x order as heatmap
build_annotation_tracks <- function(df_ann, show_facet_strips = FALSE) {
  st <- if (show_facet_strips) {
    theme()
  } else {
    theme(
      strip.text.x = element_blank(),
      strip.background.x = element_blank()
    )
  }

  # geom_tile is more reliable than geom_raster for single-row faceted strips (visible row height)
  p_prov <- ggplot(df_ann, aes(x = .data$VariantId, y = 1L, fill = .data$provenance)) +
    geom_tile(width = 1, height = 1, linewidth = 0) +
    facet_grid(. ~ LabelFacet, space = "fixed", scales = "free_x") +
    theme_ohchibi(size_panel_border = 0.3) +
    theme_facet_heatmap_x() +
    st +
    theme_track_strip() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 8, angle = 90),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.5, 1.5), breaks = NULL) +
    scale_x_discrete(expand = c(0, 0), drop = TRUE) +
    provenance_fill_scale() +
    ylab("Prov.")

  p_ver <- ggplot(df_ann, aes(x = .data$VariantId, y = 1L, fill = .data$Verdict)) +
    geom_tile(width = 1, height = 1, linewidth = 0) +
    facet_grid(. ~ LabelFacet, space = "fixed", scales = "free_x") +
    theme_ohchibi(size_panel_border = 0.3) +
    theme_facet_heatmap_x() +
    st +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 8, angle = 90),
      legend.position = "left",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.key.size = unit(0.7, "line")
    ) +
    theme_track_strip() +
    scale_y_continuous(expand = c(0, 0), limits = c(0.5, 1.5), breaks = NULL) +
    scale_x_discrete(expand = c(0, 0), drop = TRUE) +
    verdict_fill_scale() +
    ylab("Verdict")

  p_p <- ggplot(df_ann, aes(x = .data$VariantId, y = 1L, fill = .data$log10_p_else)) +
    geom_tile(width = 1, height = 1, linewidth = 0) +
    facet_grid(. ~ LabelFacet, space = "fixed", scales = "free_x") +
    theme_ohchibi(size_panel_border = 0.3) +
    theme_facet_heatmap_x() +
    st +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 8, angle = 90),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(0.8, "cm")
    ) +
    theme_track_strip() +
    scale_x_discrete(expand = c(0, 0), drop = TRUE) +
    scale_y_continuous(expand = c(0, 0), limits = c(0.5, 1.5), breaks = NULL) +
    paletteer::scale_fill_paletteer_c(
      "pals::ocean.solar",
      limits = c(-20, 0),
      oob = scales::squish,
      na.value = "#B0B0B0",
      name = "log10\np_else"
    ) +
    ylab("p_else")

  list(p_prov = p_prov, p_ver = p_ver, p_p = p_p)
}

# Fresh void panel each call (avoid reusing one ggplot object three times in ggarrange)
new_blank_panel <- function() {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"))
}

# 4 rows x 2 cols: tree | heatmap, blank | log10(p_else) raster, blank | Verdict, blank | provenance
assemble_tree_with_heatmap_column <- function(p_tree, p_heat, track_list) {
  egg::ggarrange(
    p_tree,
    p_heat,
    new_blank_panel(),
    track_list$p_p,
    new_blank_panel(),
    track_list$p_ver,
    new_blank_panel(),
    track_list$p_prov,
    nrow = 4,
    ncol = 2,
    widths = c(0.4, 1),
    heights = c(1, 0.04, 0.02, 0.02),
    padding = unit(0, "pt"),
    draw = FALSE
  )
}

# Read NR/NV matrix into a numeric matrix with VariantId as row names.
# Handles two on-disk formats:
#   1. Leading-tab header (\tSample1\tSample2...) — produced by the bash build_matrix function.
#      read.table(sep="") correctly infers row names here, but fread (used below) renames the
#      empty first column to "V1", so the explicit check below handles both cases uniformly.
#   2. Explicit named first column (V1\tSample1\tSample2...) — produced when R scripts use
#      fread+fwrite on these matrices (fwrite writes the "V1" column name literally).
#      In this case read.table(sep="") sees equal header/data field counts and does NOT infer
#      row names, causing as.matrix() to return a character matrix and value.x/value.y to be
#      non-numeric. fread with explicit sep="\t" + the check below handles this correctly.
read_nr_nv_matrix <- function(f) {
  dt <- data.table::fread(
    f, header = TRUE, sep = "\t", check.names = FALSE,
    data.table = FALSE, stringsAsFactors = FALSE
  )
  if (nrow(dt) == 0 || ncol(dt) == 0) return(matrix(numeric(0), 0, 0))
  if (!is.numeric(dt[[1]])) {
    rn <- as.character(dt[[1]])
    dt <- dt[, -1, drop = FALSE]
    rownames(dt) <- rn
  }
  m <- as.matrix(dt)
  storage.mode(m) <- "numeric"
  m
}

# --- Core data path (single read / single tree build) ---------------------------------
prepare_vaf_heatmap_bundle <- function(
    file_placed, file_nv, file_nr, file_tree, file_vep,
    heatmap_y_axis_text_size = 5, plot_tips = FALSE) {
  df_placed_norm <- normalize_placed_variants(read_placement_tsv(file_placed))
  if ("provenance" %in% names(df_placed_norm)) {
    df_placed_norm <- dplyr::filter(
      df_placed_norm,
      provenance == "phylogeny_filtered_variant_placement"
    )
    if (!nrow(df_placed_norm)) {
      warning("No rows with provenance phylogeny_filtered_variant_placement; plots may be empty.")
    }
  }

  df_nv <- read_nr_nv_matrix(file_nv) %>% reshape2::melt()
  df_nr <- read_nr_nv_matrix(file_nr) %>% reshape2::melt()

  df_vaf <- merge(df_nv, df_nr, by = c("Var1", "Var2")) %>%
    dplyr::mutate(VAF = value.x / value.y) %>%
    dplyr::rename(VariantId = Var1, biosampleName = Var2, NV = value.x, NR = value.y) %>%
    dplyr::filter(VariantId %in% df_placed_norm$VariantId)

  df_vaf$PlacedNodeId <- df_placed_norm$Branch[match(df_vaf$VariantId, df_placed_norm$VariantId)]
  df_vaf <- merge_vep_relevance_into_vaf(df_vaf, file_vep)

  rm(df_nv, df_nr)
  gc()

  tb <- table(df_placed_norm$Branch)
  df_freq <- data.frame(
    node = as.numeric(names(tb)),
    Freq = as.numeric(tb),
    stringsAsFactors = FALSE
  )
  df_freq <- df_freq[order(-df_freq$Freq), , drop = FALSE]

  tree <- ape::read.tree(file = file_tree)
  df_placed_norm <- append_descendant_tip_labels(df_placed_norm, tree)
  # tip i (1..Ntip) matches PlacedNodeId for tip placements; capture before treedata/full_join (no $ on S4).
  phy_tip_labels <- as.character(tree$tip.label)
  tree2 <- ape::ladderize(tree, right = TRUE)
  is_internal <- tree2$edge[, 2] > length(tree2$tip.label)
  order_nodes <- tree2$edge[is_internal, 2]

  df_nodes_ids <- tidytree::as_tibble(ape::makeNodeLabel(tree, method = "md5sum")) %>% as.data.frame()
  df_vaf$PlacedNodeLabel <- df_nodes_ids$label[match(df_vaf$PlacedNodeId, df_nodes_ids$node)]

  tree$edge.length <- rep(1, nrow(tree$edge))
  tree <- full_join(tree, df_freq, by = "node")

  p_tree <- ggtree(tr = tree, aes(color = Freq), ladderize = FALSE, size = 2) +
    paletteer::scale_color_paletteer_c(
      "viridis::plasma",
      na.value = "#D9D9D9",
      oob = scales::squish,
      trans = "log10",
      name = "Count\nplaced mutations"
    ) +
    geom_tiplab(size = 0, align = TRUE) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9)
    )
  p_tree$data$InternalNode <- p_tree$data$node
  p_tree <- p_tree +
    geom_label(aes(label = InternalNode), color = "black", size = 2.5) +
    geom_tiplab(size = 0, align = TRUE, color = "black") +
    coord_cartesian(clip = "off")

  order_labels_displayed <- extract_tree_data(p_tree) %>%
    subset(isTip == TRUE) %>%
    as.data.frame() %>%
    dplyr::pull(.data$label) %>%
    as.character()

  mnums <- df_placed_norm[, c("VariantId", "Branch")] %>% unique() %$% Branch %>% table()
  df_numplaced <- data.frame(Branch = names(mnums), NumPlaced = as.numeric(mnums))
  df_numplaced$Label <- paste0(df_numplaced$Branch, " (n=", df_numplaced$NumPlaced, ")")

  br_match <- match(as.character(df_vaf$PlacedNodeId), as.character(df_numplaced$Branch))
  df_vaf$NumPlacedVariantsInNodeId <- df_numplaced$NumPlaced[br_match]
  df_vaf$LabelFacet <- df_numplaced$Label[br_match]

  ntip <- ape::Ntip(tree)
  pid_v <- suppressWarnings(as.integer(df_vaf$PlacedNodeId))
  is_tip_row <- !is.na(pid_v) & pid_v >= 1L & pid_v <= ntip
  n_tip_variants <- length(unique(df_vaf$VariantId[is_tip_row]))
  tips_facet_label <- if (isTRUE(plot_tips) && n_tip_variants > 0L) {
    paste0("Tips (n=", n_tip_variants, ")")
  } else {
    NA_character_
  }
  if (!is.na(tips_facet_label)) {
    df_vaf$LabelFacet[is_tip_row] <- tips_facet_label
  }

  # Internal nodes only unless plot_tips: include tip-placed variants in their own "Tips (...)" facet.
  keep_heat <- (!is.na(pid_v) & pid_v %in% order_nodes) | (isTRUE(plot_tips) & is_tip_row)
  df_vaf_int <- df_vaf %>%
    dplyr::filter(keep_heat) %>%
    droplevels()
  df_vaf_int$biosampleName <- factor(df_vaf_int$biosampleName, levels = order_labels_displayed)

  facet_internal_lvls <- c()
  for (nd in order_nodes) {
    hit <- which(suppressWarnings(as.integer(df_numplaced$Branch)) == nd)
    if (!length(hit)) {
      hit <- which(as.character(df_numplaced$Branch) == as.character(nd))
    }
    if (!length(hit)) {
      next
    }
    lab <- df_numplaced$Label[hit[1]]
    if (!is.na(lab) && lab %in% df_vaf_int$LabelFacet) {
      facet_internal_lvls <- c(facet_internal_lvls, lab)
    }
  }
  facet_internal_lvls <- unique(facet_internal_lvls)
  order_nl <- facet_internal_lvls
  if (isTRUE(plot_tips) && !is.na(tips_facet_label) &&
        any(as.character(df_vaf_int$LabelFacet) == tips_facet_label, na.rm = TRUE)) {
    order_nl <- c(order_nl, tips_facet_label)
  }
  df_vaf_int$LabelFacet <- factor(as.character(df_vaf_int$LabelFacet), levels = order_nl)

  df_ann <- build_annotation_track_df(
    df_vaf_int, df_placed_norm, file_vep,
    heatmap_facet_levels = order_nl,
    tips_facet_label = tips_facet_label,
    tree_tip_labels = phy_tip_labels,
    order_labels_displayed = order_labels_displayed
  )
  v_levels <- levels(df_ann$VariantId)
  df_vaf_int <- apply_variant_factor_to_heatmap_df(df_vaf_int, v_levels)

  df_vaf_int$VAF[df_vaf_int$VAF == 0] <- NA

  tracks <- build_annotation_tracks(df_ann, show_facet_strips = FALSE)

  p_vaf <- ggplot(df_vaf_int, aes(VariantId, biosampleName)) +
    geom_raster(aes(fill = VAF)) +
    facet_grid(. ~ LabelFacet, space = "fixed", scales = "free") +
    theme_ohchibi(size_panel_border = 0.3) +
    theme_facet_heatmap_x() +
    theme(
      plot.margin = margin(5.5, 5.5, 0, 5.5, "pt"),
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0, "pt"),
      legend.box.spacing = unit(0.15, "line"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = heatmap_y_axis_text_size),
      axis.ticks.y = element_line(size = 0.2),
      strip.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 0),
      strip.text.y = element_text(size = 8),
      legend.text = element_text(size = 9, angle = 90, vjust = 0, hjust = 0),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 9)
    ) +
    scale_fill_paletteer_c("pals::kovesi.rainbow_bgyrm_35_85_c71", na.value = "#F5F5F5") +
    scale_x_discrete(expand = c(0, 0), drop = TRUE) +
    scale_y_discrete(expand = c(0, 0))

  list(
    p_tree = p_tree,
    p_vaf = p_vaf,
    tracks = tracks,
    df_vaf_int = df_vaf_int,
    df_ann = df_ann,
    df_placed_norm = df_placed_norm,
    order_nodes = order_nodes,
    order_labels_displayed = order_labels_displayed,
    phy_tip_labels = phy_tip_labels,
    heatmap_facet_levels = order_nl,
    tips_facet_label = tips_facet_label,
    df_numplaced = df_numplaced,
    df_node_freq = df_freq,
    tree_ggtree_input = tree,
    heatmap_y_axis_text_size = heatmap_y_axis_text_size,
    plot_tips = plot_tips
  )
}

add_genotype_and_build_digital <- function(bundle, file_gt) {
  y_ax_sz <- bundle$heatmap_y_axis_text_size
  if (is.null(y_ax_sz) || !is.finite(y_ax_sz) || y_ax_sz <= 0) {
    y_ax_sz <- 5
  }
  df <- bundle$df_vaf_int
  df_gt <- read.table(file = file_gt) %>%
    dplyr::rename(biosampleName = V1, VariantId = V2, GT = V3)
  df_gt_lookup <- df_gt
  df$VariantId <- as.character(df$VariantId)
  df <- merge(df, df_gt, by = c("biosampleName", "VariantId"), all.x = TRUE)
  df$GT[is.na(df$GT)] <- "0/0"
  df$NR[df$NR == 0] <- NA
  df$GTPlot <- df$GT
  df$GTPlot[is.na(df$NR)] <- "No coverage"
  df$GTPlot <- df$GTPlot %>%
    gsub(pattern = "1/0", replacement = "0/1") %>%
    factor()
  df$VariantId <- factor(df$VariantId, levels = levels(bundle$df_ann$VariantId))

  p_dig <- ggplot(df, aes(VariantId, biosampleName)) +
    geom_raster(aes(fill = GTPlot), alpha = 1) +
    facet_grid(. ~ LabelFacet, space = "fixed", scales = "free") +
    theme_ohchibi(size_panel_border = 0.3) +
    theme_facet_heatmap_x() +
    theme(
      plot.margin = margin(5.5, 5.5, 0, 5.5, "pt"),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = y_ax_sz),
      axis.ticks.y = element_line(size = 0.2),
      strip.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 0),
      strip.text.y = element_text(size = 8),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 9)
    ) +
    scale_fill_manual(values = paleta_gt, name = "Genotype", drop = FALSE) +
    scale_x_discrete(expand = c(0, 0), drop = TRUE) +
    scale_y_discrete(expand = c(0, 0))

  list(p_digital = p_dig, df_vaf_int = df, df_gt_lookup = df_gt_lookup)
}

# CLI: placement_tsv nv.tsv nr.tsv tree newick  gt_tsv  [vep_csq_master_relevance.tsv]  [heatmap_y_axis_text_size]  [plot_tips]
# plot_tips: "true"/"false" (argv 8, optional; default false) — include tip-placed variants in a "Tips (n=...)" facet on heatmaps and rasters.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Need at least 5 args: placement_tsv nv nr tree gt [vep_tsv [heatmap_y_axis_text_size [plot_tips]]]")
}
file_placed_variants <- args[1]
file_nv <- args[2]
file_nr <- args[3]
file_tree <- args[4]
file_gt <- args[5]
file_vep <- if (length(args) >= 6) args[6] else ""
heatmap_y_axis_text_size <- if (length(args) >= 7) {
  suppressWarnings(as.numeric(args[7]))
} else {
  5
}
if (!is.finite(heatmap_y_axis_text_size) || heatmap_y_axis_text_size <= 0) {
  heatmap_y_axis_text_size <- 5
}
plot_tips <- FALSE
if (length(args) >= 8) {
  v <- tolower(trimws(as.character(args[8])))
  plot_tips <- v %in% c("true", "1", "yes", "t")
}

bundle <- prepare_vaf_heatmap_bundle(
  file_placed_variants, file_nv, file_nr, file_tree, file_vep,
  heatmap_y_axis_text_size = heatmap_y_axis_text_size,
  plot_tips = plot_tips
)

comp_vaf <- assemble_tree_with_heatmap_column(bundle$p_tree, bundle$p_vaf, bundle$tracks)
oh.save.pdf(p = comp_vaf, outname = "res_composition.pdf", outdir = "./", width = 22, height = 18)

dig <- add_genotype_and_build_digital(bundle, file_gt)
comp_digital <- assemble_tree_with_heatmap_column(bundle$p_tree, dig$p_digital, bundle$tracks)
oh.save.pdf(p = comp_digital, outname = "res_composition_digital.pdf", outdir = "./", width = 22, height = 18)

# RDS contents:
#   composition_vaf / composition_digital — full egg::ggarrange layouts (re-print for exact PDF replay).
#   df_vaf_int — merged VAF + genotype rows used for digital geom_raster (also has VAF/NV/NR for VAF heatmap).
#   df_heatmap_vaf — same rows as used to build p_vaf (before genotype merge); geom_raster(VAF) + facets.
#   df_ann — one row per (VariantId, LabelFacet) for provenance / Verdict / log10_p_else strip geom_tile panels
#     (includes Tips facet rows; DescendantTips from placement when present).
#   heatmap_facet_levels, tips_facet_label — facet order and pooled tip panel label.
#   phy_tip_labels — as.character(read.tree()$tip.label); tip node i -> label for Tips x-ordering.
#   df_placed_norm — placement TSV after normalize_placed_variants (phylo filter); includes DescendantTips
#     (tip.label under Branch node, from the same Newick as file_tree; "; "-separated, sorted).
#   df_gt_lookup — three-column genotype table from the GT file (biosampleName, VariantId, GT).
#   df_vep_verdict — VARIANT_ID + Verdict from VEP relevance TSV (empty if no file); used when rebuilding df_ann.
#   df_node_freq, df_numplaced, order_nodes, order_labels_displayed, tree_ggtree_input — tree + heatmap ordering.
#   scale_constants — fill palettes / levels for geom_raster and geom_tile strips if rebuilding ggplot from dfs.
#   p_*, tracks — ggplot parts; reassembly needs assemble_tree_with_heatmap_column() from this script + libs.
df_heatmap_vaf <- bundle$df_vaf_int
df_vep_verdict <- read_vep_verdict_table(file_vep)
scale_constants <- list(
  paleta_gt = paleta_gt,
  verdict_fill = VERDICT_FILL,
  verdict_levels = VERDICT_LEVELS,
  vaf_fill_scale = list(package = "pals", palette = "kovesi.rainbow_bgyrm_35_85_c71", na_value = "#F5F5F5"),
  tree_count_color_scale = list(paletteer_d = "viridis::plasma", na_value = "#D9D9D9", trans = "log10"),
  log10_p_else_limits = c(-20, 0),
  log10_p_else_fill_scale = list(
    paletteer_d = "pals::ocean.solar", na_value = "#B0B0B0", limits = c(-20, 0)
  ),
  provenance_colors = c(
    phylogeny_filtered_variant_placement = "#d9a500",
    unfiltered_variant_placement = "#013b75",
    unknown = "#B8B8B8"
  )
)
res_list <- list(
  composition_vaf = comp_vaf,
  composition_digital = comp_digital,
  heatmap_y_axis_text_size = heatmap_y_axis_text_size,
  scale_constants = scale_constants,
  cli_args = args,
  p_tree = bundle$p_tree,
  p_vaf = bundle$p_vaf,
  p_vaf_digital = dig$p_digital,
  df_vaf_int = dig$df_vaf_int,
  df_heatmap_vaf = df_heatmap_vaf,
  df_ann = bundle$df_ann,
  df_placed_norm = bundle$df_placed_norm,
  df_gt_lookup = dig$df_gt_lookup,
  df_vep_verdict = df_vep_verdict,
  df_node_freq = bundle$df_node_freq,
  df_numplaced = bundle$df_numplaced,
  order_nodes = bundle$order_nodes,
  order_labels_displayed = bundle$order_labels_displayed,
  phy_tip_labels = bundle$phy_tip_labels,
  heatmap_facet_levels = bundle$heatmap_facet_levels,
  tips_facet_label = bundle$tips_facet_label,
  tree_ggtree_input = bundle$tree_ggtree_input,
  tracks = bundle$tracks,
  plot_tips = bundle$plot_tips
)
saveRDS(object = res_list, file = "res_figures.RDS")
