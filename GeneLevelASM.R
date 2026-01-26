#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(gridExtra)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
})

# ----------------------------
# CLI
# ----------------------------
opt_list <- list(
  make_option(c("--annowindow"), type="character",
              help="Annotated window file with merged ASM DMR regions (cross-annotated with GTEx eQTL/GWAS/mQTL)"),
  make_option(c("--qvalue"), type="numeric", default=0.01,
              help="FDR cutoff for windows to plot"),
  make_option(c("--delta"), type="numeric", default=25,
              help="Abs(Delta) cutoff (%) for windows to plot"),
  make_option(c("--type"), type="character",
              help="Either Normal or Cancer (used for output name)"),
  make_option(c("--out"), type="character", default=NULL,
              help="Output PDF (default: <type>.GeneLevelASM.pdf)")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.null(opt$annowindow), file.exists(opt$annowindow))
stopifnot(!is.null(opt$type), opt$type %in% c("Normal", "Cancer"))

qv   <- opt$qvalue
dcut <- opt$delta
type <- opt$type
out_pdf <- if (is.null(opt$out)) paste0(type, ".GeneLevelASM.pdf") else opt$out

# ----------------------------
# Columns + input
# ----------------------------
colnames_window <- c(
  "asId", "windowCHR", "windowL", "windowR",
  "CHR", "L", "R", "FDR", "Delta", "nCpG", "pMWU", "pKS", "PCTREF", "PCTALT",
  "SNP", "REF", "ALT", "SNPCHR", "SNPLOC", "nREF", "nALT", "DIST2L", "DIST2R",
  "RegionSize", "R2W_ratio", "nRegion", "sd_R2Wratio", "mean_R2Wratio",
  "GT", "GENE", "EVENT", "eGene", "eSNP", "dis2TSS", "maf", "maSample", "sample", "pvalue", "NES", "SE", "GWAS", "mQTL")

windows <- fread(opt$annowindow, header = FALSE, fill = TRUE)
if (ncol(windows) < length(colnames_window)) {
  stop("Input has ", ncol(windows), " cols but expected at least ", length(colnames_window),
       ". Check annowindow format.")
}
setnames(windows, colnames_window)

# Ensure key columns are numeric (prevents subtle filter issues)
windows <- windows %>%
  mutate(
    FDR    = suppressWarnings(as.numeric(FDR)),
    Delta  = suppressWarnings(as.numeric(Delta)),
    L      = suppressWarnings(as.integer(L)),
    R      = suppressWarnings(as.integer(R)),
    SNPLOC = suppressWarnings(as.integer(SNPLOC)),
    pvalue = suppressWarnings(as.numeric(pvalue))
  )

# Filter windows to plot
hits <- windows %>%
  dplyr::filter(!is.na(GENE), !is.na(FDR), !is.na(Delta), FDR <= qv, abs(Delta) >= dcut) %>%
  dplyr::arrange(dplyr::desc(abs(Delta)))

genelist <- unique(hits$GENE)
genelist <- genelist[!is.na(genelist)]

if (length(genelist) == 0) {
  stop("No genes pass filters: FDR <= ", qv, " and abs(Delta) >= ", dcut)
}

# ----------------------------
# Helpers
# ----------------------------
fmt_assoc <- function(gwas, mqtl) {
  gwas_flag <- !is.na(gwas) && gwas != "" && gwas != "NA"
  mqtl_flag <- !is.na(mqtl) && mqtl != "" && mqtl != "NA"
  
  if (gwas_flag && mqtl_flag) "GWAS + mQTL"
  else if (gwas_flag) "GWAS"
  else if (mqtl_flag) "mQTL"
  else NA_character_
}


# Build a single representative transcript per gene.
# Rule: choose transcript with largest span (deterministic).
ensg <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  keys    = genelist,
  keytype = "GENEID",
  columns = c("GENEID", "GENENAME", "GENESEQSTART", "GENESEQEND",
              "TXSEQSTART", "TXSEQEND", "TXID", "TXNAME")
) %>%
  as_tibble()

gene_symbol_map <- ensg %>%
  dplyr::filter(!is.na(GENEID)) %>%
  group_by(GENEID) %>%
  summarise(GENESYMBOL = dplyr::first(na.omit(GENENAME)), .groups = "drop")

ensg_one <- ensg %>%
  dplyr::filter(!is.na(GENEID)) %>%
  dplyr::mutate(
    gs = suppressWarnings(as.numeric(GENESEQSTART)),
    ge = suppressWarnings(as.numeric(GENESEQEND)),
    ts = suppressWarnings(as.numeric(TXSEQSTART)),
    te = suppressWarnings(as.numeric(TXSEQEND)),
    tx_span = abs(te - ts)
  ) %>%
  dplyr::group_by(GENEID) %>%
  dplyr::slice_max(tx_span, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(GENE = GENEID, gs, ge, ts, te)

# Smart label positioning:
# - label right by default
# - if SNPLOC is in the rightmost 15% of the current panel range, place label to the left
add_smart_label_pos <- function(df, right_edge_frac = 0.75, offset_frac = 0.02) {
  xr <- range(c(df$L, df$R, df$SNPLOC), na.rm = TRUE)
  span <- diff(xr)
  if (!is.finite(span) || span == 0) span <- 1
  
  threshold <- xr[1] + right_edge_frac * span
  offset <- offset_frac * span
  
  df %>%
    mutate(
      label_hjust = if_else(SNPLOC >= threshold, 1, 0),
      label_x     = if_else(SNPLOC >= threshold, SNPLOC - offset, SNPLOC + offset)
    )
}

make_block <- function(windows, gene_id, dcut) {
  windows %>%
    dplyr::filter(GENE == gene_id) %>%
    dplyr::select(-c(GT, EVENT)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(windowCHR, windowL, windowR) %>%
    dplyr::mutate(
      dist = SNPLOC - (L + R) / 2,
      ys = rank(dist, ties.method = "first"),
      ye = ys + 1
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      eQTL  = dplyr::if_else(!is.na(pvalue) & pvalue <= 0.05, "significant", "no"),
      assoc = purrr::pmap_chr(list(GWAS, mQTL), fmt_assoc),
      label_text = paste0(
        SNP, " (", CHR, ":", SNPLOC, ")",
        dplyr::if_else(!is.na(assoc), paste0(" ", assoc), ""),
        dplyr::if_else(!is.na(pvalue), paste0("|p=", formatC(pvalue, format = "e", digits = 2)), ""),"\n@",
        CHR, ":", L, "-", R)
      ) %>%
    add_smart_label_pos()
}
make_title <- function(gene_id, gene_symbol_map) {
  sym <- gene_symbol_map$GENESYMBOL[match(gene_id, gene_symbol_map$GENEID)]
  if (length(sym) == 0 || is.na(sym) || sym == "") {
    gene_id
  } else {
    paste0(sym, " (", gene_id, ")")
  }
}
plot_gene <- function(block, gene_id, coords, dcut, qv, gene_symbol_map) {
  ggplot(block, aes(xmin = L, xmax = R, ymin = ys, ymax = ye, fill = abs(Delta))) +
    geom_rect() +
    scale_fill_gradient2(
      high = "red", low = "blue",
      midpoint = dcut, limits = c(0, 100),
      name = "|Δmeth|"
    ) +
    # Gene boundaries (solid) + transcript boundaries (dashed), if available
    { if (!is.na(coords$gs)) geom_vline(xintercept = coords$gs, colour = "green",  linewidth = 0.7) } +
    { if (!is.na(coords$ge)) geom_vline(xintercept = coords$ge, colour = "orange", linewidth = 0.7) } +
    { if (!is.na(coords$ts)) geom_vline(xintercept = coords$ts, colour = "green",  linetype = "longdash", linewidth = 0.5) } +
    { if (!is.na(coords$te)) geom_vline(xintercept = coords$te, colour = "orange", linetype = "longdash", linewidth = 0.5) } +
    # SNP markers colored by eQTL status
    geom_segment(aes(x = SNPLOC, y = ys, yend = ye, colour = eQTL)) +
    scale_colour_manual(values = c(significant = "red", no = "blue"), name = "GTEx eQTL") +
    # Label ALL associations (GWAS/mQTL), with smart side-switch positioning
    geom_label(
      data = dplyr::filter(block, !is.na(label_text)),
      aes(x = label_x, y = (ys + ye) / 2, label = label_text, hjust = label_hjust),
      size = 2.2,
      label.size = 0.2,
      fill = NA,
      colour = "black"
    ) +
    facet_grid(
      rows = vars(paste(windowCHR, windowL, windowR, sep = "_")),
      scales = "free_y",
      space  = "free_y"
    ) +
    labs(
      title = make_title(gene_id, gene_symbol_map),
      x = "Genomic coordinate",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_grey(base_size = 6) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
}

# ----------------------------
# Build plots
# ----------------------------
plots <- vector("list", length(genelist))

for (i in seq_len(length(genelist))) {
  gene_id <- genelist[[i]]
  
  block <- make_block(windows, gene_id, dcut)
  if (nrow(block) == 0) next
  
  coords <- ensg_one %>% dplyr::filter(GENE == gene_id)
  if (nrow(coords) == 0) {
    coords <- tibble(gs = NA_real_, ge = NA_real_, ts = NA_real_, te = NA_real_)
  } else {
    coords <- coords[1, ]
  }
  
  plots[[i]] <- plot_gene(block, gene_id, coords, dcut, qv, gene_symbol_map)
}

plots <- plots[!vapply(plots, is.null, logical(1))]
if (length(plots) == 0) stop("No plots generated.")

ggsave(
  filename = out_pdf,
  plot = marrangeGrob(plots, nrow = 2, ncol = 1),
  width = 8.5,
  height = 11
)

message("Saved: ", out_pdf)
