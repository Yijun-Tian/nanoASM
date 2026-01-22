suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(hrbrthemes)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(gridExtra)
  library(optparse)
})

# ----------------------------
# CLI
# ----------------------------
opt_list <- list(
  make_option(c("--annoDMR"), type="character",
              help="An annotated DMR file with allele specific methylation regions, SNP was cross annotated with GWAS and mQTL for the plot"),
  make_option(c("--qvalue"), type = "numeric", default = 1e-5,
              help = "FDR criteria for regions to plot heatmap"),
  make_option(c("--delta"), type = "numeric", default = 25,
              help = "Methylation delta criteria for regions to plot heatmap"),
  make_option(c("--type"), type="character",
              help="Either Normal or Cancer")
)
opt <- parse_args(OptionParser(option_list=opt_list))

qv = opt$qvalue
d = opt$delta
type = opt$type
mat_dir <- paste0(type, "MAT")
if (!dir.exists(mat_dir)) {
  stop("Matrix directory not found: ", mat_dir)
}

# ----------------------------
# Column names
# ----------------------------
colnames_dmr = c("CHR", "L", "R", "FDR", "Delta", "nCpG", "pMWU", "pKS", "PCTREF", "PCTALT",
                 "SNP", "REF", "ALT", "SNPCHR", "SNPLOC", "nREF", "nALT", "DIST2L", "DIST2R",
                 "GWAS", "mQTL"
                 )

# ----------------------------
# Read + filter
# ----------------------------
dmr <- fread(opt$annoDMR, fill = TRUE, header = FALSE)
if (ncol(dmr) < length(colnames_dmr)) {
  stop("Input has ", ncol(dmr), " columns but expected at least ", length(colnames_dmr),
       ". Check the annoDMR format.")
}
setnames(dmr, colnames_dmr)

summary <- dmr %>%
  dplyr::filter(FDR <= qv & abs(Delta) >= d) %>%
  #dplyr::select(-c(asId, GT:SE)) %>%
  distinct() %>%
  dplyr::arrange(
    dplyr::desc((!is.na(GWAS)) + (!is.na(mQTL))),  # 2 > 1 > 0
    dplyr::desc(!is.na(GWAS)),                     # GWAS-only above mQTL-only
    dplyr::desc(abs(Delta))
  )

if (nrow(summary) == 0) {
  stop("No regions pass filters: FDR <= ", qv, " and abs(Delta) >= ", d)
}

# ----------------------------
# Helpers
# ----------------------------
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

nearest_idx <- function(vec, x) {
  which.min(abs(vec - x))
}

# Get x-index for a coordinate target; if exact match missing, use nearest
coord_to_rank <- function(map_df, target_coord) {
  idx <- which(map_df$coordinate == target_coord)
  if (length(idx) > 0) {
    return(map_df$X[idx[1]])
  }
  map_df$X[nearest_idx(map_df$coordinate, target_coord)]
}


p <- list()
for (i in 1:dim(summary)[1]){
  row <- summary[i, ]
  
  chr     <- as.character(row$CHR)
  snp     <- as.character(row$SNP)
  dmrL0   <- as.integer(row$L)
  dmrR0   <- as.integer(row$R)
  snploc  <- as.integer(row$SNPLOC)
  
  pct_ref <- as.numeric(row$PCTREF)
  pct_alt <- as.numeric(row$PCTALT)
  fdr     <- as.numeric(row$FDR)
  
  ref <- as.character(row$REF)
  alt <- as.character(row$ALT)
  
  gwas <- as.character(row$GWAS)
  mqtl <- as.character(row$mQTL)
  
  
  title <- paste0(
    snp, " (", chr, ":", snploc, ")  @", chr, ":", dmrL0, "-", dmrR0, "\n",
    "ref-mean: ", fmt_pct(pct_ref), "%; alt-mean: ", fmt_pct(pct_alt), "%; FDR = ", fdr, "\n",
    "(", gwas, " | ", mqtl, ")"
  )
  
  mat_file <- file.path(mat_dir, paste0(snp, ".matrix.txt"))
  if (!file.exists(mat_file)) {
    warning("Missing matrix file for ", snp, ": ", mat_file, " (skipping)")
    next
  }
  
  # Read header line safely (one line)
  header_line <- fread(mat_file, nrows = 0) %>% names()
  
  # Normalize column names
  header_line <- header_line %>%
    gsub("chrom", "chromosome", ., fixed = TRUE) %>%
    gsub("pos",   "coordinate", ., fixed = TRUE)
  
  # Read matrix data
  mat <- fread(mat_file, skip = 1)
  if (ncol(mat) != length(header_line)) {
    stop("Column count mismatch in ", mat_file,
         ": data has ", ncol(mat), " but header has ", length(header_line))
  }
  setnames(mat, header_line)
  
  # Basic sanity
  if (!all(c("chromosome", "coordinate") %in% names(mat))) {
    stop("Matrix missing required columns chromosome/coordinate: ", mat_file)
  }
  
  mat_tbl <- as_tibble(mat) %>%
    mutate(
      coordinate = as.integer(coordinate),
      X = rank(coordinate, ties.method = "first")
    )
  
  # infer sample counts by allele prefix in colnames (excluding chromosome/coordinate)
  sample_cols <- setdiff(names(mat_tbl), c("chromosome", "coordinate", "X"))
  # strip trailing _<id> from sample cols to get allele labels
  allele_labels <- gsub("_[A-Za-z0-9]+$", "", sample_cols)
  sizes <- table(allele_labels)
  
  if (!(ref %in% names(sizes)) || !(alt %in% names(sizes))) {
    warning("Cannot find ref/alt allele labels in matrix colnames for ", snp,
            " (ref=", ref, ", alt=", alt, "). Found: ",
            paste(names(sizes), collapse = ", "),
            " (skipping)")
    next
  }
  
  n_ref <- as.integer(sizes[[ref]])
  n_alt <- as.integer(sizes[[alt]])
  
  # map: X <-> coordinate
  map <- mat_tbl %>% dplyr::select(X, coordinate) %>% arrange(X)
  
  # Your original code uses dmr[1]+1 and dmr[2]; preserve intent but make robust
  dmrL_rank <- coord_to_rank(map, dmrL0 + 1L)
  dmrR_rank <- coord_to_rank(map, dmrR0)
  
  snp_rank_idx <- nearest_idx(map$coordinate, snploc)
  snp_rank <- map$X[snp_rank_idx]
  
  brk <- unique(round(quantile(map$X, probs = seq(0, 1, 0.1), names = FALSE)))
  brk <- brk[brk >= min(map$X) & brk <= max(map$X)]
  
  # horizontal boundary line between ref and alt rows
  # (your original logic kept; it’s basically "ref block on top vs bottom")
  y_split <- if (ref < alt) (n_ref + 0.5) else (n_alt + 0.5)
  
  long_df <- mat_tbl %>%
    dplyr::select(-chromosome, -coordinate) %>%
    pivot_longer(
      cols = -X,
      names_to = "Y",
      values_to = "Z"
    ) %>%
    mutate(
      type = if_else(str_detect(Y, paste0("^", ref, "_")), "REF", "ALT"),
      Z = as.numeric(Z)
    )
  p[[i]] <- ggplot(long_df, aes(X, Y, fill = Z)) +
    geom_tile() +
    coord_fixed(ratio = 4) +
    scale_y_discrete(name = "allelic_summary") +
    scale_x_continuous(
      name = "genomic coordinate",
      breaks = brk,
      labels = map$coordinate[match(brk, map$X)]
    ) +
    scale_fill_gradient2(
      high = "red", low = "blue", midpoint = 50,
      limits = c(0, 100),
      name = "%meth%"
    ) +
    geom_vline(xintercept = dmrL_rank, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_vline(xintercept = dmrR_rank, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_vline(xintercept = snp_rank + 0.5, linetype = "dashed", colour = "green", linewidth = 0.8) +
    geom_hline(yintercept = y_split, linetype = "dashed", colour = "black", linewidth = 0.8) +
    ggtitle(title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
}

# Remove NULLs from skipped SNPs
p <- p[!vapply(p, is.null, logical(1))]
if (length(p) == 0) stop("No plots generated (all skipped).")

out_pdf <- paste0(type, ".HeatmapASM.pdf")
ggsave(
  filename = out_pdf,
  plot = marrangeGrob(p, nrow = 2, ncol = 1),
  width = 8.5, height = 11
)
message("Saved: ", out_pdf)
