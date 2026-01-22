suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
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
              help="Annotated DMR file (ASM regions + GTEx eQTL annotations)"),
  make_option(c("--overlap"), type="character",
              help="RDS overlap dataset (minimal significance overlap between Normal/Cancer)"),
  make_option(c("--qvalue"), type="numeric", default=1e-5,
              help="FDR cutoff"),
  make_option(c("--delta"), type="numeric", default=25,
              help="Abs(Delta) cutoff"),
  make_option(c("--type"), type="character",
              help="Either Normal or Cancer (the 'main character')"),
  make_option(c("--out"), type="character", default="CancerSpecASM.dpct15.pdf",
              help="Output PDF filename")
)

opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$annoDMR), file.exists(opt$annoDMR))
stopifnot(!is.null(opt$overlap), file.exists(opt$overlap))
stopifnot(!is.null(opt$type), opt$type %in% c("Normal", "Cancer"))

qv   <- opt$qvalue
dcut <- opt$delta
type_mc <- opt$type
type_sc <- setdiff(c("Normal", "Cancer"), type_mc)[1]
out_pdf <- opt$out

# ----------------------------
# I/O + filter
# ----------------------------
overlap_C <- readRDS(opt$overlap) %>%
  dplyr::select(c(V1, V2, V3, V11, nd, td, dmrtype, snptype)) %>%
  dplyr::filter(dmrtype %in% c("onlyC", "sharedC")) %>%
  distinct()

colnames_dmr <- c(
  "CHR","L","R","FDR","Delta","nCpG","pMWU","pKS","PCTREF","PCTALT",
  "SNP","REF","ALT","SNPCHR","SNPLOC","nREF","nALT","DIST2L","DIST2R",
  "GWAS","mQTL"
)

dmr <- fread(opt$annoDMR, header=FALSE, fill=TRUE)
setnames(dmr, colnames_dmr)

collapse_unique <- function(x, sep = "|") {
  x <- x[!is.na(x)]
  ux <- unique(x)
  if (length(ux) == 0) return(NA_character_)
  if (length(ux) == 1) return(as.character(ux))
  paste(as.character(ux), collapse = sep)
}

tnplot.table <- dmr %>%
  # dplyr::select(-c(asId, GT:SE)) %>%
  dplyr::filter(FDR <= qv & abs(Delta) >= dcut) %>%
  distinct() %>%
  left_join(., overlap_C, by = c("CHR"="V1", "L"="V2", "R"="V3", "SNP"="V11")) %>%
  mutate(abs_td_nd = abs(td - nd)) %>%
  dplyr::filter(snptype == "shared") %>%
  dplyr::filter(dmrtype == "onlyC" | abs(td - nd) >= 15) %>%
  dplyr::arrange(
    dplyr::desc((!is.na(GWAS)) + (!is.na(mQTL))),  # 2 > 1 > 0
    dplyr::desc(!is.na(GWAS)),                     # GWAS-only above mQTL-only
    dplyr::desc(abs(Delta))
  )

if (nrow(tnplot.table) == 0) stop("No rows pass filters.")

# ----------------------------
# Helpers
# ----------------------------
mat_path <- function(type, snp) file.path(paste0(type, "MAT"), paste0(snp, ".matrix.txt"))

read_mat <- function(path) {
  if (!file.exists(path)) return(NULL)
  
  hdr <- names(fread(path, nrows = 0))
  hdr <- hdr %>%
    gsub("chrom", "chromosome", ., fixed = TRUE) %>%
    gsub("pos",   "coordinate", ., fixed = TRUE)
  
  mat <- fread(path, skip = 1)
  if (ncol(mat) != length(hdr)) stop("Header/data mismatch: ", path)
  setnames(mat, hdr)
  
  if (!all(c("chromosome","coordinate") %in% names(mat))) stop("Missing chromosome/coordinate in ", path)
  
  as_tibble(mat) %>% mutate(coordinate = as.integer(coordinate))
}

allele_split_line <- function(mat_tbl, ref, alt) {
  sample_cols <- setdiff(names(mat_tbl), c("chromosome","coordinate"))
  allele_labels <- gsub("_[A-Za-z0-9]+$", "", sample_cols)
  sizes <- table(allele_labels)
  
  if (!(ref %in% names(sizes)) || !(alt %in% names(sizes))) return(NA_real_)
  
  n_ref <- as.integer(sizes[[ref]])
  n_alt <- as.integer(sizes[[alt]])
  if (ref < alt) n_ref + 0.5 else n_alt + 0.5
}

nearest_coord <- function(coords, target) {
  coords[which.min(abs(coords - target))]
}

# Convert wide mat -> long, joined to global map
to_long <- function(mat_tbl, map_tbl, ref, property_label) {
  mat_tbl %>%
    dplyr::left_join(map_tbl, by = "coordinate") %>%
    dplyr::select(-chromosome) %>%
    tidyr::pivot_longer(
      cols = -c(coordinate, X),
      names_to = "Y",
      values_to = "Z"
    ) %>%
    dplyr::mutate(
      Z = as.numeric(Z),
      type = if_else(stringr::str_detect(Y, paste0("^", ref, "_")), "REF", "ALT"),
      property = property_label
    )
}
make_title <- function(row) {
  base <- paste0(
    row$SNP, " (", row$CHR, ":", row$SNPLOC, ") @",
    row$CHR, ":", row$L, "-", row$R, "\n",
    "ref-mean: ", sprintf("%.1f", as.numeric(row$PCTREF)), "%; ",
    "alt-mean: ", sprintf("%.1f", as.numeric(row$PCTALT)), "%; ",
    "FDR = ", format(as.numeric(row$FDR), digits = 3)
  )
  
  # optional nd / td info
  has_ndtd <- all(c("nd", "td") %in% names(row)) &&
    !is.na(row$nd) && !is.na(row$td)
  
  ndtd <- if (has_ndtd) {
    paste0("\nnd = ", row$nd, "; td = ", row$td,
           "; |Δ| = ", abs(row$td - row$nd))
  } else {
    ""
  }
  
  anno <- paste0(
    "\nGWAS|mQTL: ", row$GWAS, "|", row$mQTL
  )
  
  paste0(base, ndtd, anno)
}

# ----------------------------
# Main loop: one SNP per page (mc on top, sc on bottom)
# ----------------------------
plots <- vector("list", nrow(tnplot.table))

for (i in seq_len(nrow(tnplot.table))) {
  row <- tnplot.table[i, ]
  
  chr    <- as.character(row$CHR)
  snp    <- as.character(row$SNP)
  dmrL0  <- as.integer(row$L)
  dmrR0  <- as.integer(row$R)
  snploc <- as.integer(row$SNPLOC)
  ref    <- as.character(row$REF)
  alt    <- as.character(row$ALT)
  
  mc_file <- mat_path(type_mc, snp)
  sc_file <- mat_path(type_sc, snp)
  
  # If mc missing, try to generate it (your workflow)  
  if (!file.exists(mc_file)) {
    vstring <- paste0(snp, ":", ref, ":", alt, ":", chr, ":", snploc)
    inputbamdir <- paste0(type_mc, "Slice")
    outdir <- file.path("/share/lab_wangl/Yijun_Tian/RL_ASM", paste0(type_mc, "MAT"))
    cmd <- paste("bash", "makedamat.sh", shQuote(vstring), shQuote(inputbamdir), shQuote(outdir))
    system(cmd)
  }
  
  # If sc missing, try to generate it (your workflow)
  if (!file.exists(sc_file)) {
    vstring <- paste0(snp, ":", ref, ":", alt, ":", chr, ":", snploc)
    inputbamdir <- paste0(type_sc, "Slice")
    outdir <- file.path("/share/lab_wangl/Yijun_Tian/RL_ASM", paste0(type_sc, "MAT"))
    cmd <- paste("bash", "makedamat.sh", shQuote(vstring), shQuote(inputbamdir), shQuote(outdir))
    system(cmd)
  }
  
  mat_mc <- read_mat(mc_file)
  mat_sc <- read_mat(sc_file)
  
  if (is.null(mat_mc) && is.null(mat_sc)) next
  if (is.null(mat_mc) || is.null(mat_sc)) next  # for your “both panels per SNP” requirement
  
  # shared coordinate map (union) -> X rank
  map <- bind_rows(
    mat_mc %>% dplyr::select(coordinate),
    mat_sc %>% dplyr::select(coordinate)
  ) %>%
    distinct() %>%
    arrange(coordinate) %>%
    mutate(X = rank(coordinate, ties.method = "first"))
  
  coords <- map$coordinate
  
  # robust boundary positions (nearest fallback)
  dmrL_coord <- nearest_coord(coords, dmrL0 + 1L)
  dmrR_coord <- nearest_coord(coords, dmrR0)
  snp_coord  <- nearest_coord(coords, snploc)
  
  dmrL <- map$X[match(dmrL_coord, coords)]
  dmrR <- map$X[match(dmrR_coord, coords)]
  snppos <- map$X[match(snp_coord,  coords)]
  
  # split lines (per panel)
  yitc_mc <- allele_split_line(mat_mc, ref, alt)
  yitc_sc <- allele_split_line(mat_sc, ref, alt)
  
  lines_df <- tibble(
    property = factor(c("main_character", "support_character"),
                      levels = c("main_character", "support_character")),
    Lines = c(yitc_mc, yitc_sc)
  )
  
  # long-form data with shared X
  tmp <- bind_rows(
    to_long(mat_mc, map, ref, "main_character"),
    to_long(mat_sc, map, ref, "support_character")
  )
  
  brk <- unique(round(quantile(map$X, probs = seq(0, 1, 0.2), names = FALSE)))
  
  plots[[i]] <- ggplot(tmp, aes(X, Y, fill = Z)) +
    geom_tile() +
    scale_y_discrete(name = "allelic_summary") +
    scale_x_continuous(name = "genomic coordinates",
                       breaks = brk,
                       labels = map$coordinate[match(brk, map$X)]) +
    scale_fill_gradient2(high = "red", low = "blue", midpoint = 50,
                         limits = c(0, 100), name = "%meth%") +
    geom_vline(xintercept = dmrL, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_vline(xintercept = dmrR, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_vline(xintercept = snppos + 0.5, linetype = "dashed", colour = "green", linewidth = 0.8) +
    geom_hline(data = lines_df, aes(yintercept = Lines),
               linetype = "dashed", colour = "black", linewidth = 0.8) +
    facet_grid(property ~ ., scales = "free_y", space = "free_y") +
    ggtitle(make_title(row)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      strip.text.y = element_text(face = "bold")
    )
}

plots <- plots[!vapply(plots, is.null, logical(1))]
if (length(plots) == 0) stop("No plots generated.")

ggsave(out_pdf, marrangeGrob(plots, nrow = 1, ncol = 1),
       width = 8.5, height = 11)
message("Saved: ", out_pdf)
