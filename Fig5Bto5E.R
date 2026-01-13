library(ggplot2)
library(ggpubr)
library(readr)
library(tidyr)
library(patchwork)
library(ggrastr)
library(scales)
library(IRanges)
library(plotly)
library(htmlwidgets)
library(GenomicRanges)
library(dplyr)

# Load and annotate input
normal <- read_tsv("Normal.dmr.txt", col_names = FALSE) %>% mutate(condition = "Normal")
cancer <- read_tsv("Cancer.dmr.txt", col_names = FALSE) %>% mutate(condition = "Cancer")
df_all <- bind_rows(normal, cancer) %>%  mutate(dmr_id = row_number())  # Unique ID per DMR row

# Rename columns
colnames(df_all)[c(1, 4:11, 15:19)] <- c(
  "chr", "qval", "meth_diff", "n_cpg", "pval_MWU", "pval_KS",
  "meth_pct_ref", "meth_pct_alt", "rsID", "snp_pos", "n_ref", "n_alt",
  "dist_to_start", "dist_to_end"
)

# Feature engineering
df_all <- df_all %>%
  mutate(
    significant = (qval < 0.05 & abs(meth_diff) >= 10),
    abs_meth_diff = abs(meth_diff),
    dmr_size = abs(dist_to_end - dist_to_start),
    snp_to_dmr = abs(snp_pos - (X2 + X3) / 2),
    pos = snp_pos,
    AF = n_alt / (n_ref + n_alt),
    cpg_density = n_cpg / (dmr_size + 1)
  )

# Theme
theme_pub <- theme_bw(base_size = 8) +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        plot.title = element_text(size = 10, face = "bold"))



# Assign unique DMR id
df_sig <- df_all %>%
  dplyr::filter(significant)

# Construct GRanges with rsID as pseudo-chromosome
gr_normal <- GRanges(
  seqnames = df_sig$rsID[df_sig$condition == "Normal"],
  ranges = IRanges(df_sig$X2[df_sig$condition == "Normal"],
                   df_sig$X3[df_sig$condition == "Normal"]),
  dmr_id = df_sig$dmr_id[df_sig$condition == "Normal"]
)

gr_cancer <- GRanges(
  seqnames = df_sig$rsID[df_sig$condition == "Cancer"],
  ranges = IRanges(df_sig$X2[df_sig$condition == "Cancer"],
                   df_sig$X3[df_sig$condition == "Cancer"]),
  dmr_id = df_sig$dmr_id[df_sig$condition == "Cancer"]
)

# Get overlaps
hits <- findOverlaps(gr_normal, gr_cancer)

# Calculate overlap proportions for all overlapping DMRs
overlap_table <- tibble(
  normal_id = mcols(gr_normal)$dmr_id[queryHits(hits)],
  cancer_id = mcols(gr_cancer)$dmr_id[subjectHits(hits)],
  width_inter = width(pintersect(gr_normal[queryHits(hits)], gr_cancer[subjectHits(hits)])),
  width_union = width(punion(gr_normal[queryHits(hits)], gr_cancer[subjectHits(hits)]))
) %>%
  mutate(overlap_prop = width_inter / width_union)

# Keep the max proportion per DMR (normal and cancer)
normal_overlap <- overlap_table %>%
  group_by(normal_id) %>%
  summarise(overlap_prop = max(overlap_prop), .groups = "drop")

cancer_overlap <- overlap_table %>%
  group_by(cancer_id) %>%
  summarise(overlap_prop = max(overlap_prop), .groups = "drop")

# Join back to assign shared_status
df_sig <- df_sig %>%
  left_join(normal_overlap, by = c("dmr_id" = "normal_id")) %>%
  left_join(cancer_overlap, by = c("dmr_id" = "cancer_id")) %>%
  mutate(overlap_prop = coalesce(overlap_prop.x, overlap_prop.y)) %>%
  mutate(shared_status = case_when(
    !is.na(overlap_prop) ~ "Shared",
    condition == "Normal" ~ "Normal-only",
    condition == "Cancer" ~ "Cancer-only"
  )) %>%
  select(-overlap_prop.x, -overlap_prop.y)



# Step 1: Define chromosome lengths (hg38 example)
chr_lengths <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
  chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
  chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
  chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
  chr21 = 46709983, chr22 = 50818468
)

# Step 2: Add spacing between chromosomes (e.g., 5 Mb)
chr_spacing <- 5e7
chr_lengths_spaced <- chr_lengths + chr_spacing
chr_offsets_spaced <- cumsum(c(0, head(chr_lengths_spaced, -1)))
names(chr_offsets_spaced) <- names(chr_lengths)

# Step 3: Filter and mutate dataframe
df_sig <- df_sig %>% filter(chr != "chrX") %>% filter(chr %in% names(chr_lengths))
df_sig$chr <- factor(df_sig$chr, levels = names(chr_lengths))

df_plotly <- df_sig %>%
  #for mahattan plot, need to filter out more dots for better visualization, otherwise too many dots...
  #dplyr::filter(qval < 0.05 & abs(meth_diff) >= 10) %>%
  dplyr::filter(qval < 1e-5 & abs(meth_diff) >= 10) %>%
  mutate(
    global_pos = pos + chr_offsets_spaced[as.character(chr)],
    hover = paste0("rsID: ", rsID,
                   "<br>Chr: ", chr,
                   "<br>Pos: ", pos,
                   "<br>ΔMeth: ", round(meth_diff, 2),
                   "<br>|ΔMeth|: ", round(abs_meth_diff, 2),
                   "<br>Status: ", shared_status,
                   "<br>Condition: ", condition)
  )

# Step 4: Calculate center for chromosome label ticks
chr_midpoints <- chr_offsets_spaced + chr_lengths / 2

# Step 5: Generate plot
# p: Figure 5B and 5C
p <- ggplot(df_plotly, aes(x = global_pos, y = abs_meth_diff, color = shared_status, text = hover)) +
  geom_point(alpha = 0.6, size = 0.5) +  scale_x_continuous(breaks = chr_midpoints, labels = names(chr_lengths)) +
  facet_wrap(~ condition, ncol = 1, scales = "fixed") +
  labs(x = "Chromosome", y = "|ΔMeth|", color = "Shared Status") +
  theme_pub + ggtitle("Interactive Manhattan Plot with Chromosome Gaps")
#===========================================================================================================================================
check_allele_match <- function(mqtl_ref, mqtl_alt, dmr_ref, dmr_alt) {
  # reverse complement lookup
  rc <- function(base) chartr("ACGT", "TGCA", base)
  
  mqtl_pair <- c(mqtl_ref, mqtl_alt)
  dmr_pair   <- c(dmr_ref, dmr_alt)
  
  # check direct and flipped match
  if (all(mqtl_pair == dmr_pair)) {
    return("Match")
  } else if (all(rev(mqtl_pair) == dmr_pair)) {
    return("Flipped")
  } else if (all(rc(mqtl_pair) == dmr_pair)) {
    return("Strand_RC_Match")
  } else if (all(rev(rc(mqtl_pair)) == dmr_pair)) {
    return("Strand_RC_Flipped")
  } else {
    return("Mismatch")
  }
}    

library(readr)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(liftOver)
# Step 1: Read tumor mQTL file
mqtl <- read_tsv("PRAD_tumor.cis_meQTL.txt", col_types = cols())
# need liftover mqtl probe position to hg38 before overlapping
chain_path <- "hg19ToHg38.over.chain"
ch <- import.chain(chain_path)
# ---- parse hg19 "chr:pos" to GRanges ----
pos_df <- mqtl %>%
  distinct(probe_position) %>%
  transmute(
    probe_position,
    chr = sub(":.*$", "", probe_position),
    pos = as.integer(sub("^.*:", "", probe_position))
  )

gr19 <- GRanges(
  seqnames = pos_df$chr,
  ranges   = IRanges(start = pos_df$pos, end = pos_df$pos),
  probe_position = pos_df$probe_position
)

# ---- liftover ----
lifted <- liftOver(gr19, ch)
nmap <- lengths(lifted)

# build mapping table
map_tbl <- tibble(
  probe_position = mcols(gr19)$probe_position,
  liftover_n = nmap,
  liftover_status = case_when(
    liftover_n == 0 ~ "unmapped",
    liftover_n == 1 ~ "unique",
    TRUE            ~ "multi"
  ),
  probe_position_hg38 = NA_character_
)

# fill unique mappings
idx <- which(nmap == 1)
gr38_unique <- unlist(lifted[idx], use.names = FALSE)

map_tbl$probe_position_hg38[idx] <-
  paste0(as.character(seqnames(gr38_unique)), ":", start(gr38_unique))

# (optional) split chr/pos for hg38
map_tbl <- map_tbl %>%
  mutate(
    chr_hg38 = if_else(liftover_status == "unique", sub(":.*$", "", probe_position_hg38), NA_character_),
    pos_hg38 = if_else(liftover_status == "unique", as.integer(sub("^.*:", "", probe_position_hg38)), NA_integer_)
  )

# ---- join back ----
mqtl_hg38 <- mqtl %>%
  left_join(map_tbl, by = "probe_position")

# quick sanity check
table(mqtl_hg38$liftover_status, useNA = "ifany")

mqtl_hg38_unique <- mqtl_hg38 %>% filter(liftover_status == "unique")

# Step 2: Subset cancer-specific rows
df_cancer <- df_all %>%
  filter(condition == "Cancer" & qval < 0.05)

# Step 3: SNP overlap (rsID)

df_cancer_mqtl_overlap <- df_cancer %>%
  inner_join(., mqtl_hg38_unique, by = c("rsID" = "snp")) %>%
  filter(chr == sub(":.*", "", probe_position_hg38)) %>%
  mutate(
    probe_pos = as.integer(sub(".*:", "", probe_position_hg38)),
    cpg_in_dmr = probe_pos >= X2 & probe_pos <= X3
  ) %>%
  filter(cpg_in_dmr)
cat("Number of overlapping SNPs between df_cancer and mQTL:", nrow(df_cancer_mqtl_overlap), "\n")

# Combine

overlap_df <- df_cancer_mqtl_overlap %>%
  separate(alleles, into = c("mqtl_ref", "mqtl_alt"), sep = "/", remove = FALSE) %>%
  rowwise() %>%
  mutate(
    allele_match = check_allele_match(mqtl_ref, mqtl_alt, X12, X13)
  ) %>%
  ungroup() %>%
  mutate(
    meth_diff_mqtl_style = -meth_diff,
    beta_corrected = case_when(
      allele_match %in% c("Flipped", "Strand_RC_Flipped") ~ -beta,
      TRUE ~ beta
    )
  )
cat("Number of mQTL CpGs inside DMRs:", nrow(overlap_df), "\n")

# View or save
head(overlap_df)
# write_tsv(overlap_df, "cancer_dmr_mqtl_overlap.tsv")
# Scatter plot: DMR meth_diff vs mQTL beta

# plot_scatter_cancer: Figure 5D
plot_scatter_cancer <- overlap_df %>% mutate(hover = paste0(rsID,"_",probes)) %>%
  ggplot(aes(x = meth_diff_mqtl_style, y = beta_corrected)) +
  stat_cor(method = "spearman", size = 2.5, label.x.npc = "middle", label.y.npc = "top") +
  geom_point(alpha = 0.3, size = 0.5, color = "purple") +
  scale_x_continuous(limits = c(-80, 80), breaks = seq(-80, 80, by = 10)) +  
  scale_y_continuous(limits= c( -2, 2), breaks = seq(-2, 2, by = 0.5)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = "Cancer", x = "delta", y = "beta") +
  theme_pub

# Convert to Plotly and save (optional)
# ggplotly(plot_scatter, tooltip = "text")

# Step 2: Subset cancer-specific rows
df_normal <- df_all %>%
  filter(condition == "Normal" & qval < 0.05)

# Step 3: SNP overlap (rsID)

df_normal_mqtl_overlap <- df_normal %>%
  inner_join(., mqtl_hg38_unique, by = c("rsID" = "snp")) %>%
  filter(chr == sub(":.*", "", probe_position_hg38)) %>%
  mutate(
    probe_pos = as.integer(sub(".*:", "", probe_position_hg38)),
    cpg_in_dmr = probe_pos >= X2 & probe_pos <= X3
  ) %>%
  filter(cpg_in_dmr)

cat("Number of overlapping SNPs between df_normal and mQTL:", nrow(df_normal_mqtl_overlap), "\n")

# Combine

overlap_df <- df_normal_mqtl_overlap %>%
  separate(alleles, into = c("mqtl_ref", "mqtl_alt"), sep = "/", remove = FALSE) %>%
  rowwise() %>%
  mutate(
    allele_match = check_allele_match(mqtl_ref, mqtl_alt, X12, X13)
  ) %>%
  ungroup() %>%
  mutate(
    meth_diff_mqtl_style = -meth_diff,
    beta_corrected = case_when(
      allele_match %in% c("Flipped", "Strand_RC_Flipped") ~ -beta,
      TRUE ~ beta
    )
  )
cat("Number of mQTL CpGs inside DMRs:", nrow(overlap_df), "\n")

# View or save
head(overlap_df)
# write_tsv(overlap_df, "cancer_dmr_mqtl_overlap.tsv")
# Scatter plot: DMR meth_diff vs mQTL beta

# plot_scatter_normal: Figure 5E

plot_scatter_normal <- overlap_df %>% mutate(hover = paste0(rsID,"_",probes)) %>%
  ggplot(aes(x = meth_diff_mqtl_style, y = beta_corrected)) +
  stat_cor(method = "spearman", size = 2.5, label.x.npc = "middle", label.y.npc = "top") +
  geom_point(alpha = 0.3, size = 0.5, color = "purple") +
  scale_x_continuous(limits = c(-80, 80), breaks = seq(-80, 80, by = 10)) +  
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 0.5)) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = "Normal", x = "delta", y = "beta") +
  theme_pub

layout <- c(area(1,1,2,3),  area(1,4), area(2,4))

final_plot <- (p + plot_scatter_cancer + plot_scatter_normal)+  plot_layout(design = layout, guides = "collect") & theme(legend.position = "bottom")

ggsave("Figure5B-5E_jan26_fdr1e5.pdf", final_plot, width = 6.5, height = 4.5)