suppressPackageStartupMessages({
  library(venneuler)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(hrbrthemes)
  library(multidplyr)
  library(dplyr, warn.conflicts = FALSE)
  library(data.table)
  library(stringr)
  library(ggpmisc)
  library(ggpubr)
  library(annotatr)
  library(purrr)
  library(EnsDb.Hsapiens.v86)
  library(gridExtra)
  library(ggsankey)
  library(dgof)
  library(bedr)
  library(BioVenn)
})

# Figure 5F, overlapping between tested SNPs between normal and cancer samples.
# Check number from "Cancer|Normal.input.rsid.txt"
d <- venneuler(c(Normal = 4485180 - 2037091, Cancer = 2207008 - 2037091, "Normal&Cancer" = 2037091))
plot(d)
title("input SNP")

# Figure 5G, overlapping SNPs of significant DMR between normal and cancer samples.
# Check number from "Cancer|Normal.dmr.rsid.txt"
CancerDMR <- scan("Cancer.dmr.rsid.txt", what = "character")
NormalDMR <- scan("Normal.dmr.rsid.txt", what = "character")
DMRintersect <- intersect(CancerDMR, NormalDMR)
c <- venneuler(c(Normal = 174867 - 13648, Cancer = 112429 - 13648, "Normal&Cancer" = 13648))
plot(c)
title("DMR SNP")

# Figure 5H, overlapping SNPs and ASM-DMR at the same time between normal and cancer samples.
# Check number from "awk '{FS=OFS="\t"}{if ($4 <= 0.05 && ($5)^2 >= 100) print $11"_"$0}' Normal.dmr.txt | sort -k1,1 -k2,2n",
#                   "awk '{FS=OFS="\t"}{if ($4 <= 0.05 && ($5)^2 >= 100) print $11"_"$0}' Normal.dmr.txt | sort -k1,1 -k2,2n",
#                   and "NormalNCancer.dmr_intersect.txt"
i <- venneuler(c(Normal = 186163 - 11869, Cancer = 124033 - 11869, "Normal&Cancer" = 11869))
plot(i)
title("interval numbers q<=0.05 and delta >= 10%")
par(mfrow = c(1,3))
plot(d)
title("input SNP")

plot(c)
title("DMR SNP")

plot(i)
title("interval numbers q<=0.05 and delta >= 10%")
# ==================================================================================================================================================
# Figure 5I and 5J
transition_mut = c("A2G","G2A","C2T","T2C")
transversion_mut = c("A2C", "C2A", "G2T", "T2G")
shareR <- fread("NormalNCancer.dmr_intersect.txt", header = F) %>%
  dplyr::select(c(V1,V3,V4,V15,V16,V28,V29,V30)) %>%
  dplyr::rename(rsId = "V1", ns = "V3", ne = "V4", ts = "V15", te = "V16", sizeS = "V28", fracN = "V29", fracT = "V30")
normal <- fread("Normal.dmr.txt", header = F) %>%
  left_join(.,shareR, by = c("V11" = "rsId", "V2"= "ns", "V3" = "ne")) %>%
  mutate(dmrtype = case_when(!is.na(ts) ~ "sharedN", .default = "onlyN")) %>%
  dplyr::rename(PNs = "ts", PNe = "te")
cancer <- fread("Cancer.dmr.txt", header = F)%>%
  left_join(.,shareR, by = c("V11" = "rsId", "V2"= "ts", "V3" = "te")) %>%
  mutate(dmrtype = case_when(!is.na(ns) ~ "sharedC", .default = "onlyC"))  %>%
  dplyr::rename(PNs = "ns", PNe = "ne")
df <- rbind(normal, cancer) %>%
  left_join(., snpinfo, by = c("V11" = "rsId")) %>%
  mutate(qv = case_when(V4 <= 0.05 & abs(V5) >= 10 ~ "significant", .default = "not")) %>%
  mutate(delta = ntile(abs(V5), 4)) %>%
  mutate(refNalt = paste0(V12,"2",V13)) %>%
  mutate(RA = case_when(refNalt %in% transversion_mut ~ "transversion",
                        refNalt %in% transition_mut ~ "transition",
                        .default = refNalt)) %>%
  mutate(v2r = ntile(abs(V15 - 0.5*(V3+V2)),4)) %>%
  mutate(region_size = case_when(V3 - V2 <= 500 ~ "small", V3 - V2 > 500 & V3 - V2 <=1000 ~ "median", TRUE ~ "large")) %>%
  mutate(v2r_relationship = case_when(V18 * V19 < 0 ~ "reside", abs(V18) <= 200 | abs(V19) <=200 ~ "near", abs(V18) > 200 & abs(V19) >200 ~ "dist", TRUE ~ "NA")) %>%
  group_by(V11, dmrtype) %>%
  mutate(n_dmr_per_snp = n()) %>%
  ungroup() %>%
  mutate(g_dmr_per_snp = case_when(n_dmr_per_snp <= 3 ~ as.character(n_dmr_per_snp),
                                   .default = "morethan4")) %>%
  mutate(cpg_density = V6*100 / (V3 - V2)) %>% # cpg density is a metric to evaluate how many cpg sites has been included in this dmr
  mutate(dmr.direction = case_when(V5 > 0 ~ "rGTa", V5 < 0 ~ "rLTa"))
snp <- dplyr::select(df, c(V1, V15, V12, V13, V11)) %>% unique(.) %>% dplyr::rename(chr = "V1", ref_allele = "V12", alt_allele = "V13", name = "V11") %>%
  #mutate(start = V15 - 1, end = V15 + 2, score = ".", strand = "+") %>%
  mutate(start = V15 - 19, end = V15 + 20, score = ".", strand = "+") %>%
  data.frame(.)
saveRDS(snp, "snp.rds")

dmr <- dplyr::select(df, c(V1, V2, V3)) %>% unique(.) %>% dplyr::rename(chr = "V1", start = "V2", end = "V3") %>%
  mutate(score = ".", strand = "+") %>%
  data.frame(.)
saveRDS(dmr, "dmr.rds")
dmr.get <- dmr %>% arrange(chr, start) %>% get.fasta(fasta = "/share/lab_wangl/Yijun_Tian/WGS_Prostate/hg38.alignmt.fa", check.chr = FALSE, use.name.field = TRUE)
saveRDS(dmr.get, "dmr.get.rds")
snp.get <- dplyr::select(snp, c(chr, start, end, name, score, strand)) %>%
  arrange(chr, start) %>%
  get.fasta(fasta = "/share/lab_wangl/Yijun_Tian/WGS_Prostate/hg38.alignmt.fa",
            check.chr = FALSE, use.name.field = TRUE)
snp.get <- snp.get %>% separate(index, into = c("rsId", "chr", "string_start", "string_end"), sep = "::|:|-") %>%
  dplyr::select(c(rsId, sequence)) %>% mutate(sequence = toupper(sequence))
snp.seqinfo <- left_join(snp, snp.get, by = c("name" = "rsId")) %>%
  rowwise() %>%
  mutate(refstring = case_when(substr(sequence, 20, 20) == ref_allele ~ sequence, .default = NA)) %>%
  mutate(altstring = case_when(!is.na(refstring) ~ paste0(substr(sequence, 1, 19), alt_allele, substr(sequence, 21, 39)), .default = NA)) %>%
  mutate(refseq = substr(refstring, 19, 21)) %>%
  mutate(altseq = substr(altstring, 19, 21)) %>%
  mutate(cg.type = case_when(str_detect(refseq, "CG") & !str_detect(altseq, "CG") ~ "wiping",
                             !str_detect(refseq, "CG") & str_detect(altseq, "CG") ~ "forming",
                             str_detect(refseq, "CG") & str_detect(altseq, "CG") ~ "cooccur",
                             !str_detect(refseq, "CG") & !str_detect(altseq, "CG") ~ "nonexist"))
tmp <- left_join(df, snp.seqinfo, by = c("V1" = "chr", "V11" = "name", "V12" = "ref_allele", "V13" = "alt_allele", "V15" = "V15")) %>%
  # dplyr::filter(V4 <= 0.1) %>%
  # dplyr::filter(v2r_relationship %in% c("reside", "near")) %>%
  mutate(hmSEQ = case_when(V5 > 0 ~ refstring, V5 < 0 ~ altstring)) %>%
  mutate(lmSEQ = case_when(V5 > 0 ~ altstring, V5 < 0 ~ refstring))
df <- tmp
# Figure 5I
h10_data = dplyr::filter(df, qv == "significant") %>% dplyr::select(V1, V2, V3, V11, v2r_relationship, dmrtype, region_size) %>% unique(.) %>%
  mutate(sample = gsub("[a-z]+", "", dmrtype)) %>% mutate(category = gsub("[A-Z]+", "", dmrtype)) %>% mutate(size = log2(V3-V2))
h10_summary <- h10_data %>% group_by(sample, category) %>%
  summarise(Q1 = quantile(size, 0.25, prob = 0.25),
            Q2 = quantile(size, 0.50, prob = 0.50),
            Q3 = quantile(size, 0.75, prob = 0.75))
h10 <-  ggplot(h10_data, aes(x = size, color = factor(dmrtype), fill = factor(dmrtype))) +
  geom_histogram(aes(y = ..density..), bins = 100, position = "identity", alpha = 0.5) +
  geom_density(alpha = 0.6) +
  facet_grid(cols = vars(category), rows = vars(sample)) +
  geom_vline(data = h10_summary, aes(xintercept = Q1), color = "red", linetype="dashed", linewidth = 1) +
  geom_vline(data = h10_summary, aes(xintercept = Q2), color = "blue", linetype="dashed", linewidth = 1) +
  geom_vline(data = h10_summary, aes(xintercept = Q3), color = "green", linetype="dashed", linewidth = 1) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 8, face = "bold"),
        legend.position="none") +
  ggtitle("DMR type and size, KS-test significant between normal and cancer")

ggarrange(p10, h10, labels = c("A", "B"), ncol = 1, nrow = 2)

ks.test(h10_data[which(h10_data$dmrtype == "onlyC"),]$size, h10_data[which(h10_data$dmrtype == "onlyN"),]$size)
ks.test(h10_data[which(h10_data$dmrtype == "sharedC"),]$size, h10_data[which(h10_data$dmrtype == "sharedN"),]$size)

# Figure 5J
h5_data <- dplyr::filter(df, qv == "significant")  %>% dplyr::select(V1, V2, V3, V9, V10, V11, dmrtype, delta) %>% unique(.) %>%
  mutate(sample = gsub("[a-z]+", "", dmrtype)) %>% mutate(category = gsub("[A-Z]+", "", dmrtype)) %>% mutate(D = abs(V9-V10))
h5_summary <- h5_data %>% group_by(sample, category) %>%
  summarise(Q1 = quantile(D, 0.25, prob = 0.25),
            Q2 = quantile(D, 0.50, prob = 0.50),
            Q3 = quantile(D, 0.75, prob = 0.75))
h5 <-  ggplot(h5_data, aes(x = D, color = factor(dmrtype), fill = factor(dmrtype))) +
  geom_histogram(aes(y = ..density..), bins = 100, position = "identity", alpha = 0.5) +
  geom_density(alpha = 0.6) +
  facet_grid(cols = vars(category), rows = vars(sample)) +
  geom_vline(data = h5_summary, aes(xintercept = Q1), color = "red", linetype="dashed", linewidth = 1) +
  geom_vline(data = h5_summary, aes(xintercept = Q2), color = "blue", linetype="dashed", linewidth = 1) +
  geom_vline(data = h5_summary, aes(xintercept = Q3), color = "green", linetype="dashed", linewidth = 1) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 8, face = "bold"),
        legend.position="none") +
  ggtitle("DMR type and methylation delta, KS-test significant between N and C")
ggarrange(p5, h5, labels = c("A", "B"), ncol = 1, nrow = 2)
ks.test(h5_data[which(h5_data$dmrtype == "onlyC"),]$D, h5_data[which(h5_data$dmrtype == "onlyN"),]$D)
ks.test(h5_data[which(h5_data$dmrtype == "sharedC"),]$D, h5_data[which(h5_data$dmrtype == "sharedN"),]$D)

# Figure 5K
data <- df %>% dplyr::filter(qv == "significant") %>%
  dplyr::filter(cg.type %in% c("forming", "wiping")) %>%
  group_by(dmrtype, v2r_relationship) %>%
  group_modify(~make_long(.x, c(cg.type, dmr.direction)))

p13 = ggplot(data, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(aes(x = as.numeric(x) + .05,label = after_stat(paste0(node, "\nn = ", freq)))) +
  theme_sankey(base_size = 8) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 8, face = "bold"),
        legend.position="none") +
  facet_grid(rows = vars(v2r_relationship), cols = vars(dmrtype), scales = "free_y") +
  ggtitle("cpg form/wipe type")
# Below used to calculate the pct in the sankey plot
info <- df %>% dplyr::filter(qv == "significant") %>%
  dplyr::filter(cg.type %in% c("forming", "wiping"))
table(info$cg.type, info$dmr.direction, info$v2r_relationship,info$dmrtype)
