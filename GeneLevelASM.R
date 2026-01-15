suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(hrbrthemes)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggpmisc)
  library(ggpubr)
  library(annotatr)
  library(purrr)
  library(EnsDb.Hsapiens.v86)
  library(gridExtra)
  library(tidyverse)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(optparse)
})
opt_list <- list(
  make_option(c("--annowindow"), type="character",
                help="An annotated window file with merged allele specific DMR regions, it was cross annotated with GTEx eQTL for the plot"),
  make_option(c("--type"), type="character", help="Either Normal or Cancer")
)
opt <- parse_args(OptionParser(option_list=opt_list))
type = opt$type

colnames_window = c("asId","windowCHR","windowL","windowR","CHR","L","R","FDR","V8","V9","V10","V11","PCTREF","PCTALT","SNP","V15","V16","V17","SNPLOC","V19","V20","V21","V22","V23","V24","n","sd_region","mean_region","GT","GENE","EVENT","dis2TSS","maf", "maSample","sample", "pvalue","NES","SE")

windows <- fread(opt$annowindow, fill = Inf)
colnames(windows) <- colnames_window
gene <- dplyr::filter(windows, FDR<= 1e-5 & abs(V8) >= 25) %>% arrange(-abs(V8))
genelist <- unique(gene$GENE)

ensg <- AnnotationDbi:::select(EnsDb.Hsapiens.v86, keys =  genelist, keytype = "GENEID",
                               columns = c("TXID", "TXNAME", "GENEID", "GENENAME", "GENEBIOTYPE", "TXBIOTYPE",
                                           "ENTREZID", "GENESEQSTART", "GENESEQEND", "TXSEQSTART","TXSEQEND"))
ensg <- ensg[!duplicated(ensg[,1]),]
list <- genelist[!is.na(genelist)]
p <- list()
for (i in 1:length(list)){
  interest = list[[i]]
  block <- windows[which(windows$GENE == interest),] %>%
    dplyr::select(-c(GT, EVENT)) %>%
    unique(.) %>%
    group_by(windowCHR,windowL,windowR) %>%
    mutate(ys = rank(SNPLOC-(L+R)/2, ties.method = "first")) %>%
    mutate(ye = ys + 1) %>%
    ungroup() %>%
    mutate(eQTL = case_when(pvalue <= 0.05 ~ "significant", .default = "no")) %>%
    mutate(highlight = case_when(eQTL == "significant" & abs(V8) >= 25 ~ paste0(SNP," (",CHR,":",SNPLOC,") at", CHR, ":", L, "-", R),
                                 .default = "no")) 
  
  gs = unique(ensg[which(ensg$GENEID == interest),]$GENESEQSTART)
  ge = unique(ensg[which(ensg$GENEID == interest),]$GENESEQEND)
  ts = unique(ensg[which(ensg$GENEID == interest),]$TXSEQSTART)
  te = unique(ensg[which(ensg$GENEID == interest),]$TXSEQEND)
  
  p[[i]] <- ggplot(block, aes(xmin = L, xmax = R, ymin = ys, ymax = ye, fill = abs(V8))) +
    scale_fill_gradient2(high = "red", low= "blue", name="delta_meth", limits = c(0,100), midpoint = 25) +
    geom_rect() +
    geom_vline(xintercept = gs, colour = "green", linetype = "solid", linewidth = 0.75) +
    geom_vline(xintercept = ge, colour = "orange", linetype = "solid", linewidth = 0.75) +
    geom_vline(xintercept = ts, colour = "green", linetype = "longdash", linewidth = 0.5) +
    geom_vline(xintercept = te, colour = "orange", linetype = "longdash", linewidth = 0.5) +
    geom_segment(aes(x = SNPLOC, y = ys, yend = ye, color = factor(eQTL))) +
    geom_label(aes(x = SNPLOC, y = (ye + ys)/2, label = highlight), size = 1, color = "black", fill = NA, hjust = -0.5) + 
    scale_colour_manual(values = c("significant" = "red", "no" = "blue")) +
    facet_grid(rows = vars(paste(windowCHR,windowL,windowR, sep = "_")),scale= "free_y") +
    ggtitle(interest)
}
ggsave(paste0(type, ".GeneLevelASM.pdf"), marrangeGrob(p, nrow = 2, ncol = 1), width = 8.5, height = 11)