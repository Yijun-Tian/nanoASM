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
opt_list <- list(
  make_option(c("--annoDMR"), type="character",
              help="An annotated DMR file with allele specific methylation regions, it was cross annotated with GTEx eQTL for the plot"),
  make_option(c("--qvalue"), type = "numeric", default = 1e-5, help = "FDR criteria for regions to plot heatmap"),
  make_option(c("--delta"), type = "numeric", default = 25, help = "Methylation delta criteria for regions to plot heatmap"),
  make_option(c("--type"), type="character", help="Either Normal or Cancer")
)

opt <- parse_args(OptionParser(option_list=opt_list))
qv = opt$qvalue
d = opt$delta
type = opt$type

colnames_dmr = c("asId","CHR","L","R","FDR","V5","V6","V7","V8","PCTREF","PCTALT","SNP","V12","V13","V14","SNPLOC","V16","V17","V18","V19",
                 "GT","GENE","EVENT","dis2TSS","maf", "maSample","sample", "pvalue","NES","SE","GWAS","mQTL")

dmr <- fread(opt$annoDMR, fill = Inf, header = F)
colnames(dmr) <- colnames_dmr
summary <- dplyr::filter(dmr, FDR <= qv & abs(V5) >= d) %>%
  arrange(desc(abs(V5)))

p <- list()

for (i in 2:dim(summary)[1]){
  chr <- summary[i]$CHR
  snp <- summary[i]$SNP
  dmr <- c(summary[i]$L, summary[i]$R)
  
  snploc <- summary[i]$SNPLOC
  pct_ref <- summary[i]$PCTREF
  pct_alt <- summary[i]$PCTALT
  fdr <- summary[i]$FDR
  ref <- summary[i]$V12
  alt <- summary[i]$V13
  
  GWAS <- summary[i]$gwas
  MQTL <- summary[i]$mQTL
  
  title <- paste0(snp," (",chr,":",snploc,") @", chr, ":", dmr[1], "-", dmr[2], "\nref-mean: ", pct_ref, "%; alt-mean: ", pct_alt, "%; FDR = ",fdr, ".\n(", GWAS,"|",MQTL,")")
  colnames <- read.table(file = paste0(type,"MAT/",snp,".matrix.txt"), nrows = 1) %>%
    gsub('chrom', 'chromosome',.) %>%
    gsub('pos', 'coordinate',.)
  sizes <- gsub('_[a-zA-Z0-9]+$','', colnames[-1:-2]) %>% table()
  n_ref <- sizes[ref][[1]]
  n_alt <- sizes[alt][[1]]
  mat <- fread(file = paste0(type,"MAT/",snp,".matrix.txt"), skip = 1) %>%
    rename_with(., ~unlist(colnames), starts_with("V")) %>%
    as_tibble() %>%
    mutate(X=rank(coordinate))
  alleles <- names(sizes)
  map <- mat[,c("X","coordinate")]
  dmrL <- map[which(map$coordinate == dmr[1]+1),]$X
  dmrR <- map[which(map$coordinate == dmr[2]),]$X
  snppos <- which(abs(map$coordinate-snploc)==min(abs(map$coordinate-snploc)))
  brk <- round(quantile(1:dim(map)[1], probs = seq(0,1,0.1)))
  #horizontal coordinate switch caused by sample id ranking
  if (ref < alt) {
    yitc = n_ref + 0.5
  } else {
    yitc = n_alt + 0.5
  }
  
  p[[i]]  <- dplyr::select(mat, -c(chromosome, coordinate)) %>% gather(key = "Y", value = "Z", -X) %>% 
    mutate(type = case_when(str_detect(Y, paste0(ref,"_")) ~ "REF", TRUE ~ "ALT")) %>%
    ggplot(aes(X, Y, fill= Z)) + geom_tile() + coord_fixed(ratio = 4) +
    scale_y_discrete(name="allelic_summary") +
    scale_x_continuous(name="genomic coordinates", breaks= brk, labels = map$coordinate[brk]) +
    scale_fill_gradient2(high = "red", low= "blue", name="%meth%", limits = c(0,100), midpoint = 50) +
    geom_vline(xintercept = dmrL, linetype = "dashed", colour = "black",size = 1) +
    geom_vline(xintercept = dmrR, linetype = "dashed", colour = "black",size = 1) + 
    geom_vline(xintercept = snppos + 0.5, linetype = "dashed", colour = "green",size = 1) +
    geom_hline(yintercept = yitc, linetype = "dashed", color = "black", size = 1)+
    ggtitle(label=title) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
ggsave(paste0(type, ".HeatmapASM.pdf"), marrangeGrob(p, nrow = 2, ncol = 1), width = 8.5, height = 11)

