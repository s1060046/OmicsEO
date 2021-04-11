#' analyse_fmrp
#'
#' perform FMRP target analysis on DESeq result
#' @param input DESeq result
#' @param input plot choose between "comparison" or "validation". comparison for fmrp target analysis
#' @param cds_length_val cds length threshold default 1500
#' @param baseMean_val abundance threshold default 1500
#' @import ggplot2
#' @import ggpubr
#' @export
#'
#'
analyse_fmrp <- function(input, plot = c("comparison", "validation"), cds_length_val = 1500, baseMean_val = 1500){
  basic_info <- readRDS(file = "./raw_data/gene_list/detected_gene_info_unique_fmrp.rds")

  plotting_data <- left_join(x, basic_info, by = c("ensembl_gene_id"))
  message(paste("original data ", dim(x)[1], " after merging ",  dim(plotting_data)[1], sep=""))

  plotting_data <- plotting_data[!is.na(plotting_data$group),]
  plotting_data_ <- subset(plotting_data, plotting_data$group == "non-target")

  plotting_data_fmrp <- subset(plotting_data, plotting_data$group != "non-target")
  plotting_data_fmrp$group = "target"

  plotting_data_long <- subset(plotting_data_, plotting_data_$cds_length >cds_length_val)
  plotting_data_long <- subset(plotting_data_long, plotting_data_long$baseMean >baseMean_val)

  plotting_data_long$group <- "long"
  print(dim(plotting_data_long)[1])
  plotting_data <- rbind(plotting_data_fmrp, plotting_data_, plotting_data_long)
  fmr_long <- wilcox.test(plotting_data_fmrp$stat, plotting_data_long$stat)
  fmr_all <- wilcox.test(plotting_data_fmrp$stat, plotting_data_$stat)
  message("results of wilcoxon test with BH pvalue adjustment")
  message(data.frame(test = c("Fmrp target vs long", "Fmrp target vs all"),
                   p = c(fmr_long$p.value, fmr_all$p.value),
                   padj = p.adjust(c(fmr_long$p.value, fmr_all$p.value), method = "BH")))

  if(plot == "comparison"){
    ggplot(plotting_data, aes(stat, color = group, group = group)) + stat_ecdf() + theme_bw() + scale_color_manual(values = c("Dark green", "Black", "Red"))
  } else if(plot == "validation"){
    a <- ggplot(plotting_data, aes(log10(baseMean), color = group)) + geom_density()
    b <- ggplot(plotting_data, aes(log2FoldChange, color = group)) + geom_density()
    ggarrange(a,b)
  }
}
