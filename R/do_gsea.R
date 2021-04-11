#' do_gsea
#'
#' sets up gsea
#' @param data data from DEseq2
#' @param GO_terms object from getGOterms function
#' @param gsea_lim lower and upper limit for number of genes per GO term
#' @import piano dplyr
#' @export
#'
#'
do_gsea <- function(GO_terms, data, gsea_lim = c(20,500)){

  GO_terms_gsea <- GO_terms[c("ensembl_gene_id", "name_1006")]
  names(GO_terms_gsea) <- c('gene.name', "GO_terms")
  GO_terms_gsea <- loadGSC(GO_terms_gsea)


  data_gsea <- data[c('ensembl_gene_id', 'log2FoldChange')]

  row.names(data_gsea) <- data_gsea$ensembl_gene_id
  data_gsea$ensembl_gene_id <- NULL

  message("starting gsea")
  gsea_res <- runGSA(data_gsea, gsc=GO_terms_gsea,  gsSizeLim=gsea_lim, geneSetStat="gsea")

  message("finished gsea")

  gsea_data <- GSAsummaryTable(gsea_res)

  tidy_gsea <- function(x){
    pos <- subset(x, x[3] > 0)
    neg <- subset(x, x[3] < 0)

    pos <- pos[c("Name", "Genes (tot)", "Stat (dist.dir)" , "p (dist.dir.up)" , "p adj (dist.dir.up)")]
    neg <- neg[c("Name", "Genes (tot)", "Stat (dist.dir)" , "p (dist.dir.dn)" , "p adj (dist.dir.dn)")]

    names(pos) <- c("Name", "gene.NO", 'stat', "pvalue", "padj")
    names(neg) <- c("Name", "gene.NO", 'stat', "pvalue", "padj")

    exp_data <- rbind(pos,neg)

    return(exp_data)
  }

  gsea_data <- tidy_gsea(gsea_data)

  message("adding gene names")
  gene_names <- data.frame()
  loop = 0
  for(goterm in gsea_data$Name){
    go_info <- subset(GO_terms, GO_terms$name_1006 == goterm)
    go_data <- subset(data, data$ensembl_gene_id %in% go_info$ensembl_gene_id)
    go_data <- left_join(go_data, distinct(GO_terms[c("ensembl_gene_id","external_gene_name")]), by = c("ensembl_gene_id" = "ensembl_gene_id"))



    if(gsea_data[gsea_data$Name == goterm,][3] >0){
      go_data <- go_data[order(-go_data$log2FoldChange),]
    } else {
      go_data <- go_data[order(go_data$log2FoldChange),]
    }
    exp_data <- data.frame(Name = goterm,
                           external_gene_name = paste(go_data$external_gene_name, collapse = "|"))
    gene_names <- rbind(gene_names, exp_data)
    loop = loop +1
    if(loop%%100 == 0){message(paste("processed", loop, "GO terms", sep = " "))}
  }

  gsea_data <- left_join(gsea_data, gene_names, by = c("Name"= "Name"))
  return(gsea_data)
}
