#' network
#'
#'
#' @param threshold padjust threshold for considering node
#' @param direction positive terms or negative terms
#' @param selection pattern selection for TPM calculation. defult = NULL. set it to "WT" to only use WT data for
#' @import igraph
#' @export
#'
#'

networks <- function(data,full,direc, threshold){
  nodes <- get_nodes(data, direc, threshold)
  edge <- get_edge(data,full,direc, threshold)
  net_IP <- graph_from_data_frame(d=edge[!duplicated(apply(edge[1:2],1,function(x) paste(sort(x),collapse=''))),],
                                  vertices=nodes, directed=FALSE)

  V(net_IP)$size <- V(net_IP)$pvalue *5
  V(net_IP)$frame.color <- "white"

  if(direc == "pos"){V(net_IP)$color <- make_Col_pos(net_IP)}
  else if(direc == "neg"){V(net_IP)$color <- make_Col_neg(net_IP)}

  E(net_IP)$arrow.mode <- 0
  E(net_IP)$width <- E(net_IP)$weight/500
  l_IP <- layout_with_fr(net_IP)
  plot(net_IP, layout=l_IP)
}
get_sig_neg <- function(x, threshold){
  dat <- subset(x, x$padj <threshold)
  dat <- subset(dat, dat$stat <0)
  return(dat)}
get_sig_pos <- function(x, threshold){
  dat <- subset(x, x$padj <threshold)
  dat <- subset(dat, dat$stat >0)
  return(dat)}
make_Col_pos <- function(x){
  data <- data.frame(pval=V(x)$pvalue)
  pal_col <- colorRampPalette(c("yellow","red"))
  exp_table <- pal_col(20)[as.numeric(cut(data$pval,breaks = 20))]
  return(exp_table)
}
make_Col_neg <- function(x){
  data <- data.frame(pval=V(x)$pvalue)
  pal_col <- colorRampPalette(c("white","blue"))
  exp_table <- pal_col(20)[as.numeric(cut(data$pval,breaks = 20))]
  return(exp_table)
}
get_nodes <- function(x, direc, threshold){
  if(!(direc %in% c("pos", "neg"))){stop("set direction as pos or neg with quotation")}
  else if(direc == "pos"){data <- get_sig_pos(x, threshold)}
  else if(direc == "neg"){data <- get_sig_neg(x, threshold)}

  data$pvalue <- -log10(data$pvalue)
  return(data)
}
get_edge <- function(x, full,direc, threshold){
  if(!(direc %in% c("pos", "neg"))){stop("set direction as pos or neg with quotation")}
  else if(direc == "pos"){data_ <- get_sig_pos(x, threshold)}
  else if(direc == "neg"){data_ <- get_sig_neg(x, threshold)}

  full_data <- full
  name <- data_$Name


  exp <- NULL

  for(n in name){
    first <- subset(full_data, full_data$name_1006 == n)
    for(s in name){
      if(n!=s){
        second <- subset(full_data, full_data$name_1006 == s)
        ovlap <- subset(first, first$ensembl_gene_id %in% second$ensembl_gene_id)
        overlap <- dim(ovlap)[1]
        start <- first$name_1006[1]
        end <- second$name_1006[1]

        data <- data.frame(start = start, end= end, weight = overlap)

        exp <- rbind(exp,data)
      }
    }
  }
  exp <- as.data.frame(exp)
  exp <- subset(exp, exp$weight !=0)
  return(exp)
}

