totals <- function(data, genedf, plot=F, cores=1) {

  genes <- genedf$gene
  get_data <- function(gene_name) {
    gene <- get_gene(gene = gene_name, genedf=genedf, data = data, seq = NA)
    if (is.null(gene)) {
      NULL
    } else {
      coord <- genedf[genedf$gene == gene_name, ]
      gene <- gene[gene$nucleotide >= 0 & gene$nucleotide <= (coord$end - coord$start + 1), ]
      if (dim(gene)[1] == 0) {
        NULL
      } else {
        gene$frame<-as.numeric(gene$frame)
        dats <- colSums(gene, na.rm = T)
        dats<-as.data.frame(t(dats))
        dats$gene<-gene_name
        return(dats)
      }
    }
  }
  if(cores>1){
    tots <- do.call("rbind", mclapply(genes, get_data, mc.cores = cores))
  }else {
    tots <- do.call("rbind", lapply(genes, get_data))
  }
  if(is.null(data$rna)){
    totals_df<-data.frame(gene=tots$gene, sum_ribo=tots$freq, sum_coverage=NA)
  } else {
    totals_df<-data.frame(gene=tots$gene, sum_ribo=tots$freq, sum_coverage=tots$coverage)
  }

  totals<-as.data.frame(totals)
  if(plot){
    print(ggplot(totals_df, aes(x=sum_coverage, y=sum_ribo))+geom_point())
  }
  invisible(totals_df)
}

