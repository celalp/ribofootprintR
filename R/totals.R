totals <- function(data, genedf, plot=F, cores=1) {

    summer <- function(nuc, df) {
        a <- df[df$nucleotide == nuc, ]
        a <- sum(a[["freq"]], na.rm = T)
        a
    }
    genes <- genedf$gene
    get_data <- function(gene_name) {
        gene <- get_gene(gene = gene_name, genedf=genedf, data = data, seq = NA)
        if (is.null(gene)) {
            NULL
        } else {
            coord <- genedf[genedf$gene == gene_name, ]
            gene <- gene[gene$nucleotide >= 0 & gene$nucleotide <= (coord$end - coord$start + 1), ]
            counts <- sapply(unique(gene$nucleotide), summer, gene)
            gene <- as.data.frame(cbind(unique(gene[, 1:2]), counts))
            if (dim(gene)[1] == 0) {
                NULL
            } else {
                dats <- colSums(gene[, 2:3])
                gene <- data.frame(sum_coverage = dats[1], sum_ribo = dats[2])
                rownames(gene) <- gene_name
                gene
            }
        }
    }
    if(cores>1){
      totals <- do.call("rbind", mclapply(genes, get_data, mc.cores = cores))
      if (is.null(data$rna)) {
        totals$sum_coverage <- NA
      }
      }else {
        totals <- do.call("rbind", lapply(genes, get_data))
        if (is.null(data$rna)) {
          totals$sum_coverage <- NA
        }
      }
    totals<-as.data.frame(totals)
    if(plot){
      print(ggplot(totals, aes(x=sum_coverage, y=sum_ribo))+geom_point())
    }
    invisible(totals)
}

