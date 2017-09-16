metagene_bins <- function(data, genedf, bins = 100, norms = "rna", simplify = "mean", cores=1, plot=F) {

    summer <- function(nuc, df, measure) {
        a <- df[df$nucleotide == nuc, ]
        a <- sum(a[[measure]], na.rm = T)
        a
    }

    if (norms == "rna" & is.null(data[["rna"]])) {
        warning("No RNA-Seq Data Setting norms to total")
        norms = "total"
    }

    #quick and dirty to fix promise issues
    gene_list<-genedf

    norm_gene <- function(name, gene_list) {
      gene <- get_gene(gene = name, data = data, seq = NA, genedf=gene_list)
        if(is.null(gene) || dim(gene)[1]<1) {
            NULL
        } else {
            coord <- genedf[genedf$gene == name, ]
            gene <- gene[gene$nucleotide >= 0 & gene$nucleotide <= (coord$end - coord$start + 1), ]
            if(dim(gene)[1]<1){
              NULL
            } else {
            counts <- sapply(unique(gene$nucleotide), summer, gene, "freq")
            coverage <- sapply(unique(gene$nucleotide), summer, gene, "coverage")

            if (norms == "none") {
                a <- counts
                a
            } else if (norms == "rna") {
                a <- counts/coverage
                a
            } else if (norms == "total") {
                a <- counts/sum(counts, na.rm = T)
                a
            }
        }
        }
    }

    binner <- function(gene_name) {
        gene <- norm_gene(name=gene_name, gene_list=gene_list)
        gene <- noinf(gene)
        gene_bins <- 0
        for (i in 1:bins) {
            lower <- floor((length(gene)/bins) * (i - 1))
            upper <- ceiling((length(gene)/bins) * (i))
            gene_bins[i] <- mean(gene[lower:upper], na.rm = T)
        }
        if (is.null(gene_bins)) {
            gene_bins <- rep(NA, bins)
            gene_bins
        } else {
            gene_bins
        }
    }

    #noinf <- function(a) {
    #    loc <- is.infinite(a)
    #    loc2 <- is.nan(a)
    #    a <- a[!loc]
    #    a <- a[!loc2]
    #}

    #supress error when the number of nucleotides is not divisible by # of bins
    genes <- genedf$gene
    if(cores>1){
      suppressWarnings(binned <- do.call("rbind", mclapply(genes, binner, mc.cores = cores, mc.cleanup = T)))
    } else {
      suppressWarnings(binned <- do.call("rbind", lapply(genes, binner)))
    }

    rownames(binned) <- genes
    #binned <- apply(binned, 2, noinf)

    if (simplify == "mean") {
        binned <- data.frame(bin=1:bins, value=apply(binned, 2, mean, na.rm = T))
        print(ggplot(binned, aes(x=bin, y=value))+geom_line())

    } else if (simplify == "median") {
        binned <- data.frame(bin=1:bins, value=apply(binned, 2, median, na.rm = T))
        print(ggplot(binned, aes(x=bin, y=value))+geom_line())

    } else if (simplify == "none") {
        binned <- reshape2::melt(binned)
        colnames(binned) <- c("gene", "bin", "value")
        if(plot){
          print(ggplot(binned, aes(x=bin, y=value))+geom_smooth())
        }
    }
    invisible(binned)
}
