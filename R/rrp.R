rrp<-function(data, genedf, norms="total", cores=1, plot=T, percentage=50, frame=F, absolute=F, ignore_start=T){
  gene_names<-genedf$gene

  if(percentage>1){
    percentage<-percentage/100
  }

  if(is.null(data$rna) & norms=="rna"){
    warning("No RNA data setting norms to total")
    norms="total"
  }
  #add ignore start option
  get_rrp<-function(gene_name, percentage){
    gene<-get_gene(gene_name, data=data, genedf=genedf)
    if (is.null(gene)) {
      NULL
    } else {
      coord <- genedf[genedf$gene == gene_name, ]
      if(ignore_start){
        gene<-gene[gene$nucleotide > 3 & gene$nucleotide < coord$end-coord$start,]
      } else {
        gene<-gene[gene$nucleotide > 0 & gene$nucleotide < coord$end-coord$start,]
      }
      gene<-gene[order(gene$nucleotide),]
      if(dim(gene)[1]==0){
        NULL
      } else {
        if(frame){
          gene<-gene[gene$frame==frame,]
        }

        if(norms=="total"){
          sums<-cumsum(na.omit(gene$freq))
          tots<-sum(na.omit(gene$freq))
          sums<-sums/tots
          gene$cumulative<-NA
          gene$cumulative[!is.na(gene$freq)]<-sums

          #fix warning issue
          lower<-gene[max(which(gene$cumulative<percentage)), c("nucleotide", "codon", "cumulative")]
          upper<-gene[min(which(gene$cumulative>percentage)), c("nucleotide", "codon", "cumulative")]
          nuc<-weighted.mean(x=c(lower$nucleotide, upper$nucleotide),
                             w=c(percentage-lower$cumulative, upper$cumulative-percentage))
          cod<-weighted.mean(x=c(lower$codon, upper$codon),
                             w=c(1/(percentage-lower$cumulative), 1/(upper$cumulative-percentage)))
          rrp<-data.frame(nucleotide=nuc, codon=cod, gene=gene_name, lower=lower$cumulative, upper=upper$cumulative)
          if(absolute){
            invisible(rrp)
          } else {
            rrp$nucleotide<-rrp$nucleotide/coord$nnuc
            rrp$codon<-rrp$codon/coord$ncodons
            invisible(rrp)
          }
        } else {
          norms<-gene$freq/gene$coverage
          sums<-cumsum(na.omit(norms))
          tots<-sum(na.omit(gene$freq))
          sums<-sums/tots
          gene$cumulative<-NA
          gene$cumulative[!is.na(gene$freq)]<-sums
          #ix warning issue
          nuc<-weighted.mean(x=c(lower$nucleotide, upper$nucleotide),
                             w=c(percentage-lower$cumulative, upper$cumulative-percentage))
          cod<-weighted.mean(x=c(lower$codon, upper$codon),
                             w=c(percentage-lower$cumulative, upper$cumulative-percentage))
          rrp<-data.frame(nucleotide=nuc, codon=cod, gene=gene_name, lower=lower$cumulative, upper=upper$cumulative)
          if(absolute){
            rrp$gene<-gene_name
            invisible(rrp)
          } else {
            rrp$nucleotide<-rrp$nucleotide/coord$nnuc
            rrp$codon<-rrp$codon/coord$ncodons
            rrp$gene<-gene_name
            invisible(rrp)
          }
        }
      }
    }
  }

  if (cores>1){
    rel_ribo<-do.call("rbind", mclapply(gene_names, get_rrp, percentage, mc.cores = cores))
  } else {
    rel_ribo<-do.call("rbind", lapply(gene_names, get_rrp, percentage))
  }

  if(plot){
    if(absolute){
    print(ggplot(rel_ribo, aes(x=nucleotide))+geom_density())
    } else {
      print(ggplot(rel_ribo, aes(x=nucleotide))+geom_density()+xlab("fraction ORF"))

    }
  }

  invisible(rel_ribo)

}
