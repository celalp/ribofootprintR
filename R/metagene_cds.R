metagene_cds <- function(data, genedf, before_start = 50, after_start = 100, before_stop = 100,
                         after_stop = 50, norms = "total", cores=1, plot=F) {

    returner <- function(nucleotide, frame, df) {
        a <- df[df$nucleotide == nucleotide & df$frame == frame, ]
        b <- sum(a$norm, na.rm = T)
        b
    }

    if (is.null(data[["rna"]]) & norms == "rna") {
        warning("No RNA-Seq data setting norms to total")
        norms = "total"
    }

    gene_names <- genedf$gene

    get_cds <- function(gene_name, direction, norm = norms, genes=genedf) {
        gene <- get_gene(gene_name, data = data, seq = NA, genedf = genes)
        if (is.null(gene) == T) {
            NULL
        } else {
            coord <- genes[genes$gene == gene_name, ]
            if (direction == "five") {
                cdss <- data.frame(nucleotide = c((before_start * -1):after_start))
                cdss <- left_join(cdss, gene, by = "nucleotide")
            } else if (direction == "three") {
                cdss <- data.frame(nucleotide = c((coord$end - coord$start - before_stop):(coord$end - coord$start + after_stop)))
                cdss <- left_join(cdss, gene, by = "nucleotide")
                cdss$nucleotide <- cdss$nucleotide - coord$end + coord$start
            }
            if (norm == "total") {
                cdss$norm <- as.numeric(cdss$freq/sum(gene$freq, na.rm = T))
            } else if (norm == "rna") {
                cdss$norm <- as.numeric(cdss$freq/cdss$coverage)
            } else if (norm == "none") {
                cdss$norm <- as.numeric(cdss$freq)
            }
            cdss$norm[!is.finite(cdss$norm)] <- NA
            cdss
        }
    }

    if (cores>1){
    fives<-do.call("rbind", mclapply(gene_names, get_cds, "five", mc.cores = cores))
    } else {
      fives<-do.call("rbind", lapply(gene_names, get_cds, "five"))
    }
    x<-unique(fives[,c(1,3)])
    x<-na.omit(x)
    counts<-mdply(x, returner, fives)
    fives<-counts
    #fives<-dplyr::left_join(unique(fives[,c(1,2)]), counts, by="nucleotide")
    colnames(fives)[3]<-"value"
    fives$codon<-ceiling(fives$nucleotide/3)

    if(cores>1){
    threes<-do.call("rbind", mclapply(gene_names, get_cds, "three", mc.cores = cores))
    } else {
      threes<-do.call("rbind", lapply(gene_names, get_cds, "three"))
    }
    x<-unique(threes[,c(1,3)])
    x<-na.omit(x)
    counts<-mdply(x, returner, threes)
    threes<-counts
    #threes<-dplyr::left_join(unique(threes[,c(1,2)]), counts, by="nucleotide")
    colnames(threes)[3]<-"value"
    threes$codon<-ceiling(threes$nucleotide)

    fives$location<-"from_start"
    threes$location<-"from_stop"
    cds<-rbind(fives, threes)
    cds


    fives$location <- "from_start"
    threes$location <- "from_stop"
    cds <- rbind(fives, threes)

    if(plot){
      print(ggplot(cds, aes(x=codon, y=value, fill=frame))+geom_bar(stat="identity", position="dodge")+facet_grid(.~location, scales = "free_x"))
    }

    invisible(cds)
}

