phasing <- function(data, genedf, cores=1, plot=F) {
    genes <- genedf$gene
    ribo <- data[["ribo"]]
    get_phasing <- function(gene) {
        orf <- ribo[[gene]]
        if (is.null(orf)) {
            orf <- NULL
            orf
        } else {
            features <- genedf[genedf$gene == gene, ]
            orf <- orf[orf$codon > 0 & orf$codon < features$ncodon, c(1, 2, 4)]
            orf
        }
    }
    if(cores>1){
      orfs <- do.call("rbind", mclapply(genes, get_phasing, mc.cores = cores))
    } else {
      orfs <- do.call("rbind", lapply(genes, get_phasing))
    }
    phasing <- reshape2::dcast(data = orfs, width ~ frame, fun.aggregate = sum)
    colnames(phasing) <- c("read_length", "frame_0", "frame_1", "frame_2")
    phasing$read_length <- as.factor(phasing$read_length)
    phasing <- reshape2::melt(phasing)
    colnames(phasing) <- c("read_length", "frame", "number_of_reads")
    if(plot){
      print(ggplot(phasing, aes(x=read_length, y=number_of_reads, fill=frame))+geom_bar(stat="identity", position="dodge")+
              theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    }
    invisible(phasing)

}
