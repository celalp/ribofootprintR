get_gene <- function(gene_name, data, plot = F, seq = NA, genedf) {
    ribo <- data[["ribo"]][[gene_name]]
    rna <- data[["rna"]][[gene_name]]

    if (is.null(ribo) == T) {
        gene <- NULL
    } else {
        coord <- genedf[genedf$gene == gene_name, ]
        coord <- coord[coord[, 1] == gene_name, ]
        if (is.null(rna) == T) {
            gene <- ribo
            gene$frame <- as.factor(gene$frame)
            if (plot) {
                plot <- ggplot() + geom_bar(data = gene, aes(x = nucleotide, y = freq, fill = frame), stat = "identity",
                  position = "dodge")
                plot <- plot + geom_segment(aes(x = 0, y = 0, xend = (coord$end - coord$start) + 1, yend = 0), color = "blueviolet",
                  size = 2)
                print(plot)

            }
        } else {
            gene <- dplyr::full_join(rna, ribo, by = "nucleotide")
            gene$frame <- as.factor(gene$frame)

            if (plot) {
                plot <- ggplot() + geom_line(data = gene, aes(x = nucleotide, y = -coverage)) + geom_bar(data = gene[!is.na(gene$frame),
                  ], aes(x = nucleotide, y = freq, fill = frame), stat = "identity", position = "dodge")
                plot <- plot + geom_segment(aes(x = 0, y = 0, xend = coord$end - coord$start, yend = 0), color = "blueviolet", size = 2)
                print(plot)
            }
        }
    }
    invisible(gene)
}
