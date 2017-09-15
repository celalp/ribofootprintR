select_read_lengths <- function(data, ribo_length = c(27:29)) {
    ribo <- data[["ribo"]]
    rna <- data[["rna"]]
    features <- data[["features"]]
    selector <- function(gene, lengths) {
        gene <- ribo[[gene]]
        gene <- gene[gene$width >= min(lengths) & gene$width <= max(lengths), ]
        gene
    }
    ribo_selected <- lapply(names(ribo), selector, lengths = ribo_length)
    
    names(ribo_selected) <- names(ribo)
    data$ribo <- ribo_selected
    data
}
