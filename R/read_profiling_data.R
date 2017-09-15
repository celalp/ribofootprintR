read_profiling_data <- function(ribo, ribo_index, rna = NA, rna_index = NA, genedf, cores=1, A_site_adjust=T) {
    message("Loading ribo data")

    read_data <- function(data, index) {
        message("read bam file")
        gr <- readGAlignments(file = data, index = index)
        invisible(gc())
        message("Split into individual genes")
        gr <- split(gr, seqnames(gr))
        invisible(gc)
        return(gr)
    }

    get_coverage <- function(gene) {
        covs <- cov_list[[gene]]
        if (is.null(covs)) {
            NULL
        } else {
            covs <- as.data.frame(covs)
            b <- genedf[genedf$gene == gene, ]
            coverage <- data.frame(nucleotide = as.integer(1:nrow(covs)), coverage = covs$value)
            coverage$nucleotide <- coverage$nucleotide - b$start
            coverage
        }
    }

    get_frames <- function(gr) {
        if (is.null(gr)) {
            gr <- data.frame()
        } else {
            b <- genedf[genedf$gene == as.character(gr@seqnames@values), ]
            nucleotide<-gr@start-b$start
            if(A_site_adjust){
              frames<-data.frame(frame=nucleotide-b$start+1+(width(gr)-15))%%3
              frames$codon<-floor((start(gr)-b$start+(round(width(gr)/2)))/3)+1
              frames$width<-width(gr)
              } else {
              frames <- data.frame(frame = ((nucleotide)%%3))
              frames$codon <- floor(nucleotide/3)
            }
            frames$width <- width(gr)
            frames$nucleotide <- frames$codon * 3 - 2
            frames <- plyr::count(frames)
            frames <- frames[, c(4, 1, 2, 3, 5)]
            frames
        }
    }
    ribo <- read_data(data = ribo, index = ribo_index)

    message("Calling Frames")
    if(cores>1){
      ribo <- mclapply(ribo, get_frames, mc.cleanup = T, mc.cores = cores)
    } else {
      ribo <- lapply(ribo, get_frames)
    }
    invisible(gc())
    message("Done")


    if (is.na(rna)) {
        warning("No rna data")
        coverages <- NULL
    } else {
        message("Generating coverage list")
        cov_list <- coverage(rna)

        message("Done")
        message("Split into individual genes")
        genes <- genedf$gene
        coverages <- lapply(X = genes, FUN = get_coverage)
        names(coverages) <- genedf$gene
        invisible(gc())
        message("Done")
    }
    invisible(gc)
    data <- list(ribo = ribo, rna = coverages)
    invisible(data)
  }
