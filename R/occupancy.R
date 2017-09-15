occupancy <- function(data, seq, genedf, codon = T, norms = "rna", gene_list = F, frame = NA, cores=1, separate_start=T) {
    codon_summer <- function(gene, class, measure) {
        if (class == "codon") {
            if (measure == "ribo") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["counts"]])
                })
            } else if (measure == "rna") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["coverage"]])
                })
            }
        } else if (class == "aa") {
            if (measure == "ribo") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["counts"]])
                })
            } else if (measure == "rna") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["coverage"]])
                })
            }
        }
        sums
    }

    summer <- function(nuc, df) {
        a <- df[df$nucleotide == nuc, ]
        a <- sum(a[["freq"]], na.rm = T)
        a
    }

    gene_names <- genedf$gene

    get_occupancy <- function(gene_name, seq = seq, data = data) {
        coord <- genedf[genedf$gene == gene_name, ]
        gene <- get_gene(gene = gene_name, data = data, seq = seq, genedf = genedf)
        if (!is.na(frame)) {
            gene <- gene[gene$frame == frame, ]
            gene <- gene[!is.na(gene$nucleotide), ]
        }
        if (is.null(gene)) {
            NULL
        } else {
            counts <- sapply(unique(gene$nucleotide), summer, gene)
            gene <- as.data.frame(cbind(unique(gene[, 1:2]), counts))
            gene <- gene[gene$nucleotide >= 1 & gene$nucleotide <= (coord$end - coord$start + 1), ]
            if (dim(gene)[1] < 3) {
                NULL
            } else {
                gene$coverage <- rollmean(gene$coverage, 3, fill = NA, align = "left")
                gene <- gene[(1:length(gene$nucleotide/3) * 3 - 2), ]
                cods <- as.data.frame(codons(seq[[gene_name]][coord$start:coord$end]))[, 1]
                aas <- as.data.frame(translate(seq[[gene_name]][coord$start:coord$end]))[, 1]
            }

            if (dim(gene)[1] < length(cods)) {
                gene$codon <- cods[floor(gene$nucleotide/3) + 1]
                gene$aa <- aas[floor(gene$nucleotide/3) + 1]
                gene$aa <- as.character(gene$aa)
            } else if (dim(gene)[1] >= length(cods)) {
                gene <- gene[1:length(cods), ]
                gene$codon <- cods
                gene$aa <- aas
                gene$aa <- as.character(gene$aa)
            }

            if(separate_start){
            gene$codon[1]<-"start"
            }

            if (codon) {
                occup <- data.frame(codon = unique(gene$codon), ribo = codon_summer(gene, "codon", "ribo"), rna = codon_summer(gene, "codon",
                  "rna"))
                n_obs <- as.data.frame(plyr::count(gene$codon))
                colnames(n_obs) <- c("codon", "n_obs")
                occup <- plyr::join(x = occup, y = n_obs, by = "codon", type = "full")
                occup
            } else {
                occup <- data.frame(aa = unique(gene$aa), ribo = codon_summer(gene, "aa", "ribo"), rna = codon_summer(gene, "aa", "rna"))
                n_obs <- as.data.frame(plyr::count(gene$aa))
                colnames(n_obs) <- c("aa", "n_obs")
                occup <- plyr::join(x = occup, y = n_obs, by = "aa", type = "full")
                occup
            }
            occup
        }
    }

    #added supress wornings not all genes (unspliced, pseudogene, putative genes all have ORFs that are divisible by 3)
    if (gene_list) {
      if(cores>1){
        suppressWarnings(occups <- mclapply(gene_names, get_occupancy, seq, data, mc.cores = cores))
        invisible(occups)
      } else {
        suppressWarnings(occups <- lapply(gene_names, get_occupancy, seq, data))
        invisible(occups)
      }
    } else {
      if(cores>1){
        suppressWarnings(occups <- do.call("rbind", mclapply(gene_names, get_occupancy, seq, data, mc.cores = cores)))
      } else {
        suppressWarnings(occups <- do.call("rbind", lapply(gene_names, get_occupancy, seq, data)))
      }

        get_sums <- function(val, x, col) {
            a <- x[which(x[, 1] == val), col]
            a <- sum(a, na.rm = T)
        }

        if (is.null(data[["rna"]]) & norms == "rna") {
            warning("No RNA-Seq Data Setting norms to total")
            norms = "total"
        }

       if (norms == "rna") {
            occups$norm <- occups$ribo/occups$rna
            torem <- grep("Inf", occups$norm)
            torem <- c(grep("NaN", occups$norm), torem)
            torem <- unique(torem)
            if (length(torem > 0)) {
                occups <- occups[-torem, ]
            }
            unq <- unique(occups[, 1])
            value <- lapply(unq, get_sums, occups, 5)
            n_obs <- lapply(unq, get_sums, occups, 4)
            b <- data.frame(feature = unq, value = unlist(value), n_obs = unlist(n_obs))
            invisible(b)
        } else if (norms == "total") {
            unq <- unique(occups[, 1])
            value <- lapply(unq, get_sums, occups, 2)
            n_obs <- lapply(unq, get_sums, occups, 4)
            b <- data.frame(feature = unq, value = unlist(value), n_obs = unlist(n_obs))
            b$value <- b$value/sum(b$value)
            invisible(b)
        } else if (norms == "none") {
            unq <- unique(occups[, 1])
            value <- lapply(unq, get_sums, occups, 2)
            n_obs <- lapply(unq, get_sums, occups, 4)
            b <- data.frame(feature = unq, value = unlist(value), n_obs = unlist(n_obs))
            b <- data.frame(feature = unq, value = unlist(a))
            invisible(b)
        }
    }
}
