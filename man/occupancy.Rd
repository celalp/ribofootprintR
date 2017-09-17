\name{occupancy}
\alias{occupancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get ribosome occupancy at codon or amino acid level

%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
occupancy(data, seq, genedf, codon = T, norms = "rna", gene_list = F, frame = NA, cores = 1, separate_start = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  list returned by read_profiling_data

%%     ~~Describe \code{data} here~~
}
  \item{seq}{
  A DNAStringSet of the FASTA file used for alignment.

%%     ~~Describe \code{seq} here~~
}
  \item{genedf}{
  data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF.

%%     ~~Describe \code{genedf} here~~
}
  \item{codon}{
  logical, whether to return codon level or amino acid level occupancy, defaults to TRUE

%%     ~~Describe \code{codon} here~~
}
  \item{norms}{
  normalization, character either "total", "rna" or "none". Defaults to "total"

%%     ~~Describe \code{norms} here~~
}
  \item{gene_list}{
  logical, whether to aggregate all genes in genedf or return gene by gene values

%%     ~~Describe \code{gene_list} here~~
}
  \item{frame}{
  integer, 0,1,2 or NA indicating whether only reads from a specific frame should be counted. If set to NA all reads are counted.
%%     ~~Describe \code{frame} here~~
}
  \item{cores}{
  integer, number of threads to use. See details
%%     ~~Describe \code{cores} here~~
}
  \item{separate_start}{
  whether to tread start codon separately from other methionines, defaults to T, codon needs to be set to T
%%     ~~Describe \code{separate_start} here~~
}
}
\details{
Normalizations are similar to other functions. If "none" is selected then the total number of reads are returned. If "rna" is selected then the the average of the three nucleotides of the codon is used for normalization. If "total" is used then the fraction of reads that correspond to each codon or amino acid is returned.

Currently PSOCK multithreadding is not supported. For Windows set cores to 1. For operating systems that support forking the number of threads that are used can be set by the cores argument. If the cores>1 then the data is processed by mclapply instead of lapply.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
a data frame containing the number of reads (or normalized values) with the number of codons or amino acids in the selection is returned. if gene_list=T then a list is returned where each element is a data frame with occupancy values and number of codons or amino acids.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alper Celik
%%  ~~who you are~~
}
\note{
warnings are suppressed because some pseudogenes or putative proteins do not have annotated ORF lengths are not divisible by 3.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{totals}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, seq, genedf, codon = T, norms = "rna", gene_list = F,
    frame = NA, cores = 1, separate_start = T)
{
    codon_summer <- function(gene, class, measure) {
        if (class == "codon") {
            if (measure == "ribo") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["counts"]])
                })
            }
            else if (measure == "rna") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["coverage"]])
                })
            }
        }
        else if (class == "aa") {
            if (measure == "ribo") {
                sums <- sapply(unique(gene[[class]]), function(x) {
                  sum(gene[which(gene[[class]] == x), ][["counts"]])
                })
            }
            else if (measure == "rna") {
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
        gene <- get_gene(gene = gene_name, data = data, seq = seq,
            genedf = genedf)
        if (!is.na(frame)) {
            gene <- gene[gene$frame == frame, ]
            gene <- gene[!is.na(gene$nucleotide), ]
        }
        if (is.null(gene)) {
            NULL
        }
        else {
            counts <- sapply(unique(gene$nucleotide), summer,
                gene)
            gene <- as.data.frame(cbind(unique(gene[, 1:2]),
                counts))
            gene <- gene[gene$nucleotide >= 1 & gene$nucleotide <=
                (coord$end - coord$start + 1), ]
            if (dim(gene)[1] < 3) {
                NULL
            }
            else {
                gene$coverage <- rollmean(gene$coverage, 3, fill = NA,
                  align = "left")
                gene <- gene[(1:length(gene$nucleotide/3) * 3 -
                  2), ]
                cods <- as.data.frame(codons(seq[[gene_name]][coord$start:coord$end]))[,
                  1]
                aas <- as.data.frame(translate(seq[[gene_name]][coord$start:coord$end]))[,
                  1]
            }
            if (dim(gene)[1] < length(cods)) {
                gene$codon <- cods[floor(gene$nucleotide/3) +
                  1]
                gene$aa <- aas[floor(gene$nucleotide/3) + 1]
                gene$aa <- as.character(gene$aa)
            }
            else if (dim(gene)[1] >= length(cods)) {
                gene <- gene[1:length(cods), ]
                gene$codon <- cods
                gene$aa <- aas
                gene$aa <- as.character(gene$aa)
            }
            if (separate_start) {
                gene$codon[1] <- "start"
            }
            if (codon) {
                occup <- data.frame(codon = unique(gene$codon),
                  ribo = codon_summer(gene, "codon", "ribo"),
                  rna = codon_summer(gene, "codon", "rna"))
                n_obs <- as.data.frame(plyr::count(gene$codon))
                colnames(n_obs) <- c("codon", "n_obs")
                occup <- plyr::join(x = occup, y = n_obs, by = "codon",
                  type = "full")
                occup
            }
            else {
                occup <- data.frame(aa = unique(gene$aa), ribo = codon_summer(gene,
                  "aa", "ribo"), rna = codon_summer(gene, "aa",
                  "rna"))
                n_obs <- as.data.frame(plyr::count(gene$aa))
                colnames(n_obs) <- c("aa", "n_obs")
                occup <- plyr::join(x = occup, y = n_obs, by = "aa",
                  type = "full")
                occup
            }
            occup
        }
    }
    if (gene_list) {
        if (cores > 1) {
            suppressWarnings(occups <- mclapply(gene_names, get_occupancy,
                seq, data, mc.cores = cores))
            invisible(occups)
        }
        else {
            suppressWarnings(occups <- lapply(gene_names, get_occupancy,
                seq, data))
            invisible(occups)
        }
    }
    else {
        if (cores > 1) {
            suppressWarnings(occups <- do.call("rbind", mclapply(gene_names,
                get_occupancy, seq, data, mc.cores = cores)))
        }
        else {
            suppressWarnings(occups <- do.call("rbind", lapply(gene_names,
                get_occupancy, seq, data)))
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
            b <- data.frame(feature = unq, value = unlist(value),
                n_obs = unlist(n_obs))
            invisible(b)
        }
        else if (norms == "total") {
            unq <- unique(occups[, 1])
            value <- lapply(unq, get_sums, occups, 2)
            n_obs <- lapply(unq, get_sums, occups, 4)
            b <- data.frame(feature = unq, value = unlist(value),
                n_obs = unlist(n_obs))
            b$value <- b$value/sum(b$value)
            invisible(b)
        }
        else if (norms == "none") {
            unq <- unique(occups[, 1])
            value <- lapply(unq, get_sums, occups, 2)
            n_obs <- lapply(unq, get_sums, occups, 4)
            b <- data.frame(feature = unq, value = unlist(value),
                n_obs = unlist(n_obs))
            b <- data.frame(feature = unq, value = unlist(a))
            invisible(b)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
