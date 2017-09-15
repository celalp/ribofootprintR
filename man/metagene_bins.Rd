\name{metagene_bins}
\alias{metagene_bins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
divide ORF into bins and retun either gene by gene or a summary value per bin
}
\description{
Divide the ORF in the bins and return the value either gene by gene or as cumulative mean or median.
}
\usage{
metagene_bins(data, genedf, bins = 100, norms = "rna", simplify = "mean", seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  list returned by read_profiling_data() function
}
  \item{genedf}{
a data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF.
}
  \item{bins}{
number of bins to divide the ORF with defaults to 100
}
  \item{norms}{
character either "total", "rna" or "none"
}
  \item{simplify}{
character either "mean", "median", "none"
}
  \item{seq}{
DNAString set object of the FASTA file used for alignment. This is only needed for the occupancy function otherwise can be set to NA
}
}
\details{
Returns a matrix. Each ORF is divided by the number of bins specified and the mean occpancy is calculated for each ORF. If norms is selected to be "rna" the Ribosome profiling data is divided by the coverage from RNA-Seq. If norms is "total" then the fraction of reads mapping to each bin is calculted. "none" just returns average values. For the simplify argument if none is selected a molten data frame is returned where each gene and each bin is a row. This is done to make plotting using easier. For other simplify options only one value per bin is returned.
}
\value{
A data frame with values specified in details
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alper Celik
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read_profiling_data}}, \code{\link{metagene_cds}}
}
\examples{

TODO


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, genedf, bins = 100, norms = "rna", simplify = "mean",
    seq)
{
    summer <- function(nuc, df, measure) {
        a <- df[df$nucleotide == nuc, ]
        a <- sum(a[[measure]], na.rm = T)
        a
    }
    if (norms == "rna" & is.null(data[["rna"]])) {
        warning("No RNA-Seq Data Setting norms to total")
        norms = "total"
    }
    norm_gene <- function(gene_name) {
        gene <- get_gene(gene = gene_name, data = data, seq = seq)
        if (is.null(gene)) {
            NULL
        }
        else {
            coord <- genedf[genedf$gene == gene_name, ]
            gene <- gene[gene$nucleotide >= 0 & gene$nucleotide <=
                (coord$end - coord$start + 1), ]
            counts <- sapply(unique(gene$nucleotide), summer,
                gene, "freq")
            coverage <- sapply(unique(gene$nucleotide), summer,
                gene, "coverage")
            if (norms == "none") {
                a <- counts
                a
            }
            else if (norms == "rna") {
                a <- counts/coverage
                a
            }
            else if (norms == "total") {
                a <- counts/sum(counts, na.rm = T)
                a
            }
        }
    }
    binner <- function(gene_name) {
        gene <- norm_gene(gene_name)
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
        }
        else {
            gene_bins
        }
    }
    noinf <- function(a) {
        loc <- is.infinite(a)
        loc2 <- is.nan(a)
        a <- a[!loc]
        a <- a[!loc2]
    }
    genes <- genedf$gene
    binned <- do.call("rbind", mclapply(genes, binner, mc.cores = (detectCores() -
        1)))
    rownames(binned) <- genes
    binned <- apply(binned, 2, noinf)
    if (simplify == "mean") {
        binned <- colMeans(binned, na.rm = T)
        binned
    }
    else if (simplify == "median") {
        binned <- apply(binned, 2, median, na.rm = T)
        binned
    }
    else if (simplify == "none") {
        binned <- reshape2::melt(binned)
        colnames(binned) <- c("gene", "bin", "value")
        binned
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
