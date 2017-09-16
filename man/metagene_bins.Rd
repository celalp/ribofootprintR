\name{metagene_bins}
\alias{metagene_bins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Divide ORF into bins and retun either gene by gene or a summary value per bin

%%  ~~function to do ... ~~
}
\description{
Divide the ORF in the bins and return the value either gene by gene or as cumulative mean or median.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metagene_bins(data, genedf, bins = 100, norms = "rna", simplify = "mean", cores = 1, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  list returned by read_profiling_data() function

%%     ~~Describe \code{data} here~~
}
  \item{genedf}{
   data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF. This can be the whole trancriptome or a subset of the genedf used by \code{read_profiling_data}
%%     ~~Describe \code{genedf} here~~
}
  \item{bins}{
  number of bins to divide the ORF with defaults to 100

%%     ~~Describe \code{bins} here~~
}
  \item{norms}{
  character either "total", "rna" or "none"
%%     ~~Describe \code{norms} here~~
}
  \item{simplify}{
  character either "mean", "median", "none"

%%     ~~Describe \code{simplify} here~~
}
  \item{cores}{
  number of threads to use in UNIX based systems. Defaults to 1.
%%     ~~Describe \code{cores} here~~
}
  \item{plot}{
  logical, whether to plot the results. See details.
%%     ~~Describe \code{plot} here~~
}
}
\details{

%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a molten data frame. Each ORF is divided by the number of bins specified and the mean occpancy is calculated for each ORF. If norms is selected to be "rna" the Ribosome profiling data is divided by the coverage from RNA-Seq. If norms is "total" then the fraction of reads mapping to each bin is calculted. "none" just returns average values. For the simplify argument if "none" is selected a molten data frame is returned where each gene and each bin is a row. This is done to make plotting easier. For other simplify options only one value per bin is returned.

Currently PSOCK multithreadding is not supported. For Windows set cores to 1. For operating systems that support forking the number of threads that are used can be set by the cores argument. If the cores>1 then the data is processed by mclapply instead of lapply.

In the simplify argumnet if "mean" or "median" is selected a line plot is returned. If simplify is set to be "none" then plotting is done by geom_smooth which uses generalized additive models.
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{metagene_cds}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, genedf, bins = 100, norms = "rna", simplify = "mean",
    cores = 1, plot = F)
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
    gene_list <- genedf
    norm_gene <- function(name, gene_list) {
        gene <- get_gene(gene = name, data = data, seq = NA,
            genedf = gene_list)
        if (is.null(gene) || dim(gene)[1] < 1) {
            NULL
        }
        else {
            coord <- genedf[genedf$gene == name, ]
            gene <- gene[gene$nucleotide >= 0 & gene$nucleotide <=
                (coord$end - coord$start + 1), ]
            if (dim(gene)[1] < 1) {
                NULL
            }
            else {
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
    }
    binner <- function(gene_name) {
        gene <- norm_gene(name = gene_name, gene_list = gene_list)
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
    genes <- genedf$gene
    if (cores > 1) {
        suppressWarnings(binned <- do.call("rbind", mclapply(genes,
            binner, mc.cores = cores, mc.cleanup = T)))
    }
    else {
        suppressWarnings(binned <- do.call("rbind", lapply(genes,
            binner)))
    }
    rownames(binned) <- genes
    if (simplify == "mean") {
        binned <- data.frame(bin = 1:bins, value = apply(binned,
            2, mean, na.rm = T))
        print(ggplot(binned, aes(x = bin, y = value)) + geom_line())
    }
    else if (simplify == "median") {
        binned <- data.frame(bin = 1:bins, value = apply(binned,
            2, median, na.rm = T))
        print(ggplot(binned, aes(x = bin, y = value)) + geom_line())
    }
    else if (simplify == "none") {
        binned <- reshape2::melt(binned)
        colnames(binned) <- c("gene", "bin", "value")
        if (plot) {
            print(ggplot(binned, aes(x = bin, y = value)) + geom_smooth())
        }
    }
    invisible(binned)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
