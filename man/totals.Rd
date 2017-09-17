\name{totals}
\alias{totals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Returns the total number of reads that map to an ORF and the total coverage per ORF
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
totals(data, genedf, plot = F, cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  object returned by \code{read_profiling_data}

%%     ~~Describe \code{data} here~~
}
  \item{genedf}{
  data frame with the start and stop coordinated of the ORFs in the transcriptome.

%%     ~~Describe \code{genedf} here~~
}
  \item{plot}{
  logical, whether to plot a scatter plot of results, defaults to F
%%     ~~Describe \code{plot} here~~
}
  \item{cores}{
  integer, number of threads to use. See details.
%%     ~~Describe \code{cores} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a data frame where each row is a gene and sums of RNA-Seq coverage and total number of ribosome profiling reads per gene. If there is no RNA-Seq data coverages are set to NA

Currently PSOCK multithreadding is not supported. For Windows set cores to 1. For operating systems that support forking the number of threads that are used can be set by the cores argument. If the cores>1 then the data is processed by mclapply instead of lapply.

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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, genedf, plot = F, cores = 1)
{
    summer <- function(nuc, df) {
        a <- df[df$nucleotide == nuc, ]
        a <- sum(a[["freq"]], na.rm = T)
        a
    }
    genes <- genedf$gene
    get_data <- function(gene_name) {
        gene <- get_gene(gene = gene_name, genedf = genedf, data = data,
            seq = NA)
        if (is.null(gene)) {
            NULL
        }
        else {
            coord <- genedf[genedf$gene == gene_name, ]
            gene <- gene[gene$nucleotide >= 0 & gene$nucleotide <=
                (coord$end - coord$start + 1), ]
            counts <- sapply(unique(gene$nucleotide), summer,
                gene)
            gene <- as.data.frame(cbind(unique(gene[, 1:2]),
                counts))
            if (dim(gene)[1] == 0) {
                NULL
            }
            else {
                dats <- colSums(gene[, 2:3])
                gene <- data.frame(sum_coverage = dats[1], sum_ribo = dats[2])
                rownames(gene) <- gene_name
                gene
            }
        }
    }
    if (cores > 1) {
        totals <- do.call("rbind", mclapply(genes, get_data,
            mc.cores = cores))
        if (is.null(data$rna)) {
            totals$sum_coverage <- NA
        }
    }
    else {
        totals <- do.call("rbind", lapply(genes, get_data))
        if (is.null(data$rna)) {
            totals$sum_coverage <- NA
        }
    }
    totals <- as.data.frame(totals)
    if (plot) {
        print(ggplot(totals, aes(x = sum_coverage, y = sum_ribo)) +
            geom_point())
    }
    invisible(totals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
