\name{totals}
\alias{totals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
retuns the total number of reads that map to an ORF and the total coverage per ORF
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
totals(data, genedf, seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  object returned by read_profiling_data
%%     ~~Describe \code{data} here~~
}
  \item{genedf}{
  data frame with the start and stop coordinated of the ORFs in the transcriptome.

%%     ~~Describe \code{genedf} here~~
}
  \item{seq}{
  fix this
%%     ~~Describe \code{seq} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a data frame where each row is a gene and sums of RNA-Seq coverage and total number of ribosome profiling reads per gene. If there is no RNA-Seq data coverages are set to \code{NA}
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
function (data, genedf, seq)
{
    summer <- function(nuc, df) {
        a <- df[df$nucleotide == nuc, ]
        a <- sum(a[["freq"]], na.rm = T)
        a
    }
    genes <- genedf$gene
    get_data <- function(gene_name) {
        gene <- get_gene(gene = gene_name, data = data, seq = seq)
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
    totals <- do.call("rbind", mclapply(genes, get_data, mc.cores = (detectCores() -
        1)))
    if (is.null(data$rna)) {
        totals$sum_coverage <- NA
    }
    totals
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
