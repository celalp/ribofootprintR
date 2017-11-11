\name{rrp}
\alias{rrp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate relative ribosome positioning
}
\description{
Calculate where within the orf the cumulative location of ribosomal reads are.
}
\usage{
rrp(data, genedf, norms = "total", cores = 1, plot = T, percentage = 50, frame = F, absolute = F, ignore_start = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  ribodata, the output of read_profiling_data
}
  \item{genedf}{
  data.frame, stating the relative locations of start and stop locations in the transcriptome. Can be manually generated or through generate_transcriptome function
%%     ~~Describe \code{genedf} here~~
}
  \item{norms}{
  string, "total" for fraction of reads per gene or "rna" for first getting ribo/rna-seq ratio. Default "total"
%%     ~~Describe \code{norms} here~~
}
  \item{cores}{
  integer, number of threads to run. only on UNIX systems see details, defaults to 1
%%     ~~Describe \code{cores} here~~
}
  \item{plot}{
  logical, should the function produce a density plot defaults to T
%%     ~~Describe \code{plot} here~~
}
  \item{percentage}{
  float or integer, the percentage of the reads in the cumulative function. Default 50. The function converts the percentage to fraction if the value is greater than 1
%%     ~~Describe \code{percentage} here~~
}
  \item{frame}{
  F or integer 0,1,2 should the function use reads coming from a specific frame, default NA uses all reads.
%%     ~~Describe \code{frame} here~~
}
  \item{absolute}{
  Logical, should the relative or absoulte positions be returned. Defaults to F
%%     ~~Describe \code{absolute} here~~
}
  \item{ignore_start}{
  logical, should the calculations start from the first codon or the second one. This might be useful for libraries that are prepared with cycloheximide
%%     ~~Describe \code{ignore_start} here~~
}
}
\details{
Since the data is discerete the way this function goes around calculating the percent specified by the user is as follows:
The cumulative sum is calculated for each non NA data point and it's corresponding nucleotide position is found.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Invisible returns a data frame. One row per gene with nucleotide and codon location (relative or absolute depending on settings), gene name, and the percentile locations around the percent input. The weighted average (1/ distance from percent) of he largest read location that's smaller than the precent specified and smallest read location larger than the percent specified is calculated. I would be interested in suggestion for estimating an exact cdf per gene.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
This was first concieved by Richard Baker Ph.D. at UMASS Medical School. This is my version of the calculations that were done in his original perl script.
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
function (data, genedf, norms = "total", cores = 1, plot = T,
    percentage = 50, frame = F, absolute = F, ignore_start = T)
{
    gene_names <- genedf$gene
    if (percentage > 1) {
        percentage <- percentage/100
    }
    if (is.null(data$rna) & norms == "rna") {
        warning("No RNA data setting norms to total")
        norms = "total"
    }
    get_rrp <- function(gene_name, percentage) {
        gene <- get_gene(gene_name, data = data, genedf = genedf)
        if (is.null(gene)) {
            NULL
        }
        else {
            coord <- genedf[genedf$gene == gene_name, ]
            if (ignore_start) {
                gene <- gene[gene$nucleotide > 3 & gene$nucleotide <
                  coord$end - coord$start, ]
            }
            else {
                gene <- gene[gene$nucleotide > 0 & gene$nucleotide <
                  coord$end - coord$start, ]
            }
            gene <- gene[order(gene$nucleotide), ]
            if (dim(gene)[1] == 0) {
                NULL
            }
            else {
                if (frame) {
                  gene <- gene[gene$frame == frame, ]
                }
                if (norms == "total") {
                  sums <- cumsum(na.omit(gene$freq))
                  tots <- sum(na.omit(gene$freq))
                  sums <- sums/tots
                  gene$cumulative <- NA
                  gene$cumulative[!is.na(gene$freq)] <- sums
                  lower <- gene[max(which(gene$cumulative < percentage)),
                    c("nucleotide", "codon", "cumulative")]
                  upper <- gene[min(which(gene$cumulative > percentage)),
                    c("nucleotide", "codon", "cumulative")]
                  nuc <- weighted.mean(x = c(lower$nucleotide,
                    upper$nucleotide), w = c(percentage - lower$cumulative,
                    upper$cumulative - percentage))
                  cod <- weighted.mean(x = c(lower$codon, upper$codon),
                    w = c(1/(percentage - lower$cumulative),
                      1/(upper$cumulative - percentage)))
                  rrp <- data.frame(nucleotide = nuc, codon = cod,
                    gene = gene_name, lower = lower$cumulative,
                    upper = upper$cumulative)
                  if (absolute) {
                    invisible(rrp)
                  }
                  else {
                    rrp$nucleotide <- rrp$nucleotide/coord$nnuc
                    rrp$codon <- rrp$codon/coord$ncodons
                    invisible(rrp)
                  }
                }
                else {
                  norms <- gene$freq/gene$coverage
                  sums <- cumsum(na.omit(norms))
                  tots <- sum(na.omit(gene$freq))
                  sums <- sums/tots
                  gene$cumulative <- NA
                  gene$cumulative[!is.na(gene$freq)] <- sums
                  nuc <- weighted.mean(x = c(lower$nucleotide,
                    upper$nucleotide), w = c(percentage - lower$cumulative,
                    upper$cumulative - percentage))
                  cod <- weighted.mean(x = c(lower$codon, upper$codon),
                    w = c(percentage - lower$cumulative, upper$cumulative -
                      percentage))
                  rrp <- data.frame(nucleotide = nuc, codon = cod,
                    gene = gene_name, lower = lower$cumulative,
                    upper = upper$cumulative)
                  if (absolute) {
                    rrp$gene <- gene_name
                    invisible(rrp)
                  }
                  else {
                    rrp$nucleotide <- rrp$nucleotide/coord$nnuc
                    rrp$codon <- rrp$codon/coord$ncodons
                    rrp$gene <- gene_name
                    invisible(rrp)
                  }
                }
            }
        }
    }
    if (cores > 1) {
        rel_ribo <- do.call("rbind", mclapply(gene_names, get_rrp,
            percentage, mc.cores = cores))
    }
    else {
        rel_ribo <- do.call("rbind", lapply(gene_names, get_rrp,
            percentage))
    }
    if (plot) {
        if (absolute) {
            print(ggplot(rel_ribo, aes(x = nucleotide)) + geom_density())
        }
        else {
            print(ggplot(rel_ribo, aes(x = nucleotide)) + geom_density() +
                xlab("fraction ORF"))
        }
    }
    invisible(rel_ribo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
