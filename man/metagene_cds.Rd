\name{metagene_cds}
\alias{metagene_cds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get average ribosome footprint alignment values at the 5' and 3' ends of coding sequences
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metagene_cds(data, genedf, before_start = 50, after_start = 100, before_stop = 100, after_stop = 50, norms = "total", cores = 1, plot = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  list retuned by read_profiling_data

%%     ~~Describe \code{data} here~~
}
  \item{genedf}{
  data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF.

%%     ~~Describe \code{genedf} here~~
}
  \item{before_start}{
  integer, the number of nucleotide to calculate before the start nucleotide, defaults to 50

%%     ~~Describe \code{before_start} here~~
}
  \item{after_start}{
  integer, the number of nucleotide to calculate after the start nucleotide, defaults to 100

%%     ~~Describe \code{after_start} here~~
}
  \item{before_stop}{
  integer, the number of nucleotide to calculate before the stop nucleotide, defaults to 100

%%     ~~Describe \code{before_stop} here~~
}
  \item{after_stop}{
  integer, the number of nucleotide to calculate before the stop nucleotide, defaults to 50

%%     ~~Describe \code{after_stop} here~~
}
  \item{norms}{
  character, normalization, character either "total", "rna" or "none". Defaults to "total"

%%     ~~Describe \code{norms} here~~
}
  \item{cores}{
  integer, number of threads to use. See details.
%%     ~~Describe \code{cores} here~~
}
  \item{plot}{
  logical, whether to plot or not, defaults to F
%%     ~~Describe \code{plot} here~~
}
}
\details{
If norms is set to "total"" returns the average fraction of reads that map to the region specified. If norms are set to "rna" then the number of reads that map to each codon are divided by the RNA-Seq coverage. "none" returns the average number of reads per codon.

Currently PSOCK multithreadding is not supported. For Windows set cores to 1. For operating systems that support forking the number of threads that are used can be set by the cores argument. If the cores>1 then the data is processed by mclapply instead of lapply.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a data frame with normalized nucleotide values.

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
\code{\link{metagene_bins}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, genedf, before_start = 50, after_start = 100,
    before_stop = 100, after_stop = 50, norms = "total", cores = 1,
    plot = F)
{
    returner <- function(nucleotide, frame, df) {
        a <- df[df$nucleotide == nucleotide & df$frame == frame,
            ]
        b <- sum(a$norm, na.rm = T)
        b
    }
    if (is.null(data[["rna"]]) & norms == "rna") {
        warning("No RNA-Seq data setting norms to total")
        norms = "total"
    }
    gene_names <- genedf$gene
    get_cds <- function(gene_name, direction, norm = norms, genes = genedf) {
        gene <- get_gene(gene_name, data = data, seq = NA, genedf = genes)
        if (is.null(gene) == T) {
            NULL
        }
        else {
            coord <- genes[genes$gene == gene_name, ]
            if (direction == "five") {
                cdss <- data.frame(nucleotide = c((before_start *
                  -1):after_start))
                cdss <- left_join(cdss, gene, by = "nucleotide")
            }
            else if (direction == "three") {
                cdss <- data.frame(nucleotide = c((coord$end -
                  coord$start - before_stop):(coord$end - coord$start +
                  after_stop)))
                cdss <- left_join(cdss, gene, by = "nucleotide")
                cdss$nucleotide <- cdss$nucleotide - coord$end +
                  coord$start
            }
            if (norm == "total") {
                cdss$norm <- as.numeric(cdss$freq/sum(gene$freq,
                  na.rm = T))
            }
            else if (norm == "rna") {
                cdss$norm <- as.numeric(cdss$freq/cdss$coverage)
            }
            else if (norm == "none") {
                cdss$norm <- as.numeric(cdss$freq)
            }
            cdss$norm[!is.finite(cdss$norm)] <- NA
            cdss
        }
    }
    if (cores > 1) {
        fives <- do.call("rbind", mclapply(gene_names, get_cds,
            "five", mc.cores = cores))
    }
    else {
        fives <- do.call("rbind", lapply(gene_names, get_cds,
            "five"))
    }
    x <- unique(fives[, c(1, 3)])
    x <- na.omit(x)
    counts <- mdply(x, returner, fives)
    fives <- counts
    colnames(fives)[3] <- "value"
    fives$codon <- ceiling(fives$nucleotide/3)
    if (cores > 1) {
        threes <- do.call("rbind", mclapply(gene_names, get_cds,
            "three", mc.cores = cores))
    }
    else {
        threes <- do.call("rbind", lapply(gene_names, get_cds,
            "three"))
    }
    x <- unique(threes[, c(1, 3)])
    x <- na.omit(x)
    counts <- mdply(x, returner, threes)
    threes <- counts
    colnames(threes)[3] <- "value"
    threes$codon <- ceiling(threes$nucleotide)
    fives$location <- "from_start"
    threes$location <- "from_stop"
    cds <- rbind(fives, threes)
    cds
    fives$location <- "from_start"
    threes$location <- "from_stop"
    cds <- rbind(fives, threes)
    if (plot) {
        print(ggplot(cds, aes(x = codon, y = value, fill = frame)) +
            geom_bar(stat = "identity", position = "dodge") +
            facet_grid(. ~ location, scales = "free_x"))
    }
    invisible(cds)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
