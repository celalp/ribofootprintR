\name{metagene_cds}
\alias{metagene_cds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get average ribosome footprint alignment to the codin sequence at the 5' and the 3' end
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metagene_cds(data, genes, before_start = 50, after_start = 100, before_stop = 100, after_stop = 50, norms = "total", seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  list retuned by \code{read_profiling_data}
%%     ~~Describe \code{data} here~~
}
  \item{genes}{
   data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF.
%%     ~~Describe \code{genes} here~~
}
  \item{before_start}{
  the number of nucleotide to calculate before the start nucleotide, defaults to 50
%%     ~~Describe \code{before_start} here~~
}
  \item{after_start}{
   the number of nucleotide to calculate after the start nucleotide, defaults to 100
%%     ~~Describe \code{after_start} here~~
}
  \item{before_stop}{
     the number of nucleotide to calculate before the stop nucleotide, defaults to 100

%%     ~~Describe \code{before_stop} here~~
}
  \item{after_stop}{
       the number of nucleotide to calculate before the stop nucleotide, defaults to 50

%%     ~~Describe \code{after_stop} here~~
}
  \item{norms}{
  normalization, character either "total", "rna" or "none". Defaults to "total"
%%     ~~Describe \code{norms} here~~
}
  \item{seq}{
  DNAString set object of the FASTA file used for alignment. This is only needed for the occupancy function otherwise can be set to NA
%%     ~~Describe \code{seq} here~~
}
}
\details{
If norms is set to "total"" returns the average fraction of reads that map to the region specified. If norms are set to "rna" then the number of reads that map to each codon are divided by the RNA-Seq coverage. "none" returns the average number of reads per codon.
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
\code{\link{metagene_bins}}, \code{\link{totals}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, genes, before_start = 50, after_start = 100,
    before_stop = 100, after_stop = 50, norms = "total", seq)
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
    gene_names <- genes$gene
    get_cds <- function(gene_name, direction, norm = norms) {
        gene <- get_gene(gene_name, data = data, seq = seq)
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
    fives <- do.call("rbind", mclapply(gene_names, get_cds, "five",
        mc.cores = (detectCores() - 1)))
    x <- unique(fives[, c(1, 3)])
    counts <- mdply(x, returner, fives)
    fives <- dplyr::left_join(unique(fives[, c(1, 2)]), counts,
        by = "nucleotide")
    colnames(fives)[4] <- "value"
    fives$mean_value <- fives$value/length(gene_names)
    threes <- do.call("rbind", mclapply(gene_names, get_cds,
        "three", mc.cores = (detectCores() - 1)))
    x <- unique(threes[, c(1, 3)])
    counts <- mdply(x, returner, threes)
    threes <- dplyr::left_join(unique(threes[, c(1, 2)]), counts,
        by = "nucleotide")
    colnames(threes)[4] <- "value"
    threes$mean_value <- threes$value/length(gene_names)
    fives$location <- "from_start"
    threes$location <- "from_stop"
    cds <- rbind(fives, threes)
    cds
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
