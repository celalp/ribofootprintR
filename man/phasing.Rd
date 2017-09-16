\name{phasing}
\alias{phasing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count the number of aligned reads per read length per frame

%%  ~~function to do ... ~~
}
\description{
This function counts the number of reads in the [[ribo]] list for all genes and returns a data frame with the number of reads per read lenght per frame.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phasing(data, genedf, cores = 1, plot = F)
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
  \item{cores}{
  integer, how many threads to use. See details.
%%     ~~Describe \code{cores} here~~
}
  \item{plot}{
  logical whether to plot the results, defaults to F.
%%     ~~Describe \code{plot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a molten data frame that has the number of reads per read lenght per frame.

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
function (data, genedf, cores = 1, plot = F)
{
    genes <- genedf$gene
    ribo <- data[["ribo"]]
    get_phasing <- function(gene) {
        orf <- ribo[[gene]]
        if (is.null(orf)) {
            orf <- NULL
            orf
        }
        else {
            features <- genedf[genedf$gene == gene, ]
            orf <- orf[orf$codon > 0 & orf$codon < features$ncodon,
                c(1, 2, 4)]
            orf
        }
    }
    if (cores > 1) {
        orfs <- do.call("rbind", mclapply(genes, get_phasing,
            mc.cores = cores))
    }
    else {
        orfs <- do.call("rbind", lapply(genes, get_phasing))
    }
    phasing <- reshape2::dcast(data = orfs, width ~ frame, fun.aggregate = sum)
    colnames(phasing) <- c("read_length", "frame_0", "frame_1",
        "frame_2")
    phasing$read_length <- as.factor(phasing$read_length)
    phasing <- reshape2::melt(phasing)
    colnames(phasing) <- c("read_length", "frame", "number_of_reads")
    if (plot) {
        print(ggplot(phasing, aes(x = read_length, y = number_of_reads,
            fill = frame)) + geom_bar(stat = "identity", position = "dodge"))
    }
    invisible(phasing)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
