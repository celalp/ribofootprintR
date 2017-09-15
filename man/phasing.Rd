\name{phasing}
\alias{phasing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
count the number of reads per read length per frame
%%  ~~function to do ... ~~
}
\description{
this function counts the number of reads in the \code{[[ribo]]} list for all genes and returns a data frame with the number of reads per read lenght per frame.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phasing(data, genedf)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a molten data frame that has the number of reads per read lenght per frame

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
\code{\link{read_profiling_data}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, genedf)
{
    genes <- genedf$gene
    ribo <- data[["ribo"]]
    get_phasing <- function(gene) {
        orf <- ribo[[gene]]
        if (is.null(orf)) {
            orf <- NULL
        }
        else {
            features <- genedf[genedf$gene == gene, ]
            orf <- orf[orf$codon > 0 & orf$codon < features$ncodon,
                c(1, 2, 4)]
            orf
        }
    }
    orfs <- do.call("rbind", mclapply(genes, get_phasing, mc.cores = (detectCores() -
        1)))
    phasing <- reshape2::dcast(data = orfs, width ~ frame, fun.aggregate = sum)
    colnames(phasing) <- c("read_length", "frame_0", "frame_1",
        "frame_2")
    phasing$read_length <- as.factor(phasing$read_length)
    phasing <- reshape2::melt(phasing)
    colnames(phasing) <- c("read_length", "frame", "number_of_reads")
    phasing
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
