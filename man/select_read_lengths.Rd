\name{select_read_lengths}
\alias{select_read_lengths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select Relevant Read Lengths for the Profiling Library
}
\description{
Filters out reads that are too long or too short for meaningful biological analysis. For example is cycloheximide is used in the experiment read length 27-29 are usually the predominant read lengths and provide the most biologically meaningful data.
}
\usage{
select_read_lengths(data, ribo_length = c(27:29))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  list returned by code{\link{read_profiling_data}}
%%     ~~Describe \code{data} here~~
}
  \item{ribo_length}{
  read lengths that are selected. defaults to c(27:29)
%%     ~~Describe \code{ribo_length} here~~
}
}
\details{
selects desired read lengths. This is often necessary before further analysis.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a list similar to \code{\link{read_profiling_data}} except the \code{[[ribo]]} portion only contains read lengths selected. The \code{[[rna]]} portion remains untouched.
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
function (data, ribo_length = c(27:29))
{
    ribo <- data[["ribo"]]
    rna <- data[["rna"]]
    features <- data[["features"]]
    selector <- function(gene, lengths) {
        gene <- ribo[[gene]]
        gene <- gene[gene$width >= min(lengths) & gene$width <=
            max(lengths), ]
        gene
    }
    ribo_selected <- lapply(names(ribo), selector, lengths = ribo_length)
    names(ribo_selected) <- names(ribo)
    data$ribo <- ribo_selected
    data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
