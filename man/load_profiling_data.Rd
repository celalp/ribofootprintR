\name{load_profiling_data}
\alias{load_profiling_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
re-load pre-processed profiling data
%%  ~~function to do ... ~~
}
\description{
load a pre-processed profiling data as an Rlist
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load_profiling_data(file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
  name of the Rdata object
%%     ~~Describe \code{file_name} here~~
}
}
\details{
load the Rdata back to the enviroment. This is a crude function. It just loads an Rdata file back to \code{Global Environment}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
profiling list. if the Rdata file is not a list might give an error. The output is whatever is saved with \code{\link{write_profiling_data}}
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
I need to make this a profiling specific function.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{write_profiling_data}}, \code{\link{read_profiling_data}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file_name)
{
    data <- list.load(file_name)
    data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
