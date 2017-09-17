\name{write_profiling_data}
\alias{write_profiling_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write profiling library list to file.
%%  ~~function to do ... ~~
}
\description{
save the profiling library to file as and Rdata object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write_profiling_data(data, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Rlist to be written
%%     ~~Describe \code{data} here~~
}
  \item{name}{
  name of the file ".Rdata" will be pasted to the file name.
%%     ~~Describe \code{name} here~~
}
}
\details{
this is a one line function for convenience. It's intended function is to ease movement of data between different computers. The resulting Rdata object is several orders of magnitude smaller than the .bam files.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
I need to make this a profiling specific function.

%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{load_profiling_data}}, \code{\link{read_profiling_data}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, name)
{
    list.save(x = data, file = paste(name, "rdata", sep = "."))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
