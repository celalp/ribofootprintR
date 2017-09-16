\name{get_gene}
\alias{get_gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get coverage and ribosome occupancy of a gene by nucleotide
}
\description{
Returns a data.frame with RNA-Seq coverage per nucleotide and This is workhorse function in the package, it is called by most other functions.
}
\usage{
get_gene(gene_name, data, plot = F, seq = NA, genedf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_name}{
character name of the gene, must match genedf$gene entry
}
  \item{data}{
list returned by read_profiling_data() function
}
  \item{plot}{
whether to plot the data frame, default=F
}
  \item{seq}{
DNAString set object of the FASTA file used for alignment. This is only needed for the occupancy function otherwise can be set to NA. Defaults to NA
}
  \item{genedf}{
a data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame with the nucleotide number (negative if within the 5' UTR), RNA-Seq coverage (NA is no RNA-Seq library is present) and the number of reads that map to each codon (and their frame). See read_profiling_data for the frame calculations
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
function (gene_name, data, plot = F, seq = NA, genedf)
{
    ribo <- data[["ribo"]][[gene_name]]
    rna <- data[["rna"]][[gene_name]]
    if (is.null(ribo) == T) {
        gene <- NULL
    }
    else {
        coord <- genedf[genedf$gene == gene_name, ]
        coord <- coord[coord[, 1] == gene_name, ]
        if (is.null(rna) == T) {
            gene <- ribo
            gene$frame <- as.factor(gene$frame)
            if (plot) {
                plot <- ggplot() + geom_bar(data = gene, aes(x = nucleotide,
                  y = freq, fill = frame), stat = "identity",
                  position = "dodge")
                plot <- plot + geom_segment(aes(x = 0, y = 0,
                  xend = (coord$end - coord$start) + 1, yend = 0),
                  color = "blueviolet", size = 2)
                print(plot)
            }
        }
        else {
            gene <- dplyr::full_join(rna, ribo, by = "nucleotide")
            gene$frame <- as.factor(gene$frame)
            if (plot) {
                plot <- ggplot() + geom_line(data = gene, aes(x = nucleotide,
                  y = -coverage)) + geom_bar(data = gene[!is.na(gene$frame),
                  ], aes(x = nucleotide, y = freq, fill = frame),
                  stat = "identity", position = "dodge")
                plot <- plot + geom_segment(aes(x = 0, y = 0,
                  xend = coord$end - coord$start, yend = 0),
                  color = "blueviolet", size = 2)
                print(plot)
            }
        }
    }
    invisible(gene)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
