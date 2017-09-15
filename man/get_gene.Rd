\name{get_gene}
\alias{get_gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get coverage and ribosome occupancy of a gene by nucleotide
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns a data.frame with RNA-Seq coverage per nucleotide and This is workhorse function in the package, it is called by most other functions.
}
\usage{
get_gene(gene_name, data, plot = F, seq = NA, genedf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_name}{
  character name of the gene, must match genedf$gene entry
%%     ~~Describe \code{gene_name} here~~
}
  \item{data}{
  list returned by read_profiling_data() function
%%     ~~Describe \code{data} here~~
}
  \item{plot}{
  Whether to plot the data frame, default=F
%%     ~~Describe \code{plot} here~~
}
  \item{seq}{
  DNAString set object of the FASTA file used for alignment. This is only needed for the occupancy function otherwise can be set to NA
%%     ~~Describe \code{seq} here~~
}
  \item{genedf}{
  a data frame that contains the gene name, and the nucleotide coordinates of start and stop locations as well as number of nucleotites and codon of the ORF.
%%     ~~Describe \code{genedf} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame with the nucleotide number (negative if within the 5' UTR), RNA-Seq coverage (NA is no RNA-Seq library is present) and the number of reads that map to each codon (and their frame). See \code{\link{read_profiling_data}} for the frame calculations
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
\code{\link{read_profiling_data}}, \code{\link{phasing}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
TODO


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
            gene <- data.frame(nucleotide = c(-coord$start:(width(seq[names(seq) ==
                gene_name]) - coord$start)), coverage = 1)
            gene <- full_join(gene, ribo, by = "nucleotide")
            gene$frame <- as.factor(gene$frame)
            if (plot) {
                plot <- ggplot() + geom_bar(data = gene[!is.na(gene$frame),
                  ], aes(x = nucleotide, y = freq, fill = frame),
                  stat = "identity", position = "dodge")
                plot <- plot + geom_segment(aes(x = 0, y = -1,
                  xend = (coord$end - coord$start) + 1, yend = -1),
                  color = "blueviolet", size = 2)
                print(plot)
            }
        }
        else {
            gene <- dplyr::full_join(rna, ribo, by = "nucleotide")
            gene$frame <- as.factor(gene$frame)
            if (plot) {
                plot <- ggplot() + geom_line(data = gene, aes(x = nucleotide,
                  y = log2(coverage))) + geom_bar(data = gene[!is.na(gene$frame),
                  ], aes(x = nucleotide, y = freq, fill = frame),
                  stat = "identity", position = "dodge")
                plot <- plot + geom_segment(aes(x = 0, y = -1,
                  xend = coord$end - coord$start, yend = -1),
                  color = "blueviolet", size = 2)
                print(plot)
            }
        }
    }
    gene
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~gene }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~data }% __ONLY ONE__ keyword per line
