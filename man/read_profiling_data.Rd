\name{read_profiling_data}
\alias{read_profiling_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read .bam files and pre-process them for further analysis
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read_profiling_data(ribo, ribo_index, rna = NA, rna_index = NA, genedf, cores = 1, A_site_adjust = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ribo}{
  sorted and indexed .bam file for ribosome profiling libraries

%%     ~~Describe \code{ribo} here~~
}
  \item{ribo_index}{
  index file for ribosome profiling
%%     ~~Describe \code{ribo_index} here~~
}
  \item{rna}{
  sorted indexed .bam file for RNA-Seq control, can be ignored defaults to NA

%%     ~~Describe \code{rna} here~~
}
  \item{rna_index}{
  index file for RNA-Seq control, can be ignored defaults to NA

%%     ~~Describe \code{rna_index} here~~
}
  \item{genedf}{
  a data frame with the start and stop coordinated of the ORF. DO NOT USE A GTF OR A GFF FILE. This data frame can be generated through \code{generate_transcriptome}, function or manually.
%%     ~~Describe \code{genedf} here~~
}
  \item{cores}{
  integer, number of threads to use. See details.
%%     ~~Describe \code{cores} here~~
}
  \item{A_site_adjust}{
  Whether to adjust frames and start positions to the ribosome A site. See details. Defaults to T
%%     ~~Describe \code{A_site_adjust} here~~
}
}
\details{

For ribosome profiling libraries the frame and is calculated by using the formula below

\deqn{mapped nucleotide = 1 based read start - orf start (based on genedf)}

the frame is the modulus of

\deqn{frame = mapped nucleotide - ORF-start + 1 + (read length - 15)}

when divided by 3. This is done to make sure that all usable read lengths (18-30) all show predominantly 0 frame when there is good phasing. However this assumes that RNAseI is used for digestion and therefore is not compatible with MNase treated samples. For those samples please set \code{A_site_adjust=F}

For RNA-Seq libraries simple coverage is calculated by using the \code{coverage} function of Rsamtools.

Currently PSOCK multithreadding is not supported. For Windows set cores to 1. For operating systems that support forking the number of threads that are used can be set by the cores argument. If the cores>1 then the data is processed by mclapply instead of lapply.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of lists. The returned value is a list with two elements [[ribo]] and [[rna]] each of which are lists themselves. Within each sub-list each gene is a separate data frame. For profiling the columns are nucleotide # (- if in 5' UTR), frame, codon (- if in 5' UTR), read length and the number of reads that map for that frame and read_length. For RNA-Seq data, just nucleotide and coverage. If there is no RNA-Seq data then the RNA portion of the list is set to NULL
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
\code{\link{generate_transcriptome}}, \code{\link{phasing}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ribo, ribo_index, rna = NA, rna_index = NA, genedf,
    cores = 1, A_site_adjust = T)
{
    message("Loading ribo data")
    read_data <- function(data, index) {
        message("read bam file")
        gr <- readGAlignments(file = data, index = index)
        invisible(gc())
        message("Split into individual genes")
        gr <- split(gr, seqnames(gr))
        invisible(gc)
        return(gr)
    }
    get_coverage <- function(gene) {
        covs <- cov_list[[gene]]
        if (is.null(covs)) {
            NULL
        }
        else {
            covs <- as.data.frame(covs)
            b <- genedf[genedf$gene == gene, ]
            coverage <- data.frame(nucleotide = as.integer(1:nrow(covs)),
                coverage = covs$value)
            coverage$nucleotide <- coverage$nucleotide - b$start
            coverage
        }
    }
    get_frames <- function(gr) {
        if (is.null(gr)) {
            gr <- data.frame()
        }
        else {
            b <- genedf[genedf$gene == as.character(gr@seqnames@values),
                ]
            nucleotide <- gr@start - b$start
            if (A_site_adjust) {
                frames <- data.frame(frame = nucleotide - b$start +
                  1 + (width(gr) - 15))\%\%3
                frames$codon <- floor((start(gr) - b$start +
                  (round(width(gr)/2)))/3) + 1
                frames$width <- width(gr)
            }
            else {
                frames <- data.frame(frame = ((nucleotide)\%\%3))
                frames$codon <- floor(nucleotide/3)
            }
            frames$width <- width(gr)
            frames$nucleotide <- frames$codon * 3 - 2
            frames <- plyr::count(frames)
            frames <- frames[, c(4, 1, 2, 3, 5)]
            frames
        }
    }
    ribo <- read_data(data = ribo, index = ribo_index)
    message("Calling Frames")
    if (cores > 1) {
        ribo <- mclapply(ribo, get_frames, mc.cleanup = T, mc.cores = cores)
    }
    else {
        ribo <- lapply(ribo, get_frames)
    }
    invisible(gc())
    message("Done")
    if (is.na(rna)) {
        warning("No rna data")
        coverages <- NULL
    }
    else {
        message("Generating coverage list")
        cov_list <- coverage(rna)
        message("Done")
        message("Split into individual genes")
        genes <- genedf$gene
        coverages <- lapply(X = genes, FUN = get_coverage)
        names(coverages) <- genedf$gene
        invisible(gc())
        message("Done")
    }
    invisible(gc)
    data <- list(ribo = ribo, rna = coverages)
    invisible(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
