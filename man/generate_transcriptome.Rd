\name{generate_transcriptome}
\alias{generate_transcriptome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate transcriptome and genedf for alignment and analysis from BSGenome and TxDb objects
}
\description{
Returns a list whose items are the transcriptome for reads to be aligned as a DNAStringset and a dataframe that contains the start and end locations of open reading frames.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate_transcriptome(genome, TxDb, fill_utr = F, utr_fill_length = 30, cores = 1, write = T, fasta_name, genedf_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{
  A BSGenome object from bioconductor of self made from a genome fasta file
%%     ~~Describe \code{genome} here~~
}
  \item{TxDb}{
  A TxDb database, can be from bioconductor or generated from a GFF3 file.
%%     ~~Describe \code{TxDb} here~~
}
  \item{fill_utr}{
  If the genome (or a gene) does not have UTR annotations use a fixed length of sequence that if immediately after or before the start location. Defaults to F.
%%     ~~Describe \code{fill_utr} here~~
}
  \item{utr_fill_length}{
  The length of the UTR to be filled. Defaults to 30. See details for potential problems.
%%     ~~Describe \code{utr_fill_length} here~~
}
  \item{cores}{
  number of threads to use in UNIX based systems. Defaults to 1.
%%     ~~Describe \code{cores} here~~
}
  \item{write}{
  whether to write the DNAStringset and the data frame to file separately.
%%     ~~Describe \code{write} here~~
}
  \item{fasta_name}{
  the name of the fasta file (.fasta will be added).
%%     ~~Describe \code{fasta_name} here~~
}
  \item{genedf_name}{
  name of the genedf file saved as a tab delimeted plain text file (.txt will be added)
%%     ~~Describe \code{genedf_name} here~~
}
}
\details{
Returns a list where the first item is a DNAStringset object with genomic sequences for each CDS and their UTRs combined into a single string per gene. The second item is a data frame where the start and end locations of the open reading frame as well as number of codons and nucletoidies are stored. These start and stop coodrinates are relative to the fasta, they are not genomic coordinates.

For genomes without UTR annotations there is a fill_utr option and the lenght of the "UTR" is specified in the utr_fill_length. UTR filling is done when there are no 5' or 3' UTRs are specified in the genome. Selecting a large number for this can be problematic as it can bleed into other coding sequences or introns that are in UTRs (in S.Cereviase genome for example there are not UTR annotations but 24 5' UTR introns)

Currently PSOCK multi threading is not supported. For Windows set cores to 1. For operating systems that support forking the number of threads that are used can be set by the cores argument. If the cores>1 then the data is processed by mclapply instead of lapply.

For larger genomes it may take several minutes for the code to finish runnig. However this process is only needed to be performed once per genome.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
List
first item is a DNAstringset, second item is a data.frame
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
\code{\link{read_profiling_data}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (genome, TxDb, fill_utr = F, utr_fill_length = 300,
    cores = 1, write = T, fasta_name, genedf_name)
{
    if (keep_introns) {
        warning("If using a bioconductor TxDb object some introns may be ignored try creating your own from a GFF file",
            call. = F)
    }
    cds <- cdsBy(TxDb, use.names = T)
    gene_names <- names(cds)
    five_UTR <- fiveUTRsByTranscript(TxDb, use.names = T)
    three_UTR <- threeUTRsByTranscript(TxDb, use.names = T)
    chrlengths <- seqlengths(genome)
    get_ranges <- function(gene_name) {
        gene_cds <- as.data.frame(cds[[gene_name]])
        gene_cds$class <- "cds"
        if (gene_cds$strand[1] == "+") {
            gene_cds <- gene_cds[order(gene_cds$start), ]
        }
        else {
            gene_cds <- gene_cds[order(gene_cds$start, decreasing = T),
                ]
        }
        gene_threeUTR <- as.data.frame(three_UTR[[gene_name]])
        gene_fiveUTR <- as.data.frame(five_UTR[[gene_name]])
        if (dim(gene_threeUTR)[1] == 0 && fill_utr) {
            size <- dim(gene_cds)[1]
            if (gene_cds$strand[1] == "+") {
                gene_threeUTR <- gene_cds[size, ]
                gene_threeUTR$start <- gene_threeUTR$end + 1
                gene_threeUTR$end <- gene_threeUTR$start + utr_fill_length
                if (gene_threeUTR$end > chrlengths[names(chrlengths) ==
                  gene_cds$seqnames[1]]) {
                  gene_threeUTR$end <- chrlengths[names(chrlengths) ==
                    gene_cds$seqnames[1]]
                }
            }
            else {
                gene_threeUTR <- gene_cds[size, ]
                gene_threeUTR$end <- gene_threeUTR$start - 1
                gene_threeUTR$start <- gene_threeUTR$start -
                  utr_fill_length
                if (gene_threeUTR$start < 1) {
                  gene_threeUTR$start <- 1
                  gene_threeUTR$end <- gene_cds$start[1] - 1
                }
            }
        }
        if (dim(gene_fiveUTR)[1] == 0 && fill_utr) {
            if (gene_cds$strand[1] == "+") {
                gene_fiveUTR <- gene_cds[1, ]
                gene_fiveUTR$end <- gene_fiveUTR$start - 1
                gene_fiveUTR$start <- gene_fiveUTR$end - utr_fill_length
                if (gene_fiveUTR$start < 1) {
                  gene_fiveUTR$start <- 1
                  gene_fiveUTR$end <- gene_cds$start[1] - 1
                }
            }
            else {
                gene_fiveUTR <- gene_cds[1, ]
                gene_fiveUTR$start <- gene_fiveUTR$end + 1
                gene_fiveUTR$end <- gene_fiveUTR$start + utr_fill_length
                if (gene_fiveUTR$end > chrlengths[names(chrlengths) ==
                  gene_cds$seqnames[1]]) {
                  gene_fiveUTR$end <- chrlengths[names(chrlengths) ==
                    gene_cds$seqnames[1]]
                }
            }
        }
        if (dim(gene_fiveUTR)[1] > 0) {
            gene_fiveUTR$class <- "five_UTR"
            gene <- rbind.fill(gene_cds, gene_fiveUTR)
        }
        else {
            gene <- gene_cds
        }
        if (dim(gene_threeUTR)[1] > 0) {
            gene_threeUTR$class <- "three_UTR"
            gene <- rbind.fill(gene, gene_threeUTR)
        }
        else {
            gene <- gene_cds
        }
        if (gene$strand[1] == "+") {
            gene <- gene[order(gene$start), ]
        }
        else {
            gene <- gene[order(gene$start, decreasing = T), ]
        }
        gene$width <- gene$end - gene$start
        gene$cds_name <- gene_name
        gene <- gene[c(1, 2, 3, 4, 5, 7, 9)]
        gene
    }
    if (cores > 1) {
        gene_list <- mclapply(gene_names, FUN = get_ranges, mc.cores = cores)
    }
    else {
        gene_list <- lapply(gene_names, get_ranges)
    }
    names(gene_list) <- gene_names
    get_seq <- function(gene_name) {
        gene <- gene_list[[gene_name]]
        gene_gr <- makeGRangesFromDataFrame(gene)
        seq <- getSeq(x = genome, names = gene_gr)
        names(seq) <- gene$class
        seq
    }
    if (cores > 1) {
        seq_list <- mclapply(gene_names, FUN = get_seq, mc.cores = cores)
    }
    else {
        seq_list <- lapply(gene_names, get_seq)
    }
    names(seq_list) <- gene_names
    make_genedf <- function(gene_name) {
        gene <- seq_list[[gene_name]]
        five_w <- width(gene[names(gene) == "five_UTR"])
        three_w <- width(gene[names(gene) == "three_UTR"])
        cds_w <- sum(width(gene[names(gene) == "cds"]))
        gene_df <- data.frame(gene = gene_name, start = five_w +
            1, end = sum(cds_w + five_w), ncodons = cds_w/3,
            nnuc = cds_w)
        gene_df
    }
    genedf <- do.call("rbind", lapply(gene_names, make_genedf))
    fasta <- lapply(seq_list, unlist)
    fasta <- DNAStringSet(fasta)
    support <- list(fasta = fasta, genedf = genedf)
    if (write) {
        writeXStringSet(fasta, paste(fasta_name, "fasta", sep = "."),
            append = FALSE, compress = FALSE, compression_level = NA,
            format = "fasta")
        write.table(genedf, paste(genedf_name, "txt", sep = "."),
            quote = F, sep = "\t", col.names = T, row.names = F)
    }
    invisible(support)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
