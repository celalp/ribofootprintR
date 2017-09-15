\name{smith}
\alias{smith}
\docType{data}
\title{
Ribosome profiiling and RNA-Seq library from Smith et. al.
%%   ~~ data name/kind ... ~~
}
\description{
Ribosome profiling and RNA-Seq libraries from WT S.Cerevisiae.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("smith")}
\format{
 This is the results of the \code{read_profiling_data} function. In addition to the n and rna elements the genedf and DNAStringset is present in \code{smith$genedf} and \code{smith$seq}
 }
\details{
These libraries are prepared in the presence of cycloheximide. The data is published by Smith et. al. PMID:  24931603. The raw data can be accessed at SRR125839 for ribosome profiling and SRR12558470. The aligment was done to an earlier version of yeast transcriptome that was generated manualy. A similar transcriptome and genedf can be generated using \code{generate_transcriptome} function.

Alignment was done with:


\code{lib=$(find -iname "*.fastq")}

\code{for i in ${lib[*]}}

\code{    do}

\code{	   cutadapt -a CTGTAGGCA -q 10 --trim-n -m 10 -o $i-trimmed.fastq $i}

\code{       bowtie -m 100 -p 8 -n 2 -l 15 --un $i-norna transcriptome/rna $i-trimmed.fastq}

\code{       bowtie -m 4 -p 8 -n 2 -l 15 -S --best --strata transcriptome/profiling $i-norna > $i.sam}

\code{	    samtools view -b $i.sam > $i.bam}

\code{        samtools sort $i.bam > $i-sorted.bam}

\code{        samtools index -b $i-sorted.bam > $i-sorted.bam.bai}

\code{done}


%%  ~~ If necessary, more details than the __description__ above ~~
}

\source{
\link{https://www.ncbi.nlm.nih.gov/pubmed/24931603}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Cell Rep. 2014 Jun 26;7(6):1858-66
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(smith)

}
\keyword{datasets}
