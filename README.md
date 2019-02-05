# Updates to the original package

There is now a new function called rrp, short for relative ribosome positioning. This function calculates where withing the ORF the percentage of the reads specified goes up to. For example if the specified percentage is 50, the function calculates (in it's own way, see function doc for details) where within the ORF the cumulatibe 50% of the reads reside. This was originally concieved by Richard Baker at UMASS Medical School. Though conceptually similar the way these numbers are calculated are quite different from his original script. 

### Coming soon!
A function based on `plyr::mdply` as a wrapper to compare different runs of the same function. 
A wrapper function to run all the functions in the package and annotate them as `attr` for later use. 


# Introduction

This is an extremely early version of the package I'm in the process of building. The idea is to provide users with a tool to analyse ribosome footprint profiling libraries. These libraries are generated from RNase treated cell lysates. The idea is to digest single stranded sections of mRNA so that only the ribosome protected portions will be available for sequencing. This way we can measure the protein production dynamics of all mRNAs simultaneously. 

In addition to the .bam files the user also needs a FASTA file and a data.frame that contains the start and stop locations of all the genes of interest. These are NOT genomes or GFT/GFF files. In order to circumvent splicing problems one needs a "transcriptome" to align their reads to. The function to make this "transcriptome" and genedf use `generate_transcriptome()` function. 

Like most R packages this one also has a few dependencies. The R packages that ribofootprintR depends on are as follows:
GenomicAlignments, ggplot2, plyr, dplyr, Rsamtools, zoo, GenomicFeatures, reshape2, rlist and Biostrings. Please make sure that you are running a recent version of R (>=3.4) and all packages are up to date. 

To install the package please use the snippet

```{r, eval=FALSE, include=T}
source("https://bioconductor.org/biocLite.R")
biocLite(pkgs = c("GenomicFeatures", "Rsamtools", "GenomicAlignments"))
install.packages(c("rlist", "ggplot2", "plyr", "dplyr", "zoo", "reshape2"))
devtools::install_github("celalp/ribofootprintR")
```


In addition for the processing of .bam files it depends on [samtools](http://www.htslib.org/) library so make sure that samtools is in your $PATH.

## Loading the Data
Before we get to do any analysis we need to load the data. 

```{r, eval=FALSE, include=T}
library(ribofootprintR)
```

```{r, eval=FALSE, include=T}
data<-read_profiling_data(ribo = "ribo-sorted.bam", ribo_index = "ribo-sorted.bam.bai", 
                          rna="rna-sorted.bam", rna_index = "rna-sorted.bam.bai", 
                          genedf = genedf)
```

## Analysis
One important aspect of ribosome profiling data is that one first needs to determine whether profiling libraries actually have the needed three nucleotide periodicity. This can be done by using the `phasing()` function. 

```{r, eval=FALSE, include=T}
phasing<-phasing(data = data, genedf = genedf, plot = T)
```

Usually the only that are 28nt (and sometimes 27 and 29 if cycloheximide is used) long show substantial periodicity, in that they are mostly in one frame. We may need to select these read lengths before we continue with detailed analysis. For this we can use `select_read_lengths()` function

```{r, eval=FALSE, include=T}
selected<-select_read_lengths(data = data, ribo_length = 28)
```

We can of course look at this from an aggregate level as well. For example we can look at the first 100 nucleotides after the start and 100 nucleotide before the stop codon. 


```{r, eval=FALSE, include=T}
cds<-metagene_cds(data = selected, norms ="total", plot = T, genedf = genedf)
```

Considering the sizes of open reading frames can vary wildly may be looking at the first and last couple of hundred nucleotides is not enough. One way to look at aggregate data is to to divide the ORF into different bins and normalize all of the ORF. That can be done using the `metagene_bins()` function. 

```{r, eval=FALSE, include=T}
bins<-metagene_bins(data = selected, genedf = genedf, bins = 100, norms = "total", simplify = "none", plot = T)
```

Maybe we are just interested in looking at how many ribosomes are on a given mRNA and how does that compare to mRNA expression levels. `totals()` function does that.

```{r, eval=FALSE, include=T}
tots<-totals(data = selected, genedf = genedf, plot = T)
```

Finally, we may be interested in ribosome dynamics and how each codon or amino acid is translated. To look at those measures we can use the `occupancy` function. This is the longest function in the package and depending on the size of the data might take several minutes to run. 

```{r, eval=FALSE, include=T}
occ<-occupancy(data = selected, seq = seq, genedf = genedf)
```

Portability and data sharing is an important feature of research. So after some analysis the data can be written to disk with `write_profiling_data` and loaded back with `load_profiling_data` functions. These save the data as an .Rdata object and they usually are <100MB per library which makes them easy to share over dropbox or other file sharing systems. 

## Other features

There are many normalization and display options in each function. You can see the function documentation for details. For this intial iteration I left the entire R code in the package documentation for other people to see easily. 

All the calculations we did above used the entire gene set. We can just as easily re-run these calculations with any subset of genedf. For example if one wants to compare substrates of one pathway to the rest of the transcriptome we can just re-run the code twice and compare the resulting data.frames. We can also change the start and end coordinates of the genedf to focus on specific areas of specific genes. I tried to make the requirements for data analysis as little as possible so complicated comparisons can be done with little preparation. 

Almost all the functions in the package invisibly return a data.frame (except for `read_profiling_data` and `select_read_lengths`) so the end user can perform more complicated analyses if desired. 

While .bam files can be several GB in size the processing done by `read_profiling_data` condenses the data and removes a lot of redundancies. After initial processing the size of the Rlist that is produced is usually around 150MB. This allows the user to be able to load many experiments simultaneously and do rapid processing without needing to worry about running out of RAM. 

Most of the functions in the package also contain a cores argument. For operating systems that allow forking (iOS, Linux) the calculations are passed to mclappy (instead of lapply) and uses user specified number of threads I opted against using detectCores not to violate any job scheduling system settings and if one is doing the calculations on a virtual machine the detectCores function may return `NA`. To enable forking the user needs to set the cores argument >1. Currently there is no windows multicore processing because when a PSOCK cluster is initialized in R all the necessary files would also be available in that cluster. Currently there is no check for whether the operating system is suitable for forking so it is up to the user to make sure. 

### In the future

For this to package to be bioconductor compatible I need also to implement S4 methods or at the very least convert the initial Rlist that is the output of `read_profiling_data` to be an S4 object. 

Finally I would like to be able to use the posterior alignment probabilities generated by [RSEM's](https://deweylab.github.io/RSEM/) expectation maximization protocol for reads that map to multiple locations. 

Please let me know if you have other suggestions and detect any bugs or errors. 
