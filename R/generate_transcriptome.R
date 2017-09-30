generate_transcriptome<-function(genome, TxDb, fill_utr=F, utr_fill_length=300, cores=1, write=T, fasta_name, genedf_name){
  if(keep_introns){
    warning("If using a bioconductor TxDb object some introns may be ignored try creating your own from a GFF file", call. = F)
  }
  cds<-cdsBy(TxDb, use.names=T)
  gene_names<-names(cds)
  five_UTR<-fiveUTRsByTranscript(TxDb, use.names=T)
  three_UTR<-threeUTRsByTranscript(TxDb, use.names=T)
  chrlengths<-seqlengths(genome)

  #get coordinates
  #fix chromosome edges so that utr fill does not end up in the next thing
  get_ranges<-function(gene_name){
    gene_cds<-as.data.frame(cds[[gene_name]])
    gene_cds$class<-"cds"
    if(gene_cds$strand[1]=="+"){
      gene_cds<-gene_cds[order(gene_cds$start),]
    } else {
      gene_cds<-gene_cds[order(gene_cds$start, decreasing = T),]
    }
    gene_threeUTR<-as.data.frame(three_UTR[[gene_name]])
    gene_fiveUTR<-as.data.frame(five_UTR[[gene_name]])

    if(dim(gene_threeUTR)[1]==0 && fill_utr){
      size<-dim(gene_cds)[1]
      if(gene_cds$strand[1]=="+"){
        gene_threeUTR<-gene_cds[size,]
        gene_threeUTR$start<-gene_threeUTR$end+1
        gene_threeUTR$end<-gene_threeUTR$start+utr_fill_length
        if(gene_threeUTR$end>chrlengths[names(chrlengths)==gene_cds$seqnames[1]]){
          gene_threeUTR$end<-chrlengths[names(chrlengths)==gene_cds$seqnames[1]]
        }
      } else {
        gene_threeUTR<-gene_cds[size,]
        gene_threeUTR$end<-gene_threeUTR$start-1
        gene_threeUTR$start<-gene_threeUTR$start-utr_fill_length
        if(gene_threeUTR$start<1){
          gene_threeUTR$start<-1
          gene_threeUTR$end<-gene_cds$start[1]-1
        }
      }
    }

    if(dim(gene_fiveUTR)[1]==0 && fill_utr){
      if(gene_cds$strand[1]=="+"){
        gene_fiveUTR<-gene_cds[1,]
        gene_fiveUTR$end<-gene_fiveUTR$start-1
        gene_fiveUTR$start<-gene_fiveUTR$end-utr_fill_length
        if(gene_fiveUTR$start<1){
          gene_fiveUTR$start<-1
          gene_fiveUTR$end<-gene_cds$start[1]-1
        }
      } else {
        gene_fiveUTR<-gene_cds[1,]
        gene_fiveUTR$start<-gene_fiveUTR$end+1
        gene_fiveUTR$end<-gene_fiveUTR$start+utr_fill_length
        if(gene_fiveUTR$end>chrlengths[names(chrlengths)==gene_cds$seqnames[1]]){
          gene_fiveUTR$end<-chrlengths[names(chrlengths)==gene_cds$seqnames[1]]
        }
      }
    }

    if(dim(gene_fiveUTR)[1]>0){
      gene_fiveUTR$class<-"five_UTR"
      gene<-rbind.fill(gene_cds, gene_fiveUTR)
    } else {
      gene<-gene_cds
    }
    if(dim(gene_threeUTR)[1]>0){
      gene_threeUTR$class<-"three_UTR"
      gene<-rbind.fill(gene, gene_threeUTR)
    } else {
      gene<-gene_cds
    }
    if(gene$strand[1]=="+"){
      gene<-gene[order(gene$start),]
    } else {
      gene<-gene[order(gene$start, decreasing = T),]
    }
    gene$width<-gene$end-gene$start
    gene$cds_name<-gene_name
    gene<-gene[c(1,2,3,4,5,7,9)]
    gene
  }

  if(cores>1){
    gene_list<-mclapply(gene_names, FUN = get_ranges, mc.cores = cores)
  } else {
    gene_list<-lapply(gene_names, get_ranges)
  }
  names(gene_list)<-gene_names
  #get sequence
  get_seq<-function(gene_name){
    gene<-gene_list[[gene_name]]
    gene_gr<-makeGRangesFromDataFrame(gene)
    seq<-getSeq(x = genome, names=gene_gr)
    names(seq)<-gene$class
    seq
  }
  if(cores>1){
    seq_list<-mclapply(gene_names, FUN = get_seq, mc.cores = cores)
  } else {
    seq_list<-lapply(gene_names, get_seq)
  }
  names(seq_list)<-gene_names
  #make genedf
  make_genedf<-function(gene_name){
    gene<-seq_list[[gene_name]]
    five_w<-width(gene[names(gene)=="five_UTR"])
    three_w<-width(gene[names(gene)=="three_UTR"])
    cds_w<-sum(width(gene[names(gene)=="cds"]))
    gene_df<-data.frame(gene=gene_name, start=five_w+1, end=sum(cds_w+five_w), ncodons=cds_w/3, nnuc=cds_w)
    gene_df
  }
  genedf<-do.call("rbind", lapply(gene_names, make_genedf))
  fasta<-lapply(seq_list, unlist)
  fasta<-DNAStringSet(fasta)
  support<-list(fasta=fasta, genedf=genedf)
  if(write){
    writeXStringSet(fasta, paste(fasta_name, "fasta", sep="."), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
    write.table(genedf, paste(genedf_name, "txt", sep="."), quote = F, sep="\t", col.names = T, row.names = F)
  }
  invisible(support)
}


















