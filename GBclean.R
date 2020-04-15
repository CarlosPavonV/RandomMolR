GBclean <- function(infile,keep.strings,outfile){
  #Takes a FASTA file downloaded from GB and returns an alignment (as a list, as is used by seqinr) with shortened sequence names to include only relevant information
  #infile = character string containing name of file with FASTA alignment
  #keep.strings = numeric vector containing which elements of sequence name (separated by spaces) to keep for new names
  #outfile = character string containing name for FASTA output file, if argument missing output will not be written to file
  require(seqinr)
  alin <- read.fasta(infile,whole.header = T)
  names1 <- names(alin)
  names2 <- gsub("[.]1","",names1)
  names3 <- c()
  for (i in 1:length(names2)){
    split.names <- strsplit(names2,split=" ")
    names3[[i]] <- paste0(split.names[[i]][keep.strings],collapse="_")
  }
  alin2 <- alin
  names(alin2) <- names3
  for(i in 1:length(alin2)){
    attr(alin2[[i]],"name") <- names(alin2)[i]
  }
  if (missing(outfile)){
    print("Output has not been written to file, returning alignment with sequence names changed")
  } else{
    write.fasta(alin2,names=names(alin2),file.out=outfile)
  }
  return(alin2)
}
