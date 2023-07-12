removeInv <- function(alignment,output){
  #Removes invariant sites from a FASTA alignment
  #alignment = character string with the name of the file containing the original alignment
  #output = character string with the name of the output file; e.g., "./output.fas"
  seq.tmp <- readLines(alignment)
  dna.tmp <- seq.tmp[seq(2,length(seq.tmp),by=2)]
  dna.tmp <- sapply(dna.tmp,strsplit,split="")
  df.tmp <- matrix("A",ncol=length(dna.tmp[[1]]),nrow=length(dna.tmp))
  for(i in 1:nrow(df.tmp)){
    df.tmp[i,] <- dna.tmp[[i]]
  }
  fun.tmp <- function(x){
    x.tmp <- x[-which(x=="N" | x=="n")]
    l.tmp <- length(unique(x.tmp))
    if(l.tmp>1){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  inv.tmp <- apply(df.tmp,2,fun.tmp)
  dfn.tmp <- df.tmp[,-which(inv.tmp)]
  dnan.tmp <- apply(dfn.tmp,1,paste,collapse="")
  seq.tmp[seq(2,length(seq.tmp),by=2)] <- dnan.tmp
  writeLines(seq.tmp,output)
}