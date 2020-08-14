CreateIndLoc <- function(alignment,parts,prefix){
  #Takes a concatenated FASTA alignment and a partitions file in raxml style, gives files with individual locus alignments
  #alignment = character string, name of file with concatenated FASTA alignment
  #parts = character string, name of file with raxml style partitions
  #prefix = character string, appended to file names with individual alignments, e.g. "Locus"
  require(seqinr)
  require(tidyr)
  alin <- read.fasta(file=alignment,forceDNAtolower=F)
  part.df <- read.csv(parts,sep="=",header=FALSE)
  part.df <- separate(part.df,col=2,into=c("beg","end"),sep="-",convert=T,remove=T)
  for (i in 1:nrow(part.df)){
    alin.list <- list()
    for (f in 1:length(alin)){
      alin.ind <- alin[[f]][part.df$beg[[i]]:part.df$end[[i]]]
      alin.list[[f]] <- alin.ind
      names(alin.list)[f] <- names(alin)[f]
    }
    print(paste("Writing locus ",i," of ",nrow(part.df),sep=""))
    write.fasta(alin.list,names=names(alin.list),nbchar=100,file.out=paste(paste(prefix,i,sep="_"),".fas",sep=""))
  }
}
