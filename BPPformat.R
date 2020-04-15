BPPformat <- function(alignment,part.file,outfile){
  #Takes a concatenated FASTA alignment and a partition file in raxml format and gives a sequence file for BPP
  #alignment = character string, name of file containing FASTA alignment
  #part.file = character string, name of file containing partitions
  #outfile = character string, name for sequence file in BPP format
  require(seqinr)
  require(tidyr)
  require(chopper)
  alin <- read.fasta(file=alignment)
  names1 <- names(alin)
  names2 <- c()
  for (i in 1:length(names1)){
    names2[[i]] <- paste(names1[[i]],"^",i,sep="")
  }
  part.df <- read.csv(part.file,sep="=",header=FALSE)
  part.df <- separate(part.df,col=2,into=c("beg","end"),sep="-",convert=T,remove=T)
  for (i in 1:nrow(part.df)){
    print(paste("Changing sequence names in locus ",i," of ",nrow(part.df),". Step 1 of 3",sep=""))
    alin.list <- list()
    for (f in 1:length(alin)){
      alin.ind <- alin[[f]][part.df$beg[[i]]:part.df$end[[i]]]
      alin.list[[f]] <- alin.ind
      names(alin.list)[f] <- names2[f]
    }
    write.fasta(alin.list,names=names(alin.list),nbchar=100,file.out=paste(i,"temporary_delete.fas",sep="_"))
  }
  filesinf1 <- list.files()
  filesdel1 <- filesinf1[grepl("temporary_delete",filesinf1)==TRUE]
  for(i in 1:length(filesdel1)){
    print(paste("Transforming locus ",i," of ",length(filesdel1),". Step 2 of 3",sep=""))
    fas2phy(filesdel1[[i]],format="sequential",overwrite=T)
    file.remove(filesdel1[[i]])
  }
  filesinf2 <- list.files()
  filesdel2 <- filesinf2[grepl("temporary_delete.phy",filesinf2)==TRUE]
  filesdel2 <- filesdel2[order(nchar(filesdel2),filesdel2)]
  out.vec <- c()
  for(i in 1:length(filesdel2)){
    print(paste("Processing locus ",i," of ",length(filesdel2),". Step 3 of 3",sep=""))
    linesphy <- readLines(filesdel2[[i]])
    out.vec <- c(out.vec,linesphy)
    file.remove(filesdel2[[i]])
  }
  writeLines(out.vec,outfile)
}
