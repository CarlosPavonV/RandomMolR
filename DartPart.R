DartPart <- function(globject,format.part,outfile,align.prefix,resol.method){
  #Takes a genlight object and creates a file with the partitions for each locus specified in either raxml or nexus format
  #Optionally writes individual files with the alignment of each individual locus
  #globject = input genlight object
  #format.part = either "raxml" or "nexus"
  #outfile = name for the output file, e.g. "parts_nexus.txt"
  #align.prefix = optional character string. If given, individual FASTA alignments are created for each locus and the prefix appended to each file name
  #resol.method = if creating individual alignments, the method used to resolve heterozygous positions. Must be 1 or 2 (see help for gl2fasta function of dartR package)
  require(adegenet)
  require(ade4)
  or.options <- options()
  options(scipen = 999)
  trim.seqs <- globject$other$loc.metrics$TrimmedSequence
  loc.l <- nchar(trim.seqs,"char")
  par.t <- data.frame(loc.l)
  par.t$beg[[1]] <- 1
  par.t$end[[1]] <- loc.l[[1]]
  par.t$part[[1]] <- paste("DNA, Locus",1," = ",par.t$beg[[1]],"-",par.t$end[[1]],sep="")
  for (i in 2:length(loc.l)){
    par.t$beg[[i]] <- par.t$end[[i-1]]+1
    par.t$end[[i]] <- par.t$beg[[i]]-1+par.t$loc.l[[i]]
    par.t$part[[i]] <- paste("DNA, Locus",i," = ",par.t$beg[[i]],"-",par.t$end[[i]],sep="")
  }
  if (format.part == "raxml"){
    write(par.t$part,file=outfile)
  } else{
    if (format.part == "nexus"){
      temp.mat1 <- as.data.frame(matrix(1:8,nrow=2,ncol=4))
      temp.mat2 <- as.data.frame(matrix(1:4,nrow=1,ncol=4))
      colnames(temp.mat1) <- colnames(par.t)
      colnames(temp.mat2) <- colnames(par.t)
      par.t2 <- rbind(temp.mat1,par.t,temp.mat2)
      par.t2$part[[1]] <- "#nexus"
      par.t2$part[[2]] <- "begin sets;"
      par.t2$part[[length(par.t2$part)]] <- "end;"
      par.t2$part <- gsub("DNA,","charset",par.t2$part)
      for (i in 3:(length(par.t2$part)-1)){
        par.t2$part[[i]] <- paste(par.t2$part[[i]],";",sep="")
      }
      write(par.t2$part,file=outfile)
    } else{
      print("Cannot understand what output format is desired, please check arguments")
    }
  }
  options(or.options)
  if (missing(align.prefix)){
    print("Individual locus alignments not created")
  } else{
    require(dartR)
    require(seqinr)
    gl2fasta(globject,method=resol.method,outpath=getwd(),outfile="temp.fas")
    glalin <- read.fasta(file="temp.fas")
    for (i in 1:length(loc.l)){
      alin.list <- list()
      for (f in 1:length(glalin)){
        alin.ind <- glalin[[f]][par.t$beg[[i]]:par.t$end[[i]]]
        alin.list[[f]] <- alin.ind
        names(alin.list)[f] <- names(glalin)[f]
      }
      write.fasta(alin.list,names=names(alin.list),nbchar=100,file.out=paste(paste(align.prefix,i,sep="_"),".fas",sep=""))
    }
    file.remove("temp.fas")
  }
}
