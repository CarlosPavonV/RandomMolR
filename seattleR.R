seattleR <- function(alignment,part.file,outfile){
  #If you use this script, please cite it as: Pavón-Vázquez, E.A., & C.J. Pavón-Vázquez. 2019. SeATTLE: Sequence Alignment Transformation into a Table for Later Edition. Available at: https://github.com/CarlosPavonV/seattle
  #Takes a concatenated FASTA alignment and get a data frame in which columns correspond to loci
  #Output can be used, among other things, to extract SNPs using the phrynomics package
  #alignment = character string, name of file containing FASTA alignment
  #part.file = character string, name of file containing partitions
  #outfile = optional. Character string specifying name for output file, which is a tab delimited text file
  require(seqinr)
  require(tidyr)
  alin <- read.fasta(file=alignment)
  part.df <- read.csv(part.file,sep="=",header=FALSE)
  part.df <- separate(part.df,col=2,into=c("beg","end"),sep="-",convert=T,remove=T)
  alin.tab <- data.frame(matrix(nrow=length(alin),ncol=nrow(part.df)))
  for (i in 1:nrow(part.df)){
    print(paste("Processing locus ",i," of ",nrow(part.df),sep=""))
    for (f in 1:length(alin)){
      alin.ind <- alin[[f]][part.df$beg[[i]]:part.df$end[[i]]]
      seq.temp <- paste0(alin.ind,collapse="")
      seq.temp <- toupper(seq.temp)
      alin.tab[f,i] <- seq.temp
    }
  }
  rownames(alin.tab) <- names(alin)
  if (missing(outfile)){
    print("Outfile not created")
  } else{
    write.table(alin.tab,quote=F,file=outfile,row.names=T,col.names=F,sep="\t")
  }
  print("If you use this script, please cite it as: Pavón-Vázquez, E.A., & C.J. Pavón-Vázquez. 2019. SeATTLE: Sequence Alignment Transformation into a Table for Later Edition. Available at: https://github.com/CarlosPavonV/seattle")
  return(alin.tab)
}
