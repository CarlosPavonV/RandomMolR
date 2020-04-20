TerGaps2N <- function(alignment,outfile,missing.char,base.case){
  #Takes a FASTA alignment and replaces terminal gaps with N's or ?'s
  #alignment = character string, name of file containing FASTA alignment
  #outfile = character string, specifying name for output FASTA with replaced gaps
  #missing.char = character string specifying what character to use instead of terminal gaps, e.g. "?" or "n"
  #base.case = character string specifying what case to use for bases, either "lower" or "upper"
  require(seqinr)
  require(tidyr)
  alin <- read.fasta(file=alignment,as.string=T)
  len.tmp <- nchar(alin[1][1])
  for (i in 1:length(alin)){
    alin[i][1] <- gsub("^\\P{L}*","",alin[i][1],perl=T)
    while(nchar(alin[i][1])<len.tmp){
      alin[i][1] <- paste(missing.char,alin[i][1],sep="")
    }
    alin[i][1] <- gsub("\\P{L}*$","",alin[i][1],perl=T)
    while(nchar(alin[i][1])<len.tmp){
      alin[i][1] <- paste(alin[i][1],missing.char,sep="")
    }
    if(base.case=="lower"){
      for(i in 1:length(alin)){
        alin[i][1] <- tolower(alin[i][1])
      }
    } else{
      if(base.case=="upper"){
        for(i in 1:length(alin)){
          alin[i][1] <- toupper(alin[i][1])
        }
      } else{
        print("Cannot understand what case is desired for bases, please check arguments")
      }
    }
  }
  write.fasta(sequences=alin,names=names(alin),file.out=outfile)
}
