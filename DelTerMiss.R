DelTerMiss <- function(alignment,base.case,outfile){
  #Takes a FASTA alignment and deletes terminal positions that are represented solely by missing data
  #alignment = character string, name of input file containing FASTA alignment
  #base.case = one of "upper", "lower", or "gap".
    #Use "upper" if the bases in your alignment are in upper case and terminal missing data is coded as "N"
    #Use "lower" if the bases in your alignment are in lower case and terminal missing data is coded as "n"
    #Use "gap" if the terminal missing data is coded as "-"
  #outfile = character string, specifying name for output FASTA
  require(purrr)
  require(seqinr)
  if (base.case=="upper"){
    miss.char="N"
  } else {
    if (base.case=="lower"){
      miss.char="n" 
    } else {
      if (base.case=="gap"){
        miss.char="-"
      } else {
        print("Check base.case argument")
      }
    }
  } 
  alin.tmp <- read.fasta(file=alignment,forceDNAtolower = F)
  pos.del <- c()
  num.sites <- length(alin.tmp[[1]])
  for (j in 1:num.sites){
    pos.tmp <- gsub(miss.char,NA,unlist(map(alin.tmp,j)))
    if(all(is.na(pos.tmp))){
      pos.del <- c(pos.del,j)
    }
  }
  if(length(pos.del)>0){
    lead.del <- pos.del[-((which(diff(c(0,pos.del))>1)[1]):length(pos.del))]
    drag.del <- pos.del[-(1:length(which(diff(rev(c(pos.del,(num.sites+1))))<(-1))[1]:length(pos.del)))]
  }
  if(length(lead.del)>0){
    for (j in 1:length(alin.tmp)){
      alin.tmp[[j]] <- alin.tmp[[j]][-lead.del]
    }
  }
  if(length(drag.del)>0){
    for (j in 1:length(alin.tmp)){
      alin.tmp[[j]] <- alin.tmp[[j]][-drag.del]
    }
  }
  write.fasta(alin.tmp,names=names(alin.tmp),nbchar=num.sites,file.out=outfile)
}
