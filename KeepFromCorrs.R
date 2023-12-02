KeepFromCorrs <- function(input,threshold,optimize.me,verbose){
  #Often you want to keep characters/loci that are independent from each other, deleting tightly correlated characters/loci
  #Often you will want to keep as many characters as possible. Let's say trait1 and trait2 are correlated, trait 2 and trait3 are correlated, but trait1 and trait 3 are not correlated, if you delete trait2 you can keep the other two
  #Then maybe you would like to keep those characters/loci that are more informative, or have fewer missing data
  #This can be hard to do when you have many characters/loci
  #This function makes your life easier, it returns a vector with the names of the characters to keep from a set of correlated variables
  #The function tries to keep as many independent traits/characters as possible, and then it can optionally choose among pairs of variables using a criterion
  #input = a dist object corresponding to a pairwise correlation matrix,
  #        or a pairwise correlation (or p-values of correlation) matrix,
  #        or a data frame in which the first and second columns contain a unique name for each character/loci in a pairwise comparison and the third column a measure of correlation (e.g., p-value, R2)
  #        the fourth and fifth column in the data frame can contain a numeric indicator of quality for the first and second column, respectively. The criterion can be used optionally when deciding which element to keep (see under optimize.me)
  #        if no criterion is going to be used, please still include the fourth and fifth columns as NA's
  #threshold = numeric, character/loci comparisons for which a measure of correlation is BELOW the threshold will be considered to be tightly correlated
  #            if the measure of correlation is a p-value, you could indicate 0.05 for example
  #            if the measure of correlation in your dataset is R2 or a similar coefficient, indicate here {1 - the value above which you consider correlation to be tight}
  #            e.g., if you consider pairs with R2 above 0.8 to be tightly correlated, input 0.2 here and NOT 0.8
  #optimize.me = TRUE or FALSE, if TRUE a criterion (e.g., proportion of missing data, proportion of parsimony informative sites, variance) is used to decide which character/loci to keep among each tightly correlated pair
  #               keep in mind that a larger value of the criterion is considered to be a positive quality, for example for number of parsimony informative sites the larger the better
  #               if the criterion corresponds to a negative quality (e.g., proportion of missing data), please replace the criterion by {1 - criterion} or {1 / criterion} in your input (e.g. change a proportion of missing data from 0.9 to 0.1)
  #               only applicable when input is a data frame. Otherwise the decision of which character/loci to keep is arbitrary
  #verbose = TRUE or FALSE, self-explanatory
  '%notin%' <- Negate('%in%')
  if(class(input)[1]=="data.frame"){
    p.pairs <- input
    p.pairs[,1:2] <- sapply(p.pairs[,1:2],as.character)
    p.pairs[,3:5] <- sapply(p.pairs[,3:5],as.numeric)
    df.use <- p.pairs
    p.pairs <- p.pairs[p.pairs[,3]<threshold,]
    corr.all <- unique(c(p.pairs[,1],p.pairs[,2]))
    p.pairs$longcomp <- 2
    p.pairs$longcomp <- as.numeric(sapply(1:nrow(p.pairs),function(i){length(p.pairs[i,1:2][!is.na(p.pairs[i,1:2])])}))
    char.2membs <- as.character(unlist(p.pairs[,1:2]))
    counts <- as.data.frame(table(char.2membs))
    changed2NA <- as.character(counts[which(counts$Freq==max(counts$Freq)),1][1])
    p.pair.repl.max <- as.data.frame(sapply(p.pairs,gsub,pattern=as.character(counts[which(counts$Freq==max(counts$Freq)),1][1]),replacement=NA),stringsAsFactors = F)
    outway <- p.pair.repl.max[which(is.na(p.pair.repl.max[,1])|is.na(p.pair.repl.max[,2])),]
    outway$longcomp <- 1
    p.pairs <- p.pair.repl.max[which(!is.na(p.pair.repl.max[,1])&!is.na(p.pair.repl.max[,2])),]
    char.2membs <- as.character(unlist(p.pairs[,1:2]))
    counts <- as.data.frame(table(char.2membs))
    keep.track=0
    while(any(counts$Freq>1)){
      changed2NA <- c(changed2NA,as.character(counts[which(counts$Freq==max(counts$Freq)),1][1]))
      p.pair.repl.max <- as.data.frame(sapply(p.pairs,gsub,pattern=as.character(counts[which(counts$Freq==max(counts$Freq)),1][1]),replacement=NA),stringsAsFactors = F)
      outway2 <- p.pair.repl.max[which(is.na(p.pair.repl.max[,1])|is.na(p.pair.repl.max[,2])),]
      outway <- rbind(outway,outway2)
      outway$longcomp <- 1
      p.pairs <- p.pair.repl.max[which(!is.na(p.pair.repl.max[,1])&!is.na(p.pair.repl.max[,2])),]
      char.2membs <- as.character(unlist(p.pairs[,1:2]))
      counts <- as.data.frame(table(char.2membs))
      if(verbose==TRUE){
        keep.track=keep.track+1
        bad.counts <- length(counts$Freq[which(counts$Freq>1)])
        print(paste("Round ",keep.track,". ","Bad counts = ",bad.counts,": This number must be 0 in order to proceed to next step.",sep="")) 
      }
    }
    p.pairs <- rbind(outway,p.pairs)
    p.pairs[,3:6] <- sapply(p.pairs[,3:6],as.numeric)
    if(optimize.me==TRUE){
      keep1 <- ifelse(p.pairs[,4]>=p.pairs[,5],p.pairs[,1],NA)
      keep2 <- ifelse(p.pairs[,4]<p.pairs[,5],p.pairs[,2],NA)
      pos.change <- which(p.pairs$longcomp>1)
      delete.me1 <- p.pairs[pos.change,1][which(p.pairs[pos.change,1] %notin% keep1[pos.change])]
      delete.me2 <- p.pairs[pos.change,2][which(p.pairs[pos.change,2] %notin% keep2[pos.change])]
      delete.me <- unique(c(delete.me1,delete.me2))
      p.pairs[pos.change,1] <- keep1[pos.change]
      p.pairs[pos.change,2] <- keep2[pos.change]
      keep.me <- unique(c(p.pairs[,1],p.pairs[,2]))
      keep.me <- keep.me[which(!is.na(keep.me))]
      unique.input <- unique(c(df.use[,1],df.use[,2]))
      unique.input <- unique.input[which(unique.input %notin% corr.all)]
      keep.me <- c(keep.me,unique.input)
      keep.me <- keep.me[which(keep.me %notin% delete.me)]
      keep.me <- keep.me[which(keep.me %notin% changed2NA)]
      keep.me <- sort(keep.me)
    } else{
      pos.change <- which(p.pairs$longcomp>1)
      delete.me <- p.pairs[pos.change,2]
      keep.me <- unique(c(p.pairs[,1],p.pairs[,2]))
      keep.me <- keep.me[which(!is.na(keep.me))]
      unique.input <- unique(c(df.use[,1],df.use[,2]))
      unique.input <- unique.input[which(unique.input %notin% corr.all)]
      keep.me <- c(keep.me,unique.input)
      keep.me <- keep.me[which(keep.me %notin% delete.me)]
      keep.me <- keep.me[which(keep.me %notin% changed2NA)]
      keep.me <- sort(keep.me)
      }
  } else{
    if(class(input)[1]=="dist"|class(input)[1]=="matrix"){
      p.pairs0 <- as.matrix(input)
      p.pairs <- t(combn(colnames(p.pairs0),2))
      p.pairs <- data.frame(p.pairs, dist=p.pairs0[p.pairs])
      p.pairs[,1:2] <- sapply(p.pairs[,1:2],as.character)
      p.pairs[,3] <- sapply(p.pairs[,3],as.numeric)
      df.use <- p.pairs
      p.pairs <- p.pairs[p.pairs[,3]<threshold,]
      corr.all <- unique(c(p.pairs[,1],p.pairs[,2]))
      p.pairs$longcomp <- 2
      p.pairs$longcomp <- as.numeric(sapply(1:nrow(p.pairs),function(i){length(p.pairs[i,1:2][!is.na(p.pairs[i,1:2])])}))
      char.2membs <- as.character(unlist(p.pairs[,1:2]))
      counts <- as.data.frame(table(char.2membs))
      changed2NA <- as.character(counts[which(counts$Freq==max(counts$Freq)),1][1])
      p.pair.repl.max <- as.data.frame(sapply(p.pairs,gsub,pattern=as.character(counts[which(counts$Freq==max(counts$Freq)),1][1]),replacement=NA),stringsAsFactors = F)
      outway <- p.pair.repl.max[which(is.na(p.pair.repl.max[,1])|is.na(p.pair.repl.max[,2])),]
      outway$longcomp <- 1
      p.pairs <- p.pair.repl.max[which(!is.na(p.pair.repl.max[,1])&!is.na(p.pair.repl.max[,2])),]
      char.2membs <- as.character(unlist(p.pairs[,1:2]))
      counts <- as.data.frame(table(char.2membs))
      keep.track=0
      while(any(counts$Freq>1)){
        changed2NA <- c(changed2NA,as.character(counts[which(counts$Freq==max(counts$Freq)),1][1]))
        p.pair.repl.max <- as.data.frame(sapply(p.pairs,gsub,pattern=as.character(counts[which(counts$Freq==max(counts$Freq)),1][1]),replacement=NA),stringsAsFactors = F)
        outway2 <- p.pair.repl.max[which(is.na(p.pair.repl.max[,1])|is.na(p.pair.repl.max[,2])),]
        outway <- rbind(outway,outway2)
        outway$longcomp <- 1
        p.pairs <- p.pair.repl.max[which(!is.na(p.pair.repl.max[,1])&!is.na(p.pair.repl.max[,2])),]
        char.2membs <- as.character(unlist(p.pairs[,1:2]))
        counts <- as.data.frame(table(char.2membs))
        if(verbose==TRUE){
          keep.track=keep.track+1
          bad.counts <- length(counts$Freq[which(counts$Freq>1)])
          print(paste("Round ",keep.track,". ","Bad counts = ",bad.counts,": This number must be 0 in order to proceed to next step.",sep="")) 
        }
      }
      p.pairs <- rbind(outway,p.pairs)
      pos.change <- which(p.pairs$longcomp>1)
      delete.me <- p.pairs[pos.change,2]
      p.pairs[pos.change,2] <- NA
      keep.me <- unique(c(p.pairs[,1],p.pairs[,2]))
      keep.me <- keep.me[which(!is.na(keep.me))]
      unique.input <- unique(c(df.use[,1],df.use[,2]))
      unique.input <- unique.input[which(unique.input %notin% corr.all)]
      keep.me <- c(keep.me,unique.input)
      keep.me <- keep.me[which(keep.me %notin% delete.me)]
      keep.me <- keep.me[which(keep.me %notin% changed2NA)]
      keep.me <- sort(keep.me)    } else {
      cat("Couldn't understand input, please check.\n")
    }
  }
  return(keep.me)
}
