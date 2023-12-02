KeepFromGraphs <- function(input,threshold,optimize.me,return.groups){
  #Often you want to keep characters/loci that are independent from each other, deleting tightly correlated characters/loci
  #This function is similar to KeepFromCorrs but more stringent, it uses graphs to delete groups of correlated variables.
  #Let's say trait1 is correlated to trait2, trait2 is correlated to trait3, and trait1 and trait3 are uncorrelated. Because 1 is connected to 3 through their mutual correlation with 2, the three loci will be considered as part of a group of correlated variables
  #Among groups, you would probably like to keep the characters/loci that is more informative, or has fewer missing data
  #This function returns a vector with the names of the characters to keep from a set of correlated variables and optionally a list with the groups of correlated variables
  #input = a dist object corresponding to a pairwise correlation matrix,
  #        or a pairwise correlation (or p-values of correlation) matrix,
  #        or a data frame in which the first and second columns contain a unique name for each character/loci in a pairwise comparison and the third column a measure of correlation (e.g., p-value, R2)
  #        the fourth and fifth column in the data frame can contain a numeric indicator of quality for the first and second column, respectively. The criterion can be used optionally when deciding which element to keep (see under optimize.me)
  #        if no criterion is going to be used, please still include the fourth and fifth columns as NA's
  #threshold = numeric, character/loci comparisons for which a measure of correlation is BELOW the threshold will be considered to be tightly correlated
  #            if the measure of correlation is a p-value, you could indicate 0.05 for example
  #            if the measure of correlation in your dataset is R2 or a similar coefficient, indicate here {1 - the value above which you consider correlation to be tight}
  #            e.g., if you consider pairs with R2 above 0.8 to be tightly correlated, input 0.2 here and NOT 0.8
  #optimize.me = TRUE or FALSE, if TRUE a criterion (e.g., proportion of missing data, proportion of parsimony informative sites, variance) is used to decide which character/loci to keep among each group of correlated variables
  #               keep in mind that a larger value of the criterion is considered to be a positive quality, for example for number of parsimony informative sites the larger the better
  #               if the criterion corresponds to a negative quality (e.g., proportion of missing data), please replace the criterion by {1 - criterion} or {1 / criterion} in your input (e.g. change a proportion of missing data from 0.9 to 0.1)
  #               only applicable when input is a data frame. Otherwise the decision of which character/loci to keep is arbitrary
  #return.groups = TRUE or FALSE. If TRUE, the function will return a list with two elements: a vector with the characters/loci to be kept and a list with the groups of correlated variables.
  #                If FALSE, it will return just the vector of characters/loci to be kept.
  require(igraph)
  '%notin%' <- Negate('%in%')
  if(class(input)[1]=="data.frame"){
    input[,1:2] <- sapply(input[,1:2],as.character)
    input[,3:5] <- sapply(input[,3:5],as.numeric)
    input.all.tr <- unique(c(input[,1],input[,2]))
    sigld.input <- input[which(input[,3]<threshold),]
    corr.unique <- unique(c(sigld.input[,1],sigld.input[,2]))
    conect.list <- lapply(1:length(corr.unique),function(i){c(corr.unique[i],unique(c(sigld.input[which(sigld.input[,1]==corr.unique[i]),2],sigld.input[which(sigld.input[,2]==corr.unique[i]),1])))})
    mat.con <-  sapply(conect.list,function(x) sapply(conect.list,function(y) length(intersect(x,y))>0))
    dimnames(mat.con) <- list(corr.unique,corr.unique)
    groups.con <-  groups(components(graph_from_adjacency_matrix(mat.con)))
    if(optimize.me==TRUE){
      score.tr <- data.frame(trait=c(input[,1],input[,2]),score=c(input[,4],input[,5]))
      score.tr <- unique(score.tr)
      score.list <- lapply(1:length(groups.con),function(i){score.tr[which(score.tr$trait %in% groups.con[[i]]),2]})
      pos.max.sc <- lapply(1:length(groups.con),function(i){which(score.list[[i]]==max(score.list[[i]]))[1]})
      keep1 <- sapply(1:length(groups.con),function(i){groups.con[[i]][pos.max.sc[[i]]]})
    } else{
      keep1 <- sapply(1:length(groups.con),function(i){groups.con[[i]][1]})
    }
    keep2 <- input.all.tr[which(input.all.tr %notin% corr.unique)]
    keep.fin <- c(keep1,keep2)
    keep.fin <- sort(keep.fin)
  } else{
    if(class(input)[1]=="dist"|class(input)[1]=="matrix"){
      input0 <- as.matrix(input)
      input <- t(combn(colnames(input0),2))
      input <- data.frame(input, dist=input0[input])
      input[,1:2] <- sapply(input[,1:2],as.character)
      input[,3] <- sapply(input[,3],as.numeric)
      input.all.tr <- unique(c(input[,1],input[,2]))
      sigld.input <- input[which(input[,3]<threshold),]
      corr.unique <- unique(c(sigld.input[,1],sigld.input[,2]))
      conect.list <- lapply(1:length(corr.unique),function(i){c(corr.unique[i],unique(c(sigld.input[which(sigld.input[,1]==corr.unique[i]),2],sigld.input[which(sigld.input[,2]==corr.unique[i]),1])))})
      mat.con <-  sapply(conect.list,function(x) sapply(conect.list,function(y) length(intersect(x,y))>0))
      dimnames(mat.con) <- list(corr.unique,corr.unique)
      groups.con <-  groups(components(graph_from_adjacency_matrix(mat.con)))
      keep1 <- sapply(1:length(groups.con),function(i){groups.con[[i]][1]})
      keep2 <- input.all.tr[which(input.all.tr %notin% corr.unique)]
      keep.fin <- c(keep1,keep2)
      keep.fin <- sort(keep.fin)
    } else{
      cat("Couldn't understand input, please check.\n")
    }
  }
  if(return.groups==TRUE){
    ret.list <- list()
    ret.list$kept.traits <- keep.fin
    ret.list$trait.groups <- groups.con
    return(ret.list)
  } else{
    return(keep.fin)
  }
}
