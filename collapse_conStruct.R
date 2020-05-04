collapse_conStruct <- function(freqs,coords,miles,threshold,logdist,filterlocs,filtergen){
  #Collapses localities that are geographically close to each other based on a distance threshold
  #freqs= a matrix of frequencies in conStruct format
  #coords = matrix of coordinates (rownames=populations or individuals,column1=longitude,colum2=latitude)
  #miles = logical, if TRUE the output units are miles, if FALSE kilometers
  #threshold = distance threshold for collapsing localities, in same units as decided by miles argument
  #logdist = should the geographic distance matrix be log-transformed?
  #filterlocs = optional, localities with a proportion of missing loci higher than the threshold (after collapsing) will be deleted
  #filtergen = optional, loci with a proportion of missing localities higher than the threshold (after collapsing) will be deleted
  #Value:a list containing:
  #freqs.fin = the frequency matrix after collapsing nearby localities, optionally filtered by proportion of missing data
  #collapsed.locs = data frame indicating the original locality name, a tag indicating what name was given to the locality after collapsing, and the mean longitude and latitiude for each tag
  #coords.coll.fin = the collapsed coordinate matrix
  #geo.dist.fin = the collapsed geographic distance matrix, optionally log-transformed
  #Note: ploidy is assumed to be the same for all collapsed individuals!
  require(fields)
  geo.distan <- rdist.earth(coords,miles=miles)
  dimnames(geo.distan) <- list(rownames(coords),rownames(coords))
  names.tax <- t(combn(colnames(geo.distan), 2))
  geo.distan.table <- data.frame(names.tax, dist=geo.distan[names.tax])
  small.dist <- geo.distan.table[which(geo.distan.table$dist<threshold),]
  row.names(small.dist) <- 1:nrow(small.dist)
  small.dist$X1 <- as.character(small.dist$X1)
  small.dist$X2 <- as.character(small.dist$X2)
  coords.coll <- as.data.frame(coords)
  coords.coll$Group <- rownames(coords.coll)
  while(any(coords.coll$Group %in% small.dist[,2])){
    for (i in 1:nrow(coords.coll)){
      if(coords.coll$Group[i] %in% small.dist[,2]){
        coords.coll$Group[i] <- small.dist[which(small.dist[,2]==coords.coll$Group[i]),1][1]
      }
    }
  }
  coords.coll.fin <- aggregate(coords.coll[,1:2], list(coords.coll$Group), mean)
  rownames(coords.coll.fin) <- coords.coll.fin$Group.1
  coords.coll.fin <- coords.coll.fin[,2:3]
  collapsed.locs <- data.frame(Sample=rownames(coords.coll),TagCollapsed=coords.coll$Group,meanLon=coords.coll$Lon,meanLat=coords.coll$Lat)
  if(logdist==TRUE){
    geo.dist.fin <- log(rdist.earth(coords.coll.fin,miles=miles))
  } else{
    geo.dist.fin <- rdist.earth(coords.coll.fin,miles=miles)
  }
  diag(geo.dist.fin) <- 0
  geo.dist.fin[which(!is.finite(geo.dist.fin))] <- 0
  geo.dist.fin[which(geo.dist.fin<0)] <- 0
  freqs.mod <- merge(freqs,coords.coll,by="row.names")
  rownames(freqs.mod) <- freqs.mod$Row.names
  freqs.mod <- freqs.mod[,c(2:(ncol(freqs)+1),ncol(freqs.mod))]
  colnames(freqs.mod)[ncol(freqs.mod)] <- "Group"
  freqs.fin <- as.matrix(aggregate(freqs.mod[,1:(length(freqs.mod)-1)], list(freqs.mod$Group), mean, na.rm=T))
  rownames(freqs.fin) <- freqs.fin[,1]
  freqs.fin <- freqs.fin[,-1]
  if(!missing(filterlocs)&!missing(filtergen)){
    freqs.fin <- freqs.fin[which(rowMeans(is.na(freqs.fin)) <= filterlocs), which(colMeans(is.na(freqs.fin)) <= filtergen)]
  } else{
    if(!missing(filterlocs)&missing(filtergen)){
      freqs.fin <- freqs.fin[which(rowMeans(is.na(freqs.fin)) <= filterlocs),]
    } else{
      if(missing(filterlocs)&!missing(filtergen)){
        freqs.fin <- freqs.fin[,which(colMeans(is.na(freqs.fin)) <= filtergen)]
      } else{
        freqs.fin <- freqs.fin
      }
    }
  }
  freqs.fin.locs <- rownames(freqs.fin)
  freqs.fin.gen <- colnames(freqs.fin)
  freqs.fin <- sapply(freqs.fin,as.numeric)
  freqs.fin <- matrix(freqs.fin,ncol=length(freqs.fin.gen),nrow=length(freqs.fin.locs))
  dimnames(freqs.fin) <- list(freqs.fin.locs,freqs.fin.gen)
  list.out <- list()
  list.out$freqs.fin <- freqs.fin
  list.out$collapsed.locs <- collapsed.locs
  list.out$coords.coll.fin <- as.matrix(coords.coll.fin)
  list.out$geo.dist.fin <- geo.dist.fin
  return(list.out)
}
