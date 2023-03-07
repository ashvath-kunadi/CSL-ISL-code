evalsim <- function(dat, gdat){
  res <- 0
  index <- 1
  for(i in 6:7){
    thrfdat <- gdat[[i]]
      for(j in 1:dim(thrfdat)[2]){
        res[index] <- sum(dat-thrfdat[,j], na.rm = T)
        index <- index+1
      }
  }
  return(res)
}
