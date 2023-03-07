numderd <- function(dudz, Et0dat, n, time, ntr){
  d0 <- 0
  comb <- dim(dudz)[1]/n
  tim <- rep(time, comb)
  for(i in 1:length(tim)){
    if(tim[i] %in% ntr){
      d0[i] <- dudz$zeval[i]-(rep(Et0dat$ustar,comb)[i]/dudz$dscdz[i]/0.4)
      if(!is.na(d0[i])){
        if(d0[i]<0){
        d0[i] <- NA
        }
      }
    }
    else{
      d0[i] <- NA
    }
  }
  ans <- data.frame(d0, tim, cmb = rep(1:comb, each = n))
  return(ans)
}
