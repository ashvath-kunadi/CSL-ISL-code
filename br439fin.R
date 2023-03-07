br439prf <- function(n, canh, pdat, coln = 2){
  library(ggplot2)
  library(gtools)
  pdatab <- pdat[pdat$H > canh,]
  pdatab <- pdatab[order(pdatab$H, pdatab$time),]
  senn <- dim(pdatab)[1]/n
  if(senn > 2){
    senpr <- combinations(senn, 2, v=1:senn, set=TRUE, repeats.allowed=FALSE)
  }else{
    senpr <- matrix(data = c(1,2), nrow=1)
  }
  dscdzdf <- data.frame()
  for(j in 1:(dim(senpr)[1])){
    arns <- 0
    sub1 <- pdatab[(1+((senpr[j,1]-1)*n)):(senpr[j,1]*n),]
    sub2 <- pdatab[(1+((senpr[j,2]-1)*n)):(senpr[j,2]*n),]
    for(i in 1:n){
      arns[i] <- (sub1[i,coln]-sub2[i,coln])/(sub1$H[i]-sub2$H[i])
    }
    ans <- data.frame(dscdz = arns, sen1 = sub1$H[i], sen2 = sub2$H[i], time = sub1$time)
    ans$zeval <- round(min(sub1$H[i], sub2$H[i])+(abs(sub1$H[i]-sub2$H[i])/2),1)
    dscdzdf <- rbind(dscdzdf, ans)
  }
  return(dscdzdf)
}
