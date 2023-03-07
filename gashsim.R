gashsim <- function(gdat){
  Pg <- gdat[[1]]
  Rbar <- gdat[[2]]
  Ebar <- gdat[[3]]
  recsim <- data.frame()
  for(sim in 1:1000){
    S <- runif(1, min = 0.1, max = 3.5)
    p <- runif(1, min = 0.01, max = 0.8)
    Pgdash <- (-gdat[[2]]*S/gdat[[3]])*log(1-(gdat[[3]]/(gdat[[2]]*(1-p))))
    Pgdash <- par[[3]]
    Int <- 0
    for(i in 1:length(Pg)){
      if((Pg[i])<Pgdash){
        Int[i] <- (1-p)*Pg[i]
      }else{
        Int[i] <- (1-p)*Pgdash + Ebar/Rbar*(Pg[i]-Pgdash)
      }
      recsim <- rbind(recsim, cbind(i, S, p, Pgdash, Int))
    }
  }
  return(recsim)
}
