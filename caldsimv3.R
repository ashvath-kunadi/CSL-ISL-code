caldsim <- function(usedat, gpfldat){
  library(deSolve)
  rfdat <- usedat[[1]]
  ET0dat <- usedat[[2]]
  L <- 0
  q <- 0 
  v <- 0
  rec <- 1
  nfun <- function(m,q){
    if(q > 1.0000001){
      pdcs <- 0
      for(x in 1:floor(q-0.0000001)){
        pdcs <- pdcs+(x-q)*m^(x)/factorial(x)
      }
      ans <- (q*(1-exp(-m))+exp(-m)*pdcs)
      return(ans)
    }else{
      ans <- q*(1-exp(-m))
      return(ans)
    }
  }
  dndm <- function(m,q){
    if(q > 1.0000001){
      pdcs <- 0
      for(x in 1:floor(q-0.0000001)){
        pdcs <- pdcs+((x-q)*(x-m)*m^(x-1)/factorial(x))
      }
      ans <- q*exp(-m)+exp(-m)*pdcs
      return(ans)
    }else{
      ans <- q*exp(-m)
      return(ans)
    }
  }
  recsim <- data.frame()
  for(sim in 1:100){
    L[sim] <- runif(1, min = 1, max = 500)*10^-3
    q[sim] <- runif(1, min = 0.25, max = 5)
    v[sim] <- 4/3*pi*abs(rnorm(1, mean = 1, sd = 0.25))^3
    a <- 1/(v[sim]*L[sim])
    for(k in 1:length(unique(rfdat$ID))){
      IDen <- unique(rfdat$ID)[k]
      subrf <- rfdat[(rfdat$ID == IDen),]
      subrf <- subrf[-c(1:59),2]
      subev <- ET0dat[(ET0dat$ID == IDen),]
      subev <- subev[-c(1:59),2]
      if(sum(is.na(subev)) == 0){
        qsim <- q[sim]
        caldhal <- function(t, y, parms){
            with(as.list(y), {
              dPdt <- subrf[t]
              dIdt <- subev[t]
              dmdt <- a*dPdt-a*dIdt/dndm(m, qsim)
              list(dmdt)
            })
        }
        yini <- c(m=0)
        times <- seq(from = 1, to = length(subrf), by = 1)
        out <- ode(func = caldhal, y = yini, times = times, parms = NULL, method = "bdf")
        n <- nfun(out[,2], q[sim])
        C <- n/a
        recsim <- rbind(recsim, cbind(subrf, subev, C, out, n, IDen, sim, L = L[sim], q = q[sim], v = v[sim]))
      }
    }
  }
  return(recsim)
}
