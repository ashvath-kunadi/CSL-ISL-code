rutsim <- function(usedat, gpfldat, dt = 5){
  library(data.table)
  library(deSolve)
  Ssim <- 0
  psim <- 0 
  bsim <- 0
  asim <- 0 
  recsim <- data.frame()
  rfdat <- usedat[[1]]
  ET0dat <- usedat[[2]]
  drainage <- function(C, a, b){
      ans <- exp(a + b*C)
      if(ans > C){
        ans <- C
      }
      return(ans)
  }
  evap <- function(S, C, ev){
    if(C > S){
      ans <- ev
    }
    else{
      ans <- ev*C/S
    }
    return(ans)
  }
  for(sim in 1:500){
      Ssim[sim] <- runif(1, min = 0.1, max = 3.5)
      psim[sim] <- runif(1, min = 0.01, max = 0.8)
      bsim[sim] <- runif(1, min = 1, max = 6)
      asim[sim] <- log(0.002)-bsim[sim]*Ssim[sim]
      for(k in 1:length(unique(rfdat$ID))){
        IDen <- unique(rfdat$ID)[k]
        subrf <- rfdat[(rfdat$ID == IDen),]
        subrf <- subrf[-c(1:59),2]
        subev <- ET0dat[(ET0dat$ID == IDen),]
        subev <- subev[-c(1:59),2]
            if(sum(is.na(subev)) == 0){
              S <- Ssim[sim]
              p <- psim[sim]
              b <- bsim[sim]
              a <- asim[sim]
              caldhal <- function(t, y, parms){
                with(as.list(y), {
                  dPdt <- subrf[t]
                  dIdt <- subev[t]
                  dCdt <- (1-p)*dPdt-drainage(C, a, b)-evap(S, C, dIdt)
                  dThdt <- p*dPdt + drainage(C, a, b)
                  list(dCdt, dThdt)
                })
              }
              yini <- c(C=0)
              times <- seq(from = 1, to = length(subrf), by = 1)
              out <- ode(func = caldhal, y = yini, times = times, parms = NULL, method = "bdf")
              recsim <- rbind(recsim, cbind(subrf, subev, S, p, b, a, out, IDen, sim))
            }
      }
  }
  return(recsim)
}

