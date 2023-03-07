z0ndnlsnslp <- function(wsp, ntr, tim, canh) {
    z0 <- 0
    d <- 0
    err <- 0
    wspabl <- wsp[wsp$H >= canh,]
    for(i in 1:length(tim)){
      if(tim[i] %in% ntr){
        usethis <- wspabl[wspabl$time == tim[i],-3]
        if(sum(is.na(usethis$ustar)) == 0 & sum(is.na(usethis$WS)) < 2 ){
          model1 <- nls(WS ~ (ustar)/0.4*log((H-b)/a), data = usethis,
                        start = list(a = 0.1, b = 0.1), trace = F, na.action(na.exclude(1)),
                        algorithm = "port", lower = c(10^-6, 0.005), upper = c(4, (canh-0.01)),
                        control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
          z0[i] <- summary(model1)$coefficients[1]
          d[i] <- summary(model1)$coefficients[2]
          err[i] <- summary(model1)$sigma
          # if(z0[i] == 10^-6 | d[i] == 0.005 | z0[i] == 4 | d[i] == (canh-0.01)){
          #   z0[i] <- NA
          #   d[i] <- NA
          # }
          usethis <- usethis[usethis$H== max(usethis$H),]
        }
        else{
          z0[i] <- NA
          d[i] <- NA
          err[i] <- NA
        }
      }
      else{
        z0[i] <- NA
        d[i] <- NA
        err[i] <- NA
      }
    }
    return(list(z0, d, err))
}
