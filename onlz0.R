onlz0 <- function(ntr, fixedd, canh, wsp, tim){
  z0fd <- 0
  err <- 0
  wspabl <- wsp[wsp$H >= canh,]
  for(i in 1:length(tim)){
    if(tim[i] %in% ntr){
      usethis <- wspabl[wspabl$time == tim[i],-3]
      if(sum(is.na(usethis$ustar)) == 0 & sum(is.na(usethis$WS)) < 2 ){
        model1 <- nls(WS ~ (ustar)/0.4*log((H-fixedd)/a), data = usethis,
                      start = list(a = 0.1), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(10^-6), upper = c(4),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0fd[i] <- summary(model1)$coefficients[1]
        err[i] <- summary(model1)$sigma
      }
      else{
        z0fd[i] <- NA
        err[i] <- NA
      }
      
    }
    else{
      z0fd[i] <- NA
      err[i] <- NA
    }
  }
  top2 <- data.frame(z0fd, err)
}
