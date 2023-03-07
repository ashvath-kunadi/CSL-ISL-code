z0hnlsfin <- function(khcon, khvar, tap, d, Et0dat, tim, IR = F){
  z0h1n <- 0
  errh1 <- 0
  z0h2n <- 0
  errh2 <- 0
  z0h3n <- 0
  errh3 <- 0
  dh3 <- 0
  z0h4n <- 0
  errh4 <- 0
  z0h5n <- 0
  errh5 <- 0
  z0h6n <- 0
  errh6 <- 0
  dh6 <- 0
  z0h7n <- 0
  errh7 <- 0
  z0h8n <- 0
  errh8 <- 0
  z0h9n <- 0
  errh9 <- 0
  dh9 <- 0
  tapabl <- tap[tap$H > 5,c(1,2,7)]
  tapabl$kHcon <- khcon
  tapabl$kHvar <- khvar
  tapabl$vard <- d
  cd <- mean(d, na.rm = T)
  if(IR == T){
    sur <- data.frame(time = Et0dat[,5], Ta = Et0dat[,19])
  }else{
    sur <- tap[tap$H == 1.25, c(2,7)]
  }
  for(i in 1:length(tim)){
    if(tim[i] %in% ntr){
      usethis <- tapabl[tapabl$time == tim[i],-3]
      usethis$cd <- cd
      usethis$FH <- Et0dat$H[Et0dat$time == tim[i]]
      usethis$Cpm <- Et0dat$Cpm[Et0dat$time == tim[i]]
      usethis$rhom <- Et0dat$rhom[Et0dat$time == tim[i]]
      usethis$ustar <- Et0dat$ustar[Et0dat$time == tim[i]]
      usethis$sur <- sur$Ta[sur$time == tim[i]]
      if(sum(complete.cases(usethis)) == 3){
        model1 <- nls((sur-Ta) ~ (FH)/0.4/ustar/Cpm/rhom*log((H-cd)/a), data = usethis,
                      start = list(a = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24), upper = c(8),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h1n[i] <- summary(model1)$coefficients[1]
        errh1[i] <- summary(model1)$sigma
        model2 <- nls((sur-Ta) ~ (FH)/0.4/ustar/Cpm/rhom*log((H-vard)/a), data = usethis,
                      start = list(a = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24), upper = c(8),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h2n[i] <- summary(model2)$coefficients[1]
        errh2[i] <- summary(model2)$sigma
        model3 <- nls((sur-Ta) ~ (FH)/0.4/ustar/Cpm/rhom*log((H-b)/a), data = usethis,
                      start = list(a = 0.01, b = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24, 10^-3), upper = c(8, 7),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h3n[i] <- summary(model3)$coefficients[1]
        dh3[i] <- summary(model3)$coefficients[2]
        errh3[i] <- summary(model3)$sigma
        model4 <- nls((sur-Ta) ~ (FH)/kHcon/ustar/Cpm/rhom*log((H-cd)/a), data = usethis,
                      start = list(a = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24), upper = c(8),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h4n[i] <- summary(model4)$coefficients[1]
        errh4[i] <- summary(model4)$sigma
        model5 <- nls((sur-Ta) ~ (FH)/kHcon/ustar/Cpm/rhom*log((H-vard)/a), data = usethis,
                      start = list(a = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24), upper = c(8),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h5n[i] <- summary(model5)$coefficients[1]
        errh5[i] <- summary(model5)$sigma
        model6 <- nls((sur-Ta) ~ (FH)/kHcon/ustar/Cpm/rhom*log((H-b)/a), data = usethis,
                      start = list(a = 0.01, b = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24, 10^-3), upper = c(8, 7),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h6n[i] <- summary(model6)$coefficients[1]
        dh6[i] <- summary(model6)$coefficients[2]
        errh6[i] <- summary(model6)$sigma
        model7 <- nls((sur-Ta) ~ (FH)/kHvar/ustar/Cpm/rhom*log((H-cd)/a), data = usethis,
                      start = list(a = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24), upper = c(8),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h7n[i] <- summary(model7)$coefficients[1]
        errh7[i] <- summary(model7)$sigma
        model8 <- nls((sur-Ta) ~ (FH)/kHvar/ustar/Cpm/rhom*log((H-vard)/a), data = usethis,
                      start = list(a = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24), upper = c(8),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h8n[i] <- summary(model8)$coefficients[1]
        errh8[i] <- summary(model8)$sigma
        model9 <- nls((sur-Ta) ~ (FH)/kHvar/ustar/Cpm/rhom*log((H-b)/a), data = usethis,
                      start = list(a = 0.01, b = 0.01), trace = F, na.action(na.exclude(1)),
                      algorithm = "port", lower = c(8.432798e-24, 10^-3), upper = c(8, 7),
                      control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
        z0h9n[i] <- summary(model9)$coefficients[1]
        dh9[i] <- summary(model9)$coefficients[2]
        errh9[i] <- summary(model9)$sigma
      }
      else{
        z0h1n[i] <- NA
        errh1[i] <- NA
        z0h2n[i] <- NA
        errh2[i] <- NA
        z0h3n[i] <- NA
        errh3[i] <- NA
        dh3[i] <- NA
        z0h4n[i] <- NA
        errh4[i] <- NA
        z0h5n[i] <- NA
        errh5[i] <- NA
        z0h6n[i] <- NA
        errh6[i] <- NA
        dh6[i] <- NA
        z0h7n[i] <- NA
        errh7[i] <- NA
        z0h8n[i] <- NA
        errh8[i] <- NA
        z0h9n[i] <- NA
        errh9[i] <- NA
        dh9[i] <- NA
      }
    }
    else{
      z0h1n[i] <- NA
      errh1[i] <- NA
      z0h2n[i] <- NA
      errh2[i] <- NA
      z0h3n[i] <- NA
      errh3[i] <- NA
      dh3[i] <- NA
      z0h4n[i] <- NA
      errh4[i] <- NA
      z0h5n[i] <- NA
      errh5[i] <- NA
      z0h6n[i] <- NA
      errh6[i] <- NA
      dh6[i] <- NA
      z0h7n[i] <- NA
      errh7[i] <- NA
      z0h8n[i] <- NA
      errh8[i] <- NA
      z0h9n[i] <- NA
      errh9[i] <- NA
      dh9[i] <- NA
    }
  }
  z0h1n[z0h1n == 8.432798e-24 | z0h1n == 8] <- NA
  z0h2n[z0h2n == 8.432798e-24 | z0h2n == 8] <- NA
  z0h3n[z0h3n == 8.432798e-24 | z0h3n == 8] <- NA
  z0h4n[z0h4n == 8.432798e-24 | z0h4n == 8] <- NA
  z0h5n[z0h5n == 8.432798e-24 | z0h5n == 8] <- NA
  z0h6n[z0h6n == 8.432798e-24 | z0h6n == 8] <- NA
  z0h7n[z0h7n == 8.432798e-24 | z0h7n == 8] <- NA
  z0h8n[z0h8n == 8.432798e-24 | z0h8n == 8] <- NA
  z0h9n[z0h9n == 8.432798e-24 | z0h9n == 8] <- NA
  ans <- data.frame(z0h = c(z0h1n, z0h2n, z0h3n, z0h4n, z0h5n,z0h6n,z0h7n,z0h8n,z0h9n), 
                    d = c(rep(cd, i), d, dh3, rep(cd, i), d, dh6, rep(cd, i), d, dh9), 
                    kscheme = c(rep(c("k", "khc", "khv"), each = i*3)))
  return(ans)
}
