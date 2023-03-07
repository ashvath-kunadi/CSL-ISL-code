momentum <- function(Et0dat, wdp,wsp,tap,ahp, stb, dudz, dTdz, dahdz, senh, canh, plotcond = F, plotloc = "H:/My Documents/Project and related work/Code/Papier Code", plotprefix = "Site", ntrlim = 0.05){
  library(corrplot)
  library(robustbase)
  library(plyr)
  library(dplyr)
  library(MASS)
  library(ggplot2)
  library(ggExtra)
  tim <- Et0dat$time
  n <- dim(Et0dat)[1]
  ntr <- tim[stb == "D"]
  #ntr <- tim[dtdzpsqstb == "D"]
  #ntrcond <- (dtdzpsqstbas == "D")
  # monobhL <- with(Et0dat, -(ustar^3*rhom)/(0.4*9.80665*((H/(Tair+273)/Cpm)+(0.61*FLh/Lheat))))
  monobhL <- with(Et0dat, monobhl)
  ##first estimate stull 382 wind speed
  est1 <- z0ndnlsnslp(wsp, ntr, tim, canh)
  
  ##second brutsaert 4.5
  # est2 <- numderd(dudz, Et0dat, n, tim, ntr)
  
  ##looks like the estimates match
  cd1 <- mean(est1[[2]], na.rm = T)
  cz01 <- mean(est1[[1]], na.rm = T)
  vard1 <- est1[[2]]
  varz01 <- est1[[1]]
  
  ##testing whether the errors were due to any other factors
  ##stability checks
  cd2 <- cd1 - 0.020
  first <- T
  count <- 1
  while(abs(cd1-cd2) > 0.01 & count < 10){
    zeta1 <- (senh-cd1)/monobhL
    fkn2017cat <- rep("VUST", length(zeta1))
    fkn2017cat[zeta1 >= -1] <- "UST"
    fkn2017cat[zeta1 >= -ntrlim] <- "NTR"
    fkn2017cat[zeta1 >= ntrlim] <- "STB"
    fkn2017cat[zeta1 >= 1] <- "VSTB"
    fkn2017cat[zeta1 >= 7] <- NA
    fkn2017cat[zeta1 <= -5] <- NA
    ##new stability based on L
    ntr <- tim[fkn2017cat == "NTR"]
    ntrcond <- (fkn2017cat == "NTR")
    
    ##first estimate stull 382 wind speed
    est3 <- z0ndnlsnslp(wsp, ntr, tim, canh)
    cd1 <- cd2
    vard2 <- est3[[2]]
    cd2 <- mean(vard2, na.rm = T)
    if(first){
      first <- F
    }
    count <- count+1
  }
  if(plotcond){
    x <- getwd()
    top <- data.frame(fkn2017cat, stb, zeta1)
    top$fkn2017cat <- factor(fkn2017cat, levels = c("VUST", "UST", "NTR", "STB", "VSTB"))
    g1 <- ggplot(top, aes(x = fkn2017cat, fill = stb)) + geom_bar(position = "stack", stat = "count") + xlab(expression(paste(zeta, " Stability Criteria"))) + theme_bw() + labs(fill = paste("Temperature Profile", "\n", "Stability Criteria", "\n", "Unstable to Stable"))
    g2 <- ggplot(top, aes(x = fkn2017cat, fill = stb)) + geom_bar(position = "fill", stat = "count") + theme_bw()
    setwd(plotloc)
    jpeg(paste(plotprefix, 'Stability comparison 1.jpg'))
    print(g1)
    dev.off()
    jpeg(paste(plotprefix, 'Stability comparison 2.jpg'))
    print(g2)
    dev.off()
    setwd(x)
  }
  ##second brutsaert 4.5
  est4 <- numderd(dudz, Et0dat, n, tim, ntr)
  if(plotcond){
    x <- getwd()
    setwd(plotloc)
    top <- data.frame(dog = vard2, est4)
    g1 <- ggplot(top, aes(dog, d0)) + facet_grid(.~cmb) + geom_point(alpha = 0.2) + ylim(0,10) + geom_smooth() + coord_flip() + theme_bw()   +labs(x=expression(Old~~d[0]), y=expression(New~~d[0]))
    # ctop <- top[complete.cases(top),]
    # ctop <- ctop[ctop$d0 < 10 & ctop$d0 >= 0,]
    # kd <- with(ctop[ctop$cmb == 1,1:2], MASS::kde2d(dog, d0, n = 50))
    # plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface()
    # setwd(plotloc)
    jpeg(paste(plotprefix, 'd0 check 1.jpg'))
    print(g1)
    dev.off()
    setwd(x)
  }
  varz02 <- est3[[1]]
  cz02 <- mean(varz02, na.rm = T)
  
  # beta <- wsp$ustar[wsp$H == 7]/wsp$WS[wsp$H == 7]
  # tcond <- !is.na(vard2)
  # try1 <- data.frame(ws7 = wsp$WS[wsp$H == 7], ws19 = wsp$WS[wsp$H == 1.9], beta)
  # try2 <- try1[tcond,]
  # model2 <- nls(ws7 ~ ((ws7/(besselI(2*sqrt(a/beta/0.4),0)-(besselI(2*sqrt(z0g*a/canh/beta/0.4),0)*besselK(2*sqrt(a/beta/0.4),0)/besselK(2*sqrt(z0g*a/canh/beta/0.4),0))))*besselI(2*sqrt(a/beta/0.4), 0)+((ws7/(besselI(2*sqrt(a/beta/0.4),0)-(besselI(2*sqrt(z0g*a/canh/beta/0.4),0)*besselK(2*sqrt(a/beta/0.4),0)/besselK(2*sqrt(z0g*a/canh/beta/0.4),0))))*besselI(2*sqrt(z0g*a/canh/beta/0.4),0)/besselK(2*sqrt(z0g*a/canh/beta/0.4),0))*besselK(2*sqrt(a/beta/0.4),0)), data = try2,
  #               start = list(z0g = 0.01, a = 0.1), trace = F, na.action(na.exclude(1)),
  #               algorithm = "port", lower = c(10^-6, 10^-6, 10^-6), upper = c(1000, 1000, 1000),
  #               control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))  ##testing whether the errors were due to any other factors
  # model1 <- nls(ws19 ~ ((ws7/(besselI(2*sqrt(a/beta/0.4),0)-(besselI(2*sqrt(z0g*a/canh/beta/0.4),0)*besselK(2*sqrt(a/beta/0.4),0)/besselK(2*sqrt(z0g*a/canh/beta/0.4),0))))*besselI(2*sqrt(1.9/7*a/beta/0.4), 0)+((ws7/(besselI(2*sqrt(a/beta/0.4),0)-(besselI(2*sqrt(z0g*a/canh/beta/0.4),0)*besselK(2*sqrt(a/beta/0.4),0)/besselK(2*sqrt(z0g*a/canh/beta/0.4),0))))*besselI(2*sqrt(z0g*a/canh/beta/0.4),0)/besselK(2*sqrt(z0g*a/canh/beta/0.4),0))*besselK(2*sqrt(1.9/7*a/beta/0.4),0)), data = try2,
  #               start = list(z0g = 0.01, a = 0.1), trace = F, na.action(na.exclude(1)),
  #               algorithm = "port", lower = c(10^-6, 10^-6, 10^-6), upper = c(1000, 1000, 1000),
  #               control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))  ##testing whether the errors were due to any other factors
  # 
  # topcor <- data.frame(d0 = vard2, z0 = varz02, logratio = log((senh-vard2)/varz02), u = Et0dat$u, ustar = Et0dat$ustar, stability = zeta1, error = est3[[3]], bta = Et0dat$ustar/Et0dat$u)
  # if(plotcond){
  #   x <- getwd()
  #   setwd(plotloc)
  #   g1 <- ggplot(topcor, aes(d0)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +geom_vline(xintercept = (2/3*canh)) + theme_bw() + labs(x =  expression(d[0]))
  #   jpeg(paste(plotprefix, 'd0 density.jpg'))
  #   print(g1)
  #   dev.off()
  #   g1 <- ggplot(topcor, aes(z0)) + geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +geom_vline(xintercept = (0.1*canh)) + theme_bw() + labs(x =  expression(z[0]))
  #   jpeg(paste(plotprefix, 'z0 density.jpg'))
  #   print(g1)
  #   dev.off()
  #   
  #   topcorg <- topcor[complete.cases(topcor),]
  #   cor_matrix <- cor(topcorg, use = "complete.obs")
  #   g1 <- corrplot.mixed(cor_matrix, lower = "circle", upper = "number", tl.pos = "lt", diag = "u")
  #   jpeg(paste(plotprefix, 'd0 z0 correlation.jpg'))
  #   print(g1)
  #   dev.off()
  #   mtop <- melt(topcorg, id = "d0")
  #   g1 <- ggplot(mtop, aes(value, d0)) + facet_wrap(.~variable, scales = "free") + geom_point(alpha = 0.01) + theme_bw()  +labs(y=expression(paste(d[0], " m")))
  #   jpeg(paste(plotprefix, 'd0 variation.jpg'))
  #   print(g1)
  #   dev.off()
  #   setwd(x)
  #   try <- topcorg[,c(1,6)]
  #   try <- try[order(abs(try$stability)),]
  #   try$AbsoluteAverage <- cumsum(try$d0)/(1:length(try$d0))
  #   try <- try[order(try$stability),]
  #   try$GrossAverage <- cumsum(try$d0)/(1:length(try$d0))
  #   try2 <- melt(try, id = c("stability", "d0"))
  #   g1 <- ggplot(try2, aes(stability, value)) + facet_grid(.~variable) + geom_line() + theme_bw() + labs(x =  expression(zeta), y = expression(d[0]))
  #   jpeg(paste(plotprefix, 'd0 variation w zeta.jpg'))
  #   print(g1)
  #   dev.off()
  #   
  # top2 <- fundf

   #   jpeg(paste(plotprefix, 'd0 variation w u.jpg'))
  #   print(g1)
  #   dev.off()
  #   
  #   
  #   jpeg(paste(plotprefix, 'z0 variation w ustar.jpg'))
  #   print(g1)
  #   dev.off()
  #   setwd(x)
  #   
  # }
  # 
  # sampz0 <- cbind(topcor[,-c(1,3,4,5)], Et0dat)
  # sampz0g <- sampz0[complete.cases(sampz0),]
  # compmod <-  lm(z0~., data = sampz0g)
  # step.model <- stepAIC(compmod, direction = "both")
  # 
  # sampd0 <- cbind(topcor[,-c(2,3,4,5)], Et0dat)
  # sampd0g <- sampd0[complete.cases(sampd0),]
  # compmod <-  lm(d0~., data = sampd0g)
  # step.model <- stepAIC(compmod, direction = "both")
  # 
  # dmodel1 <- lm(vard2~u, data = topcor)
  # dmodel2 <- lm(vard2~u+ustar, data = topcor)
  # dmodel3 <- nls(vard2 ~ (1-(2*(bta)^3*a/0.4/canh))*canh, data = topcor,
  #                start = list(a = 0.1), trace = F, na.action(na.exclude(1)),
  #                algorithm = "port", lower = c(10^-6), upper = c(10000),
  #                control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
  # dmodel4 <- nls(vard2 ~ 7 - (ustar*2*7/(0.4*2*
  #                                          sqrt (a/bta/
  #                                                  0.4)*((ws7/(besselI (2*sqrt (a/bta/0.4), 
  #                                                                       0) - (besselI (2*sqrt (z0g*a/canh/bta/0.4), 0)*
  #                                                                               besselK (2*sqrt (a/bta/0.4), 0)/
  #                                                                               besselK (2*sqrt (z0g*a/canh/bta/0.4), 0))))*
  #                                                          besselI (2*sqrt (a/bta/0.4), 
  #                                                                   1) - ((ws7/(besselI (2*sqrt (a/bta/0.4), 
  #                                                                                        0) - (besselI (2*sqrt (z0g*a/canh/bta/0.4), 0)*
  #                                                                                                besselK (2*sqrt (a/bta/0.4), 0)/
  #                                                                                                besselK (2*sqrt (z0g*a/canh/bta/0.4), 0))))*
  #                                                                           besselI (2*sqrt (z0g*a/canh/bta/0.4), 0)/
  #                                                                           besselK (2*sqrt (z0g*a/canh/bta/0.4), 0))*
  #                                                          besselK (2*sqrt (a/bta/0.4), 1)))), data = topcor,
  #                start = list(z0g = 0.01, a = 0.1), trace = F, na.action(na.exclude(1)),
  #                algorithm = "port", lower = c(10^-6, 10^-6), upper = c(1000, 1000),
  #                control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))  ##testing whether the errors were due to any other factors
  # topcor4 <- topcor
  # topcor4$bta[(topcor$bta < 0.2)] <- NA
  # dmodel4 <- nls((7-d0) ~ (ustar*2*7/(0.4*2*
  #                                          sqrt(a/bta/0.4)*(c1*besselI (2*sqrt (a/bta/0.4),0)) + 
  #                                                          (c2*besselK (2*sqrt (a/bta/0.4),0)))), data = topcor4,
  #                start = list(a = 0.01, c1 = 0.1, c2 = 0.1), trace = F, na.action(na.exclude(1)),
  #                algorithm = "port", lower = c(10^-6, 10^-6, 10^-6), upper = c(100, 10000, 10000),
  #                control = nls.control(maxiter = 5000, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))  ##testing whether the errors were due to any other factors
  # 
  # ##new model2
  # #read dinghal and naraynmurthy 2019 and wang 2012 to come up with the forcing and response necessary
  # trydd <- data.frame(d0 = vard2, conv = (2/3*canh), cd = cd2, mod1 = predict(dmodel1, newdata = data.frame(u = Et0dat$u)), 
  #                     mod2 = predict(dmodel2, newdata = data.frame(u = Et0dat$u, ustar = Et0dat$ustar)),
  #                     mod3 = predict(dmodel3, newdata = data.frame(bta = Et0dat$ustar/Et0dat$u)))
  # # mdd <- melt(trydd, id = "d0")
  # names(mdd) <- c("d0", "d0_scheme", "d0_est")
  # mdd$err <- mdd$d0-mdd$d0_est
  # mddtop <- melt(mdd, id = c("d0", "d0_scheme"))
  
  ##try fixing d0 and then see the values that you get for z0
  newz0 <-  onlz0(ntr, cd2, canh, wsp, tim)
  varz03 <- newz0$z0fd

  if(plotcond){
    x <- getwd()
    setwd(plotloc)
    topcor2 <- data.frame(z0 = varz02, z0fd = varz03, Thmr = log((senh-cd2)/varz03), u = Et0dat$u, ustar = Et0dat$ustar, zeta = zeta1, err = est3[[3]])
    topcorg2 <- topcor2[complete.cases(topcor2),]
    comp <- data.frame(z0wd = topcor$logratio, z0wod = topcor2$Thmr)
    g1 <- ggplot(comp, aes(z0wd, z0wod)) + geom_point(alpha = 0.1, aes(col = Et0dat$ustar)) + theme_minimal() + labs(x = expression(paste("Thompson ratio with variable ", d[0])), y = expression(paste("Thompson ratio with constant ", d[0])), main = expression(paste("Comparison of ", z[0]))) + geom_abline(slope = 1, intercept = 0) + geom_smooth(method = "lm", formula = y~0+x) + geom_smooth(method = "lm", col = "red")
    jpeg(paste(plotprefix, 'z0 check 1.jpg'))
    print(g1)
    dev.off()
    cor_matrix <- cor(topcorg2, use = "complete.obs")
    g2 <- corrplot.mixed(cor_matrix, lower = "circle", upper = "number", tl.pos = "lt", diag = "u")
    jpeg(paste(plotprefix, 'z0 correlation.jpg'))
    print(g2)
    dev.off()
    mtop <- melt(topcorg2, id = "z0")
    g1 <- ggplot(mtop, aes(value, z0)) + facet_grid(.~variable, scales = "free") + geom_point(alpha = 0.01) + theme_bw()
    jpeg(paste(plotprefix, 'z0 variation.jpg'))
    print(g1)
    dev.off()
    setwd(x)
  }
  
  # z0model1 <- lm(z0~u, data = topcor)
  # z0model2 <- lm(z0~u+ustar, data = topcor)
  # z0model3 <- nls(z0 ~ (2*bta^3*a/0.4*exp(-0.4/bta)), data = topcor,
  #                           start = list(a = 0.1), trace = F, na.action(na.exclude(1)),
  #                           algorithm = "port", lower = c(10^-6), upper = c(10000),
  #                           control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
  # z0model4 <- nls(log((7-d0)/z0) ~ (0.4/ustar*((c1*besselI (2*sqrt (a/bta/0.4),1)) - 
  #                                          (c2*besselK (2*sqrt (a/bta/0.4),1)))), data = topcor4,
  #                start = list(a = 0.01, c1 = 0.1, c2 = 0.1), trace = F, na.action(na.exclude(1)),
  #                algorithm = "port", lower = c(10^-6, 10^-6, 10^-6), upper = c(100, 10000, 10000),
  #                control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))  ##testing whether the errors were due to any other factors
  # 
  # ##new models
  # #usm1 <- 1/topcor$ustar
  # #us2 <- (topcor$ustar)^2
  # #z0model3 <- lm(z0)
  # tryzd <- data.frame(z0 = varz02, conv = (0.1*canh), cz0 = cz02, mod1 = predict(z0model1, newdata = data.frame(u = Et0dat$u)), 
  #                     mod2 = predict(z0model2, newdata = data.frame(u = Et0dat$u, ustar = Et0dat$ustar)),
  #                     mod3 = predict(z0model3, newdata = data.frame(bta = Et0dat$ustar/Et0dat$u)))
  # kvnh <- khnkvcalc(Et0dat, vard2, wsp, tap, ahp,  dudz, dTdz, dahdz, fkn2017cat, senh, canh)
  # ##phim calc
  # ansdf <- phimfnc(Et0dat, trydd, dudz, senh, monobhL, plotcond, plotloc, plotprefix)
  ##less go
  ##getting to phism
  
  #return(list(fkn2017cat, cz01, cd1, cz02,cd2, varz02, vard2, z0model1, z0model2, dmodel1, dmodel2, ansdf, kvnh))
  return(list(fkn2017cat, cz01, cd1, cz02,cd2, varz02, vard2, varz03, est4))
  
}  
    