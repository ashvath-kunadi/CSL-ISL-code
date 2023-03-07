phimfnc <- function(Et0dat, trydd, dudz, senh, monobhL, plotcond = F, plotloc = "H:/My Documents/Project and related work/Code/Papier Code", plotprefix = "Site"){
  dscn <- length(unique(dudz$zeval))
  betaus <- 0
  errus <- 0 
  betast <- 0
  errst <- 0
  algo <- c("lm","M","MM", "tau", "CM", "mtl")
  for(i in 1:5){
    dsub <- trydd[,(i+1)]
    zetafi <-   (senh-dsub)/monobhL
    for(j in 1:dscn){
      dudzsub <- dudz[dudz$zeval == unique(dudz$zeval)[j],]
      phism <- 0.4/Et0dat$ustar*(dudzsub$zeval-dsub)*dudzsub$dscdz
      sub <- data.frame(zeta = zetafi, phism)
      sub <- sub[complete.cases(sub),]
      for(k in 2:10){
        unstb <- sub[sub$zeta <= 0 & sub$zeta >= -k,]
        stb <- sub[sub$zeta > 0 & sub$zeta <= k,]
        if(k == 10){
          unstb <- sub[sub$zeta < 0,]
          stb <- sub[sub$zeta > 0,]
        }
        for(l in 1:6){
          if(l == 1){
            model1 <- nls(phism ~ (1-(b*zeta))^(-1/4), data = unstb, start = list(b = 16), algorithm = "port", lower = c(0.05), upper = c(100))
            model2 <- nls(phism ~ (1+(b*zeta)), data = stb, start = list(b = 6), algorithm = "port", lower = c(0.05), upper = c(100))
            betaus[((i-1)*54*dscn)+((j-1)*54)+((k-2)*6)+l] <- summary(model1)$coefficients[1]
            betast[((i-1)*54*dscn)+((j-1)*54)+((k-2)*6)+l] <- summary(model2)$coefficients[1]
          }else if(l == 2){
            model1 <- nlrob(phism ~ (1-(b*zeta))^(-1/4), data = unstb, method = "M", start = list(b = 16), algorithm = "port", lower = c(0.05), upper = c(100))
            model2 <- nlrob(phism ~ 1+(b*zeta), data = stb, method = "M", start = list(b = 6), algorithm = "port", lower = c(0.05), upper = c(100))
            betaus[((i-1)*54*dscn)+((j-1)*54)+((k-2)*6)+l] <- model1$coefficients
            betast[((i-1)*54*dscn)+((j-1)*54)+((k-2)*6)+l] <- model2$coefficients
          }else{
            model1 <- nlrob(phism ~ (1-(b*zeta))^(-1/4), data = unstb, method = algo[l], lower = c(b = 0.05), upper = c(b = 100))
            model2 <- nlrob(phism ~ 1+(b*zeta), data = stb, method = algo[l], lower = c(b = 0.05), upper = c(b = 100))
            betaus[((i-1)*54*dscn)+((j-1)*54)+((k-2)*6)+l] <- model1$coefficients
            betast[((i-1)*54*dscn)+((j-1)*54)+((k-2)*6)+l] <- model2$coefficients
          }
        }
        if(plotcond){
          x <- getwd()
          phismbr <- function(zee){
            len <- length(zee)
            ans <- rep(NA, len)
            for(p in 1:len){
              if(!is.na(zee[p])){
                if(zee[p] <= 0){
                  ans[p] <- (1-16*zee[p])^-0.25
                }else if(zee[p] > 0 & zee[p] <= 1){
                  ans[p] <- 1 + 5*zee[p]
                }else{
                  ans[p] <- 6
                }
              }else{
                ans[p] <- NA
              }
            }
            return(ans)
          }
          phismfk <- function(zee){
            len <- length(zee)
            ans <- rep(NA, len)
            for(p in 1:len){
              if(zee[p] <= 0){
                ans[p] <- (1-19.3*zee[p])^-0.25
              }else{
                ans[p] <- 1 + 6*zee[p]
              }
            }
            return(ans)
          }
          phismfn <- function(zee, bs, bus){
            len <- length(zee)
            ans <- rep(NA, len)
            for(p in 1:len){
              if(zee[p] <= 0){
                ans[p] <- (1-bus*zee[p])^-0.25
              }else if(zee[p] > 0){
                ans[p] <- 1 + bs*zee[p]
              }
            }
            return(ans)
          }
          schn <- c("nls", "m", "mm", "cm", "tau", "bustmtl")
          schc <- c( "coral","brown", "blue", "green", "orange", "purple")
          g <- ggplot(sub, aes(zeta, phism)) + xlim(c(-5,5)) +  ylim(c(0, 10)) + 
              geom_point(alpha = 0.01) + theme(legend.position='none') + 
              ggtitle(paste("phi and zeta with ", as.character(names(trydd)[i+1]), " at ", unique(dudz$zeval)[j], " zeta limit at ", k))
          for(l in 1:length(unique(schn))){
              g <- g + stat_function(fun = phismfn, n = 1000, args = list(betast[((i-1)*162)+((j-1)*54)+((k-2)*6)+l], betaus[((i-1)*162)+((j-1)*54)+((k-2)*6)+l]), colour = schc[l]) + annotate("text", x = -1, y = 10-l/2, label = paste(schn[l], " Stb b -", round(betast[((i-1)*162)+((j-1)*54)+((k-2)*6)+l],1), " Unstb b -",round(betaus[((i-1)*162)+((j-1)*54)+((k-2)*6)+l],1)), col = schc[l])
          }
          g <- g + theme_bw() + labs(x = expression(zeta), y = expression(phi[m])) + stat_function(fun = phismfk, col = "red") + stat_function(fun = phismbr, col = "darkolivegreen") + theme_minimal() + labs(x = expression(zeta), y = expression(phi[m]))
          setwd(plotloc)
          jpeg(paste(plotprefix, " phi and zeta with ", as.character(names(trydd)[i+1]), " at ", unique(dudz$zeval)[j], " zeta limit at ", k,".jpg"))
          print(g)
          dev.off()
          setwd(x)
        }
      }
    }
  }
  
  ansdf <- data.frame(bu = betaus, bs = betast, dsch = rep(names(trydd)[2:6], each = dscn*9*6), 
                      zeval = rep(rep(unique(dudz$zeval), each = 9*6), 5), ztlim = rep(rep(2:10, each = 6), dscn*5),
                      algo)
  return(ansdf)
}
