khnkvcalc <- function(Et0dat, vd, wsp, tap, ahp, dudz, dTdz, dahdz, stb, senh, canh){
  ##estimating value of av and ah
  kdf <- data.frame()
  kvdf <- data.frame()
  khdf <- data.frame()
  cd <- mean(vd, na.rm = T)
  uapab <- dudz[dudz$zeval > canh, ]
  tapab <- dTdz[dTdz$zeval > canh, ]
  ahpab <- dahdz[dahdz$zeval > canh,]
  udux <- length(unique(uapab$zeval))
  udta <- length(unique(tapab$zeval))
  udah <- length(unique(ahpab$zeval))
  dlab <- c("const", "var")
  dmat <- cbind(cd, vd)
  ##argument for using only the middl eone is that the d0 values from the middle one gave the best values, so yeah
  for(j in 1:udux){
    subdz <- uapab[uapab$zeval==unique(uapab$zeval)[j],]
    for(i in 1:2){
      dsub <- dmat[,i]
      subsen1 <- wsp[wsp$H == subdz$sen1[1],]
      subsen2 <- wsp[wsp$H == subdz$sen2[1],]
      k45 <- Et0dat$ustar/(senh-dsub)/subdz$dscdz
      k46 <- Et0dat$ustar/(subsen2[,2]-subsen1[,2])*log((subsen2[,1]-dsub)/(subsen1[,1]-dsub))
      ans <-  data.frame(k = c(k45, k46), d = dsub, dlab = dlab[i], z = subdz$zeval, frm = rep(c("4.5","4.6"), each = length(k45)))
      kdf <- rbind(kdf, ans)
    }
  }
  for(j in 1:udah){
    subdz <- ahpab[ahpab$zeval==unique(ahpab$zeval)[j],]
    for(i in 1:2){
      dsub <- dmat[,i]
      subsen1 <- ahp[ahp$H == subdz$sen1[1],]
      subsen2 <- ahp[ahp$H == subdz$sen2[1],]
      x <- cbind(subsen1[,6],subsen1[,6])
      subrhom <- apply(x, 1, mean)
      kv411 <- -Et0dat$FLh/Et0dat$Lheat/Et0dat$ustar/subrhom/(senh-dsub)/subdz$dscdz
      kv413 <- Et0dat$FLh/Et0dat$Lheat/Et0dat$ustar/subrhom/(subsen1[,5]-subsen2[,5])*log((subsen2[,1]-dsub)/(subsen1[,1]-dsub))
      ans <-  data.frame(kv = c(kv411, kv413), d = dsub, dlab = dlab[i], z = subdz$zeval, frm = rep(c("4.11","4.13"), each = length(kv411)))
      kvdf <- rbind(kvdf, ans)
    }
  }
  for(j in 1:udta){
    subdz <- tapab[tapab$zeval==unique(tapab$zeval)[j],]
    for(i in 1:2){
      dsub <- dmat[,i]
      kh415 <- -Et0dat$H/Et0dat$Cpm/Et0dat$rhom/Et0dat$ustar/(senh-dsub)/subdz$dscdz
      subsen1 <- tap[tap$H == subdz$sen1[1],]
      subsen2 <- tap[tap$H == subdz$sen2[1],]
      kh416 <- Et0dat$H/Et0dat$Cpm/Et0dat$rhom/Et0dat$ustar/(subsen1[,3]-subsen2[,3])*log((subsen2[,1]-dsub)/(subsen1[,1]-dsub))
      ans <-  data.frame(kh = c(kh415, kh416), d = dsub, dlab = dlab[i], z = subdz$zeval, frm = rep(c("4.15","4.16"), each = length(kh415)))
      khdf <- rbind(khdf, ans)
    }
  }
  kvdf$kv[abs(kvdf$kv) > 0.025] <- NA
  kvdf$kv[stb != "NTR"] <- NA
  khdf$kh[abs(khdf$kh) > 1] <- NA
  kvdf$kv[stb != "NTR"] <- NA
  return(list(kvdf, khdf,kdf))
}
