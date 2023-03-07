prrprsc <- function(csvfile = "H:/My Documents/Project and related work/Cool gingin/L3data/AMF_US-Prr_BASE_HH_3-5.csv", canh = 3, ntrlim = 0.01){
  library(lubridate)
  rawdat <- read.csv(csvfile , skip = 2, na.strings = -9999)
  prdat <- data.frame(time = strptime(as.character(rawdat$TIMESTAMP_START), format = "%Y%m%d%H%M", tz = "GMT"))
  prdat$ustar <- rawdat[,names(rawdat) == "USTAR_1_1_1"]
  rawdat <- rawdat[!is.na(prdat$ustar),]
  prdat <- prdat[!is.na(prdat$ustar),]
  wsn <- names(rawdat)[grepl("WS_1_[0123456789]{1,2}_1$", names(rawdat))]
  wsh <- c(16,13,11,9,7.5,6,4.5,3,1.9,1.5)
  wspprr <- data.frame()
  for(i in 1:length(wsn)){
    varnam <- paste("WS_1_", as.character(i), "_1", sep = "")
    ans <- rawdat[,names(rawdat) == varnam]
    ansdf <- cbind(WS = ans, H = wsh[i])
    if(i == 1){
      wspprr <- ansdf
    }else{
      wspprr <- rbind(wspprr, ansdf)
    }
  }
  
  
  wspprr <- as.data.frame(wspprr)
  ggplot(wspprr, aes(WS, H, group = H)) + geom_boxplot() + theme_bw() + ggtitle("Poker Flat Research Range")
  
  tan <- names(rawdat)[grepl("TA_1_[0123456789]{1,2}_1$", names(rawdat))]
  tah <- c(16,13,11,9,6,4.5,3,1.9,1.5)
  tapprr <- data.frame()
  for(i in 1:length(tan)){
    varnam <- paste("TA_1_", as.character(i), "_1", sep = "")
    ans <- rawdat[,names(rawdat) == varnam]
    ansdf <- cbind(TA = ans, H = tah[i])
    if(i == 1){
      tapprr <- ansdf
    }else{
      tapprr <- rbind(tapprr, ansdf)
    }
  }
 
  
  tapprr <- as.data.frame(tapprr)
  # ggplot(tapprr, aes(TA, H, group = H)) + geom_boxplot() + theme_bw() + ggtitle("Poker Flat Research Range")

  rhn <- names(rawdat)[grepl("RH_1_[0123456789]{1,2}_1$", names(rawdat))]
  rhh <- c(16,13,11,9,6,4.5,3,1.9,1.5)
  rhpprr <- data.frame()
  for(i in 1:length(rhn)){
    varnam <- paste("RH_1_", as.character(i), "_1", sep = "")
    ans <- rawdat[,names(rawdat) == varnam]
    ansdf <- cbind(RH = ans, H = tah[i])
    if(i == 1){
      rhpprr <- ansdf
    }else{
      rhpprr <- rbind(rhpprr, ansdf)
    }
  }
  
  
  rhpprr <- as.data.frame(rhpprr)
  # ggplot(rhpprr, aes(RH, H, group = H)) + geom_boxplot() + theme_bw() + ggtitle("Poker Flat Research Range")
  
  prdat$UH <- wspprr[wspprr$H == 3,1]
  prdat$U <- wspprr[wspprr$H == 11,1]
  prdat$WD <- rawdat[,names(rawdat) == "WD_1_1_1"]
  prdat$P <- rawdat[,names(rawdat) == "PA"]
  prdat$Tair <- tapprr[tapprr$H == 11,1]
  prdat$H <- rawdat[,names(rawdat) == "H_1_1_1"]
  prdat$Flh <- rawdat[,names(rawdat) == "LE_1_1_1"]
  # rawdat <- rawdat[complete.cases(prdat),]
  # prdat <- prdat[complete.cases(prdat),]
  prdat$x <- rawdat$H2O_1_1_1
  prdat$x[prdat$x == -9999] <- NA
  ##getting x which is the specific humidity
  prdat$x <- (prdat$x*18.02/1000)/((1-(prdat$x/1000))*28.97)
  prdat$rhom <- with(prdat, (P*1000/
                               286.9/
                               (Tair+273.15)*
                               (1+x))/
                       (1+(x*286.9/461.5)))
  prdat$Cpm <- 1005 + 1820*prdat$x
  prdat$Lheat <- (2500.8-(2.36*(prdat$Tair))+(0.0016*(prdat$Tair^2))-(0.00006*(prdat$Tair^3)))*1000
  prdat$L <- with(prdat, -(ustar^3*rhom)/(0.4*9.80665*((H/(Tair+273)/Cpm)+(0.61*Flh/Lheat))))
  ##initial estimate of d0
  cd1 <- 2/3*canh
  cd2 <- cd1-0.2
  count <- 1
  
  wspprr$time <- prdat$time
  wspprr$ustar <- prdat$ustar
  
  while(abs(cd1-cd2) > 0.01 & count < 10){
    zeta1 <- (15-cd1)/prdat$L
    fkn2017cat <- rep("VUST", length(zeta1))
    fkn2017cat[zeta1 >= -1] <- "UST"
    fkn2017cat[zeta1 >= -ntrlim] <- "NTR"
    fkn2017cat[zeta1 >= ntrlim] <- "STB"
    fkn2017cat[zeta1 >= 1] <- "VSTB"
    fkn2017cat[zeta1 >= 7] <- NA
    fkn2017cat[zeta1 <= -5] <- NA
    ##new stability based on L
    ntr <- prdat$time[fkn2017cat == "NTR"]
    ntrcond <- (fkn2017cat == "NTR")
    est <- z0ndnlsnslp(wspprr, ntr, prdat$time, 3)
    ##first estimate stull 382 wind speed
    cd1 <- cd2
    vard2 <- est[[2]]
    cd2 <- mean(vard2, na.rm = T)
    count <- count+1
  }
  gngpup <- data.frame(d0 = vard2, z0 = est[[1]], err = est[[3]], snsmx = 16)
  for(i in unique(wspprr$H)[unique(wspprr$H)>4.5 & unique(wspprr$H)<16]){
    cd1rf <- 2/3*canh
    cd2rf <- cd1rf-0.2
    count <- 1
    rfwspr <- wspprr[wspprr$H <= i,]
    while(abs(cd1rf-cd2rf) > 0.01 & count < 10){
      zeta1 <- (15-cd1)/prdat$L
      fkn2017cat <- rep("VUST", length(zeta1))
      fkn2017cat[zeta1 >= -1] <- "UST"
      fkn2017cat[zeta1 >= -ntrlim] <- "NTR"
      fkn2017cat[zeta1 >= ntrlim] <- "STB"
      fkn2017cat[zeta1 >= 1] <- "VSTB"
      fkn2017cat[zeta1 >= 7] <- NA
      fkn2017cat[zeta1 <= -5] <- NA
      ##new stability based on L
      ntr <- prdat$time[fkn2017cat == "NTR"]
      ntrcond <- (fkn2017cat == "NTR")
      estrf <- z0ndnlsnslp(rfwspr, ntr, prdat$time, 3)
      ##first estimate stull 382 wind speed
      cd1rf <- cd2rf
      vard2rf <- estrf[[2]]
      cd2rf <- mean(vard2rf, na.rm = T)
      count <- count+1
    }
    ans <- data.frame(d0 = vard2rf, z0 = estrf[[1]], err = estrf[[3]], snsmx = i)
    gngpup <- rbind(gngpup, ans)
  }
  gngpup <- cbind(gngpup, data.frame(d0ref = vard2, z0ref = est[[1]], errref = est[[3]]))
  z0lmwht <- tapply(gngpup$z0, gngpup$snsmx, FUN = function(x) lm(x~est[[1]]))
  d0lmwht <- tapply(gngpup$d0, gngpup$snsmx, FUN = function(x) lm(x~vard2))
  gngpdn <- data.frame(d0 = vard2, z0 = est[[1]], err = est[[3]], snsmn = 3)
  for(i in unique(wspprr$H)[unique(wspprr$H)>3 & unique(wspprr$H)<13]){
    cd1isl <- 2/3*canh
    cd2isl <- cd1isl-0.2
    count <- 1
    islwspr <- wspprr[wspprr$H >= i,]
    while(abs(cd1isl-cd2isl) > 0.01 & count < 10){
      zeta1 <- (15-cd1)/prdat$L
      fkn2017cat <- rep("VUST", length(zeta1))
      fkn2017cat[zeta1 >= -1] <- "UST"
      fkn2017cat[zeta1 >= -ntrlim] <- "NTR"
      fkn2017cat[zeta1 >= ntrlim] <- "STB"
      fkn2017cat[zeta1 >= 1] <- "VSTB"
      fkn2017cat[zeta1 >= 7] <- NA
      fkn2017cat[zeta1 <= -5] <- NA
      ##new stability based on L
      ntr <- prdat$time[fkn2017cat == "NTR"]
      ntrcond <- (fkn2017cat == "NTR")
      estisl <- z0ndnlsnslp(islwspr, ntr, prdat$time, 3)
      ##first estimate stull 382 wind speed
      cd1isl <- cd2isl
      vard2isl <- estisl[[2]]
      cd2isl <- mean(vard2isl, na.rm = T)
      count <- count+1
    }
    ans <- data.frame(d0 = vard2isl, z0 = estisl[[1]], err = estisl[[3]], snsmn = i)
    gngpdn <- rbind(gngpdn, ans)
  }
  
  z0lmwht2 <- tapply(gngpdn$z0, gngpdn$snsmn, FUN = function(x) lm(x~est[[1]]))
  d0lmwht2 <- tapply(gngpdn$d0, gngpdn$snsmn, FUN = function(x) lm(x~vard2))
  
  prdat$d0 <- vard2/canh
  prdat$beta <- with(prdat, ustar/UH)
  prdat$mnth <- month(prdat$time)
  prdat$season <- "winter"
  prdat$season[(prdat$mnth >= 5 & prdat$mnth < 11)] <- "summer"
  prdat$z0 <- est[[1]]/canh
  return(list(prdat, wspprr, gngpdn, gngpup, z0lmwht, z0lmwht2, d0lmwht, d0lmwht2, fkn2017cat))
}
