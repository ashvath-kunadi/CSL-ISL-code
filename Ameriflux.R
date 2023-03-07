amriflxprsc <- function(csvfile, canh, ws1z, ws2z, ntrlim, wscan = NA){
  library(lubridate)
  library(dplyr)
  rawdat <- read.csv(csvfile , skip = 2, na.strings = -9999)
  prdat <- data.frame(time = strptime(as.character(rawdat$TIMESTAMP_START), format = "%Y%m%d%H%M", tz = "GMT"))
  if("USTAR" %in% names(rawdat)){
    prdat$ustar <- rawdat[,names(rawdat) == "USTAR"]
  }else{
    prdat$ustar <- rawdat[,names(rawdat) == "USTAR_PI_F_1_1_1"]
  }
  rawdat <- rawdat[!is.na(prdat$ustar),]
  prdat <- prdat[!is.na(prdat$ustar),]
  prdat$UH <- rawdat[,names(rawdat) == "WS_1_2_1"]
  prdat$U <- rawdat[,names(rawdat) == "WS_1_1_1"]
  prdat$WD <- rawdat[,names(rawdat) == "WD_1_1_1"]
  
  if("PA" %in% names(rawdat)){
    prdat$P <- rawdat[,names(rawdat) == "PA"]
  }else{
    prdat$P <- rawdat[,names(rawdat) == "PA_1_1_1"]
  }
  prdat$Tair <- rawdat[,names(rawdat) == "TA_1_1_1"]
  if("H" %in% names(rawdat)){
    prdat$H <- rawdat$H
  }else{
    prdat$H <- rawdat[,names(rawdat) == "H_1_1_1"]
  }  
  if("LE" %in% names(rawdat)){
    prdat$Flh <- rawdat$LE
  }else{
    prdat$Flh <- rawdat[,names(rawdat) == "LE_1_1_1"]
  }
  rawdat <- rawdat[complete.cases(prdat),]
  prdat <- prdat[complete.cases(prdat),]
  if(is.na(wscan)){
    wscanchck <- NA
  }else{
    wscanchck <- data.frame(time = prdat$time, UH = prdat$UH, ustar = prdat$ustar, zcan = wscan[1], uz = rawdat[,names(rawdat) == "WS_1_3_1"])
  }
  if(length(wscan)>1){
    ans <- wscanchck
    for(i in 2:length(wscan)){
      ans$uz <- rawdat[,names(rawdat) == paste("WS_1_",(i+2),"_1", sep = "")]
      ans$zcan <- wscan[2]
      wscanchck <- rbind(wscanchck, ans)
    }
  }
  for(i in 1:ncol(prdat)){
    cond <- (prdat[,i] == -9999)
    prdat[cond,i] <- NA
  }
  
  if("H2O" %in% names(rawdat)){
    prdat$x <- rawdat$H2O
    prdat$x[prdat$x == -9999] <- NA
    ##getting x which is the specific humidity
    prdat$x <- (prdat$x*18.02/1000)/((1-(prdat$x/1000))*28.97)
    prdat$rhom <- with(prdat, (P*1000/286.9/(Tair+273.15)*(1+x))/(1+(x*286.9/461.5)))
  }else if("H2O_1_1_1" %in% names(rawdat)){
    prdat$x <- rawdat$H2O_1_1_1
    prdat$x[prdat$x == -9999] <- NA
    ##getting x which is the specific humidity
    prdat$x <- (prdat$x*18.02/1000)/((1-(prdat$x/1000))*28.97)
    prdat$rhom <- with(prdat, (P*1000/286.9/(Tair+273.15)*(1+x))/(1+(x*286.9/461.5)))
  }else{
    prdat$rhom <- NA
  }
  prdat$Cpm <- 1005 + 1820*prdat$x
  prdat$Lheat <- (2500.8-(2.36*(prdat$Tair))+(0.0016*(prdat$Tair^2))-(0.00006*(prdat$Tair^3)))*1000
  if("MO_LENGTH" %in% names(rawdat)){
    prdat$L <- rawdat$MO_LENGTH
  }else{
    prdat$L <- with(prdat, -(ustar^3*rhom)/(0.4*9.80665*((H/(Tair+273)/Cpm)+(0.61*Flh/Lheat))))
  }
  wscanchck <- wscanchck[complete.cases(prdat),]
  prdat <- prdat[complete.cases(prdat),]
  ##initial estimate of d0
  cd1 <- 2/3*canh
  cd2 <- cd1-0.2
  count <- 1
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
    
    ##first estimate stull 382 wind speed
    d0 <- 0
    zeval <- mean(c(ws1z,ws2z))
    limit <- NA
    for(i in 1:length(prdat$time)){
      if(prdat$time[i] %in% ntr){
        d0[i] <- with(prdat, zeval-(ustar[i]/0.4/((U[i]-UH[i])/(ws1z-ws2z))))
        if(d0[i]>=canh){
          limit[i] <- "Upper"
        }else if(d0[i]<=0){
          limit[i] <- "Lower"
        }else{
          limit[i] <- "Good"
        }
      }
      else{
        d0[i] <- NA
        limit[i] <- NA
      }

    }
    cd1 <- cd2
    vard2 <- d0
    cd2 <- mean(vard2, na.rm = T)
    count <- count+1
  }
  prdat$d0 <- d0/canh
  prdat$beta <- with(prdat, ustar/UH)
  prdat$mnth <- month(prdat$time)
  prdat$season <- "winter"
  prdat$season[(prdat$mnth >= 5 & prdat$mnth < 11)] <- "summer"
  return(list(prdat, wscanchck, limit))
}
