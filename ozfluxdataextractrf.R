ozfluxdataextractrf <- function(fileoz, textdestination, sttim, endtim, z0deter){
  library(ncdf4)
  library(lubridate)
  library(reshape2)
  nc_data <- nc_open(fileoz)
  if(sum(textdestination %in% list.files())!=1){
    sink(textdestination)
    print(nc_data)
    sink()
  }
  time <- paste0(ncvar_get(nc_data, "Year"), "-", ncvar_get(nc_data, "Month"),  "-", ncvar_get(nc_data, "Day"), " ", ncvar_get(nc_data, "Hour"), "-", ncvar_get(nc_data, "Minute"))
  time <- ymd_hm(time)
  Rs <- ncvar_get(nc_data, "Fld")+ncvar_get(nc_data, "Fsd")
  Rs <- aclnrintp(Rs, 0.7)
  Rn <- ncvar_get(nc_data, "Fn")
  Rn <- aclnrintp(Rn, 0.7)
  Rn_fil <- ncvar_get(nc_data, "Fn_4cmpt")
  if("Ta_HMP_15m" %in% names(nc_data$var)){
    Tair <- ncvar_get(nc_data, "Ta_HMP_15m")
  }else{
    Tair <- ncvar_get(nc_data, "Ta")
  }
  ta15 <- ncvar_get(nc_data, "Ta_HMP_15m")
  Tair <- aclnrintp(Tair, 0.7)
  Tair[is.na(Tair)] <- ta15[is.na(Tair)]
  rh15 <- ncvar_get(nc_data, "RH")
  P <- ncvar_get(nc_data, "ps")
  P <- aclnrintp(P, 0.7)
  u <- ncvar_get(nc_data, "Ws")
  u <- aclnrintp(u, 0.7)
  ea <- ncvar_get(nc_data, "e")
  ea <- aclnrintp(ea, 0.7)
  es <- ncvar_get(nc_data, "esat")
  es <- aclnrintp(es, 0.7)
  Tair[is.na(Tair)] <- 237.3*log(es[is.na(Tair)]/.6108)/(17.27-log(es[is.na(Tair)]/.6108))
  es[is.na(es)] <- 0.6108*exp(17.27*Tair[is.na(es)]/(237.3+Tair[is.na(es)]))
  G <- ncvar_get(nc_data, "Fg")
  G <- aclnrintp(G, 0.7)
  FLh <- ncvar_get(nc_data, "Fe")
  rhom <- ncvar_get(nc_data, "rhom")
  rhod <- ncvar_get(nc_data, "rhod")
  Pr <- ncvar_get(nc_data, "Precip")
  Cpm <- ncvar_get(nc_data, "Cpm")
  Cpd <- ncvar_get(nc_data, "Cpd")
  if("Ts_4cm" %in% names(nc_data$var)){
    Ts1 <- ncvar_get(nc_data, "Ts_4cm")
  }else if("Ts_10cm" %in% names(nc_data$var)){
    Ts1 <- ncvar_get(nc_data, "Ts_10cm")
  }else{
    Ts1 <- ncvar_get(nc_data, "Ts")
  }
  wsn <- names(nc_data$var)[grepl("Ws_WS4_[0123456789]{1,2}m$", names(nc_data$var))]
  Tan <- names(nc_data$var)[grepl("Ta_HMP_[0123456789]{1,2}m$", names(nc_data$var))]
  Ahn <- names(nc_data$var)[grepl("Ah_HMP_[0123456789]{1,2}m$", names(nc_data$var))]
  #plot(Ts, Ts4)
  #abline(coef = c(0,1))
  ##assuming that the temperature and velocity sensors are placed at the same height
  if(length(wsn) >1){
      for(i in  1:length(wsn)){
        ws1 <- ncvar_get(nc_data, wsn[i])
        if(i == 1){
          ws <- ws1
        }else{
          ws <- cbind(ws, ws1)
        }
      }
      for(i in  1:length(Tan)){
        Ta1 <- ncvar_get(nc_data, Tan[i])
        Ah1 <- ncvar_get(nc_data, Ahn[i])
        if(i == 1){
          Ta <- Ta1
          Ah <- Ah1
        }else{
          Ta <- cbind(Ta, Ta1)
          Ah <- cbind(Ah, Ah1)
        }
      }
      wsp <- melt(data.frame(ws))
      variable <- gsub("Ws_WS4_",'',wsn)
      variable <- gsub("m", '', variable)
      variable <- as.integer(variable)
      variable[variable == 10] <- 10.6
      variable[variable == 16] <- 15.5
      variable[variable == 2] <- 1.9
      variable[variable == 7] <- 7
      wsp$variable <- rep(variable, each = length(ws1))
      Tap <- melt(data.frame(Ta))
      Ahp <- melt(data.frame(Ah))
      variable <- gsub("Ta_HMP_",'',Tan)
      variable <- gsub("m", '', variable)
      variable <- as.integer(variable)
      variable[variable == 10] <- 10.6
      variable[variable == 15] <- 14.7
      variable[variable == 2] <- 1.25
      variable[variable == 4] <- 4.2
      variable[variable == 7] <- 7.4
      Tap$variable <- rep(variable, each = length(Ta1))
      Ahp$variable <- rep(variable, each = length(Ah1))
      names(wsp) <- c("H", "WS")
      names(Tap) <- c("H", "Ta")
      names(Ahp) <- c("H", "Ah")
      #to show Richie
      ##Make the virtual temperature kelvin
      Ahp$e <- 461.52*(Ahp$Ah/1000)*(Tap$Ta+273.15)/1000
      Ahp$P <- rep(P,5)+(rep(rhod,5)*9.80665*(15-Ahp$H)/1000)
      Ahp$q <- 0.622*Ahp$e*Ahp$P/(Ahp$P-(0.378*Ahp$e))/1000
      Ahp$rhom <- Ahp$P*1000/287.04/(Tap$Ta+273.15)*(1-Ahp$e/Ahp$P)
      Tap$Tavir1 <- ((Tap$Ta+273.15)*(1+(Ahp$q*0.61)))-273.15
      # Tap$Tapot <- Tap$Ta*(1000/rep(P, 5))^(287.04*(1-(0.23*Ahp$q1))/rep(Cpd, 5))
      # Tap$Tapotvir <- Tap$Tapot*(1+(Ahp$q*0.61))
      # Tap$Tavirpot <- Tap$Tavir*(1000/rep(P, 5))^(287.04*(1-(0.23*Ahp$q))/rep(Cpd, 5))
      adilpsrt <- round((9.80665/Cpd),3)
      dtvdzrec <- 0
      for(i in 1:length(Ta1)){
          dz <- Tap[c(i, (i+length(Ta1)), (i+2*length(Ta1)),(i+3*length(Ta1))), 1]
          dTv <- Tap[c(i, (i+length(Ta1)), (i+2*length(Ta1)),(i+3*length(Ta1))), 3]
          if(sum(is.na(dTv))<2){
            dtvdz <- lm(dTv~dz)
            dtvdzrec[i] <- round(dtvdz$coefficients[2], 4)
          }else{
            dtvdzrec[i] <- NA
          }
      }
      wsp$time <- time
      Tap$time <- time
      Ahp$time <- time
      ustar <- ncvar_get(nc_data, "ustar")
      wsp$ustar <- ustar
    
    ##wind direction
    wdn <- names(nc_data$var)[grepl("Wd_WS4_[0123456789]{1,2}m$", names(nc_data$var))]
    for(i in  1:length(wdn)){
      wd1 <- ncvar_get(nc_data, wdn[i])
      if(i == 1){
        wd <- wd1
      }else{
        wd <- cbind(wd, wd1)
      }
    }
    wdp <- melt(data.frame(wd))
    variable <- gsub("Wd_WS4_",'',wdn)
    variable <- gsub("m", '', variable)
    variable <- as.integer(variable)
    variable[variable == 10] <- 10.6
    variable[variable == 16] <- 15.5
    variable[variable == 2] <- 1.9
    variable[variable == 7] <- 7
    wdp$variable <- rep(variable, each = length(wd1))
    names(wdp) <- c("H", "WD")
    wdp$time <- time 
  }else{
    wsp <- NA
    wdp <- NA
    Tap <- NA
    Ahp <- NA 
    dtvdzrec <- NA
  }
  ustar <- ncvar_get(nc_data, "ustar")
  H <- ncvar_get(nc_data, "Fh")
  monobhl <- ncvar_get(nc_data, "L")
  if("IR_TargT" %in% names(nc_data$var)){
    IRtemp <- ncvar_get(nc_data, "IR_TargT")
  }else{
    IRtemp <- NA
  }
  ##Potential ET calculations
  Lheat <- ncvar_get(nc_data, "Lv")
  rhow <- 1000 #kg/m3
  psyc <- 0.665*10^-3*P #kPa/C
  D <- 4098*(0.6108*exp(17.27*Tair/(Tair+237.3)))/(Tair+237.3)^2 ##kPa/C
  dt <- 48
  ETvar <-as.data.frame(cbind( u, ustar,Tair, H, time, D, Rn, dt, psyc, rhow, Lheat,  es, ea, rhom, Cpm, monobhl,  G, FLh, IRtemp))
  
  ##estimating the soil moisture profile
  soilwaternames <- names(nc_data$var)[grepl("Sws_[0123456789]{2,3}cm$", names(nc_data$var))]
  for(i in 1:length(soilwaternames)){
    d1 <- ncvar_get(nc_data, soilwaternames[i])
    if(i == 1){
      d <- d1
    }else{
      d <- cbind(d, d1)
    }
  }
  smp <- melt(as.data.frame(d))
  variable <- gsub("Sws_", '', soilwaternames)
  variable <- gsub("cm", '', variable)
  variable <- as.integer(variable)
  smp$variable <- rep(variable, each = length(d1))
  names(smp) <- c("D", "Theta")
  smp$time <- time
  
  soiltempnames <- names(nc_data$var)[grepl("Ts_[0123456789]{2,3}cm$", names(nc_data$var))]
  for(i in 1:length(soiltempnames)){
    d1 <- ncvar_get(nc_data, soiltempnames[i])
    if(i == 1){
      d <- d1
    }else{
      d <- cbind(d, d1)
    }
  }
  stp <- melt(as.data.frame(d))
  variable <- gsub("Ts_", '', soiltempnames)
  variable <- gsub("cm", '', variable)
  variable <- as.integer(variable)
  stp$variable <- rep(variable, each = length(d1))
  names(stp) <- c("sen", "temp")
  stp$time <- time
  
  soiltempnames2 <- names(nc_data$var)[grepl("Ts_CS650_[0123456789]{2,3}cm[123]{0,1}$", names(nc_data$var))]
  if(length(soiltempnames2) != 0){
    for(i in 1:length(soiltempnames2)){
      d1 <- ncvar_get(nc_data, soiltempnames2[i])
      if(i == 1){
        d <- d1
      }else{
        d <- cbind(d, d1)
      }
    }
    stp2 <- melt(as.data.frame(d))
    variable <- gsub("Ts_CS650_", '', soiltempnames2)
    variable <- gsub("cm", '', variable)
    stp2$variable <- rep(variable, each = length(d1))
    names(stp) <- c("sen", "temp")
    stp2$time <- time
  }else{
    stp2 <- NA
  }
  
  soiltempnames3 <- names(nc_data$var)[grepl("Ts_T108_[0123456789]{2}cm[ab]{0,1}$", names(nc_data$var))]
  if(length(soiltempnames3) != 0){
    for(i in 1:length(soiltempnames3)){
      d1 <- ncvar_get(nc_data, soiltempnames3[i])
      if(i == 1){
        d <- d1
      }else{
        d <- cbind(d, d1)
      }
    }
    stp3 <- melt(as.data.frame(d))
    variable <- gsub("Ts_T108_", '', soiltempnames3)
    variable <- gsub("cm", '', variable)
    stp3$variable <- rep(variable, each = length(d1))
    names(stp) <- c("sen", "temp")
    stp3$time <- time
  }else{
    stp3 <- NA
  }
  
  surtemp <- names(nc_data$var)[grepl("Ts_4cm[abcd]{0,1}$", names(nc_data$var))]
  if(length(surtemp) != 0){
    for(i in 1:length(surtemp)){
      d1 <- ncvar_get(nc_data, surtemp[i])
      if(i == 1){
        d <- d1
      }else{
        d <- cbind(d, d1)
      }
    }
    surt <- melt(as.data.frame(d))
    variable <- gsub("Ts_4cm", '', surtemp)
    surt$variable <- rep(variable, each = length(d1))
    names(stp) <- c("sen", "temp")
    surt$time <- time
  }else{
    surt <- NA
  }
    
  stdat <- list(stp, stp2, stp3, surt)
  nc_close(nc_data)
  return(list(time, ETvar, Pr, smp, wsp, wdp, Tap, Ahp, dtvdzrec, stdat, Rs))
}
