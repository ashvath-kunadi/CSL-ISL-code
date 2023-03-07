ozfluxdataextractrfyanc <- function(fileoz, textdestination, sttim, endtim, z0deter){
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
  Rn <- ncvar_get(nc_data, "Fn")
  Tair <- ncvar_get(nc_data, "Ta")
  rh <- ncvar_get(nc_data, "RH")
  P <- ncvar_get(nc_data, "ps")
  u <- ncvar_get(nc_data, "Ws_CSAT")
  Wd <- ncvar_get(nc_data, "Wd")
  ea <- ncvar_get(nc_data, "e")
  es <- ncvar_get(nc_data, "esat")
  G <- ncvar_get(nc_data, "Fg")
  rhod <- ncvar_get(nc_data, "rhod")
  rhom <- ncvar_get(nc_data, "rhom")
  FLh <- ncvar_get(nc_data, "Fe")
  Pr <- ncvar_get(nc_data, "Precip")
  Cpm <- ncvar_get(nc_data, "Cpm")
  Ts <- ncvar_get(nc_data, "Ts")
  wsn <- names(nc_data$var)[grepl("Ws_[0123456789]{1,2}m$", names(nc_data$var))]
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
        if(i == 1){
          Ta <- Ta1
        }else{
          Ta <- cbind(Ta, Ta1)
        }
      }
      for(i in  1:length(Ahn)){
        Ah1 <- ncvar_get(nc_data, Ahn[i])
        if(i == 1){
          Ah <- Ah1
        }else{
          Ah <- cbind(Ah, Ah1)
        }
      }
      wsp <- melt(data.frame(ws))
      variable <- gsub("Ws_",'',wsn)
      variable <- gsub("m", '', variable)
      variable <- as.integer(variable)
      wsp$variable <- rep(variable, each = length(ws1))
      Tap <- melt(data.frame(Ta))
      variable <- gsub("Ta_HMP_",'',Tan)
      variable <- gsub("m", '', variable)
      variable <- as.integer(variable)
      Tap$variable <- rep(variable, each = length(Ta1))
      Ahp <- melt(data.frame(Ah))
      variable <- gsub("Ah_HMP_",'',Ahn)
      variable <- gsub("m", '', variable)
      variable <- as.integer(variable)
      Ahp$variable <- rep(variable, each = length(Ah1))
      names(wsp) <- c("H", "WS")
      names(Tap) <- c("H", "Ta")
      names(Ahp) <- c("H", "Ah")
      #to show Richie
      
      Ahp$e <- 461.52*(Ahp$Ah/1000)*(Tap$Ta+273.15)/1000
      Ahp$P <- rep(P,length(unique(Ahp$H)))+(rep(rhod,length(unique(Ahp$H)))*9.80665*(15-Ahp$H)/1000)
      Ahp$q <- 0.622*Ahp$e*Ahp$P/(Ahp$P-(0.378*Ahp$e))/1000
      Ahp$rhom <- Ahp$P*1000/287.04/(Tap$Ta+273.15)*(1-Ahp$e/Ahp$P)
      Tap$Tavir1 <- Tap$Ta*(1+(Ahp$q*0.61))
      dtvdzrec <- 0
      for(i in 1:length(Ta1)){
          dz <- Tap[c(i, (i+length(Ta1)), (i+2*length(Ta1)),(i+3*length(Ta1))), 1]
          dTv <- Tap[c(i, (i+length(Ta1)), (i+2*length(Ta1)),(i+3*length(Ta1))), 2]
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
    
  }else{
    wsp <- NA
    Tap <- NA
    Ahp <- NA 
    dtvdzrec <- NA
  }
  rhom <- ncvar_get(nc_data, "rhom")
  rhod <- ncvar_get(nc_data, "rhod")
  H <- ncvar_get(nc_data, "Fh")
  monobhl <- ncvar_get(nc_data, "L")
  surtemp <- ncvar_get(nc_data, "Ts_3cm")
  
  ##Potential ET calculations
  psyc <- 0.665*10^-3*P #kPa/C
  D <- 4098*(0.6108*exp(17.27*Tair/(Tair+237.3)))/(Tair+237.3)^2 ##kPa/C
  dt <- 48
  ETvar <- data.frame(u, ustar,Tair, H, time, D, Rn, dt, psyc, Lheat = 2458280,  es, ea, rhom, Cpm, monobhl,  G, FLh)
  nc_close(nc_data)
  return(list(time, ETvar, Pr, wsp, Tap, Ahp, dtvdzrec, surtemp,Wd))
}