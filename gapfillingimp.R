gapfillingimp <- function(taggeddata, evtim, ET0, timint = 60){
  library(lubridate)
  ETgpfil <- data.frame()
  for(i in 6:9){
    dat <- taggeddata[[i]]
    gpfilledrgdat <- data.frame()
    if(i == 7){
      refdat <- taggeddata[[10]]
    }
    for(j in 1:length(unique(dat$i))){
        ID <- unique(dat$i)[j]
        sub <- dat[dat$i == ID, ]
        if(i == 7){
           subref <- refdat[refdat$ID == ID,]
        }
        lowtime <- evtim[[1]][ID]
        hightime <- evtim[[2]][ID]
        # Set the minute and second to the nearest 10 minute value
        month(lowtime) <- floor(month(lowtime)/timint*60*60*24*30) * timint/60/60/24/30
        month(hightime) <- ceiling(month(hightime)/timint*60*60*24*30) * timint/60/60/24/30
        day(lowtime) <- floor(day(lowtime)/timint*60*60*24) * timint/60/60/24
        day(hightime) <- round(ceiling(day(hightime)/timint*60*60*24)/(60*60*24/timint))
        hour(lowtime) <- floor(hour(lowtime)/timint*60*60) * timint/60/60
        hour(hightime) <- ceiling(hour(hightime)/timint*60*60) * timint/60/60
        minute(lowtime) <- floor(minute(lowtime)/timint*60) * timint/60
        minute(hightime) <- ceiling(minute(hightime)/timint*60) * timint/60
        second(lowtime) <- 0
        second(hightime) <- 0
        #  Set the breakpoints at 10 minute intervals
        breakpoints <- seq.POSIXt(lowtime, hightime, by = timint)
        if(i == 6){
          output <- aggregate(mm ~ cut(time, breaks = breakpoints), sub, sum)
          names(output) <- c("time", "mm")
          output$time <- as.POSIXct(as.character(output$time), origin = "1970-01-01 00:00:00", tz = "GMT")
          output$mm <- output$mm/length(unique(sub$name))
          nodat <- data.frame(time = breakpoints, nd = 0)
          comb <- merge(nodat, output, by = "time", all=T)
          comb$mm[is.na(comb$mm)] <- 0
          comb <- comb[,-2]
          gpfilledrgdat <- rbind(gpfilledrgdat,cbind(comb, ID))
          ##Evaporation treated
          ETgpfilsub <- ET0[ET0$time>=(lowtime-(29*60)) & ET0$time<=(hightime+(29*60)),]
          combET <- merge(nodat, ETgpfilsub, by = "time", all=T)
          cond <- (as.numeric(combET$time-min(combET$time))%%1800) == 0
          x <- rle(cond)$lengths
          cx <- cumsum(x)
          lx <- rle(cond)$values
          combET$try <- combET$mm/30
          for(k in 1:(length(x)-1)){
            if(lx[k] | k == 1){
              lval <- combET$try[cx[k]]
              hval <- combET$try[cx[k+2]]
              valdif <- hval - lval
              ltim <- combET$time[cx[k]]
              htim <- combET$time[cx[k+2]]
            }else{
              combET$try[(cx[k-1]+1):(cx[k])] <- lval + as.numeric(combET$time[(cx[k-1]+1):(cx[k])]-ltim)/as.numeric(htim-ltim)*valdif
            }
          }
          combET <- combET[combET$time>=(lowtime) & combET$time<=(hightime),c(1,4)]
          names(combET) <- c("time", "mm")
          combET$i <- ID
          ETgpfil <- rbind(ETgpfil, combET)
        }else if(i == 7){
          sfl <- sum(sub$mm, na.rm = T)
          srf <- sum(subref$mm, na.rm = T)
          subref$mm <- subref$mm*sfl/srf
          gpfilledrgdat <- rbind(gpfilledrgdat, subref)
        }else{
          for(k in 1:length(unique(sub$name))){
            namdat <- sub[sub$name==unique(sub$name)[k],]
            output <- aggregate(mm ~ cut(time, breaks = breakpoints), namdat, sum)
            names(output) <- c("time", "mm")
            output$time <- as.POSIXct(as.character(output$time), origin = "1970-01-01 00:00:00", tz = "GMT")
            nodat <- data.frame(time = breakpoints, nd = 0)
            comb <- merge(nodat, output, by = "time", all=T)
            comb$mm[is.na(comb$mm)] <- 0
            comb <- comb[,-2]
            comb$name <- unique(sub$name)[k]
            gpfilledrgdat <- rbind(gpfilledrgdat, cbind(comb, ID))
          }
        }
    }
   taggeddata[[4+i]] <-  gpfilledrgdat
  }
  togive <- list(taggeddata[[10]],taggeddata[[11]], taggeddata[[12]], taggeddata[[13]], ETgpfil)
  return(togive)
}
