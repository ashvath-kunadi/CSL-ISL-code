strmsepev <- function(rf, lagtimes = 2*60*60){
    library(lubridate)
    crf <- rf[order(rf$time),-1]
    ##finding the beginning and ending time of rainfall events
    tmdif <- diff(crf$time)
    stormID <- 0
    i <- 1
    strmbgntim <- NA
    strmendtim <- NA
    check <- T
    crf <- crf[crf$mm > 0,]
    cond <- tmdif<(1.5*lagtimes)
    while(i <= dim(crf)[1]){
        stormID <- stormID+1
        check <- T
        if(i < dim(crf)[1]){
            while(cond[i]){
                if(check){
                  strmbgntim[stormID] <- crf[i,1]
                  check <- F
                }
                i <- i+1
                if(i > length(cond)){
                  break()
                }
            }
            strmendtim[stormID] <- crf[i,1]
            if(is.na(strmbgntim[stormID])){
                strmbgntim[stormID] <- crf[i,1]
            }
            strmbgntim[stormID] <- strmbgntim[stormID] - (lagtimes/2)
            strmendtim[stormID] <- strmendtim[stormID] + lagtimes
        }else{
            if(tmdif[i-1]>(1.5*lagtimes)){
              strmbgntim[stormID] <- crf[i,1] - (lagtimes/2)
              strmendtim[stormID] <- crf[i,1] + lagtimes
            }else{
              strmendtim[stormID-1] <- crf[i,1] + lagtimes
            }
        }
        i <- i+1
    }
    strmbgntim <- as.POSIXct(strmbgntim, tz = "GMT", origin = "1970-01-01")
    strmendtim <- as.POSIXct(strmendtim, tz = "GMT", origin = "1970-01-01")
    timdiff <- as.numeric(strmendtim-strmbgntim-(3/2*lagtimes/60/60))
    return(list(strmbgntim, strmendtim,timdiff))
}
