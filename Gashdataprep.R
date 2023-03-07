prepdat <- function(evtim, impst, ET0){
    fnflrf <- impst[[7]]
    fnrf <- impst[[6]]
    fnthrfN <- impst[[8]]
    fnthrfS <- impst[[9]]
    stt <- evtim[[1]]
    ent <- evtim[[2]]
    sid <- unique(fnflrf$i)
    td <- evtim[[3]]
    td <- td[sid] ##in hours hopefully
    ##getting event and average rainfall
    fltot <- tapply(fnflrf$mm, fnflrf$i, sum)
    rftot <- tapply(fnrf$mm, list(fnrf$i, fnrf$name), sum)
    rfplfl <- cbind(rftot, fltot)
    Pg <- apply(rfplfl, 1, FUN = function(x) mean(x, na.rm = T)) #in mm
    Rbar <- mean(Pg/td)
    dfrfplfl <- data.frame(rfplfl, ID = sid)
    ##getting average evaporation rate
    ETev <- 0
    for(i in 1:length(unique(stt))){
        subev <- ET0[(ET0$time > stt[i]) & (ET0$time < ent[i]) ,]
        tdev <- subev$time[length(subev$time)]-subev$time[1]
        totev <- sum(subev$mm)
        ETev[i] <- totev/as.numeric(tdev)
    }
    ETsel <- ETev[sid]
    Ebar <- mean(ETsel, na.rm = T)
    
    ##getting reference interception
    totthrfN <- tapply(fnthrfN$mm, list(fnthrfN$i, fnthrfN$name), sum)
    meanthrfN <- apply(totthrfN, 1, FUN = function(x) mean(x, na.rm=T))
    intNpl <- cbind(totthrfN, meanthrfN)
    dfintNpl <- data.frame(intNpl, ID = as.integer(rownames(intNpl)))
    totthrfS <- tapply(fnthrfS$mm, list(fnthrfS$i, fnthrfS$name), sum)
    meanthrfS <- apply(totthrfS, 1, FUN = function(x) mean(x, na.rm=T))
    intSpl <- cbind(totthrfS, meanthrfS)
    dfintSpl <- data.frame(intSpl, ID = as.integer(rownames(intSpl)))
    tryN <- merge(dfintNpl, dfrfplfl, by = "ID", all = T)
    tryS <- merge(dfintSpl, dfrfplfl, by = "ID", all = T)
    totintN <- Pg-tryN[,2:7] 
    totintS <- Pg-tryS[,2:7] 
    ##incase you want total interception loss
    INL <- apply(totintN, 2, FUN = function(x) sum(x,na.rm = T))
    ISL <- apply(totintS, 2, FUN = function(x) sum(x,na.rm = T))
    return(list(Pg, Rbar, Ebar, tryN, tryS, totintN, totintS, INL, ISL))
}
