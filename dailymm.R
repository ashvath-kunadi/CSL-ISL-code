dailymm <- function(dat){
    prodat <- data.frame()
    if(sum(names(dat) == "date") != 1){
      dat$date <- as.Date(dat$time)
    }
    for(i in 1:length(unique(dat$name))){
      subdat <- dat[unique(dat$name)[i] == dat$name, ]
      tp <- tapply(subdat$mm, subdat$date, FUN = function(x) sum(x, na.rm = T))
      date <- names(tp)
      prodat <- rbind(prodat, cbind(name = as.character(unique(dat$name)[i]), date, tp))
    }
    prodat <- as.data.frame(prodat)
    prodat$date <- as.Date(as.character(prodat$date))
    prodat$tp <- as.numeric(as.character(prodat$tp))
    return(prodat)  
}
