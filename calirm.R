calirm <- function(rf, thrfN, thrfS, refrg, dates, refthrsh = 1, calitips = 50, mrntime = 9, evntime = 20, timbwtips = 60, thresh = 95, img = F){
    library(ggplot2)
    library(lubridate)
    rmdat <- dmy(as.character(dates$calidates))
    rf <- rf[!(as.Date(rf$time) %in% rmdat),]
    thrfN <- thrfN[!(as.Date(thrfN$time) %in% rmdat),]
    thrfS <- thrfS[!(as.Date(thrfS$time) %in% rmdat),]
    refrg <- refrg[!(as.Date(refrg$time) %in% rmdat),]
    if(img){
        for(i in 1:length(unique(rf$name))){
            sub <- rf[rf$name == unique(rf$name)[i],]
            jpeg(file = paste(unique(rf$name)[i], '.jpeg', sep = ''))
            g1 <- ggplot(sub, aes(time, cumsum(mm)))+geom_point()
            g1 <- g1 + theme_minimal() + xlab("Time") + ylab("Cumilative Rainfall in mm")
            print(g1)
            dev.off()
        }
        for(i in 1:length(unique(thrfN$name))){
          sub <- thrfN[thrfN$name == unique(thrfN$name)[i],]
          jpeg(file = paste(unique(thrfN$name)[i], '.jpeg', sep = ''))
          g1 <- ggplot(sub, aes(time, cumsum(mm)))+geom_point()
          g1 <- g1 + theme_minimal() + xlab("Time") + ylab("Cumilative Rainfall in mm")
          print(g1)
          dev.off()
        }
        for(i in 1:length(unique(thrfS$name))){
          sub <- thrfS[thrfS$name == unique(thrfS$name)[i],]
          jpeg(file = paste(unique(thrfS$name)[i], '.jpeg', sep = ''))
          g1 <- ggplot(sub, aes(time, cumsum(mm)))+geom_point()
          g1 <- g1 + theme_minimal() + xlab("Time") + ylab("Cumilative Rainfall in mm")
          print(g1)
          dev.off()
        }
    }
    drefrg <- dailymm(refrg)
    refdate <- unique(drefrg$date[drefrg$tp > refthrsh])
    dthrfN <- dailymm(thrfN)
    thrfdateN <- unique(dthrfN$date[dthrfN$tp > 0])
    dthrfS <- dailymm(thrfS)
    thrfdateS <- unique(dthrfS$date[dthrfS$tp > 0])
    dtrf <- dailymm(rf)
    rfdate <- unique(dtrf$date[dtrf$tp > 0])
    
    funkydates <- unique(c(thrfdateN[!(thrfdateN %in% refdate)], thrfdateS[!(thrfdateS %in% refdate)], rfdate[!(rfdate %in% refdate)]))
    
    for(i in 1:length(funkydates)){
      for(j in 1:length(unique(rf$name))){
        sub <- rf[rf$name == unique(rf$name)[j],]
        sub <- sub[as.Date(sub$time) == funkydates[i],]
        sub <- sub[order(sub$time),]
        cond1 <- dim(sub)[1] > calitips
        if(cond1){
          if(hour(median(sub$time)) > mrntime & hour(median(sub$time)) < evntime){
            if(img){
                jpeg(file = paste(unique(rf$name)[j], as.character(funkydates[i]), '.jpeg', sep = ''))
                g1 <- ggplot(sub, aes(time, cumsum(mm)))+geom_point()
                g1 <- g1 + theme_minimal() + xlab("Time") + ylab("Cumilative Rainfall in mm") + ggtitle(as.character(funkydates[i]))
                print(g1)
                dev.off()
                td <- as.numeric(diff(sub$time))
                jpeg(file = paste(unique(rf$name)[j], as.character(funkydates[i]), 'hist.jpeg', sep = ''))
                r <- hist(td, breaks = c(timbwtips*0:3, max(td)), xlim = range(0, min(1200,max(td))), col = "blue1")
                text(r$mids, r$density, r$counts, adj = c(.5, -.5), col = "blue3")
                lines(r, lty = 3, border = "purple") # -> lines.histogram(*)
                dev.off()
            }else{
              td <- as.numeric(diff(sub$time))
              r <- hist(td, breaks = c(timbwtips*0:3, max(td)), xlim = range(0, min(1200,max(td))), col = "blue1")
            }
            if((r$counts/sum(r$counts))[1] > thresh/100){
              rf <- rf[!(rf$time %in% sub$time),]
            }
          }
        }
      }
      for(j in 1:length(unique(thrfN$name))){
        sub <- thrfN[thrfN$name == unique(thrfN$name)[j],]
        sub <- sub[as.Date(sub$time) == funkydates[i],]
        sub <- sub[order(sub$time),]
        cond1 <- dim(sub)[1] > calitips
        if(cond1){
          if(hour(median(sub$time)) > mrntime & hour(median(sub$time)) < evntime){
            if(img){
                jpeg(file = paste(unique(thrfN$name)[j], as.character(funkydates[i]), '.jpeg', sep = ''))
                g1 <- ggplot(sub, aes(time, cumsum(mm)))+geom_point()
                g1 <- g1 + theme_minimal() + xlab("Time") + ylab("Cumilative Throughfall in mm")
                print(g1)
                dev.off()
                td <- as.numeric(diff(sub$time))
                jpeg(file = paste(unique(thrfN$name)[j], as.character(funkydates[i]), 'hist.jpeg', sep = ''))
                r <- hist(td, breaks = c(timbwtips*0:3, max(td)), xlim = range(0, min(1200,max(td))), col = "blue1")
                text(r$mids, r$density, r$counts, adj = c(.5, -.5), col = "blue3")
                lines(r, lty = 3, border = "purple") # -> lines.histogram(*)
                dev.off()
            }else{
              td <- as.numeric(diff(sub$time))
              r <- hist(td, breaks = c(timbwtips*0:3, max(td)), xlim = range(0, min(1200,max(td))), col = "blue1")
            }
            if((r$counts/sum(r$counts))[1] > thresh/100){
              thrfN <- thrfN[!(thrfN$time %in% sub$time),]
            }
          }
        }
      }
      for(j in 1:length(unique(thrfS$name))){
        sub <- thrfS[thrfS$name == unique(thrfS$name)[j],]
        sub <- sub[as.Date(sub$time) == funkydates[i],]
        sub <- sub[order(sub$time),]
        cond1 <- dim(sub)[1] > calitips
        if(cond1){
          if(hour(median(sub$time)) > mrntime & hour(median(sub$time)) < evntime){
            if(img){
                jpeg(file = paste(unique(thrfS$name)[j], as.character(funkydates[i]), '.jpeg', sep = ''))
                g1 <- ggplot(sub, aes(time, cumsum(mm)))+geom_point()
                g1 <- g1 + theme_minimal() + xlab("Time") + ylab("Cumilative Throughfall in mm")
                print(g1)
                dev.off()
                td <- as.numeric(diff(sub$time))
                jpeg(file = paste(unique(thrfS$name)[j], as.character(funkydates[i]), 'hist.jpeg', sep = ''))
                r <- hist(td, breaks = c(timbwtips*0:3, max(td)), xlim = range(0, min(1200,max(td))), col = "blue1")
                text(r$mids, r$density, r$counts, adj = c(.5, -.5), col = "blue3")
                lines(r, lty = 3, border = "purple") # -> lines.histogram(*)
                dev.off()
            }else{
              td <- as.numeric(diff(sub$time))
              r <- hist(td, breaks = c(timbwtips*0:3, max(td)), xlim = range(0, min(1200,max(td))), col = "blue1")
            }
            if((r$counts/sum(r$counts))[1] > thresh/100){
              thrfS <- thrfS[!(thrfS$time %in% sub$time),]
            }
          }
        }
      }
    }
    return(list(rf, thrfN, thrfS))
}
