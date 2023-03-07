formatrg <- function(rgdatfr){
    library(lubridate)
    rmna <- is.na(rgdatfr[,4])
    rgdatfr <- rgdatfr[!rmna, ]
    time <- dmy_hms(paste(rgdatfr[,2], rgdatfr[,3], sep = " "))
    tips <- c(1, diff(rgdatfr[, 4], 1))
    return(cbind(name = rgdatfr[1, 7],time, tips))
}