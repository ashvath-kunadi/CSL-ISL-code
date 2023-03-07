kdiagn <- function(k, lim){
  k[abs(kh) > lim] <- NA
  topcor <- data.frame(k, d, z0, wd, Et0dat$u, Et0dat$ustar, monobhL, monobhLwoE, Et0dat$IRtemp, err)  
  topcorg <- topcor[complete.cases(topcor),]
  cor_matrix <- cor(topcorg, use = "complete.obs")
  corrplot.mixed(cor_matrix, lower = "circle", upper = "number", tl.pos = "lt", diag = "u")
  topcor$dnn <- "night"
  topcor$dnn[Et0dat$Rn > 0] <- "day"
  topcor$WDfac <- "N"
  topcor$WDfac[topcor$wd > 225 & topcor$wd < 315] <- "W"
  topcor$WDfac[topcor$wd > 135 & topcor$wd < 225] <- "S"
  topcor$WDfac[topcor$wd > 45 & topcor$wd < 135] <- "E"
  topcor$bintmp <- "high"
  topcor$bintmp[Et0dat$IRtemp > 17 & Et0dat$IRtemp < 23] <- "mid"
  topcor$bintmp[Et0dat$IRtemp <= 17] <- "low"
  ggplot(topcor, aes(kh, Et0dat.ustar)) + facet_grid(dnn~WDfac) + geom_point(aes(col = bintmp), alpha = 0.1) + theme_bw()
  ggplot(topcor, aes(kh, Et0dat.ustar)) + facet_grid(bintmp~WDfac) + geom_point(aes(col = dnn), alpha = 0.1) + theme_bw()
  tapply(topcor$kh, list(topcor$WDfac, topcor$dnn), FUN = function(x)mean(x, na.rm = T))
}