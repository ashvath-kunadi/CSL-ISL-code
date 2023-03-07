evalsimcald <- function(sim, gpfldat, img = F){
  try1 <- sim
  cond <- (try1$subrf != 0)
  condfin <- c(FALSE, cond[-length(cond)])
  try1$thsim <- 0
  wazza <- try1$C[condfin]-try1$C[cond]
  try1$thsim[condfin] <- try1$subrf[cond]-try1$subev[condfin]-wazza
  try2 <- data.frame()
  err <- data.frame()
  for(i in 1:length(unique(try1$IDen))){
    sub <- try1[try1$IDen == unique(try1$IDen)[i],]
    for(j in 1:length(unique(try1$sim))){
      sub2 <- sub[sub$sim == unique(sub$sim)[j], ]
      sub2$crf <- cumsum(sub2$subrf)
      sub2$cev <- cumsum(sub2$subev)
      sub2$cthsim <- cumsum(sub2$thsim)
      sub3 <- sub2[,c(3,4,10,11,12,13,14,15)]
      msub3 <- melt(sub3, id = c("time"))
      if(img){
        jpeg(paste("storm", as.character(unique(sub2$IDen)[1]), "q -", (sub2$q)[1], "L -", (sub2$L)[1], ".jpg"), width = 1000, height = 800)
        g1 <- ggplot(msub3, aes(time, value)) + geom_point(aes(colour = variable)) + theme_classic()
        print(g1)
        dev.off()
      }
      try2 <- rbind(try2, sub3)
      err <- rbind(err, c(sub2$cthsim[length(sub2$cthsim)],sub2$crf[length(sub2$cthsim)],sub2$IDen[1],sub2$q[1], sub2$v[1], sub2$L[1]))
    }
  }
  names(err) <- c("Thsim", "rf","ID", "q", "v", "L")
  norplot <- (gpfldat[[3]])
  norcheck <- data.frame(tapply(norplot$mm, list(norplot$ID, norplot$name), sum))
  norcheck$mean <- apply(norcheck, 1, FUN = function(x)mean(x, na.rm = T))
  norcheck$ID <- as.numeric(rownames(norcheck))
  sthplot <- (gpfldat[[4]])
  sthcheck <- data.frame(tapply(sthplot$mm, list(sthplot$ID, sthplot$name), sum))
  sthcheck$mean <- apply(sthcheck, 1, FUN = function(x)mean(x, na.rm = T))
  sthcheck$ID <- as.numeric(rownames(sthcheck))
  library(dplyr)
  masternor <- norcheck %>% left_join(err, by = "ID")
  mastersth <- sthcheck %>% left_join(err, by = "ID")
  nordifsim <- data.frame()
  sthdifsim <- data.frame()
  for(i in 1:6){
      diff <- masternor[,c(10,11,12)]
      diff$simdif <- masternor[,i]-masternor[,8]
      simis <- tapply(diff$simdif, diff$q, FUN = function(x)sum(x, na.rm = T))
      topl <- data.frame(cbind(simis, as.numeric(names(simis))))
      names(topl) <- c("simdif", "q")
      par <- diff[1:100,]
      dat <- cbind(topl, par[order(par$q), 1:3])
      dat$ID <- names(masternor)[i]
      nordifsim <- rbind(nordifsim, dat)
      
      diff <- mastersth[,c(10,11,12)]
      diff$simdif <- mastersth[,i]-mastersth[,8]
      simis <- tapply(diff$simdif, diff$q, FUN = function(x)sum(x, na.rm = T))
      topl <- data.frame(cbind(simis, as.numeric(names(simis))))
      names(topl) <- c("simdif", "q")
      par <- diff[1:100,]
      dat <- cbind(topl, par[order(par$q), 1:3])
      dat$ID <- names(mastersth)[i]
      sthdifsim <- rbind(sthdifsim, dat)
  }
  topl <- nordifsim[,-3]
  topl <- melt(topl, id = c("simdif", "ID"))
  g2 <- ggplot(topl, aes(value, simdif))+facet_grid(ID~variable, scales = "free")+geom_point(alpha = 0.5)+geom_smooth()+theme_classic()+geom_hline(yintercept = 0)
  jpeg("caldsim nor.jpg", width = 800, height = 1000)
  print(g2)
  dev.off()
  topl <- sthdifsim[,-3]
  topl <- melt(topl, id = c("simdif", "ID"))
  g2 <- ggplot(topl, aes(value, simdif))+facet_grid(ID~variable, scales = "free")+geom_point(alpha = 0.5)+geom_smooth()+theme_classic()+geom_hline(yintercept = 0)
  jpeg("caldsim sth.jpg", width = 800, height = 1000)
  print(g2)
  dev.off()
  return(list(masternor, mastersth, nordifsim, sthdifsim))
}
