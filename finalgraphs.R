library(reshape2)
library(ggplot2)
library(viridis)
library(ggExtra)
library(scales)
library(MASS)
zeta <- seq(0.01, 3, by = 0.01)
n <- length(zeta)
Inoue <- function(zeta, beta, cdahc) {exp(cdahc/2/beta^2*(zeta-1))}
Cowanorg <-  function(zeta, beta, cdah,d0h) {sqrt(sinh(sqrt(cdah/beta/0.4/(1-d0h))*zeta)/sinh(sqrt(cdah/beta/0.4/(1-d0h))))}
Thom <- function(zeta, beta, cdah, d0h) {1/(1+(sqrt(cdah/6/beta/0.4/(1-d0h))*(1-zeta)))^2}
Wang <- function(zeta, beta, cdahc,z0gbyhc){
  (besselI (2*sqrt (cdahc/beta/0.4*zeta), 0)-
     (besselI (2*sqrt (z0gbyhc*cdahc/beta/0.4), 0)*
        besselK (2*sqrt (cdahc/beta/0.4*zeta), 0)/besselK (2*sqrt (z0gbyhc*cdahc/beta/0.4), 0)))/
    (besselI (2*sqrt (cdahc/beta/0.4), 0) - 
       (besselI (2*sqrt (z0gbyhc*cdahc/beta/0.4), 0)*besselK (2*sqrt (cdahc/beta/0.4), 0)/
          besselK (2*sqrt (z0gbyhc*cdahc/beta/0.4), 0)))}

candens <- c(0.004,0.04,0.4,4)
betafr <- c(0.05, 0.15, 0.33, 0.33)
canlabs <- c("Sparse-    =0.05,       =0.004", "Transitional-    =0.15,       =0.04", "Dense-    =0.33,       =0.4", "Denser-    =0.33,       =0.4")
canlabs <- factor(canlabs, levels = c("Sparse-    =0.05,       =0.004", "Transitional-    =0.15,       =0.04", "Dense-    =0.33,       =0.4", "Denser-    =0.33,       =0.4"))
mod <- c("Inoue", "Cowan", "Thom", "Wang")
top <- data.frame(zeta = rep(zeta, length(candens)*4), uzuh = NA,
                  Model = rep(mod, each = n*length(candens)), 
                  candens = rep(candens, each = n),
                  canlabs = rep(canlabs, each = n),
                  beta = rep(betafr, each = n), d0 = NA, z0 = NA)

cw68 <- function(x,y,par){
  ((2*x/(1 - par[1])/0.4)-(sqrt(y/x/0.4/(1 - par[1]))/tanh(sqrt(y/x/0.4/(1 - par[1])))))
}
for(i in 1:nrow(top)){
  if(i <= 4*n){
    top$d0[i] <- (1 - (2*top$beta[i]^3/top$candens[i]/0.4))
  }else if(i > 4*n & i <= (8*n)+1){
    result <- optim(par = 0.999999, fn = cw68, x = top$beta[i], y =top$candens[i],
                    method = "Brent", lower = 0.0000000000001, upper = 0.9999999999999)
    top$d0[i] <- result$par
  }else if(i > 8*n & i <= (12*n)+1){
    top$d0[i] <- (1 - (6/4*top$beta[i]^3/top$candens[i]/0.4))
  }else{
    top$d0[i] <- (1 - (top$beta[i]*(besselI(2*sqrt(top$candens[i]/0.4/top$beta[i]),0) + 
                            (besselI(2*sqrt(0.00001/top$beta[i]/0.4*top$candens[i]),0)*
                               besselK(2*sqrt(top$candens[i]/top$beta[i]/0.4),0)/
                               besselK(2*sqrt(0.00001*top$candens[i]/top$beta[i]/0.4),0)))/
                       (0.4*sqrt(top$candens[i]/top$beta[i]/0.4)*
                          (besselI(2*sqrt(top$candens[i]/top$beta[i]/0.4),1) - 
                             (besselI(2*sqrt(0.00001/top$beta[i]/0.4*top$candens[i]),0)*
                                besselK(2*sqrt(top$candens[i]/top$beta[i]/0.4),1)/
                                besselK(2*sqrt(0.00001/top$beta[i]/0.4*top$candens[i]),0))))))
  }
}
top$d0[top$d0<=0] <- 0.0000000000001
top$d0[top$d0>=1] <- 0.9999999999999

top$z0 <- (1-top$d0)/exp(0.4/top$beta)

for(i in 1:nrow(top)[1]){
  if(top$zeta[i]>=1){
    top$uzuh[i] <- top$beta[i]/0.4*log((top$zeta[i]-top$d0[i])/top$z0[i])
  }else{
    if(top$Model[i] == ("Inoue")){
      top$uzuh[i] <- Inoue(top$zeta[i], top$beta[i], top$candens[i])
    }else if(top$Model[i] == ("Cowan")){
      top$uzuh[i] <- Cowanorg(top$zeta[i], top$beta[i], top$candens[i], top$d0[i])
    }else if(top$Model[i] == ("Thom")){
      top$uzuh[i] <- Thom(top$zeta[i], top$beta[i], top$candens[i], top$d0[i])
    }else if(top$Model[i] == ("Wang")){
      top$uzuh[i] <- Wang(top$zeta[i], top$beta[i], top$candens[i], 0.00001)
    }
  }
}

top2 <- top
g1 <- ggplot(top, aes(zeta, uzuh)) + facet_wrap(.~canlabs,scales = "free") + geom_line(aes(col = Model), size = 1)
g1+theme_classic()+coord_flip() + ylab(expression(frac(U(z),U[H]))) + xlab(expression(frac(z,H[c]))) + 
  scale_color_manual(values = c("limegreen","springgreen", "forestgreen", "darkgreen"))
g1 <- ggplot(top, aes(zeta, uzuh)) + facet_wrap(.~Model,scales = "free") + geom_line(aes(col = as.factor(candens))) + 
  ylab(expression(frac(u[z],u[h]))) + xlab(expression(frac(z,H[c])))
g1+theme_classic()+coord_flip() 

my_breaks = unique(top$cdahc)[-3]
g1 <- ggplot(top, aes(zeta, uzuh)) + facet_wrap(.~Model) + geom_point(aes(col = candens)) + 
  ylab(expression(frac(u[z],u[h]))) + xlab(expression(frac(z,H[c])))
g1+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+coord_flip() + scale_fill_gradient(name = "count", trans = "log",
                      breaks = my_breaks, labels = my_breaks)
##changes for Sally
#Panel A
cbet <- 0.1
cbetustr <- 0.5
varcdah <- c(0.01,0.07,0.5,10)
panela <- data.frame(zeta = rep(zeta, length(varcdah)),beta = cbet, ustar = cbetustr, cdah = rep(varcdah, each = length(zeta)))
panela$uzuh <- NA
panela$d0 <- NA
for(i in 1:nrow(panela)){
  panela$uzuh[i] <- Wang(panela$zeta[i], panela$beta[i], panela$cdah[i], 0.00001)
  panela$d0[i] <- (1 - (panela$beta[i]*(besselI(2*sqrt(panela$cdah[i]/0.4/panela$beta[i]),0) + 
                                    (besselI(2*sqrt(0.00001/panela$beta[i]/0.4*panela$cdah[i]),0)*
                                       besselK(2*sqrt(panela$cdah[i]/panela$beta[i]/0.4),0)/
                                       besselK(2*sqrt(0.00001*panela$cdah[i]/panela$beta[i]/0.4),0)))/
                       (0.4*sqrt(panela$cdah[i]/panela$beta[i]/0.4)*
                          (besselI(2*sqrt(panela$cdah[i]/panela$beta[i]/0.4),1) - 
                             (besselI(2*sqrt(0.00001/panela$beta[i]/0.4*panela$cdah[i]),0)*
                                besselK(2*sqrt(panela$cdah[i]/panela$beta[i]/0.4),1)/
                                besselK(2*sqrt(0.00001/panela$beta[i]/0.4*panela$cdah[i]),0))))))
}
panela$z0m <- with(panela, (1-d0)/exp(0.4/beta))
panela[panela$zeta>1,5] <- with(panela[panela$zeta>1,], beta/0.4*log((zeta-d0)/z0m))
ndlogprf <- function(beta,zeta,d0,z0m){
  ans <- NA
  for(i in 1:length(zeta)){
    if(zeta[i] > d0[i]){
      ans[i] <- beta[i]/0.4*log((zeta[i]-d0[i])/z0m[i])
    }else{
      ans[i] <- 0
    }
  }
  ans[ans < 0] <- 0
  return(ans)
}
panela$lgprf <- with(panela, ndlogprf(beta,zeta,d0,z0m))
ggplot(panela, aes(uzuh, zeta, group = cdah)) + geom_line(aes(x = lgprf, y = zeta, group = cdah, col = cdah), linetype = 2, size = 1) + 
  geom_line(aes(col = cdah), size = 1) +  
  scale_color_gradient(name = expression(paste(c[d],aH[c])), trans = "log",breaks = varcdah, labels = varcdah, low = "darkgoldenrod", high = "forestgreen") + 
  theme_classic() + xlab(expression(frac(U[z],U[H]))) + ylab(expression(frac(z, H[c]))) + ggtitle(expression(paste(beta,"=0.1")))

#Panel B
varbet <- c(0.125,0.25,0.375,0.5)
cbetustr <- 0.5
ccdah <- 0.8
panelb <- data.frame(zeta = rep(zeta, length(varbet)),beta = rep(varbet, each = length(zeta)), ustar = cbetustr, cdah = ccdah)
panelb$uzuh <- NA
panelb$d0 <- NA
for(i in 1:nrow(panelb)){
  panelb$uzuh[i] <- Wang(panelb$zeta[i], panelb$beta[i], panelb$cdah[i], 0.00001)
  panelb$d0[i] <- (1 - (panelb$beta[i]*(besselI(2*sqrt(panelb$cdah[i]/0.4/panelb$beta[i]),0) + 
                                          (besselI(2*sqrt(0.00001/panelb$beta[i]/0.4*panelb$cdah[i]),0)*
                                             besselK(2*sqrt(panelb$cdah[i]/panelb$beta[i]/0.4),0)/
                                             besselK(2*sqrt(0.00001*panelb$cdah[i]/panelb$beta[i]/0.4),0)))/
                          (0.4*sqrt(panelb$cdah[i]/panelb$beta[i]/0.4)*
                             (besselI(2*sqrt(panelb$cdah[i]/panelb$beta[i]/0.4),1) - 
                                (besselI(2*sqrt(0.00001/panelb$beta[i]/0.4*panelb$cdah[i]),0)*
                                   besselK(2*sqrt(panelb$cdah[i]/panelb$beta[i]/0.4),1)/
                                   besselK(2*sqrt(0.00001/panelb$beta[i]/0.4*panelb$cdah[i]),0))))))
}
panelb$z0m <- with(panelb, (1-d0)/exp(0.4/beta))
panelb[panelb$zeta>1,5] <- with(panelb[panelb$zeta>1,], beta/0.4*log((zeta-d0)/z0m))
ndlogprf <- function(beta,zeta,d0,z0m){
  ans <- NA
  for(i in 1:length(zeta)){
    if(zeta[i] > d0[i]){
      ans[i] <- beta[i]/0.4*log((zeta[i]-d0[i])/z0m[i])
    }else{
      ans[i] <- 0
    }
  }
  ans[ans < 0] <- 0
  return(ans)
}
panelb$lgprf <- with(panelb, ndlogprf(beta,zeta,d0,z0m))
ggplot(panelb, aes(uzuh, zeta, group = beta)) + geom_line(aes(x = lgprf, y = zeta, group = beta, col = beta), linetype = 2, size = 1) + 
  geom_line(aes(col = beta), size = 1) +  scale_color_gradient(name = expression(beta), breaks = varbet, labels = varbet) + theme_classic() +
  xlab(expression(frac(U[z],U[H]))) + ylab(expression(frac(z, H[c]))) + ggtitle(expression(paste(c[d],aH[c],"=0.8")))

#Panel C
cbet <- 0.3
ccdah <- 0.4
panelc <- data.frame(zeta = rep(zeta, 4),beta = cbet, ustar = cbetustr, cdah = ccdah, Model = rep(c("Inoue", "Cowan", "Thom", "Wang"), each = n))
panelc$uzuh <- NA
panelc$d0 <- NA
for(i in 1:nrow(panelc)){
  if(i <= n){
    panelc$uzuh[i] <- Inoue(panelc$zeta[i], panelc$beta[i], panelc$cdah[i])
    panelc$d0[i] <- (1 - (2*panelc$beta[i]^3/panelc$cdah[i]/0.4))
  }else if(i > n & i <= 2*n){
    result <- optim(par = 0.999999, fn = cw68, x = panelc$beta[i], y =panelc$cdah[i],
                    method = "Brent", lower = 0.0000000000001, upper = 0.9999999999999)
    panelc$d0[i] <- result$par
    panelc$uzuh[i] <- Cowanorg(panelc$zeta[i], panelc$beta[i], panelc$cdah[i], panelc$d0[i])
  }else if(i > 2*n & i <= 3*n){
    panelc$d0[i] <- (1 - (6/4*panelc$beta[i]^3/panelc$cdah[i]/0.4))
    panelc$uzuh[i] <- Thom(panelc$zeta[i], panelc$beta[i], panelc$cdah[i], panelc$d0[i])
  }else{
    panelc$uzuh[i] <- Wang(panelc$zeta[i], panelc$beta[i], panelc$cdah[i], 0.00001)
    panelc$d0[i] <- (1 - (panelc$beta[i]*(besselI(2*sqrt(panelc$cdah[i]/0.4/panelc$beta[i]),0) + 
                                      (besselI(2*sqrt(0.00001/panelc$beta[i]/0.4*panelc$cdah[i]),0)*
                                         besselK(2*sqrt(panelc$cdah[i]/panelc$beta[i]/0.4),0)/
                                         besselK(2*sqrt(0.00001*panelc$cdah[i]/panelc$beta[i]/0.4),0)))/
                         (0.4*sqrt(panelc$cdah[i]/panelc$beta[i]/0.4)*
                            (besselI(2*sqrt(panelc$cdah[i]/panelc$beta[i]/0.4),1) - 
                               (besselI(2*sqrt(0.00001/panelc$beta[i]/0.4*panelc$cdah[i]),0)*
                                  besselK(2*sqrt(panelc$cdah[i]/panelc$beta[i]/0.4),1)/
                                  besselK(2*sqrt(0.00001/panelc$beta[i]/0.4*panelc$cdah[i]),0))))))
  }
}

panelc$z0m <- with(panelc, (1-d0)/exp(0.4/beta))
panelc[panelc$zeta>1,6] <- with(panelc[panelc$zeta>1,], beta/0.4*log((zeta-d0)/z0m))
panelc$lgprf <- with(panelc, ndlogprf(beta,zeta,d0,z0m))
ggplot(panelc, aes(uzuh, zeta, group = Model)) + geom_line(aes(x = lgprf, y = zeta, group = Model, col = Model), linetype = 2, size = 1) + 
  geom_line(aes(col = Model), size = 1)  + theme_classic() + scale_color_manual(values= c("darkorange4", "goldenrod4", "goldenrod", "lightgoldenrod3")) +
  xlab(expression(frac(U[z],U[H]))) + ylab(expression(frac(z, H[c]))) + ggtitle(expression(paste(c[d],aH[c],"=0.1,",beta,"=1"))) 


cw68 <- function(x,y,par){
  ((2*x/(1 - par[1])/0.4)-(sqrt(y/x/0.4/(1 - par[1]))/tanh(sqrt(y/x/0.4/(1 - par[1])))))
}
for(i in 1:nrow(top)){
  
}
top$d0[top$d0<=0] <- 0.0000000000001
top$d0[top$d0>=1] <- 0.9999999999999

top$z0 <- (1-top$d0)/exp(0.4/top$beta)

for(i in 1:nrow(top)[1]){
  if(top$zeta[i]>=1){
    top$uzuh[i] <- top$beta[i]/0.4*log((top$zeta[i]-top$d0[i])/top$z0[i])
  }else{
    if(top$Model[i] == ("Inoue")){
      top$uzuh[i] <- Inoue(top$zeta[i], top$beta[i], top$candens[i])
    }else if(top$Model[i] == ("Cowan")){
      top$uzuh[i] <- Cowanorg(top$zeta[i], top$beta[i], top$candens[i], top$d0[i])
    }else if(top$Model[i] == ("Thom")){
      top$uzuh[i] <- Thom(top$zeta[i], top$beta[i], top$candens[i], top$d0[i])
    }else if(top$Model[i] == ("Wang")){
      top$uzuh[i] <- Wang(top$zeta[i], top$beta[i], top$candens[i], 0.00001)
    }
  }
}

top2 <- top
g1 <- ggplot(top, aes(zeta, uzuh)) + facet_wrap(.~canlabs,scales = "free") + geom_line(aes(col = Model), size = 1)
g1+theme_classic()+coord_flip() + ylab(expression(frac(U(z),U[H]))) + xlab(expression(frac(z,H[c]))) + 
  scale_color_manual(values = c("limegreen","springgreen", "forestgreen", "darkgreen"))
g1 <- ggplot(top, aes(zeta, uzuh)) + facet_wrap(.~Model,scales = "free") + geom_line(aes(col = as.factor(cdahc))) + 
  ylab(expression(frac(u[z],u[h]))) + xlab(expression(frac(z,H[c])))
g1+theme_classic()+coord_flip() 


betacontr <- seq(0.01, 1, by = 0.01)
cdahcontr <- unique(c(c(2:10 %o% 10^(-2:1)),seq(0.1,0.5, by = 0.01), seq(1,5, by = 0.1), seq(10, 50, by = 1), 0.018, 0.015,0.013, 0.011, 0.01, 0.009,0.007,0.006,0.005,0.004))
cntr <- data.frame(beta = rep(betacontr, each=length(cdahcontr)), cdah = rep(cdahcontr, length(betacontr)))
cntr$Inoue <- with(cntr, (1 - (2*beta^3/cdah/0.4)))
cntr$Cowan <- NA
cw68 <- function(x,y,par){
  ((2*x/(1 - par[1])/0.4)-(sqrt(y/x/0.4/(1 - par[1]))/tanh(sqrt(y/x/0.4/(1 - par[1])))))
}
cw682 <- function(x,y,par){
  (x^2-((sqrt(y*x*0.4*(1-par[1])))/sinh(sqrt(y/x/0.4/(1-par[1])))*(cosh(sqrt(y/x/0.4/(1-par[1])))-1)))
}
checkcw2 <- NA
for(i in 1:dim(cntr)[1]){
  result <- optim(par = 0.999999, fn = cw68, x = cntr[i,1], y =cntr[i,2],
                  method = "Brent", lower = 0.0000000000001, upper = 0.9999999999999)
  cntr[i,4] <- result$par
  result <- optim(par = 0.999999, fn = cw682, x = cntr[i,1], y =cntr[i,2],
                  method = "Brent", lower = 0.0000000000001, upper = 0.9999999999999)
  checkcw2[i] <- result$par
}
cntr$Thom <- with(cntr, 1 - 6*beta^3/cdah/(4*0.4))
cntr$Wang <- with(cntr, 1 - (beta*(besselI(2*sqrt(cdah/0.4/beta),0) + 
                                       (besselI(2*sqrt(0.001/beta/0.4*cdah),0)*
                                          besselK(2*sqrt(cdah/beta/0.4),0)/
                                          besselK(2*sqrt(0.001*cdah/beta/0.4),0)))/
                                 (0.4*sqrt(cdah/beta/0.4)*
                                    (besselI(2*sqrt(cdah/beta/0.4),1) - 
                                       (besselI(2*sqrt(0.001/beta/0.4*cdah),0)*
                                          besselK(2*sqrt(cdah/beta/0.4),1)/
                                          besselK(2*sqrt(0.001/beta/0.4*cdah),0))))))

mcntrt <- melt(cntr, id = c("beta", "cdah"))
mcntrt$z0 <- with(mcntrt, (1-value)/exp(0.4/beta))
mcntrt$value[mcntrt$value<0.001] <- NA
mcntrt$value[mcntrt$value>1] <- NA
mcntrt$z0[mcntrt$z0<0.00001] <- NA
mcntrt$z0[mcntrt$z0>1] <- NA
library(viridis)

ggplot(mcntrt, aes(x = cdah, y = beta, z = value)) +facet_grid(variable~.)+ geom_contour_filled(bins = 10) + 
  scale_x_continuous(trans = 'log10',breaks=c(0.01, 0.1, 1, 10, 100), 
                     labels=c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2))) + 
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_classic() + ylab(expression(beta)) + xlab(expression(paste(c[d],aH[c]))) + 
  geom_vline(xintercept = 0.004, col = "#2C7BB6", size = 1) + #geom_hline(yintercept = 0.3333, ltype = 2) +
  geom_vline(xintercept = 0.04, col = "#ABD9E9", size = 1) + geom_vline(xintercept = 0.4, col = "#FDAE61", size = 1) + 
  geom_vline(xintercept = 4, col = "#D7191C", size = 1) + theme(strip.background = element_blank()) +
  scale_fill_viridis(discrete = TRUE,name = expression(frac(d[0],H[c])),labels = c("<0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"))
  

ggplot(mcntrt[complete.cases(mcntrt),], aes(x = cdah, y = beta, z = z0)) +facet_grid(variable~.)+ geom_contour_filled(bins = 10) + 
  scale_x_continuous(trans = 'log10',breaks=c(0.01, 0.1, 1, 10, 100), 
                     labels=c(expression(10^-2), expression(10^-1), expression(10^0), expression(10^1), expression(10^2))) + 
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_classic() + ylab(expression(beta)) + xlab(expression(paste(c[d],aH[c]))) + 
  geom_vline(xintercept = 0.004, col = "#2C7BB6", size = 1) + #geom_hline(yintercept = 0.3333, ltype = 2) + 
  geom_vline(xintercept = 0.04, col = "#ABD9E9", size = 1) + geom_vline(xintercept = 0.4, col = "#FDAE61", size = 1) + 
  geom_vline(xintercept = 4, col = "#D7191C", size = 1) + theme(strip.background = element_blank()) +
  scale_fill_viridis(discrete = TRUE,name = expression(frac(z["0m"],H[c])),
                     labels = c("<0.07", "0.14", "0.21", "0.28", "0.35", "0.42", "0.49", "0.56", "0.63", "0.7"))


topl <- cntr
topl$density <- NA
topl$density[topl$cdah == 0.004] <- "Sparse"
topl$density[topl$cdah == 0.04] <- "Transitional"
topl$density[topl$cdah == 0.4] <- "Dense"
topl$density[topl$cdah == 4] <- "Denser"
topl <- topl[!is.na(topl$density),]
topl$density <- factor(topl$density, levels = c("Denser", "Dense", "Transitional", "Sparse"))
mtopl <- melt(topl[,-2], id = c("beta", "density"))
mtopl$z0 <- with(mtopl, (1-value)/exp(0.4/beta))
mtopl$value[mtopl$value>1] <- NA
mtopl$value[mtopl$value<0] <- NA
ggplot(mtopl, aes(beta, value)) + geom_line(aes(col = density), size = 1) + 
  facet_grid(variable~.) + theme_classic() + ylab(expression(frac(d[0],H[c]))) +   
  scale_color_brewer(palette = "RdYlBu", name = "Canopy\n Density") + xlab(expression(beta)) + 
  scale_y_continuous(breaks = c(0,0.5,1)) + theme(strip.background = element_blank())


mtopl$z0[mtopl$z0>1] <- NA
mtopl$z0[mtopl$z0<0] <- NA
mtopl$z0[mtopl$value<10^-8] <- NA
ggplot(mtopl[complete.cases(mtopl),], aes(beta, z0)) + geom_line(aes(col = density), size = 1) + 
  facet_grid(variable~.) + theme_classic() + ylab(expression(frac(z["0m"],H[c]))) +     scale_y_continuous(breaks = c(0,0.3,0.6)) +
  scale_color_brewer(palette = "RdYlBu", name = "Canopy\n Density") + xlab(expression(beta)) + theme(strip.background = element_blank())
  


limin <- function(cdahc) {return((0.4*cdahc/2)^(1/3))}
limth <- function(cdahc) {return((0.4*4*cdahc/6)^(1/3))}
limcwfn <- function(x,par){  sum((2*par[1]/0.4)-(sqrt(x/par[1]/0.4)/tanh(sqrt(x/par[1]/0.4))))^2}
limcw <- function(cdahc) {

    results <- optim(par = 1, fn = limcwfn, x = cdahc, method = "Brent", lower = 0.0000000000001, upper = 5000)
    ans <-results$par
  
  return(ans)
}
limwgfn <- function(x,par){  sum(1 -  (par[1]*(besselI(2*sqrt(x/0.4/par[1]),0) + 
                                               (besselI(2*sqrt(0.001/par[1]/0.4*x),0)*
                                                  besselK(2*sqrt(x/par[1]/0.4),0)/
                                                  besselK(2*sqrt(0.001/par[1]/0.4*x),0)))/
                                            (0.4*sqrt(x/par[1]/0.4)*
                                               (besselI(2*sqrt(x/par[1]/0.4),1) - 
                                                  (besselI(2*sqrt(0.001/par[1]/0.4*x),0)*
                                                     besselK(2*sqrt(x/par[1]/0.4),1)/
                                                     besselK(2*sqrt(0.001/par[1]/0.4*x),0))))))^2}
limwg <- function(cdahc) {

    results <- optim(par = 1, fn = limwgfn, x = cdahc, method = "Brent", lower = 0.0000000000001, upper = 5000)
    ans <-results$par
  
  return(ans)
}
cdahcontr <- unique(c(c(2:10 %o% 10^(-2:1)),seq(0.1,0.5, by = 0.01), seq(1,5, by = 0.1), seq(10, 50, by = 1), 0.004, 0.0004))

ans <- 0
lim <- data.frame(cdahc = rep(sort(cdahcontr),4), model = NA, beta = NA)
lim$model <- rep(c("Inoue", "Cowan", "Thom", "Wang"), each = nrow(lim)/4)
lim$beta[1:(nrow(lim)/4)] <-limin(sort(cdahcontr))
for(i in 1:(nrow(lim)/4)){  ans[i] <- limcw(lim$cdahc[i])  }
lim$beta[(nrow(lim)/4+1):(nrow(lim)*2/4)] <- ans
lim$beta[(nrow(lim)*2/4+1):(nrow(lim)*3/4)] <-limth(lim$cdahc[(nrow(lim)*2/4+1):(nrow(lim)*3/4)])
ans <- 0
for(i in 1:(nrow(lim)/4)){  ans[i] <- limwg(lim$cdahc[i])  }
lim$beta[(nrow(lim)*3/4+1):nrow(lim)] <-ans
ggplot(lim, aes(cdahc, beta)) + geom_line(aes(col = model)) + theme_classic() + ylab(expression(beta)) + xlab(expression(c[d]~aH[c])) + scale_x_continuous(trans = 'log10') + ylim(0,5)

##setting up working directory as the place with codes 
setwd("H:/My Documents/Project and related work/Code/Papier Code")

##reading the fuctions for this code
source("fileread2.R")
source("formatrg.R")
source("daily.R")
source("dailymm.R")
source("ozfluxdataextractrf.R")
source("ozfluxdataextractrfas.R")
source("ozfluxdataextractrfyanc.R")
source("calirm.R")
source("gapfillingimp.R")
source("strmsepev.R")
source("aclnrintp.R")
source("estimateplots.R")
source("rutsprsim.R")
source("suggestedchanges.R")
source("plotstroms.R")
source("evalsimrt.R")
source("evalsimrtsp.R")
source("rutsim.R")
source("caldsimv3.R")
source("evalsimcald.R")
source("Gashdataprep.R")
source("gashsim.R")
source("evalsimgash.R")
source("figureoutET.R")
source("pasqstability.R")
source("momentum.R")
source("enbal.R")
source("br439fin.R")
source("zondnlsnslp.R")
source("onlz0.R")
source("numderdmod.R")
source("khnkvcalcfin.R")
source("z0hnlsfin.R")
source("kdiagn.R")
source("fixkod.R")
source("rmna.R")
source("phim.R")
source("Ameriflux.R")
source("Ameriflux - Prr.R")

##reading in all the data
path2dat <- "H:/My Documents/Project and related work/Rainfall-Throughfall"
dat <- fileread2(path2dat)

##making the data usable
thrfN <- data.frame()
thrfS <- data.frame()
rf <- data.frame()
for(i in 1:length(dat)){
  if(grepl("9092", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE) | grepl("9080", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9047", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9093", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9078", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)){
    tempdat <- formatrg(dat[[i]])
    thrfN <- rbind(thrfN, tempdat)
  }
  if(grepl("9086", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE) | grepl("9079", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9076", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9082", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9077", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)){
    tempdat <- formatrg(dat[[i]])
    thrfS <- rbind(thrfS, tempdat)
  }
  if(grepl("T*RG", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9097", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE) | grepl("9098", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)| grepl("9099", substr(dat[[i]]$Name[1], 1, 4), ignore.case = TRUE)){
    tempdat <- formatrg(dat[[i]])
    rf <- rbind(rf, tempdat)
  } 
}

##formatting the data correctly
thrfN$time <- as.POSIXct(as.POSIXlt(as.integer(as.character(thrfN$time)), origin = "1970-01-01", tz = "GMT"))
thrfS$time <- as.POSIXct(as.POSIXlt(as.integer(as.character(thrfS$time)), origin = "1970-01-01", tz = "GMT"))
rf$time <- as.POSIXct(as.POSIXlt(as.integer(as.character(rf$time)), origin = "1970-01-01", tz = "GMT"))
thrfN$tips <- as.integer(as.character(thrfN$tips))*0.2
thrfS$tips <- as.integer(as.character(thrfS$tips))*0.2
rf$tips <- as.integer(as.character(rf$tips))*0.2
thrfN$name <- substr(as.character(thrfN$name), 3,4)
thrfS$name <- substr(as.character(thrfS$name), 3,4)
rf$name <- substr(as.character(rf$name), 1,4)
names(rf) <- c("name", "time", "mm")
names(thrfN) <- c("name", "time", "mm")
names(thrfS) <- c("name", "time", "mm")

##making room in the memory
rm(dat, tempdat, fileread2, formatrg, i)

##data before 2017 seems quite irregular so we'll remove that
thrfN <- thrfN[thrfN$time > ymd("2017-01-01"),]
thrfS <- thrfS[thrfS$time > ymd("2017-01-01"),]
mintime <- min(c(thrfN$time, thrfS$time))
rf <- rf[rf$time>mintime,]
maxtime <- max(c(thrfN$time, thrfS$time, rf$time))

##obtaining relavant flux tower data
setwd("H:/My Documents/Project and related work/Cool gingin/L3data")
fluxdat11 <- ozfluxdataextractrf("Gingin_2011_L3.nc", "2011.txt")
fluxdat12 <- ozfluxdataextractrf("Gingin_2012_L3.nc", "2012.txt")
fluxdat13 <- ozfluxdataextractrf("Gingin_2013_L3.nc", "2013.txt")
fluxdat14 <- ozfluxdataextractrf("Gingin_2014_L3.nc", "2014.txt")
fluxdat15 <- ozfluxdataextractrf("Gingin_2015_L3.nc", "2015.txt")
fluxdat16 <- ozfluxdataextractrf("Gingin_2016_L3.nc", "2016.txt")
fluxdat17 <- ozfluxdataextractrf("Gingin_2017_L3.nc", "2017.txt")
fluxdat18 <- ozfluxdataextractrf("Gingin_2018_L3.nc", "2018.txt")
fluxdat19 <- ozfluxdataextractrf("Gingin_2019_L3.nc", "2019.txt")

##zero plane displacement (d) and roughness length (z0)
wdp <- rbind(fluxdat15[[6]], fluxdat16[[6]], fluxdat17[[6]], fluxdat18[[6]], fluxdat19[[6]])
wsp <- rbind(fluxdat15[[5]], fluxdat16[[5]], fluxdat17[[5]], fluxdat18[[5]], fluxdat19[[5]])
tap <- rbind(fluxdat15[[7]], fluxdat16[[7]], fluxdat17[[7]], fluxdat18[[7]], fluxdat19[[7]])
ahp <- rbind(fluxdat15[[8]], fluxdat16[[8]], fluxdat17[[8]], fluxdat18[[8]], fluxdat19[[8]])
stp <- list(fluxdat15[[10]], fluxdat16[[10]], fluxdat17[[10]], fluxdat18[[10]], fluxdat19[[10]])
Et0dat <- rbind(fluxdat15[[2]], fluxdat16[[2]], fluxdat17[[2]], fluxdat18[[2]], fluxdat19[[2]])
n <- dim(Et0dat)[1]
##Passquill stability
dtdzslp <- c(fluxdat15[[9]], fluxdat16[[9]], fluxdat17[[9]], fluxdat18[[9]], fluxdat19[[9]])
dudz <- br439prf(n, 6.5, wsp)
dTdz <- br439prf(n, 6.5, tap,3)
dahdz <- br439prf(n, 6.5, ahp,5)
drhodz <- br439prf(n, 6.5, ahp, 6)
Ridat <- merge(dudz, dTdz,  by = c("zeval", "time"))
Ridat <- merge(Ridat, dahdz,  by = c("zeval", "time"))
Ridat <- Ridat[,-c(4,5,7,8,10,11)]
if(dim(Ridat)[1] == 0){
  rws <- unique(dudz$zeval)
  cols <- unique(dahdz$zeval)
  pick <- matrix(data = NA, nrow = length(rws), ncol = length(cols))
  for(i in 1:length(rws)){
    for(j in 1:length(cols)){
      pick[i,j] <- abs(rws[i]-cols[j])
    }
  }
  loc <- which(pick == min(pick), arr.ind = TRUE)
  subuz <- dudz[dudz$zeval == rws[loc[1]],c(1,4)]
  subahz <- dahdz[dahdz$zeval == cols[loc[2]],c(1,4)]
  subtz <- dTdz[dTdz$zeval == cols[loc[2]],c(1,4)]
  Ridat <- merge(subuz, subtz,  by = "time")
  Ridat <- merge(Ridat, subahz,  by = "time")
  Ridat$zeval <- mean(c(rws[loc[1]], cols[loc[2]]))
  Ridat <- Ridat[,c(5,1:4)]
}
names(Ridat) <- c("zeval", "time", "dudz", "dTdz", "dahdz")
##there are two 2019-01-01 points so need to correct for that mannualy in Ridat
Ridat <- Ridat[-(70130:70135),]
Richno <- with(Ridat, 9.806/Et0dat$Tair/dudz^2*(dTdz+(0.61*Et0dat$Tair*(dahdz))))
dtdzpsqstb <- pasqstability((dtdzslp*100), Et0dat, Richno)

##energy balance
fltim <- c(fluxdat15[[1]], fluxdat16[[1]], fluxdat17[[1]], fluxdat18[[1]], fluxdat19[[1]])
ebal <- enbal(Et0dat, fltim, dtdzpsqstb)

##whats up with L and dtvdz
ginmom <- momentum(Et0dat, wdp,wsp,tap,ahp, dtdzpsqstb, dudz, dTdz, dahdz, 15.5, 7, F, plotloc = "H:/My Documents/Project and related work/Code/Papier Code", "Gingin", 0.01)

##obtaining relavant flux tower data
setwd("H:/My Documents/Project and related work/Cool gingin/L3data")
asfluxdat12 <- ozfluxdataextractrfas("ASM_EC_2012_L3_Corrected.nc", "as2012.txt")
asfluxdat13 <- ozfluxdataextractrfas("ASM_EC_2013_L3_Corrected.nc", "as2013.txt")
asfluxdat14 <- ozfluxdataextractrfas("ASM_EC_2014_L3_Corrected.nc", "as2014.txt")
asfluxdat15 <- ozfluxdataextractrfas("ASM_EC_2015_L3_Corrected.nc", "as2015.txt")
asfluxdat16 <- ozfluxdataextractrfas("ASM_EC_2016_L3_Corrected.nc", "as2016.txt")
asfluxdat17 <- ozfluxdataextractrfas("ASM_EC_2017_L3_Corrected.nc", "as2017.txt")
asfluxdat18 <- ozfluxdataextractrfas("ASM_EC_2018_L3_Corrected.nc", "as2018.txt")
asfluxdat19 <- ozfluxdataextractrfas("ASM_EC_2019_L3_Corrected.nc", "as2019.txt")

yancfluxdat14 <- ozfluxdataextractrfyanc("Yanco_2014_L3.nc", "yanc2014.txt")
yancfluxdat15 <- ozfluxdataextractrfyanc("Yanco_2015_L3.nc", "yanc2015.txt")
yancfluxdat16 <- ozfluxdataextractrfyanc("Yanco_2016_L3.nc", "yanc2016.txt")
yancfluxdat17 <- ozfluxdataextractrfyanc("Yanco_2017_L3.nc", "yanc2017.txt")
yancfluxdat18 <- ozfluxdataextractrfyanc("Yanco_2018_L3.nc", "yanc2018.txt")
yancfluxdat19 <- ozfluxdataextractrfyanc("Yanco_2019_L3.nc", "yanc2019.txt")
yancfluxdat20 <- ozfluxdataextractrfyanc("Yanco_2020_L3.nc", "yanc2020.txt")

##zero plane displacement (d) and roughness length (z0)
wspas <- rbind(asfluxdat12[[4]], asfluxdat13[[4]], asfluxdat14[[4]], asfluxdat15[[4]], asfluxdat16[[4]], asfluxdat17[[4]], asfluxdat18[[4]], asfluxdat19[[4]])
tapas <- rbind(asfluxdat12[[5]], asfluxdat13[[5]], asfluxdat14[[5]], asfluxdat15[[5]], asfluxdat16[[5]], asfluxdat17[[5]], asfluxdat18[[5]], asfluxdat19[[5]])
ahpas <- rbind(asfluxdat12[[6]], asfluxdat13[[6]], asfluxdat14[[6]], asfluxdat15[[6]], asfluxdat16[[6]], asfluxdat17[[6]], asfluxdat18[[6]], asfluxdat19[[6]])
Et0datas <- rbind(asfluxdat12[[2]], asfluxdat13[[2]], asfluxdat14[[2]], asfluxdat15[[2]], asfluxdat16[[2]], asfluxdat17[[2]], asfluxdat18[[2]], asfluxdat19[[2]])
dtdzslpas <- c(asfluxdat12[[7]], asfluxdat13[[7]], asfluxdat14[[7]], asfluxdat15[[7]], asfluxdat16[[7]], asfluxdat17[[7]], asfluxdat18[[7]], asfluxdat19[[7]])
wspas <- rmNA(-9999, wspas)
tapas <- rmNA(-9999, tapas)
ahpas <- rmNA(-9999, ahpas)
Et0datas <- rmNA(-9999, Et0datas)
dtdzslpas <- rmNA(-9999, dtdzslpas)
nas <- length(dtdzslpas)
dudzas <- br439prf(nas, 6.5, wspas)
dTdzas <- br439prf(nas, 6.5, tapas)
dahdzas <- br439prf(nas, 6.5, ahpas)

Ridatas <- merge(dudzas, dTdzas,  by = c("zeval", "time"))
Ridatas <- merge(Ridatas, dahdzas,  by = c("zeval", "time"))
Ridatas <- Ridatas[,-c(4,5,7,8,10,11)]
names(Ridatas) <- c("zeval", "time", "dudz", "dTdz", "dahdz")
Richno <- with(Ridatas, 9.806/Et0datas$Tair/dudz^2*(dTdz+(0.61*Et0datas$Tair*(dahdz))))
hist(Richno[abs(Richno)<50])
dtdzpsqstbas <- pasqstability((dtdzslpas*100), Et0datas, Richno)

##energy balance
fltimas <- c(asfluxdat12[[1]], asfluxdat13[[1]], asfluxdat14[[1]], asfluxdat15[[1]], asfluxdat16[[1]], asfluxdat17[[1]], asfluxdat18[[1]], asfluxdat19[[1]])
ebalas <- enbal(Et0datas, fltimas, dtdzpsqstbas)
senh <- 11.6
canh <- 6.5
surtemp <- tapas[tapas$H == min(tapas$H, na.rm = T),]
surtemp <- surtemp[order(surtemp$time),2]
flrfas <- c(asfluxdat12[[3]], asfluxdat13[[3]], asfluxdat14[[3]], asfluxdat15[[3]], asfluxdat16[[3]], asfluxdat17[[3]], asfluxdat18[[3]], asfluxdat19[[3]])
flrfas <- rmNA(-9999, flrfas)
ra <- with(Et0datas, rhom*Cpm*(surtemp-Tair)/H)
dowas <- rep("dry", nas)
dowas[flrfas!=0] <- "wet"
tapply(ra, dowas, FUN= function(x)sum(is.na(x))/length(x)*100)

##Yanco
Et0daty <- rbind(yancfluxdat14[[2]], yancfluxdat15[[2]], yancfluxdat16[[2]], yancfluxdat17[[2]], yancfluxdat18[[2]], yancfluxdat19[[2]], yancfluxdat20[[2]])
wspy <- rbind(yancfluxdat14[[4]], yancfluxdat15[[4]], yancfluxdat16[[4]], yancfluxdat17[[4]], yancfluxdat18[[4]], yancfluxdat19[[4]], yancfluxdat20[[4]])
tapy <- rbind(yancfluxdat14[[5]], yancfluxdat15[[5]], yancfluxdat16[[5]], yancfluxdat17[[5]], yancfluxdat18[[5]], yancfluxdat19[[5]], yancfluxdat20[[5]])
ahpy <- rbind(yancfluxdat14[[6]], yancfluxdat15[[6]], yancfluxdat16[[6]], yancfluxdat17[[6]], yancfluxdat18[[6]], yancfluxdat19[[6]], yancfluxdat20[[6]])
dtdzslpy <- c(yancfluxdat14[[7]], yancfluxdat15[[7]], yancfluxdat16[[7]], yancfluxdat17[[7]], yancfluxdat18[[7]], yancfluxdat19[[7]], yancfluxdat20[[7]])
surtempy <- c(yancfluxdat14[[8]], yancfluxdat15[[8]], yancfluxdat16[[8]], yancfluxdat17[[8]], yancfluxdat18[[8]], yancfluxdat19[[8]], yancfluxdat20[[8]])
wspy <- rmNA(-9999, wspy)
tapy <- rmNA(-9999, tapy)
ahpy <- rmNA(-9999, ahpy)
Et0daty <- rmNA(-9999, Et0daty)
dtdzslpy <- rmNA(-9999, dtdzslpy)
ny <- length(dtdzslpy)
dudzy <- br439prf(ny, 0.3, wspy)
dTdzy <- br439prf(ny, 0.3, tapy)
dahdzy <- br439prf(ny, 0.3, ahpy)
Ridaty <- merge(dudzy, dTdzy,  by = c("zeval", "time"))
Ridaty <- merge(Ridaty, dahdzy,  by = c("zeval", "time"))
Ridaty <- Ridaty[,-c(4,5,7,8,10,11)]
names(Ridaty) <- c("zeval", "time", "dudz", "dTdz", "dahdz")
Richno <- with(Ridaty, 9.806/Et0daty$Tair/dudz^2*(dTdz+(0.61*Et0daty$Tair*(dahdz))))
hist(Richno[abs(Richno)<50])
dtdzpsqstby <- pasqstability((dtdzslpy*100), Et0daty, Richno)

##energy balance
fltimy <- c(yancfluxdat14[[1]], yancfluxdat15[[1]], yancfluxdat16[[1]], yancfluxdat17[[1]], yancfluxdat18[[1]], yancfluxdat19[[1]], yancfluxdat20[[1]])
ebaly <- enbal(Et0daty, fltimy, dtdzpsqstby)
senhy <- 8
canhy <- 0.3
flrfy <- c(yancfluxdat14[[3]], yancfluxdat15[[3]], yancfluxdat16[[3]], yancfluxdat17[[3]], yancfluxdat18[[3]], yancfluxdat19[[3]], yancfluxdat20[[3]])
flrfy <- rmNA(-9999, flrfy)
ra <- with(Et0daty, rhom*Cpm*(surtempy-Tair)/H)
dowy <- rep("dry", ny)
dowy[flrfy!=0] <- "wet"
tapply(ra, dowy, FUN= function(x)sum(is.na(x))/length(x)*100)

##whats up with L and dtvdz
alspmom <- momentum(Et0datas, wdpas,wspas,tapas,ahpas, dtdzpsqstbas, dudzas, dTdzas, dahdzas, 11.6, 6.5, F, plotloc = "H:/My Documents/Project and related work/Code/Papier Code", "Alice Springs", 0.01)
yancmom <- momentum(Et0daty, wdpy,wspy,tapy,ahpy, dtdzpsqstby, dudzy, dTdzy, dahdzy, 8, 0.3, F, plotloc = "H:/My Documents/Project and related work/Code/Papier Code", "Yanco", 0.01)

findfhdiv <- data.frame(d0 = c(ginmom[[7]]/7,alspmom[[7]]/6.5,yancmom[[7]]/0.3),z0=c(ginmom[[6]]/7,alspmom[[6]]/6.5,yancmom[[6]]/0.3), Site = c(rep("Gingin", length(ginmom[[7]])),rep("Alice Springs", length(alspmom[[7]])),rep("Yanco", length(yancmom[[7]]))))
findfhdiv$time <- as.POSIXct(c(Et0dat$time, Et0datas$time, Et0daty$time), origin = "1970-01-01 00:00:00")
findfhdiv$U <- c(Et0dat$u, Et0datas$u, Et0daty$u)
findfhdiv$ustar <- c(Et0dat$ustar, Et0datas$ustar, Et0daty$ustar)
findfhdiv$beta <- c(Et0dat$ustar/wsp[wsp$H == 7,2], Et0datas$ustar/wspas[wspas$H == 6.62,2], rep(NA, length(Et0daty$ustar)))
findfhdiv$WD <- c(wdp[wdp$H == 15.5,2], 
                  c(asfluxdat12[[9]], asfluxdat13[[9]], asfluxdat14[[9]], asfluxdat15[[9]], 
                    asfluxdat16[[9]], asfluxdat17[[9]], asfluxdat18[[9]], asfluxdat19[[9]]),
                  c(yancfluxdat14[[9]], yancfluxdat15[[9]], yancfluxdat16[[9]], 
                    yancfluxdat17[[9]], yancfluxdat18[[9]], yancfluxdat19[[9]], 
                    yancfluxdat20[[9]]))

setwd("H:/My Documents/Project and related work/Cool gingin/L3data")
USCMW <- amriflxprsc("AMF_US-CMW_BASE_HH_1-5.csv", 7, 14, 8, 0.01, 3)
USCMWcan <- USCMW[[2]]
USCMWlim <- USCMW[[3]]
USCMW <- USCMW[[1]]
USWCr <- amriflxprsc("AMF_US-WCr_BASE_HH_20-5.csv", 25, 29.6, 24.4, 0.01, c(12.2,2))
USWCrcan <- USWCr[[2]]
USWCrlim <- USWCr[[3]]
USWCr <- USWCr[[1]]
USxBr <- amriflxprsc("AMF_US-xBR_BASE_HH_4-5.csv", 23, 35.68, 25.48, 0.01, c(19.2,15.5))
USxBrcan <- USxBr[[2]]
USxBrlim <- USxBr[[3]]
USxBr <- USxBr[[1]]
USCMW$Site <- "Charleston"
USWCr$Site <- "Willow Creek"
USWCrog <- USWCr[USWCrlim=="Good",]
USWCrcan <- USWCrcan[USWCrlim=="Good",]
USWCrcan <- USWCrcan[USWCrog$ustar > 0.3,]
USWCr <- USWCrog[USWCrog$ustar > 0.3,]
USWCrcan <- USWCrcan[USWCr$WD < 90 | USWCr$WD > 180,]
USWCr <- USWCr[USWCr$WD < 90 | USWCr$WD > 180,]


USxBr$Site <- "Bartlett"
Amcomb <- rbind(USCMW, USWCr, USxBr)
Amcomb$z0 <- (1-Amcomb$d0)/exp(0.4/Amcomb$beta)
USPrrlst <- prrprsc()
USPrr <- USPrrlst[[1]]
USPrr$Site <- "Poker Flat"
USPrrcan <- USPrrlst[[2]]
USPrr <- USPrr[,c(1:18, 20, 19)]
Amcomb <- rbind(Amcomb, USPrr)
fincomb <- rbind(findfhdiv[,c(4,3,5,6,7,1,2,8)],Amcomb[,c(1,19,4,2,16,15,20,5)])

tags <- c("[337.5,22.5) North", "[22.5,67.5) North East","[67.5,112.5) East", "[112.5,157.5) South East",
          "[157.5,202.5) South", "[202.5,247.5) South West","[247.5,292.5) West", 
          "[292.5,337.5) North West")
try2grp <- fincomb %>%
  mutate(tag = case_when(
    WD >= 337.5 | WD < 22.5 ~ tags[1],
    WD >= 22.5 & WD < 67.5 ~ tags[2],
    WD >= 67.5 & WD < 112.5 ~ tags[3],
    WD >= 112.5 & WD < 157.5 ~ tags[4],
    WD >= 157.5 & WD < 202.5 ~ tags[5],
    WD >= 202.5 & WD < 247.5 ~ tags[6],
    WD >= 247.5 & WD < 292.5 ~ tags[7],
    WD >= 292.5 & WD < 337.5 ~ tags[8]))

try2grp$tag <- factor(try2grp$tag, levels = tags, ordered = FALSE)
try1grp <- try2grp[complete.cases(try2grp[!(is.na(try2grp$d0))&!(is.na(try2grp$WD)),]),]
g1 <- ggplot(data = try1grp, mapping = aes(x=tag,y=d0)) + facet_grid(Site~., scales = "free") +
      geom_jitter(color='coral',alpha=0.2) +
      geom_boxplot(fill="coral",color="black",alpha=0.3) +
      labs(y=expression(d[0]/H[c]), x = 'Direction') + 
      ylim(0,1) + 
      # geom_hline(data=cmpr[cmpr$estimator == "Mean",], aes(yintercept = d0))+
      # guides(color=FALSE) +
      theme_classic()
g1

try1grp <- try2grp[complete.cases(try2grp[!(is.na(try2grp$z0))&!(is.na(try2grp$WD)),]),]
g2 <- ggplot(data = try2grp, mapping = aes(x=tag,y=z0)) + facet_grid(Site~., scales = "free") +
  geom_jitter(color='cadetblue',alpha=0.2) +
  geom_boxplot(fill="cadetblue",color="black",alpha=0.3) +
  labs(y=expression(z[0]/H[c]), x = 'Direction') + ylim(0,1) + 
  # geom_hline(data=cmpr[cmpr$estimator == "Mean",], aes(yintercept = d0))+
  # guides(color=FALSE) +
  theme_classic()
g2

fincomb <- fincomb %>%
  mutate(cnh = case_when(
    Site == "Gingin" ~ 7,
    Site == "Alice Springs" ~ 6.5,
    Site == "Yanco" ~ 0.3,
    Site == "Charleston" ~ 7,#https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek" ~ 24, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 23, #https://www.neonscience.org/field-sites/bart
    Site == "Poker Flat" ~ 3 #https://www.sciencedirect.com/science/article/pii/S187396521300008X#fig2
  ))
##get bartlett info
fincomb <- fincomb %>%
  mutate(LAI = case_when(
    Site == "Gingin" ~ 0.7,
    Site == "Alice Springs" ~ 0.3,
    Site == "Yanco" ~ 0,
    Site == "Charleston" ~ 1.3, #https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek" ~ 5.3, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 3,
    Site == "Poker Flat" ~ 0.73 #https://ieeexplore.ieee.org/document/6589947
  ))
fincomb <- fincomb <- fincomb %>%
  mutate(techn = case_when(
    Site == "Gingin" ~ "nls",
    Site == "Alice Springs" ~ "nls",
    Site == "Yanco" ~ "nls",
    Site == "Charleston" ~ "numerd", #https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek" ~ "numerd", #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ "numerd",
    Site == "Poker Flat" ~ "nls" #https://ieeexplore.ieee.org/document/6589947
  ))

##density plot
fincomb$d0t <- with(fincomb,d0*cnh)
fincomb$z0t <- with(fincomb,z0*cnh)
fincomb$limit <- NA
fincomb$limit[fincomb$d0t == (fincomb$cnh-0.01)] <- "Upper"
fincomb$limit[fincomb$d0t == 0.005] <- "Lower"
fincomb$limit[fincomb$Site == "Charleston"] <- USCMWlim
fincomb$limit[fincomb$Site == "Willow Creek"] <- NA
fincomb$limit[fincomb$Site == "Bartlett"] <- USxBrlim
fincomb$Season <- "Leaf on"
fincomb$Season[month(fincomb$time) > 11 | month(fincomb$time) < 6] <- "Leaf off"
fincomb$Season[fincomb$Site == "Gingin" | fincomb$Site == "Alice Springs" | fincomb$Site == "Poker Flat" | fincomb$Site == "Yanco"] <- "Evergreen"
fincomb$Site[fincomb$Site == "Willow Creek" & year(fincomb$time) < 2007] <- "Willow Creek < 2007"
fincomb$Site[fincomb$Site == "Willow Creek" & year(fincomb$time) > 2007] <- "Willow Creek > 2010"
ggplot(fincomb[!is.na(fincomb$d0),], aes(U)) + facet_wrap(.~Site, scales = "free") + geom_density(aes(col = Season, fill = Season), alpha = 0.3) + scale_fill_manual(values = c("green", "burlywood4", "forestgreen")) + scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + theme_classic() + xlab("U [m/s]")
ggplot(fincomb[!is.na(fincomb$d0),], aes(ustar)) + facet_wrap(.~Site, scales = "free") + geom_density(aes(col = Season, fill = Season), alpha = 0.3) + scale_fill_manual(values = c("green", "burlywood4", "forestgreen")) + scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + theme_classic() + xlab("u* [m/s]")
ggplot(fincomb[!is.na(fincomb$d0),], aes(beta)) + facet_wrap(.~Site, scales = "free") + geom_density(aes(col = Season, fill = Season), alpha = 0.3) + scale_fill_manual(values = c("green", "burlywood4", "forestgreen")) + scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + theme_classic() + xlab(expression(beta)) + geom_vline(xintercept = 1/3)

load("//uniwa.uwa.edu.au/userhome/students9/22789059/My Documents/Project and related work/Cool gingin/L3data/.RData")
canhsit <- tapply(fincomb$cnh, fincomb$Site, FUN = function(x) mean(x, na.rm = T))
LAIsit <- tapply(fincomb$LAI, fincomb$Site, FUN = function(x) mean(x, na.rm = T))
d0stnd <- as.data.frame(tapply(fincomb$d0[is.na(fincomb$limit) | fincomb$limit == "Good"], list(fincomb$Site[is.na(fincomb$limit) | fincomb$limit == "Good"], fincomb$Season[is.na(fincomb$limit) | fincomb$limit == "Good"]), 
                        FUN = function(x) mean(x, na.rm = T)))
d0stnd$Site <- row.names(d0stnd)
d0stnd$cnh <- canhsit
d0stnd <- melt(d0stnd, id = c("Site", "cnh"))
d0stnd <- d0stnd[complete.cases(d0stnd),]

z0stnd <- as.data.frame(tapply(fincomb$z0[is.na(fincomb$limit) | fincomb$limit == "Good"], list(fincomb$Site[is.na(fincomb$limit) | fincomb$limit == "Good"], fincomb$Season[is.na(fincomb$limit) | fincomb$limit == "Good"]), 
                               FUN = function(x) mean(x, na.rm = T)))
z0stnd$Site <- row.names(z0stnd)
z0stnd$cnh <- canhsit
z0stnd$LAI <- c(0.3,0.7,0.73,0.5,3.6,1,5.3,5.3,4.5,1.6,5.3,5.3)
z0stnd <- melt(z0stnd, id = c("Site", "cnh"))
z0stnd <- z0stnd[complete.cases(z0stnd),]

cmpr <- data.frame(Site = c(rep(names(canhsit), 2), d0stnd$Site), 
                   estimator = c(rep(c("Convention", "Raupach"),each = length(canhsit)), paste("Mean-", d0stnd$variable)),
                   d0=c((2/3*canhsit), (1-(1-exp(-sqrt(7.5*LAIsit*2)))/sqrt(7.5*LAIsit*2))*canhsit, d0stnd$value*d0stnd$cnh),
                   z0=c((0.12*canhsit), ((1-exp(-sqrt(7.5*LAIsit*2)))/sqrt(7.5*LAIsit*2))*exp(-(0.4/0.3)-0.193)*canhsit,z0stnd$value*z0stnd$cnh))

write.csv(cmpr, "d0z0vals.csv")
cmprhdiv <- data.frame(Site = rep(names(d0stnd), 3), estimator = rep(c("Convention", "Raupach", "Mean"),each = length(d0stnd)),  
                       d0=c(rep((2/3),length(d0stnd)), (1-(1-exp(-sqrt(7.5*LAIsit*2)))/sqrt(7.5*LAIsit*2)), d0stnd),
                       z0=c(rep((0.12),length(d0stnd)), ((1-exp(-sqrt(7.5*LAIsit*2)))/sqrt(7.5*LAIsit*2))*exp(-(0.4/0.3)-0.193), z0stnd))

gtop <- ggplot(fincomb[(is.na(fincomb$limit) | fincomb$limit == "Good")& (fincomb$Site == "Gingin" | fincomb$Site == "Alice Springs" | fincomb$Site == "Charleston"),], aes(d0t)) + 
  geom_density(aes(fill = Season), alpha = 0.4) + 
  facet_wrap(.~Site, scales = "free")+ theme_classic() + geom_vline(aes(xintercept = cnh), col = "black", alpha = 0.2) + 
  xlab(expression(d[0] ~~ "[m]")) + coord_flip()  +
  scale_fill_manual(values = c("forestgreen", "burlywood4", "darkgreen"))
gtop + geom_rug(data = cmpr[(cmpr$Site == "Gingin" | cmpr$Site == "Alice Springs" | cmpr$Site == "Charleston"),], aes(x = d0, col = estimator), size = 1, length = unit(0.1, "npc")) +
  scale_color_manual(name = element_blank(), values = c("deepskyblue", "forestgreen", "burlywood4", "darkgreen", "slateblue4"))  + theme(strip.background = element_blank())

gbtm <- ggplot(fincomb[(is.na(fincomb$limit) | fincomb$limit == "Good")& (fincomb$Site == "Willow Creek < 2007" | fincomb$Site == "Willow Creek > 2010" | fincomb$Site == "Bartlett"),], aes(d0t)) + 
  geom_density(aes(fill = Season), alpha = 0.4) + 
  facet_wrap(.~Site, scales = "free")+ theme_classic() + geom_vline(aes(xintercept = cnh), col = "black", alpha = 0.2) + 
  xlab(expression(d[0] ~~ "[m]")) + coord_flip() + 
  scale_fill_manual(values = c("burlywood4", "darkgreen"))
gbtm + geom_rug(data = cmpr[(cmpr$Site == "Willow Creek < 2007" | cmpr$Site == "Willow Creek > 2010" | cmpr$Site == "Bartlett"),], aes(x = d0, col = estimator), size = 1, length = unit(0.1, "npc")) +
  scale_color_manual(name = element_blank(), values = c("deepskyblue", "burlywood4","darkgreen",  "slateblue4"))  + theme(strip.background = element_blank())

gmid <- ggplot(fincomb[is.na(fincomb$limit) & (fincomb$Site == "Poker Flat"),], aes(d0t)) + geom_density(aes(fill = Season), alpha = 0.5) +
  facet_wrap(.~Site)+ theme_classic() + geom_vline(aes(xintercept = cnh), col = "black", alpha = 0.2) +
  xlab(expression(d[0] ~~ "[m]")) + coord_flip() +
  scale_fill_manual(values = c("forestgreen"))
gmid + geom_rug(data = cmpr[(cmpr$Site == "Poker Flat"),], aes(x = d0, col = estimator), size = 1, length = unit(0.1, "npc")) +
  scale_color_manual(name = element_blank(), values = c("deepskyblue", "forestgreen","slateblue4"))  + theme(strip.background = element_blank())

gmid2 <- ggplot(fincomb[is.na(fincomb$limit) & (fincomb$Site == "Yanco"),], aes(d0t)) + geom_density(aes(fill = Season), alpha = 0.5) +
  facet_wrap(.~Site)+ theme_classic() + geom_vline(aes(xintercept = cnh), col = "black", alpha = 0.2) +
  xlab(expression(d[0] ~~ "[m]")) + coord_flip() +
  scale_fill_manual(values = c("forestgreen"))
gmid2 + geom_rug(data = cmpr[(cmpr$Site == "Yanco"),], aes(x = d0, col = estimator), size = 1, length = unit(0.1, "npc")) +
  scale_y_continuous(breaks = c(0,(7.5*10^13)), labels = c(0, expression(paste(7.5,x,10^13)))) +
  scale_color_manual(name = element_blank(), values = c("deepskyblue", "forestgreen","slateblue4"))  + theme(strip.background = element_blank())

ggplot(z0stnd, aes(LAI, value, label = Site))   + geom_point(mapping = aes(x=0,y=0)) + ylim(0,.3) + theme_classic() + xlim(0,4/0.6) + 
  stat_function(fun = function(x)(((20^-(LAI))*(6-(LAI)))*exp(-(0.5*0.3/0.16*((20^-(LAI))*(6-(LAI))*20))^-0.5))) + 
  ylab(expression(frac(z[0],H[c]))) + geom_text(vjust = 0.1, hjust = 0.05, nudge_x = 0.01,angle = 45)+ geom_point(aes(col = variable)) + 
  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) 


ggplot(z0stnd, aes(LAI, value, label = Site))   + geom_point(mapping = aes(x=0,y=0)) + ylim(0,.3) + theme_classic() + xlim(0,4/0.6) + 
  #stat_function(fun = function(x)(((20^-(LAI))*(6-(LAI)))*exp(-(0.5*0.3/0.16*((20^-(LAI))*(6-(LAI))*20))^-0.5))) + 
  ylab(expression(frac(z[0],H[c]))) + #geom_text(vjust = 0.1, hjust = 0.05, nudge_x = 0.01,angle = 45)+ 
  geom_point(aes(col = variable)) + 
  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) 

dectry2 <- nls(value ~ ((a^-(LAI))*((LAI)))*exp(-1/sqrt(0.5*b/0.16*((a^-(LAI))*(LAI)))), data = z0stnd[z0stnd$variable!="Evergreen",],
  start = list(a = 2, b = 1.6), trace = T, na.action(na.exclude(1)),
  algorithm = "port", lower = c(1, 0.03), upper = c(20, 4),
  control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
evgrtry2 <- nls(value ~ ((a^-(LAI))*((LAI)))*exp(-1/sqrt(0.5*b/0.16*((a^-(LAI))*(LAI)))), data = z0stnd[z0stnd$variable=="Evergreen",],
  start = list(a = 2, b = 1.6), trace = T, na.action(na.exclude(1)),
  algorithm = "port", lower = c(1, 0.03), upper = c(20, 4),
  control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
try2 <- nls(value ~ ((a^-(LAI))*((LAI)))*exp(-1/sqrt(0.5*b/0.16*((a^-(LAI))*(LAI)))), data = z0stnd,
  start = list(a = 2, b = 1.6), trace = T, na.action(na.exclude(1)),
  algorithm = "port", lower = c(1, 0.03), upper = c(20, 4),
  control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = T))


ggplot(z0stnd, aes(LAI, value, label = Site))   + geom_point(mapping = aes(x=0,y=0)) + ylim(0,.3) + theme_classic() + xlim(0,4/0.6) + 
    stat_function(fun = function(x)((2^-(x))*((x))*exp(-1/sqrt(0.5*1.2/(0.4^2)*((2^-(x))*((x)))))), colour = "darkgreen", linetype = 2, size = 1.2) +  
  stat_function(fun = function(x)((1.46396^-(x))*((x))*exp(-1/sqrt(0.5*0.0888711/(0.4^2)*((1.46396^-(x))*((x)))))), colour = "burlywood4", linetype = 2, size = 1.2) +  
  ylab(expression(frac(z["0m"],H[c]))) + geom_point(aes(fill = variable), shape =21, size = 4) + 
  scale_fill_manual(name = "Season",values = c("forestgreen", "burlywood4", "darkgreen")) 

ggplot(fincomb[is.na(fincomb$limit),], aes(z0t)) + geom_density(aes(fill = Season), alpha = 0.5) + xlim(0,5) + 
  facet_wrap(.~Site, scales = "free")+ theme_classic() + geom_vline(aes(xintercept = cnh), col = "black", alpha = 0.2) + 
  xlab(expression(z[0] ~~ "[m]")) + geom_vline(data = cmpr, aes(xintercept = z0, col = estimator), size = 1 ) + coord_flip() + scale_color_manual(values = c("deepskyblue", "darkgreen", "burlywood4", "forestgreen", "slateblue4")) + scale_fill_manual(values = c("green", "burlywood4", "forestgreen"))
ggplot(fincomb[is.na(fincomb$limit),], aes(z0t)) + geom_density(col = "darkslategray",fill = "darkslategray", alpha = 0.2) + 
  facet_wrap(.~Site, scales = "free")+ theme_classic() +  
  xlab(expression(z[0])) + geom_vline(data = cmpr, aes(xintercept = z0, col = estimator)) + coord_flip()
# ggplot(fincomb, aes(d0)) + geom_density(col = "red",fill = "red", alpha = 0.2) + facet_wrap(.~Site, scales = "free_x")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab(paste(expression(d[0]), "/", expression(H[c]))) + geom_vline(data = cmprhdiv, aes(xintercept = d0, col = estimator))
# ggplot(fincomb, aes(z0t)) + geom_density(col = "red", fill = "red", alpha = 0.2) + facet_wrap(.~Site, scales = "free")+ 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   xlab(paste(expression(z[0]), "m"))+ geom_vline(data = cmpr, aes(xintercept = z0, col = estimator))
# ggplot(fincomb, aes(z0)) + geom_density(col = "red",fill = "red", alpha = 0.2) + facet_wrap(.~Site, scales = "free_x")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab(paste(expression(z[0]), "/", expression(H[c])))  + geom_vline(data = cmprhdiv, aes(xintercept = z0, col = estimator))

prf <- fincomb[fincomb$techn == "nls",]
prf <- prf[complete.cases(prf),]
prf <- prf[prf$d0t != 0.005,]
prf$z0oth <- (1-prf$d0)/exp(0.4/prf$beta)
ggplot(prf, aes(z0, z0oth)) + geom_point(shape = 1, col = "grey", alpha = 0.2) + facet_wrap(.~Site) + theme_bw() + geom_abline(intercept = 0, slope = 1) + geom_smooth(method = "lm") + xlim(0,0.5) + ylim(0,0.5)
ggplot(fincomb[fincomb$Site!="Yanco" & fincomb$d0t != 0.005,], aes(d0, z0)) + geom_point(col = "grey", alpha = 0.1) + facet_wrap(techn~Site, scales = "free")+ theme_bw() + geom_smooth(method = "lm")

ggplot(fincomb[is.na(fincomb$limit) & !is.na(fincomb$d0),], aes(U)) + geom_density(col = "forestgreen",fill = "forestgreen", alpha = 0.2) + 
  facet_wrap(Site~., scales = "free")+ theme_classic() 

ggplot(fincomb[is.na(fincomb$limit) & !is.na(fincomb$d0),], aes(ustar)) + geom_density(col = "firebrick4",fill = "firebrick4", alpha = 0.2) + 
  facet_wrap(Site~., scales = "free")+ theme_classic() 

top <- data.frame()
for(i in 1:length(unique(fincomb$Site))){
  sub <- fincomb[fincomb$Site == unique(fincomb$Site)[i],]
  ans <- quantile(sub$beta[sub$beta != Inf], probs = 0.98, na.rm=T)
  sub$top <- T
  sub$top[sub$beta > ans] <- F
  sub$ovrlim <- T
  sub$ovrlim[sub$beta > 0.3333] <- F
  top <- rbind(top, sub)
}
top <- top[top$Site != "Yanco", ]

ggplot(top[is.na(top$limit) & top$top,], aes(beta)) + geom_density(aes(fill = Season), alpha = 0.5) + 
  facet_wrap(.~Site, scales = "free")+ theme_classic() + geom_vline(xintercept = 0.333) + 
  xlab(expression(beta)) + scale_fill_manual(values = c("green", "burlywood4", "forestgreen"))

##conventional results table
set.seed(338)
calval <- sample(rep(c("Cal", "Val"), nrow(fincomb) * c(0.5, 0.5)))
fincomb2 <- fincomb[is.na(fincomb$limit) & (calval == "Cal"), -14]
dzlst2 <- list()
udlst2 <- list()
ustrdlst2 <- list()
uzlst2 <- list()
ustrzlst2 <- list()
count <- 1
st <- NA
sn <- NA
for(i in 1:length(unique(fincomb2$Site))){
  sub <- fincomb2[complete.cases(fincomb2) & (fincomb2$Site == unique(fincomb2$Site)[i]) ,]
  for(j in 1:length(unique(sub$Season))){
    sub2 <- sub[(sub$Season == unique(sub$Season)[j]) ,]
    print(sub2[1,])
    if(nrow(sub2) != 0){
      dzlst2[[count]] <- with(sub2, lm(d0t~z0t))
      udlst2[[count]] <- with(sub2, lm(d0t~U))
      ustrdlst2[[count]] <- with(sub2, lm(d0t~ustar))
      uzlst2[[count]] <- with(sub2, lm(z0t~U))
      ustrzlst2[[count]] <- with(sub2, lm(z0t~ustar))
      st[count] <- unique(fincomb2$Site)[i]
      sn[count] <- unique(sub$Season)[j]
      count <- count+1
    }else{
      st[count] <- NA
      sn[count] <- NA
      count <- count+1
    }
  }
}
cmblst2 <- list(dzlst2, udlst2, ustrdlst2, uzlst2, ustrzlst2)
tbfstf2 <- data.frame(prvr = c(rep("d0", (count-1)*3), 
                               rep("z0", (count-1)*2)), 
                     exvr = c(rep("z0", (count-1)), 
                              rep(rep(c("U", "ustar"), each = (count-1)),2)),
                     Site = rep(st, 5), Season = rep(sn, 5),
                     Slope = NA, R2 = NA, RMSE = NA , Intercept = NA)
for(i in 1:5){
  for(j in 1:(count-1)){
    if(j != 3 & j != 4){
      tbfstf2$Slope[(i-1)*(count-1)+j] <- cmblst2[[i]][[j]]$coefficients[2]
      tbfstf2$Intercept[(i-1)*(count-1)+j] <- cmblst2[[i]][[j]]$coefficients[1]
      tbfstf2$R2[(i-1)*(count-1)+j] <- summary(cmblst2[[i]][[j]])$r.squared
      tbfstf2$RMSE[(i-1)*(count-1)+j] <- summary(cmblst2[[i]][[j]])$sigma
    }
  }
}
write.csv(tbfstf2, "conventionalexplanation.csv")
##finally putting an lc thing
Ind0lc <- with(fincomb, 2*beta^3/0.4/(1-d0))
Inz0lc <- with(fincomb, 2*beta^3/0.4/z0/exp(0.4/beta))
cw68d0lc <- function(x,y,par){
  sum((2*x/(1 - y)/0.4)-(sqrt(par[1]/x/0.4/(1 - y))/tanh(sqrt(par[1]/x/0.4/(1 - y)))))^2
}
cw68z0lc <- function(x,y,par){
  sum((2*x/y/0.4/exp(0.4/x))-(sqrt(par[1]/x/0.4/y/exp(0.4/x))/tanh(sqrt(par[1]/x/0.4/y/exp(0.4/x)))))^2
}
Cwd0lc <- 0
for(i in 1:length(fincomb$d0)){
  if(sum(is.na(c(fincomb$beta[i], fincomb$d0[i])))==0){
    results <- optim(par = 1, fn = cw68d0lc, x = fincomb$beta[i], y =fincomb$d0[i], method = "Brent", lower = 0.0000000000001, upper = 500000)
    Cwd0lc[i] <- results$par
  }else{
    Cwd0lc[i] <- NA
  }
}
Cwz0lc <- 0
for(i in 1:length(fincomb$z0)){
  if(sum(is.na(c(fincomb$beta[i], fincomb$z0[i])))==0){
    results <- optim(par = 0.4, fn = cw68z0lc, x = fincomb$beta[i], y =fincomb$z0[i], method = "Brent", lower = 0.0000000000001, upper = 500000)
    Cwz0lc[i] <- results$par
  }else{
    Cwz0lc[i] <- NA
  }
}
Thd0lc <- with(fincomb, 6*beta^3/(4*0.4*(1-d0)))
Thz0lc <- with(fincomb, 6*beta^3/z0/(4*0.4)/exp(0.4/beta))
wg12d0lc <- function(x,y,par){
  sum(1 - y - (x*(besselI(2*sqrt(par[1]/0.4/x),0) + 
                    (besselI(2*sqrt(0.000001/x/0.4*par[1]),0)*
                       besselK(2*sqrt(par[1]/x/0.4),0)/
                       besselK(2*sqrt(0.000001*par[1]/x/0.4),0)))/
                 (0.4*sqrt(par[1]/x/0.4)*
                    (besselI(2*sqrt(par[1]/x/0.4),1) - 
                       (besselI(2*sqrt(0.000001/x/0.4*par[1]),0)*
                          besselK(2*sqrt(par[1]/x/0.4),1)/
                          besselK(2*sqrt(0.000001/x/0.4*par[1]),0))))))^2}
wg12z0lc <- function(x,y,par){
  sum(((x*(besselI(2*sqrt(par[1]/0.4/x),0) + 
             (besselI(2*sqrt(0.000001/x/0.4*par[1]),0)*
                besselK(2*sqrt(par[1]/x/0.4),0)/
                besselK(2*sqrt(0.000001*par[1]/x/0.4),0)))/
          (0.4*sqrt(par[1]/x/0.4)*
             (besselI(2*sqrt(par[1]/x/0.4),1) - 
                (besselI(2*sqrt(0.000001/x/0.4*par[1]),0)*
                   besselK(2*sqrt(par[1]/x/0.4),1)/
                   besselK(2*sqrt(0.000001/x/0.4*par[1]),0)))))/exp(0.4/x))-y)^2}
Wgd0lc <- 0
for(i in 1:length(fincomb$d0)){
  if(sum(is.na(c(fincomb$beta[i], fincomb$d0[i])))==0){
    results <- optim(par = 1, fn = wg12d0lc, x = fincomb$beta[i], y =fincomb$d0[i], method = "Brent", lower = 0.0000000000001, upper = 10000)
    Wgd0lc[i] <- results$par
  }else{
    Wgd0lc[i] <- NA
  }
}
Wgz0lc <- 0
for(i in 1:length(fincomb$z0)){
  if(sum(is.na(c(fincomb$beta[i], fincomb$z0[i])))==0){
    results <- optim(par = 0.4, fn = wg12z0lc, x = fincomb$beta[i], y =fincomb$z0[i], method = "Brent", lower = 0.0000000000001, upper = 10000)
    Wgz0lc[i] <- results$par
  }else{
    Wgz0lc[i] <- NA
  }
}
toplcdist <- as.data.frame(cbind(fincomb[,c(6,7,2,5,1,15)], calval))
mtoplcdist <- melt(toplcdist, id = c("Site", "beta", "time", "Season", "calval"))
mtoplcdist$Inoue <- c(Ind0lc, Inz0lc)
mtoplcdist$Cowan <- c(Cwd0lc, Cwz0lc)
mtoplcdist$Thom <- c(Thd0lc, Thz0lc)
mtoplcdist$Wang <- c(Wgd0lc, Wgz0lc)
names(mtoplcdist)[6:7] <- c("advarbyhc","valadvar")
mmtoplcdist <- melt(mtoplcdist, id = c("Site", "beta", "time", "Season", "calval", "advarbyhc","valadvar"))
mmtoplcdistd0 <- mmtoplcdist[mmtoplcdist$advarbyhc == "d0",]

betovr <- as.data.frame(tapply(toplcdist$beta[!is.na(toplcdist$d0)], list(toplcdist$Site[!is.na(toplcdist$d0)], toplcdist$Season[!is.na(toplcdist$d0)]), 
                               FUN = function(x) sum(x<1/3, na.rm = T)/sum(!is.na(x))))
betovr$Site <- row.names(betovr)
betovr <- melt(betovr, id = "Site")
betovr <- betovr[complete.cases(betovr),]
betmn <- as.data.frame(tapply(toplcdist$beta[!is.na(toplcdist$d0)], list(toplcdist$Site[!is.na(toplcdist$d0)], toplcdist$Season[!is.na(toplcdist$d0)]), 
                               FUN = function(x) mean(x[x<13.7],na.rm = T)))
betmn$Site <- row.names(betmn)
betmn <- melt(betmn, id = "Site")
betmn <- betmn[complete.cases(betmn),]
betsd <- as.data.frame(tapply(toplcdist$beta[!is.na(toplcdist$d0)], list(toplcdist$Site[!is.na(toplcdist$d0)], toplcdist$Season[!is.na(toplcdist$d0)]), 
                              FUN = function(x) sd(x[x<13.7],na.rm = T)))
betsd$Site <- row.names(betsd)
betsd <- melt(betsd, id = "Site")
betsd <- betsd[complete.cases(betsd),]

d0sd <- as.data.frame(tapply((fincomb$d0*fincomb$cnh)[is.na(fincomb$limit)], list(fincomb$Site[is.na(fincomb$limit)], fincomb$Season[is.na(fincomb$limit)]), 
                             FUN = function(x) sd(x, na.rm = T)))
d0sd$Site <- row.names(d0sd)
d0sd <- melt(d0sd, "Site")
d0sd <- d0sd[complete.cases(d0sd),]

z0sd <- as.data.frame(tapply((fincomb$z0*fincomb$cnh)[is.na(fincomb$limit)], list(fincomb$Site[is.na(fincomb$limit)], fincomb$Season[is.na(fincomb$limit)]), 
                             FUN = function(x) sd(x, na.rm = T)))
z0sd$Site <- row.names(z0sd)
z0sd <- melt(z0sd, "Site")
z0sd <- z0sd[complete.cases(z0sd),]

sdvals <- data.frame(Site = d0sd$Site, Season = d0sd$variable, d0sd = d0sd$value, z0sd = z0sd$value)
betvals<- data.frame(Site = betmn$Site, Season = betmn$variable, betmn = betmn$value, betsd = betsd$value, betovr = betovr$value)
write.csv(sdvals,"d0z0sd.csv")
write.csv(betvals,"betvals.csv")

z0err <- z0stnd
z0err$lower <- z0err$value - z0sd$value/2/z0err$cnh
z0err$upper <- z0err$value + z0sd$value/2/z0err$cnh
ggplot(z0err, aes(LAI, value, label = Site))   + geom_point(mapping = aes(x=0,y=0)) + ylim(0,.3) + theme_classic() + xlim(0,4/0.6) + 
  stat_function(fun = function(x)((2^-(x))*((x))*exp(-1/sqrt(0.5*1.2/(0.4^2)*((2^-(x))*((x)))))), colour = "darkgreen", linetype = 2, size = 1.2) +  
  stat_function(fun = function(x)((1.46396^-(x))*((x))*exp(-1/sqrt(0.5*0.0888711/(0.4^2)*((1.46396^-(x))*((x)))))), colour = "burlywood4", linetype = 2, size = 1.2) +  
  ylab(expression(frac(z["0m"],H[c]))) + geom_point(aes(fill = variable), shape =21, size = 4) + 
  scale_fill_manual(name = "Season",values = c("forestgreen", "burlywood4", "darkgreen")) + geom_errorbar(aes(ymin = lower, ymax = upper, col = variable), size = 1) + 
  scale_color_manual(name = "Season",values = c("forestgreen", "burlywood4", "darkgreen"))


ggplot(mmtoplcdist[complete.cases(mmtoplcdist),], aes(value)) + facet_grid(Site~variable, scales = "free") + 
  geom_density(aes(col = advarbyhc, fill = advarbyhc), alpha = 0.2) + scale_x_continuous(trans = 'log10') + 
  theme_classic() + xlab(expression(c[d]~aH[c]))

wndmods <- c("Inoue", "Cowan", "Thom", "Wang")
mmtoplcdistd0 <- mmtoplcdistd0[(mmtoplcdistd0$Site == "Gingin" | mmtoplcdistd0$Site == "Charleston") & 
                               (mmtoplcdistd0$variable == "Inoue" | mmtoplcdistd0$variable == "Wang") & complete.cases(mmtoplcdistd0),]
library(scales)
library(viridis)
##Figure 7 final
ggplot(mmtoplcdistd0[complete.cases(mmtoplcdistd0),], aes(value)) + facet_grid(Site~variable, scales = "free_y") + 
  geom_density(aes(fill = Season), alpha = 0.4) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + xlab(expression(c[d]~aH[c])) + scale_fill_manual(values = c("forestgreen", "burlywood4", "darkgreen")) + theme(strip.background = element_blank())
##Figure 10
ggplot(mmtoplcdist[mmtoplcdist$value >10^-4 & mmtoplcdist$value < 10^4 & complete.cases(mmtoplcdist) & mmtoplcdist$advarbyhc == "d0" &
                     (mmtoplcdist$variable == "Inoue" | mmtoplcdist$variable == "Wang"),], aes(value, beta)) + 
  facet_grid(variable~.) + geom_point(size = 0.6, shape = 1, aes(col = valadvar)) + theme_classic() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + 
  ylim(0,2) + scale_color_viridis(name = expression(frac(d[0],H[c]))) + xlab(expression(c[d]~aH[c])) + ylab(expression(beta))   + theme(strip.background = element_blank())

ggplot(mmtoplcdistd0[complete.cases(mmtoplcdistd0),], aes(value, beta)) + facet_grid(Site~variable, scales = "free") + 
  geom_point(aes(col = Season, fill = Season), alpha = 0.2, size=0.2) + scale_x_continuous(trans = 'log10') + 
  theme_classic() + xlab(expression(c[d]~aH[c])) + ylab(expression(beta))
trmdn <- tapply(mmtoplcdist$value[mmtoplcdist$calval == "Cal" & mmtoplcdist$advarbyhc == "d0"], 
                list(mmtoplcdist$Site[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"], 
                     mmtoplcdist$variable[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"], 
                     mmtoplcdist$Season[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"]), 
                FUN = function(x)(median(x,na.rm = T)))
tag <- as.data.frame(rbind(trmdn[,,1],trmdn[,,2], trmdn[,,3]))
tag$Site <- rownames(trmdn[,,1])
tag$Season <- rep(c("Evergreen","Leaf off", "Leaf on"), each= 8)
mtag <- melt(tag, id = c("Site", "Season"))
mtag <- mtag[complete.cases(mtag),]

trmn <- tapply(mmtoplcdist$value[mmtoplcdist$calval == "Cal" & mmtoplcdist$advarbyhc == "d0"], 
                list(mmtoplcdist$Site[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"], 
                     mmtoplcdist$variable[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"], 
                     mmtoplcdist$Season[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"]), 
                FUN = function(x)(mean(x[x != Inf],na.rm = T)))
tag2 <- as.data.frame(rbind(trmn[,,1],trmn[,,2], trmn[,,3]))
tag2$Site <- rownames(trmn[,,1])
tag2$Season <- rep(c("Evergreen","Leaf off", "Leaf on"), each= 8)
mtag2 <- melt(tag2, id = c("Site", "Season"))
mtag2 <- mtag2[complete.cases(mtag2),]

trsd <- tapply(mmtoplcdist$value[mmtoplcdist$calval == "Cal" & mmtoplcdist$advarbyhc == "d0"], 
               list(mmtoplcdist$Site[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"], 
                    mmtoplcdist$variable[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"], 
                    mmtoplcdist$Season[mmtoplcdist$calval == "Cal"& mmtoplcdist$advarbyhc == "d0"]), 
               FUN = function(x)(sd(x[x != Inf],na.rm = T)))
tag3 <- as.data.frame(rbind(trsd[,,1],trsd[,,2], trsd[,,3]))
tag3$Site <- rownames(trsd[,,1])
tag3$Season <- rep(c("Evergreen","Leaf off", "Leaf on"), each= 8)
mtag3 <- melt(tag3, id = c("Site", "Season"))
mtag3 <- mtag3[complete.cases(mtag3),]

lctab <- data.frame(Site = mtag$Site, Season = mtag$Season, Model = mtag$variable,
                    lcmdn = mtag$value, lcmn = mtag2$value, lcsd = mtag3$value)
write.csv(lctab, "lctbl.csv")

mtagtop <- mtag[(mtag$Site == "Gingin" |  mtag$Site == "Charleston") & (mtag$variable == "Inoue" |  mtag$variable == "Wang"),]
ggplot(mmtoplcdistd0, aes(value)) + facet_grid(Site~variable, scales = "free") + 
  geom_density(aes(fill = Season), alpha = 0.5) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(values = c("green", "burlywood4", "forestgreen")) + 
  theme(strip.background = element_blank()) + 
  geom_vline(data = mtagtop, aes(xintercept = value, col = Season)) + scale_color_manual(values = c("green", "burlywood4", "forestgreen"))


mtagtop2 <- mtag[mtag$Site == "Willow Creek 2" | mtag$Site == "Bartlett",]
top2 <-  mmtoplcdist[mmtoplcdist$Site == "Bartlett",]
top2 <- top2[year(top2$time) == 2019,]
top2$plot <- F
top2$plot[(top2$value > 10^-3) & (top2$value < 10)] <- T
top2$bimon <- "Jan-Feb"
top2$bimon[month(top2$time)>=3 & month(top2$time)<=4] <- "Mar-Apr"
top2$bimon[month(top2$time)>=5 & month(top2$time)<=6] <- "May-Jun"
top2$bimon[month(top2$time)>=7 & month(top2$time)<=8] <- "Jul-Aug"
top2$bimon[month(top2$time)>=9 & month(top2$time)<=10] <- "Sep-Oct"
top2$bimon[month(top2$time)>=11 & month(top2$time)<=12] <- "Nov-Dec"
top2$bimon <- factor(top2$bimon, levels = c("Jan-Feb", "Mar-Apr", "May-Jun", "Jul-Aug", "Sep-Oct", "Nov-Dec"))
top2$plot[(top2$variable == "Cowan" | top2$variable == "Thom") & top2$plot] <- F
##Figure 8
ggplot(top2[complete.cases(top2) & top2$plot,], aes(value)) + facet_grid(bimon~.) + ggtitle("Bartlett 2019") +
  geom_density(aes(fill = variable), alpha = 0.5) + scale_x_continuous(trans = 'log10') + 
  theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(name = "Model", values = c("darkorange4", "lightgoldenrod3")) + 
  theme(strip.background = element_blank()) + scale_y_continuous(breaks = c(0,1,2))



top2 <-  mmtoplcdist[mmtoplcdist$Site == "Willow Creek > 2010",]
top2 <- top2[year(top2$time) == 2013,]
top2$plot <- F
top2$plot[(top2$value > 0.1) & (top2$value < 10^5)] <- T
top2$bimon <- "Jan-Feb"
top2$bimon[month(top2$time)>=3 & month(top2$time)<=4] <- "Mar-Apr"
top2$bimon[month(top2$time)>=5 & month(top2$time)<=6] <- "May-Jun"
top2$bimon[month(top2$time)>=7 & month(top2$time)<=8] <- "Jul-Aug"
top2$bimon[month(top2$time)>=9 & month(top2$time)<=10] <- "Sep-Oct"
top2$bimon[month(top2$time)>=11 & month(top2$time)<=12] <- "Nov-Dec"
top2$bimon <- factor(top2$bimon, levels = c("Jan-Feb", "Mar-Apr", "May-Jun", "Jul-Aug", "Sep-Oct", "Nov-Dec"))
top2$plot[(top2$variable == "Cowan" | top2$variable == "Thom") & top2$plot] <- F
ggplot(top2[complete.cases(top2) & top2$plot,], aes(value)) + facet_grid(bimon~.) + ggtitle("Willow Creek 2013") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  geom_density(aes(fill = variable), alpha = 0.5) +
  theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(values = c("darkorange4", "lightgoldenrod3"))

top2 <-  mmtoplcdist[mmtoplcdist$Site == "Willow Creek > 2010",]
top2 <- top2[year(top2$time) == 2014,]
top2$plot <- F
top2$plot[(top2$value > 0.1) & (top2$value < 10^5)] <- T
top2$bimon <- "Jan-Feb"
top2$bimon[month(top2$time)>=3 & month(top2$time)<=4] <- "Mar-Apr"
top2$bimon[month(top2$time)>=5 & month(top2$time)<=6] <- "May-Jun"
top2$bimon[month(top2$time)>=7 & month(top2$time)<=8] <- "Jul-Aug"
top2$bimon[month(top2$time)>=9 & month(top2$time)<=10] <- "Sep-Oct"
top2$bimon[month(top2$time)>=11 & month(top2$time)<=12] <- "Nov-Dec"
top2$bimon <- factor(top2$bimon, levels = c("Jan-Feb", "Mar-Apr", "May-Jun", "Jul-Aug", "Sep-Oct", "Nov-Dec"))
top2$plot[(top2$variable == "Cowan" | top2$variable == "Thom") & top2$plot] <- F
ggplot(top2[complete.cases(top2) & top2$plot,], aes(value)) + facet_grid(bimon~.) + ggtitle("Willow Creek 2014") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  geom_density(aes(fill = variable), alpha = 0.5) + 
  theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(values = c("darkorange4", "lightgoldenrod3"))

top2 <-  mmtoplcdist[mmtoplcdist$Site == "Willow Creek > 2010",]
top2 <- top2[year(top2$time) == 2015,]
top2$plot <- F
top2$plot[(top2$value > 0.1) & (top2$value < 10^5)] <- T
top2$bimon <- "Jan-Feb"
top2$bimon[month(top2$time)>=3 & month(top2$time)<=4] <- "Mar-Apr"
top2$bimon[month(top2$time)>=5 & month(top2$time)<=6] <- "May-Jun"
top2$bimon[month(top2$time)>=7 & month(top2$time)<=8] <- "Jul-Aug"
top2$bimon[month(top2$time)>=9 & month(top2$time)<=10] <- "Sep-Oct"
top2$bimon[month(top2$time)>=11 & month(top2$time)<=12] <- "Nov-Dec"
top2$bimon <- factor(top2$bimon, levels = c("Jan-Feb", "Mar-Apr", "May-Jun", "Jul-Aug", "Sep-Oct", "Nov-Dec"))
top2$plot[(top2$variable == "Cowan" | top2$variable == "Thom") & top2$plot] <- F
ggplot(top2[complete.cases(top2) & top2$plot,], aes(value)) + facet_grid(bimon~.) + ggtitle("Willow Creek 2015") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  geom_density(aes(fill = variable), alpha = 0.5) + 
  theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(values = c("darkorange4", "lightgoldenrod3"))

top2 <-  mmtoplcdist[mmtoplcdist$Site == "Willow Creek > 2010",]
top2 <- top2[year(top2$time) == 2016,]
top2$plot <- F
top2$plot[(top2$value > 0.1) & (top2$value < 10^8)] <- T
top2$bimon <- "Jan-Feb"
top2$bimon[month(top2$time)>=3 & month(top2$time)<=4] <- "Mar-Apr"
top2$bimon[month(top2$time)>=5 & month(top2$time)<=6] <- "May-Jun"
top2$bimon[month(top2$time)>=7 & month(top2$time)<=8] <- "Jul-Aug"
top2$bimon[month(top2$time)>=9 & month(top2$time)<=10] <- "Sep-Oct"
top2$bimon[month(top2$time)>=11 & month(top2$time)<=12] <- "Nov-Dec"
top2$bimon <- factor(top2$bimon, levels = c("Jan-Feb", "Mar-Apr", "May-Jun", "Jul-Aug", "Sep-Oct", "Nov-Dec"))
top2$plot[(top2$variable == "Cowan" | top2$variable == "Thom") & top2$plot] <- F
ggplot(top2[complete.cases(top2) & top2$plot,], aes(value)) + facet_grid(bimon~.) + ggtitle("Willow Creek 2016") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  geom_density(aes(fill = variable), alpha = 0.5) + 
  theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(values = c("darkorange4", "lightgoldenrod3"))
# ggplot(top[complete.cases(top),], aes(value, beta)) + facet_grid(Site~variable, scales = "free") + 
#   geom_point(aes(col = Season, fill = Season), alpha = 0.2, size=0.2) + scale_x_continuous(trans = 'log10') + 
#   theme_classic() + xlab(expression(c[d]~aH[c])) + ylab(expression(beta)) + geom_hline(yintercept = 0.33)

trmdn <- tapply(mmtoplcdistd0$value[rep(calval == "Cal",4)], list(mmtoplcdistd0$Site[rep(calval == "Cal",4)], 
                                                                  mmtoplcdistd0$variable[rep(calval == "Cal",4)],
                                                                  mmtoplcdistd0$Season[rep(calval == "Cal",4)]), 
                FUN = function(x)(median(x,na.rm = T)))
tag <- as.data.frame(rbind(trmdn[,,1],trmdn[,,2], trmdn[,,3]))
tag$Site <- row.names(tag)[1:8]
tag$Site <- gsub("...2007", " < 2007", tag$Site)
tag$Site <- gsub("...2010", " > 2010", tag$Site)
tag$Site <- gsub("\\.", " ", tag$Site)
tag$Season <- rep(c("Evergreen","Leaf off", "Leaf on"), each= 8)
mtag <- melt(tag, id = c("Site", "Season"))
mtag <- mtag %>%
  mutate(cnh = case_when(
    Site == "Gingin" ~ 7,
    Site == "Alice Springs" ~ 6.5,
    Site == "Yanco" ~ 0.3,
    Site == "Charleston" ~ 7,#https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek < 2007" ~ 24, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Willow Creek > 2010" ~ 24, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 23, #https://www.neonscience.org/field-sites/bart
    Site == "Poker Flat" ~ 3 #https://www.sciencedirect.com/science/article/pii/S187396521300008X#fig2
  ))
##get bartlett info
mtag <- mtag %>%
  mutate(LAI = case_when(
    Site == "Gingin" ~ 0.7,
    Site == "Alice Springs" ~ 0.3,
    Site == "Yanco" ~ 0,
    Site == "Charleston" ~ 1.3, #https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek < 2007" ~ 5.3, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Willow Creek > 2010" ~ 5.3, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 3,
    Site == "Poker Flat" ~ 0.73 #https://ieeexplore.ieee.org/document/6589947
  ))
mtag <- mtag[complete.cases(mtag),]
mtag$cd <- mtag$value/mtag$LAI/mtag$cnh

ggplot(mmtoplcdist[complete.cases(mmtoplcdist) &  (mmtoplcdist$value > 10^-2) & (mmtoplcdist$value < (5*10^3)),], aes(value)) + 
  facet_grid(Site~variable, scales = "free") + geom_hline(yintercept = 0.333) +
  geom_density(aes(col = Season, fill = Season), alpha = 0.2) + scale_x_continuous(trans = 'log10') + 
  theme_classic() + xlab(expression(c[d]~aH[c])) + geom_vline(data = mtag, aes(xintercept = value, col = Season)) 
  #geom_label(data = mcmpr[mcmpr$Estimate == "Median",], aes(x=10, y = 1, label = as.character(round(value,2)))) 

mmtoplcdistz0 <- mmtoplcdist[mmtoplcdist$advarbyhc == "z0",]
ggplot(mmtoplcdistd0[mmtoplcdistd0$value >10^-4 & mmtoplcdistd0$value < 10^4 & complete.cases(mmtoplcdistd0),], aes(value, beta)) +
  facet_grid(variable~.) + geom_point(size = 0.6, shape = 1, aes(col = valadvar)) + theme_classic() + scale_x_continuous(trans = "log10") +
  ylim(0,2) + scale_color_viridis() + xlab(expression(c[d]~aH[c])) + ylab(expression(beta)) + geom_hline(yintercept = 1/3)
ggplot(mmtoplcdistz0[mmtoplcdistz0$value >10^-4 & mmtoplcdistz0$value < 10^4 & complete.cases(mmtoplcdistz0) & mmtoplcdistz0$valadvar < 0.7,], aes(value, beta)) +
  facet_grid(variable~.) + geom_point(size = 0.6, shape = 1, aes(col = valadvar)) + theme_classic() + scale_x_continuous(trans = "log10") +
  ylim(0,1) + scale_color_viridis() + xlab(expression(c[d]~aH[c])) + ylab(expression(beta)) + geom_hline(yintercept = 1/3)


cw68 <- function(x,y,par){
  (2*x/(1 - par[1])/0.4)-(sqrt(y/x/0.4/(1 - par[1]))/tanh(sqrt(y/x/0.4/(1 - par[1]))))
}
##Supplementary material figure
ggplot(mmtoplcdist[mmtoplcdist$advarbyhc == "d0",], aes(value)) + facet_grid(Site~variable, scales = "free") + 
     geom_density(aes(fill = Season), alpha = 0.5) + scale_x_continuous(trans = 'log10') + 
     theme_classic() + xlab(expression(c[d]~aH[c]))  + scale_fill_manual(values = c("green", "burlywood4", "forestgreen")) + 
     geom_vline(data = mtagtop, aes(xintercept = value, col = Season)) + scale_color_manual(values = c("green", "burlywood4", "forestgreen"))

##Figure 9 + Supplementary
for(i in 1:length(unique(mmtoplcdist$Site))){
    if(unique(mmtoplcdist$Site)[i] != "Yanco"){
    sub <- mmtoplcdist[complete.cases(mmtoplcdist) & (mmtoplcdist$Site == unique(mmtoplcdist$Site)[i]) & 
                           mmtoplcdist$calval == "Cal" & mmtoplcdist$advarbyhc == "d0",]
    toplc <- data.frame()
    for(j in c(.1,.5,.9)){
        toplc <- rbind(toplc, as.data.frame(tapply(sub$value, list(sub$variable, sub$Season) , FUN = function(x)quantile(x, probs = j))))
    }
    toplc$Model <- row.names(toplc)[1:4]
    mtoplc <- melt(toplc, id = "Model")
    llbt <- quantile(sub$beta, probs = .05, na.rm = T)
    ulbt <- quantile(sub$beta, probs = .95, na.rm = T)
    betacontr <- seq(llbt, ulbt, length.out = 100)
    cntrfin <- data.frame(beta = rep(betacontr, nrow(mtoplc)), Model = rep(mtoplc$Model,each = length(betacontr)), 
                       Season = rep(mtoplc$variable,each = length(betacontr)),
                       cdah = rep(mtoplc$value,each = length(betacontr)), 
                       cdahcat = rep(c("low", "med", "high"), each = 400), d0 = NA)
    cntrfin$d0[cntrfin$Model == "Inoue"] <- with(cntrfin[cntrfin$Model == "Inoue",], (1 - (2*beta^3/cdah/0.4)))
    for(k in 1:length(cntrfin$d0[cntrfin$Model == "Cowan"])){
      result <- optim(par = 0.999999, fn = cw68, x = cntrfin$beta[cntrfin$Model == "Cowan"][k], y =cntrfin$cdah[cntrfin$Model == "Cowan"][k],
                      method = "Brent", lower = 0.0000000000001, upper = 0.9999999999999)
      cntrfin$d0[cntrfin$Model == "Cowan"][k] <- result$par
    }
    cntrfin$d0[cntrfin$Model == "Thom"] <- with(cntrfin[cntrfin$Model == "Thom",], (1 - 6*beta^3/cdah/(4*0.4)))
    cntrfin$d0[cntrfin$Model == "Wang"] <- with(cntrfin[cntrfin$Model == "Wang",], (1 - (beta*(besselI(2*sqrt(cdah/0.4/beta),0) + 
                                                                                                 (besselI(2*sqrt(0.00001/beta/0.4*cdah),0)*
                                                                                                    besselK(2*sqrt(cdah/beta/0.4),0)/
                                                                                                    besselK(2*sqrt(0.00001*cdah/beta/0.4),0)))/
                                                                                           (0.4*sqrt(cdah/beta/0.4)*
                                                                                              (besselI(2*sqrt(cdah/beta/0.4),1) - 
                                                                                                 (besselI(2*sqrt(0.00001/beta/0.4*cdah),0)*
                                                                                                    besselK(2*sqrt(cdah/beta/0.4),1)/
                                                                                                    besselK(2*sqrt(0.00001/beta/0.4*cdah),0)))))))
    names(sub)[8] <- "Model"
    names(sub)[7] <- "d0"
    sub <- sub[sub$beta != Inf,]
    # df_ell <- data.frame()
    # for(j in 1:length(unique(sub$Season))){
    #   sub2 <- sub[sub$Season == unique(sub$Season)[j],]
    #   for(k in 1:length(unique(sub$Model))){
    #     df_ell <- rbind(df_ell, cbind(as.data.frame(with(sub2[sub2$Model==unique(sub2$Model)[k],], ellipse(cor(beta, d0), 
    #                                                                                                   scale=c(sd(beta),sd(d0)), 
    #                                                                                                   centre=c(mean(beta),mean(d0))))),
    #                                   Season=unique(sub$Season)[j], Model = unique(sub2$Model)[k]))
    #   }
    # }
    # 
    ribdat <- as.data.frame(cntrfin[cntrfin$cdahcat == "low",c(1,2,3,6)])
    ribdat$high <- cntrfin[cntrfin$cdahcat == "high",6]
    g2 <- ggplot(cntrfin[cntrfin$cdahcat == "med",], aes(beta, d0)) + 
      geom_ribbon(data= ribdat,aes(x = beta, ymin = d0, ymax = high, fill = Season), alpha=0.1)+ 
      geom_area(data = ribdat[ribdat$d0 <0,],aes(x=beta,y= high, fill = Season), alpha=0.1) +
      geom_point(data = sub, size = 0.6, shape = 1, aes(x = beta, y = d0, col = Season), alpha = 0.5) + 
      geom_line(aes(col = Season),size = 1) + facet_grid(.~Model, scales = "free") + 
      ylim (0,1) + theme_classic() + xlim(0,ulbt) +  ggtitle(sub$Site[1])  + xlab(expression(beta)) + 
      ylab(expression(frac(d[0],H[c]))) +  theme(strip.background = element_blank()) 
    if(length(unique(sub$Season))==1){
        g2 <- g2 + scale_color_manual(values = "forestgreen") + scale_fill_manual(values = "palegreen3")
    }else{
        g2 <- g2 + scale_color_manual(values = c("burlywood4", "darkgreen")) + scale_fill_manual(values = c("tan3", "seagreen3"))
    }
      # geom_path(data = df_ell, aes(x=x, y= y, col = Season))
    print(g2)
    }
  if(unique(mmtoplcdist$Site)[i] == "Gingin" | unique(mmtoplcdist$Site)[i] == "Charleston"){
    g2 <- ggplot(cntrfin[cntrfin$cdahcat == "med" & (cntrfin$Model == "Inoue" | cntrfin$Model == "Wang"),], aes(beta, d0)) + 
      geom_ribbon(data= ribdat[(ribdat$Model == "Inoue" | ribdat$Model == "Wang"),],aes(x = beta, ymin = d0, ymax = high, fill = Season), alpha=0.1)+ 
      geom_area(data = ribdat[ribdat$d0 <0 & (ribdat$Model == "Inoue" | ribdat$Model == "Wang"),],aes(x=beta,y= high, fill = Season), alpha=0.1) +
      geom_point(data = sub[(sub$Model == "Inoue" | sub$Model == "Wang"),], size = 0.6, shape = 1, aes(x = beta, y = d0, col = Season), alpha = 0.5) + 
      geom_line(aes(col = Season),size = 1) + facet_grid(.~Model, scales = "free") + 
      ylim (0,1) + theme_classic() + xlim(0,0.9) +  ggtitle(sub$Site[1])  + xlab(expression(beta)) + 
      ylab(expression(frac(d[0],H[c]))) +   theme(strip.background = element_blank())
    if(length(unique(sub$Season))==1){
      g2 <- g2 + scale_color_manual(values = "forestgreen") + scale_fill_manual(values = "palegreen3")
    }else{
      g2 <- g2 + scale_color_manual(values = c("burlywood4", "darkgreen")) + scale_fill_manual(values = c("tan3", "seagreen3"))
    }
    # geom_path(data = df_ell, aes(x=x, y= y, col = Season))
    print(g2)
  }
}

 ##new set
toplcdist <- fincomb[,c(6,7,2,5,1)]
mtoplcdist <- melt(toplcdist, id = c("Site", "beta", "time"))
names(mtoplcdist)[4:5] <- c("advarbyhc","valadvar")
mtoplcdist$Wglc <- c(Wgd0lc, Wgz0lc)
sub <- mtoplcdist[complete.cases(mtoplcdist) & (mtoplcdist$Site == unique(mtoplcdist$Site)[4]),]
sub <- sub[sub$advarbyhc == "d0",]
xmx <- quantile(sub$Wglc, probs = .85, na.rm = T)
g1 <- ggplot(sub, aes(Wglc)) + geom_density(aes(fill = as.factor(quarter(time))), alpha = 0.2) + theme_bw() + 
      ggtitle(unique(mtoplcdist$Site)[4])  + xlim(0,xmx)
g1 +  scale_fill_viridis(discrete = TRUE)


sub <- mtoplcdist[complete.cases(mtoplcdist) & (mtoplcdist$Site == unique(mtoplcdist$Site)[5]),]
sub <- sub[sub$advarbyhc == "d0",]
xmx <- quantile(sub$Wglc, probs = .85, na.rm = T)
g1 <- ggplot(sub, aes(Wglc)) +  geom_density(aes(fill = as.factor(quarter(time))), alpha = 0.2) + theme_bw() + 
      ggtitle(unique(mtoplcdist$Site)[5])  + xlim(0,xmx)
g1 +  scale_fill_viridis(discrete = TRUE)

sub <- mtoplcdist[complete.cases(mtoplcdist) & (mtoplcdist$Site == unique(mtoplcdist$Site)[6]),]
sub <- sub[sub$advarbyhc == "d0",]
xmx <- quantile(sub$Wglc, probs = .95, na.rm = T)
g1 <- ggplot(sub, aes(Wglc)) +  geom_density(aes(fill = as.factor(quarter(time))), alpha = 0.2) + theme_bw() + 
        ggtitle(unique(mtoplcdist$Site)[6])  + xlim(0,xmx)
g1 +  scale_fill_viridis(discrete = TRUE)

# sub <- mtoplcdist[complete.cases(mtoplcdist) & (mtoplcdist$Site == unique(mtoplcdist$Site)[7]),]
# sub <- sub[sub$advarbyhc == "d0",]
# xmx <- quantile(sub$Wglc, probs = .85, na.rm = T)
# g1 <- ggplot(sub, aes(Wglc)) +  geom_density(aes(fill = as.factor(quarter(time))), alpha = 0.2) + theme_bw() + ggtitle(unique(mtoplcdist$Site)[7])+ xlim(0,xmx)
# g1+  scale_fill_viridis(discrete = TRUE)
sub <- mtoplcdist[complete.cases(mtoplcdist) & (mtoplcdist$Site == unique(mtoplcdist$Site)[1] | mtoplcdist$Site == unique(mtoplcdist$Site)[2] | mtoplcdist$Site == unique(mtoplcdist$Site)[7]),]
sub <- sub[sub$advarbyhc == "d0",]
xmx <- quantile(sub$Wglc, probs = .9, na.rm = T)
g1 <- ggplot(sub, aes(Wglc)) + facet_wrap(.~Site, scales = "free") + geom_density(aes(fill = as.factor(quarter(time))), alpha = 0.2) + theme_bw() + xlim(0,xmx) + ggtitle("Evergreen Sites")
print(g1+  scale_fill_viridis(discrete = TRUE))

sub <- mtoplcdist[complete.cases(mtoplcdist) & (mtoplcdist$Site == unique(mtoplcdist$Site)[1] | mtoplcdist$Site == unique(mtoplcdist$Site)[4] | mtoplcdist$Site == unique(mtoplcdist$Site)[6] | mtoplcdist$Site == unique(mtoplcdist$Site)[5]),]
sub <- sub[sub$advarbyhc == "d0",]
sub$Season <- "Leaf off"
sub$Season[month(sub$time) > 5 & month(sub$time) < 11] <- "Leaf on"
sub$Season[sub$Site == "Gingin"] <- ""
sub$cat <- paste(sub$Site, sub$Season)
toplc <- tapply(sub$Wglc, sub$cat, FUN = function(x)quantile(x, probs = c(.1,.5,.9)))
betacontr <- seq(0.01, 1, length.out = 100)
cntr <- data.frame(beta = rep(betacontr, 21), cat = rep(names(toplc),each = 300), cdah = NA, 
                   cdahcat = rep(rep(c("low", "med", "high"), each = 100), 7), d0 = NA)
for(i in 1:7){
  for(j in 1:3){
    cntr[(((i-1)*300)+((j-1)*100)+1):(((i-1)*300)+((j-1)*100)+100),3] <- toplc[[i]][j]
    cntr[(((i-1)*300)+((j-1)*100)+1):(((i-1)*300)+((j-1)*100)+100),5] <- (1 - (betacontr*(besselI(2*sqrt(toplc[[i]][j]/0.4/betacontr),0) + 
                                                                      (besselI(2*sqrt(0.000001/betacontr/0.4*toplc[[i]][j]),0)*
                                                                         besselK(2*sqrt(toplc[[i]][j]/betacontr/0.4),0)/
                                                                         besselK(2*sqrt(0.000001*toplc[[i]][j]/betacontr/0.4),0)))/
                                                           (0.4*sqrt(toplc[[i]][j]/betacontr/0.4)*(besselI(2*sqrt(toplc[[i]][j]/betacontr/0.4),1) - 
                                                            (besselI(2*sqrt(0.000001/betacontr/0.4*toplc[[i]][j]),0)*
                                                            besselK(2*sqrt(toplc[[i]][j]/betacontr/0.4),1)/
                                                            besselK(2*sqrt(0.000001/betacontr/0.4*toplc[[i]][j]),0))))))
  }
}

g2 <- ggplot(cntr, aes(beta, d0)) + geom_point(data = sub, aes(x = beta, y = valadvar), col = "grey", alpha = 0.1) + 
      geom_line(aes(col = cdahcat)) + facet_wrap(.~cat) + ylim (0,1) + xlim(0,1)+ 
      theme_classic() 
print(g2)

ggplot(mmtoplcdist[complete.cases(mmtoplcdist),], aes(value)) + facet_grid(Site~variable, scales = "free") + 
  geom_density(aes(col = advarbyhc, fill = advarbyhc), alpha = 0.2) + scale_x_continuous(trans = 'log10') + theme_classic() + xlab(expression(c[d]~aH[c]))

##last graph
count <- 1
fintop <- data.frame()
fintop2 <- data.frame()
for(i in 1:length(unique(fincomb$Site))){
  sub <- fincomb[fincomb$Site == unique(fincomb$Site)[i] & calval == "Val" & is.na(fincomb$limit),-14]
  if(count == 4){
    count <- 5
  }
  for(j in 1:length(unique(sub$Season))){
    sub2 <- sub[(sub$Season == unique(sub$Season)[j]) & complete.cases(sub),]
    if(nrow(sub2) != 0){
      insb <- mtag[mtag$Site == unique(fincomb2$Site)[i] & mtag$Season == unique(sub$Season)[j] & mtag$variable == "Inoue",4]
      cwsb <- mtag[mtag$Site == unique(fincomb2$Site)[i] & mtag$Season == unique(sub$Season)[j] & mtag$variable == "Cowan",4]
      cwans <- 0
      for(k in 1:length(sub2$beta)){
        result <- optim(par = 0.999999, fn = cw68, x = sub2$beta[k], y =cwsb,
                        method = "Brent", lower = 0.0000000000001, upper = 0.9999999999999)
        cwans[k] <- result$par*sub2$cnh[k]
      }
      thsb <- mtag[mtag$Site == unique(fincomb2$Site)[i] & mtag$Season == unique(sub$Season)[j] & mtag$variable == "Thom",4]
      wgsb <- mtag[mtag$Site == unique(fincomb2$Site)[i] & mtag$Season == unique(sub$Season)[j] & mtag$variable == "Wang",4]
      ans <- data.frame(Site = sub$Site[1], Season = sub2$Season[1],
                        Emperical = sub2$d0*sub2$cnh, U_LM = udlst2[[count]]$coefficients[1]+udlst2[[count]]$coefficients[2]*sub2$U,
                        ustar_LM = ustrdlst2[[count]]$coefficients[1]+ustrdlst2[[count]]$coefficients[2]*sub2$ustar,
                        Inoue = (1 - (2*sub2$beta^3/insb/0.4))*sub2$cnh,
                        Cowan = cwans,
                        Thom = (1 - (6/4*sub2$beta^3/thsb/0.4))*sub2$cnh,
                        Wang = (1 - (sub2$beta*(besselI(2*sqrt(wgsb/0.4/sub2$beta),0) + 
                                                 (besselI(2*sqrt(0.000001/sub2$beta/0.4*wgsb),0)*
                                                    besselK(2*sqrt(wgsb/sub2$beta/0.4),0)/
                                                    besselK(2*sqrt(0.000001*wgsb/sub2$beta/0.4),0)))/
                                       (0.4*sqrt(wgsb/sub2$beta/0.4)*(besselI(2*sqrt(wgsb/sub2$beta/0.4),1) - 
                                                                          (besselI(2*sqrt(0.000001/sub2$beta/0.4*wgsb),0)*
                                                                             besselK(2*sqrt(wgsb/sub2$beta/0.4),1)/
                                                                             besselK(2*sqrt(0.000001/sub2$beta/0.4*wgsb),0))))))*sub2$cnh)
      fintop <- rbind(fintop, ans)
      ans2 <- data.frame(Site = sub$Site[1], Season = sub2$Season[1],
                        Emperical = sub2$z0*sub2$cnh, U_LM = uzlst2[[count]]$coefficients[1]+uzlst2[[count]]$coefficients[2]*sub2$U,
                        ustar_LM = ustrzlst2[[count]]$coefficients[1]+ustrzlst2[[count]]$coefficients[2]*sub2$ustar,
                        Inoue = (2*sub2$beta^3/insb/0.4)/exp(0.4/sub2$beta)*sub2$cnh,
                        Cowan = (1-cwans/sub2$cnh)/exp(0.4/sub2$beta)*sub2$cnh,
                        Thom = (6/4*sub2$beta^3/thsb/0.4)/exp(0.4/sub2$beta)*sub2$cnh,
                        Wang = ((sub2$beta*(besselI(2*sqrt(wgsb/0.4/sub2$beta),0) + 
                                                  (besselI(2*sqrt(0.000001/sub2$beta/0.4*wgsb),0)*
                                                     besselK(2*sqrt(wgsb/sub2$beta/0.4),0)/
                                                     besselK(2*sqrt(0.000001*wgsb/sub2$beta/0.4),0)))/
                                       (0.4*sqrt(wgsb/sub2$beta/0.4)*(besselI(2*sqrt(wgsb/sub2$beta/0.4),1) - 
                                                                        (besselI(2*sqrt(0.000001/sub2$beta/0.4*wgsb),0)*
                                                                           besselK(2*sqrt(wgsb/sub2$beta/0.4),1)/
                                                                           besselK(2*sqrt(0.000001/sub2$beta/0.4*wgsb),0))))))/
                          exp(0.4/sub2$beta)*sub2$cnh)
      fintop2 <- rbind(fintop2, ans2)
      count <- count+1
    }else{
      count <- count+1
    }
  }
}

mans <- melt(fintop, id = c("Site","Season","Emperical"))
mans <- mans[complete.cases(mans),]
mans <- mans %>%
  mutate(cnh = case_when(
    Site == "Gingin" ~ 7,
    Site == "Alice Springs" ~ 6.5,
    Site == "Yanco" ~ 0.3,
    Site == "Charleston" ~ 7,#https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek < 2007" ~ 24,
    Site == "Willow Creek > 2010" ~ 24,#https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 23, #https://www.neonscience.org/field-sites/bart
    Site == "Poker Flat" ~ 3 #https://www.sciencedirect.com/science/article/pii/S187396521300008X#fig2
  ))
mans2 <- mans[mans$value > 0 & mans$value < mans$cnh,]

whywesuck <- data.frame()
for(i in 1:length(unique(mans2$Site))){
  sub <- mans2[mans2$Site == unique(mans2$Site)[i],]
  for(j in 1:length(unique(sub$Season))){
    sub2 <- sub[sub$Season == unique(sub$Season)[j],]
    for(k in 1:length(unique(sub2$variable))){
      sub3 <- sub2[sub2$variable == unique(sub2$variable)[k],]
      mod <- lm(value~Emperical, data= sub3)
      ans <- data.frame(Site = unique(mans2$Site)[i],
                        Season = unique(sub$Season)[j],
                        variable = unique(sub2$variable)[k],
                        r2 = summary(mod)$r.square,
                        adjr2 = summary(mod)$adj.r.square,
                        rmse = summary(mod)$sigma,
                        int = summary(mod)$coefficient[1],
                        slp = summary(mod)$coefficient[2])
      whywesuck <- rbind(whywesuck,ans)
    }
  }
}
write.csv(whywesuck, "whywesuck.csv")

mans <- melt(fintop2, id = c("Site","Season","Emperical"))
mans <- mans[complete.cases(mans),]
mans <- mans %>%
  mutate(cnh = case_when(
    Site == "Gingin" ~ 7,
    Site == "Alice Springs" ~ 6.5,
    Site == "Yanco" ~ 0.3,
    Site == "Charleston" ~ 7,#https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek < 2007" ~ 24,
    Site == "Willow Creek > 2010" ~ 24,#https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 23, #https://www.neonscience.org/field-sites/bart
    Site == "Poker Flat" ~ 3 #https://www.sciencedirect.com/science/article/pii/S187396521300008X#fig2
  ))
mans2 <- mans[mans$value > 0 & mans$value < mans$cnh,]

whywesuck <- data.frame()
for(i in 1:length(unique(mans2$Site))){
  sub <- mans2[mans2$Site == unique(mans2$Site)[i],]
  for(j in 1:length(unique(sub$Season))){
    sub2 <- sub[sub$Season == unique(sub$Season)[j],]
    for(k in 1:length(unique(sub2$variable))){
      sub3 <- sub2[sub2$variable == unique(sub2$variable)[k],]
      mod <- lm(value~Emperical, data= sub3)
      ans <- data.frame(Site = unique(mans2$Site)[i],
                        Season = unique(sub$Season)[j],
                        variable = unique(sub2$variable)[k],
                        r2 = summary(mod)$r.square,
                        adjr2 = summary(mod)$adj.r.square,
                        rmse = summary(mod)$sigma,
                        int = summary(mod)$coefficient[1],
                        slp = summary(mod)$coefficient[2])
      whywesuck <- rbind(whywesuck,ans)
    }
  }
}
write.csv(whywesuck, "whywesuckless.csv")

ggplot(mans2, aes(value)) + geom_density(aes(x = Emperical, col = Season, fill = Season), alpha = 0.2) + theme_classic() + 
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") + xlab(expression(d[0]~~"[m]")) +
  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + scale_fill_manual(values = c("green", "burlywood4", "forestgreen"))
ggplot(mans2[mans2$Site == "Gingin",], aes(value)) + geom_density(aes(x = Emperical, col = Season, fill = Season), alpha = 0.2) + theme_classic() +
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") + 
  scale_color_manual(values = c("green")) + scale_fill_manual(values = c( "green")) + xlab(expression(d[0]~~"[m]")) 
ggplot(mans2[mans2$Site == "Charleston",], aes(value)) + geom_density(aes(x = Emperical, col = Season, fill = Season), alpha = 0.2) + theme_classic() +
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") +
  scale_color_manual(values = c("burlywood4", "forestgreen")) + scale_fill_manual(values = c( "burlywood4", "forestgreen")) + xlab(expression(d[0]~~"[m]"))

print(g1)
g2 <- ggplot(mans, aes(value)) + geom_density(aes(x = NLS), col = "grey", fill = "grey", alpha = 0.2) + theme_classic() + facet_wrap(.~variable, scales = "free") + geom_density(aes(col = variable))+ coord_flip() + ggtitle(unique(fincomb$Site)[i])
print(g2)
st[count] <- unique(fincomb2$Site)[i]
sn[count] <- unique(sub$Season)[j]
ggplot(top[complete.cases(top),], aes(U, wlc)) + geom_point(shape = 1, alpha = 0.3, col = "grey") + facet_wrap(.~Site) + ylim(0,3) + theme_classic() + geom_smooth()

mans2 <- melt(fintop2, id = c("Site","Season","Emperical"))
mans2 <- mans2[complete.cases(mans2),]
mans2 <- mans2 %>%
  mutate(cnh = case_when(
    Site == "Gingin" ~ 7,
    Site == "Alice Springs" ~ 6.5,
    Site == "Yanco" ~ 0.3,
    Site == "Charleston" ~ 7,#https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek 1" ~ 24,
    Site == "Willow Creek 2" ~ 24,#https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 23, #https://www.neonscience.org/field-sites/bart
    Site == "Poker Flat" ~ 3 #https://www.sciencedirect.com/science/article/pii/S187396521300008X#fig2
  ))
mans2 <- mans2[mans2$value >0 & mans2$value < 7,]
ggplot(mans2, aes(value)) + geom_density(aes(x = Emperical, col = Season, fill = Season), alpha = 0.2) + theme_classic() + 
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") + 
  xlab(expression(z[0]~~"[m]")) +  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + 
  scale_fill_manual(values = c("green", "burlywood4", "forestgreen"))
ggplot(mans2[mans2$Site == "Gingin" & mans2$variable != "Cowan" & mans2$variable != "Thom",], aes(value)) + geom_density(aes(x = Emperical, fill = Season), col = "white",alpha = 0.2) + theme_classic() +
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") +
  scale_color_manual(values = c("forestgreen")) + scale_fill_manual(values = c( "forestgreen")) + xlab(expression(d[0]~~"[m]")) +   theme(strip.background = element_blank())

ggplot(mans[mans$Site == "Bartlett",], aes(value)) + geom_density(aes(x = Emperical, col = Season, fill = Season), alpha = 0.2) + theme_classic() +
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") +  theme(strip.background = element_blank()) +
  scale_color_manual(values = c("burlywood4", "forestgreen")) + scale_fill_manual(values = c( "burlywood4", "forestgreen")) + xlab(expression(d[0]~~"[m]"))

ggplot(mans2[mans2$Site == "Charleston" & mans2$variable != "Cowan" & mans2$variable != "Thom",], aes(value)) + geom_density(aes(x = Emperical, fill = Season), col = "white", alpha = 0.2) + theme_classic() +
  facet_wrap(.~variable) + geom_density(aes(col = Season))+ coord_flip() + facet_grid(Site~variable, scales = "free") +  theme(strip.background = element_blank()) +
  scale_color_manual(values = c("burlywood4", "darkgreen")) + scale_fill_manual(values = c( "burlywood4", "darkgreen")) + xlab(expression(d[0]~~"[m]"))

print(g1)
g2 <- ggplot(mans, aes(value)) + geom_density(aes(x = NLS), col = "grey", fill = "grey", alpha = 0.2) + theme_classic() + facet_wrap(.~variable, scales = "free") + geom_density(aes(col = variable))+ coord_flip() + ggtitle(unique(fincomb$Site)[i])
print(g2)
st[count] <- unique(fincomb2$Site)[i]
sn[count] <- unique(sub$Season)[j]
ggplot(top[complete.cases(top),], aes(U, wlc)) + geom_point(shape = 1, alpha = 0.3, col = "grey") + facet_wrap(.~Site) + ylim(0,3) + theme_classic() + geom_smooth()
library(mltools)
ansd0r2 <- data.frame()
ansz0r2 <- data.frame()
ansd0rm <- data.frame()
ansz0rm <- data.frame()
for(i in 1:length(unique(fintop$Site))){
  subd0 <- fintop[fintop$Site == unique(fintop$Site)[i],]
  subz0 <- fintop2[fintop2$Site == unique(fintop2$Site)[i],]
  subd0r2 <- as.data.frame(cbind(data.frame(Site = unique(fintop2$Site)[i]), subd0[1,4:9]))
  subz0r2 <- as.data.frame(cbind(data.frame(Site = unique(fintop2$Site)[i]), subd0[1,4:9]))
  subd0rm <- as.data.frame(cbind(data.frame(Site = unique(fintop2$Site)[i]), subd0[1,4:9]))
  subz0rm <- as.data.frame(cbind(data.frame(Site = unique(fintop2$Site)[i]), subd0[1,4:9]))
  d0bins <- (0:20)/20*max(subd0$Emperical, na.rm = T)
  z0bins <- seq(from = as.numeric(quantile(subz0$Emperical, probs = .05)), to = as.numeric(quantile(subz0$Emperical, probs = .95)), length.out = 20)
  xd0 <- apply(subd0[,3:9], 2, FUN = function(x) table(bin_data(x, bins = d0bins)))
  xz0 <- apply(subz0[,3:9], 2, FUN = function(x) table(bin_data(x, bins = z0bins)))
  for(j in 2:7){
    mod <- lm(xd0[,j]~xd0[,1])
    subd0r2[1,j] <- summary(mod)$adj.r.squared
    subd0rm[1,j] <- summary(mod)$sigma
    mod <- lm(xz0[,j]~xz0[,1])
    subz0r2[1,j] <- summary(mod)$adj.r.squared
    subz0rm[1,j] <- summary(mod)$sigma
  }
  ansd0r2 <- rbind(ansd0r2, subd0r2)
  ansz0r2 <- rbind(ansz0r2, subz0r2)
  ansd0rm <- rbind(ansd0rm, subd0rm)
  ansz0rm <- rbind(ansz0rm, subz0rm)
}
ansd0r2$var <- "d0"
ansz0r2$var <- "z0"
anstbl <- rbind(ansd0r2, ansz0r2)
write.csv(anstbl, "fintbl.csv")
##checking velocity
USCMWcan$Site <- "Charleston"
USCMWcan$Inoue <- Ind0lc[fincomb$Site=="Charleston"]
USCMWcan$Cowan <- Cwd0lc[fincomb$Site=="Charleston"]
USCMWcan$Thom <- Thd0lc[fincomb$Site=="Charleston"]
USCMWcan$Wang <- Wgd0lc[fincomb$Site=="Charleston"]
USCMWcan$beta <- fincomb$beta[fincomb$Site=="Charleston"]
USCMWcan$d0h <- fincomb$d0[fincomb$Site=="Charleston"]
USWCrcan$Site <- "Willow Creek"
USWCrcan$Inoue <- Ind0lc[fincomb$Site=="Willow Creek"]
USWCrcan$Cowan <- Cwd0lc[fincomb$Site=="Willow Creek"]
USWCrcan$Thom <- Thd0lc[fincomb$Site=="Willow Creek"]
USWCrcan$Wang <- Wgd0lc[fincomb$Site=="Willow Creek"]
USWCrcan$beta <- fincomb$beta[fincomb$Site=="Willow Creek"]
USWCrcan$d0h <- fincomb$d0[fincomb$Site=="Willow Creek"]
USxBrcan$Site <- "Bartlett"
USxBrcan$Inoue <- Ind0lc[fincomb$Site=="Bartlett"]
USxBrcan$Cowan <- Cwd0lc[fincomb$Site=="Bartlett"]
USxBrcan$Thom <- Thd0lc[fincomb$Site=="Bartlett"]
USxBrcan$Wang <- Wgd0lc[fincomb$Site=="Bartlett"]
USxBrcan$beta <- fincomb$beta[fincomb$Site=="Bartlett"]
USxBrcan$d0h <- fincomb$d0[fincomb$Site=="Bartlett"]
USPrrcan <- USPrrlst[[2]]
USPrrcan <- USPrrcan[USPrrcan$H < 3,]
USPrrcan$Site <- "Poker Flat"
USPrrcan$Inoue <- Ind0lc[fincomb$Site=="Poker Flat"]
USPrrcan$Cowan <- Cwd0lc[fincomb$Site=="Poker Flat"]
USPrrcan$Thom <- Thd0lc[fincomb$Site=="Poker Flat"]
USPrrcan$Wang <- Wgd0lc[fincomb$Site=="Poker Flat"]
USPrrcan$beta <- fincomb$beta[fincomb$Site=="Poker Flat"]
USPrrcan$d0h <- fincomb$d0[fincomb$Site=="Poker Flat"]
OZGin <- wsp[wsp$H < 7,]
OZGin$Site <- "Gingin"
OZGin$Inoue <- Ind0lc[fincomb$Site=="Gingin"]
OZGin$Cowan <- Cwd0lc[fincomb$Site=="Gingin"]
OZGin$Thom <- Thd0lc[fincomb$Site=="Gingin"]
OZGin$Wang <- Wgd0lc[fincomb$Site=="Gingin"]
OZGin$beta <- fincomb$beta[fincomb$Site=="Gingin"]
OZGin$d0h <- fincomb$d0[fincomb$Site=="Gingin"]
OZAls <- wspas[wspas$H < 6.5,]
OZAls$Site <- "Alice Springs"
OZAls$Inoue <- Ind0lc[fincomb$Site=="Alice Springs"]
OZAls$Cowan <- Cwd0lc[fincomb$Site=="Alice Springs"]
OZAls$Thom <- Thd0lc[fincomb$Site=="Alice Springs"]
OZAls$Wang <- Wgd0lc[fincomb$Site=="Alice Springs"]
OZAls$beta <- fincomb$beta[fincomb$Site=="Alice Springs"]
OZAls$d0h <- fincomb$d0[fincomb$Site=="Alice Springs"]
setcan <- rbind(OZGin, OZAls, USPrrcan[,c(2,1,3:11)])
setcan$UH <- setcan$ustar/setcan$beta
setcan2 <- rbind(USCMWcan,USWCrcan,USxBrcan)
names(setcan2)[4:5] <- c("H", "WS")
canchck <- rbind(setcan, setcan2)
rm(OZGin, USPrrcan, OZAls, USCMWcan,USWCrcan,USxBrcan)
canchck$WS[canchck$WS==-9999] <- NA
canchck <- canchck[complete.cases(canchck),]
canchck$WSUH <- canchck$WS/canchck$UH
canchck <- canchck[canchck$WSUH < 1,]
canchck <- canchck %>%
  mutate(cnh = case_when(
    Site == "Gingin" ~ 7,
    Site == "Alice Springs" ~ 6.5,
    Site == "Yanco" ~ 0.3,
    Site == "Charleston" ~ 7,#https://pubag.nal.usda.gov/download/6822/PDF
    Site == "Willow Creek" ~ 24, #https://www.sciencedirect.com/science/article/pii/S0168192304001601
    Site == "Bartlett" ~ 23, #https://www.neonscience.org/field-sites/bart
    Site == "Poker Flat" ~ 3 #https://www.sciencedirect.com/science/article/pii/S187396521300008X#fig2
  ))
canchck$zeta <- canchck$H/canchck$cnh
canchck$limit <- NA
canchck$limit[canchck$d0t == (canchck$cnh-0.01)] <- "Upper"
canchck$limit[canchck$d0t == 0.005] <- "Lower"
canchck <- canchck[is.na(canchck$limit), -16]
canchck$Season <- "Leaf on"
canchck$Season[month(canchck$time) > 11 | month(canchck$time) < 6] <- "Leaf off"
canchck$Season[canchck$Site == "Gingin" | canchck$Site == "Alice Springs" | canchck$Site == "Poker Flat"] <- "Evergreen"
canchck$Site[canchck$Site == "Willow Creek" & year(canchck$time) < 2007] <- "Willow Creek 1"
canchck$Site[canchck$Site == "Willow Creek" & year(canchck$time) > 2007] <- "Willow Creek 2"

mcanchck <- melt(canchck[,c(15,13,10,11,3,5,6,7,8,9,16)], id= names(canchck)[c(15,13,10,11,3,5,16)])
mcanchck$valuevarlc[1:(nrow(canchck))] <- Inoue(canchck$zeta, canchck$beta, canchck$Inoue)
mcanchck$valuevarlc[(nrow(canchck)+1):(2*nrow(canchck))] <- Cowanorg(canchck$zeta, canchck$beta, canchck$Cowan, canchck$d0h)
mcanchck$valuevarlc[(2*nrow(canchck)+1):(3*nrow(canchck))] <- Thom(canchck$zeta, canchck$beta, canchck$Thom, canchck$d0h)
mcanchck$valuevarlc[(3*nrow(canchck)+1):(4*nrow(canchck))] <- Wang(canchck$zeta, canchck$beta, canchck$Wang, 0.000001)
fitness <- data.frame()
for(i in 1:length(unique(mcanchck$Site))){
  for(j in 1:length(unique(mcanchck$variable))){
    sub <- mcanchck[mcanchck$Site == unique(mcanchck$Site)[i] & mcanchck$variable == unique(mcanchck$variable)[j],]
    x <- lm(sub$valuevarlc ~ sub$value)
    fitness <- rbind(fitness, data.frame(Site = unique(mcanchck$Site)[i], model = unique(mcanchck$variable)[j], slp = x$coefficients[2], int = x$coefficients[1], rmse = summary(x)$sigma, adjr2 = summary(x)$adj.r.squared, n = nrow(sub)))
  }
}
tomrg <- mtag
names(tomrg)[4] <- "CnstLc"
mmcanchck <- merge(mcanchck, tomrg)
mmcanchck$valuecnstlc <- NA
mmcanchck$valuecnstlc[mmcanchck$variable == "Cowan"] <- Cowanorg(mmcanchck$zeta[mmcanchck$variable == "Cowan"], 
                                                     mmcanchck$beta[mmcanchck$variable == "Cowan"], 
                                                     mmcanchck$CnstLc[mmcanchck$variable == "Cowan"], 
                                                     mmcanchck$d0h[mmcanchck$variable == "Cowan"])
mmcanchck$valuecnstlc[mmcanchck$variable == "Inoue"] <- Inoue(mmcanchck$zeta[mmcanchck$variable == "Inoue"], 
                                                                    mmcanchck$beta[mmcanchck$variable == "Inoue"], 
                                                                    mmcanchck$CnstLc[mmcanchck$variable == "Inoue"])
mmcanchck$valuecnstlc[mmcanchck$variable == "Thom"] <- Thom(mmcanchck$zeta[mmcanchck$variable == "Thom"], 
                                                                     mmcanchck$beta[mmcanchck$variable == "Thom"], 
                                                                     mmcanchck$CnstLc[mmcanchck$variable == "Thom"],
                                                                     mmcanchck$d0h[mmcanchck$variable == "Thom"])
mmcanchck$valuecnstlc[mmcanchck$variable == "Wang"] <- Wang(mmcanchck$zeta[mmcanchck$variable == "Wang"], 
                                                            mmcanchck$beta[mmcanchck$variable == "Wang"], 
                                                            mmcanchck$CnstLc[mmcanchck$variable == "Wang"], 
                                                            0.000001)

fitness <- data.frame()
fitnesscnst <- data.frame()
for(i in 1:length(unique(mmcanchck$Site))){
  for(j in 1:length(unique(mmcanchck$variable))){
    sub <- mmcanchck[mmcanchck$Site == unique(mmcanchck$Site)[i] & mmcanchck$variable == unique(mmcanchck$variable)[j],]
    for(k in 1:length(unique(sub$Season))){
      x <- lm(sub$valuevarlc[sub$Season == unique(sub$Season)[k]] ~ sub$value[sub$Season == unique(sub$Season)[k]])
      fitness <- rbind(fitness, data.frame(Site = unique(mmcanchck$Site)[i], model = unique(mmcanchck$variable)[j], 
                                           Season = unique(sub$Season)[k], slp = x$coefficients[2], int = x$coefficients[1], 
                                           rmse = summary(x)$sigma, adjr2 = summary(x)$adj.r.squared, n = nrow(sub)))
      x <- lm(sub$valuecnstlc[sub$Season == unique(sub$Season)[k]] ~ sub$value[sub$Season == unique(sub$Season)[k]])
      fitnesscnst <- rbind(fitness, data.frame(Site = unique(mmcanchck$Site)[i], model = unique(mmcanchck$variable)[j], 
                                               Season = unique(sub$Season)[k], slp = x$coefficients[2], int = x$coefficients[1], 
                                               rmse = summary(x)$sigma, adjr2 = summary(x)$adj.r.squared, n = nrow(sub)))
    }
  }
}
write.csv(fitnesscnst, "worstperforming.csv")
ggplot(mmcanchck, aes(valuecnstlc)) + geom_density(aes(col = Season)) + facet_wrap(Site~variable, scales = "free") + theme_classic() + 
  geom_density(aes(x = WSUH), fill= "season", col = "black", alpha = 0.3)
ggplot(mmcanchck, aes(WSUH, valuecnstlc)) + geom_point(aes(col = Season), size = 0.4, shape = 1, alpha = 0.4) + 
  facet_grid(paste(Site,round(zeta,2))~variable) + theme_classic() + geom_abline(intercept = 0, slope = 1)  + 
  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + xlab(expression(U[obs])) + ylab(expression(U[mod]))
ggplot(mmcanchck, aes(WSUH, valuevarlc)) + geom_point(aes(col = Season), size = 0.4, shape = 1, alpha = 0.4) + 
  facet_grid(paste(Site,round(zeta,2))~variable) + theme_classic() + geom_abline(intercept = 0, slope = 1)  + 
  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + xlab(expression(U[obs])) + ylab(expression(U[mod]))

top <- mmcanchck[(mmcanchck$Site == "Gingin" | mmcanchck$Site == "Charleston") & (mmcanchck$variable == "Inoue" | mmcanchck$variable == "Wang"),]
#Figure 11
ggplot(top, aes(WSUH, valuecnstlc)) + geom_point(aes(col = Season), size = 0.4, shape = 1, alpha = 0.4) + 
  facet_grid(paste(Site,round(zeta,2))~variable) + theme_classic() + geom_abline(intercept = 0, slope = 1)  + 
  scale_color_manual(values = c("green", "burlywood4", "forestgreen")) + xlab(expression(frac(U[obs], U[H]))) + ylab(expression(frac(U[mod], U[H])))

##Supplementary material
numderd <- function(dudz, Et0dat, n, time, ntr){
  d0 <- 0
  comb <- dim(dudz)[1]/n
  tim <- rep(time, comb)
  for(i in 1:length(tim)){
    if(tim[i] %in% ntr){
      d0[i] <- dudz$zeval[i]-(rep(Et0dat$ustar,comb)[i]/dudz$dscdz[i]/0.4)
      if(!is.na(d0[i])){
        if(d0[i]<0){
          d0[i] <- NA
        }
      }
    }
    else{
      d0[i] <- NA
    }
  }
  ans <- data.frame(d0, tim, cmb = rep(1:comb, each = n))
  return(ans)
}
ntrgin <- Et0dat$time[ginmom[[1]]=="NTR"]
ginpr1 <- numderd(dudz, Et0dat, n, Et0dat$time, ntrgin)
ginpr1$d0[ginpr1$d0 >7] <- NA
ginpr1$nls <- ginmom[[7]]
ggplot(ginpr1, aes(d0, nls)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(.~cmb) + theme_bw()
x <- summary(with(ginpr1[ginpr1$cmb == 2,], lm(nls~d0)))

ntrals <- Et0datas$time[alspmom[[1]]=="NTR"]
alspr1 <- numderd(dudzas, Et0datas, nas, Et0datas$time, ntrals)
alspr1$d0[alspr1$d0 >6.5] <- NA
alspr1$nls <- alspmom[[7]]
ggplot(alspr1, aes(d0, nls)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(.~cmb)
x <- summary(with(alspr1[alspr1$cmb == 2,], lm(nls~d0)))

ntry <- Et0daty$time[yancmom[[1]]=="NTR"]
yncpr1 <- numderd(dudzy, Et0daty, ny, Et0daty$time, ntry)
yncpr1$d0[yncpr1$d0 >.3] <- NA
yncpr1$nls <- yancmom[[7]]
yncpr1$nls[yncpr1$nls == 0.05 | yncpr1$nls == 0.3] <- NA
ggplot(yncpr1, aes(d0, nls)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(.~cmb)
x <- summary(with(yncpr1[yncpr1$cmb == 2,], lm(nls~d0)))

ntrpr <- USPrrlst[[1]]$time[USPrrlst[[9]] == "NTR"]
dudzprr <- br439prf(length(USPrrlst[[9]]), 3, USPrrlst[[2]][,c(2,1,3,4)])
prrpr1 <- numderd(dudzprr, USPrrlst[[1]], length(USPrrlst[[9]]), USPrrlst[[1]][,1], ntrpr)
prrpr1$d0[prrpr1$d0 > 3] <- NA
prrpr1$nls <- USPrrlst[[1]]$d0*3
prrpr1$nls[prrpr1$nls <= 0.0051] <- NA
ggplot(prrpr1, aes(d0, nls)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(.~cmb)
x <- summary(with(yncpr1[yncpr1$cmb == 2,], lm(nls~d0)))

r2pr1 <- 0
adr2pr1 <- 0
rmpr1 <- 0
slppr1 <- 0
for(i in 1:max(prrpr1$cmb)){
  sub <- prrpr1[complete.cases(prrpr1) & prrpr1$cmb == i,-3]
  if(sum(!is.na(sub$d0))>2){
    mod <- summary(with(sub, lm(nls~d0)))
    r2pr1[i] <- mod$r.squared
    adr2pr1[i] <- mod$adj.r.squared
    rmpr1[i] <- mod$sigma
    slppr1[i] <- mod$coefficients[2]
  }
}
r2pr1 
adr2pr1
rmpr1 
slppr1
