  fixkod <- function(richstb, fkn2017cat, Et0dat){
  top <- data.frame(dudzfwd = richstb[[1]][,1],dudzcd = richstb[[1]][,3],dudzbwd = richstb[[1]][,2], ustar = Et0dat$ustar)
  top <- top[fkn2017cat == "NTR",]
  mtop <- melt(top, id = "ustar")
  ggplot(mtop, aes(ustar, value)) + facet_grid(.~variable) + geom_point(alpha = 0.1) + geom_smooth(method = "lm", formula = y~0+x) + geom_smooth(method = "lm", col = "red", formula = y~x) + theme_bw() + xlab(expression(paste(u["*"]," m/s"))) + ylab(expression(paste(partialdiff,"u", "/", partialdiff,"z"," 1/s")))
  
  lmdudzufwd <- lm(dudzfwd~ 0 + ustar, data = top)
  lmdudzubwd <- lm(dudzbwd~ 0 +ustar, data = top)
  lmdudzucd <- lm(dudzcd~ 0 +ustar, data = top)
  mtop$slp <- lmdudzubwd$coefficients[1]
  mtop$slp[mtop$variable == "dudzfwd"] <- lmdudzufwd$coefficients[1]
  mtop$slp[mtop$variable == "dudzcd"] <- lmdudzucd$coefficients[1]
  mtop$slp <- round(mtop$slp, 2)
  mtop$z <- 13.05
  mtop$z[mtop$variable == "dudzfwd"] <- 8.8
  mtop$z[mtop$variable == "dudzcd"] <- 11.25
  slps <- unique(mtop$slp)
  zees <- unique(mtop$z)
  fun1 <- function(slp, heigh, d0){
    k <- (slp*(heigh-d0))^-1
    return(k)
  }
  d0top <- seq(0.1,5, by = 0.1)
  topkc <- data.frame(d0 = rep(d0top, 3), sch = rep(c("Forward", "Central", "Backward"), each = length(d0top)), 
                      k = c(fun1(slps[1], zees[1],d0top), fun1(slps[2], zees[2],d0top), fun1(slps[3], zees[3],d0top)))
  ggplot(topkc, aes(d0, k))  + annotate("rect", xmin=c(0), xmax=c(5), ymin=c(0.35) , ymax=c(0.43), color="black", fill="grey", alpha = 0.2) +
    geom_text(aes(x  = 1, y = 0.375, label = "Foken Table 2.11 range"), color = "black") + geom_line(aes(col = sch)) + 
    theme_bw() + ggtitle("k and d0 relationship") + geom_hline(yintercept = 0.4) + geom_vline(xintercept = 2*6.8/3)
  mtop <- melt(top, id = "ustar")
  lab <- paste("At z = 13.05\n, partialdiffu%/%partialdiffz ==", as.character(round(lmdudzubwd$coefficients[2]),2), "%*%u[*] + ",  as.character(round(lmdudzubwd$coefficients[1]),2), "\nRMSE == ", round(summary(lmdudzubwd)$r.squared,2), "; R2 == ", round(summary(lmdudzubwd)$sigma,2))
  ann_text <- data.frame(ustar = 0, value = 0.6,lab = lab, variable = factor("dudzbwd",levels = c("dudzfwd","dudzcd","dudzbwd")))
  
  mtop$int <- lmdudzubwd$coefficients[1]
  mtop$int[mtop$variable == "dudzfwd"] <- lmdudzufwd$coefficients[1]
  mtop$int[mtop$variable == "dudzcd"] <- lmdudzucd$coefficients[1]
  mtop$int <- round(mtop$int, 2)
  
  mtop$r2 <- summary(lmdudzubwd)$r.squared
  mtop$r2[mtop$variable == "dudzfwd"] <- summary(lmdudzufwd)$r.squared
  mtop$r2[mtop$variable == "dudzcd"] <- summary(lmdudzucd)$r.squared
  mtop$r2 <- round(mtop$r2, 2)
  mtop$rse <- summary(lmdudzubwd)$sigma
  mtop$rse[mtop$variable == "dudzfwd"] <- summary(lmdudzufwd)$sigma
  mtop$rse[mtop$variable == "dudzcd"] <- summary(lmdudzucd)$sigma
  mtop$rse <- round(mtop$rse, 2)
  }