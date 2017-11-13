spaov <-
function(resp, nestfct, crossfct, block, plots = TRUE){
  #create factors
  nestfct <- as.factor(nestfct)
  crossfct <- as.factor(crossfct)
  block <- as.factor(block)
  
  #calculate degrees of freedom
  dfnest <- length(levels(nestfct)) - 1
  dffct <- length(levels(crossfct)) - 1
  dfblock <- length(levels(block)) - (1 + dfnest)
  dfinteract <- dfnest * dffct
  dfT <- length(resp)
  dfE <- dfT - sum(dfnest, dffct, dfblock, dfinteract, 1)
  
  
  #calculate SS
  lm1 <- lm(resp ~ nestfct + block + nestfct*crossfct)
  SSE <- sum(lm1$residuals^2)
  SST <- sum(resp^2)
  aov1 <- aov(resp ~ nestfct + block + nestfct*crossfct)
  SSnest <- summary(aov1)[[1]]$`Sum Sq`[1]
  SSblock <- summary(aov1)[[1]]$`Sum Sq`[2]
  SSfct <- summary(aov1)[[1]]$`Sum Sq`[3]
  SSinteract <- summary(aov1)[[1]]$`Sum Sq`[4]
  
  #calculate MS
  MSE <- SSE/dfE
  MSnest <- SSnest/dfnest
  MSblock <- SSblock/dfblock
  MSfct <- SSfct/dffct
  MSinteract <- SSinteract/dfinteract
  
  #F stat
  Fvnest <- MSnest/MSblock
  Fvblock <- MSblock/MSE
  Fvfct <- MSfct/MSE
  Fvinteract <- MSinteract/MSE
  
  #P value
  Pvnest <- pf(Fvnest, dfnest, dfblock, lower.tail = FALSE)
  Pvblock <- pf(Fvblock, dfblock, dfE, lower.tail = FALSE)
  Pvfct <- pf(Fvfct, dffct, dfE, lower.tail = FALSE)
  Pvinteract <- pf(Fvinteract, dfinteract, dfE, lower.tail = FALSE)
  
  #put it all in a table
  source <- c("Nest Factor", "Block", "Crossed Factor", "Interaction", "Residuals", "Total")
  DF <- c(dfnest, dfblock, dffct, dfinteract, dfE, dfT)
  SS <- c(SSnest, SSblock, SSfct, SSinteract, SSE, SST)
  MS <- c(round(MSnest, digits = 6),
          round(MSblock, digits = 6),
          round(MSfct, digits = 6),
          round(MSinteract, digits = 6),
          round(MSE, digits = 6), "")
  Fv <- c(round(Fvnest, digits = 6),
          round(Fvblock, digits = 6),
          round(Fvfct, digits = 6),
          round(Fvinteract, digits = 6),
          "", "")
  Pv <- c(round(Pvnest, digits = 6),
          round(Pvblock, digits = 6),
          round(Pvfct, digits = 6),
          round(Pvinteract, digits = 6),
          "", "")
  
  tabeel <- data.frame(row.names = source,
                       "Df" = DF,
                       "Sum Sq" = SS,
                       "Mean Sq" = MS,
                       "F value" = Fv,
                       "P value" = Pv)
  
  if (plots == TRUE) {
    print(tabeel)
    #Autobuild diagnostic plots
    par(mfrow = c(1, 2))
    plot(aov1, which = c(1:2))
  }
  else if (plots == FALSE) {
    print(tabeel)
  }

}
