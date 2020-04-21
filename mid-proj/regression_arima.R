aic_reg_table <- function(P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      try(table[p+1,q+1] <- arima(x=log(conv$Total.Volume), 
                                  order=c(p,0,q),
                                  xreg=log(conv$AveragePrice))$aic)
    }
  }
  dimnames(table) <- list(paste("AR",0:P, sep=""),paste("MA",0:Q,sep=""))
  table
}

temp = aic_reg_table(5,5)
kable(temp,digits=2)

arma11 <- arima(x=log(conv$Total.Volume), 
                order=c(0,1,1),
                xreg=log(conv$AveragePrice))

acf1 <- acf(arma11$residuals,main="ACF of ARMA(1,1) residuals", lag.max = 169)

acf1$acf[53]

sarma11 <- arima(x=log(conv$Total.Volume), 
                 order=c(1,0,1),
                 seasonal=list(order=c(1,0,0),period=52), 
                 xreg=log(conv$AveragePrice))
sarma11

acf2 <- acf(sarma11$residuals,main="ACF of SARMA(1,1) residuals", lag.max = 169)

ggplot() + 
  geom_point(aes(x=conv$Date, y=c(sarma11$residuals)), color="#5e7987") +
  xlab('Date') +
  ylab('Residuals')


xe = sarma11$residuals

aic_table <- function(data,P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      try(table[p+1,q+1] <- arima(data,order=c(p,0,q),method = "ML")$aic)
    }
  }
  dimnames(table) <- list(paste("AR",0:P, sep=""),paste("MA",0:Q,sep=""))
  table
}

temp_aic_table = aic_table(xe,5,5) 




