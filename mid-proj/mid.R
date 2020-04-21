# save(conv, file = "conv.RData")
# save(org, file = "org.RData")

ggplot(conv, aes(x=Date, y=Total.Volume/1e4)) + 
  geom_line(aes(color="#2983bb")) +
  geom_line(data=org, aes(x=Date, y=Total.Volume/1e3, color="#ff9900")) +
  xlab('Date') +
  ylab('Conventional total Vol. [K]')  +
  scale_colour_manual(name = 'type', 
                      values =c("#2983bb"="#2983bb","#ff9900"="#ff9900"), 
                      labels = c('conventional','organic')) + 
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Organic total Vol. [K]"))

date <- conv$Date
u1 <- log(conv$AveragePrice)
# u1_ts <- ts(u1,start=2015,frequency=1)
# spec1 <- spectrum(u1_ts, spans=c(3,5,3), 
#                   main="Smoothed periodogram of average price")
# date <- seq(from=2015,length=length(u1),by=1/52)
u1_ts <- ts(u1,start=2015,frequency=52)
u_low <- ts(loess(u1~date,span=0.5)$fitted,start=2015,frequency=52)
u_hi <- ts(u1 - loess(u1~date,span=0.1)$fitted,start=2015,frequency=52)
u_cycles <- u1 - u_hi - u_low
plot(ts.union(u1, u_low,u_hi,u_cycles),
     main="Decomposition of average price as trend + noise + cycles")


which.max(spec1$spec)

ggplot() + 
  geom_line(aes(x=spec1$freq, y=spec1$spec), color="#74787a") + 
  labs(title="Smoothed periodogram of average price", 
       # subtitle="Year 2015 - 2018", 
       x="frequency [cycles per week]",
       y="spectrum") + 
  scale_y_continuous(trans = 'log10') + 
  geom_vline(xintercept = spec1$freq[41], colour="#2983bb", linetype=2) + 
  geom_vline(xintercept = spec1$freq[3], colour="#ff9900", linetype=2)


date <- conv$Date
u2 <- conv$Total.Volume
u2_ts <- ts(u2,start=2015,frequency=1)
spec2 <- spectrum(u2_ts,spans=c(3,5, 3), plot=FALSE,
                  main="Smoothed periodogram of average price")
ggplot() + 
  geom_line(aes(x=spec2$freq, y=spec2$spec/1e13), color="#c04851") + 
  labs(title="Spectrum comparison", 
       subtitle="Total volume (red line) vs. Average price (blue line)", 
       x="frequency [cycles per week]",
       y="Average price") + 
  # scale_y_continuous(trans = 'log10') + 
  geom_vline(xintercept = spec2$freq[41], colour="#f1939c", linetype=2) + 
  geom_vline(xintercept = spec2$freq[4], colour="#f1939c", linetype=2)  + 
  geom_line(aes(x=spec1$freq, y=spec1$spec), color="#144a74") + 
  scale_y_continuous(trans = 'log10') + 
  geom_vline(xintercept = spec1$freq[41], colour="#2983bb", linetype=3) + 
  geom_vline(xintercept = spec1$freq[3], colour="#2983bb", linetype=3) + 
  scale_y_continuous(trans = 'log10', 
                     sec.axis = sec_axis(~.*1e13, name = "Total volume"))

cor.test(conv$AveragePrice, conv$Total.Volume, 
         method=c("pearson", "kendall", "spearman"))


acf(conv$AveragePrice,main="ACF of ARMA(3,2) residuals", lag.max = 169)

conv$residual = c(conv$residual)
ggplot(conv, aes(x=Date, y=residual)) + 
  geom_point()


aic_reg_table <- function(P,Q){
  table <- matrix(NA,(P+1),(Q+1))
  for(p in 0:P) {
    for(q in 0:Q) {
      try(table[p+1,q+1] <- arima(x=log(conv$AveragePrice), 
                                  order=c(p,0,q),
                                  seasonal=list(order=c(1,0,0),period=52), 
                                  xreg=log(conv$Total.Volume))$aic)
    }
  }
  dimnames(table) <- list(paste("AR",0:P, sep=""),paste("MA",0:Q,sep=""))
  table
}

temp = aic_reg_table(5,5)
kable(temp,digits=2)

arma10 <- arima(x=log(conv$Total.Volume), 
                order=c(1,0,0),
                seasonal=list(order=c(1,0,0),period=52), 
                xreg=log(conv$AveragePrice))

acf0 <- acf(arma10$residuals,main="ACF of ARMA(1,1) residuals", lag.max = 169)




