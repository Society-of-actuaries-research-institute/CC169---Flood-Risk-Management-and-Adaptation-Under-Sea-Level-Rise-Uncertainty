# losses to all buildings in study region by Wang 2015
x=c(1873,2479,4161,8871) 
loss = c(0.012, 0.12, 0.45, 1.925) # from Figure 3 in Wang 2015 and Figure 10 in Wang 2014.
plot(x,loss,type='l')
x2 = x^2
model=lm(loss ~ x + x2)
summary(model)

#find natural flood threshold u
xx=seq(1600,3000,1)
y=coef(model)[1] + coef(model)[2]*xx + coef(model)[3]*xx^2
y2=y^2
xx[which(y2 == min(y2))]
y2[which(y2 == min(y2))]
plot(xx,y2)
# plot(xx,y)


# plot damage curve
png(file = paste(fig, "damagecurve_qld.png", sep=""), res = 200, width = 1200, height = 800)
par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
plot(x, loss, xlab = 'Sea level (mm)', ylab = 'Loss (billion AUD)')
lines(x, predict(model, type='response'))
dev.off()











