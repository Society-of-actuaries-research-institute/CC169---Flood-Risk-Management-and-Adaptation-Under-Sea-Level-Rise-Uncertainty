# Damage curve without adaptation
loss <- c(6.5, 30.2, 54.4, 66.1, 74, 108, 135, 181.2)
rp <- c(100, 250, 500,750, 1000, 2000, 5000, 10000) #return period
# storm surge corresponding to each return period
#x=c(2752.479, 3130.493, 3468.840, 3690.924, 3860.431, 4313.775, 5024.713, 5661.912)
x=c(2852.829, 3337.924, 3794.566, 4104.750, 4346.709, 5014.500, 6116.775, 7155.684)
plot(x,loss,type='l')
x2 = x^2
model=lm(loss ~ x + x2)
summary(model)

# fit linear model
model=lm(loss ~ x)
summary(model)
coef(model) 
# a=0; b=0.0393

#find natural flood threshold u
xx=seq(1000,3000,1)
y=coef(model)[1] + coef(model)[2]*xx 
y2=y^2
u = xx[which(y2 == min(y2))]  
u #u = 2506
y2[which(y2 == min(y2))]
plot(xx,y2)
# plot(xx,y)

### FIND flood level for top cover limit

ustar = max(xx[which(y2<= 1)])
ustar # ustar =2531


# # plot damage curve
# png(file = paste(fig, "damagecurve_nyc.png", sep=""), res = 200, width = 1200, height = 800)
# par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
# plot(x, loss, xlab = 'Sea level (mm)', ylab = 'Loss (billion USD)')
# lines(x, predict(model, type='response'))
# dev.off()

#=== Damage curve for water proofing ===
# cost of wet flood proofing existing buildings: SM page 29.
# +2 feet , high effectiveness.
BCR = 0.88; NPV = -29525931/10^9

# # +4 feet = 1.219 meter, high effectiveness.
# BCR = 0.68; NPV = -161919690/10^9
out=proofing_cost(BCR, NPV, kappa,k, a, u, alpha0, scale, xi )
at2=0
bt2=out$bt
bt2
C = out$C
C


#=== Damage curve plot for 3 cases: no adaptation, dyke, elevation ===
x = seq(2000, 4500, 10)
z=x-u
lna = coef(model)[2]*z # loss from no adaptation, based on regression line
lna[lna<0] = 0
zd=x-u-k2  #effective flood height for dyke
ldk=lna  #loss under dyke project
ldk[zd<0]=0

#wetproofing
lnaw=bt2*z # loss from no adaptation for wetproofing project
lnaw[lnaw<0] = 0

ze=x-u-k1 #effective flood height for water proofing
lwp=lnaw  #loss under wetproofing project
lwp[(ze<0) & (z>0)]=kappa1*lnaw[(ze<0) & (z>0)]

png(file = paste(fig, "damagecurves_nyc.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(x[1], x[length(x)]), c(0, max(lna)), cex=0, 
     xlab='Sea level (mm)', ylab='Loss ($US billion)')
colour <- c('black', 'blue', 'black','red')
style <- c(2,1,2,1)
weight=c(1,2,1,1.5)
lines(x,lna, col=colour[1], lty=style[1],lwd=weight[1])
lines(x,ldk, col=colour[2], lty=style[2],lwd=weight[2])

lines(x,lnaw, col=colour[3], lty=style[3],lwd=weight[1])
lines(x,lwp, col=colour[4], lty=style[4],lwd=weight[4])
legend('topleft',legend=c('no adaptation', 'barrier and dyke', 'flood-proofing'),
       lty=style[-3],col=colour[-3],lwd=weight[-3], bty = "n") 
dev.off()
#--------------------------






