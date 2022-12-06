
#=== Dike Stage 1 ===
I <- I1
kappa =kappa2; k1=k1;
#mus <- c(0,mu0,0,mu0)
#mus = c(0, mu0,0,mu0)
#gammas = c(0,0,gamma0,gamma0)
mus = c(0, mu0,2*mu0,mu0)
gammas <- c(gamma0,gamma0,gamma0,2*gamma0)

result <- data.frame(matrix(NA, 3, 8))
result[,1] <- c('NPV','Total value with optionality','Difference')



for (i in 1:4){
  mu=mus[i]; gamma=gammas[i]
  
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u, #Threhsold for the damage function Before adpatation
                  a,b,u, # damage Curve 
                  u+k1, #Threhsold for the damage function After adpatation
                  a,b,u,#Damage curve
                  delta,kappa,
                  scale,xi, 
                  f=matrix(0,nt,ns-1))
  
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,5) #NPV
  
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  
  check = out$f[1,]
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  ov <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  ov$f
  
  result[2,2*i] <- round(ov$f,5)  #option value
  result[3,2*i] <- result[2,2*i] - result[1,2*i] #difference
  if(i==4){
    #write.table(data, file = paste(mpath, "gvalue", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
    png(file = paste(fig, "projectvalue_cph_dike_stage1.png", sep=""),
        res = 200,
        width = 1200,
        height = 800)
    par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
    plot(data[,1], data[,2], xlab = expression(alpha), ylab = 'Project value (Euro billion)', type='l')
    dev.off()
  }
  if(i==1){B0 <- data.frame(ov$B)}
  if(i>1){B0 <- cbind(B0, ov$B[,2])}
}

names(B0) <- c('time', 'c1', 'c2', 'c3', 'c4')
write.table(B0,
            file = paste(opath, "threshold_cph_dike_stage1", ".csv", sep=""),
            sep=",",
            row.names = FALSE,
            col.names = TRUE)

#=== Make Latex table ====
head0 <- paste("\\hline &", "$\\mu=0, \\gamma = 1.31\\%$","&&", "$\\mu=6.5, \\gamma = 1.31\\%$","&&", "$\\mu=13, \\gamma = 1.31\\%$","&&",
               "$\\mu=6.5, \\gamma = 1.31\\%$", 
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8}" ,sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)
result[,1] <- c('NPV','Total value with optionality','Difference')

library(xtable)
x.width <- xtable(result,digit = 5, label=paste('analysis_cph_dike',sep=''),
                  caption=paste("Investment analysis for Barrier and Dike Project using 
                  NPV rule and real options methods. We report (in billion Euro) the NPV, 
                  the total value with optionality and the difference of the two values. 
                  Results are reported for four cases: i) no climatic change and zero exposure growth, 
                  ii) climatic change and zero exposure growth, iii) no climatic change and exposure growth, 
                                iv) climatic change and exposure growth. 
                                All cases assume $\\sigma = 7.21 $mm.",sep=""),
                  align ="llrrrrrrr")
align(x.width) <-"llrrrrrrr"
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "analysis_cph_dike_stage1", ".tex", sep=""), scalebox = 0.80)


#== Exercise boundary ===
data <- read.csv(paste(opath, "threshold_cph_dike_stage1.csv", sep=""), sep=",", header=TRUE)

#idx = seq(14,
#          114,
#          by = 1)

#rownames(data) = data$time
#data = data[idx*12,]
data <- data[data$time>10*12,]
data <- data[data$time<10*100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'green', 'blue', 'red')
style <- seq(1,4,1)

png(file = paste(fig, "exerciseboundary_cph_dike_stage1.png", sep=""),
    res = 200,
    width = 1400,
    height = 1000)

xtick = seq(20,100,by =20)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:5]), max(data[,2:5])+300), cex=0, 
     xlab='Year',
     ylab='Investment threshold (mm)',
     xaxt = "n",
     ylim = c(-400,1200) )

axis(sid = 1,
     at = seq(200, 1000, by = 200),
     labels = xtick)
for (i in 1:4){
  lines(data[,1], data[,1+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c(expression(paste(mu, "=0    ", gamma, "=1.31%" )),
                           expression(paste(mu, "=6.5 ", gamma, "=1.31%" )),
                           expression(paste(mu, "=13    ", gamma, "=1.31%" )),
                           expression(paste(mu, "=6.5 ", gamma, "=2.6%" ))),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
dev.off()
#----------------------



#== Dike Second Stage ===
I <- I2; k=k2
kappa =kappa1; 
a0=a1; b0=b1; 

result <- data.frame(matrix(NA, 3, 8))
result[,1] <- c('NPV','Total value with optionality','Difference')

#=== Investment Analysis ===
#mus <- c(0,mu0,0,mu0)
mus = c(0, mu0,2*mu0,mu0)
gammas <- c(gamma0,gamma0,gamma0,2*gamma0)


for (i in 1:4){
#i=4
  mu=mus[i]; gamma=gammas[i]
  
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u+k1,
                  a,b, u,
                  u+ k1 + k2,
                  a,b,u,
                  delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds  #interpolating for current value of alpha0
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,5) #NPV
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  ov <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  ov$f
  
  result[2,2*i] <- round(ov$f,5) #option value
  result[3,2*i] <- result[2,2*i] - result[1,2*i] #difference
  if(i==4){
    #write.table(data, file = paste(mpath, "gvalue", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
    png(file = paste(fig, "projectvalue_cph_dike_stage_2.png", sep=""), res = 200, width = 1200, height = 800)
    par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
    plot(data[,1], data[,2], xlab = expression(alpha), ylab = 'Project value (Euro billion)', type='l')
    dev.off()
  }
  if(i==1){B0 <- data.frame(ov$B)}
  if(i>1){B0 <- cbind(B0, ov$B[,2])}
}

names(B0) <- c('time', 'c1', 'c2', 'c3', 'c4')
write.table(B0,
            file = paste(opath, "threshold_cph_stage2", ".csv", sep=""),
            sep=",",
            row.names = FALSE,
            col.names = TRUE)

#=== Make Latex table ====
head0 <- paste("\\hline &", "$\\mu=0, \\gamma = 1.31\\%$","&&", "$\\mu=6.5, \\gamma = 1.31\\%$","&&", "$\\mu=13, \\gamma = 1.31\\%$","&&",
               "$\\mu=6.5, \\gamma = 1.31\\%$", 
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8}" ,sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 5, label=paste('analysis_cph_dike_stage2',sep=''),
                  caption=paste("Investment analysis for dike  using 
                  NPV rule and real options methods. We report (in billion Euro) the NPV, 
                  the total value with optionality and the difference of the two values. 
                  Results are reported for four cases: i) no climatic change and zero exposure growth, 
                  ii) climatic change and zero exposure growth, iii) no climatic change and exposure growth, 
                                iv) climatic change and exposure growth. 
                                All cases assume $\\sigma = 7.21 $mm.",sep=""),
                  align ="llrrrrrrr")
align(x.width) <-"llrrrrrrr"
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "analysis_cph_dike_stage2", ".tex", sep=""), scalebox = 0.80)


#== Exercise boundary ===
data <- read.csv(paste(opath, "threshold_cph_stage2.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>10*12,]
data <- data[data$time<10*100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'green', 'blue', 'red')
style <- seq(1,4,1)
xtick = seq(20,100,by = 20 )
png(file = paste(fig, "exerciseboundary_cph_dike_stage2.png", sep=""),
    res = 200,
    width = 1400,
    height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:5]), max(data[,2:5])+350), cex=0, 
     xlab='Year',
     ylab='Investment threshold (mm)',
     xaxt = "n",
     ylim = c(-400,1200))
axis(sid = 1,
     at = seq(200, 1000, by = 200),
     labels = xtick)
for (i in 1:4){
  lines(data[,1], data[,1+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c(expression(paste(mu, "=0    ", gamma, "=1.31%" )),
                           expression(paste(mu, "=6.5 ", gamma, "=1.31%" )),
                           expression(paste(mu, "=13    ", gamma, "=1.31%" )),
                           expression(paste(mu, "=6.5 ", gamma, "=2.6%" ))),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
dev.off()
#----------------------





#== Dike All in one stage ==#

I <- I1and2; k=k2
kappa =kappa1; 
a0=a1; b0=b1; 

result <- data.frame(matrix(NA, 3, 8))
result[,1] <- c('NPV','Total value with optionality','Difference')

#=== Investment Analysis ===
#mus <- c(0,mu0,0,mu0)
#mus = c(0.79/12, mu0,0.79/12,mu0)

#gammas <- c(0,0,gamma0,gamma0)


mus = c(0, mu0,2*mu0,mu0)
gammas <- c(gamma0,gamma0,gamma0,2*gamma0)

for (i in 1:4){
  #i=4
  mu=mus[i]; gamma=gammas[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,
                  a0, b0, u,
                  u + k1 + k2,
                  a1,b1,u,
                  delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds  #interpolating for current value of alpha0
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,5) #NPV
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  ov <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  ov$f
  result[2,2*i] <- round(ov$f,5) #option value
  result[3,2*i] <- result[2,2*i] - result[1,2*i] #difference
  if(i==4){
    #write.table(data, file = paste(mpath, "gvalue", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
    png(file = paste(fig, "projectvalue_cph_dike_all_in_one.png", sep=""),
        res = 200,
        width = 1200,
        height = 800)
    par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
    plot(data[,1], data[,2], xlab = expression(alpha), ylab = 'Project value (Euro billion)', type='l')
    dev.off()
  }
  if(i==1){B0 <- data.frame(ov$B)}
  if(i>1){B0 <- cbind(B0, ov$B[,2])}
}

names(B0) <- c('time', 'c1', 'c2', 'c3', 'c4')
write.table(B0,
            file = paste(opath, "threshold_cph_all_in_one", ".csv", sep=""),
            sep=",",
            row.names = FALSE,
            col.names = TRUE)


#=== Make Latex table ====
head0 <- paste("\\hline &", "$\\mu=0, \\gamma = 1.31\\%$","&&", "$\\mu=6.5, \\gamma = 1.31\\%$","&&", "$\\mu=13, \\gamma = 1.31\\%$","&&",
               "$\\mu=6.5, \\gamma = 1.31\\%$",  
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8}" ,sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 5, label=paste('analysis_cph_all_in_one',sep=''),
                  caption=paste("Investment analysis for Dike using 
                  NPV rule and real options methods. We report (in billion Euros) the NPV, 
                  the total value with optionality and the difference of the two values. 
                  Results are reported for four cases: i) no climatic change and zero exposure growth, 
                  ii) climatic change and zero exposure growth, iii) no climatic change and exposure growth, 
                                iv) climatic change and exposure growth. 
                                All cases assume $\\sigma = 7.21 $mm.",sep=""),
                  align ="llrrrrrrr")
align(x.width) <-"llrrrrrrr"
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),
      add.to.row = addtorow, 
      sanitize.text=function(x){x},
      caption.placement = 'top',
      table.placement='H',
      floating = TRUE,
      include.rownames=FALSE,
      include.colnames=FALSE,
      file = paste(tab, "analysis_cph_dike_all_in_one", ".tex", sep=""),
      scalebox = 0.80)


#== Exercise boundary ===
data <- read.csv(paste(opath, "threshold_cph_all_in_one.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>10*12,]
data <- data[data$time<10*100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'green', 'blue', 'red')
style <- seq(1,4,1)
xtick = seq(20,100,by = 20)
png(file = paste(fig, "exerciseboundary_cph_all_in_one.png", sep=""),
    res = 200,
    width = 1400,
    height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:5]), max(data[,2:5])+350), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)',
     xaxt = "n",
     ylim = c(-400,1200) )
axis(sid = 1,
     at = seq(200, 1000, by = 200),
     labels = xtick)
for (i in 1:4){
  lines(data[,1], data[,1+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c(expression(paste(mu, "=0    ", gamma, "=1.31%" )),
                           expression(paste(mu, "=6.5 ", gamma, "=1.31%" )),
                           expression(paste(mu, "=13    ", gamma, "=1.31%" )),
                           expression(paste(mu, "=6.5 ", gamma, "=2.6%" ))),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
dev.off()



#== Sequencing: proofing-dike and dike-proofing
# proofing is indexed 1, dike is indexed 2.

result <- data.frame(matrix(NA, 3, 12))
result[,1] <- c('NPV','Total value with optionality','Difference')
#=== Investment Analysis ===
# dike after proofing:proportion (kappa1*b1)/b2 of total exposure has been protected by proofing, 
# so the relevant damage curve for dike is [1-(kappa1*b1)/b2]b2 = b2-kappa1*b1
# proofing after dike: benefit of proofing = D(k2,u,L) - D(k1,u,L)
# so put flood threshold before investment as k2
# Investment cost for sequencing
I1 = (160000 +  1371 * (k1^2))/1000000000
I12 = (160000 +  1371 * (k2)^2 +   1371 * (2*k1 * k2))/1000000000
I2 = (160000 +  1371 * (k2^2))/1000000000
I21 = (160000 +  1371 * (k1)^2 + 1371 * (2*k1 * k2))/1000000000

# INvestment Cost of building the k1+k2 heihgt dike in one go ( saving on fixed cost)

#I12 = (1600000 +  1371 * (k1^2) + 1371 * (2*k1^2 ) + 1371 * (k1^2)) /1000000000 




# order of vector below: c(proofing, dike after proofing, dike, proofing after dike)
Is <- c(I1, I12, I2, I12)
# flood threshold before investment
x0 <- c(u,
        u+k1,
        u,
        u+k2) 
# flood threshold after investment
x1 <- c(u+k1,
        u+k1 + k2,
        u+k2,
        u+k1 +k2)
# effectiveness of investment
kappas=c(kappa1,kappa2,kappa2,kappa1)
# intercept of damage curve before investment
u0 <- c(u,u,u,u)
# intercept of damage curve after investment, see Figure 1.
u1 <- c(u,u,u,u)
# quadratic term before investment
as0 <- c(a1,a2,a2,a1)
# slope term before investment
#bs0 <- c(b1,b2-kappa1*b1,b2,b1)
bs0 = c(b1,b2,b2,b1)
# quadratic term after investment
as1 <- c(a1,a2,a2,a1)
# slope term after investment
bs1 <- c(b1,b2,b2,b1)
for (i in 1:4){
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu0, sigma, r,gamma0, nt, ns, Smin, Smax,T,
                  x0[i],as0[i],bs0[i],u0[i], x1[i],as1[i],bs1[i],u1[i],
                  delta,kappas[i], scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - Is[i]  # NPV
  if(i<3){
    result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - Is[i],5)
  }else{
    result[1,2+2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - Is[i],5)
  }
  
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  out1 <- .Fortran("binop", mu0, sigma,gamma0, alpha0, r, nto, T, N, V,Is[i], f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out1$f
  if(i<3){
    result[2,2*i] <- round(out1$f,5)
    result[3,2*i] <- result[2,2*i] - result[1,2*i]
  }else{
    result[2,2+2*i] <- round(out1$f,5)
    result[3,2+2*i] <- result[2,2+2*i] - result[1,2+2*i]
  }
  if(i==1){B0 <- data.frame(out1$B)}
  if(i>1){B0 <- cbind(B0, out1$B[,2])}
}
names(B0) <- c('time', 'c1', 'c12', 'c2', 'c21')
write.table(B0, file = paste(opath, "threshold_portfolio_cph", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
result[,6]=unlist(mapply(function(i) result[i,2]+result[i,4], 1:3,SIMPLIFY=FALSE))
result[,12]=unlist(mapply(function(i) result[i,8]+result[i,10], 1:3,SIMPLIFY=FALSE))


#==== Investment threshold ===
data <- read.csv(paste(opath, "threshold_portfolio_cph.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>12*10,]
data <- data[data$time<12*100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'red')
style <- seq(1,2,1)
png(file = paste(fig, "exerciseboundary_cph_portfolio_12.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:3]), max(data[,2:3])+500), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)')
for (i in 1:2){
  lines(data[,1], data[,1+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c("Dike Stage 1", "Dike Stage 2"),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
dev.off()

#==== Investment threshold ===
data <- read.csv(paste(opath, "threshold_portfolio_cph.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>12*10,]
data <- data[data$time<12*100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('red', 'black')
style <- seq(1,2,1)
png(file = paste(fig, "exerciseboundary_cph_portfolio_21.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,4:5]), max(data[,4:5])+1000), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)')
for (i in 1:2){
  lines(data[,1], data[,3+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c("Dike Stage 2","Dike Stage 1"),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
dev.off()
#----------------------

#=== Make Latex table ====
head0 <- paste("\\hline &", "\\multicolumn{5}{c}{Dike Stage 1 then Dike Stage 2}", "&&",
               "\\multicolumn{5}{c}{Dike Stage 2 then Dike Stage 1}", "\\\\","\n",
               "\\cmidrule{2-6} \\cmidrule{8-12}" ,sep="")

head1 <- paste("&", "Dike Stage 1", "&&","Dike Stage 2", '&&', "Total",'&&', "Dike Stage 2", "&&","Dike Stage 1", '&&', "Total",
               "\\\\","\n","\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10} \\cmidrule{12-12}", sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- -1
addtorow$pos[[2]] <- 0
addtorow$command <- c(head0,head1)

library(xtable)
x.width <- xtable(result,digit = 5, label=paste('analysis_cph_portfolio_de',sep=''),
                  caption=paste("Investment analysis for investment sequences using 
                  NPV rule and real options methods. We report (in billion euro) the NPV, 
                  the total value with optionality and the difference of the two values.",sep=""),
                  align ="llrrrrrrrrrrr")
align(x.width) <-"llrrrrrrrrrrr"
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "analysis_cph_portfolio_de", ".tex", sep=""), scalebox = 0.65)

I1 = (160000 +  1371 * (k1^2))/1000000000
#I1 = (1371 * (k1^2))/1000000000
#I1 = 0.973
I2 = (160000 +  1371 * (k2^2) + 1371 * (2*k1*k2 ) )/1000000000
#I2 = 2.919




I12 = (160000 +  1371 *(k1+k2)^2)/1000000000




