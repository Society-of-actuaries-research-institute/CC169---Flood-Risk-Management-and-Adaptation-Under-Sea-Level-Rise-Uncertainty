

#=== Barrier and Dike ===
I <- I2
kappa =kappa2; k=k2;
mus <- c(0,mu0,2*mu0,mu0)
gammas <- c(gamma0,gamma0,gamma0,2*gamma0)
result <- data.frame(matrix(NA, 3, 8))
result[,1] <- c('NPV','Total value with optionality','Difference')
for (i in 1:4){
  mu=mus[i]; gamma=gammas[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,a,b,u,u+k,a,b,u, delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,3) #NPV
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  ov <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, To, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  ov$f
  result[2,2*i] <- round(ov$f,3)  #option value
  result[3,2*i] <- result[2,2*i] - result[1,2*i] #difference
  if(i==4){
    #write.table(data, file = paste(mpath, "gvalue", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
    png(file = paste(fig, "projectvalue_nyc_dike.png", sep=""), res = 200, width = 1200, height = 800)
    par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
    plot(data[,1], data[,2], xlab = expression(alpha), ylab = 'Project value ($US billion)', type='l')
    dev.off()
  }
  if(i==1){B0 <- data.frame(ov$B)}
  if(i>1){B0 <- cbind(B0, ov$B[,2])}
}
names(B0) <- c('time', 'c1', 'c2', 'c3', 'c4')
write.table(B0, file = paste(opath, "threshold_nyc_dike", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)

#=== Make Latex table ====
head0 <- paste("\\hline &", "$\\mu=0, \\gamma = 1.31\\%$","&&", "$\\mu=6.5, \\gamma =  1.31\\%$","&&", "$\\mu=13, \\gamma =  1.31\\%\\%$","&&",
               "$\\mu=6.5, \\gamma = 2.6\\%$", 
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8}" ,sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)
result[,1] <- c('NPV','Total value with optionality','Difference')

library(xtable)
x.width <- xtable(result,digit = 2, label=paste('analysis_nyc_dike',sep=''),
                  caption=paste("Investment analysis for Barrier and Dike Project using 
                  NPV rule and real options methods. We report (in billion USD) the NPV, 
                  the total value with optionality and the difference of the two values. 
                  Results are reported for four cases: i) no climatic change and zero exposure growth, 
                  ii) climatic change and zero exposure growth, iii) no climatic change and exposure growth, 
                                iv) climatic change and exposure growth. 
                                All cases assume $\\sigma = 25 $mm.",sep=""),
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
      file = paste(tab, "analysis_nyc_ike", ".tex", sep=""), scalebox = 0.80)


#== Exercise boundary ===
data <- read.csv(paste(opath, "threshold_nyc_dike.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>3,]
data <- data[data$time<100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'green', 'blue', 'red')
style <- seq(1,4,1)
png(file = paste(fig, "exerciseboundary_nyc_dike.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:5]), max(data[,2:5])+300), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)')
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



#== Wet-proofing project ===
I <- I1; k=k1
kappa =kappa1; 
a0=a1; b0=b1; 

result <- data.frame(matrix(NA, 3, 8))
result[,1] <- c('NPV','Total value with optionality','Difference')

#=== Investment Analysis ===
mus <- c(0,mu0,2*mu0,mu0)
gammas <- c(gamma0,gamma0,gamma0,2*gamma0)

for (i in 1:4){
#i=4
  mu=mus[i]; gamma=gammas[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,a0,b0,u,u+k,a1,b1,u, delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds  #interpolating for current value of alpha0
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,3) #NPV
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  ov <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  ov$f
  result[2,2*i] <- round(ov$f,3) #option value
  result[3,2*i] <- result[2,2*i] - result[1,2*i] #difference
  if(i==4){
    #write.table(data, file = paste(mpath, "gvalue", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
    png(file = paste(fig, "projectvalue_nyc_proofing.png", sep=""), res = 200, width = 1200, height = 800)
    par(mfrow=c(1,1), mar=c(4,4,2,2)+0.1) #Margin count from xaxis, yaxis and so on
    plot(data[,1], data[,2], xlab = expression(alpha), ylab = 'Project value (USD billion)', type='l')
    dev.off()
  }
  if(i==1){B0 <- data.frame(ov$B)}
  if(i>1){B0 <- cbind(B0, ov$B[,2])}
}

names(B0) <- c('time', 'c1', 'c2', 'c3', 'c4')
write.table(B0, file = paste(opath, "threshold_nyc_proofing", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)

#=== Make Latex table ====
head0 <- paste("\\hline &", "$\\mu=0, \\gamma = 1.31\\%$","&&", "$\\mu=6.5, \\gamma =  1.31\\%$","&&", "$\\mu=13, \\gamma =  1.31\\%\\%$","&&",
               "$\\mu=6.5, \\gamma = 2.6\\%$", 
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8}" ,sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 2, label=paste('analysis_nyc_proofing',sep=''),
                  caption=paste("Investment analysis for wet flood-proofing project using 
                  NPV rule and real options methods. We report (in billion USD) the NPV, 
                  the total value with optionality and the difference of the two values. 
                  Results are reported for four cases: i) no climatic change and zero exposure growth, 
                  ii) climatic change and zero exposure growth, iii) no climatic change and exposure growth, 
                                iv) climatic change and exposure growth. 
                                All cases assume $\\sigma = 25 $mm.",sep=""),
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
      file = paste(tab, "analysis_nyc_proofing", ".tex", sep=""), scalebox = 0.80)


#== Exercise boundary ===
data <- read.csv(paste(opath, "threshold_nyc_proofing.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>6,]
data <- data[data$time<100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'green', 'blue', 'red')
style <- seq(1,4,1)
png(file = paste(fig, "exerciseboundary_nyc_proofing.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:5]), max(data[,2:5])+350), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)')
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


#== Sequencing: proofing-dike and dike-proofing
# proofing is indexed 1, dike is indexed 2.

result <- data.frame(matrix(NA, 3, 12))
result[,1] <- c('NPV','Total value with optionality','Difference')
#=== Investment Analysis ===
# dike after proofing:proportion (kappa1*b1)/b2 of total exposure has been protected by proofing, 
# so the relevant damage curve for dike is [1-(kappa1*b1)/b2]b2 = b2-kappa1*b1
# proofing after dike: benefit of proofing = D(k2,u,L) - D(k1,u,L)
# so put flood threshold before investment as k2
# order of vector below: c(proofing, dike after proofing, dike, proofing after dike)
Is <- c(I1, I2, I2, I1)
# flood threshold before investment
x0 <- c(u,u,u,u+k2) 
# flood threshold after investment
x1 <- c(u+k1,u+k2,u+k2,u+k1)
# effectiveness of investment
kappas=c(kappa1,kappa2,kappa2,kappa1)
# intercept of damage curve before investment
u0 <- c(u,u,u,u)
# intercept of damage curve after investment, see Figure 1.
u1 <- c(u,u,u,u)
# quadratic term before investment
as0 <- c(a1,a2,a2,a1)
# slope term before investment
bs0 <- c(b1,b2-kappa1*b1,b2,b1)
# quadratic term after investment
as1 <- c(a1,a2,a2,a1)
# slope term after investment
bs1 <- c(b1,b2-kappa1*b1,b2,b1)
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
    result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - Is[i],3)
  }else{
    result[1,2+2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - Is[i],3)
  }
  
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  out1 <- .Fortran("binop", mu0, sigma,gamma0, alpha0, r, nto, T, N, V,Is[i], f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out1$f
  if(i<3){
    result[2,2*i] <- round(out1$f,3)
    result[3,2*i] <- result[2,2*i] - result[1,2*i]
  }else{
    result[2,2+2*i] <- round(out1$f,3)
    result[3,2+2*i] <- result[2,2+2*i] - result[1,2+2*i]
  }
  if(i==1){B0 <- data.frame(out1$B)}
  if(i>1){B0 <- cbind(B0, out1$B[,2])}
}
names(B0) <- c('time', 'c1', 'c12', 'c2', 'c21')
write.table(B0, file = paste(opath, "threshold_portfolio", ".csv", sep=""), sep=",", row.names = FALSE, col.names = TRUE)
result[,6]=unlist(mapply(function(i) result[i,2]+result[i,4], 1:3,SIMPLIFY=FALSE))
result[,12]=unlist(mapply(function(i) result[i,8]+result[i,10], 1:3,SIMPLIFY=FALSE))


#==== Investment threshold ===
data <- read.csv(paste(opath, "threshold_portfolio.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>9,]
data <- data[data$time<100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'red')
style <- seq(1,2,1)
png(file = paste(fig, "exerciseboundary_nyc_portfolio_ed.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,2:3]), max(data[,2:3])+500), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)')
for (i in 1:2){
  lines(data[,1], data[,1+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c("Flood-proofing", "Barrier and dike"),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
dev.off()

#==== Investment threshold ===
data <- read.csv(paste(opath, "threshold_portfolio.csv", sep=""), sep=",", header=TRUE)
data <- data[data$time>6,]
data <- data[data$time<100,]
data$id <- NA
for (i in 1:(nrow(data)/2)){data$id[(i-1)*2] <- 1}
data <- na.omit(data)
colour <- c('black', 'red')
style <- seq(1,2,1)
#png(file = paste(fig, "exerciseboundary_nyc_portfolio_de.png", sep=""), res = 200, width = 1400, height = 1000)
plot(c(data$time[1], data$time[nrow(data)]), c(min(data[,4:5]), max(data[,4:5])+1000), cex=0, 
     xlab='Year', ylab='Investment threshold (mm)')
for (i in 1:2){
  lines(data[,1], data[,3+i], col=colour[i], lty=style[i])
}
legend('topright',legend=c("Barrier and dike","Flood-proofing"),
       lty=style,col=colour, bty = "n") 
#title(main="Frequency mitigation", font.main=1)
#dev.off()
#----------------------

#=== Make Latex table ====
head0 <- paste("\\hline &", "\\multicolumn{5}{c}{Flood-proofing then dike}", "&&",
               "\\multicolumn{5}{c}{Dike then flood-proofing}", "\\\\","\n",
               "\\cmidrule{2-6} \\cmidrule{8-12}" ,sep="")

head1 <- paste("&", "Flood-proofing", "&&","Barrier and Dike", '&&', "Total",'&&', "Barrier and Dike", "&&","Flood-proofing", '&&', "Total",
               "\\\\","\n","\\cmidrule{2-2} \\cmidrule{4-4} \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10} \\cmidrule{12-12}", sep="")
addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- -1
addtorow$pos[[2]] <- 0
addtorow$command <- c(head0,head1)

library(xtable)
x.width <- xtable(result,digit = 2, label=paste('analysis_nyc_portfolio_de',sep=''),
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
      file = paste(tab, "analysis_nyc_portfolio_de", ".tex", sep=""), scalebox = 0.65)








