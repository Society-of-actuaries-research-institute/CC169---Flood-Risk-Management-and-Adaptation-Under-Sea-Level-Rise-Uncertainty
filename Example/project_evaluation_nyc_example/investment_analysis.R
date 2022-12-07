

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
      file = paste(tab, "analysis_nyc_dike", ".tex", sep=""), scalebox = 0.80)


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








