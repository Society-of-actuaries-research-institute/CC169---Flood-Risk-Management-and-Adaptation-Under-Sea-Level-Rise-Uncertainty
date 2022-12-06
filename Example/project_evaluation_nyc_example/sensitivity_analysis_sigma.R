
sigmas = c(7 , 15, 25 , 30, 45 )

I1s <- rep(I1,length(sigmas))
I2s <- rep(13+0.118/r,length(sigmas))
#=== Barrier and Dike ===
Is <- I2s;k=k2;
kappa =kappa2; 
mu <- mu0
gamma <- gamma0
result <- data.frame(matrix(NA, 3, length(sigmas)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')

for (i in 1:length(sigmas)){
  sigma = sigmas[i]
  I = Is[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,a,b,u,u+k,a,b,u, delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,3)
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  out1 <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out1$f
  result[2,2*i] <- round(out1$f,3) #option value
  result[3,2*i] <- result[2,2*i] -result[1,2*i] # value added by real options method
}


#=== Make Latex table ====
head0 <- paste("\\hline &", paste("r=",sigmas[-length(sigmas)],"\\%", "&&",collapse='',sep=''),
               paste("r=",sigmas[length(sigmas)], "\\%",sep=''),
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} 
               \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10}" ,sep="")

addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 2, label=paste('sensitivity_dike_nyc_discount',sep=''),
                  caption=paste("Investment analysis for Barrier and Dike Project using 
                  NPV rule and real options methods under different discount rates. We report (in billion USD) the NPV, 
                  the total value with optionality and the difference between the two values. 
                  ",sep=""),
                  align =paste('ll',paste(rep('r',length(sigmas)*2-1),collapse=''),sep=''))
align(x.width) <-paste('ll',paste(rep('r',length(sigmas)*2-1),collapse=''),sep='')
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "sensitivity_dike_nyc_sigmas", ".tex", sep=""), scalebox = 0.90)





























