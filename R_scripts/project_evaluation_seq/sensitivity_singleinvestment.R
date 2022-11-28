
nt <- as.integer(10000); ns <- as.integer(10000)
Smin = -2000; Smax = 10000; T=2000
nto=as.integer(100000); 
dt = T/nt
ds = (Smax-Smin)/ns
r0=r
#=== Sensitivity analysis - Discount rate ===
rs=seq(0.02,0.06,0.01); 
mu <- mu0;gamma <- gamma0
kappa=1
result <- data.frame(matrix(NA, 3, length(rs)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')
for (i in 1:length(rs)){
  r=rs[i]; 
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,a,b,u,u+k,a,b,u+k, delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
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
head0 <- paste("\\hline &", paste("r=",100*rs[-length(rs)],"\\%", "&&",collapse='',sep=''),
               paste("r=",100*rs[length(rs)], "\\%",sep=''),
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} 
               \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10}" ,sep="")

addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 2, label=paste('sensitivity_singleinvestment_discount_qld',sep=''),
                  caption=paste("Sensitivity analysis for SEQ adaptation with respect to discount rate. We report (in billion AUD) the NPV, 
                  the total value with optionality and the difference between the two values. 
                  ",sep=""),
                  align =paste('ll',paste(rep('r',length(rs)*2-1),collapse=''),sep=''))
align(x.width) <-paste('ll',paste(rep('r',length(rs)*2-1),collapse=''),sep='')
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "sensitivity_singleinvestment_discount_qld", ".tex", sep=""), scalebox = 0.90)


#=== Sensitivity analysis - Sea Level Rise ===
r <- r0; gamma <- gamma0
kappa=1
mus <- c(0,0.5,1,1.5,2)*8.7 + 0.38
result <- data.frame(matrix(NA, 3, length(mus)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')
for (i in 1:length(mus)){
  mu = mus[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,a,b,u,u+k,a,b,u+k, delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
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
head0 <- paste("\\hline &", paste("$\\mu$=",(mus-0.38)[-length(mus)], "&&",collapse='',sep=''),
               paste("$\\mu$=",(mus-0.38)[length(mus)], sep=''),
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} 
               \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10}" ,sep="")

addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 2, label=paste('sensitivity_singleinvestment_slr_qld',sep=''),
                  caption=paste("Sensitivity analysis for SEQ adaptation with respect to the expected sea level rise. We report (in billion AUD) the NPV, 
                  the total value with optionality and the difference between the two values. 
                  ",sep=""),
                  align =paste('ll',paste(rep('r',length(mus)*2-1),collapse=''),sep=''))
align(x.width) <-paste('ll',paste(rep('r',length(mus)*2-1),collapse=''),sep='')
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result))

print(x.width,   hline.after = c(nrow(result)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "sensitivity_singleinvestment_slr_qld", ".tex", sep=""), scalebox = 0.90)




































