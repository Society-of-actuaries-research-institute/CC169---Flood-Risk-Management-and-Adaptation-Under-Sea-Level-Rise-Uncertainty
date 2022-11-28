
r1=r
rs=c(0.02,0.03,r1,0.05,0.06) 


I1 = (160000 +  1371 * (k1^2))/1000000000
#I1 = (1371 * (k1^2))/1000000000
#I1 = 0.973
I2 = (160000 +  1371 * (k2^2) + 1371 * (2*k1*k2 ) )/1000000000
#I2 = 2.919

# INvestment Cost of building the k1+k2 heihgt dyke in one go ( saving on fixed cost)

#I12 = (1600000 +  1371 * (k1^2) + 1371 * (2*k1^2 ) + 1371 * (k1^2)) /1000000000 


I1and2 = (160000 +  1371 *(k1+k2)^2)/1000000000




I1s <- rep(I1,length(rs))
I2s <- rep(I2,length(rs)) 
I12s = rep(I1and2, length(rs))

#=== Barrier Stage 1 ===
Is <- I1s;k=k1;
kappa =kappa2; 
mu <- mu0
gamma <- gamma0
result <- data.frame(matrix(NA, 3, length(rs)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')
for (i in 1:length(rs)){
  r=rs[i]; I=Is[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,
                  a,b,u,
                  u+k1,
                  a,b,u,
                  delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,5)
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  out1 <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out1$f
  result[2,2*i] <- round(out1$f,5) #option value
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
x.width <- xtable(result,digit = 5, label=paste('sensitivity_dyke_cph_discount_stage1',sep=''),
                  caption=paste("Investment analysis for Dyke Project Stage 1 (k=750 mm) using 
                  NPV rule and real options methods under different discount rates. We report (in billion Euros) the NPV, 
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
      file = paste(tab, "sensitivity_dyke_cph_discount_stage1", ".tex", sep=""), scalebox = 0.90)



#==Barrier Stage 2 ===
Is <- I2s; k=k2
kappa =kappa1; 
mu <- mu0
gamma <- gamma0
a0=a1; b0=b1; 

result <- data.frame(matrix(NA, 3, length(rs)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')

#=== Investment Analysis ===
for (i in 1:length(rs)){
  r=rs[i]; I=Is[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u+k1,
                  a0,b0,u,
                  u+k1+k2,
                  a1,b1,u,
                  delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,5)
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  out <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out$f
  result[2,2*i] <- round(out$f,5)
  result[3,2*i] <- result[2,2*i] - max(result[1,2*i],0)
}
#=== Make Latex table ==
head0 <- paste("\\hline &", paste("r=",100*rs[-length(rs)],"\\%", "&&",collapse='',sep=''),
               paste("r=",100*rs[length(rs)], "\\%",sep=''),
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} 
               \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10}" ,sep="")

addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 5, label=paste('sensitivity_dyke_cph_discount_stage2',sep=''),
                  caption=paste("Investment analysis for Dyke project Stage 2 (k = 500mm) using 
                  NPV rule and real options methods under different discount rates. We report (in billion USD) the NPV, 
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
      file = paste(tab, "sensitivity_dyke_cph_discount_stage2", ".tex", sep=""), scalebox = 0.90)






#==Barrier One Shot ===
Is <- I12s; k=k2
kappa =kappa1; 
mu <- mu0
gamma <- gamma0
a0=a1; b0=b1; 

result <- data.frame(matrix(NA, 3, length(rs)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')

#=== Investment Analysis ===
for (i in 1:length(rs)){
  r=rs[i]; I=Is[i]
  dyn.load(paste(ftpath, "crank.dll", sep=""))
  out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
                  u,
                  a0,b0,u,
                  u+k1+k2,
                  a1,b1,u,
                  delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
  
  dyn.unload(paste(ftpath, "crank.dll", sep=""))
  j=floor((alpha0-Smin)/ds)
  Sj=Smin + j*ds
  dx=(alpha0-Sj)/ds
  out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I  # NPV
  result[1,2*i] <- round(out$f[1,j]+dx*(out$f[1,j+1]-out$f[1,j]) - I,5)
  x <- seq(Smin+ds, Smax-ds, ds)
  data <- data.frame(x, out$f[1,])
  V <- as.matrix(data)
  N <- as.integer(nrow(V))
  dyn.load(paste(ftpath, "binop.dll", sep=""))
  out <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out$f
  result[2,2*i] <- round(out$f,5)
  result[3,2*i] <- result[2,2*i] - max(result[1,2*i],0)
}
#=== Make Latex table ==
head0 <- paste("\\hline &", paste("r=",100*rs[-length(rs)],"\\%", "&&",collapse='',sep=''),
               paste("r=",100*rs[length(rs)], "\\%",sep=''),
               "\\\\","\n",        "\\cmidrule{2-2} \\cmidrule{4-4} 
               \\cmidrule{6-6} \\cmidrule{8-8} \\cmidrule{10-10}" ,sep="")

addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- 0
addtorow$command <- c(head0)

library(xtable)
x.width <- xtable(result,digit = 5, label=paste('sensitivity_dyke_cph_discount_one_shot',sep=''),
                  caption=paste("Investment analysis for Dyke project Stage 2 (k = 500mm) using 
                  NPV rule and real options methods under different discount rates. We report (in billion USD) the NPV, 
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
      file = paste(tab, "sensitivity_dyke_cph_discount_one_shot", ".tex", sep=""), scalebox = 0.90)






















