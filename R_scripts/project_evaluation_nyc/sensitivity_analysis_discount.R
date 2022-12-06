
rs=seq(0.02,0.06,0.01); 
I1s <- rep(I1,length(rs)); I2s <- 13+0.118/rs
#=== Barrier and Dike ===
Is <- I2s;k=k2;
kappa =kappa2; 
mu <- mu0
gamma <- gamma0
result <- data.frame(matrix(NA, 3, length(rs)*2))
result[,1] <- c('NPV','Total value with optionality','Difference')
for (i in 1:length(rs)){
  r=rs[i]; I=Is[i]
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
head0 <- paste("\\hline &", paste("r=",100*rs[-length(rs)],"\\%", "&&",collapse='',sep=''),
               paste("r=",100*rs[length(rs)], "\\%",sep=''),
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
      file = paste(tab, "sensitivity_dike_nyc_discount", ".tex", sep=""), scalebox = 0.90)



#== Wet proofing project ===
Is <- I1s; k=k1
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
                  u,a0,b0,u,u+k,a1,b1,u, delta,kappa, scale,xi,f=matrix(0,nt,ns-1))
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
  out <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,I, f=0, B = matrix(0,nto+1,2))
  dyn.unload(paste(ftpath, "binop.dll", sep=""))
  out$f
  result[2,2*i] <- round(out$f,3)
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
x.width <- xtable(result,digit = 2, label=paste('sensitivity_elevation_nyc_discount',sep=''),
                  caption=paste("Investment analysis for wet flood-proofing project using 
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
      file = paste(tab, "sensitivity_elevation_nyc_discount", ".tex", sep=""), scalebox = 0.90)


#== Sequencing: proofing-dike and dike-proofing
mu <- mu0
gamma <- gamma0

#=== Investment Analysis ===
result0 <- data.frame(matrix(NA, 3, length(rs)*4))
result0[,1] <- c('NPV','Total value with optionality','Difference')
for(ii in 1:length(rs)){
  result <- data.frame(matrix(NA, 3, 12))
  result[,1] <- c('NPV','Total value with optionality','Difference')
  Is <- c(I1s[ii], I2s[ii], I2s[ii], I1s[ii])
  r=rs[ii]
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
    out <- .Fortran("crank", mu, sigma, r, gamma, nt, ns, Smin, Smax,T,
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
    out1 <- .Fortran("binop", mu, sigma,gamma, alpha0, r, nto, T, N, V,Is[i], f=0, B = matrix(0,nto+1,2))
    dyn.unload(paste(ftpath, "binop.dll", sep=""))
    out1$f
    if(i<3){
      result[2,2*i] <- round(out1$f,3) #option value
      result[3,2*i] <- result[2,2*i] - result[1,2*i] # value added
    }else{
      result[2,2+2*i] <- round(out1$f,3)
      result[3,2+2*i] <- result[2,2+2*i] - result[1,2+2*i]
    }
  }
  result[,6]=unlist(mapply(function(i) result[i,2]+result[i,4], 1:3,SIMPLIFY=FALSE))
  result[,12]=unlist(mapply(function(i) result[i,8]+result[i,10], 1:3,SIMPLIFY=FALSE))
  result0[,2+2*(ii-1)]=result[,6]
  result0[,2+2*length(rs)+2*(ii-1)]=result[,12]
}


#=== Make Latex table ====
head0 <- paste("\\hline &", paste("\\multicolumn{", length(rs)*2-1, "}{c}{Flood-proofing then dike}",sep=''),
               "&&",paste("\\multicolumn{", length(rs)*2-1, "}{c}{Dike then flood-proofing}",sep=''),
               "\\\\","\n",
               paste(
                 paste("\\cmidrule{2-",length(rs)*2, "}",sep=''),
                 paste("\\cmidrule{",length(rs)*2+2 ,"-",length(rs)*4, "}",sep='')
                 ,collapse=''),
               sep='')

idl = 2*seq(1,2*length(rs),1)
head1 <- paste("&", paste("r=",100*rs,"\\%", "&&",collapse='',sep=''),
               paste("r=",100*rs[-length(rs)],"\\%", "&&",collapse='',sep=""),
               paste("r=",100*rs[length(rs)], "\\%",sep=''),
               "\\\\","\n",        
               paste(unlist(mapply(function(i) paste("\\cmidrule{", idl[i],"-",idl[i], "}",sep=''), 1:length(idl), SIMPLIFY=FALSE)),collapse=''),
               sep="")

addtorow <- list()
addtorow$pos <- list()
addtorow$pos[[1]] <- -1
addtorow$pos[[2]] <- 0
addtorow$command <- c(head0,head1)

library(xtable)
x.width <- xtable(result0,digit = 2, label=paste('sensitivity_portfolio_de_nyc_discount',sep=''),
                  caption=paste("Investment sequences under different discount rates. We report (in billion euro) the NPV, 
                  the total value with optionality and the difference between the two values.",sep=""),
                  align =paste('ll',paste(rep('r',length(rs)*4-1),collapse=''),sep=''))
align(x.width) <-paste('ll',paste(rep('r',length(rs)*4-1),collapse=''),sep='')
comment <- list()
comment$pos <- list()
comment$pos[[1]] <- c(nrow(result0))

print(x.width,   hline.after = c(nrow(result0)),add.to.row = addtorow,
      sanitize.text=function(x){x}, caption.placement = 'top',
      floating = FALSE, include.rownames=FALSE,include.colnames=FALSE)

print(x.width,hline.after = c(nrow(result0)),add.to.row = addtorow, 
      sanitize.text=function(x){x}, caption.placement = 'top',table.placement='H',
      floating = TRUE,include.rownames=FALSE, include.colnames=FALSE,
      file = paste(tab, "sensitivity_portfolio_de_nyc_discount", ".tex", sep=""), scalebox = 0.65)









































