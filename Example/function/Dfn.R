
#Expected annual damage
# This function calculate the flood damage when flood threshold is x; and damage curve is at,bt,u, see Figure 1 in paper.
Dfn <- function(x,a,b,u,alpha,s,xi){
  c=a*u*u -b*u
  hx <- pevd(x, loc = alpha, scale = s, shape = xi, type = "GEV")
  if(xi>0 & (x>alpha -s/xi)){
    out <- a*((s/xi)^2)*(
      pgamma((1+xi*(x-alpha)/s)^(-1/xi), 1-2*xi)*gamma(1-2*xi)  # see?pgamma
    )+
      2*(s/xi)*(a*alpha-a*s/xi + b/2)*(
        pgamma((1+xi*(x-alpha)/s)^(-1/xi), 1-xi)*gamma(1-xi)
      )+
      (a*(alpha-s/xi)^2+b*(alpha-s/xi)+c)*(1-hx)  
  }
  if(xi>0 & (x<alpha -s/xi)){
    out <- a*(alpha-s/xi)*(alpha-s/xi) +
      (s*a/xi)*((s/xi)*gamma(1-2*xi)-2*((s/xi)-alpha)*gamma(1-xi))+
      b*alpha + (s*b/xi)*(gamma(1-xi)-1)+c
  }
  if(xi<0 & (x<alpha -s/xi)){
    out <- a*((s/xi)^2)*(
      pgamma((1+xi*(x-alpha)/s)^(-1/xi), 1-2*xi)*gamma(1-2*xi)  # see?pgamma
    )+
      2*(s/xi)*(a*alpha-a*s/xi + b/2)*(
        pgamma((1+xi*(x-alpha)/s)^(-1/xi), 1-xi)*gamma(1-xi)
      )+
      (a*(alpha-s/xi)^2+b*(alpha-s/xi)+c)*(1-hx)  
  }
  if(xi<0 & (x>alpha -s/xi))  {
    out <- 0
  }
  return(out)
}
# Dfn(u,u, alpha,s,xi, at,bt,measure=2)  # Annual expected flood damage without adaptation, see Figure 1 in paper.
# # Behaviour of damage function as alpha changes
# alphat <- seq(3000, 20000, 10)
# y <- mapply(function(i) Dfn(u,u, alpha=alphat[i],s,xi, at,bt,measure=2), 1:length(alphat), SIMPLIFY = FALSE)
# plot(alphat, y)
#----------