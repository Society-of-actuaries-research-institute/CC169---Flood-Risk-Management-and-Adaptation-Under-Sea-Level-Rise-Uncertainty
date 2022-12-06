r=0.04;  
mu0=6 - 0.15
sigma= 25;gamma0 = 0.01; 
alpha0 = 1642; scale = 131; xi = 0.2747
delta = 0.03;
u = 2506;  a=0; b=0.0393

ustar =2800


# proofing is indexed 1, dike is indexed 2.
I2 <- 13+0.118/r
k2=1000
a2=0; b2=0.0393;
kappa2 =1;

nt <- as.integer(5000); ns <- as.integer(5000)
Smin = -2000; Smax = 10000; T=2000
To = 100; nto=as.integer(100000); 
dt = T/nt
ds = (Smax-Smin)/ns