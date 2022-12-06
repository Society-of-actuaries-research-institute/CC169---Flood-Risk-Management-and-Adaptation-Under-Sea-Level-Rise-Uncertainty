r = 12*( (1+0.04)^(1/12) -1);
#r = 0.04
## Mean See Level Brownian Motion
mu0 = (6.5 - (- 0.79) )/12
#mu0 = 8.7 - (- 0.79)
sigma = 25/sqrt(12) 
#sigma = 25
## Loss Exposure Growth parameter
gamma0 =12*( (1+0.0131)^(1/12)-1 ) 
#gamma0 = 0.0131
## Safety Loading
delta = 0.03
## Starting Values of Max Tide Residual
# = 402.188 1180  1000
alpha0 =  402.188

scale = 254.56 #254.56#183.32  #254.56
#scale = 183.32
# xi = -0.05, 1.5, -0.0554
xi = -0.05

## yearly parameters
#r = 0.04
#gamma0 = 0.0131
#mu0=8.7 -(- 0.79)
#sigma = 25

## Current Threshold (Check with CHI)
u = 330.92
ustar = 600
## Loss Curve Parameters
a = 0.0000007453
b = 0.001562

# Dike is indexed 1, and then second dike is indexed 2.
k1= 1000  #750
k2= 1000 #500

I1 = (160000 +  1371 * (k1^2))/1000000000
#I1 = (1371 * (k1^2))/1000000000
#I1 = 0.973
I2 = (160000 +  1371 * (k2^2) + 1371 * (2*k1*k2 ) )/1000000000
#I2 = 2.919

# INvestment Cost of building the k1+k2 heihgt dike in one go ( saving on fixed cost)

#I12 = (1600000 +  1371 * (k1^2) + 1371 * (2*k1^2 ) + 1371 * (k1^2)) /1000000000 


I1and2 = (160000 +  1371 *(k1+k2)^2)/1000000000

#(1600000 +  1371 * (k1 + k2)^2) /1000000000 

a1=0.0000007453; b1=0.001562;
a2=0.0000007453; b2=0.001562;

kappa1 = 1
kappa2 = 1


nt <- as.integer(5000); ns <- as.integer(5000)
Smin = -2000; Smax = 10000; T=2000

To = 100
nto=as.integer(100000) 
dt = T/nt
ds = (Smax-Smin)/ns