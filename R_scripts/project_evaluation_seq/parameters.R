r=0.04; 
mu0=8.7 + 0.38 #This is net amount: 3.6 for mu; 0.15 from theta. 
sigma=25;
gamma0 = 0.01
delta = 0.03; 

alpha0=1666;scale=53.23; xi=0.45;k=1000;
I1 =  0.0294*980000*43000/10^9  #cost is in billion AUD
# cost of adapting new buildings: 
g=0.01
N = 29400*(1+g)^(seq(1,100,1))
Nl = c(29400,N[-100])
AF = (1+r)^(-seq(1,100,1))
I2 = sum(0.002*10*0.01*150*1285*(N-Nl)*AF)/10^9  #cost is in billion AUD
I = I1 + I2

# Loss curve
u=1765
ustar = 2000

a=1.744*(10^(-8)); b=8.535*(10^(-5)); c = -0.205

# loss curve
M = seq(0, 5000,10)
L = a*M^2 + b*M + c
plot(M,L)

