library(ggplot2)

############# FUNCTION #####################################################
LocationLevel_simulation <- function (N, T, mu.x, sigma.x, x0, dt) {
  
  n=N
  N <- ifelse(N%%2 == 0, N, N + 1) # make N an even integer
  nsteps <- length(seq(0, T, dt)) - 1
  
  shock <- matrix(rnorm((nsteps) * N/2,mean=0, sd = 1),ncol = N/2)
  
  output <- matrix(NA, nrow = nsteps + 1, ncol = N)
  output[1, ] <- x0
  
  time_length = seq(1, nsteps,1)*dt
  
  output[2:(nsteps+1), seq(1, N, 2)] = 
    mapply(function(j) x0 + mu.x*time_length + shock[,j]*sigma.x* sqrt(time_length), 1:(N/2), SIMPLIFY = TRUE)
  output[2:(nsteps+1),seq(2, N, 2)] = 
    mapply(function(j) x0 + mu.x*time_length - shock[,j]*sigma.x* sqrt(time_length), 1:(N/2), SIMPLIFY = TRUE)
  return(output[, 1:n])
  
  
  
  
  return(output[, 1:n])
}




########################## LOW SEA LEVEL RISE#########################
mu0 = 1.8
g = ( (1+gamma0)^(1/12)-1 ) 

################ LOSS DISTRIBUTION AT 100 YEARS with NO ADAPTATIOn ############

dt = T/(12*T)

set.seed(100)
alpha_simulated = LocationLevel_simulation(N = 2000,
                                           T = 100,
                                           mu.x = mu0/12,
                                           sigma.x = sigma/sqrt(12),
                                           x0 = alpha0,
                                           dt = dt)

alpha_simulated_no_policy = alpha_simulated[(12*3):(12*100),]

no_policy = matrix(0,
                   nrow = nrow(alpha_simulated_no_policy),
                   ncol = ncol(alpha_simulated_no_policy))
no_policy_top_cover = matrix(0,
                             nrow = nrow(alpha_simulated_no_policy),
                             ncol = ncol(alpha_simulated_no_policy))
measure = 2
for(i in 1: nrow(alpha_simulated_no_policy)){
  for( j in 1: ncol(alpha_simulated_no_policy)){
    
    no_policy[i,j] = (1+delta) * exp(12* g) * Dfn(x = u,
                                                  a = a,
                                                  b = b,
                                                  u = u,
                                                  alpha  = alpha_simulated_no_policy[i,j],
                                                  s = scale,
                                                  xi = xi)
    
    no_policy_top_cover[i,j] = no_policy[i,j] - (1+delta) * exp(12* g) * Dfn(x = ustar,
                                                                             a = a,
                                                                             b = b,
                                                                             u = u,
                                                                             alpha  = alpha_simulated_no_policy[i,j],
                                                                             s = scale,
                                                                             xi = xi)
    
    
  }
}



####################### DYKE ################################

# Stage 1 Boundary
boundary = read.csv(paste(opath, "threshold_cph_dyke_stage1.csv", sep=""), sep=",", header=TRUE)
# monthly comparison

idx_monthly =round( seq(3,100, by = 1/12)*1000)

boundary = boundary[idx_monthly,]

boundary_check1  = boundary$c4

# Stage 2 Boundary
boundary = read.csv(paste(opath, "threshold_cph_stage2.csv", sep=""), sep=",", header=TRUE)
# monthly comparison

idx_monthly =round( seq(3,100, by = 1/12)*1000)

boundary = boundary[idx_monthly,]

boundary_check2  = boundary$c4


alpha_simulated_dyke = alpha_simulated[(12*3):(12*100),]

dyke = matrix(0,
              ncol = ncol(alpha_simulated_dyke),
              nrow = nrow(alpha_simulated_dyke))

top_cover = matrix(0,
                   ncol = ncol(alpha_simulated_dyke),
                   nrow = nrow(alpha_simulated_dyke))


stopping_time1 = matrix(0,
                        nrow = ncol(alpha_simulated_dyke))
stopping_time2 =  matrix(0,
                         nrow = ncol(alpha_simulated_dyke))


results = data.frame(matrix(0,
                            nrow = 2,
                            ncol = 5))

colnames(results) = c("0%",
                      "25%",
                      "50%",
                      "75%",
                      "100%")

for ( j in 1:ncol(alpha_simulated_dyke)){
  check1 = 0
  check2 = 0
  for ( i in 1:nrow(alpha_simulated_dyke)){
    
    if ((alpha_simulated_dyke[i,j]<=boundary_check1[i]) && (check1 == 0) ){
      
      dyke[i,j] = (1+delta) * exp(12* g) * Dfn(x = u,
                                               a = a,
                                               b = b,
                                               u = u,
                                               alpha  = alpha_simulated_dyke[i,j],
                                               s = scale,
                                               xi = xi)
      
      top_cover[i,j] =  dyke[i,j] - (1+delta) * exp(12* g) * Dfn(x = ustar,
                                                                 a = a,
                                                                 b = b,
                                                                 u = u,
                                                                 alpha  = alpha_simulated_dyke[i,j],
                                                                 s = scale,
                                                                 xi = xi)
      
      
    } else if ( (alpha_simulated_dyke[i,j]<=boundary_check2[i]) && (check2 == 0) ) {
      
      dyke[i,j] = (1+delta) * exp(12* g) * Dfn(x = u +k1,
                                               a = a,
                                               b = b,
                                               u = u,
                                               alpha  = alpha_simulated_dyke[i,j],
                                               s = scale,
                                               xi = xi)
      
      top_cover[i,j] =  dyke[i,j] - (1+delta) * exp(12* g) * Dfn(x = ustar + k1,
                                                                 a = a,
                                                                 b = b,
                                                                 u = u,
                                                                 alpha  = alpha_simulated_dyke[i,j],
                                                                 s = scale,
                                                                 xi = xi)
      
      if (check1 == 0){
        stopping_time1[j] = i
        check1 = 1
      }
      
      
    } else {
      dyke[i,j] = (1+delta) * exp(12* g) *  Dfn(x = u+k1+k2 ,
                                                a = a,
                                                b = b,
                                                u = u+k1,
                                                alpha  = alpha_simulated_dyke[i,j],
                                                s = scale,
                                                xi = xi)
      
      top_cover[i,j] =  dyke[i,j]  - (1+delta) * exp(12* g) *  Dfn(x = ustar+k1 +k2,
                                                                   a = a,
                                                                   b = b,
                                                                   u = u+k1,
                                                                   alpha  = alpha_simulated_dyke[i,j],
                                                                   s = scale,
                                                                   xi = xi)
      
      if (check2 == 0){
        stopping_time2[j] = i
        check2 = 1
      }
    }
    
    
  }
  
}

LM_nopolicy = function(alpha, a, b,u){
  c=a*u*u -b*u
  L = a * alpha^2 + b*alpha + c
  return(L)
}

LM_dyke = function(alpha, a, b,u,k){
  u = u + k
  c=a*u*u -b*u
  L =max( a * (alpha)^2 + b*(alpha) + c,0)
  return(L)
}

LossesNoPolicy = unlist(mapply(function(i) LM_nopolicy(alpha = alpha_simulated_dyke[1165,i],
                                                       a = a,
                                                       b= b,
                                                       u =u),
                               1:length(alpha_simulated_dyke[1165,]),
                               SIMPLIFY = FALSE))

LossesPolicy = unlist(mapply(function(i) LM_dyke(alpha = alpha_simulated_dyke[1165,i],
                                                 a = a,
                                                 b = b,
                                                 u = u,
                                                 k = k2),
                             1:length(alpha_simulated_dyke[1165,]),
                             SIMPLIFY = FALSE))

df = data.frame(Losses = c(LossesNoPolicy,
                           LossesPolicy),
                Adaptation = c(rep("No Adaptation",
                                   length(LossesNoPolicy)),
                               rep("Dyke",
                                   length(LossesPolicy))),
                MeanNoPolicy = mean(LossesNoPolicy),
                MeanPolicy = mean(LossesPolicy))


f = ggplot(df, aes( x = Losses,
                    color = Adaptation,
                    fill = Adaptation)) 

f = f + geom_histogram(binwidth=0.0001,
                       alpha=0.5,
                       position="identity")


f = f + geom_vline(data=df, aes(xintercept = MeanNoPolicy), colour="blue")

f = f + geom_vline(data=df, aes(xintercept = MeanPolicy), colour="red")


f = f + theme_bw()

f = f + xlab("Billion Dollars") + ylab("Loss Disitrubtion")

#f

ggsave(filename = "Loss_Distribution_dyke_CPH_lowmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)


## Premiums at 100 Years

idx_nopolicy= 12 * c(10,
                     20,
                     30,
                     40,
                     50,
                     60,
                     70,
                     80,
                     90,
                     100)+ 1


idx_nopolicy = seq(7*12,1165,
                   by =12*10 )+1


idx_policy = seq(7*12,1165,
                 by =12*10 )+1





df = data.frame(value = c( no_policy[idx_nopolicy[1],],
                           dyke[idx_policy[1],],
                           no_policy[idx_nopolicy[2],],
                           dyke[idx_policy[2],],
                           no_policy[idx_nopolicy[3],],
                           dyke[idx_policy[3],],
                           no_policy[idx_nopolicy[4],],
                           dyke[idx_policy[4],],
                           no_policy[idx_nopolicy[5],],
                           dyke[idx_policy[5],],
                           no_policy[idx_nopolicy[6],],
                           dyke[idx_policy[6],],
                           no_policy[idx_nopolicy[7],],
                           dyke[idx_policy[7],],
                           no_policy[idx_nopolicy[8],],
                           dyke[idx_policy[8],],
                           no_policy[idx_nopolicy[9],],
                           dyke[idx_policy[9],],
                           no_policy[idx_nopolicy[10],],
                           dyke[idx_policy[10],]),
                name = factor(c( rep(10, length(no_policy[idx_nopolicy[1],]) + length(dyke[idx_policy[1],]) ),
                                 rep(20, length(no_policy[idx_nopolicy[2],]) + length(dyke[idx_policy[2],]) ),
                                 rep(30, length(no_policy[idx_nopolicy[3],]) + length(dyke[idx_policy[3],]) ),
                                 rep(40, length(no_policy[idx_nopolicy[4],]) + length(dyke[idx_policy[4],]) ),
                                 rep(50, length(no_policy[idx_nopolicy[5],]) + length(dyke[idx_policy[5],]) ),
                                 rep(60, length(no_policy[idx_nopolicy[6],]) + length(dyke[idx_policy[6],]) ),
                                 rep(70, length(no_policy[idx_nopolicy[7],]) + length(dyke[idx_policy[7],]) ),
                                 rep(80, length(no_policy[idx_nopolicy[8],]) + length(dyke[idx_policy[8],]) ),
                                 rep(90, length(no_policy[idx_nopolicy[9],]) + length(dyke[idx_policy[9],]) ),
                                 rep(100, length(no_policy[idx_nopolicy[10],]) + length(dyke[idx_policy[10],]) ))),
                Adaptation =  rep( c(rep("No Adaptation", length(no_policy[idx_nopolicy[1],])),
                                     rep("Dyke", length(dyke[idx_nopolicy[1],]))),10 )
                
)

df$Adaptation = factor(df$Adaptation,
                       levels = c("No Adaptation",
                                  "Dyke"))





group.colors <- c("No Adaptation" = "#FF3333",
                  "Dyke" = "#3361FF")

f = ggplot(df, aes( x = name,
                    y = value,
                    fill  = Adaptation))

f  = f + geom_boxplot()


f = f + theme_bw()

f = f + xlab("Years") + ylab("Premium Distribution ($Billion)")

f = f + scale_fill_manual(values=group.colors)

f = f + theme(axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              axis.title.x = element_text(size = 17),
              axis.title.y = element_text(size = 17),
              legend.text = element_text(size= 17),
              legend.title = element_text(size= 17))
f




ggsave(filename = "Premium_Distribution_dyke_time_CPH_lowmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)



### TOP COVER LIMIT

df = data.frame(value = c( no_policy_top_cover[idx_nopolicy[1],],
                           top_cover[idx_policy[1],],
                           no_policy_top_cover[idx_nopolicy[2],],
                           top_cover[idx_policy[2],],
                           no_policy_top_cover[idx_nopolicy[3],],
                           top_cover[idx_policy[3],],
                           no_policy_top_cover[idx_nopolicy[4],],
                           top_cover[idx_policy[4],],
                           no_policy_top_cover[idx_nopolicy[5],],
                           top_cover[idx_policy[5],],
                           no_policy_top_cover[idx_nopolicy[6],],
                           top_cover[idx_policy[6],],
                           no_policy_top_cover[idx_nopolicy[7],],
                           top_cover[idx_policy[7],],
                           no_policy_top_cover[idx_nopolicy[8],],
                           top_cover[idx_policy[8],],
                           no_policy_top_cover[idx_nopolicy[9],],
                           top_cover[idx_policy[9],],
                           no_policy_top_cover[idx_nopolicy[10],],
                           top_cover[idx_policy[10],]),
                name = factor(c( rep(10, length(no_policy_top_cover[idx_nopolicy[1],]) + length(top_cover[idx_policy[1],]) ),
                                 rep(20, length(no_policy_top_cover[idx_nopolicy[2],]) + length(top_cover[idx_policy[2],]) ),
                                 rep(30, length(no_policy_top_cover[idx_nopolicy[3],]) + length(top_cover[idx_policy[3],]) ),
                                 rep(40, length(no_policy_top_cover[idx_nopolicy[4],]) + length(top_cover[idx_policy[4],]) ),
                                 rep(50, length(no_policy_top_cover[idx_nopolicy[5],]) + length(top_cover[idx_policy[5],]) ),
                                 rep(60, length(no_policy_top_cover[idx_nopolicy[6],]) + length(top_cover[idx_policy[6],]) ),
                                 rep(70, length(no_policy_top_cover[idx_nopolicy[7],]) + length(top_cover[idx_policy[7],]) ),
                                 rep(80, length(no_policy_top_cover[idx_nopolicy[8],]) + length(top_cover[idx_policy[8],]) ),
                                 rep(90, length(no_policy_top_cover[idx_nopolicy[9],]) + length(top_cover[idx_policy[9],]) ),
                                 rep(100, length(no_policy_top_cover[idx_nopolicy[10],]) + length(top_cover[idx_policy[10],]) ))),
                Adaptation =  rep( c(rep("No Adaptation", length(no_policy_top_cover[idx_nopolicy[1],])),
                                     rep("Dyke", length(top_cover[idx_nopolicy[1],]))),10 )
                
)

df$Adaptation = factor(df$Adaptation,
                       levels = c("No Adaptation",
                                  "Dyke"))


group.colors <- c("No Adaptation" = "#FF3333",
                  "Dyke" = "#3361FF")


f = ggplot(df, aes( x = name,
                    y = value,
                    fill  = Adaptation))

f  = f + geom_boxplot()


f = f + theme_bw()

f = f + xlab("Years") + ylab("Premium Distribution ($Billion)")

f = f + scale_fill_manual(values=group.colors)

f = f + theme(axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              axis.title.x = element_text(size = 17),
              axis.title.y = element_text(size = 17),
              legend.text = element_text(size= 17),
              legend.title = element_text(size= 17))
f



ggsave(filename = "Premium_Distribution_Top_cover_dyke_time_CPH_lowmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)













df = data.frame(value = c(  no_policy[1165,],
                            dyke[1165,]),
                Adaptation = c(rep("No Adaptation",
                                   length(no_policy[1165,])),
                               rep("Dyke",
                                   length(dyke[1165,]))),
                MeanNoPolicy = mean(no_policy[1165,]),
                MeanPolicy = mean(dyke[1165,]))

f = ggplot(df, aes( x = value,
                    color = Adaptation,
                    fill = Adaptation)) 

f = f + geom_histogram(binwidth=0.002,
                       alpha=0.5,
                       position="identity")


f = f + geom_vline(data=df, aes(xintercept = MeanNoPolicy), colour="blue")

f = f + geom_vline(data=df, aes(xintercept = MeanPolicy), colour="red")


f = f + theme_bw()

f = f + xlab("Billion Dollars") + ylab("Premium Disitrubtion")

f


ggsave(filename = "Premium_Distribution_Dyke_CPH_lowmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)


## Stopping time
df = data.frame(value = stopping_time1)

f = ggplot(df, aes(x = value))

f = f + geom_histogram()
f = f + theme_bw()

f = f + xlab("Time") + ylab("Stopping Time Disitrubtion (stage1")


f




ggsave(filename = "Stopping_time_dyk_CPH_stage1_lowmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)

df = data.frame(value = stopping_time2[stopping_time2>0])

f = ggplot(df, aes(x = value))

f = f + geom_histogram()
f = f + theme_bw()

f = f + xlab("Time") + ylab("Stopping Time Disitrubtion (stage2)")


f




ggsave(filename = "Stopping_time_dyk_CPH_stage2_lowmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)



## Print Table quantile

results[1,] = as.numeric(quantile(no_policy[1165,]))
results[2,] = as.numeric(quantile(dyke[1165,]))

write.table(results, "Quantile_Dyke_CPH_lowmslvl.csv", sep = "&")




########################## HIGH SEA LEVEL RISE#########################
mu0 = 19
g = ( (1+gamma0)^(1/12)-1 ) 

################ LOSS DISTRIBUTION AT 100 YEARS with NO ADAPTATIOn ############

dt = T/(12*T)

set.seed(100)
alpha_simulated = LocationLevel_simulation(N = 2000,
                                           T = 100,
                                           mu.x = mu0/12,
                                           sigma.x = sigma/sqrt(12),
                                           x0 = alpha0,
                                           dt = dt)

alpha_simulated_no_policy = alpha_simulated[(12*3):(12*100),]

no_policy = matrix(0,
                   nrow = nrow(alpha_simulated_no_policy),
                   ncol = ncol(alpha_simulated_no_policy))
no_policy_top_cover = matrix(0,
                             nrow = nrow(alpha_simulated_no_policy),
                             ncol = ncol(alpha_simulated_no_policy))
measure = 2
for(i in 1: nrow(alpha_simulated_no_policy)){
  for( j in 1: ncol(alpha_simulated_no_policy)){
    
    no_policy[i,j] = (1+delta) * exp(12* g) * Dfn(x = u,
                                                  a = a,
                                                  b = b,
                                                  u = u,
                                                  alpha  = alpha_simulated_no_policy[i,j],
                                                  s = scale,
                                                  xi = xi)
    
    no_policy_top_cover[i,j] = no_policy[i,j] - (1+delta) * exp(12* g) * Dfn(x = ustar,
                                                                             a = a,
                                                                             b = b,
                                                                             u = u,
                                                                             alpha  = alpha_simulated_no_policy[i,j],
                                                                             s = scale,
                                                                             xi = xi)
    
    
  }
}



####################### DYKE ################################

# Stage 1 Boundary
boundary = read.csv(paste(opath, "threshold_cph_dyke_stage1.csv", sep=""), sep=",", header=TRUE)
# monthly comparison

idx_monthly =round( seq(3,100, by = 1/12)*1000)

boundary = boundary[idx_monthly,]

boundary_check1  = boundary$c4

# Stage 2 Boundary
boundary = read.csv(paste(opath, "threshold_cph_stage2.csv", sep=""), sep=",", header=TRUE)
# monthly comparison

idx_monthly =round( seq(3,100, by = 1/12)*1000)

boundary = boundary[idx_monthly,]

boundary_check2  = boundary$c4


alpha_simulated_dyke = alpha_simulated[(12*3):(12*100),]

dyke = matrix(0,
              ncol = ncol(alpha_simulated_dyke),
              nrow = nrow(alpha_simulated_dyke))

top_cover = matrix(0,
                   ncol = ncol(alpha_simulated_dyke),
                   nrow = nrow(alpha_simulated_dyke))


stopping_time1 = matrix(0,
                        nrow = ncol(alpha_simulated_dyke))
stopping_time2 =  matrix(0,
                         nrow = ncol(alpha_simulated_dyke))


results = data.frame(matrix(0,
                            nrow = 2,
                            ncol = 5))

colnames(results) = c("0%",
                      "25%",
                      "50%",
                      "75%",
                      "100%")

for ( j in 1:ncol(alpha_simulated_dyke)){
  check1 = 0
  check2 = 0
  for ( i in 1:nrow(alpha_simulated_dyke)){
    
    if ((alpha_simulated_dyke[i,j]<=boundary_check1[i]) && (check1 == 0) ){
      
      dyke[i,j] = (1+delta) * exp(12* g) * Dfn(x = u,
                                               a = a,
                                               b = b,
                                               u = u,
                                               alpha  = alpha_simulated_dyke[i,j],
                                               s = scale,
                                               xi = xi)
      
      top_cover[i,j] =  dyke[i,j] - (1+delta) * exp(12* g) * Dfn(x = ustar,
                                                                 a = a,
                                                                 b = b,
                                                                 u = u,
                                                                 alpha  = alpha_simulated_dyke[i,j],
                                                                 s = scale,
                                                                 xi = xi)
      
      
    } else if ( (alpha_simulated_dyke[i,j]<=boundary_check2[i]) && (check2 == 0) ) {
      
      dyke[i,j] = (1+delta) * exp(12* g) * Dfn(x = u +k1,
                                               a = a,
                                               b = b,
                                               u = u,
                                               alpha  = alpha_simulated_dyke[i,j],
                                               s = scale,
                                               xi = xi)
      
      top_cover[i,j] =  dyke[i,j] - (1+delta) * exp(12* g) * Dfn(x = ustar + k1,
                                                                 a = a,
                                                                 b = b,
                                                                 u = u,
                                                                 alpha  = alpha_simulated_dyke[i,j],
                                                                 s = scale,
                                                                 xi = xi)
      
      if (check1 == 0){
        stopping_time1[j] = i
        check1 = 1
      }
      
      
    } else {
      dyke[i,j] = (1+delta) * exp(12* g) *  Dfn(x = u+k1+k2 ,
                                                a = a,
                                                b = b,
                                                u = u+k1,
                                                alpha  = alpha_simulated_dyke[i,j],
                                                s = scale,
                                                xi = xi)
      
      top_cover[i,j] =  dyke[i,j]  - (1+delta) * exp(12* g) *  Dfn(x = ustar+k1 +k2,
                                                                   a = a,
                                                                   b = b,
                                                                   u = u+k1,
                                                                   alpha  = alpha_simulated_dyke[i,j],
                                                                   s = scale,
                                                                   xi = xi)
      
      if (check2 == 0){
        stopping_time2[j] = i
        check2 = 1
      }
    }
    
    
  }
  
}

LM_nopolicy = function(alpha, a, b,u){
  c=a*u*u -b*u
  L = a * alpha^2 + b*alpha + c
  return(L)
}

LM_dyke = function(alpha, a, b,u,k){
  u = u + k
  c=a*u*u -b*u
  L =max( a * (alpha)^2 + b*(alpha) + c,0)
  return(L)
}

LossesNoPolicy = unlist(mapply(function(i) LM_nopolicy(alpha = alpha_simulated_dyke[1165,i],
                                                       a = a,
                                                       b= b,
                                                       u =u),
                               1:length(alpha_simulated_dyke[1165,]),
                               SIMPLIFY = FALSE))

LossesPolicy = unlist(mapply(function(i) LM_dyke(alpha = alpha_simulated_dyke[1165,i],
                                                 a = a,
                                                 b = b,
                                                 u = u,
                                                 k = k2),
                             1:length(alpha_simulated_dyke[1165,]),
                             SIMPLIFY = FALSE))

df = data.frame(Losses = c(LossesNoPolicy,
                           LossesPolicy),
                Adaptation = c(rep("No Adaptation",
                                   length(LossesNoPolicy)),
                               rep("Dyke",
                                   length(LossesPolicy))),
                MeanNoPolicy = mean(LossesNoPolicy),
                MeanPolicy = mean(LossesPolicy))


f = ggplot(df, aes( x = Losses,
                    color = Adaptation,
                    fill = Adaptation)) 

f = f + geom_histogram(binwidth=0.0001,
                       alpha=0.5,
                       position="identity")


f = f + geom_vline(data=df, aes(xintercept = MeanNoPolicy), colour="blue")

f = f + geom_vline(data=df, aes(xintercept = MeanPolicy), colour="red")


f = f + theme_bw()

f = f + xlab("Billion Dollars") + ylab("Loss Distribution")

#f

ggsave(filename = "Loss_Distribution_dyke_CPH_highmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)


## Premiums at 100 Years

idx_nopolicy= 12 * c(10,
                     20,
                     30,
                     40,
                     50,
                     60,
                     70,
                     80,
                     90,
                     100)+ 1


idx_nopolicy = seq(7*12,1165,
                   by =12*10 )+1


idx_policy = seq(7*12,1165,
                 by =12*10 )+1





df = data.frame(value = c( no_policy[idx_nopolicy[1],],
                           dyke[idx_policy[1],],
                           no_policy[idx_nopolicy[2],],
                           dyke[idx_policy[2],],
                           no_policy[idx_nopolicy[3],],
                           dyke[idx_policy[3],],
                           no_policy[idx_nopolicy[4],],
                           dyke[idx_policy[4],],
                           no_policy[idx_nopolicy[5],],
                           dyke[idx_policy[5],],
                           no_policy[idx_nopolicy[6],],
                           dyke[idx_policy[6],],
                           no_policy[idx_nopolicy[7],],
                           dyke[idx_policy[7],],
                           no_policy[idx_nopolicy[8],],
                           dyke[idx_policy[8],],
                           no_policy[idx_nopolicy[9],],
                           dyke[idx_policy[9],],
                           no_policy[idx_nopolicy[10],],
                           dyke[idx_policy[10],]),
                name = factor(c( rep(10, length(no_policy[idx_nopolicy[1],]) + length(dyke[idx_policy[1],]) ),
                                 rep(20, length(no_policy[idx_nopolicy[2],]) + length(dyke[idx_policy[2],]) ),
                                 rep(30, length(no_policy[idx_nopolicy[3],]) + length(dyke[idx_policy[3],]) ),
                                 rep(40, length(no_policy[idx_nopolicy[4],]) + length(dyke[idx_policy[4],]) ),
                                 rep(50, length(no_policy[idx_nopolicy[5],]) + length(dyke[idx_policy[5],]) ),
                                 rep(60, length(no_policy[idx_nopolicy[6],]) + length(dyke[idx_policy[6],]) ),
                                 rep(70, length(no_policy[idx_nopolicy[7],]) + length(dyke[idx_policy[7],]) ),
                                 rep(80, length(no_policy[idx_nopolicy[8],]) + length(dyke[idx_policy[8],]) ),
                                 rep(90, length(no_policy[idx_nopolicy[9],]) + length(dyke[idx_policy[9],]) ),
                                 rep(100, length(no_policy[idx_nopolicy[10],]) + length(dyke[idx_policy[10],]) ))),
                Adaptation =  rep( c(rep("No Adaptation", length(no_policy[idx_nopolicy[1],])),
                                     rep("Dyke", length(dyke[idx_nopolicy[1],]))),10 )
                
)

df$Adaptation = factor(df$Adaptation,
                       levels = c("No Adaptation",
                                  "Dyke"))





group.colors <- c("No Adaptation" = "#FF3333",
                  "Dyke" = "#3361FF")

f = ggplot(df, aes( x = name,
                    y = value,
                    fill  = Adaptation))

f  = f + geom_boxplot()


f = f + theme_bw()

f = f + xlab("Years") + ylab("Premium Distribution ($Billion)")

f = f + scale_fill_manual(values=group.colors)

f = f + theme(axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              axis.title.x = element_text(size = 17),
              axis.title.y = element_text(size = 17),
              legend.text = element_text(size= 17),
              legend.title = element_text(size= 17))
f




ggsave(filename = "Premium_Distribution_dyke_time_CPH_highmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)



### TOP COVER LIMIT

df = data.frame(value = c( no_policy_top_cover[idx_nopolicy[1],],
                           top_cover[idx_policy[1],],
                           no_policy_top_cover[idx_nopolicy[2],],
                           top_cover[idx_policy[2],],
                           no_policy_top_cover[idx_nopolicy[3],],
                           top_cover[idx_policy[3],],
                           no_policy_top_cover[idx_nopolicy[4],],
                           top_cover[idx_policy[4],],
                           no_policy_top_cover[idx_nopolicy[5],],
                           top_cover[idx_policy[5],],
                           no_policy_top_cover[idx_nopolicy[6],],
                           top_cover[idx_policy[6],],
                           no_policy_top_cover[idx_nopolicy[7],],
                           top_cover[idx_policy[7],],
                           no_policy_top_cover[idx_nopolicy[8],],
                           top_cover[idx_policy[8],],
                           no_policy_top_cover[idx_nopolicy[9],],
                           top_cover[idx_policy[9],],
                           no_policy_top_cover[idx_nopolicy[10],],
                           top_cover[idx_policy[10],]),
                name = factor(c( rep(10, length(no_policy_top_cover[idx_nopolicy[1],]) + length(top_cover[idx_policy[1],]) ),
                                 rep(20, length(no_policy_top_cover[idx_nopolicy[2],]) + length(top_cover[idx_policy[2],]) ),
                                 rep(30, length(no_policy_top_cover[idx_nopolicy[3],]) + length(top_cover[idx_policy[3],]) ),
                                 rep(40, length(no_policy_top_cover[idx_nopolicy[4],]) + length(top_cover[idx_policy[4],]) ),
                                 rep(50, length(no_policy_top_cover[idx_nopolicy[5],]) + length(top_cover[idx_policy[5],]) ),
                                 rep(60, length(no_policy_top_cover[idx_nopolicy[6],]) + length(top_cover[idx_policy[6],]) ),
                                 rep(70, length(no_policy_top_cover[idx_nopolicy[7],]) + length(top_cover[idx_policy[7],]) ),
                                 rep(80, length(no_policy_top_cover[idx_nopolicy[8],]) + length(top_cover[idx_policy[8],]) ),
                                 rep(90, length(no_policy_top_cover[idx_nopolicy[9],]) + length(top_cover[idx_policy[9],]) ),
                                 rep(100, length(no_policy_top_cover[idx_nopolicy[10],]) + length(top_cover[idx_policy[10],]) ))),
                Adaptation =  rep( c(rep("No Adaptation", length(no_policy_top_cover[idx_nopolicy[1],])),
                                     rep("Dyke", length(top_cover[idx_nopolicy[1],]))),10 )
                
)

df$Adaptation = factor(df$Adaptation,
                       levels = c("No Adaptation",
                                  "Dyke"))


group.colors <- c("No Adaptation" = "#FF3333",
                  "Dyke" = "#3361FF")


f = ggplot(df, aes( x = name,
                    y = value,
                    fill  = Adaptation))

f  = f + geom_boxplot()


f = f + theme_bw()

f = f + xlab("Years") + ylab("Premium Distribution ($Billion)")

f = f + scale_fill_manual(values=group.colors)

f = f + theme(axis.text.x = element_text(size = 17),
              axis.text.y = element_text(size = 17),
              axis.title.x = element_text(size = 17),
              axis.title.y = element_text(size = 17),
              legend.text = element_text(size= 17),
              legend.title = element_text(size= 17))
f



ggsave(filename = "Premium_Distribution_Top_cover_dyke_time_CPH_highmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)













df = data.frame(value = c(  no_policy[1165,],
                            dyke[1165,]),
                Adaptation = c(rep("No Adaptation",
                                   length(no_policy[1165,])),
                               rep("Dyke",
                                   length(dyke[1165,]))),
                MeanNoPolicy = mean(no_policy[1165,]),
                MeanPolicy = mean(dyke[1165,]))

f = ggplot(df, aes( x = value,
                    color = Adaptation,
                    fill = Adaptation)) 

f = f + geom_histogram(binwidth=0.002,
                       alpha=0.5,
                       position="identity")


f = f + geom_vline(data=df, aes(xintercept = MeanNoPolicy), colour="blue")

f = f + geom_vline(data=df, aes(xintercept = MeanPolicy), colour="red")


f = f + theme_bw()

f = f + xlab("Billion Dollars") + ylab("Premium Disitrubtion")

f


ggsave(filename = "Premium_Distribution_Dyke_CPH_highmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)


## Stopping time
df = data.frame(value = stopping_time1)

f = ggplot(df, aes(x = value))

f = f + geom_histogram()
f = f + theme_bw()

f = f + xlab("Time") + ylab("Stopping Time Disitrubtion (stage1")


f




ggsave(filename = "Stopping_time_dyk_CPH_stage1_highmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)

df = data.frame(value = stopping_time2[stopping_time2>0])

f = ggplot(df, aes(x = value))

f = f + geom_histogram()
f = f + theme_bw()

f = f + xlab("Time") + ylab("Stopping Time Disitrubtion (stage2)")


f




ggsave(filename = "Stopping_time_dyk_CPH_stage2_highmslvl.eps",
       plot = f,
       device="eps",
       dpi = 600)



## Print Table quantile

results[1,] = as.numeric(quantile(no_policy[1165,]))
results[2,] = as.numeric(quantile(dyke[1165,]))

write.table(results, "Quantile_Dyke_CPH_highmslvl.csv", sep = "&")



################### THE END ################ 


