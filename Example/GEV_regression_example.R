require('extRemes')
require('gsl')
require(zoo)
require(lubridate)
require(dplyr)

############################ FUNCTIONS #######################################
LR_test_MaxTideRes = function(covariates, covariates.list.location, covariates.list.scale, 
                              use.phi = FALSE, type = "GEV",method = "GMLE" ){
	# This function produces the LR test for each coefficient in the GEV regression.
  
  

  formula.loc = as.formula(paste("~",  paste(covariates.list.location, collapse = "+")))
  formula.scale = as.formula(paste("~",  paste(covariates.list.scale, collapse = "+")))
  
  model1 = fevd(x = MaxTideResDeTrend,
                data = covariates,
                location.fun = formula.loc,
                scale.fun = formula.scale,
                shape = ~1,
                use.phi = use.phi,
                type = type,
                method = method
  )
  
  
  pvalue.location = data.frame(matrix(0,
                             nrow = 1,
                             ncol = length(covariates.list.location)))
  
  llik1 = model1$results$value
  
  for ( i in 1:length(covariates.list.location)){
    formula.loc = as.formula(paste("~",  paste(covariates.list.location[-i], collapse = "+")))
    formula.scale = as.formula(paste("~",  paste(covariates.list.scale, collapse = "+")))
    
    model0 = fevd(x = MaxTideResDeTrend,
                  data = covariates,
                  location.fun = formula.loc,
                  scale.fun = formula.scale,
                  shape = ~1,
                  use.phi = use.phi,
                  type = type,
                  method = method
    )
    
    llik0 = model0$results$value
    
    pvalue.location[1,i] = lr.test(model0,model1)$p.value
    
  }
  
  colnames(pvalue.location) = covariates.list.location
  rownames(pvalue.location) = "location"
  
  pvalue.scale = data.frame(matrix(0,
                                      nrow = 1,
                                      ncol = length(covariates.list.scale)))
  
  for ( i in 1:length(covariates.list.scale)){
    formula.loc = as.formula(paste("~",  paste(covariates.list.location, collapse = "+")))
    formula.scale = as.formula(paste("~",  paste(covariates.list.scale[-i], collapse = "+")))
    
    model0 = fevd(x = MaxTideResDeTrend,
                  data = covariates,
                  location.fun = formula.loc,
                  scale.fun = formula.scale,
                  shape = ~1,
                  use.phi = use.phi,
                  type = type,
                  method = method
    )
    
    llik0 = model0$results$value
    
    pvalue.scale[1,i] = lr.test(model0,model1)$p.value
    
  } 
  
  
  colnames(pvalue.scale) = covariates.list.scale
  rownames(pvalue.scale) ="scale"
  out = list(pvalue.location, pvalue.scale)
  names(out) = c("location", "scale")
  return(out)
  
}
################################ NYC ##########################################
setwd("./data/nyc")


data_hourly =  read.csv(file = paste("hourly_series_predict", ".csv", sep=""),
                        sep=",",
                        header=TRUE)

data_hourly$TideRes = data_hourly$wl-data_hourly$wlh



max_TideRes = data.frame(aggregate(data_hourly$TideRes,
                                    data.frame(as.yearmon(data_hourly$time)),
                                    FUN="max"))

names(max_TideRes) = c('Time','MaxTideRes')

max_WaterLvl = data.frame(aggregate(data_hourly$wl,
                                    data.frame(as.yearmon(data_hourly$time)),
                                    FUN="max"))

names(max_WaterLvl) = c('Time','MaxWaterLvl')

yearly_median = aggregate(data_hourly$TideRes,
                          data.frame(year(data_hourly$time)),
                          FUN = "median")
names(yearly_median) = c("Year", "Median")

yearly_mean = aggregate(data_hourly$TideRes,
                        data.frame(year(data_hourly$time)),
                        FUN = "mean")
names(yearly_mean) = c("Year", "Mean")


monthly_obs = aggregate(data_hourly$TideRes,
                         data.frame(as.yearmon(data_hourly$time)),
                         FUN = "length")

colnames(monthly_obs) = c("Year", "Observations")

monthly_obs$Observations = monthly_obs$Observations/max(monthly_obs$Observations)



data_monthly = merge(x = max_WaterLvl, 
                    y = max_TideRes,
                    by = "Time",
                    all.y = TRUE)

data_monthly$YearlyMedian = 0
data_monthly$YearlyMean = 0

for ( t in 1:length(data_monthly$Time)){
  
  data_monthly$YearlyMedian[t] = yearly_median$Median[yearly_median$Year == year(data_monthly$Time[t])]
  data_monthly$YearlyMean[t] = yearly_mean$Mean[yearly_mean$Year == year(data_monthly$Time[t])]
  
}

data_monthly$MaxDeTrend = data_monthly$MaxWaterLvl-data_monthly$YearlyMedian

data_monthly$Observation = monthly_obs$Observations

data_monthly = subset(data_monthly,data_monthly$Time>=1950)

data_monthly = subset(data_monthly, data_monthly$Time<2019)

data_monthly$MaxTideResDetrend = data_monthly$MaxTideRes - abs(data_monthly$YearlyMean)

data_monthly = subset(data_monthly, data_monthly$Observation>0.6)




### Covariates
# From Menéndez1 and Philip L. Woodworth (2010) three set up for Sea Level

covariates = data.frame(Time = data_monthly$Time,
                        LinearTrend = 1:length(data_monthly$Time)/12,
                        ExpTrend = exp((1:length(data_monthly$Time))/12),
                        AnnCos = cos( 2*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnSin = sin( 2*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnCos2 = cos( 4*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnSin2 = sin( 4*pi*c(1:length(data_monthly$Time))/12 ),
                        NodalCos = cos( 2*pi*c(1:length(data_monthly$Time))/(12*18.6)),
                        NodalSin = sin( 2*pi*c(1:length(data_monthly$Time))/(12*18.6)),
                        PerigeanCos =cos( 2*pi*c(1:length(data_monthly$Time))/(12*4.4)),
                        PerigeanSin = sin(2*pi*c(1:length(data_monthly$Time))/(12*4.4))
                        )






# Climate indexes


nao = read.csv("nao.csv", header = F)
nao = subset(nao, nao$V1<=2018)

nino = read.csv("Nino 34.csv")
nino = subset(nino, nino$X.1>=1950 & nino$X.1<2019)

time_range = format(seq(ymd('1950-01-01'),ymd('2018-12-31'),by='month'), "%Y-%m")


missing_months = setdiff(as.yearmon(time_range), as.yearmon(data_monthly$Time
) )


index_keep = is.na(match(as.yearmon(time_range), as.yearmon(missing_months)))


nao_values = unlist(nao[,2:13])

nao_values = nao_values[index_keep]

nino_values = unlist(nino[,3:14])
nino_values = nino_values[index_keep]

covariates$NAO = nao_values

covariates$NINO = nino_values

covariates$MaxWaterLvl = data_monthly$MaxWaterLvl
covariates$MaxTideRes = data_monthly$MaxTideRes
covariates$MaxDeTrend = data_monthly$MaxDeTrend
covariates$MaxTideResDeTrend = data_monthly$MaxTideResDetrend

covariates$LagNAO = lag(covariates$NAO)
covariates$LagNINO = lag(covariates$NINO)


covariates= covariates[-1,]

# Stationary model. The estimates are used as starting values
model_stationary = fevd(x =MaxTideResDeTrend,
              data = covariates,
              location.fun = ~1,
              scale.fun = ~1,
              shape = ~1,
              use.phi = T,
              type = "GEV",
              method = "GMLE",
              units = "MM"
)


# Model selected for NYC

covariates.list.location = c("LinearTrend",
                             #"ExpTrend")
                             "AnnCos",
                             "AnnSin",
                             #"AnnCos2",
                             #"AnnSin2")
                             "NodalCos",
                             #"NodalSin")
                             "PerigeanCos",
                             "NAO",
                             "NINO",
                             "NAO:NINO")


covariates.list.scale = c(#"LinearTrend",
						  #"ExpTrend")
						  #"AnnCos",
						  "AnnSin",
						  #"AnnCos2")
						  "AnnSin2",
						  "NodalCos",
						  #"NodalSin")
						  "PerigeanCos",
						  "NAO",
						  "NINO",
						  "NAO:NINO")




formula.loc = as.formula(paste("~",  paste(covariates.list.location, collapse = "+")))
formula.scale = as.formula(paste("~",  paste(covariates.list.scale, collapse = "+")))


model_nyc = fevd(x = MaxTideResDeTrend,
                 data = covariates,
                 location.fun = formula.loc,
                 scale.fun = formula.scale,
                 shape = ~1,
                 use.phi = T,
                 type = "GEV",
                 method = "GMLE"
)


pvalue_results_nyc = LR_test_MaxTideRes(covariates = covariates,
                                    covariates.list.location = covariates.list.location,
                                    covariates.list.scale = covariates.list.scale,
                                    use.phi = T,
                                    type = "GEV",
                                    method = "GMLE")


################################ CPH ##########################################

setwd("./data/cph")


data_hourly =  read.csv(file = paste("hourly_series_predict_den", ".csv", sep=""),
                        sep=",",
                        header=TRUE)

data_hourly$TideRes = data_hourly$wl-data_hourly$wlh



max_TideRes = data.frame(aggregate(data_hourly$TideRes,
                                   data.frame(as.yearmon(data_hourly$time)),
                                   FUN="max"))

names(max_TideRes) = c('Time','MaxTideRes')

max_WaterLvl = data.frame(aggregate(data_hourly$wl,
                                    data.frame(as.yearmon(data_hourly$time)),
                                    FUN="max"))

names(max_WaterLvl) = c('Time','MaxWaterLvl')

yearly_median = aggregate(data_hourly$TideRes,
                          data.frame(year(data_hourly$time)),
                          FUN = "median")
names(yearly_median) = c("Year", "Median")

yearly_mean = aggregate(data_hourly$TideRes,
                        data.frame(year(data_hourly$time)),
                        FUN = "mean")
names(yearly_mean) = c("Year", "Mean")


monthly_obs = aggregate(data_hourly$TideRes,
                        data.frame(as.yearmon(data_hourly$time)),
                        FUN = "length")

colnames(monthly_obs) = c("Year", "Observations")

monthly_obs$Observations = monthly_obs$Observations/max(monthly_obs$Observations)



data_monthly = merge(x = max_WaterLvl, 
                     y = max_TideRes,
                     by = "Time",
                     all.y = TRUE)

data_monthly$YearlyMedian = 0
data_monthly$YearlyMean = 0

for ( t in 1:length(data_monthly$Time)){
  
  data_monthly$YearlyMedian[t] = yearly_median$Median[yearly_median$Year == year(data_monthly$Time[t])]
  data_monthly$YearlyMean[t] = yearly_mean$Mean[yearly_mean$Year == year(data_monthly$Time[t])]
  
}

data_monthly$MaxDeTrend = data_monthly$MaxWaterLvl-data_monthly$YearlyMedian

data_monthly$Observation = monthly_obs$Observations

data_monthly = subset(data_monthly,data_monthly$Time>=1950)

data_monthly = subset(data_monthly, data_monthly$Time<2013)

data_monthly$MaxTideResDetrend = data_monthly$MaxTideRes - abs(data_monthly$YearlyMean)

data_monthly = subset(data_monthly, data_monthly$Observation>0.6)

## COVARIATES FOR SWEDEN
covariates = data.frame(Time = data_monthly$Time,
                        LinearTrend = 1:length(data_monthly$Time)/12,
                        ExpTrend = exp((1:length(data_monthly$Time))/12),
                        AnnCos = cos( 2*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnSin = sin( 2*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnCos2 = cos( 4*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnSin2 = sin( 4*pi*c(1:length(data_monthly$Time))/12 ),
                        
                        
                        
                        NodalCos = cos( 2*pi*c(1:length(data_monthly$Time))/(12*18.6)),
                        NodalSin = sin( 2*pi*c(1:length(data_monthly$Time))/(12*18.6)),
                        PerigeanCos =cos( 2*pi*c(1:length(data_monthly$Time))/(12*4.4)),
                        PerigeanSin = sin(2*pi*c(1:length(data_monthly$Time))/(12*4.4))
)


## Climate index: SCAND

scand = read.csv("scand.csv")

scand = subset(scand, scand$YEAR<2013)
time_range = format(seq(ymd('1950-01-01'),ymd('2012-12-31'),by='month'), "%Y-%m")


missing_months = setdiff(as.yearmon(time_range), as.yearmon(data_monthly$Time
) )


index_keep = is.na(match(as.yearmon(time_range), as.yearmon(missing_months)))

length(index_keep[index_keep==F])

scand_values = scand$INDEX

scand_values = scand_values[index_keep]


covariates$SCAND = scand_values

covariates$MaxWaterLvl = data_monthly$MaxWaterLvl
covariates$MaxTideRes = data_monthly$MaxTideRes
covariates$MaxDeTrend = data_monthly$MaxDeTrend
covariates$MaxTideResDeTrend = data_monthly$MaxTideResDetrend

covariates$LagSCAND = lag(covariates$SCAND)

############# GEV REGRESSION ON COPENHAGEN##################################

covariates = covariates[-1,]



initial = c(mean(covariates$MaxTideResDeTrend),
            log(var(covariates$MaxTideResDeTrend)), 0.02)


names(initial) =c("location", "scale", "shape") 



model_stationary = fevd(x =  MaxTideResDeTrend,
                        data = covariates,
                        location.fun = ~1,
                        scale.fun = ~1,
                        shape = ~1,
                        use.phi = T,
                        type = "GEV",
                        method = "GMLE",
                        units = "MM",
                        initial = as.list(initial)
)


###################### MODEL SELECTED  FOR COPENAGHEN ##############################

initial.location = model_stationary$results$par[1]
initial.scale = model_stationary$results$par[2]
initial.shape = model_stationary$results$par[3]

covariates.list.location = c("LinearTrend",
                    #"ExpTrend")
                    "AnnCos",
                    "AnnSin",
                    #"AnnCos2"
                    "AnnSin2",
                    "NodalCos",
                    #"NodalSin")
                    #"PerigeanCos")
                    #"PerigeanSin")
                    "SCAND")

covariates.list.scale = c("LinearTrend",
                    #"ExpTrend")
                    #"AnnCos",
                    #"AnnSin",
                    #"AnnCos2",
                    #"AnnSin2",
                    #"NodalCos",
                    #"NodalSin",
                    #"PerigeanCos",
                    #"PerigeanSin")
                    "SCAND")


# INITIAL VALUES
formula.loc.ols  = as.formula(paste("MaxTideResDeTrend~",  paste(covariates.list.location, collapse = "+")))



ols.location = lm(formula.loc.ols,
                  data = covariates)


formula.scale.ols = as.formula(paste("log(MaxTideResDeTrend^2)~",  paste(covariates.list.scale, collapse = "+")))



ols.scale = lm(formula.scale.ols,
                  data = covariates)


initial = c(as.numeric(ols.location$coefficients),
            as.numeric(ols.scale$coefficients),
            initial.shape)

names(initial) = c( unlist(lapply(1:length(ols.location$coefficients), function(i){
                   paste("mu",i-1, sep = "")})),
                   unlist(lapply(1:length(ols.scale$coefficients), function(i){
                     paste("phi",i-1, sep = "")})),
                   "shape")





formula.scale = as.formula(paste("~",  paste(covariates.list.scale, collapse = "+")))
formula.loc = as.formula(paste("~",  paste(covariates.list.location, collapse = "+")))

model_den = fevd(x = MaxTideResDeTrend,
                 data = covariates,
                 location.fun = formula.loc,
                 scale.fun = formula.scale,
                 shape = ~1,
                 use.phi = T,
                 type = "GEV",
                 method = "GMLE",
                 initial = as.list(initial),
                 verbose = T,
                 iter = 1000
)

#################################### SEQ ###########################################


setwd("./data/seq")


data_hourly =  read.csv(file = paste("hourly_series_predict_bris", ".csv", sep=""),
                        sep=",",
                        header=TRUE)

data_hourly$TideRes = data_hourly$wl-data_hourly$wlh



max_TideRes = data.frame(aggregate(data_hourly$TideRes,
                                   data.frame(as.yearmon(data_hourly$time)),
                                   FUN="max"))

names(max_TideRes) = c('Time','MaxTideRes')

max_WaterLvl = data.frame(aggregate(data_hourly$wl,
                                    data.frame(as.yearmon(data_hourly$time)),
                                    FUN="max"))

names(max_WaterLvl) = c('Time','MaxWaterLvl')

yearly_median = aggregate(data_hourly$TideRes,
                          data.frame(year(data_hourly$time)),
                          FUN = "median")
names(yearly_median) = c("Year", "Median")

yearly_mean = aggregate(data_hourly$TideRes,
                        data.frame(year(data_hourly$time)),
                        FUN = "mean")
names(yearly_mean) = c("Year", "Mean")


monthly_obs = aggregate(data_hourly$TideRes,
                        data.frame(as.yearmon(data_hourly$time)),
                        FUN = "length")

colnames(monthly_obs) = c("Year", "Observations")

monthly_obs$Observations = monthly_obs$Observations/max(monthly_obs$Observations)



data_monthly = merge(x = max_WaterLvl, 
                     y = max_TideRes,
                     by = "Time",
                     all.y = TRUE)

data_monthly$YearlyMedian = 0
data_monthly$YearlyMean = 0

for ( t in 1:length(data_monthly$Time)){
  
  data_monthly$YearlyMedian[t] = yearly_median$Median[yearly_median$Year == year(data_monthly$Time[t])]
  data_monthly$YearlyMean[t] = yearly_mean$Mean[yearly_mean$Year == year(data_monthly$Time[t])]
  
}

data_monthly$MaxDeTrend = data_monthly$MaxWaterLvl-data_monthly$YearlyMedian

data_monthly$Observation = monthly_obs$Observations

data_monthly = subset(data_monthly,data_monthly$Time>=1984)

data_monthly = subset(data_monthly, data_monthly$Time<2019)

data_monthly$MaxTideResDetrend = data_monthly$MaxTideRes - abs(data_monthly$YearlyMean)

data_monthly = subset(data_monthly, data_monthly$Observation>0.6)

# Covariates
covariates = data.frame(Time = data_monthly$Time,
                        LinearTrend = 1:length(data_monthly$Time)/12,
                        ExpTrend = exp((1:length(data_monthly$Time))/12),
                        AnnCos = cos( 2*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnSin = sin( 2*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnCos2 = cos( 4*pi*c(1:length(data_monthly$Time))/12 ),
                        AnnSin2 = sin( 4*pi*c(1:length(data_monthly$Time))/12 ),
                        
                        
                        
                        NodalCos = cos( 2*pi*c(1:length(data_monthly$Time))/(12*18.6)),
                        NodalSin = sin( 2*pi*c(1:length(data_monthly$Time))/(12*18.6)),
                        PerigeanCos =cos( 2*pi*c(1:length(data_monthly$Time))/(12*4.4)),
                        PerigeanSin = sin(2*pi*c(1:length(data_monthly$Time))/(12*4.4))
)






# Climate indexes


iod = read.csv("iod.csv", header = F)

aao = read.csv("aao.csv", header = F)
soi = read.csv("soi.csv")

soi = soi[1:71,]


time_range = format(seq(ymd('1984-01-01'),ymd('2018-12-31'),by='month'), "%Y-%m")


missing_months = setdiff(as.yearmon(time_range), as.yearmon(data_monthly$Time))

index_keep = is.na(match(as.yearmon(time_range), as.yearmon(missing_months)))

iod = subset(iod, iod$V1 >= 1984 &iod$V1< 2019)

iod_values = unlist(iod[,2:13])

iod_values = iod_values[index_keep]

aao = subset(aao, aao$V1 >= 1984 &aao$V1< 2019)

aao_values = unlist(aao[,2:13])

aao_values = aao_values[index_keep]

soi = subset(soi,soi$YEAR>=1984 & soi$YEAR< 2019)

soi_values = as.numeric(unlist(soi[,2:13]))

soi_values = soi_values[index_keep]

covariates$IOD = iod_values

covariates$AAO = aao_values

covariates$SOI = soi_values

covariates$MaxWaterLvl = data_monthly$MaxWaterLvl
covariates$MaxTideRes = data_monthly$MaxTideRes
covariates$MaxDeTrend = data_monthly$MaxDeTrend
covariates$MaxTideResDeTrend = data_monthly$MaxTideResDetrend

covariates = subset(covariates, covariates$MaxTideResDeTrend>0)



########################### GEV REGRESSION BRISBANE ###################

initial = c(mean(covariates$MaxTideResDeTrend),
            log(var(covariates$MaxTideResDeTrend)), 0.002)


names(initial) =c("location", "scale", "shape") 



model_stationary = fevd(x =  MaxTideResDeTrend,
                        data = covariates,
                        location.fun = ~1,
                        scale.fun = ~1,
                        shape = ~1,
                        use.phi = T,
                        type = "GEV",
                        method = "GMLE",
                        units = "MM",
                        initial = as.list(initial),
                        verbose = T
)

initial.location = model_stationary$results$par[1]
initial.scale = model_stationary$results$par[2]
initial.shape = model_stationary$results$par[3]


###################### MODEL SELECTED BRISBANE ################################
covariates.list.location = c("LinearTrend",
                    #"ExpTrend")
                    #"AnnCos")
                    #"AnnSin")
                    #"AnnCos2")
                    "AnnSin2",
                    #"NodalCos")
                    #"NodalSin")
                    #"PerigeanCos")
                    #"PerigeanSin")
                    "IOD",
                    "AAO",
                    "SOI")

covariates.list.scale = c("LinearTrend",
                    #"ExpTrend")
                    #"AnnCos",
                    #"AnnSin",
                    #"AnnCos2",
                    #"AnnSin2",
                    #"NodalCos",
                    #"NodalSin",
                    #"PerigeanCos",
                    #"PerigeanSin")
                    "IOD",
                    "AAO",
                    "SOI")

# INITIAL VALUES
formula.loc.ols  = as.formula(paste("MaxTideResDeTrend~",  paste(covariates.list.location, collapse = "+")))



ols.location = lm(formula.loc.ols,
                  data = covariates)


formula.scale.ols = as.formula(paste("log(MaxTideResDeTrend^2)~",  paste(covariates.list.scale, collapse = "+")))



ols.scale = lm(formula.scale.ols,
               data = covariates)


initial = c(as.numeric(ols.location$coefficients),
            as.numeric(ols.scale$coefficients),
            initial.shape)

names(initial) = c( unlist(lapply(1:length(ols.location$coefficients), function(i){
  paste("mu",i-1, sep = "")})),
  unlist(lapply(1:length(ols.scale$coefficients), function(i){
    paste("phi",i-1, sep = "")})),
  "shape")





formula.scale = as.formula(paste("~",  paste(covariates.list.scale, collapse = "+")))
formula.loc = as.formula(paste("~",  paste(covariates.list.location, collapse = "+")))

model_bris = fevd(x = MaxTideResDeTrend,
                 data = covariates,
                 location.fun = formula.loc,
                 scale.fun = formula.scale,
                 shape = ~1,
                 use.phi = T,
                 type = "GEV",
                 method = "GMLE",
                 initial = as.list(initial),
                 verbose = T,
                 iter = 1000
)

