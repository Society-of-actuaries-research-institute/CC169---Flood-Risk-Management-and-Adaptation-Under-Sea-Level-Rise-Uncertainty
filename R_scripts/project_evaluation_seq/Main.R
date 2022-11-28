rm(list=ls(all=TRUE)) 
#setwd("C:/E drive/Flood/SOA/sent/Report/empirical/qld")
setwd("C:/Users/matte/Desktop/SOA project/Multiple Investments/empirical/qld")

opath= './output/'
fpath= './../function/'
ftpath='./../fortran/'
tab=  './tables/'
fig=  './figures/'
file.sources = list.files(fpath, pattern="*.R$",full.names=TRUE, ignore.case=TRUE)
invisible(sapply(file.sources,source))
source(paste("parameters.R",sep=''))

#== 1. De-tide ===
#source(paste("detide.R",sep=''))

#== 2. Estimate GEV ====
#source(paste("estimation_gev.R",sep=''))

#== 3 Estimate loss curve ===
#source(paste("damage_curve.R",sep=''))

#== 4. Market risk premium for mean sea level process ===
#source(paste("risk_premium.R",sep=''))

#=== 5. Single Investment===
source(paste("investment_analysis.R",sep=''))
source(paste("sensitivity_singleinvestment.R",sep=''))
source(paste("sensitivity_singleinvestment_sigma.R",sep=''))

source(paste("loss_distribution_QLD.R",sep=''))

source(paste("sensitivity_analysis_loss_distribution.R", sep =''))

































