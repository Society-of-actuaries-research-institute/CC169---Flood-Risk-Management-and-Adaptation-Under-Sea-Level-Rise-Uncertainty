rm(list=ls(all=TRUE)) 
# add your working path here:
setwd("...")
##
setwd("./project_evaluation_nyc")
opath= './output/'
fpath= './../function/'
ftpath='./../fortran/'
tab=  './output/'
fig=  './output/'
file.sources = list.files(fpath, pattern="*.R$",full.names=TRUE, ignore.case=TRUE)
invisible(sapply(file.sources,source))


#== 1. Get parameters for NYC case study
source(paste("parameters.R",sep=''))

#== 2. Project evaluation using binomial tree
source(paste("investment_analysis.R",sep=''))

#== 3. Sensitivity Analysis Project Evaluation
source(paste("sensitivity_analysis_discount.R",sep=''))
source(paste("sensitivity_analysis_slr.R",sep=''))
source(paste("sensitivity_analysis_sigma.R",sep=''))


#== 4. Loss Distribution in 100 Years
source(paste("loss_distribution_NYC.R",sep=''))
#== 5. Sensistivity Analysis Loss distribution
source(paste("sensitivity_analysis_loss_distribution.R",sep=''))






























