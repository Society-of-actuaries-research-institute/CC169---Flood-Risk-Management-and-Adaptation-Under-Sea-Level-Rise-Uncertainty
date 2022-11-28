 # Code for flood risk project using real option analysis

 This repository contains the relevant code for the project on flood risk and real option analysis. There are several parts. 

 # Part 1: Data Collection
 The folder data contains the relevant data for the three case studies in New York City, South East Queensland, and Copenhagen. For sea level data we use hourly data available from the archives of the University of Hawaii Sea Level Centre (http://uhslc.soest.hawaii.edu). Climate indices used in this study are:  NAO and the Ni ̃no 3.4 indices for the US, the SCA for Europe, and  the AAO, the IOD and the SOI for Australia. These indices are collected from different sources, and stored in the relevant data folders.

 # Part 2: Deseasonalization

  Sea water level are deseasonalized using the matlab function "deseasonalization_matlab.m" in the folder matlab_scripts. This function is based on the Utide package in R by Codiga (2011). The function takes the sea level as input and estimates tide parameters using 4 year hourly data collected from the University of Hawaii ("hourly_series_2014.csv" for NYC, "hourly_series_2014_bris.csv" for SEQ, and "hourly_series_2014.csv" for CPH), and then using estimated tide parameters, predict the tide component for the whole hourly data set. The output ("hourly_series_predict.csv", "hourly_series_predict_bris.csv", or "hourly_series_predict_den.csv") now has two variables: wl is the original water level with trend and tide; wh is the predicted tide (seasonality component)

  When running this function, be sure to:
  1. indicate correctly the input file corresponding to the case study;
  2. use the correct latitude values;
  3. indicate the correct output file.

  # Part 3: Statistical Analysis of Extreme Sea Water Level Events
  The maxima of the non-tidal residual (the time series of the sea level deseasonalized) can then be modelled using extreme value theory and the GAMLSS approach: the monthly maxima of the non tidal residual are assume to follow a generalzied extreme value theory distribution with location and scale parameters depending on covariates. The shape parameter is assumed to be stationary for stability issues. The file called "GEV_regression.R" in R_scripts contains the codes for replicate what is presented in the report in section 3.5 and based on the framework discussed in sections 3.4, 3.4.1, and 3.4.2. The initial values for the maxiumum likelihood estimation are taken from the stationary model, i.e. the model without covariates. Covariates selection is done following the discussion in section 3.4.2 of the report. The significativity of each coeffcient is investigated using a likelihood ratio test via the function "LR_test_MaxTideRes" in the R file. Results are stored in model_nyc, model_bris, and model_den. 

  # Part 4: Climate Adaptation Project Evaluation
  This part evaluates climate adaptation projects for the three case studies. It is based on Section 4 of the report, and can be used to replicate the numerical results of Section 5. The relevant code can be found in R_scripts and it is divided in three folders "project_evalution_nyc", "project_evalution_seq", and "project_evalution_den" corresponding to the NYC, SEQ, and CPH case studies respectively. In the case study of New York City, we consider single investments of i) a barrier and dyke project and ii) a water proofing project, as well as a multiple investment when the two projects are considered together. Consideration of multiple investments is necessary for the formulation of dynamic investment pathways. In South East Queensland, we investigate the optimal time to invest in a house elevation project. In the Copenhagen,we consider the problem of investing in a dyke system where the ultimate height of the dyke is achieved in one go (a single investment) or achieved in two stages (a multistage investment).

  In each case, the file "Main.R" is used to execute the relevant section of the code:
  1. "parameters.R": loads the parameters needed for the project evalaution. The choice of the paramters depends on the estimation results for the generalized extreme value distribution in each case, and on previous findings in the literature.
  2. "investment_analysis.R": performs the project evaluation, considering the decision to invest into the climate adaptation policy as an american option. 4 scenarios are considered: no climate change impact and growth exposure; climate change impact and growth exposure; high climate change impact and growth exposure; climate change impact and high growth exposure. For each climate adaptation policy analyzed, the code outputs the project value measured both as the Net Present Value and as an american option, stored in the variable results, and the investment boundaries. The variable results is used to construct the Tables in Section 5.The numerical approximation is done via binomial tree (see Section 4 of the report). 
  3. "sensitivity_analysis_discount.R", "sensitivity_analysis_slr.R", "sensitivity_analysis_sigma.R": these three files contains the code relevant to the sensitivity analysis with repect to the dscount rate, the sea level rise, and the sea level rise uncertainty. Each fine can be run separately, and the variable result is used to construct the Tables in Section 5 of the report.

  # Part 5: Premium Distribution Evolution

  








