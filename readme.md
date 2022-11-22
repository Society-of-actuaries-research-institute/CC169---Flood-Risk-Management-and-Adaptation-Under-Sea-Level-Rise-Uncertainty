 # Code for flood risk project using real option analysis

 # This repository contains the relevant code for the project on flood risk and real option analysis
 # There are several parts. 

 # Part 1: Data Collection
 The folder data contains the relevant data for the three case studies in New York City, South East Queensland, and Copenhagen. For sea level data we use hourly data available from the archives of the University of Hawaii Sea Level Centre (http://uhslc.soest.hawaii.edu). Climate indices used in this study are:  NAO and the Ni ̃no 3.4 indices for the US, the SCA for Europe, and  the AAO, the IOD and the SOI for Australia. These indices are collected from different sources, and stored in the relevant data folders.

 # Part 1: Deseasonalization

  Sea water level are deseasonalized using the matlab function "deseasonalization_matlab.m" in the folder matlab_scripts. This function is based on the Utide package in R by Codiga (2011). The function takes the sea level as input and preduces

   When running this function, be sure to:
  1. indicate correctly the input file corresponding to the case study;
  2. use the correct latitude values;
  3. indicate the correct output file.