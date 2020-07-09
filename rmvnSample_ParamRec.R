# This script was written by Alaina Pearce in 2020 
# to generate the data used to test parameter recovery
# when sampling from the multivariate normal distribution for
# number of bites, total intake (grams), and model parameters
#     FMP model includes r and theta as parameter
#     Kissileff model includes intercept, linear, and quadratic coefficients
#
# Input Arguments:
#     nSample: number of samples to be randomly pulled from multivariate normal distribution, default = 100
#     nSims: number of simulations per sample (right now programed for just basic recovery with procNoise); only needed if 
#            datONLY = FALSE. default = 100
#     model: which model to use, default = "FPM"
#     write.dat: logical, want to write data to .csv. default = TRUE - will write to Data directory
#     datOnly: logical, only generate data set and skip parameter recover, default is TRUE. if FALSE< function will
#                return both the sampled data set and the parameter dataset as separate items in a list
#     keepBites: logical, only relevat if datONLY = FALSE. if TRUE will return bite data from parameter simulation as
#                 another item in the list returned
#
#
#     Copyright (C) 20120 Alaina L Pearce
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

#### Set up ####
library(bitemodelr)
library(MASS)


#if running outside of a script that loads the data, need to load data here
#SimDat_Fogel2017 = read.csv('Data/ParamDat_Fogel2017.csv')

# Simulate and extract parameter estimates ####
# Choose different values based on multivariate distribution and compare 
# FMP and Kissileff curves. Fit bite timing with FPM first and then use to get
# cumulative intake from quadratic to avoide the end timing issue

mvnSample_ParamRec = function(nSample = 100, nSim = 100, model = "FPM", write.dat = TRUE, datOnly = TRUE, keepBites = FALSE){
  set.seed(1234)
  
  TotalSamples = 0
  while(TotalSamples < nSample){
    rowsNeed = 100 - TotalSamples
    
    if(rowsNeed < 10){
      newsample = 10
    } else {
      newsample = rowsNeed
    }
    
    if(model == "FPM"){
      rmv_dat = as.data.frame(mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 10:11)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 10:11)]), empirical = TRUE))
      
      
      rmv_dat$r_check = (-1*rmv_dat$theta)/rmv_dat$TotalIntake_g
      
      rmv_datKeep = rmv_dat[rmv_dat$r_check < rmv_dat$r, ]
      rmv_datKeep = rmv_datKeep[rmv_datKeep$r >= min(SimDat_Fogel2017$r) & rmv_datKeep$r <= max(SimDat_Fogel2017$r), ]
      rmv_datKeep = rmv_datKeep[rmv_datKeep$theta >= min(SimDat_Fogel2017$theta) & rmv_datKeep$theta <= max(SimDat_Fogel2017$theta), ]
    } else if(model == "Kissilef"){
      rmv_dat = as.data.frame(mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 17:19)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 17:19)]), empirical = TRUE))
      
      rmv_dat$sqrt_term = rmv_dat$linear^2 - 4 * (rmv_dat$int - rmv_dat$TotalIntake_g) * rmv_dat$quad
      
      rmv_datKeep = rmv_dat[rmv_dat$sqrt_term > 0, ]
      rmv_datKeep = rmv_datKeep[rmv_datKeep$int >= min(SimDat_Fogel2017$int) & rmv_datKeep$int <= max(SimDat_Fogel2017$int), ]
      rmv_datKeep = rmv_datKeep[rmv_datKeep$linear >= min(SimDat_Fogel2017$linear) & rmv_datKeep$linear <= max(SimDat_Fogel2017$linear), ]
      rmv_datKeep = rmv_datKeep[rmv_datKeep$quad >= min(SimDat_Fogel2017$quad) & rmv_datKeep$linear <= max(SimDat_Fogel2017$quad), ]
    }
    
    rmv_datKeep = rmv_datKeep[rmv_datKeep$TotalIntake_g >= min(SimDat_Fogel2017$TotalIntake_g) & rmv_datKeep$TotalIntake_g <= max(SimDat_Fogel2017$TotalIntake_g), ]
    rmv_datKeep = rmv_datKeep[rmv_datKeep$nBites >= min(SimDat_Fogel2017$nBites) & rmv_datKeep$nBites <= max(SimDat_Fogel2017$nBites), ]
    
    if(TotalSamples == 0){
      SimDat_rmvn = rmv_datKeep
    } else if(nrow(rmv_datKeep) > rowsNeed){
      rmv_datKeep = rmv_datKeep[sample(nrow(rmv_datKeep), rowsNeed), ]
      SimDat_rmvn = rbind(SimDat_rmvn, rmv_datKeep)
    } else if(nrow(rmv_datKeep) <= TotalSamples){
      SimDat_rmvn = rbind(SimDat_rmvn, rmv_datKeep)
    }
    
    TotalSamples = nrow(SimDat_rmvn)
  }
  
  
  if(isTRUE(datOnly)){
    if(isTRUE(write.dat)){
      write.csv(SimDat_rmvn, 'Data/SimDat_rmvn.csv', row.names = FALSE)
    }
    
    return(SimDat_rmvn)
    
  } else if(isFALSE(datOnly)){
    #get standard parameter recovery with process noise
    for(l in 1:nrow(SimDat_rmvn)) {
      
      if(model == "FPM"){
        paramRec = ParamRecovery(nBites = SimDat_rmvn$nBites[l], Emax = SimDat_rmvn$TotalIntake_g[l], parameters = c(SimDat_rmvn$theta[l], SimDat_rmvn$r[l]), time_fn = FPM_Time, fit_fn = FPM_Fit, keepBites = FALSE, intake_fn = FPM_Intake, CI = TRUE, nSims = 100)
      }
      
      paramRec$ID = l
      paramRec$r_fit = ifelse(paramRec$initial_r < paramRec$u95CI_r & paramRec$initial_r > paramRec$l95CI_r, TRUE, FALSE)
      paramRec$theta_fit = ifelse(paramRec$initial_theta < paramRec$u95CI_theta & paramRec$initial_theta > paramRec$l95CI_theta, TRUE, FALSE)
      
      SimDat_rmvn$rFit_n[l] = nrow(paramRec[paramRec$r_fit, ])
      SimDat_rmvn$rFit_p[l] = SimDat_rmvn$rFit_n[l]/100
      
      SimDat_rmvn$thetaFit_n[l] = nrow(paramRec[paramRec$theta_fit, ])
      SimDat_rmvn$thetaFit_p[l] = SimDat_rmvn$thetaFit_n[l]/100
      
      if(l == 1){
        FPM_paramRecDat_rmvn = paramRec
      } else {
        FPM_paramRecDat_rmvn = rbind(FPM_paramRecDat_rmvn, paramRec)
      }
    }
    
    if(isTRUE(write.dat)){
      write.csv(SimDat_rmvn, 'Data/SimDat_rmvn.csv', row.names = FALSE)
      if(isTRUE(keepBites)){
        write.csv(FPM_paramRecDat_rmvn$biteDat_paramRecov, 'Data/SimDat_rmvnBiteDat.csv', row.names = FALSE)
        write.csv(FPM_paramRecDat_rmvn$paramDat, 'Data/SimDat_rmvnParamRec.csv', row.names = FALSE)
      } else {
        write.csv(FPM_paramRecDat_rmvn, 'Data/SimDat_rmvnParamRec.csv', row.names = FALSE)
      }
    }
    
    if(isTRUE(keepBites)){
      SimDat_rmvn_list = list(SimDat_rmvnDat = SimDat_rmvn,
                             SimDat_rmvnParamDat = FPM_paramRecDat_rmvn$paramDat,
                             SSimDat_rmvnBiteDat = FPM_paramRecDat_rmvn$biteDat_paramRecov)
      
    } else {
      SimDat_rmvn_list = list(SimDat_rmvnDat = SimDat_rmvn,
                             SimDat_rmvnParamDat = FPM_paramRecDat_rmvn)
      
    }
    return(SimDat_rmvn_list)
  }
}
