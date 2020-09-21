# This script was written by Alaina Pearce in 2020 
# to generate the data used to test parameter recovery
# when sampling from the multivariate normal distribution for
# number of bites, total intake (grams), and model parameters
#     FMP model includes r and theta as parameter
#     Kissileff model includes intercept, linear, and quadratic coefficients
#
#     Copyright (C) 2020 Alaina L Pearce
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
#
# Input Arguments:
#     nSample: number of samples to be randomly pulled from multivariate normal distribution, default <- 100
#     nSims: number of simulations per sample (right now programed for just basic recovery with procNoise); only needed if 
#            datONLY <- FALSE. default <- 100
#     model: which model to use, default <- "FPM"
#     write.dat: logical, want to write data to .csv. default <- TRUE - will write to Data directory
#     datOnly: logical, only generate data set and skip parameter recovery, default is TRUE. if FALSE, function will
#                return both the sampled data set and the parameter dataset as separate items in a list
#     paramCI: A list of strings with the names of the parameters to compute CIs for. Optional. If none specified, no CI will be computed
#     bound: A string indicating which confidence bound to return: 'upper, 'lower', or 'both'. Default <- 'both'
#     data_str: (optional) A string you want to use to name ouput dataset - the model used and number of samples will automatically be part of the name. Default is 'simDat' which results in the dataset name(s) if you use 100 samples: data - model_simDat_rmvnDat_rmvnDat100.csv and parameter recovery - data - model_simDat_rmvnDat_rmvnParamRec100.csv

#### Set up ####
library(bitemodelr)
library(MASS)


#if running outside of a script that loads the data, need to load data here
#SimDat_Fogel2017 <- read.csv('Data/ParamDat_Fogel2017.csv')

# Simulate and extract parameter estimates ####
# Choose different values based on multivariate distribution and compare 
# FMP and Kissileff curves. Fit bite timing with FPM first and then use to get
# cumulative intake from quadratic to avoide the end timing issue

rmvnSample_ParamRec = function(nSample = 100, nSim = 100, procNoise = TRUE, model = "FPM", write.dat = TRUE, datOnly = TRUE, paramCI, bound = 'both', data_str = 'simDat'){
  set.seed(1234)
  
  TotalSamples <- 0
  while(TotalSamples < nSample){
    rowsNeed <- 100 - TotalSamples
    
    if(rowsNeed < 10){
      newsample <- 10
    } else {
      newsample <- rowsNeed
    }
    
    if(model == "FPM"){
      rmvn_dat <- as.data.frame(mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 10:11)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 10:11)]), empirical = TRUE))
      
      
      rmvn_dat$r_check <- (-1*rmvn_dat$theta)/rmvn_dat$TotalIntake_g
      
      rmvn_datKeep <- rmvn_dat[rmvn_dat$r_check < rmvn_dat$r, ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$r >= min(SimDat_Fogel2017$r) & rmvn_datKeep$r <= max(SimDat_Fogel2017$r), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$theta >= min(SimDat_Fogel2017$theta) & rmvn_datKeep$theta <= max(SimDat_Fogel2017$theta), ]
    
      
      } else if(model == "Kissileff"){
      rmvn_dat <- as.data.frame(mvrnorm(newsample, mu = colMeans(SimDat_Fogel2017[c(2, 7, 17:19)]), Sigma = cov(SimDat_Fogel2017[c(2, 7, 17:19)]), empirical = TRUE))   
      
      rmvn_datKeep <- rmvn_dat[rmvn_dat$int >= min(SimDat_Fogel2017$int) & rmvn_dat$int <= max(SimDat_Fogel2017$int), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$linear >= min(SimDat_Fogel2017$linear) & rmvn_datKeep$linear <= max(SimDat_Fogel2017$linear), ]
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$quad >= min(SimDat_Fogel2017$quad) & rmvn_datKeep$linear <= max(SimDat_Fogel2017$quad), ]
      }
    
    
    rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$nBites >= min(SimDat_Fogel2017$nBites) & rmvn_datKeep$nBites <= max(SimDat_Fogel2017$nBites), ]
    rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$TotalIntake_g >= min(SimDat_Fogel2017$TotalIntake_g) & rmvn_datKeep$TotalIntake_g <= max(SimDat_Fogel2017$TotalIntake_g), ]
    
    if (nrow(rmvn_datKeep) > 0){
      rmvn_datKeep$time_calc <- 'Y'
      
      for (r in 1:nrow(rmvn_datKeep)){
        
        grams.bite_avg <- rep(rmvn_datKeep$TotalIntake_g[r]/round(rmvn_datKeep$nBites[r]), round(rmvn_datKeep$nBites[r]))
        # get cumulative intake
        grams.cumulative_avg <- cumsum(grams.bite_avg)
        
        message_long <- rep(FALSE, round(rmvn_datKeep$nBites[r]))
        
        # get long list of parameters
        if (model == 'FPM') {
          params_long <- rep(list(c(rmvn_datKeep$theta[r], rmvn_datKeep$r[r])), round(rmvn_datKeep$nBites[r]))
          simTime <- mapply(FPM_Time, intake = grams.cumulative_avg,
                            parameters = params_long, Emax = rmvn_datKeep$TotalIntake_g[r], message = message_long)
        } else if (model == 'Kissileff'){
          params_long <- rep(list(c(rmvn_datKeep$int[r], rmvn_datKeep$linear[r], rmvn_datKeep$quad[r])), round(rmvn_datKeep$nBites[r]))
          simTime <- mapply(Kissileff_Time, intake <- grams.cumulative_avg,
                            parameters <- params_long, message <- message_long)
        }
        
        if(length(unlist(simTime)) != round(rmvn_datKeep$nBites[r])){
          rmvn_datKeep$time_calc[r] <- 'N'
        } 
        
        n_negTime = sum(unlist(simTime) < 0)
        
        if (n_negTime > 0){
          rmvn_datKeep$time_calc[r] <- 'N'
        }
      }
      
      rmvn_datKeep <- rmvn_datKeep[rmvn_datKeep$time_calc == 'Y', ]

    }
    

    if(TotalSamples == 0){
      SimDat_rmvn <- rmvn_datKeep
    } else if(nrow(rmvn_datKeep) > rowsNeed){
      rmvn_datKeep <- rmvn_datKeep[sample(nrow(rmvn_datKeep), rowsNeed), ]
      SimDat_rmvn <- rbind(SimDat_rmvn, rmvn_datKeep)
    } else if(nrow(rmvn_datKeep) <= TotalSamples){
      SimDat_rmvn <- rbind(SimDat_rmvn, rmvn_datKeep)
    }
    
    TotalSamples <- nrow(SimDat_rmvn)
  }
  
  if(isTRUE(write.dat)){
    if(isTRUE(procNoise)){
      write.csv(SimDat_rmvn[1:ncol(SimDat_rmvn)-1], paste0('Data/', model, '_', data_str, '_rmvnDat', nSample, '.csv'), row.names = FALSE)
  } else {
    write.csv(SimDat_rmvn[1:ncol(SimDat_rmvn)-1], paste0('Data/', model, '_', data_str, '_rmvnDat', nSample, '.csv'), row.names = FALSE)
    }
  }

  if(isTRUE(datOnly)){
    return(SimDat_rmvn)
    
  } else if(isFALSE(datOnly)){
    #get standard parameter recovery with process noise
    
    for(l in 1:nrow(SimDat_rmvn)) {
      
      SimDat_rmvn$nBites_round = round(SimDat_rmvn$nBites)
      
      if(model == "FPM"){
        if(hasArg(paramCI)){
          paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], parameters = c(SimDat_rmvn$theta[l], SimDat_rmvn$r[l]), time_fn = FPM_Time, fit_fn = FPM_Fit, nSims = 1, procNoise = procNoise, intake_fn = FPM_Intake, paramCI = paramCI, bound = bound)
        } else {
          paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], parameters = c(SimDat_rmvn$theta[l], SimDat_rmvn$r[l]), time_fn = FPM_Time, fit_fn = FPM_Fit, nSims = 1, procNoise = procNoise, intake_fn = FPM_Intake)
        }
        
      } else if (model == "Kissileff"){
        if(hasArg(paramCI)){
          paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], parameters = c(SimDat_rmvn$int[l], SimDat_rmvn$linear[l], SimDat_rmvn$quad[l]), time_fn = Kissileff_Time, fit_fn = Kissileff_Fit, nSims = 1, procNoise = procNoise, intake_fn = Kissileff_Intake, paramCI = paramCI, bound = bound)
        } else {
          paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], parameters = c(SimDat_rmvn$int[l], SimDat_rmvn$linear[l], SimDat_rmvn$quad[l]), time_fn = Kissileff_Time, fit_fn = Kissileff_Fit, nSims = 1, procNoise = procNoise, intake_fn = Kissileff_Intake)
        }
        
        paramRec$ID <- l
      }
      
      if(l == 1){
        FPM_paramRecDat_rmvn <- paramRec
      } else {
        FPM_paramRecDat_rmvn <- rbind(FPM_paramRecDat_rmvn, paramRec)
      }
    }
    
    #check if param in fitted
    if(hasArg(paramCI)){
      if(bound == 'both'){
        if(model == "FPM"){
          FPM_paramRecDat_rmvn$r_fit <- ifelse(FPM_paramRecDat_rmvn$initial_r < FPM_paramRecDat_rmvn$u95CI_r & FPM_paramRecDat_rmvn$initial_r > FPM_paramRecDat_rmvn$l95CI_r, TRUE, FALSE)
          FPM_paramRecDat_rmvn$theta_fit <- ifelse(FPM_paramRecDat_rmvn$initial_theta < FPM_paramRecDat_rmvn$u95CI_theta & FPM_paramRecDat_rmvn$initial_theta > FPM_paramRecDat_rmvn$l95CI_theta, TRUE, FALSE)
        } else if (model == "Kissileff"){
          FPM_paramRecDat_rmvn$int_fit <- ifelse(FPM_paramRecDat_rmvn$initial_int < FPM_paramRecDat_rmvn$u95CI_int & FPM_paramRecDat_rmvn$initial_int > FPM_paramRecDat_rmvn$l95CI_int, TRUE, FALSE)
          FPM_paramRecDat_rmvn$linear_fit <- ifelse(FPM_paramRecDat_rmvn$initial_linear < FPM_paramRecDat_rmvn$u95CI_linear & FPM_paramRecDat_rmvn$initial_linear > FPM_paramRecDat_rmvn$l95CI_linear, TRUE, FALSE)
          FPM_paramRecDat_rmvn$quad_fit <- ifelse(FPM_paramRecDat_rmvn$initial_quad < FPM_paramRecDat_rmvn$u95CI_quad & FPM_paramRecDat_rmvn$initial_quad > FPM_paramRecDat_rmvn$l95CI_quad, TRUE, FALSE)
        }
      }
    }
    
    
    if(isTRUE(write.dat)){
      if(isTRUE(procNoise)){
        write.csv(FPM_paramRecDat_rmvn, paste0('Data/', model, '_', data_str, '_rmvnParamRec', nSample, '.csv'), row.names=FALSE)
      } else {
        write.csv(FPM_paramRecDat_rmvn, paste0('Data/', model, '_', data_str, '_rmvnParamRec', nSample, '.csv'), row.names=FALSE)
      }
    }  
    
    SimDat_rmvn_list <- list(SimDat_rmvnDat <- SimDat_rmvn[1:4],
                              SimDat_rmvnParamDat <- FPM_paramRecDat_rmvn)
      
    return(SimDat_rmvn_list)
  }
}
