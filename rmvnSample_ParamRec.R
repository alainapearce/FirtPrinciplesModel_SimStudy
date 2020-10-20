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
#     model_str: The base model to use--'FPM' for the first principles model and 'Kissileff' for the quadratic model. Default is 'FPM'.
#     procNoise (optional) A logical indicator for adding random process noise to the bite data by jittering bite size (and thus estimated timing). Default value is TRUE; if FALSE will use average bite size to estimate initail bite timing.
#     measureNoise (optional) A string indicating they type of measurement noise to add. The options include:
#     BiteSize' - will use average bite size for parameter recovery; 'BiteTiming' - add noise to bite timing (jittered);
#     or 'Both' - will apply both types of measurement noise. This noise is applied to bite data after initial
#     parameterization and before parameter recovery. Default is no measurement error.
#     Noise_biteSizeSD (optional) This allows you to enter the standard deviation of individuals bites sizes and will replace the default procNoise routine (jittered bite sizes). Bite sizes will be randomly chosen from a normal distribution truncated at min = 0 with mean = Emax/nBites and standard deviation equal to the entered value. procNoise must be set to TRUE, otherwise this argument will be ignored.
#     @param mNoise_biteTimeSD (optional) This allows you to enter the standard deviation for adjusting bite timing and will
#     replace the default (jittered bite timing). The noise add to each timepoint will be chosen from a normal distribution
#     with mean = 0 and standard deviation entered. measureNoise must be set to to 'BiteTiming' or 'Both' otherwise this
#     argument will be ignored. Note: the normal distribution will be truncated at at each timepoint so that the time for
#     timepoint t is not less than timepoint t-1.
#     mNoise_biteSizeCat (option) This allows you to alter the default for bite size error (average bite size) by
#     entering category cut points or NA to skip this measurement error. Cut points must equal n - 1 categories (e.g., if
#     want three categories you would enter the small-medium and medium-large large cut/boundry points). Cut points will
#     be left/lower inclusive but exclude upper boundary. Bite sizes within each category will be set to the average bite
#     size for that category. This will replace the default measureNoise routine (all bites = average bite size).
#     measureNoise must be set to to 'BiteSize' or 'Both' otherwise this argument will be ignored.write.dat: logical, want to write data to .csv. default <- TRUE - will write to Data directory
#     datOnly: logical, only generate data set and skip parameter recovery, default is TRUE. if FALSE, function will
#                return both the sampled data set and the parameter dataset as separate items in a list
#     paramCI: (optional) A list of strings with the names of the parameters to compute CIs for. Optional. If none specified, no CI will be computed
#     bound: (only required if paramCI is used) A string indicating which confidence bound to return: 'upper, 'lower', or 'both'. Default <- 'both'
#     rmse (optional) A string indicating which measure to compute root mean square error for. Options include: 'timing' for bite timing, 'intake' for cumulative intake', 'both' to compute for both timing and intake. If not specified, rmse is not computed. Default is to not compute rmse.
#     data_str: (optional) A string you want to use to name output dataset - the model used and number of samples will automatically be part of the name. Default is 'simDat' which results in the dataset name(s) if you use 100 samples: data - model_simDat_rmvnDat_rmvnDat100.csv and parameter recovery - data - model_simDat_rmvnDat_rmvnParamRec100.csv
# scaleFactor (optional) A scaling factor to adjust the standard deviation of the multivariate normal distribution. Will be applied to all sampled variables. E.g., a value of 0.5 will scale the standard deviation by half and the variance by a quarter using pre- (S) and post-matrix (S transpose - ST) multiplication of the covariance matrix (C) with the scaling factor on the diagonal (S x C x ST)

#### Set up ####
library(bitemodelr)

# Simulate and extract parameter estimates ####
# Choose different values based on multivariate distribution and compare 
# FMP and Kissileff curves. Fit bite timing with FPM first and then use to get
# cumulative intake from quadratic to avoide the end timing issue

rmvnSample_ParamRec = function(nSample = 100, nSim = 100, model_str = "FPM", procNoise = TRUE, measureNoise = FALSE, pNoise_biteSizeSD = NA, mNoise_biteTimeSD = NA, mNoise_biteSizeCat = "mean", write.dat = TRUE, datOnly = TRUE, paramCI = NA, bound = 'both', rmse = NA, data_str = 'simDat', scaleFactor = NA){
  
  #get rmvn sample
  SimDat_rmvn = rmvnSample(nSample, model_str, write.dat, data_str)
  
  for(l in 1:nrow(SimDat_rmvn)) {
    
    SimDat_rmvn$nBites_round = round(SimDat_rmvn$nBites)
    
    if(model_str == "FPM" | model_str == "Both" | model_str == "both"){
      paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], 
                                parameters = c(SimDat_rmvn$theta[l], SimDat_rmvn$r[l]), 
                                model_str = model_str, nSims = 1, procNoise = procNoise, 
                                measureNoise = measureNoise, pNoise_biteSizeSD = pNoise_biteSizeSD, 
                                mNoise_biteTimeSD = mNoise_biteTimeSD, mNoise_biteSizeCat = mNoise_biteSizeCat, 
                                paramCI = paramCI, bound = bound, rmse = rmse)
    } 
    
    if (model_str == "Kissileff" | model_str == "Both" | model_str == "both"){
      paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], 
                                parameters = c(SimDat_rmvn$int[l], SimDat_rmvn$linear[l], SimDat_rmvn$quad[l]), 
                                model_str = model_str, nSims = 1, procNoise = procNoise, 
                                measureNoise = measureNoise, pNoise_biteSizeSD = pNoise_biteSizeSD, 
                                mNoise_biteTimeSD = mNoise_biteTimeSD, mNoise_biteSizeCat = mNoise_biteSizeCat, 
                                paramCI = paramCI, bound = bound, rmse = rmse)
    }
    
    paramRec$ID <- l
    
    if(l == 1){
      FPM_paramRecDat_rmvn <- paramRec
    } else {
      FPM_paramRecDat_rmvn <- rbind(FPM_paramRecDat_rmvn, paramRec)
    }
  }
  
  #check if param in fitted
  if(!is.na(paramCI)){
    if(bound == 'both'){
      if(model_str == "FPM"){
        FPM_paramRecDat_rmvn$r_fit <- ifelse(FPM_paramRecDat_rmvn$initial_r < FPM_paramRecDat_rmvn$u95CI_r & FPM_paramRecDat_rmvn$initial_r > FPM_paramRecDat_rmvn$l95CI_r, TRUE, FALSE)
        FPM_paramRecDat_rmvn$theta_fit <- ifelse(FPM_paramRecDat_rmvn$initial_theta < FPM_paramRecDat_rmvn$u95CI_theta & FPM_paramRecDat_rmvn$initial_theta > FPM_paramRecDat_rmvn$l95CI_theta, TRUE, FALSE)
      } else if (model_str == "Kissileff"){
        FPM_paramRecDat_rmvn$int_fit <- ifelse(FPM_paramRecDat_rmvn$initial_int < FPM_paramRecDat_rmvn$u95CI_int & FPM_paramRecDat_rmvn$initial_int > FPM_paramRecDat_rmvn$l95CI_int, TRUE, FALSE)
        FPM_paramRecDat_rmvn$linear_fit <- ifelse(FPM_paramRecDat_rmvn$initial_linear < FPM_paramRecDat_rmvn$u95CI_linear & FPM_paramRecDat_rmvn$initial_linear > FPM_paramRecDat_rmvn$l95CI_linear, TRUE, FALSE)
        FPM_paramRecDat_rmvn$quad_fit <- ifelse(FPM_paramRecDat_rmvn$initial_quad < FPM_paramRecDat_rmvn$u95CI_quad & FPM_paramRecDat_rmvn$initial_quad > FPM_paramRecDat_rmvn$l95CI_quad, TRUE, FALSE)
      }
    }
  }
  
  
  if(isTRUE(write.dat)){
    if(!is.na(paramCI)){
      if(!is.na(rmse)){
        write.csv(FPM_paramRecDat_rmvn, paste0('Data/', model_str, '_', data_str, '_rmvnParamRecCI_RMSE', nSample, '.csv'), row.names=FALSE)
      } else {
        write.csv(FPM_paramRecDat_rmvn, paste0('Data/', model_str, '_', data_str, '_rmvnParamRecCI', nSample, '.csv'), row.names=FALSE)
      }
    } else if (!is.na(rmse)){
      write.csv(FPM_paramRecDat_rmvn, paste0('Data/', model_str, '_', data_str, '_rmvnParamRecRMSE', nSample, '.csv'), row.names=FALSE)
    }
  }  
  
  SimDat_rmvn_list <- list(SimDat_rmvnDat <- SimDat_rmvn[1:4],
                           SimDat_rmvnParamDat <- FPM_paramRecDat_rmvn)
  
  return(SimDat_rmvn_list)
}

