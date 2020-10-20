# This script was written by Alaina Pearce in 2020 
# to generate distributions of RMSE and/or parameter recovery 
# using multiple models (e.g., Quadratic to generate bite data and 
# FPM to recover parameters)
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
library(ggplot2)

source('rmvnSample_ParamRec.R')
# Get random sample ####
rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'BiteSize',
                                          model_str = "FPM", datOnly = FALSE, write.dat = FALSE,
                                          data_str = 'simDat_procNoise_measureNoise', rmse = 'both')

# Estimate bite data with model 1 ####
for(l in 1:nrow(SimDat_rmvn)) {
  
  SimDat_rmvn$nBites_round = round(SimDat_rmvn$nBites)
  
  if(model_str == "FPM"){
    paramRec <- ParamRecovery(nBites = SimDat_rmvn$nBites_round[l], Emax = SimDat_rmvn$TotalIntake_g[l], 
                              parameters = c(SimDat_rmvn$theta[l], SimDat_rmvn$r[l]), 
                              model_str = model_str, nSims = 1, procNoise = procNoise, 
                              measureNoise = measureNoise, pNoise_biteSizeSD = pNoise_biteSizeSD, 
                              mNoise_biteTimeSD = mNoise_biteTimeSD, mNoise_biteSizeCat = mNoise_biteSizeCat, 
                              paramCI = paramCI, bound = bound, rmse = rmse)
  } else if (model_str == "Kissileff"){
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

# Recover/estimate RMSE ####
#Test RMSE - estimation per 500 randomly sampled microstructure sets ####
##loop to get 500 samples -- only has to be done once

## FPM ##

##no noise
# for (l in 1:5){
#   FPM_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = FALSE, model_str = "FPM",
#                                          datOnly = FALSE, write.dat = FALSE, data_str = 'simDat', rmse = 'both')
# 
#   if (l == 1){
#     FPM_rmvnDat100_RMSE = FPM_rmvn100_rmse_subset[[1]]
#     FPM_rmvnParamRec_RMSE = FPM_rmvn100_rmse_subset[[2]]
#   } else {
#     FPM_rmvnDat100_RMSE = rbind(FPM_rmvnDat100_RMSE, FPM_rmvn100_rmse_subset[[1]])
#     FPM_rmvnParamRec_RMSE = rbind(FPM_rmvnParamRec_RMSE, FPM_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(FPM_rmvnDat100_RMSE, 'Data/FPM_simDat_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(FPM_rmvnParamRec_RMSE, 'Data/FPM_simDat_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

# process noise
# for (l in 1:5){
#   FPM_procNoise_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, model_str = "FPM",
#                                                    datOnly = FALSE, write.dat = FALSE, data_str = 'simDat_procNoise', rmse = 'both')
# 
#   if (l == 1){
#     FPM_procNoise_rmvnDat100_RMSE = FPM_procNoise_rmvn100_rmse_subset[[1]]
#     FPM_procNoise_rmvnParamRec_RMSE = FPM_procNoise_rmvn100_rmse_subset[[2]]
#   } else {
#     FPM_procNoise_rmvnDat100_RMSE = rbind(FPM_procNoise_rmvnDat100_RMSE, FPM_procNoise_rmvn100_rmse_subset[[1]])
#     FPM_procNoise_rmvnParamRec_RMSE = rbind(FPM_procNoise_rmvnParamRec_RMSE, FPM_procNoise_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(FPM_procNoise_rmvnDat100_RMSE, 'Data/FPM_simDat_procNoise_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(FPM_procNoise_rmvnParamRec_RMSE, 'Data/FPM_simDat_procNoise_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process and measurement noise (bite size)
# for (l in 1:5){
#   FPM_procNoise_measureNoiseSize_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'BiteSize',
#                                                                 model_str = "FPM", datOnly = FALSE, write.dat = FALSE,
#                                                                 data_str = 'simDat_procNoise_measureNoise', rmse = 'both')
#   if (l == 1){
#     FPM_procNoise_measureNoiseSize_rmvnDat100_RMSE = FPM_procNoise_measureNoiseSize_rmvn100_rmse_subset[[1]]
#     FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE = FPM_procNoise_measureNoiseSize_rmvn100_rmse_subset[[2]]
#   } else {
#     FPM_procNoise_measureNoiseSize_rmvnDat100_RMSE = rbind(FPM_procNoise_measureNoiseSize_rmvnDat100_RMSE, FPM_procNoise_measureNoiseSize_rmvn100_rmse_subset[[1]])
#     FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE = rbind(FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE, FPM_procNoise_measureNoiseSize_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(FPM_procNoise_measureNoiseSize_rmvnDat100_RMSE, 'Data/FPM_simDat_procNoise_measureNoiseSize_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE, 'Data/FPM_simDat_procNoise_measureNoiseSize_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process and measurement noise (jittered timing)
# for (l in 1:5){
#   FPM_procNoise_measureNoiseTime_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'BiteTiming',
#                                                                        model_str = "FPM", datOnly = FALSE, write.dat = FALSE,
#                                                                        data_str = 'simDat_procNoise_measureNoise', rmse = 'both')
#   if (l == 1){
#     FPM_procNoise_measureNoiseTime_rmvnDat100_RMSE = FPM_procNoise_measureNoiseTime_rmvn100_rmse_subset[[1]]
#     FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE = FPM_procNoise_measureNoiseTime_rmvn100_rmse_subset[[2]]
#   } else {
#     FPM_procNoise_measureNoiseTime_rmvnDat100_RMSE = rbind(FPM_procNoise_measureNoiseTime_rmvnDat100_RMSE, FPM_procNoise_measureNoiseTime_rmvn100_rmse_subset[[1]])
#     FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE = rbind(FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE, FPM_procNoise_measureNoiseTime_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(FPM_procNoise_measureNoiseTime_rmvnDat100_RMSE, 'Data/FPM_simDat_procNoise_measureNoiseTime_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE, 'Data/FPM_simDat_procNoise_measureNoiseTime_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process and measurement noise (bite size and jittered timing)
# for (l in 1:5){
#   FPM_procNoise_measureNoise_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'both',
#                                                                        model_str = "FPM", datOnly = FALSE, write.dat = FALSE,
#                                                                        data_str = 'simDat_procNoise_measureNoise', rmse = 'both')
#   if (l == 1){
#     FPM_procNoise_measureNoise_rmvnDat100_RMSE = FPM_procNoise_measureNoise_rmvn100_rmse_subset[[1]]
#     FPM_procNoise_measureNoise_rmvnParamRec_RMSE = FPM_procNoise_measureNoise_rmvn100_rmse_subset[[2]]
#   } else {
#     FPM_procNoise_measureNoise_rmvnDat100_RMSE = rbind(FPM_procNoise_measureNoise_rmvnDat100_RMSE, FPM_procNoise_measureNoise_rmvn100_rmse_subset[[1]])
#     FPM_procNoise_measureNoise_rmvnParamRec_RMSE = rbind(FPM_procNoise_measureNoise_rmvnParamRec_RMSE, FPM_procNoise_measureNoise_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(FPM_procNoise_measureNoise_rmvnDat100_RMSE, 'Data/FPM_simDat_procNoise_measureNoise_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(FPM_procNoise_measureNoise_rmvnParamRec_RMSE, 'Data/FPM_simDat_procNoise_measureNoise_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## Kissileff
#no noise
# for (l in 1:5){
#   Kissileff_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = FALSE, model_str = "Kissileff",
#                                          datOnly = FALSE, write.dat = FALSE, data_str = 'simDat', rmse = 'both')
# 
#   if (l == 1){
#     Kissileff_rmvnDat100_RMSE = Kissileff_rmvn100_rmse_subset[[1]]
#     Kissileff_rmvnParamRec_RMSE = Kissileff_rmvn100_rmse_subset[[2]]
#   } else {
#     Kissileff_rmvnDat100_RMSE = rbind(Kissileff_rmvnDat100_RMSE, Kissileff_rmvn100_rmse_subset[[1]])
#     Kissileff_rmvnParamRec_RMSE = rbind(Kissileff_rmvnParamRec_RMSE, Kissileff_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(Kissileff_rmvnDat100_RMSE, 'Data/Kissileff_simDat_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(Kissileff_rmvnParamRec_RMSE, 'Data/Kissileff_simDat_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process noise
# for (l in 1:5){
#   Kissileff_procNoise_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, model_str = "Kissileff",
#                                                    datOnly = FALSE, write.dat = FALSE, data_str = 'simDat_procNoise', rmse = 'both')
# 
#   if (l == 1){
#     Kissileff_procNoise_rmvnDat100_RMSE = Kissileff_procNoise_rmvn100_rmse_subset[[1]]
#     Kissileff_procNoise_rmvnParamRec_RMSE = Kissileff_procNoise_rmvn100_rmse_subset[[2]]
#   } else {
#     Kissileff_procNoise_rmvnDat100_RMSE = rbind(Kissileff_procNoise_rmvnDat100_RMSE, Kissileff_procNoise_rmvn100_rmse_subset[[1]])
#     Kissileff_procNoise_rmvnParamRec_RMSE = rbind(Kissileff_procNoise_rmvnParamRec_RMSE, Kissileff_procNoise_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(Kissileff_procNoise_rmvnDat100_RMSE, 'Data/Kissileff_simDat_procNoise_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(Kissileff_procNoise_rmvnParamRec_RMSE, 'Data/Kissileff_simDat_procNoise_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process and measurement noise (bite size)
# for (l in 1:5){
#   Kissileff_procNoise_measureNoiseSize_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'BiteSize',
#                                                                 model_str = "Kissileff", datOnly = FALSE, write.dat = FALSE,
#                                                                 data_str = 'simDat_procNoise_measureNoise', rmse = 'both')
#   if (l == 1){
#     Kissileff_procNoise_measureNoiseSize_rmvnDat100_RMSE = Kissileff_procNoise_measureNoiseSize_rmvn100_rmse_subset[[1]]
#     Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE = Kissileff_procNoise_measureNoiseSize_rmvn100_rmse_subset[[2]]
#   } else {
#     Kissileff_procNoise_measureNoiseSize_rmvnDat100_RMSE = rbind(Kissileff_procNoise_measureNoiseSize_rmvnDat100_RMSE, Kissileff_procNoise_measureNoiseSize_rmvn100_rmse_subset[[1]])
#     Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE = rbind(Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE, Kissileff_procNoise_measureNoiseSize_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(Kissileff_procNoise_measureNoiseSize_rmvnDat100_RMSE, 'Data/Kissileff_simDat_procNoise_measureNoiseSize_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE, 'Data/Kissileff_simDat_procNoise_measureNoiseSize_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process and measurement noise (bite timing)
# for (l in 1:5){
#   Kissileff_procNoise_measureNoiseTime_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'BiteTiming',
#                                                                              model_str = "Kissileff", datOnly = FALSE, write.dat = FALSE,
#                                                                              data_str = 'simDat_procNoise_measureNoise', rmse = 'both')
#   if (l == 1){
#     Kissileff_procNoise_measureNoiseTime_rmvnDat100_RMSE = Kissileff_procNoise_measureNoiseTime_rmvn100_rmse_subset[[1]]
#     Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE = Kissileff_procNoise_measureNoiseTime_rmvn100_rmse_subset[[2]]
#   } else {
#     Kissileff_procNoise_measureNoiseTime_rmvnDat100_RMSE = rbind(Kissileff_procNoise_measureNoiseTime_rmvnDat100_RMSE, Kissileff_procNoise_measureNoiseTime_rmvn100_rmse_subset[[1]])
#     Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE = rbind(Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE, Kissileff_procNoise_measureNoiseTime_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(Kissileff_procNoise_measureNoiseTime_rmvnDat100_RMSE, 'Data/Kissileff_simDat_procNoise_measureNoiseTime_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE, 'Data/Kissileff_simDat_procNoise_measureNoiseTime_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

## process and measurement noise (bite size and timing)
# for (l in 1:5){
#   Kissileff_procNoise_measureNoise_rmvn100_rmse_subset = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'both',
#                                                                              model_str = "Kissileff", datOnly = FALSE, write.dat = FALSE,
#                                                                              data_str = 'simDat_procNoise_measureNoise', rmse = 'both')
#   if (l == 1){
#     Kissileff_procNoise_measureNoise_rmvnDat100_RMSE = Kissileff_procNoise_measureNoise_rmvn100_rmse_subset[[1]]
#     Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE = Kissileff_procNoise_measureNoise_rmvn100_rmse_subset[[2]]
#   } else {
#     Kissileff_procNoise_measureNoise_rmvnDat100_RMSE = rbind(Kissileff_procNoise_measureNoise_rmvnDat100_RMSE, Kissileff_procNoise_measureNoise_rmvn100_rmse_subset[[1]])
#     Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE = rbind(Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE, Kissileff_procNoise_measureNoise_rmvn100_rmse_subset[[2]])
#   }
# }
# 
# write.csv(Kissileff_procNoise_measureNoise_rmvnDat100_RMSE, 'Data/Kissileff_simDat_procNoise_measureNoise_rmvnDatRMSE100t5.csv', row.names = FALSE)
# write.csv(Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE, 'Data/Kissileff_simDat_procNoise_measureNoise_rmvnParamRecRMSE100t5.csv', row.names = FALSE)

#### Data Load ####
##load FPM data
FPM_rmvnParamRec_RMSE = read.csv('Data/FPM_simDat_rmvnParamRecRMSE100t5.csv')
FPM_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
FPM_rmvnParamRec_RMSE$noise = 'NoNoise'
FPM_rmvnDat_RMSE = read.csv('Data/FPM_simDat_rmvnDatRMSE100t5.csv')

FPM_procNoise_rmvnParamRec_RMSE = read.csv('Data/FPM_simDat_procNoise_rmvnParamRecRMSE100t5.csv')
FPM_procNoise_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
FPM_procNoise_rmvnParamRec_RMSE$noise = 'ProcNoise'
FPM_procNoise_rmvnDat_RMSE = read.csv('Data/FPM_simDat_procNoise_rmvnDatRMSE100t5.csv')

FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE = read.csv('Data/FPM_simDat_procNoise_measureNoiseSize_rmvnParamRecRMSE100t5.csv')
FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE$noise = 'ProcNoise_MeasureNoiseSize'
FPM_procNoise_measureNoiseSize_rmvnDat_RMSE = read.csv('Data/FPM_simDat_procNoise_measureNoiseSize_rmvnDatRMSE100t5.csv')

FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE = read.csv('Data/FPM_simDat_procNoise_measureNoiseTime_rmvnParamRecRMSE100t5.csv')
FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE$noise = 'ProcNoise_MeasureNoiseTime'
FPM_procNoise_measureNoiseTime_rmvnDat_RMSE = read.csv('Data/FPM_simDat_procNoise_measureNoiseTime_rmvnDatRMSE100t5.csv')

FPM_procNoise_measureNoise_rmvnParamRec_RMSE = read.csv('Data/FPM_simDat_procNoise_measureNoise_rmvnParamRecRMSE100t5.csv')
FPM_procNoise_measureNoise_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
FPM_procNoise_measureNoise_rmvnParamRec_RMSE$noise = 'ProcNoise_MeasureNoise'
FPM_procNoise_measureNoise_rmvnDat_RMSE = read.csv('Data/FPM_simDat_procNoise_measureNoise_rmvnDatRMSE100t5.csv')

##load Kissileff data
Kissileff_rmvnParamRec_RMSE = read.csv('Data/Kissileff_simDat_rmvnParamRecRMSE100t5.csv')
Kissileff_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
Kissileff_rmvnParamRec_RMSE$noise = 'NoNoise'
Kissileff_rmvnDat_RMSE = read.csv('Data/Kissileff_simDat_rmvnDatRMSE100t5.csv')

Kissileff_procNoise_rmvnParamRec_RMSE = read.csv('Data/Kissileff_simDat_procNoise_rmvnParamRecRMSE100t5.csv')
Kissileff_procNoise_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
Kissileff_procNoise_rmvnParamRec_RMSE$noise = 'ProcNoise'
Kissileff_procNoise_rmvnDat_RMSE = read.csv('Data/Kissileff_simDat_procNoise_rmvnDatRMSE100t5.csv')

Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE = read.csv('Data/Kissileff_simDat_procNoise_measureNoiseSize_rmvnParamRecRMSE100t5.csv')
Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE$noise = 'ProcNoise_MeasureNoiseSize'
Kissileff_procNoise_measureNoiseSize_rmvnDat_RMSE = read.csv('Data/Kissileff_simDat_procNoise_measureNoiseSize_rmvnDatRMSE100t5.csv')

Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE = read.csv('Data/Kissileff_simDat_procNoise_measureNoiseTime_rmvnParamRecRMSE100t5.csv')
Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE$noise = 'ProcNoise_MeasureNoiseTime'
Kissileff_procNoise_measureNoiseTime_rmvnDat_RMSE = read.csv('Data/Kissileff_simDat_procNoise_measureNoiseTime_rmvnDatRMSE100t5.csv')

Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE = read.csv('Data/Kissileff_simDat_procNoise_measureNoise_rmvnParamRecRMSE100t5.csv')
Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE$ID = seq(1, 500, 1)
Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE$noise = 'ProcNoise_MeasureNoise'
Kissileff_procNoise_measureNoise_rmvnDat_RMSE = read.csv('Data/Kissileff_simDat_procNoise_measureNoise_rmvnDatRMSE100t5.csv')

#### RMSE Distribution Graphs - FPM ####
#combine for graphing
FPM_allNoise_rmvnDat_RMSE = rbind(FPM_rmvnParamRec_RMSE, FPM_procNoise_rmvnParamRec_RMSE, FPM_procNoise_measureNoiseSize_rmvnParamRec_RMSE, FPM_procNoise_measureNoiseTime_rmvnParamRec_RMSE, FPM_procNoise_measureNoise_rmvnParamRec_RMSE)
FPM_allNoise_rmvnDat_RMSE$noise = factor(FPM_allNoise_rmvnDat_RMSE$noise)

#intake
FPM_rmvnParamRec_RMSEtiming_density = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity')

FPM_rmvnParamRec_RMSEtiming_density_facet = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

FPM_rmvnParamRec_RMSEtiming_hist = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity')

FPM_rmvnParamRec_RMSEtiming_hist_facet = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

#intake
FPM_rmvnParamRec_RMSEintake_density = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity')

FPM_rmvnParamRec_RMSEintake_density_facet = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

FPM_rmvnParamRec_RMSEintake_hist = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity')

FPM_rmvnParamRec_RMSEintake_hist_facet = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

#### RMSE Distribution Graphs - Kissileff ####

#combine for graphing
Kissileff_allNoise_rmvnDat_RMSE = rbind(Kissileff_rmvnParamRec_RMSE, Kissileff_procNoise_rmvnParamRec_RMSE, Kissileff_procNoise_measureNoiseSize_rmvnParamRec_RMSE, Kissileff_procNoise_measureNoiseTime_rmvnParamRec_RMSE, Kissileff_procNoise_measureNoise_rmvnParamRec_RMSE)
Kissileff_allNoise_rmvnDat_RMSE$noise = factor(Kissileff_allNoise_rmvnDat_RMSE$noise)

#timing
Kissileff_rmvnParamRec_RMSEtiming_density = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity')

Kissileff_rmvnParamRec_RMSEtiming_density_facet = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

Kissileff_rmvnParamRec_RMSEtiming_hist = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity')

Kissileff_rmvnParamRec_RMSEtiming_hist_facet = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_timing, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

#intake
Kissileff_rmvnParamRec_RMSEintake_density = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity')

Kissileff_rmvnParamRec_RMSEintake_density_facet = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_density(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

Kissileff_rmvnParamRec_RMSEintake_hist = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity')

Kissileff_rmvnParamRec_RMSEintake_hist_facet = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(x = rmse_intake, group = noise, fill = noise)) +
  geom_histogram(alpha=0.6, position = 'identity') + facet_wrap(~noise, scales = "free")

## FPM ##
#number NA
FPM_allNoise_rmvnDat_RMSE$timing_nNA_any = ifelse(FPM_allNoise_rmvnDat_RMSE$rmse_timing_nNA > 0, 'Y', 'N')
FPM_nNAtime_tab = xtabs(~timing_nNA_any + noise, data = FPM_allNoise_rmvnDat_RMSE)

FPM_allNoise_rmvnDat_RMSE$intake_nNA_any = ifelse(FPM_allNoise_rmvnDat_RMSE$rmse_timing_nNA > 0, 'Y', 'N')
FPM_nNAintake_tab = xtabs(~rmse_intake_nNA + noise, data = FPM_allNoise_rmvnDat_RMSE)

## Kissileff ##
#number NA for timing
Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any = ifelse(Kissileff_allNoise_rmvnDat_RMSE$rmse_timing_nNA > 0, 'Y', 'N')
Kissileff_nNAtime_tab = xtabs(~timing_nNA_any + noise, data = Kissileff_allNoise_rmvnDat_RMSE)
Kissileff_nNAtime_mean = means.function(Kissileff_allNoise_rmvnDat_RMSE[Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any == 'Y', ], 
                                        Kissileff_allNoise_rmvnDat_RMSE[Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any == 'Y', ]$rmse_timing_nNA, 
                                        Kissileff_allNoise_rmvnDat_RMSE[Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any == 'Y', ]$noise)
Kissileff_nNAtime_range = range.function(Kissileff_allNoise_rmvnDat_RMSE[Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any == 'Y', ], 
                                         Kissileff_allNoise_rmvnDat_RMSE[Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any == 'Y', ]$rmse_timing_nNA, 
                                         Kissileff_allNoise_rmvnDat_RMSE[Kissileff_allNoise_rmvnDat_RMSE$timing_nNA_any == 'Y', ]$noise)

Kissileff_allNoise_rmvnDat_RMSE$intake_nNA_any = ifelse(Kissileff_allNoise_rmvnDat_RMSE$rmse_intake_nNA > 0, 'Y', 'N')
Kissileff_nNAintake_tab = xtabs(~intake_nNA_any + noise, data = Kissileff_allNoise_rmvnDat_RMSE)

## Correlation with parameters
##FPM - timing
FPM_allNoise_rmvnDat_RMSE$r_dif = FPM_allNoise_rmvnDat_RMSE$initial_r - FPM_allNoise_rmvnDat_RMSE$r
FPM_allNoise_rmvnDat_RMSE$theta_dif = FPM_allNoise_rmvnDat_RMSE$initial_theta - FPM_allNoise_rmvnDat_RMSE$theta
FPM_allNoise_rmvnDat_RMSE$rmse_timing_relAvgBiteSize = FPM_allNoise_rmvnDat_RMSE$rmse_timing/(FPM_allNoise_rmvnDat_RMSE$Emax/FPM_allNoise_rmvnDat_RMSE$nBites)

FPM_corRMSE_rDif_timing = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = r_dif, x = rmse_timing)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and r - r recovered') + 
  scale_y_continuous(name='r - r recovered') +
  scale_x_continuous(name='RMSE timing') +
  facet_grid( ~ noise)

FPM_corRMSE.relAvgBite_rDif_timing = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = r_dif, x = rmse_timing_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and r - r recovered') + 
  scale_y_continuous(name='r - r recovered') +
  scale_x_continuous(name='RMSE timing/Average Bite Size') +
  facet_grid( ~ noise)

FPM_corRMSE_thetaDif_timing = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = theta_dif, x = rmse_timing)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and theta - theta recovered') + 
  scale_y_continuous(name='theta - theta recovered') +
  scale_x_continuous(name='RMSE timing') +
  facet_grid( ~ noise)

FPM_corRMSE.relAvgBite_thetaDif_timing = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = theta_dif, x = rmse_timing_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and theta - theta recovered') + 
  scale_y_continuous(name='theta - theta recovered') +
  scale_x_continuous(name='RMSE timing/Average Bite Size') +
  facet_grid( ~ noise)

##FPM - size
FPM_allNoise_rmvnDat_RMSE$rmse_intake_relAvgBiteSize = FPM_allNoise_rmvnDat_RMSE$rmse_intake/(FPM_allNoise_rmvnDat_RMSE$Emax/FPM_allNoise_rmvnDat_RMSE$nBites)

FPM_corRMSE_rDif_size = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = r_dif, x = rmse_intake)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and r - r recovered') + 
  scale_y_continuous(name='r - r recovered') +
  scale_x_continuous(name='RMSE size') +
  facet_grid( ~ noise)

FPM_corRMSE.relAvgBite_rDif_size = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = r_dif, x = rmse_intake_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and r - r recovered') + 
  scale_y_continuous(name='r - r recovered') +
  scale_x_continuous(name='RMSE size/Average Bite Size') +
  facet_grid( ~ noise)

FPM_corRMSE_thetaDif_size = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = theta_dif, x = rmse_intake)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and theta - theta recovered') + 
  scale_y_continuous(name='theta - theta recovered') +
  scale_x_continuous(name='RMSE size') +
  facet_grid( ~ noise)

FPM_corRMSE.relAvgBite_thetaDif_size = ggplot(FPM_allNoise_rmvnDat_RMSE, aes(y = theta_dif, x = rmse_intake_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and theta - theta recovered') + 
  scale_y_continuous(name='theta - theta recovered') +
  scale_x_continuous(name='RMSE size/Average Bite Size') +
  facet_grid( ~ noise)

##Kissileff - Timing
Kissileff_allNoise_rmvnDat_RMSE$int_dif = Kissileff_allNoise_rmvnDat_RMSE$initial_int - Kissileff_allNoise_rmvnDat_RMSE$int
Kissileff_allNoise_rmvnDat_RMSE$linear_dif = Kissileff_allNoise_rmvnDat_RMSE$initial_linear - Kissileff_allNoise_rmvnDat_RMSE$linear
Kissileff_allNoise_rmvnDat_RMSE$quad_dif = Kissileff_allNoise_rmvnDat_RMSE$initial_quad - Kissileff_allNoise_rmvnDat_RMSE$quad

Kissileff_allNoise_rmvnDat_RMSE$rmse_timing_relAvgBiteSize = Kissileff_allNoise_rmvnDat_RMSE$rmse_timing/(Kissileff_allNoise_rmvnDat_RMSE$Emax/Kissileff_allNoise_rmvnDat_RMSE$nBites)

Kissileff_corRMSE_intDif_timing = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = int_dif, x = rmse_timing)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and int - int recovered') + 
  scale_y_continuous(name='int - int recovered') +
  scale_x_continuous(name='RMSE timing') +
  facet_grid( ~ noise)

Kissileff_corRMSE.relAvgBite_intDif_timing = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = int_dif, x = rmse_timing_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and int - int recovered') + 
  scale_y_continuous(name='int - int recovered') +
  scale_x_continuous(name='RMSE timing/Average Bite Size') +
  facet_grid( ~ noise)

Kissileff_corRMSE_linearDif_timing = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = linear_dif, x = rmse_timing)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and linear - linear recovered') + 
  scale_y_continuous(name='linear - linear recovered') +
  scale_x_continuous(name='RMSE timing') +
  facet_grid( ~ noise)

Kissileff_corRMSE.relAvgBite_linearDif_timing = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = linear_dif, x = rmse_timing_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and linear - linear recovered') + 
  scale_y_continuous(name='linear - linear recovered') +
  scale_x_continuous(name='RMSE timing/Average Bite Size') +
  facet_grid( ~ noise)

Kissileff_corRMSE_quadDif_timing = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = quad_dif, x = rmse_timing)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and quad - quad recovered') + 
  scale_y_continuous(name='quad - quad recovered') +
  scale_x_continuous(name='RMSE timing') +
  facet_grid( ~ noise)

Kissileff_corRMSE.relAvgBite_quadDif_timing = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = quad_dif, x = rmse_timing_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE Timing and quad - quad recovered') + 
  scale_y_continuous(name='quad - quad recovered') +
  scale_x_continuous(name='RMSE timing/Average Bite Size') +
  facet_grid( ~ noise)

##Kissileff - size
Kissileff_allNoise_rmvnDat_RMSE$rmse_intake_relAvgBiteSize = Kissileff_allNoise_rmvnDat_RMSE$rmse_intake/(Kissileff_allNoise_rmvnDat_RMSE$Emax/Kissileff_allNoise_rmvnDat_RMSE$nBites)

Kissileff_corRMSE_intDif_size = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = int_dif, x = rmse_intake)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and int - int recovered') + 
  scale_y_continuous(name='int - int recovered') +
  scale_x_continuous(name='RMSE size') +
  facet_grid( ~ noise)

Kissileff_corRMSE.relAvgBite_intDif_size = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = int_dif, x = rmse_intake_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and int - int recovered') + 
  scale_y_continuous(name='int - int recovered') +
  scale_x_continuous(name='RMSE size/Average Bite Size') +
  facet_grid( ~ noise)

Kissileff_corRMSE_linearDif_size = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = linear_dif, x = rmse_intake)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and linear - linear recovered') + 
  scale_y_continuous(name='linear - linear recovered') +
  scale_x_continuous(name='RMSE size') +
  facet_grid( ~ noise)

Kissileff_corRMSE.relAvgBite_linearDif_size = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = linear_dif, x = rmse_intake_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and linear - linear recovered') + 
  scale_y_continuous(name='linear - linear recovered') +
  scale_x_continuous(name='RMSE size/Average Bite Size') +
  facet_grid( ~ noise)

Kissileff_corRMSE_quadDif_size = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = quad_dif, x = rmse_intake)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and quad - quad recovered') + 
  scale_y_continuous(name='quad - quad recovered') +
  scale_x_continuous(name='RMSE size') +
  facet_grid( ~ noise)

Kissileff_corRMSE.relAvgBite_quadDif_size = ggplot(Kissileff_allNoise_rmvnDat_RMSE, aes(y = quad_dif, x = rmse_intake_relAvgBiteSize)) +
  geom_smooth(method = lm, linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Correlation between RMSE size and quad - quad recovered') + 
  scale_y_continuous(name='quad - quad recovered') +
  scale_x_continuous(name='RMSE size/Average Bite Size') +
  facet_grid( ~ noise)
