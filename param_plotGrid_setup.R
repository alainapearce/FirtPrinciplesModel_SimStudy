# This script was written by Alaina Pearce in 2020 
# to generate the bite data used to create the param_plotGrid
# figures that combine parameters at the lowest, 25th percentile,
# mean, 75th percentile, and highest values from the distributions of 
# parameters based on the Fogel2017 simulated dataset
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

#if running outside of a script that loads the data, need to load data here
#SimDat_Fogel2017 = read.csv('Data/ParamDat_Fogel2017.csv')

#### Extreme Paratemters ###

Emax_quant = quantile(SimDat_Fogel2017$TotalIntake_g)
Emax_mean = round(mean(SimDat_Fogel2017$TotalIntake_g), 2)
nBites_quant = round(quantile(SimDat_Fogel2017$nBites))
nBites_mean = round(mean(SimDat_Fogel2017$nBites))
BiteSize_quant = quantile(SimDat_Fogel2017$BiteSize_g)

## Kissileff Model
linear_quants = quantile(SimDat_Fogel2017$linear)
linear_mean = mean(SimDat_Fogel2017$linear)
quad_quants = quantile(SimDat_Fogel2017$quad)
quad_mean = mean(SimDat_Fogel2017$quad)
int_quant = quantile(SimDat_Fogel2017$int)
int_mean = mean(SimDat_Fogel2017$int)

# not all values can be together as the sqrt term in the equation to derived
# time must be positive for all values
# > linear_quants
# 0%       25%       50%       75%      100% 
# -3.628937  8.353430 13.576919 20.853112 66.185516 
# > quad_quants
# 0%        25%        50%        75%       100% 
# -7.6062035 -0.8673148 -0.4568702 -0.2564256  2.3671355 
# > int_quant
# 0%         25%         50%         75%        100% 
# -43.4252209   0.2196949   2.3052960   6.2698686  41.2033124 
# 
# q1 Int can do: linear q1-5, quad q5
# linear q3-5, quad q4
# linear q5, quad q3
# linear q5, quad q2
# 
# q2 Int can do: linear q1-5, quad q5
# linear q3-5, quad q4
# linear q4-5, quad q3
# linear q4-5, quad q2
# linear q5, quad q1
# 
# q3 Int can do: linear q2-5, quad q5
# linear q3-5, quad q4
# linear q4-5, quad q3
# linear q4-5, quad q2
# linear q5, quad q1
# 
# q4 Int can do: linear q2-5, quad q5
# linear q3-5, quad q4
# linear q4-5, quad q3
# linear q4-5, quad q2
# linear q5, quad q1
# 
# q5 Int can do: linear q4-5, quad q5
# linear q2-5, quad q4
# linear q3-5, quad q3
# linear q4-5, quad q2
# linear q5, quad q1
#                
## All Thetas, highest R
counter = 0
mod_run = 'N'
for(i in 1:5){
  if(i == 3){
    int = int_mean
  } else {
    int = int_quant[[i]]
  }
  
  for(q in 1:5){
    if(q == 3){
      quad = quad_mean 
    } else {
      quad = quad_quants[[q]]
    }
    
    for(l in 1:5){
      if(l == 3){
        linear = linear_mean 
      } else {
        linear = linear_quants[[l]]
      }
      
      if(q == 5){
        if(i < 3){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } else if(i >= 3 & i < 5 & l >= 2){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } else if(i == 5 & l >= 4){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } 
      } else if(q == 4){
        if(i <= 4 & l >= 3){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } else if(i == 5 & l >= 2){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        }
      } else if(q == 3){
        if(i == 1 & l == 5){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } else if(i == 5 & l >= 3){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } else if(i >= 2 & l >= 4){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } 
      } else if(q == 2){
        if(i >=2 & l >= 4){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } else if(i == 1 & l == 5){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        } 
      } else if(q == 1){
        if(i > 1 & l == 5){
          procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE)
          procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(int, linear, quad), time_fn = Kissileff_Time, procNoise = TRUE, maxDur = 30)
          mod_run = 'Y'
        }
      }
      
      if(mod_run == 'Y'){
        procNoise$int = paste0("i = ", round(int, 2))
        procNoise_30min$int = paste0("i = ", round(int, 2))
        
        procNoise$quad = paste0("q = ", round(quad, 2))
        procNoise_30min$quad = paste0("q = ", round(quad, 2))
        
        procNoise$linear = paste0("l = ", round(linear, 2))
        procNoise_30min$linear = paste0("l = ", round(linear, 2))
        
        counter = counter + 1
        if(counter == 1){
          Kissileff_paramGrid_procNoiseDat = procNoise
          Kissileff_paramGrid_procNoise_30minDat = procNoise_30min
        } else {
          Kissileff_paramGrid_procNoiseDat = rbind(Kissileff_paramGrid_procNoiseDat, procNoise)
          Kissileff_paramGrid_procNoise_30minDat = rbind(Kissileff_paramGrid_procNoise_30minDat, procNoise_30min)
        }
        
        #reset
        mod_run = 'N'
      }
    }
  }
}


Kissileff_paramGrid_procNoiseDat$int = factor(Kissileff_paramGrid_procNoiseDat$int, levels = c("i = -43.43", "i = 0.22", "i = 3.13", "i = 6.27", "i = 41.2"))
Kissileff_paramGrid_procNoiseDat$linear = factor(Kissileff_paramGrid_procNoiseDat$linear, levels = c("l = -3.63", "l = 8.35", "l = 15.83", "l = 20.85", "l = 66.19"))
Kissileff_paramGrid_procNoiseDat$quad = factor(Kissileff_paramGrid_procNoiseDat$quad, levels = c("q = -7.61", "q = -0.87", "q = -0.73", "q = -0.26", "q = 2.37"))


Kissileff_paramGrid_procNoise_30minDat$int = factor(Kissileff_paramGrid_procNoise_30minDat$int, levels = c("i = -43.43", "i = 0.22", "i = 3.13", "i = 6.27", "i = 41.2"))
Kissileff_paramGrid_procNoise_30minDat$linear = factor(Kissileff_paramGrid_procNoise_30minDat$linear, levels = c("l = -3.63", "l = 8.35", "l = 15.83", "l = 20.85", "l = 66.19"))
Kissileff_paramGrid_procNoise_30minDat$quad = factor(Kissileff_paramGrid_procNoise_30minDat$quad, levels = c("q = -7.61", "q = -0.87", "q = -0.73", "q = -0.26", "q = 2.37"))

write.csv(Kissileff_paramGrid_procNoiseDat, 'Data/Kissileff_paramGrid_All_procNoiseDat.csv', row.names = FALSE)
write.csv(Kissileff_paramGrid_procNoise_30minDat, 'Data/Kissileff_paramGrid_All_procNoise_30minDat.csv', row.names = FALSE)

## FPM

theta_quants = quantile(SimDat_Fogel2017$theta)
theta_mean = mean(SimDat_Fogel2017$theta)
r_quants = quantile(SimDat_Fogel2017$r)
r_mean = mean(SimDat_Fogel2017$r)


## All Thetas, highest R
counter = 0
for(t in 1:5){
  for(r in 1:5){
    if(t == 3 & r == 3){
      procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_mean, r_mean), time_fn = FPM_Time, procNoise = TRUE)
      procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_mean, r_mean), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
    } else if(t == 3 & r != 3){
      if(r == 1){
        procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], 0), time_fn = FPM_Time, procNoise = TRUE)
        procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], 0), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
      } else {
        procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_mean, r_quants[[r]]), time_fn = FPM_Time, procNoise = TRUE)
        procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_mean, r_quants[[r]]), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
      }
    } else if(t != 3 & r == 3){
      procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], r_mean), time_fn = FPM_Time, procNoise = TRUE)
      procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], r_mean), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
    } else if(t != 3 & r != 3){
      if(r == 1){
        procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], 0), time_fn = FPM_Time, procNoise = TRUE)
        procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], 0), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
      } else {
        procNoise = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], r_quants[[r]]), time_fn = FPM_Time, procNoise = TRUE)
        procNoise_30min = simBites(nBites = nBites_mean, Emax = Emax_mean, parameters = c(theta_quants[[t]], r_quants[[r]]), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
      }
    }
    
    if (r == 1){
      procNoise$r = paste0("r = ", 0)
      procNoise_30min$r = paste0("r = ", 0)
    } else {
      procNoise$r = paste0("r = ", round(r_quants[[r]], 2))
      procNoise_30min$r = paste0("r = ", round(r_quants[[r]], 2))
    }
    
    procNoise$theta = paste0("t = ", round(theta_quants[[t]], 2))
    procNoise_30min$theta = paste0("t = ", round(theta_quants[[t]], 2))
    
    counter = counter + 1
    if(counter == 1){
      FPM_paramGrid_procNoiseDat = procNoise
      FPM_paramGrid_procNoise_30minDat = procNoise_30min
    } else {
      FPM_paramGrid_procNoiseDat = rbind(FPM_paramGrid_procNoiseDat, procNoise)
      FPM_paramGrid_procNoise_30minDat = rbind(FPM_paramGrid_procNoise_30minDat, procNoise_30min)
    }
  }
}

FPM_paramGrid_procNoiseDat$r = factor(FPM_paramGrid_procNoiseDat$r, levels = c("r = 0", "r = 0.06", "r = 0.14", "r = 0.25", "r = 1.87"))
FPM_paramGrid_procNoise_30minDat$r = factor(FPM_paramGrid_procNoise_30minDat$r, levels = c("r = 0", "r = 0.06", "r = 0.14", "r = 0.25", "r = 1.87"))

FPM_paramGrid_procNoiseDat$theta = factor(FPM_paramGrid_procNoiseDat$theta, levels = c("t = 0.46", "t = 7.74", "t = 13.65", "t = 23.03", "t = 109.4"))
FPM_paramGrid_procNoise_30minDat$theta = factor(FPM_paramGrid_procNoise_30minDat$theta, levels = c("t = 0.46", "t = 7.74", "t = 13.65", "t = 23.03", "t = 109.4"))

write.csv(FPM_paramGrid_procNoiseDat, 'Data/FPM_paramGrid_All_procNoiseDat.csv', row.names = FALSE)
write.csv(FPM_paramGrid_procNoise_30minDat, 'Data/FPM_paramGrid_All_procNoise_30minDat.csv', row.names = FALSE)

##Small and negative R
###theta set to 70 (Very high) and intake to the 25th percentile of intake; nBites calculated so bite size is at 25th percentile (1.25 g/bite)
###theta set to 75th percentile (23.0276) and intake to the 25th percentile of intake to 17g (very low); nBites calculated so bite size is at 25th (1.25 g/bite) percentile
r_near0set = c(r_quants[[1]], -.5, -1*r_quants[[4]], -1*r_quants[[3]], -1*r_quants[[2]], 0, r_quants[[2]], r_quants[[3]], r_quants[[4]])
counter = 0
for(r in 1:length(r_near0set)){
  for(t in 1:2){
    if(t == 1){
      procNoise = simBites(nBites = round(17/BiteSize_quant[[2]]), Emax = 17, parameters = c(70, r_near0set[[r]]), time_fn = FPM_Time, procNoise = TRUE)
      procNoise_30min = simBites(nBites = round(17/BiteSize_quant[[2]]), Emax = 17, parameters = c(70, r_near0set[[r]]), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
      
      procNoise$theta = "t = 70"
      procNoise_30min$theta = "t = 70"
      
    } else if(t == 2){
      procNoise = simBites(nBites = round(17/BiteSize_quant[[2]]), Emax = 17, parameters = c(theta_quants[[4]], r_near0set[[r]]), time_fn = FPM_Time, procNoise = TRUE)
      procNoise_30min = simBites(nBites = round(17/BiteSize_quant[[2]]), Emax = 17, parameters = c(theta_quants[[4]], r_near0set[[r]]), time_fn = FPM_Time, procNoise = TRUE, maxDur = 30)
      
      procNoise$theta = paste0("t = ", round(theta_quants[[4]], 2))
      procNoise_30min$theta = paste0("t = ", round(theta_quants[[4]], 2))
    }
    
    procNoise$r = paste0("r = ", round(r_near0set[[r]], 2))
    procNoise_30min$r = paste0("r = ", round(r_near0set[[r]], 2))
    
    counter = counter + 1
    if(counter == 1){
      FPM_paramGrid_r_procNoiseDat = procNoise
      FPM_paramGrid_r_procNoise_30minDat = procNoise_30min
    } else {
      FPM_paramGrid_r_procNoiseDat = rbind(FPM_paramGrid_r_procNoiseDat, procNoise)
      FPM_paramGrid_r_procNoise_30minDat = rbind(FPM_paramGrid_r_procNoise_30minDat, procNoise_30min)
    }
  }
}

FPM_paramGrid_r_procNoiseDat$r = factor(FPM_paramGrid_r_procNoiseDat$r, levels = c("r = -1.29", "r = -0.5", "r = -0.25", "r = -0.14", "r = -0.06", "r = 0", "r = 0.06", "r = 0.14", "r = 0.25"))
FPM_paramGrid_r_procNoise_30minDat$r = factor(FPM_paramGrid_r_procNoise_30minDat$r, levels = c("r = -1.29", "r = -0.5", "r = -0.25", "r = -0.14", "r = -0.06", "r = 0", "r = 0.06", "r = 0.14", "r = 0.25"))

FPM_paramGrid_r_procNoiseDat$theta = factor(FPM_paramGrid_r_procNoiseDat$theta, levels = c("t = 23.03", "t = 70"))
FPM_paramGrid_r_procNoise_30minDat$theta = factor(FPM_paramGrid_r_procNoise_30minDat$theta, levels = c("t = 23.03", "t = 70"))

write.csv(FPM_paramGrid_r_procNoiseDat, 'Data/FPM_paramGrid_r_procNoiseDat.csv', row.names = FALSE)
write.csv(FPM_paramGrid_r_procNoise_30minDat, 'Data/FPM_paramGrid_r_procNoise_30minDat.csv', row.names = FALSE)