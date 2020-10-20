# This script was written by Alaina Pearce in 2020 
# to generate distributions of parameters at the quartiles 
# and r very small r values
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

#source(param_plotGrid_setup.R)
Emax_mean = round(mean(SimDat_Fogel2017$TotalIntake_g), 2)

Kissileff_paramGrid_procNoiseDat = read.csv('Data/Kissileff_paramGrid_All_procNoiseDat.csv')
Kissileff_paramGrid_procNoise_30minDat = read.csv('Data/Kissileff_paramGrid_All_procNoise_30minDat.csv')

FPM_paramGrid_procNoiseDat = read.csv('Data/FPM_paramGrid_All_procNoiseDat.csv')
FPM_paramGrid_procNoise_30minDat = read.csv('Data/FPM_paramGrid_All_procNoise_30minDat.csv')

FPM_paramGrid_r_procNoiseDat = read.csv('Data/FPM_paramGrid_r_procNoiseDat.csv')
FPM_paramGrid_r_procNoise_30minDat = read.csv('Data/FPM_paramGrid_r_procNoise_30minDat.csv')

## Kissileff Model
Kissileff_paramGrid_procNoise_int1 = ggplot(Kissileff_paramGrid_procNoiseDat[Kissileff_paramGrid_procNoiseDat$int == "i = -43.43", ], aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Linear and Quadratic Density Distributions at Lowest Intercept of -43.43 (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_vline(xintercept = 30, linetype = 3, color = 'darkgrey') +
  facet_grid(linear ~ quad)

Kissileff_paramGrid_procNoise_intmean = ggplot(Kissileff_paramGrid_procNoiseDat[Kissileff_paramGrid_procNoiseDat$int == "i = 3.13", ], aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Linear and Quadratic Density Distributions at Lowest Intercept of -43.43 (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_vline(xintercept = 30, linetype = 3, color = 'darkgrey') +
  facet_grid(linear ~ quad)

Kissileff_paramGrid_procNoise_int5 = ggplot(Kissileff_paramGrid_procNoiseDat[Kissileff_paramGrid_procNoiseDat$int == "i = 41.2", ], aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Linear and Quadratic Density Distributions at Lowest Intercept of -43.43 (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_vline(xintercept = 30, linetype = 3, color = 'darkgrey') +
  facet_grid(linear ~ quad)

Kissileff_paramGrid_procNoise_30min_int1 = ggplot(Kissileff_paramGrid_procNoise_30minDat[Kissileff_paramGrid_procNoise_30minDat$int == "i = -43.43", ], aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions - 30 min Time Limit (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_hline(yintercept=Emax_mean, linetype = 3, color = 'darkgrey') +
  facet_grid(linear ~ quad)

Kissileff_paramGrid_procNoise_30min_intmean = ggplot(Kissileff_paramGrid_procNoise_30minDat[Kissileff_paramGrid_procNoise_30minDat$int == "i = 3.13", ], aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions - 30 min Time Limit (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_hline(yintercept=Emax_mean, linetype = 3, color = 'darkgrey') +
  facet_grid(linear ~ quad)

Kissileff_paramGrid_procNoise_30min_int5 = ggplot(Kissileff_paramGrid_procNoise_30minDat[Kissileff_paramGrid_procNoise_30minDat$int == "i = 41.2", ], aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions - 30 min Time Limit (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_hline(yintercept=Emax_mean, linetype = 3, color = 'darkgrey') +
  facet_grid(linear ~ quad)

## FPM

###all combinations
FPM_paramGrid_procNoise = ggplot(FPM_paramGrid_procNoiseDat, aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_vline(xintercept = 30, linetype = 3, color = 'darkgrey') +
  facet_grid(theta ~ r)

FPM_paramGrid_procNoise_30min = ggplot(FPM_paramGrid_procNoise_30minDat, aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions - 30 min Time Limit (55 bites, Emax = 107.62 grams)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_hline(yintercept=Emax_mean, linetype = 3, color = 'darkgrey') +
  facet_grid(theta ~ r)

###extreme r values - negative and small pos
FPM_paramGrid_r_procNoise = ggplot(FPM_paramGrid_r_procNoiseDat, aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions ') + 
  ggtitle('Emax calculated to work at lowest r and Bites Size set to 25th percentile: 14 bites, Emax = 17 grams') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_vline(xintercept = 30, linetype = 3, color = 'darkgrey') +
  facet_grid(theta ~ r)

FPM_paramGrid_r_procNoise_30min = ggplot(FPM_paramGrid_r_procNoise_30minDat, aes(x = EstimatedTime_procNoise, y = CumulativeGrams_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  geom_point(color = 'blue', shape = 3)+
  ggtitle('Intake Curves From Parameter Density Distributions ') + 
  ggtitle('Emax calculated to work at lowest r and Bites Size set to 25th percentile: 14 bites, Emax = 17 grams') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  geom_hline(yintercept=17, linetype = 3, color = 'darkgrey') +
  facet_grid(theta ~ r)