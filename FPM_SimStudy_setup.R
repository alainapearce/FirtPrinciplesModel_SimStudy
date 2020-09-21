############ Basic Data Load/Setup ############
library(reporttools)
library(xtable)
library(car)
library(ggplot2)
library(reshape2)
library(stats)
#library(rstudioapi)
library(psych)
library(RColorBrewer)
library(bitemodelr)
library(GGally)
library(MASS)


#### set up ####

source('functions.R')

## ggplot theme
theme_set(theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
      panel.background = element_blank()))
###################################################
####                                           
#### Create Densityributions for model parameters ####
####
###################################################

##Create a simulated dataset
# Fogel et al., 2017: A description of an ‘obesogenic’ eating style
# that promotes higher energy intake and is associated with greater
# adiposity in 4.5 year-old children: Results from the GUSTO cohort
# Slow vs Faster Eaters (median = 6.41 g/min) ####
# N = 386
#table 2                      Slow - 192    Fast - 194    t         p
#  Bites (#)	                57.7 ± 2.5	  68.4 ± 2.5	  3.04	  0.003
#  Oral exposure per bite (s)	20.1 ± 0.9	  15.6 ± 0.5	  4.11	  < 0.0001
#  Bite size (g/bite)	        1.4 ± 0.1	    2.4 ± 0.1	    9.17	  < 0.0001
#  Chews per gram	            13.9 ± 0.5	  6.7 ± 0.1	    13.30	  < 0.0001
#  Sips (#)	                  8.6 ± 0.5	    8.2 ± 0.5	    0.54	  0.58
#  mealtime (%)	              75.0 ± 1.0	  76.0 ± 1.0	  0.56	  0.57
#  Total oral exposure (min)	15.1 ± 0.4	  15.2 ± 0.4	  0.08	  0.93
#  kcal	                      175.3 ± 6.09	306.8 ± 9.9	  11.28	  < 0.0001
#NOTE: oral exposure and mealtime(%) correlated (r=0.33)
#NOTE: bite size and mealtime(%) correlated (r=0.17)
#NOTE: Bites(#) and bite size correlated (r=-0.42)
#NOTE: Bites(#) and oral exposure correlated (r=0.54)

#NOTE: 7/9 foods had eating rates from 9.7-11.6, 1 was 5.8, and 1 was 15.1

#### Simulate Data and check characteristics ####

SimDat_Fogel2017 = read.csv('Data/ParamDat_Fogel2017.csv')

##Correlation Matrices
#original Fogel2017
Fogel2017_ReportedCorMat = matrix(c(NA, NA, NA, NA, NA,
                                    '0.15**', NA, NA, NA, NA,
                                    '0.54*', '-0.05', NA, NA, NA,
                                     '0.11*', '-0.02', '0.33***', NA, NA,
                                     '-0.58****', '-0.25***', '0.02', '0.16**', NA,
                                     NA, NA, NA, NA, NA,
                                     '-0.42***', '0.55*', '-0.01', '0.17**', '0.54***'), byrow = TRUE, nrow = 7)

rownames(Fogel2017_ReportedCorMat) = c('nBites', 'EatRate', 'TotalOE_min', 'ActiveMeal_pcent', 'BiteOE_sec', 'TotalIntake_g', 'BiteSize_g')
colnames(Fogel2017_ReportedCorMat) = c('nBites', 'EatRate', 'TotalOE_min', 'ActiveMeal_pcent', 'BiteOE_sec')


#original simulated

Fogel2017_corVars = SimDat_Fogel2017[c(2, 8, 6, 4:5, 7, 3)]
Fogel2017_corVarNames = names(SimDat_Fogel2017)[c(2, 8, 6, 4:5, 7, 3)]
Fogel2017_corMat = cor.matrix(Fogel2017_corVars, Fogel2017_corVarNames)
  

# Correlations between behaviors
SimDat_Fogel2017$EatRate_group = ifelse(SimDat_Fogel2017$EatRate < median(SimDat_Fogel2017$EatRate), 'Slow', 'Fast')

microstructure_corPlot = ggpairs(SimDat_Fogel2017, columns = c(2, 8, 6, 4:5, 7, 3), mapping = ggplot2::aes(colour=EatRate_group), legend = 1,
                                 upper = list(continuous = wrap("cor", size=3)),
                                 lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.2))) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.position = "bottom")


##Means(SD) for Fast vs Slow

#bites
bites_ttest = t.test(nBites~EatRate_group, data = SimDat_Fogel2017)
bites_se = se.function(SimDat_Fogel2017, SimDat_Fogel2017$nBites, SimDat_Fogel2017$EatRate_group)
bites_mean = means.function(SimDat_Fogel2017, SimDat_Fogel2017$nBites, SimDat_Fogel2017$EatRate_group)
bites_range = range.function(SimDat_Fogel2017, SimDat_Fogel2017$nBites, SimDat_Fogel2017$EatRate_group)

#bite size
bitesize_ttest = t.test(BiteSize_g~EatRate_group, data = SimDat_Fogel2017)
bitesize_se = se.function(SimDat_Fogel2017, SimDat_Fogel2017$BiteSize_g, SimDat_Fogel2017$EatRate_group)
bitesize_mean = means.function(SimDat_Fogel2017, SimDat_Fogel2017$BiteSize_g, SimDat_Fogel2017$EatRate_group)
bitesize_range = range.function(SimDat_Fogel2017, SimDat_Fogel2017$BiteSize_g, SimDat_Fogel2017$EatRate_group)

#oral exposure per bite
OE.bite_ttest = t.test(BiteOE_sec~EatRate_group, data = SimDat_Fogel2017)
OE.bite_se = se.function(SimDat_Fogel2017, SimDat_Fogel2017$BiteOE_sec, SimDat_Fogel2017$EatRate_group)
OE.bite_mean = means.function(SimDat_Fogel2017, SimDat_Fogel2017$BiteOE_sec, SimDat_Fogel2017$EatRate_group)
OE.bite_range = range.function(SimDat_Fogel2017, SimDat_Fogel2017$BiteOE_sec, SimDat_Fogel2017$EatRate_group)

#active mealtime
percActive_ttest = t.test(ActiveMeal_pcent~EatRate_group, data = SimDat_Fogel2017)
percActive_se = se.function(SimDat_Fogel2017, SimDat_Fogel2017$ActiveMeal_pcent, SimDat_Fogel2017$EatRate_group)
percActive_mean = means.function(SimDat_Fogel2017, SimDat_Fogel2017$ActiveMeal_pcent, SimDat_Fogel2017$EatRate_group)
percActive_range = range.function(SimDat_Fogel2017, SimDat_Fogel2017$ActiveMeal_pcent, SimDat_Fogel2017$EatRate_group)

#total oral exposure
OE.total_ttest = t.test(TotalOE_min~EatRate_group, data = SimDat_Fogel2017)
OE.total_se = se.function(SimDat_Fogel2017, SimDat_Fogel2017$TotalOE_min, SimDat_Fogel2017$EatRate_group)
OE.total_mean = means.function(SimDat_Fogel2017, SimDat_Fogel2017$TotalOE_min, SimDat_Fogel2017$EatRate_group)
OE.total_range = range.function(SimDat_Fogel2017, SimDat_Fogel2017$TotalOE_min, SimDat_Fogel2017$EatRate_group)

#total gram intake
total_g_ttest = t.test(TotalIntake_g~EatRate_group, data = SimDat_Fogel2017)
total_g_se = se.function(SimDat_Fogel2017, SimDat_Fogel2017$TotalIntake_g, SimDat_Fogel2017$EatRate_group)
total_g_mean = means.function(SimDat_Fogel2017, SimDat_Fogel2017$TotalIntake_g, SimDat_Fogel2017$EatRate_group)
total_g_range = range.function(SimDat_Fogel2017, SimDat_Fogel2017$TotalIntake_g, SimDat_Fogel2017$EatRate_group)

#### Densityribution Plots ####

## Eating Rate
SimDat_Fogel2017_EatRateHist = ggplot(SimDat_Fogel2017, aes(x = EatRate, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_EatRateHist_group = ggplot(SimDat_Fogel2017, aes(x = EatRate, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EatRateDensity_group = ggplot(SimDat_Fogel2017, aes(x = EatRate, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

## Kissileff Model

SimDat_Fogel2017_intHist = ggplot(SimDat_Fogel2017, aes(x = int, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_intHist_group = ggplot(SimDat_Fogel2017, aes(x = int, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_linearHist = ggplot(SimDat_Fogel2017, aes(x = linear, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_linearHist_group = ggplot(SimDat_Fogel2017, aes(x = linear, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_quadHist = ggplot(SimDat_Fogel2017, aes(x = quad, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_quadHist_group = ggplot(SimDat_Fogel2017, aes(x = quad, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_intDensity_group = ggplot(SimDat_Fogel2017, aes(x = int, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_linearDensity_group = ggplot(SimDat_Fogel2017, aes(x = linear, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_quadDensity_group = ggplot(SimDat_Fogel2017, aes(x = quad, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

## FPM Model
SimDat_Fogel2017_thetaHist = ggplot(SimDat_Fogel2017, aes(x = theta, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_thetaHist_group = ggplot(SimDat_Fogel2017, aes(x = theta, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_rHist = ggplot(SimDat_Fogel2017, aes(x = r, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_rHist_group = ggplot(SimDat_Fogel2017, aes(x = r, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EmaxHist = ggplot(SimDat_Fogel2017, aes(x = TotalIntake_g, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_EmaxHist_group = ggplot(SimDat_Fogel2017, aes(x = TotalIntake_g, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_thetaDensity_group = ggplot(SimDat_Fogel2017, aes(x = theta, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_rDensity_group = ggplot(SimDat_Fogel2017, aes(x = r, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EmaxDensity_group = ggplot(SimDat_Fogel2017, aes(x = TotalIntake_g, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

### Correlations between parameters ####
param_corPlot = ggpairs(SimDat_Fogel2017, columns = c(7, 10:11, 17:19), mapping = ggplot2::aes(colour=EatRate_group), legend = 1, 
                        upper = list(continuous = wrap("cor", size=3)),
                        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.2))) +
  theme(legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())


#### Extreme Paratemters ####
#no need to source after datasets are created and saved
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

# Draw random sample from multivariate distribution and test parameter recovery with only procNoise ####
source('rmvnSample_ParamRec.R')

#Test CI estimation with no error added - estimation per 100 randomly sampled microstructure sets
# only need to run if want to re-generate
# FPM_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = FALSE, model = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat')
# FPM_procNoise_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, model = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise')
# FPM_procNoise_measureNoise_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = TRUE, model = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise')

FPM_procNoise_rmvnParamRec = read.csv('Data/FPM_simDat_procNoise_rmvnParamRec100.csv')
FPM_procNoise_rmvnParamRec$ID = seq(1, 100, 1)
FPM_procNoise_rmvnDat = read.csv('Data/FPM_simDat_procNoise_rmvnDat100.csv')

FPM_rmvnParamRec = read.csv('Data/FPM_simDat_rmvnParamRec100.csv')
FPM_rmvnParamRec$ID = seq(1, 100, 1)
FPM_rmvnDat = read.csv('Data/FPM_simDat_rmvnDat100.csv')

##re-parameterization of r to e^r
# FPM_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = FALSE, model = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_rReParam')
# FPM_procNoise_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = TRUE, model = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise_rReParam')

FPM_rReParam_procNoise_rmvnParamRec = read.csv('Data/FPM_simDat_procNoise_rReParam_rmvnParamRec100.csv')
FPM_rReParam_procNoise_rmvnParamRec$ID = seq(1, 100, 1)
FPM_rReParam_procNoise_rmvnDat = read.csv('Data/FPM_simDat_procNoise_rReParam_rmvnDat100.csv')

FPM_rReParam_rmvnParamRec = read.csv('Data/FPM_simDat_rReParam_rmvnParamRec100.csv')
FPM_rReParam_rmvnParamRec$ID = seq(1, 100, 1)
FPM_rReParam_rmvnDat = read.csv('Data/FPM_simDat_rReParam_rmvnDat100.csv')

# Kissileff_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = FALSE, model = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat')
# Kissileff_procNoise_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = TRUE, model = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat_procNoise')

Kissileff_procNoise_rmvnParamRec = read.csv('Data/Kissileff_simDat_procNoise_rmvnParamRec100.csv')
Kissileff_procNoise_rmvnParamRec$ID = seq(1, 100, 1)
Kissileff_rmvnDat = read.csv('Data/Kissileff_simDat_rmvnDat100.csv')
Kissileff_rmvnParamRec = read.csv('Data/Kissileff_simDat_rmvnParamRec100.csv')
Kissileff_rmvnParamRec$ID = seq(1, 100, 1)
Kissileff_rmvnDat = read.csv('Data/Kissileff_simDat_rmvnDat100.csv')


#plot CIs in forestplot-like plot
FPM_rmvn100_CIplot_r <- ggplot(data=FPM_rReParam_rmvnParamRec, aes(y=ID, x=initial_r, xmin=l95CI_r, xmax=u95CI_r)) +
  geom_point(aes(col=r_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_r, xmax=u95CI_r, col=r_fit),height=0.5,cex=1) +
  xlab("r(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_procNoiseCIplot_r <- ggplot(data=FPM_rReParam_procNoise_rmvnParamRec, aes(y=ID, x=initial_r, xmin=l95CI_r, xmax=u95CI_r)) +
  geom_point(aes(col=r_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_r, xmax=u95CI_r, col=r_fit),height=0.5,cex=1) +
  xlab("r(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_CIplot_theta <- ggplot(data=FPM_rmvnParamRec, aes(y=ID, x=initial_theta, xmin=l95CI_theta, xmax=u95CI_theta)) +
  geom_point(aes(col=theta_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_theta, xmax=u95CI_theta, col=theta_fit),height=0.5,cex=1) +
  xlab("theta(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_procNoiseCIplot_theta <- ggplot(data=FPM_procNoise_rmvnParamRec, aes(y=ID, x=initial_theta, xmin=l95CI_theta, xmax=u95CI_theta)) +
  geom_point(aes(col=theta_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_theta, xmax=u95CI_theta, col=theta_fit),height=0.5,cex=1) +
  xlab("theta(95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_CIplot_int <- ggplot(data=Kissileff_rmvnParamRec, aes(y=ID, x=initial_int, xmin=l95CI_int, xmax=u95CI_int)) +
  geom_point(aes(col=int_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_int, xmax=u95CI_int, col=int_fit),height=0.5,cex=1) +
  xlab("Intercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_CIplot_int <- ggplot(data=Kissileff_procNoise_rmvnParamRec, aes(y=ID, x=initial_int, xmin=l95CI_int, xmax=u95CI_int)) +
  geom_point(aes(col=int_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_int, xmax=u95CI_int, col=int_fit),height=0.5,cex=1) +
  xlab("Intercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_CIplot_linear <- ggplot(data=Kissileff_rmvnParamRec, aes(y=ID, x=initial_linear, xmin=l95CI_linear, xmax=u95CI_linear)) +
  geom_point(aes(col=linear_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_linear, xmax=u95CI_linear, col=linear_fit),height=0.5,cex=1) +
  xlab("linearercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_CIplot_linear <- ggplot(data=Kissileff_procNoise_rmvnParamRec, aes(y=ID, x=initial_linear, xmin=l95CI_linear, xmax=u95CI_linear)) +
  geom_point(aes(col=linear_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_linear, xmax=u95CI_linear, col=linear_fit),height=0.5,cex=1) +
  xlab("linearercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_CIplot_quad <- ggplot(data=Kissileff_rmvnParamRec, aes(y=ID, x=initial_quad, xmin=l95CI_quad, xmax=u95CI_quad)) +
  geom_point(aes(col=quad_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_quad, xmax=u95CI_quad, col=quad_fit),height=0.5,cex=1) +
  xlab("quadercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_CIplot_quad <- ggplot(data=Kissileff_procNoise_rmvnParamRec, aes(y=ID, x=initial_quad, xmin=l95CI_quad, xmax=u95CI_quad)) +
  geom_point(aes(col=quad_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_quad, xmax=u95CI_quad, col=quad_fit),height=0.5,cex=1) +
  xlab("quadercept (95% CI)") + ylab("Random Sample Number")

#generate bite for 1
FPM_procNoise_rmvnParamRec_bothFit = FPM_procNoise_rmvnParamRec[FPM_procNoise_rmvnParamRec$r_fit == TRUE & FPM_procNoise_rmvnParamRec$theta_fit == TRUE, ]

FPM_procNoise_rmvnParamRec_true_Bites1 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$initial_theta[1], FPM_procNoise_rmvnParamRec_bothFit$initial_r[1]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_fit_Bites1 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$theta[1], FPM_procNoise_rmvnParamRec_bothFit$r[1]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_upper95CI_Bites1 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$u95CI_theta[1], FPM_procNoise_rmvnParamRec_bothFit$u95CI_r[1]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_lower95CI_Bites1 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$l95CI_theta[1], FPM_procNoise_rmvnParamRec_bothFit$l95CI_r[1]), time_fn = FPM_Time, procNoise = TRUE)


FPM_procNoise_rmvnParamRec_95CI_CumulativeIntake1 = ggplot(FPM_procNoise_rmvnParamRec_true_Bites1, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_rmvnParamRec_fit_Bites1, color = 'blue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_upper95CI_Bites1, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_lower95CI_Bites1, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI: heta = 9.13, r = -0.028') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')
  
FPM_procNoise_rmvnParamRec_true_Bites2 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[2], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[2], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$initial_theta[2], FPM_procNoise_rmvnParamRec_bothFit$initial_r[2]), time_fn = FPM_Time, procNoise = TRUE)
  
FPM_procNoise_rmvnParamRec_fit_Bites2 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[2], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[2], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$theta[2], FPM_procNoise_rmvnParamRec_bothFit$r[2]), time_fn = FPM_Time, procNoise = TRUE)
  
FPM_procNoise_rmvnParamRec_upper95CI_Bites2 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[2], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[2], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$u95CI_theta[2], FPM_procNoise_rmvnParamRec_bothFit$u95CI_r[2]), time_fn = FPM_Time, procNoise = TRUE)
  
FPM_procNoise_rmvnParamRec_lower95CI_Bites2 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[2], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[2], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$l95CI_theta[2], FPM_procNoise_rmvnParamRec_bothFit$l95CI_r[2]), time_fn = FPM_Time, procNoise = TRUE)
  
FPM_procNoise_rmvnParamRec_95CI_CumulativeIntake2 = ggplot(FPM_procNoise_rmvnParamRec_lower95CI_Bites2, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_rmvnParamRec_fit_Bites2, color = 'blue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_upper95CI_Bites2, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
    stat_smooth(data = FPM_procNoise_rmvnParamRec_lower95CI_Bites2, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI: theta = 28.37, r = 0.83') +
    scale_y_continuous(name='E(t)') +
    scale_x_continuous(name='Time (min)') 

FPM_procNoise_rmvnParamRec_true_Bites3 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[3], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[3], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$initial_theta[3], FPM_procNoise_rmvnParamRec_bothFit$initial_r[3]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_fit_Bites3 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[3], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[3], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$theta[3], FPM_procNoise_rmvnParamRec_bothFit$r[3]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_upper95CI_Bites3 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[3], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[3], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$u95CI_theta[3], FPM_procNoise_rmvnParamRec_bothFit$u95CI_r[3]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_lower95CI_Bites3 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[3], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[3], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$l95CI_theta[3], FPM_procNoise_rmvnParamRec_bothFit$l95CI_r[3]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_95CI_CumulativeIntake3 = ggplot(FPM_procNoise_rmvnParamRec_lower95CI_Bites3, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_rmvnParamRec_fit_Bites3, color = 'blue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_upper95CI_Bites3, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_lower95CI_Bites3, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI: theta = 16.47, r = 0.24') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') 

FPM_procNoise_rmvnParamRec_true_Bites4 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[4], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[4], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$initial_theta[4], FPM_procNoise_rmvnParamRec_bothFit$initial_r[4]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_fit_Bites4 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[4], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[4], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$theta[4], FPM_procNoise_rmvnParamRec_bothFit$r[4]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_upper95CI_Bites4 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[4], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[4], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$u95CI_theta[4], FPM_procNoise_rmvnParamRec_bothFit$u95CI_r[4]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_lower95CI_Bites4 = simBites(nBites = FPM_procNoise_rmvnParamRec_bothFit$nBites[4], Emax = FPM_procNoise_rmvnParamRec_bothFit$Emax[4], parameters = c(FPM_procNoise_rmvnParamRec_bothFit$l95CI_theta[4], FPM_procNoise_rmvnParamRec_bothFit$l95CI_r[4]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_rmvnParamRec_95CI_CumulativeIntake4 = ggplot(FPM_procNoise_rmvnParamRec_lower95CI_Bites4, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_rmvnParamRec_fit_Bites4, color = 'blue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_upper95CI_Bites4, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = FPM_procNoise_rmvnParamRec_lower95CI_Bites4, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI: theta = 5.41, r = 0.10') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') 
  
# Measurement Error ####
# Discretization of bite size
paramRec1_2catMean = ParamRecovery(nBites = FPM_rmvnDat$nBites[1], Emax = FPM_rmvnDat$TotalIntake_g[1], parameters = c(FPM_rmvnDat$theta[1], FPM_rmvnDat$r[1]), time_fn = FPM_Time, fit_fn = FPM_Fit, keepBites = TRUE, intake_fn = FPM_Intake, paramCI = list(r = TRUE, theta = TRUE), nSims = 100, simVar = 'biteSize', simValue = FPM_rmvnDat$TotalIntake_g[1]/FPM_rmvnDat$nBites[1])

paramRec1_2catMean_params = paramRec1_2catMean$paramDat
paramRec1_2catMean_params$CIfit_r = ifelse(round(paramRec1_2catMean_params$initial_r, 5) <= round(paramRec1_2catMean_params$u95CI_r, 5) & round(paramRec1_2catMean_params$initial_r, 5) >= round(paramRec1_2catMean_params$l95CI_r, 5), 'Y', 'N')


#get parameter recovery with bite timing error





