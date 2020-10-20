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
#### Create Density distributions for model parameters ####
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

#SimDat_Fogel2017  -  data is now part of the package

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
param_corPlot = ggpairs(SimDat_Fogel2017, columns = c(7, 10:11, 17:19), 
                        upper = list(continuous = wrap("cor", size=3)),
                        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.2))) +
  theme(legend.position = "bottom", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())


#### Extreme Paratemters ####
source('ParmDist_graphs.R')

# Draw random sample from multivariate distribution and test parameter recovery with log-likehood CI estimation ####
source('rmvnSample_ParamRec_CIgraphs.R')

##recovery tables

###FPM
FPM_rmvnParamRec$noise = 'No noise'
FPM_procNoise_rmvnParamRec$noise = 'Process noise'
FPM_procNoise_measureNoise_rmvnParamRec$noise = 'Measurement noise'

FPM_rmvnParamRec_all = rbind(FPM_rmvnParamRec, FPM_procNoise_rmvnParamRec, FPM_procNoise_measureNoise_rmvnParamRec)
FPM_rmvnParamRec_all$noise = factor(FPM_rmvnParamRec_all$noise)

FPM_rmvnParamRec_rfit_tab = xtabs(~r_fit + noise, data = FPM_rmvnParamRec_all)
FPM_rmvnParamRec_thetafit_tab = xtabs(~theta_fit + noise, data = FPM_rmvnParamRec_all)

###Kissileff
Kissileff_rmvnParamRec$noise = 'No noise'
Kissileff_procNoise_rmvnParamRec$noise = 'Process noise'
Kissileff_procNoise_measureNoise_rmvnParamRec$noise = 'Measurement noise'

Kissileff_rmvnParamRec_all = rbind(Kissileff_rmvnParamRec, Kissileff_procNoise_rmvnParamRec, Kissileff_procNoise_measureNoise_rmvnParamRec)
Kissileff_rmvnParamRec_all$noise = factor(Kissileff_rmvnParamRec_all$noise)

Kissileff_rmvnParamRec_intfit_tab = xtabs(~int_fit + noise, data = Kissileff_rmvnParamRec_all)
Kissileff_rmvnParamRec_linearfit_tab = xtabs(~linear_fit + noise, data = Kissileff_rmvnParamRec_all)
Kissileff_rmvnParamRec_quadfit_tab = xtabs(~quad_fit + noise, data = Kissileff_rmvnParamRec_all)

# Draw random sample from multivariate distribution and get distributions of RMSE ####
source('rmvnSample_RMSEdist_graphs.R')

# Model to model RMSE ####
