############ Basic Data Load/Setup ############
library(reporttools)
library(xtable)
library(car)
library(ggplot2)
library(reshape2)
library(stats)
#library(rstudioapi)
library(psych)
library(faux)
#installed from devtools::install_github("debruine/faux")
library(RColorBrewer)
library(bitemodelr)


#### set up ####

#set piorking directory to location of script--not needed pihen called 
#through Rmarkdopin doc. Uncomment belopi if running locally/manually
# this.dir = getActiveDocumentContext()$path
# setwd(dirname(this.dir))

source('functions.R')
###################################################
####                                           
#### Create distributions for model parameters ####
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

#### Simulate Data and check characteristics
source('Fogel2017_SimDat.R')
SimDat_Fogel2017 = Fogel2017_simDat(500)

##Correlation Matrices

#original Fogel2017
Fogel2017_ReportedCorMat = matrix(c(NA, NA, NA, NA, NA, NA, 
                                    '-0.42*', NA, NA, NA, NA, NA,
                                    '-0.58*', '0.54*', NA, NA, NA, NA,
                                     '0.11*', '0.17*', '0.16*', NA, NA, NA,
                                     '0.54*', '-0.01', '0.02', '0.33', NA, NA, 
                                     '0.15*', '0.55*', '-0.25*', '-0.02', '-0.05', NA), byrow = TRUE, nrow = 6)

rownames(Fogel2017_ReportedCorMat) = c('nBites', 'BiteSize_g', 'BiteOralExposure_sec', 'ActiveMeal_pcent', 'TotalOralExposure_min', 'EatRate_g.min')
colnames(Fogel2017_ReportedCorMat) = c('nBites', 'BiteSize_g', 'BiteOralExposure_sec', 'ActiveMeal_pcent', 'TotalOralExposure_min', 'EatRate_g.min')


#orignial simulated
Fogel2017_corVars = SimDat_Fogel2017[c(2:5, 7:8, 6, 9:10)]
Fogel2017_corVarNames = names(SimDat_Fogel2017)[c(2:5, 7:8, 6, 9:10)]
Fogel2017_corMat = cor.matrix(Fogel2017_corVars, Fogel2017_corVarNames)
  
#rounded bites 
Fogel2017_BiteRound_corVars = SimDat_Fogel2017[c(11, 3:5, 13:14, 12, 15:16)]
Fogel2017_BiteRound_corMat = cor.matrix(Fogel2017_BiteRound_corVars, Fogel2017_corVarNames)

##Means(SD) for Fast vs Slow
SimDat_Fogel2017$EatRate_group = ifelse(SimDat_Fogel2017$EatRate_g.min < median(SimDat_Fogel2017$EatRate_g.min), 'Slow', 'Fast')
SimDat_Fogel2017$BiteRound_EatRate_group = ifelse(SimDat_Fogel2017$BiteRound_EatRate_g.min < median(SimDat_Fogel2017$BiteRound_EatRate_g.min), 'Slow', 'Fast')

#### Generate parameter distributions from Fogel_simDat ####

#get bite data set using random time sampling from a logistic curve
source('simBitesLogit.R')
SimBites_Fogel2017_list = t(mapply(simBitesLogit, mealdur = SimDat_Fogel2017$BiteRound_MealDur_min, nBites = SimDat_Fogel2017$BiteRound_nBites, Emax = SimDat_Fogel2017$BiteRound_TotalIntake_g, id = SimDat_Fogel2017$ID))

SimBites_Fogel2017 = data.frame(matrix(c(unlist(SimBites_Fogel2017_list)), byrow = FALSE, ncol = 4))
names(SimBites_Fogel2017) = c('ID', 'Bite', 'SampledTime', 'EstimatedCumulativeIntake')

#fit parameters to the bite datasets
SimBites_Fogel2017_params = IntakeModelParams(data = SimBites_Fogel2017, timeVar = 'SampledTime', intakeVar = 'EstimatedCumulativeIntake', fit_fn = FPM_Fit, idVar = 'ID')

SimDat_Fogel2017 = merge(SimDat_Fogel2017, SimBites_Fogel2017_params, by = 'ID')

#### Make Distribution Graphs ####
#histograms
SimDat_Fogel2017_thetaHist = ggplot(SimDat_Fogel2017, aes(x = theta, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_thetaHist_group = ggplot(SimDat_Fogel2017, aes(x = theta, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_rHist = ggplot(SimDat_Fogel2017, aes(x = r, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_rHist_group = ggplot(SimDat_Fogel2017, aes(x = r, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EmaxHist = ggplot(SimDat_Fogel2017, aes(x = BiteRound_TotalIntake_g, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_EmaxHist_group = ggplot(SimDat_Fogel2017, aes(x = BiteRound_TotalIntake_g, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EatRateHist = ggplot(SimDat_Fogel2017, aes(x = BiteRound_EatRate_g.min, group = EatRate_group)) +
  geom_histogram(bins = 20, position = 'identity')

SimDat_Fogel2017_EatRateHist_group = ggplot(SimDat_Fogel2017, aes(x = BiteRound_EatRate_g.min, group = EatRate_group)) +
  geom_histogram(bins = 20, aes(fill = EatRate_group), alpha=0.6, position = 'identity')

#Density
SimDat_Fogel2017_thetaDist_group = ggplot(SimDat_Fogel2017, aes(x = theta, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_rDist_group = ggplot(SimDat_Fogel2017, aes(x = r, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EmaxDist_group = ggplot(SimDat_Fogel2017, aes(x = BiteRound_TotalIntake_g, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')

SimDat_Fogel2017_EatRateDist_group = ggplot(SimDat_Fogel2017, aes(x = BiteRound_EatRate_g.min, group = EatRate_group)) +
  geom_density(aes(fill = EatRate_group), alpha=0.6, position = 'identity')
