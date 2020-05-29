# This function was written by Alaina Pearce in 2020 
# to simulate data based on the mean and correlational
# structure of average meal microstructure in children
# reported in Fogel et al., 2017.
# 
#     Copyright (C) 2015 Alaina L Pearce
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


############ Basic Data Load/Setup ############
library(stats)
library(faux)
#installed from devtools::install_github("debruine/faux")

#### Fogel et al., 2017: A description of an ‘obesogenic’ eating style ####
#### that promotes higher energy intake and is associated with greater ####
#### adiposity in 4.5 year-old children: Results from the GUSTO cohort ####
#### Slow vs Faster Eaters (median = 6.41 g/min) ####
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

#this function will simulate a dataset of 
#n children
Fogel2017_simDat = function(n){
  #### Simuluate Data ####
  Fogel2017_r = matrix(c(1, -0.42, 0.54, 0.11, 
                         -0.42, 1, -0.01, 0.17, 
                         0.54, -0.01, 1, 0.33,
                         0.11, 0.17, 0.33, 1), byrow = TRUE, nrow = 4)
  
  Fogel2017_means_slow = c(57.7, 1.4, 15.1, 75, 175.3)
  Fogel2017_means_fast = c(68.4, 2.4, 15.2, 76.0, 306.8)
  Fogel2017_sds_slow = c(2.5, 0.1, 0.4, 1.0, 6.09)
  Fogel2017_sds_fast = c(2.5, 0.1, 0.4, 1.0, 9.9)
  
  SimDat_slow <- rnorm_multi(n = n/2, vars = 4,
                                       mu = Fogel2017_means_slow[1:4],
                                       sd = Fogel2017_sds_slow[1:4],
                                       r = Fogel2017_r, 
                                       varnames = c("nBites", "BiteSize", "OralExposure", "ActiveMeal"),
                                       empirical = TRUE)
  
  SimDat_slow$nBites = round(SimDat_slow$nBites)
  SimDat_slow$TotalGrams = SimDat_slow$BiteSize*SimDat_slow$nBites
  SimDat_slow$TotalOralExposure = (SimDat_slow$OralExposure*SimDat_slow$nBites)/60
  SimDat_slow$EatRate = SimDat_slow$TotalGrams/SimDat_slow$TotalOralExposure
  SimDat_slow$kcal = rnorm_pre(SimDat_slow$EatRate, mu = Fogel2017_means_slow[5], sd = Fogel2017_sds_slow[5], r = 0.61)
  
  SimDat_fast <- rnorm_multi(n = n/2, vars = 4,
                                       mu = Fogel2017_means_fast[1:4],
                                       sd = Fogel2017_sds_fast[1:4],
                                       r = Fogel2017_r, 
                                       varnames = c("nBites", "BiteSize", "OralExposure", "ActiveMeal"),
                                       empirical = TRUE)
  SimDat_fast$nBites = round(SimDat_fast$nBites)
  SimDat_fast$TotalGrams = SimDat_fast$BiteSize*SimDat_fast$nBites
  SimDat_fast$TotalOralExposure = (SimDat_fast$OralExposure*SimDat_fast$nBites)/60
  SimDat_fast$EatRate = SimDat_fast$TotalGrams/SimDat_fast$TotalOralExposure
  SimDat_fast$kcal = rnorm_pre(SimDat_fast$EatRate, mu = Fogel2017_means_fast[5], sd = Fogel2017_sds_fast[5], r = 0.61)
  
  
  SimDat = rbind(SimDat_slow, SimDat_fast)
  SimDat$MealDur = SimDat$OralExposure/(SimDat$ActiveMeal/100)
  
  SimDat_cormat = cor.matrix(SimDat[c(1:4, 6:9)], names(SimDat)[c(1:4, 6:9)])
  
  SimDat$id = seq(1, n, by=1)
  return(SimDat)
}
