# This script was written by Alaina Pearce in 2020 
# to create:
# 1) a databased based on the mean and correlational
# structure of average meal microstructure in children
# reported in Fogel et al., 2017. 
# 2) generate and save parameter distributes from generated
# database
# 
#     Copyright (C) 2012 Alaina L Pearce
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


###now a script in raw-data called DataGen_SimDat_Fogel2017

############ Basic Data Load/Setup ############
library(bitemodelr)
source('Fogel2017_SimDat.R')
#only need to do once
SimDat_Fogel2017 = Fogel2017_simDat(500)
write.csv(SimDat_Fogel2017, 'Data/SimDat_Fogel2017.csv', row.names = FALSE)

#### Generate bite data ####
#sample bite timing from a logistic curve and use average bite size to get cumulative intake
source('simBitesLogit.R')

SimBites_Fogel2017_list = t(mapply(simBitesLogit, mealdur = SimDat_Fogel2017$MealDur_min, nBites = SimDat_Fogel2017$nBites, Emax = SimDat_Fogel2017$TotalIntake_g, id = SimDat_Fogel2017$ID))

SimBites_Fogel2017 = data.frame(matrix(c(unlist(SimBites_Fogel2017_list)), byrow = FALSE, ncol = 4))
names(SimBites_Fogel2017) = c('ID', 'Bite', 'SampledTime', 'EstimatedCumulativeIntake')

write.csv(SimBites_Fogel2017, 'Data/SimBites_Fogel2017.csv', row.names = FALSE)

#### FPM Model ####
#fit parameters to the bite datasets
FPM_SimBites_Fogel2017_params = IntakeModelParams(data = SimBites_Fogel2017, timeVar = 'SampledTime', intakeVar = 'EstimatedCumulativeIntake', fit_fn = FPM_Fit, idVar = 'ID', CI = FALSE)

#### Kissileff Model ####
#fit parameters to the bite datasets
Kissileff_SimBites_Fogel2017_params = IntakeModelParams(data = SimBites_Fogel2017, timeVar = 'SampledTime', intakeVar = 'EstimatedCumulativeIntake', fit_fn = Kissileff_Fit, idVar = 'ID', CI = FALSE)

#### Add parameters to data ####
SimDat_Fogel2017 = merge(SimDat_Fogel2017, FPM_SimBites_Fogel2017_params, by = 'ID')
names(SimDat_Fogel2017)[12:16] = c('FPM_value', 'FPM_counts', 'FPM_counts_gradiant', 'FPM_convergence', 'FPM_method')

SimDat_Fogel2017 = merge(SimDat_Fogel2017, Kissileff_SimBites_Fogel2017_params, by = 'ID')
names(SimDat_Fogel2017)[20:24] = c('Kissileff_value', 'Kissileff_counts', 'Kissileff_counts_gradiant', 'Kissileff_convergence', 'Kissileff_method')
write.csv(SimDat_Fogel2017, 'Data/ParamDat_Fogel2017.csv', row.names = FALSE)
