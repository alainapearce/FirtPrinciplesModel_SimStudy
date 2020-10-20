# This function was written by Alaina
# Pearce in 2020 to estimate intake using the incorrect formula in Thomas
# et al., 2017
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

#Estimate times using the correct version of the
#Thomas et al., 2017 First-Principles dynamic model

#input arguments
#E_t:cumulative intake at a time point - single value
#params: values for theta and r
#Emax: total intake in grams

FPMincorrect_Intake = function(intake, parameters, Emax){
  #params = c(theta, r)

  #since it is a logistic function, theoretically E_t will never be Emax. If
  #E_t = Emax, use 99% of Emax to get estimate for last timepoint
  if(intake == Emax){
    intake = intake*.9999
  }

  (Emax/(Emax*parameters[2]+parameters[1]))*log((Emax*((intake*parameters[2])+1))/(Emax - intake))
}

#Estimate times using the correct version of the
#Thomas et al., 2017 First-Principles dynamic model

#input arguments
#E_t:cumulative intake at a time point - single value
#params: values for theta and r
#Emax: total intake in grams

FPMincorrect_Time = function(intake, parameters, Emax){
  #params = c(theta, r)

  #since it is a logistic function, theoretically E_t will never be Emax. If
  #E_t = Emax, use 99% of Emax to get estimate for last timepoint
  if(intake == Emax){
    intake = intake*.9999
  }

  (Emax/(Emax*parameters[2]+parameters[1]))*log((Emax*((intake*parameters[2])+1))/(Emax - intake))
}
