# This script was written by Alaina Pearce in 2020 
# to generate the bite data using random sampling from
# the multivariate normal distribution generated from 
# Fogel et al., 2017 means and a scaled covariance structure. Additionally,
# parameter recovery and confidence interval estimation using the 
# log-likelihood space is used and graphed using forest-like plots.
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
#### Generate Data ####
#only needs to be done once - can just load data after that
# Test param recovery CI - estimation per 100 randomly sampled microstructure sets ####
# only need to run if want to re-generate

# FPM_scaled0.5_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = FALSE, model_str = "FPM",
#                                      datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both',
#                                      data_str = 'simDat', scaleFactor = 0.5)
# FPM_procNoise_scaled0.5_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, model_str = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise', scaleFactor = 0.5)
# FPM_procNoise_measureNoise_scaled0.5_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'both', model_str = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise_measureNoise', scaleFactor = 0.5)

# Kissileff_scaled0.5_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = FALSE, model_str = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat', scaleFactor = 0.5)
# Kissileff_procNoise_scaled0.5_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = TRUE, model_str = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat_procNoise', scaleFactor = 0.5)
# Kissileff_procNoise_measureNoise_scaled0.5_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = TRUE, measureNoise = 'both', model_str = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat_procNoise_measureNoise', scaleFactor = 0.5)


# FPM_scaled0.25_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = FALSE, model_str = "FPM",
#                                      datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both',
#                                      data_str = 'simDat', scaleFactor = 0.25)
# FPM_procNoise_scaled0.25_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, model_str = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise', scaleFactor = 0.25)
# FPM_procNoise_measureNoise_scaled0.25_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, procNoise = TRUE, measureNoise = 'both', model_str = "FPM", datOnly = FALSE, paramCI = c('r', 'theta'), bound = 'both', data_str = 'simDat_procNoise_measureNoise', scaleFactor = 0.25)

# Kissileff_scaled0.25_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = FALSE, model_str = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat', scaleFactor = 0.25)
# Kissileff_procNoise_scaled0.25_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = TRUE, model_str = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat_procNoise', scaleFactor = 0.25)
# Kissileff_procNoise_measureNoise_scaled0.25_rmvn100_CI = rmvnSample_ParamRec(nSample = 100, nSim = 1, procNoise = TRUE, measureNoise = 'both', model_str = "Kissileff", datOnly = FALSE, paramCI = c('int', 'linear', 'quad'), bound = 'both', data_str = 'simDat_procNoise_measureNoise', scaleFactor = 0.25)


#### Data Load ####
##load FPM data
FPM_scaled_rmvnParamRec = read.csv('Data/FPM_simDat_scaled_rmvnParamRecCI100.csv')
FPM_scaled_rmvnParamRec$ID = seq(1, 100, 1)
FPM_rmvnDat = read.csv('Data/FPM_simDat_rmvnDatCI100.csv')

FPM_procNoise_scaled_rmvnParamRec = read.csv('Data/FPM_simDat_procNoise_scaled_rmvnParamRecCI100.csv')
FPM_procNoise_scaled_rmvnParamRec$ID = seq(1, 100, 1)
FPM_procNoise_rmvnDat = read.csv('Data/FPM_simDat_procNoise_rmvnDatCI100.csv')

FPM_procNoise_measureNoise_scaled_rmvnParamRec = read.csv('Data/FPM_simDat_procNoise_measureNoise_scaled_rmvnParamRecCI100.csv')
FPM_procNoise_measureNoise_scaled_rmvnParamRec$ID = seq(1, 100, 1)
FPM_procNoise_measureNoise_rmvnDat = read.csv('Data/FPM_simDat_procNoise_measureNoise_rmvnDatCI100.csv')

## reParam of R
FPM_procNoise_rReParam_scaled_rmvnParamRec = read.csv('Data/FPM_simDat_rReParam_scaled_rmvnParamRecCI100.csv')
FPM_procNoise_rReParam_scaled_rmvnParamRec$ID = seq(1, 100, 1)
FPM_procNoise_rReParam_rmvnDat = read.csv('Data/FPM_simDat_rReParam_rmvnDatCI100.csv')

##load Kissileff data
Kissileff_scaled_rmvnParamRec = read.csv('Data/Kissileff_simDat_scaled_rmvnParamRecCI100.csv')
Kissileff_scaled_rmvnParamRec$ID = seq(1, 100, 1)
Kissileff_rmvnDat = read.csv('Data/Kissileff_simDat_rmvnDatCI100.csv')

Kissileff_procNoise_scaled_rmvnParamRec = read.csv('Data/Kissileff_simDat_procNoise_scaled_rmvnParamRecCI100.csv')
Kissileff_procNoise_scaled_rmvnParamRec$ID = seq(1, 100, 1)
Kissileff_procNoise_rmvnDat = read.csv('Data/Kissileff_simDat_procNoise_rmvnDatCI100.csv')

Kissileff_procNoise_measureNoise_scaled_rmvnParamRec = read.csv('Data/Kissileff_simDat_procNoise_measureNoise_scaled_rmvnParamRecCI100.csv')
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$ID = seq(1, 100, 1)
Kissileff_procNoise_measureNoise_rmvnDat = read.csv('Data/Kissileff_simDat_procNoise_measureNoise_rmvnDatCI100.csv')

#### Parameter Recovery CI graphs - FPM ####
#plot CIs in forestplot-like plot 
FPM_scaled_rmvn100_CIplot_r <- ggplot(data=FPM_scaled_rmvnParamRec, aes(y=ID, x=initial_r, xmin=l95CI_r, xmax=u95CI_r)) +
  geom_point(aes(col=r_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_r, xmax=u95CI_r, col=r_fit),height=0.5,cex=1) +
  xlab("r(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_procNoise_CIplot_r <- ggplot(data=FPM_procNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_r, xmin=l95CI_r, xmax=u95CI_r)) +
  geom_point(aes(col=r_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_r, xmax=u95CI_r, col=r_fit),height=0.5,cex=1) +
  xlab("r(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_procNoise_measureNoise_CIplot_r <- ggplot(data=FPM_procNoise_measureNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_r, xmin=l95CI_r, xmax=u95CI_r)) +
  geom_point(aes(col=r_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_r, xmax=u95CI_r, col=r_fit),height=0.5,cex=1) +
  xlab("r(95% CI)") + ylab("Random Sample Number")

FPM_scaled_rmvn100_CIplot_theta <- ggplot(data=FPM_scaled_rmvnParamRec, aes(y=ID, x=initial_theta, xmin=l95CI_theta, xmax=u95CI_theta)) +
  geom_point(aes(col=theta_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_theta, xmax=u95CI_theta, col=theta_fit),height=0.5,cex=1) +
  xlab("theta(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_procNoise_CIplot_theta <- ggplot(data=FPM_procNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_theta, xmin=l95CI_theta, xmax=u95CI_theta)) +
  geom_point(aes(col=theta_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_theta, xmax=u95CI_theta, col=theta_fit),height=0.5,cex=1) +
  xlab("theta(95% CI)") + ylab("Random Sample Number")

FPM_rmvn100_procNoise_measureNoise_CIplot_theta <- ggplot(data=FPM_procNoise_measureNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_theta, xmin=l95CI_theta, xmax=u95CI_theta)) +
  geom_point(aes(col=theta_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_theta, xmax=u95CI_theta, col=theta_fit),height=0.5,cex=1) +
  xlab("theta(95% CI)") + ylab("Random Sample Number")

##r ReParam
FPM_rReParam_rmvn100_procNoise_CIplot_r <- ggplot(data=FPM_procNoise_rReParam_scaled_rmvnParamRec, aes(y=ID, x=initial_r, xmin=l95CI_r, xmax=u95CI_r)) +
  geom_point(aes(col=r_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_r, xmax=u95CI_r, col=r_fit),height=0.5,cex=1) +
  xlab("Re-parameterization of r as e^r; r(95% CI)") + ylab("Random Sample Number")

#### Parameter Recovery CI graphs - FPM ####
#plot CIs in forestplot-like plot 
Kissileff_scaled_rmvn100_CIplot_int <- ggplot(data=Kissileff_scaled_rmvnParamRec, aes(y=ID, x=initial_int, xmin=l95CI_int, xmax=u95CI_int)) +
  geom_point(aes(col=int_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_int, xmax=u95CI_int, col=int_fit),height=0.5,cex=1) +
  xlab("Intercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_CIplot_int <- ggplot(data=Kissileff_procNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_int, xmin=l95CI_int, xmax=u95CI_int)) +
  geom_point(aes(col=int_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_int, xmax=u95CI_int, col=int_fit),height=0.5,cex=1) +
  xlab("Intercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_measureNoise_CIplot_int <- ggplot(data=Kissileff_procNoise_measureNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_int, xmin=l95CI_int, xmax=u95CI_int)) +
  geom_point(aes(col=int_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_int, xmax=u95CI_int, col=int_fit),height=0.5,cex=1) +
  xlab("Intercept (95% CI)") + ylab("Random Sample Number")

Kissileff_scaled_rmvn100_CIplot_linear <- ggplot(data=Kissileff_scaled_rmvnParamRec, aes(y=ID, x=initial_linear, xmin=l95CI_linear, xmax=u95CI_linear)) +
  geom_point(aes(col=linear_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_linear, xmax=u95CI_linear, col=linear_fit),height=0.5,cex=1) +
  xlab("linearercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_CIplot_linear <- ggplot(data=Kissileff_procNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_linear, xmin=l95CI_linear, xmax=u95CI_linear)) +
  geom_point(aes(col=linear_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_linear, xmax=u95CI_linear, col=linear_fit),height=0.5,cex=1) +
  xlab("linearercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_measureNoise_CIplot_linear <- ggplot(data=Kissileff_procNoise_measureNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_linear, xmin=l95CI_linear, xmax=u95CI_linear)) +
  geom_point(aes(col=linear_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_linear, xmax=u95CI_linear, col=linear_fit),height=0.5,cex=1) +
  xlab("linearercept (95% CI)") + ylab("Random Sample Number")

Kissileff_scaled_rmvn100_CIplot_quad <- ggplot(data=Kissileff_scaled_rmvnParamRec, aes(y=ID, x=initial_quad, xmin=l95CI_quad, xmax=u95CI_quad)) +
  geom_point(aes(col=quad_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_quad, xmax=u95CI_quad, col=quad_fit),height=0.5,cex=1) +
  xlab("quadercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_CIplot_quad <- ggplot(data=Kissileff_procNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_quad, xmin=l95CI_quad, xmax=u95CI_quad)) +
  geom_point(aes(col=quad_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_quad, xmax=u95CI_quad, col=quad_fit),height=0.5,cex=1) +
  xlab("quadercept (95% CI)") + ylab("Random Sample Number")

Kissileff_rmvn100_procNoise_measureNoise_CIplot_quad <- ggplot(data=Kissileff_procNoise_measureNoise_scaled_rmvnParamRec, aes(y=ID, x=initial_quad, xmin=l95CI_quad, xmax=u95CI_quad)) +
  geom_point(aes(col=quad_fit)) + 
  geom_errorbarh(aes(xmin=l95CI_quad, xmax=u95CI_quad, col=quad_fit),height=0.5,cex=1) +
  xlab("quadercept (95% CI)") + ylab("Random Sample Number")

#### Get bites for 1 random sample ####
#FPM 
FPM_scaled_rmvnParamRec_bothFit = FPM_scaled_rmvnParamRec[FPM_scaled_rmvnParamRec$r_fit == TRUE & FPM_scaled_rmvnParamRec$theta_fit == TRUE, ]
FPM_scaled_rmvnParamRec_true_Bites = simBites(nBites = FPM_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_scaled_rmvnParamRec_bothFit$initial_theta[1], FPM_scaled_rmvnParamRec_bothFit$initial_r[1]), time_fn = FPM_Time, procNoise = FALSE)
FPM_scaled_rmvnParamRec_fit_Bites = simBites(nBites = FPM_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_scaled_rmvnParamRec_bothFit$theta[1], FPM_scaled_rmvnParamRec_bothFit$r[1]), time_fn = FPM_Time, procNoise = FALSE)
FPM_scaled_rmvnParamRec_upper95CI_Bites = simBites(nBites = FPM_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_scaled_rmvnParamRec_bothFit$u95CI_theta[1], FPM_scaled_rmvnParamRec_bothFit$u95CI_r[1]), time_fn = FPM_Time, procNoise = FALSE)
FPM_scaled_rmvnParamRec_lower95CI_Bites = simBites(nBites = FPM_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_scaled_rmvnParamRec_bothFit$l95CI_theta[1], FPM_scaled_rmvnParamRec_bothFit$l95CI_r[1]), time_fn = FPM_Time, procNoise = FALSE)

FPM_procNoise_scaled_rmvnParamRec_bothFit = FPM_procNoise_scaled_rmvnParamRec[FPM_procNoise_scaled_rmvnParamRec$r_fit == TRUE & FPM_procNoise_scaled_rmvnParamRec$theta_fit == TRUE, ]
FPM_procNoise_scaled_rmvnParamRec_true_Bites = simBites(nBites = FPM_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_scaled_rmvnParamRec_bothFit$initial_theta[1], FPM_procNoise_scaled_rmvnParamRec_bothFit$initial_r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_scaled_rmvnParamRec_fit_Bites = simBites(nBites = FPM_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_scaled_rmvnParamRec_bothFit$theta[1], FPM_procNoise_scaled_rmvnParamRec_bothFit$r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_scaled_rmvnParamRec_upper95CI_Bites = simBites(nBites = FPM_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_scaled_rmvnParamRec_bothFit$u95CI_theta[1], FPM_procNoise_scaled_rmvnParamRec_bothFit$u95CI_r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_scaled_rmvnParamRec_lower95CI_Bites = simBites(nBites = FPM_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_scaled_rmvnParamRec_bothFit$l95CI_theta[1], FPM_procNoise_scaled_rmvnParamRec_bothFit$l95CI_r[1]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit = FPM_procNoise_measureNoise_scaled_rmvnParamRec[FPM_procNoise_scaled_rmvnParamRec$r_fit == TRUE & FPM_procNoise_measureNoise_scaled_rmvnParamRec$theta_fit == TRUE, ]
FPM_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites = simBites(nBites = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$initial_theta[1], FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$initial_r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites = simBites(nBites = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$theta[1], FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_measureNoise_scaled_rmvnParamRec_upper95CI_Bites = simBites(nBites = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$u95CI_theta[1], FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$u95CI_r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_measureNoise_scaled_rmvnParamRec_lower95CI_Bites = simBites(nBites = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$l95CI_theta[1], FPM_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$l95CI_r[1]), time_fn = FPM_Time, procNoise = TRUE)

## Kissileff
Kissileff_scaled_rmvnParamRec_bothFit = Kissileff_scaled_rmvnParamRec[Kissileff_scaled_rmvnParamRec$int_fit == TRUE & Kissileff_scaled_rmvnParamRec$linear_fit == TRUE & Kissileff_scaled_rmvnParamRec$quad_fit == TRUE, ]
Kissileff_scaled_rmvnParamRec_true_Bites = simBites(nBites = Kissileff_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_scaled_rmvnParamRec_bothFit$initial_int[1], Kissileff_scaled_rmvnParamRec_bothFit$initial_linear[1], Kissileff_scaled_rmvnParamRec_bothFit$initial_quad[1]), time_fn = Kissileff_Time, procNoise = FALSE)
Kissileff_scaled_rmvnParamRec_fit_Bites = simBites(nBites = Kissileff_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_scaled_rmvnParamRec_bothFit$int[1], Kissileff_scaled_rmvnParamRec_bothFit$linear[1], Kissileff_scaled_rmvnParamRec_bothFit$quad[1]), time_fn = Kissileff_Time, procNoise = FALSE)
Kissileff_scaled_rmvnParamRec_upper95CI_Bites = simBites(nBites = Kissileff_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_scaled_rmvnParamRec_bothFit$u95CI_int[1], Kissileff_scaled_rmvnParamRec_bothFit$u95CI_linear[1], Kissileff_scaled_rmvnParamRec_bothFit$u95CI_quad[1]), time_fn = Kissileff_Time, procNoise = FALSE)
Kissileff_scaled_rmvnParamRec_lower95CI_Bites = simBites(nBites = Kissileff_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_scaled_rmvnParamRec_bothFit$l95CI_int[1], Kissileff_scaled_rmvnParamRec_bothFit$l95CI_linear[1], Kissileff_scaled_rmvnParamRec_bothFit$l95CI_quad[1]), time_fn = Kissileff_Time, procNoise = FALSE)

Kissileff_procNoise_scaled_rmvnParamRec_bothFit = Kissileff_procNoise_scaled_rmvnParamRec[Kissileff_procNoise_scaled_rmvnParamRec$int_fit == TRUE & Kissileff_procNoise_scaled_rmvnParamRec$linear_fit == TRUE & Kissileff_procNoise_scaled_rmvnParamRec$quad_fit == TRUE, ]
Kissileff_procNoise_scaled_rmvnParamRec_true_Bites = simBites(nBites = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_scaled_rmvnParamRec_bothFit$initial_int[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$initial_linear[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$initial_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_scaled_rmvnParamRec_fit_Bites = simBites(nBites = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_scaled_rmvnParamRec_bothFit$int[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$linear[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_scaled_rmvnParamRec_upper95CI_Bites = simBites(nBites = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_scaled_rmvnParamRec_bothFit$u95CI_int[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$u95CI_linear[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$u95CI_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_scaled_rmvnParamRec_lower95CI_Bites = simBites(nBites = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_scaled_rmvnParamRec_bothFit$l95CI_int[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$l95CI_linear[1], Kissileff_procNoise_scaled_rmvnParamRec_bothFit$l95CI_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)

Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec[Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$int_fit == TRUE & Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$linear_fit == TRUE & Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$quad_fit == TRUE, ]
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites = simBites(nBites = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$initial_int[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$initial_linear[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$initial_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites = simBites(nBites = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$int[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$linear[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_upper95CI_Bites = simBites(nBites = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$u95CI_int[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$u95CI_linear[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$u95CI_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_lower95CI_Bites = simBites(nBites = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$nBites[1], Emax = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$Emax[1], parameters = c(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$l95CI_int[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$l95CI_linear[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_bothFit$l95CI_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)

#### Plot Btes with CI bounds ####
FPM_scaled_rmvnParamRec_95CI_CumulativeIntake = ggplot(FPM_scaled_rmvnParamRec_true_Bites, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avgBite)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_scaled_rmvnParamRec_fit_Bites, color = 'blue') +
  stat_smooth(data = FPM_scaled_rmvnParamRec_upper95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = FPM_scaled_rmvnParamRec_lower95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI - No Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

FPM_procNoise_scaled_rmvnParamRec_95CI_CumulativeIntake = ggplot(FPM_procNoise_scaled_rmvnParamRec_true_Bites, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_scaled_rmvnParamRec_fit_Bites, color = 'blue') +
  stat_smooth(data = FPM_procNoise_scaled_rmvnParamRec_upper95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = FPM_procNoise_scaled_rmvnParamRec_lower95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI - Process Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

FPM_procNoise_measureNoise_scaled_rmvnParamRec_95CI_CumulativeIntake = ggplot(FPM_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites, color = 'blue') +
  stat_smooth(data = FPM_procNoise_measureNoise_scaled_rmvnParamRec_upper95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = FPM_procNoise_measureNoise_scaled_rmvnParamRec_lower95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI - Process Noise and Measurement Error') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

Kissileff_scaled_rmvnParamRec_95CI_CumulativeIntake = ggplot(Kissileff_scaled_rmvnParamRec_true_Bites, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avgBite)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = Kissileff_scaled_rmvnParamRec_fit_Bites, color = 'blue') +
  stat_smooth(data = Kissileff_scaled_rmvnParamRec_upper95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = Kissileff_scaled_rmvnParamRec_lower95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI - No Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

Kissileff_procNoise_scaled_rmvnParamRec_95CI_CumulativeIntake = ggplot(Kissileff_procNoise_scaled_rmvnParamRec_true_Bites, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = Kissileff_procNoise_scaled_rmvnParamRec_fit_Bites, color = 'blue') +
  stat_smooth(data = Kissileff_procNoise_scaled_rmvnParamRec_upper95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = Kissileff_procNoise_scaled_rmvnParamRec_lower95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI - Process Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_95CI_CumulativeIntake = ggplot(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites, color = 'blue') +
  stat_smooth(data = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_upper95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  stat_smooth(data = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_lower95CI_Bites, method = 'gam', formula = y~s(x), linetype = 2, color = 'cornflowerblue') +
  ggtitle('Parameter Recovery with 95%CI - Process Noise and Measurement Error') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

#### Generate Cumulative Intake Curves Poor Fit of Parameters ####
##FPM
# FPM_scaled_rmvnParamRec_noFit = FPM_scaled_rmvnParamRec[FPM_scaled_rmvnParamRec$r_fit == FALSE & FPM_scaled_rmvnParamRec$theta_fit == FALSE, ]
FPM_scaled_rmvnParamRec_noFit = FPM_scaled_rmvnParamRec[FPM_scaled_rmvnParamRec$r_fit == FALSE, ]

FPM_scaled_rmvnParamRec_true_Bites_noFit = simBites(nBites = FPM_scaled_rmvnParamRec_noFit$nBites[1], Emax = FPM_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(FPM_scaled_rmvnParamRec_noFit$initial_theta[1], FPM_scaled_rmvnParamRec_noFit$initial_r[1]), time_fn = FPM_Time, procNoise = FALSE)
FPM_scaled_rmvnParamRec_fit_Bites_noFit = simBites(nBites = FPM_scaled_rmvnParamRec_noFit$nBites[1], Emax = FPM_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(FPM_scaled_rmvnParamRec_noFit$theta[1], FPM_scaled_rmvnParamRec_noFit$r[1]), time_fn = FPM_Time, procNoise = FALSE)

FPM_procNoise_scaled_rmvnParamRec_noFit = FPM_procNoise_scaled_rmvnParamRec[FPM_procNoise_scaled_rmvnParamRec$r_fit == FALSE & FPM_procNoise_scaled_rmvnParamRec$theta_fit == FALSE, ]
FPM_procNoise_scaled_rmvnParamRec_true_Bites_noFit = simBites(nBites = FPM_procNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = FPM_procNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(FPM_procNoise_scaled_rmvnParamRec_noFit$initial_theta[1], FPM_procNoise_scaled_rmvnParamRec_noFit$initial_r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_scaled_rmvnParamRec_fit_Bites_noFit = simBites(nBites = FPM_procNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = FPM_procNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(FPM_procNoise_scaled_rmvnParamRec_noFit$theta[1], FPM_procNoise_scaled_rmvnParamRec_noFit$r[1]), time_fn = FPM_Time, procNoise = TRUE)

FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit = FPM_procNoise_measureNoise_scaled_rmvnParamRec[FPM_procNoise_scaled_rmvnParamRec$r_fit == FALSE & FPM_procNoise_measureNoise_scaled_rmvnParamRec$theta_fit == FALSE, ]
FPM_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites_noFit = simBites(nBites = FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$initial_theta[1], FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$initial_r[1]), time_fn = FPM_Time, procNoise = TRUE)
FPM_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites_noFit = simBites(nBites = FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$theta[1], FPM_procNoise_measureNoise_scaled_rmvnParamRec_noFit$r[1]), time_fn = FPM_Time, procNoise = TRUE)

##Kissileff
Kissileff_scaled_rmvnParamRec_noFit = Kissileff_scaled_rmvnParamRec[Kissileff_scaled_rmvnParamRec$int_fit == FALSE & Kissileff_scaled_rmvnParamRec$linear_fit == FALSE & Kissileff_scaled_rmvnParamRec$quad_fit == FALSE, ]
Kissileff_scaled_rmvnParamRec_true_Bites_noFit = simBites(nBites = Kissileff_scaled_rmvnParamRec_noFit$nBites[1], Emax = Kissileff_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(Kissileff_scaled_rmvnParamRec_noFit$initial_int[1], Kissileff_scaled_rmvnParamRec_noFit$initial_linear[1], Kissileff_scaled_rmvnParamRec_noFit$initial_quad[1]), time_fn = Kissileff_Time, procNoise = FALSE)
Kissileff_scaled_rmvnParamRec_fit_Bites_noFit = simBites(nBites = Kissileff_scaled_rmvnParamRec_noFit$nBites[1], Emax = Kissileff_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(Kissileff_scaled_rmvnParamRec_noFit$int[1], Kissileff_scaled_rmvnParamRec_noFit$linear[1], Kissileff_scaled_rmvnParamRec_noFit$quad[1]), time_fn = Kissileff_Time, procNoise = FALSE)

Kissileff_procNoise_scaled_rmvnParamRec_noFit = Kissileff_procNoise_scaled_rmvnParamRec[Kissileff_procNoise_scaled_rmvnParamRec$int_fit == FALSE & Kissileff_procNoise_scaled_rmvnParamRec$linear_fit == FALSE & Kissileff_procNoise_scaled_rmvnParamRec$quad_fit == FALSE, ]
Kissileff_procNoise_scaled_rmvnParamRec_true_Bites_noFit = simBites(nBites = Kissileff_procNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = Kissileff_procNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(Kissileff_procNoise_scaled_rmvnParamRec_noFit$initial_int[1], Kissileff_procNoise_scaled_rmvnParamRec_noFit$initial_linear[1], Kissileff_procNoise_scaled_rmvnParamRec_noFit$initial_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_scaled_rmvnParamRec_fit_Bites_noFit = simBites(nBites = Kissileff_procNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = Kissileff_procNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(Kissileff_procNoise_scaled_rmvnParamRec_noFit$int[1], Kissileff_procNoise_scaled_rmvnParamRec_noFit$linear[1], Kissileff_procNoise_scaled_rmvnParamRec_noFit$quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)

Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec[Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$int_fit == FALSE & Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$linear_fit == FALSE & Kissileff_procNoise_measureNoise_scaled_rmvnParamRec$quad_fit == FALSE, ]
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites_noFit = simBites(nBites = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$initial_int[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$initial_linear[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$initial_quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)
Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites_noFit = simBites(nBites = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$nBites[1], Emax = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$Emax[1], parameters = c(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$int[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$linear[1], Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_noFit$quad[1]), time_fn = Kissileff_Time, procNoise = TRUE)

#### Compare Cumulative Intake Curves: poor param Rec vs 'True' ####
##FPM
FPM_scaled_rmvnParamRec_CumulativeIntake_noFit = ggplot(FPM_scaled_rmvnParamRec_true_Bites_noFit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avgBite)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_scaled_rmvnParamRec_fit_Bites_noFit, color = 'blue') +
  ggtitle('Parameter Recovery - No Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

FPM_procNoise_scaled_rmvnParamRec_95CI_CumulativeIntake_noFit = ggplot(FPM_procNoise_scaled_rmvnParamRec_true_Bites_noFit, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_scaled_rmvnParamRec_fit_Bites_noFit, color = 'blue') +
  ggtitle('Parameter Recovery - Process Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

FPM_procNoise_measureNoise_scaled_rmvnParamRec_95CI_CumulativeIntake_noFit = ggplot(FPM_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites_noFit, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = FPM_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites_noFit, color = 'blue') +
  ggtitle('Parameter Recovery - Process Noise and Measurement Error') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

##Kissileff 
Kissileff_scaled_rmvnParamRec_CumulativeIntake_noFit = ggplot(Kissileff_scaled_rmvnParamRec_true_Bites_noFit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avgBite)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = Kissileff_scaled_rmvnParamRec_fit_Bites_noFit, color = 'blue') +
  ggtitle('Parameter Recovery - No Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

Kissileff_procNoise_scaled_rmvnParamRec_95CI_CumulativeIntake_noFit = ggplot(Kissileff_procNoise_scaled_rmvnParamRec_true_Bites_noFit, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = Kissileff_procNoise_scaled_rmvnParamRec_fit_Bites_noFit, color = 'blue') +
  ggtitle('Parameter Recovery - Process Noise') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')

Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_95CI_CumulativeIntake_noFit = ggplot(Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_true_Bites_noFit, aes(y = CumulativeGrams_procNoise, x = EstimatedTime_procNoise)) +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(data = Kissileff_procNoise_measureNoise_scaled_rmvnParamRec_fit_Bites_noFit, color = 'blue') +
  ggtitle('Parameter Recovery - Process Noise and Measurement Error') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)')