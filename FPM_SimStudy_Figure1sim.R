############ Basic Data Load/Setup ############
library(reporttools)
library(xtable)
library(car)
library(ggplot2)
library(reshape2)
library(stats)
library(rstudioapi)
library(psych)
library(bitemodelr)

#### set up ####

#set piorking directory to location of script--not needed pihen called
#through Rmarkdopin doc. Uncomment belopi if running locally/manually
# this.dir = getActiveDocumentContext()$path
# setwd(dirname(this.dir))

source('functions.R')

#### Source Setup Scripts ####
source('FPM_incorrect.R')
#### Generate curves from Figure 1 with bite level data ####
# Figure 1A: theta = 30g/min, r = 0.17 1/min, Emax = 400
# Figure 1B: theta = 30g/min, r = 0.17 1/min, Emax = 600
# Figure 1C: theta = 30g/min, r = 0.5 1/min, Emax = 400
# Figure 1C: theta = 50g/min, r = 0.17 1/min, Emax = 400

Figure1Params = data.frame(matrix("", nrow = 4, ncol = 5))
names(Figure1Params) = c('figure', 'theta', 'r', 'Emax', 'MealDuration')

Figure1Params$figure = c('A', 'B', 'C', 'D')
Figure1Params$theta = 30
Figure1Params$r = 0.17
Figure1Params$Emax = 400
Figure1Params$MealDuration = 30

Figure1Params$Emax[2] = 600
Figure1Params$r[3] = 0.5
Figure1Params$theta[4] = 50

#### Replication with sampling every 250 ms ####
##orignial equation
cont_CumulativeIntake_Figure1A_orig = simBites(time_fn = FPMincorrect_Time, nBites = 7200, parameters = c(Figure1Params$theta[1], Figure1Params$r[1]), Emax = Figure1Params$Emax[1])
cont_CumulativeIntake_Figure1A_orig$Figure = '1A'
cont_CumulativeIntake_Figure1B_orig = simBites(time_fn = FPMincorrect_Time, nBites = 7200, parameters = c(Figure1Params$theta[2], Figure1Params$r[2]), Emax = Figure1Params$Emax[2])
cont_CumulativeIntake_Figure1B_orig$Figure = '1B'
cont_CumulativeIntake_Figure1C_orig = simBites(time_fn = FPMincorrect_Time, nBites = 7200, parameters = c(Figure1Params$theta[3], Figure1Params$r[3]), Emax = Figure1Params$Emax[3])
cont_CumulativeIntake_Figure1C_orig$Figure = '1C'
cont_CumulativeIntake_Figure1D_orig = simBites(time_fn = FPMincorrect_Time, nBites = 7200, parameters = c(Figure1Params$theta[4], Figure1Params$r[4]), Emax = Figure1Params$Emax[4])
cont_CumulativeIntake_Figure1D_orig$Figure = '1D'

cont_CumulativeIntake_Figure1_orig_long = rbind.data.frame(cont_CumulativeIntake_Figure1A_orig, cont_CumulativeIntake_Figure1B_orig, cont_CumulativeIntake_Figure1C_orig, cont_CumulativeIntake_Figure1D_orig)
cont_CumulativeIntake_Figure1_orig_long$group = 'cont'

#corrected equation
cont_CumulativeIntake_Figure1A_cor = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(Figure1Params$theta[1], Figure1Params$r[1]), Emax = Figure1Params$Emax[1])
cont_CumulativeIntake_Figure1A_cor$Figure = '1A'
cont_CumulativeIntake_Figure1B_cor = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(Figure1Params$theta[2], Figure1Params$r[2]), Emax = Figure1Params$Emax[2])
cont_CumulativeIntake_Figure1B_cor$Figure = '1B'
cont_CumulativeIntake_Figure1C_cor = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(Figure1Params$theta[3], Figure1Params$r[3]), Emax = Figure1Params$Emax[3])
cont_CumulativeIntake_Figure1C_cor$Figure = '1C'
cont_CumulativeIntake_Figure1D_cor = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(Figure1Params$theta[4], Figure1Params$r[4]), Emax = Figure1Params$Emax[4])
cont_CumulativeIntake_Figure1D_cor$Figure = '1D'

cont_CumulativeIntake_Figure1_cor_long = rbind.data.frame(cont_CumulativeIntake_Figure1A_cor, cont_CumulativeIntake_Figure1B_cor, cont_CumulativeIntake_Figure1C_cor, cont_CumulativeIntake_Figure1D_cor)
cont_CumulativeIntake_Figure1_cor_long$group = 'cont'

cont_CumulativeIntake_Figure1_plotA = ggplot(cont_CumulativeIntake_Figure1A_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  stat_smooth(data = cont_CumulativeIntake_Figure1A_orig, method = 'gam', formula = y~s(x), linetype = 2, color = 'black') +
  ggtitle('Replication of Figure 1A with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

cont_CumulativeIntake_Figure1_plotB = ggplot(cont_CumulativeIntake_Figure1B_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  stat_smooth(data = cont_CumulativeIntake_Figure1B_orig, method = 'gam', formula = y~s(x), linetype = 2, color = 'black') +
  ggtitle('Replication of Figure 1B with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

cont_CumulativeIntake_Figure1_plotC = ggplot(cont_CumulativeIntake_Figure1C_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  stat_smooth(data = cont_CumulativeIntake_Figure1C_orig, method = 'gam', formula = y~s(x), linetype = 2, color = 'black') +
  ggtitle('Replication of Figure 1C with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

cont_CumulativeIntake_Figure1_plotD = ggplot(cont_CumulativeIntake_Figure1D_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  stat_smooth(data = cont_CumulativeIntake_Figure1D_orig, method = 'gam', formula = y~s(x), linetype = 2, color = 'black') +
  ggtitle('Replication of Figure 1D with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

#facet wrap
cont_CumulativeIntake_Figure1_grid = ggplot(cont_CumulativeIntake_Figure1_cor_long, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x),linetype = 1, color = 'black') +
  stat_smooth(data = cont_CumulativeIntake_Figure1_orig_long, method = 'gam', formula = y~s(x), linetype = 2, color = 'black') +
  ggtitle('Replication of Figure 1 with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~Figure)

#### Replication with 100 bite sampling ####

#corrected equation
bite100_CumulativeIntake_Figure1A_cor = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(Figure1Params$theta[1], Figure1Params$r[1]), Emax = Figure1Params$Emax[1])
bite100_CumulativeIntake_Figure1A_cor$Figure = '1A'
bite100_CumulativeIntake_Figure1B_cor = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(Figure1Params$theta[2], Figure1Params$r[2]), Emax = Figure1Params$Emax[2])
bite100_CumulativeIntake_Figure1B_cor$Figure = '1B'
bite100_CumulativeIntake_Figure1C_cor = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(Figure1Params$theta[3], Figure1Params$r[3]), Emax = Figure1Params$Emax[3])
bite100_CumulativeIntake_Figure1C_cor$Figure = '1C'
bite100_CumulativeIntake_Figure1D_cor = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(Figure1Params$theta[4], Figure1Params$r[4]), Emax = Figure1Params$Emax[4])
bite100_CumulativeIntake_Figure1D_cor$Figure = '1D'

bite100_CumulativeIntake_Figure1_cor_long = rbind.data.frame(bite100_CumulativeIntake_Figure1A_cor, bite100_CumulativeIntake_Figure1B_cor, bite100_CumulativeIntake_Figure1C_cor, bite100_CumulativeIntake_Figure1D_cor)
bite100_CumulativeIntake_Figure1_cor_long$group = 'bite100'

bite100_CumulativeIntake_Figure1_plotA = ggplot(bite100_CumulativeIntake_Figure1A_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1A_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1A with 100 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

bite100_CumulativeIntake_Figure1_plotB = ggplot(bite100_CumulativeIntake_Figure1B_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1B_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1B with 100 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

bite100_CumulativeIntake_Figure1_plotC = ggplot(bite100_CumulativeIntake_Figure1C_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1C_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1C with 100 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

bite100_CumulativeIntake_Figure1_plotD = ggplot(bite100_CumulativeIntake_Figure1D_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1D_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1D with 100 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

#facet wrap
bite100_CumulativeIntake_Figure1_grid = ggplot(bite100_CumulativeIntake_Figure1_cor_long, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1_cor_long, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1 with 100 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~Figure)

#### Replication with 33 bite sampling ####

#corrected equation
bite33_CumulativeIntake_Figure1A_cor = ParamRecovery(nBites = 30, Emax = Figure1Params$Emax[1], parameters = c(Figure1Params$theta[1], Figure1Params$r[1]), time_fn = FPM_Time, fit_fn = FPM_Fit, simVar = 'bitesSampled', simValue = 5, keepBites = TRUE)
bite33_CumulativeIntake_Figure1A_Et = bite33_CumulativeIntake_Figure1A_cor$biteDat_paramRecov
bite33_CumulativeIntake_Figure1A_Et$Figure = '1A'
bite33_CumulativeIntake_Figure1A_paramsRec = bite33_CumulativeIntake_Figure1A_cor$paramDat

bite33_CumulativeIntake_Figure1B_cor = ParamRecovery(nBites = 30, Emax = Figure1Params$Emax[2], parameters = c(Figure1Params$theta[2], Figure1Params$r[2]), time_fn = FPM_Time, fit_fn = FPM_Fit, simVar = 'bitesSampled', simValue = 5, keepBites = TRUE)
bite33_CumulativeIntake_Figure1B_Et = bite33_CumulativeIntake_Figure1B_cor$biteDat_paramRecov
bite33_CumulativeIntake_Figure1B_Et$Figure = '1B'
bite33_CumulativeIntake_Figure1B_paramsRec = bite33_CumulativeIntake_Figure1B_cor$paramDat

bite33_CumulativeIntake_Figure1C_cor = ParamRecovery(nBites = 30, Emax = Figure1Params$Emax[3], parameters = c(Figure1Params$theta[3], Figure1Params$r[3]), time_fn = FPM_Time, fit_fn = FPM_Fit, simVar = 'bitesSampled', simValue = 5, keepBites = TRUE)
bite33_CumulativeIntake_Figure1C_Et = bite33_CumulativeIntake_Figure1C_cor$biteDat_paramRecov
bite33_CumulativeIntake_Figure1C_Et$Figure = '1C'
bite33_CumulativeIntake_Figure1C_paramsRec = bite33_CumulativeIntake_Figure1C_cor$paramDat

bite33_CumulativeIntake_Figure1D_cor = ParamRecovery(nBites = 30, Emax = Figure1Params$Emax[4], parameters = c(Figure1Params$theta[4], Figure1Params$r[4]), time_fn = FPM_Time, fit_fn = FPM_Fit, simVar = 'bitesSampled', simValue = 5, keepBites = TRUE)
bite33_CumulativeIntake_Figure1D_Et = bite33_CumulativeIntake_Figure1D_cor$biteDat_paramRecov
bite33_CumulativeIntake_Figure1D_Et$Figure = '1D'
bite33_CumulativeIntake_Figure1D_paramsRec = bite33_CumulativeIntake_Figure1D_cor$paramDat

#Graphs
bite33_CumulativeIntake_Figure1_plotA = ggplot(bite33_CumulativeIntake_Figure1A_Et[bite33_CumulativeIntake_Figure1A_Et$simNum == 1, ], aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1A_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1A with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite33_CumulativeIntake_Figure1_plotB = ggplot(bite33_CumulativeIntake_Figure1B_Et[bite33_CumulativeIntake_Figure1B_Et$simNum == 1, ], aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1B_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1B with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite33_CumulativeIntake_Figure1_plotC = ggplot(bite33_CumulativeIntake_Figure1C_Et[bite33_CumulativeIntake_Figure1C_Et$simNum == 1, ], aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1C_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1C with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite33_CumulativeIntake_Figure1_plotD = ggplot(bite33_CumulativeIntake_Figure1D_Et[bite33_CumulativeIntake_Figure1D_Et$simNum == 1, ], aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1D_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1D with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

#facet wrap
bite33_CumulativeIntake_Figure1A_Et_33bite = bite33_CumulativeIntake_Figure1A_Et[bite33_CumulativeIntake_Figure1A_Et$nBites == 33, ]
bite33_CumulativeIntake_Figure1B_Et_33bite = bite33_CumulativeIntake_Figure1B_Et[bite33_CumulativeIntake_Figure1B_Et$nBites == 33, ]
bite33_CumulativeIntake_Figure1C_Et_33bite = bite33_CumulativeIntake_Figure1C_Et[bite33_CumulativeIntake_Figure1C_Et$nBites == 33, ]
bite33_CumulativeIntake_Figure1D_Et_33bite = bite33_CumulativeIntake_Figure1D_Et[bite33_CumulativeIntake_Figure1D_Et$nBites == 33, ]

bite33_CumulativeIntake_Figure1A_Et_33bite = bite33_CumulativeIntake_Figure1A_Et_33bite[bite33_CumulativeIntake_Figure1A_Et_33bite$simNum == min(bite33_CumulativeIntake_Figure1A_Et_33bite$simNum), ]
bite33_CumulativeIntake_Figure1B_Et_33bite = bite33_CumulativeIntake_Figure1B_Et_33bite[bite33_CumulativeIntake_Figure1B_Et_33bite$simNum == min(bite33_CumulativeIntake_Figure1B_Et_33bite$simNum), ]
bite33_CumulativeIntake_Figure1C_Et_33bite = bite33_CumulativeIntake_Figure1C_Et_33bite[bite33_CumulativeIntake_Figure1C_Et_33bite$simNum == min(bite33_CumulativeIntake_Figure1C_Et_33bite$simNum), ]
bite33_CumulativeIntake_Figure1D_Et_33bite = bite33_CumulativeIntake_Figure1D_Et_33bite[bite33_CumulativeIntake_Figure1D_Et_33bite$simNum == min(bite33_CumulativeIntake_Figure1D_Et_33bite$simNum), ]

bite33_CumulativeIntake_Figure1_Et_long = rbind.data.frame(bite33_CumulativeIntake_Figure1A_Et_33bite, bite33_CumulativeIntake_Figure1B_Et_33bite, bite33_CumulativeIntake_Figure1C_Et_33bite, bite33_CumulativeIntake_Figure1D_Et_33bite)

bite33_CumulativeIntake_Figure1_grid = ggplot(bite33_CumulativeIntake_Figure1_Et_long, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1_cor_long, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_point(shape = 23, color = 'red', size = 3) +
  ggtitle('Replication of Figure 1 with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~Figure)

#### Plot all data by sampling rate - corrected only ####
bite100_CumulativeIntake_Figure1A_cor$nBites = 100
cont_CumulativeIntake_Figure1A_cor$nBites = 7200
all_CumulativeIntake_Figure1_plotA = ggplot(cont_CumulativeIntake_Figure1A_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x), se = FALSE, color = 'black') +
  geom_point(data = bite100_CumulativeIntake_Figure1A_cor, shape = 23, aes(color = factor(nBites)), size = 3) +
  geom_point(data = bite33_CumulativeIntake_Figure1A_Et_33bite, shape = 23, aes(color = factor(nBites)), size = 3) +
  scale_color_manual(name = 'Sampling Rate',
                     breaks = c( '7200', '100', '33'),
                     values = c('black', 'blue', 'red'),
                     labels = c('250ms', '100 bites/30 min', '33 bites/30min'))+
  ggtitle('Replication of Figure 1A with All Sampling') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

bite100_CumulativeIntake_Figure1B_cor$nBites = 100
cont_CumulativeIntake_Figure1B_cor$nBites = 7200
all_CumulativeIntake_Figure1_plotB = ggplot(cont_CumulativeIntake_Figure1B_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x), se = FALSE, color = 'black') +
  geom_point(data = bite100_CumulativeIntake_Figure1B_cor, shape = 23, aes(color = factor(nBites)), size = 3) +
  geom_point(data = bite33_CumulativeIntake_Figure1B_Et_33bite, shape = 23, aes(color = factor(nBites)), size = 3) +
  scale_color_manual(name = 'Sampling Rate',
    breaks = c( '7200', '100', '33'),
    values = c('black', 'blue', 'red'),
    labels = c('250ms', '100 bites/30 min', '33 bites/30min'))+
  ggtitle('Replication of Figure 1B with All Sampling') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite100_CumulativeIntake_Figure1C_cor$nBites = 100
cont_CumulativeIntake_Figure1C_cor$nBites = 7200
all_CumulativeIntake_Figure1_plotC = ggplot(cont_CumulativeIntake_Figure1C_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x), se = FALSE, color = 'black') +
  geom_point(data = bite100_CumulativeIntake_Figure1C_cor, shape = 23, aes(color = factor(nBites)), size = 3) +
  geom_point(data = bite33_CumulativeIntake_Figure1C_Et_33bite, shape = 23, aes(color = factor(nBites)), size = 3) +
  scale_color_manual(name = 'Sampling Rate',
    breaks = c( '7200', '100', '33'),
    values = c('black', 'blue', 'red'),
    labels = c('250ms', '100 bites/30 min', '33 bites/30min'))+
  ggtitle('Replication of Figure 1C with All Sampling') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite100_CumulativeIntake_Figure1D_cor$nBites = 100
cont_CumulativeIntake_Figure1D_cor$nBites = 7200
all_CumulativeIntake_Figure1_plotD = ggplot(cont_CumulativeIntake_Figure1D_cor, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  geom_smooth(method = 'gam', formula = y~s(x), se = FALSE, color = 'black') +
  geom_point(data = bite100_CumulativeIntake_Figure1D_cor, shape = 23, aes(color = factor(nBites)), size = 3) +
  geom_point(data = bite33_CumulativeIntake_Figure1D_Et_33bite, shape = 23, aes(color = factor(nBites)), size = 3) +
  scale_color_manual(name = 'Sampling Rate',
    breaks = c( '7200', '100', '33'),
    values = c('black', 'blue', 'red'),
    labels = c('250ms', '100 bites/30 min', '33 bites/30min'))+
  ggtitle('Replication of Figure 1D with All Sampling') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

#### Try to Recover Parameters from 250ms Sampling E(T) ####
#Get fitted parameters based on bite data
cont_optimize = IntakeModelParams(data = cont_CumulativeIntake_Figure1_cor_long, time = 'EstimatedTime_avg', intake = 'CumulativeGrams_avgBite', fit_fn = FPM_Fit, idVar = 'Figure')
cont_optimize$Emax = c(max(cont_CumulativeIntake_Figure1A_cor$CumulativeGrams_avgBite), max(cont_CumulativeIntake_Figure1B_cor$CumulativeGrams_avgBite),
  max(cont_CumulativeIntake_Figure1C_cor$CumulativeGrams_avgBite), max(cont_CumulativeIntake_Figure1D_cor$CumulativeGrams_avgBite))

cont_optimize

#get bite data from fitted parameters
cont_CumulativeIntake_Figure1A_fit = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(cont_optimize$theta[1], cont_optimize$r[1]), Emax = cont_optimize$Emax[1])
cont_CumulativeIntake_Figure1A_fit$Figure = '1A'
cont_CumulativeIntake_Figure1B_fit = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(cont_optimize$theta[2], cont_optimize$r[2]), Emax = cont_optimize$Emax[2])
cont_CumulativeIntake_Figure1B_fit$Figure = '1B'
cont_CumulativeIntake_Figure1C_fit = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(cont_optimize$theta[3], cont_optimize$r[3]), Emax = cont_optimize$Emax[3])
cont_CumulativeIntake_Figure1C_fit$Figure = '1C'
cont_CumulativeIntake_Figure1D_fit = simBites(time_fn = FPM_Time, nBites = 7200, parameters = c(cont_optimize$theta[4], cont_optimize$r[4]), Emax = cont_optimize$Emax[4])
cont_CumulativeIntake_Figure1D_fit$Figure = '1D'


#correct
cont_CumulativeIntake_Figure1_fit_long = rbind.data.frame(cont_CumulativeIntake_Figure1A_fit, cont_CumulativeIntake_Figure1B_fit, cont_CumulativeIntake_Figure1C_fit, cont_CumulativeIntake_Figure1D_fit)
cont_CumulativeIntake_Figure1_fit_long$Fit = 'FitParams'

#Plot orignial values vs fitted values
cont_CumulativeIntake_Figure1A_fit$Fit = 'FitParams'
cont_CumulativeIntake_Figure1A_cor$Fit = 'Original'

cont_CumulativeIntake_Figure1A_fit.orig = ggplot(cont_CumulativeIntake_Figure1A_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1A_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1A with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

cont_CumulativeIntake_Figure1B_fit$Fit = 'FitParams'
cont_CumulativeIntake_Figure1B_cor$Fit = 'Original'

cont_CumulativeIntake_Figure1B_fit.orig = ggplot(cont_CumulativeIntake_Figure1B_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1B_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1B with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

cont_CumulativeIntake_Figure1C_fit$Fit = 'FitParams'
cont_CumulativeIntake_Figure1C_cor$Fit = 'Original'

cont_CumulativeIntake_Figure1C_fit.orig = ggplot(cont_CumulativeIntake_Figure1C_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1C_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1C with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

cont_CumulativeIntake_Figure1D_fit$Fit = 'FitParams'
cont_CumulativeIntake_Figure1D_cor$Fit = 'Original'

cont_CumulativeIntake_Figure1D_fit.orig = ggplot(cont_CumulativeIntake_Figure1D_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1D_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1D with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())


#Facet wrap
cont_CumulativeIntake_Figure1_cor_long$Fit = 'Original'
cont_CumulativeIntake_Figure1_fit.orig_grid = ggplot(cont_CumulativeIntake_Figure1_fit_long, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1_cor_long, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1A with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~Figure)


#### Try to Recover Parameters from 100 bite Sampling E(T) ####
bite100_optimize = IntakeModelParams(data = bite100_CumulativeIntake_Figure1_cor_long, time = 'EstimatedTime_avg', intake = 'CumulativeGrams_avgBite', fit_fn = FPM_Fit, idVar = 'Figure')
bite100_optimize$Emax = c(max(bite100_CumulativeIntake_Figure1A_cor$CumulativeGrams_avgBite), max(bite100_CumulativeIntake_Figure1B_cor$CumulativeGrams_avgBite),
  max(bite100_CumulativeIntake_Figure1C_cor$CumulativeGrams_avgBite), max(bite100_CumulativeIntake_Figure1D_cor$CumulativeGrams_avgBite))

#get bite data back
bite100_CumulativeIntake_Figure1A_fit = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(bite100_optimize$theta[1], bite100_optimize$r[1]), Emax = bite100_optimize$Emax[1])
bite100_CumulativeIntake_Figure1A_fit$Figure = '1A'
bite100_CumulativeIntake_Figure1B_fit = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(bite100_optimize$theta[2], bite100_optimize$r[2]), Emax = bite100_optimize$Emax[2])
bite100_CumulativeIntake_Figure1B_fit$Figure = '1B'
bite100_CumulativeIntake_Figure1C_fit = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(bite100_optimize$theta[3], bite100_optimize$r[3]), Emax = bite100_optimize$Emax[3])
bite100_CumulativeIntake_Figure1C_fit$Figure = '1C'
bite100_CumulativeIntake_Figure1D_fit = simBites(time_fn = FPM_Time, nBites = 100, parameters = c(bite100_optimize$theta[4], bite100_optimize$r[4]), Emax = bite100_optimize$Emax[4])
bite100_CumulativeIntake_Figure1D_fit$Figure = '1D'

#correct
bite100_CumulativeIntake_Figure1_fit_long = rbind.data.frame(bite100_CumulativeIntake_Figure1A_fit, bite100_CumulativeIntake_Figure1B_fit, bite100_CumulativeIntake_Figure1C_fit, bite100_CumulativeIntake_Figure1D_fit)
bite100_CumulativeIntake_Figure1_fit_long$Fit = 'FitParams'


#individual graphs
bite100_CumulativeIntake_Figure1A_fit$Fit = 'FitParams'
bite100_CumulativeIntake_Figure1A_cor$Fit = 'Original'

bite100_CumulativeIntake_Figure1A_fit.orig = ggplot(bite100_CumulativeIntake_Figure1A_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = bite100_CumulativeIntake_Figure1A_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1A with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite100_CumulativeIntake_Figure1B_fit$Fit = 'FitParams'
bite100_CumulativeIntake_Figure1B_cor$Fit = 'Original'

bite100_CumulativeIntake_Figure1B_fit.orig = ggplot(bite100_CumulativeIntake_Figure1B_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = bite100_CumulativeIntake_Figure1B_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1B with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite100_CumulativeIntake_Figure1C_fit$Fit = 'FitParams'
bite100_CumulativeIntake_Figure1C_cor$Fit = 'Original'

bite100_CumulativeIntake_Figure1C_fit.orig = ggplot(bite100_CumulativeIntake_Figure1C_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = bite100_CumulativeIntake_Figure1C_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1C with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

bite100_CumulativeIntake_Figure1D_fit$Fit = 'FitParams'
bite100_CumulativeIntake_Figure1D_cor$Fit = 'Original'

bite100_CumulativeIntake_Figure1D_fit.orig = ggplot(bite100_CumulativeIntake_Figure1D_fit, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = bite100_CumulativeIntake_Figure1D_cor, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1D with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank())

#facet wrap
bite100_CumulativeIntake_Figure1_cor_long$Fit = 'Original'

bite100_CumulativeIntake_Figure1_fit.orig_grid = ggplot(bite100_CumulativeIntake_Figure1_fit_long, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = bite100_CumulativeIntake_Figure1_cor_long, method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  geom_smooth(method = 'gam', formula = y~s(x), aes(linetype = Fit, color = Fit)) +
  ggtitle('Replication of Figure 1 with 250 ms sampling (for 30 min meal)') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  scale_linetype_manual(values=c(1, 3),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  scale_color_manual(values=c('cornflowerblue', 'black'),
    name = 'Model Parameters',
    breaks = c( 'Original', 'FitParams'),
    labels = c('Original', 'FitParams')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~Figure)

#### Try to Recover Parameters from 33 bite Sampling E(T) ####
bite33_CumulativeIntake_Figure1A_Et$nBiteGroup = ifelse(bite33_CumulativeIntake_Figure1A_Et$nBites < 26, '15-25',
  ifelse(bite33_CumulativeIntake_Figure1A_Et$nBites < 36, '26-35', ifelse(bite33_CumulativeIntake_Figure1A_Et$nBites < 46,'36-45', '45+')))

bite33_CumulativeIntake_Figure1_plotA_fit.orig = ggplot(bite33_CumulativeIntake_Figure1A_Et, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1A_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 3, se = TRUE) +
  geom_jitter(aes(color = nBites)) +
  scale_colour_gradientn(colors = rainbow(10)) +
  #scale_color_viridis_c(option = 'magma') +
  #guides(colour=FALSE) +
  ggtitle('Replication of Figure 1A with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~factor(nBiteGroup))

bite33_CumulativeIntake_Figure1B_Et$nBiteGroup = ifelse(bite33_CumulativeIntake_Figure1B_Et$nBites < 26, '15-25',
  ifelse(bite33_CumulativeIntake_Figure1B_Et$nBites < 36, '26-35', ifelse(bite33_CumulativeIntake_Figure1B_Et$nBites < 46,'36-45', '45+')))

bite33_CumulativeIntake_Figure1_plotB_fit.orig = ggplot(bite33_CumulativeIntake_Figure1B_Et, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1B_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 3, se = TRUE) +
  geom_jitter(aes(color = nBites)) +
  scale_colour_gradientn(colors = rainbow(10)) +
  #scale_color_viridis_c(option = 'magma') +
  #guides(colour=FALSE) +
  ggtitle('Replication of Figure 1B with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~factor(nBiteGroup))

bite33_CumulativeIntake_Figure1C_Et$nBiteGroup = ifelse(bite33_CumulativeIntake_Figure1C_Et$nBites < 26, '15-25',
  ifelse(bite33_CumulativeIntake_Figure1C_Et$nBites < 36, '26-35', ifelse(bite33_CumulativeIntake_Figure1C_Et$nBites < 46,'36-45', '45+')))

bite33_CumulativeIntake_Figure1_plotC_fit.orig = ggplot(bite33_CumulativeIntake_Figure1C_Et, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1C_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 3, se = TRUE) +
  geom_jitter(aes(color = nBites)) +
  scale_colour_gradientn(colors = rainbow(10)) +
  #scale_color_viridis_c(option = 'magma') +
  #guides(colour=FALSE) +
  ggtitle('Replication of Figure 1C with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~factor(nBiteGroup))

bite33_CumulativeIntake_Figure1D_Et$nBiteGroup = ifelse(bite33_CumulativeIntake_Figure1D_Et$nBites < 26, '15-25',
  ifelse(bite33_CumulativeIntake_Figure1D_Et$nBites < 36, '26-35', ifelse(bite33_CumulativeIntake_Figure1D_Et$nBites < 46,'36-45', '45+')))

bite33_CumulativeIntake_Figure1_plotD_fit.orig = ggplot(bite33_CumulativeIntake_Figure1D_Et, aes(y = CumulativeGrams_avgBite, x = EstimatedTime_avg)) +
  stat_smooth(data = cont_CumulativeIntake_Figure1D_cor, method = 'gam', formula = y~s(x), linetype = 1, color = 'black') +
  geom_smooth(method = 'gam', formula = y~s(x), linetype = 3, se = TRUE) +
  geom_jitter(aes(color = nBites)) +
  scale_colour_gradientn(colors = rainbow(10)) +
  #scale_color_viridis_c(option = 'magma') +
  #guides(colour=FALSE) +
  ggtitle('Replication of Figure 1D with 33 Bites') +
  scale_y_continuous(name='E(t)') +
  scale_x_continuous(name='Time (min)') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    panel.background = element_blank()) + facet_wrap(~factor(nBiteGroup))
