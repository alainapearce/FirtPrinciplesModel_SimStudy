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
library(RColorBrewer)
#installed from devtools::install_github("debruine/faux")

#### set up ####

#set piorking directory to location of script--not needed pihen called 
#through Rmarkdopin doc. Uncomment belopi if running locally/manually
# this.dir = getActiveDocumentContext()$path
# setwd(dirname(this.dir))

source('functions.R')

#### Source Setup Scripts ####
source('CurveFit.R')
source('BiteSimulated.R')
source('CumulativeIntakeCurve.R')

#### Generate curves from Figure 1A with bite level data ####
# Figure 1A: theta = 30g/min, r = 0.17 1/min, Emax = 400
# Figure 1B: theta = 30g/min, r = 0.17 1/min, Emax = 600
# Figure 1C: theta = 30g/min, r = 0.5 1/min, Emax = 400
# Figure 1C: theta = 50g/min, r = 0.17 1/min, Emax = 400

source('BiteRateModel_Figure1sim.R')

#### Fogel et al., 2017: A description of an ‘obesogenic’ eating style ####
#### that promotes higher energy intake and is associated with greater ####
#### adiposity in 4.5 year-old children: Results from the GUSTO cohort ####
#### Slow vs Faster Eaters (median = 6.41 g/min) ####
#### N = 386 ####
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

source('Fogel2017_SimDat.R')
SimDat_Fogel2017 = Fogel2017_simDat(500)
SimDat_Fogel2017$Emax = SimDat_Fogel2017$TotalGrams
SimDat_Fogel2017$EatRate_group = ifelse(SimDat_Fogel2017$EatRate <=6.41, 'Slow', 'Fast')

#### Generate parameter distributions from Fogel_simDat ####

SimBites_Fogel2017_bitesize0.25 = CIC_simulate(SimDat_Fogel2017, bitesize_sd = 0.25)





SimBites_Fogel2017_CumulativeIntake_BiteSizeVar = ggplot(SimBites_Fogel2017_bitesize_long, 
                                                      aes(y = EstimatedIntake, x = Time,
                                                          color = BiteSizeEstimate, linetype = AverageYN)) +
  geom_smooth(method = 'loess', formula = y~x, se = F) +
  scale_color_manual(values=c(brewer.pal(8, "Dark2"), 'black'))+
  guides(linetype=FALSE)+
  ggtitle('Cumulative Intake Curves from Simulated Bite Data: Avg Bite Size vs differing levels of variability') +
  scale_y_continuous(name='Estimated Intake (E(t))') +
  scale_x_continuous(name='Time (min)') +
  labs(colour = "BiteSizeEstimate",linetype = "AverageYN") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank())

SimBites_Fogel2017_CumulativeIntake_BiteSizeVar_byIDslow = ggplot(SimBites_Fogel2017_bitesize_long[SimBites_Fogel2017_bitesize_long$EatRate_group == 'Slow', ], 
                                                      aes(y = EstimatedIntake, x = Time,
                                                          color = BiteSizeEstimate, linetype = AverageYN)) +
  geom_smooth(method = 'loess', formula = y~x, se = F) +
  scale_color_manual(values=c(brewer.pal(8, "Dark2"), 'black'))+
  guides(linetype=FALSE)+
  ggtitle('Cumulative Intake Curves from Simulated Bite Data: Avg Bite Size vs differing levels of variability') +
  scale_y_continuous(name='Estimated Intake (E(t))') +
  scale_x_continuous(name='Time (min)') +
  labs(colour = "BiteSizeEstimate",linetype = "AverageYN") + 
  theme(strip.text = element_text(size = 6, margin = margin(.05, 0, .05, 0, "cm")), legend.position="bottom", 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        panel.background = element_blank()) + facet_wrap(~id)

SimBites_Fogel2017_CumulativeIntake_BiteSizeVar_byIDfast = ggplot(SimBites_Fogel2017_bitesize_long[SimBites_Fogel2017_bitesize_long$EatRate_group == 'Fast', ], 
                                                                  aes(y = EstimatedIntake, x = Time,
                                                                      color = BiteSizeEstimate, linetype = AverageYN)) +
  geom_smooth(method = 'loess', formula = y~x, se = F) +
  scale_color_manual(values=c(brewer.pal(8, "Dark2"), 'black'))+
  guides(linetype=FALSE)+
  ggtitle('Cumulative Intake Curves from Simulated Bite Data: Avg Bite Size vs differing levels of variability') +
  scale_y_continuous(name='Estimated Intake (E(t))') +
  scale_x_continuous(name='Time (min)') +
  labs(colour = "BiteSizeEstimate",linetype = "AverageYN") + 
  theme(strip.text = element_text(size = 6, margin = margin(.05, 0, .05, 0, "cm")), legend.position="bottom", 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        panel.background = element_blank()) + facet_wrap(~id)


##run optimization
SimFit_Fogel2017_bitesizeAvg = CIC_fit(SimBites_Fogel2017_bitesize0.25, intake = 'CumulativeGrams_avgBite', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesizeAvg$BiteSizeEstimate = 'Average'
SimFit_Fogel2017_bitesize0.25 = CIC_fit(SimBites_Fogel2017_bitesize0.25, intake = 'CumulativeGrams_0.25', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize0.25$BiteSizeEstimate = 'sd_0.25g'
SimFit_Fogel2017_bitesize0.5 = CIC_fit(SimBites_Fogel2017_bitesize0.5, intake = 'CumulativeGrams_0.5', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize0.5$BiteSizeEstimate = 'sd_0.5g'
SimFit_Fogel2017_bitesize0.75 = CIC_fit(SimBites_Fogel2017_bitesize0.75, intake = 'CumulativeGrams_0.75', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize0.75$BiteSizeEstimate = 'sd_0.75g'
SimFit_Fogel2017_bitesize1 = CIC_fit(SimBites_Fogel2017_bitesize1, intake = 'CumulativeGrams_1', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize1$BiteSizeEstimate = 'sd_1g'
SimFit_Fogel2017_bitesize1.5 = CIC_fit(SimBites_Fogel2017_bitesize1.5, intake = 'CumulativeGrams_1.5', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize1.5$BiteSizeEstimate = 'sd_1.5g'
SimFit_Fogel2017_bitesize2 = CIC_fit(SimBites_Fogel2017_bitesize2, intake = 'CumulativeGrams_2', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize2$BiteSizeEstimate = 'sd_2g'
SimFit_Fogel2017_bitesize2.5 = CIC_fit(SimBites_Fogel2017_bitesize2.5, intake = 'CumulativeGrams_2.5', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize2.5$BiteSizeEstimate = 'sd_2.5g'
SimFit_Fogel2017_bitesize3 = CIC_fit(SimBites_Fogel2017_bitesize3, intake = 'CumulativeGrams_3', parameters = c(20, .1), method = 'Thomas2017')
SimFit_Fogel2017_bitesize3$BiteSizeEstimate = 'sd_3g'

SimFit_Fogel2017_bitesize_long = rbind(SimFit_Fogel2017_bitesize0.25, SimFit_Fogel2017_bitesize0.5, SimFit_Fogel2017_bitesize0.75,
                                       SimFit_Fogel2017_bitesize1, SimFit_Fogel2017_bitesize1.5, SimFit_Fogel2017_bitesize2,
                                       SimFit_Fogel2017_bitesize2.5, SimFit_Fogel2017_bitesize3, SimFit_Fogel2017_bitesizeAvg)
SimFit_Fogel2017_bitesize_long$id = factor(SimFit_Fogel2017_bitesize_long$id)

SimFit_Fogel2017_bitesize_listbyid = t(sapply(levels(SimFit_Fogel2017_bitesize_long$id), function(x){
  SimFit_Fogel2017_bitesize_long[SimFit_Fogel2017_bitesize_long$id == x, ]}))

SimFit_Fogel2017_tabnames = c('id', 'theta', 'r', 'neg2ll', 'counts.function')

SimFit_Fogel2017_means = data.frame(SimFit_Fogel2017_bitesizeAvg$id, matrix(mapply(mean, SimFit_Fogel2017_bitesize_listbyid[, 2:5]), nrow = 40, byrow = FALSE))
SimFit_Fogel2017_sd = data.frame(SimFit_Fogel2017_bitesizeAvg$id, matrix(mapply(sd, SimFit_Fogel2017_bitesize_listbyid[, 2:5]), nrow = 40, byrow = FALSE))
SimFit_Fogel2017_min = data.frame(SimFit_Fogel2017_bitesizeAvg$id, matrix(mapply(range, SimFit_Fogel2017_bitesize_listbyid[, 2:5])[1, ], nrow = 40, byrow = FALSE))
SimFit_Fogel2017_max = data.frame(SimFit_Fogel2017_bitesizeAvg$id, matrix(mapply(range, SimFit_Fogel2017_bitesize_listbyid[, 2:5])[2, ], nrow = 40, byrow = FALSE))

SimFit_Fogel2017_tab = cbind(SimFit_Fogel2017_means[1:2], SimFit_Fogel2017_sd[2], SimFit_Fogel2017_min[2], SimFit_Fogel2017_max[2],
                             SimFit_Fogel2017_means[3], SimFit_Fogel2017_sd[3], SimFit_Fogel2017_min[3], SimFit_Fogel2017_max[3],
                             SimFit_Fogel2017_means[4], SimFit_Fogel2017_sd[4], SimFit_Fogel2017_min[4], SimFit_Fogel2017_max[4],
                             SimFit_Fogel2017_means[5], SimFit_Fogel2017_sd[5], SimFit_Fogel2017_min[5], SimFit_Fogel2017_max[5])
names(SimFit_Fogel2017_tab) = c('id', 'theta_mean', 'theta_sd', 'theta_min', 'theta_max',
                                'r_mean', 'r_sd', 'r_min', 'r_max',
                                'n2ll_mean', 'n2ll_sd', 'n2ll_min', 'n2ll_max',
                                'count_mean', 'count_sd', 'count_min', 'count_max')


SimFit_Fogel2017_tab = merge(SimFit_Fogel2017_tab, SimDat_Fogel2017[c(10, 12)], by = 'id')

SimFit_Fogel2017_tab_slow = SimFit_Fogel2017_tab[SimFit_Fogel2017_tab$EatRate_group == 'Slow', ]
SimFit_Fogel2017_tab_fast = SimFit_Fogel2017_tab[SimFit_Fogel2017_tab$EatRate_group == 'Fast', ]

# ##recover intake curves
# SimFit_Fogel2017_bitesizeAvg_merge = merge(SimFit_Fogel2017_bitesizeAvg, SimDat_Fogel2017, by = 'id')
# SimFit_Fogel2017_bitesizeAvg_merge = merge(SimBites_Fogel2017_bitesize[c(1, 3)], SimFit_Fogel2017_bitesizeAvg_merge, by = 'id')
# 
# EstBites_Fogel2017_bitesizeAvg = CIC_estimateBites(SimFit_Fogel2017_bitesizeAvg_merge, method = 'Thomas2017')
# EstBites_Fogel2017_bitesizeAvg$BiteSizeEstimate = 'Average'
# EstBites_Fogel2017_bitesize0.25 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_0.25g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize0.25$BiteSizeEstimate = 'sd_0.25g'
# EstBites_Fogel2017_bitesize0.5 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_0.5g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize0.5$BiteSizeEstimate = 'sd_0.5g'
# EstBites_Fogel2017_bitesize0.75 =CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_0.75g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize0.75$BiteSizeEstimate = 'sd_0.75g'
# EstBites_Fogel2017_bitesize1 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_1g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize1$BiteSizeEstimate = 'sd_1g'
# EstBites_Fogel2017_bitesize1.5 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_1.5g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize1.5$BiteSizeEstimate = 'sd_1.5g'
# EstBites_Fogel2017_bitesize2 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_2g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize2$BiteSizeEstimate = 'sd_2g'
# EstBites_Fogel2017_bitesize2.5 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_2.5g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize2.5$BiteSizeEstimate = 'sd_2.5g'
# EstBites_Fogel2017_bitesize3 = CIC_estimateBites(SimFit_Fogel2017_bitesize_long_merge[SimFit_Fogel2017_bitesize_long_merge$BiteSizeEstimate == 'sd_3g', ], method = 'Thomas2017')
# EstBites_Fogel2017_bitesize3$BiteSizeEstimate = 'sd_3g'
# 

# 
# SimBites_Fogel2017_bitesize = cbind(SimBites_Fogel2017_bitesize0.25, SimBites_Fogel2017_bitesize0.5[4],
#                                     SimBites_Fogel2017_bitesize0.75[4], SimBites_Fogel2017_bitesize1[4],
#                                     SimBites_Fogel2017_bitesize1.5[4], SimBites_Fogel2017_bitesize2[4],
#                                     SimBites_Fogel2017_bitesize2.5[4], SimBites_Fogel2017_bitesize3[4])
# 
# 
# SimBites_Fogel2017_bitesize_long = melt(SimBites_Fogel2017_bitesize, id.vars = names(SimBites_Fogel2017_bitesize)[c(1:3)])
# SimBites_Fogel2017_bitesize_long$BiteSizeEstimate = ifelse(SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_0.25', 'sd_0.25g', ifelse(
#   SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_0.5', 'sd_0.5g', ifelse(
#     SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_0.75', 'sd_0.75g', ifelse(
#       SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_1', 'sd_1g', ifelse(
#         SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_1.5', 'sd_1.5g', ifelse(
#           SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_2', 'sd_2g', ifelse(
#             SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_2.5', 'sd_2.5g', ifelse(
#               SimBites_Fogel2017_bitesize_long$variable == 'CumulativeGrams_3', 'sd_3g', 'Average'
#             ))))))))
# SimBites_Fogel2017_bitesize_long$BiteSizeEstimate = factor(SimBites_Fogel2017_bitesize_long$BiteSizeEstimate,
#                                                            levels = c('sd_0.25g', 'sd_0.5g', 'sd_0.75g', 
#                                                                       'sd_1g', 'sd_1.5g', 'sd_2g', 'sd_2.5g', 
#                                                                       'sd_3g', 'Average'))
# BiteSizeEstimateLevels = levels(SimBites_Fogel2017_bitesize_long$BiteSizeEstimate)
# SimBites_Fogel2017_bitesize_long$EstimatedIntake = SimBites_Fogel2017_bitesize_long$value
# SimBites_Fogel2017_bitesize_long = SimBites_Fogel2017_bitesize_long[c(1:3, 6:7)]
# SimBites_Fogel2017_bitesize_long$AverageYN = ifelse(SimBites_Fogel2017_bitesize_long$BiteSizeEstimate == 'Average', 'Y', 'N')
# SimBites_Fogel2017_bitesize_long$AverageYN = factor(SimBites_Fogel2017_bitesize_long$AverageYN, levels = c('Y', 'N'))
# 
# SimBites_Fogel2017_bitesize_long = merge(SimBites_Fogel2017_bitesize_long, SimDat_Fogel2017[c(10, 12)], by = 'id')
# 
# 
# SimBites_Fogel2017_CumulativeIntake_BiteSizeVar = ggplot(SimBites_Fogel2017_bitesize_long, 
#                                                          aes(y = EstimatedIntake, x = Time,
#                                                              color = BiteSizeEstimate, linetype = AverageYN)) +
#   geom_smooth(method = 'loess', formula = y~x, se = F) +
#   scale_color_manual(values=c(brewer.pal(8, "Dark2"), 'black'))+
#   guides(linetype=FALSE)+
#   ggtitle('Cumulative Intake Curves from Simulated Bite Data: Avg Bite Size vs differing levels of variability') +
#   scale_y_continuous(name='Estimated Intake (E(t))') +
#   scale_x_continuous(name='Time (min)') +
#   labs(colour = "BiteSizeEstimate",linetype = "AverageYN") + 
#   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
#         panel.background = element_blank())
# 
# SimBites_Fogel2017_CumulativeIntake_BiteSizeVar_byIDslow = ggplot(SimBites_Fogel2017_bitesize_long[SimBites_Fogel2017_bitesize_long$EatRate_group == 'Slow', ], 
#                                                                   aes(y = EstimatedIntake, x = Time,
#                                                                       color = BiteSizeEstimate, linetype = AverageYN)) +
#   geom_smooth(method = 'loess', formula = y~x, se = F) +
#   scale_color_manual(values=c(brewer.pal(8, "Dark2"), 'black'))+
#   guides(linetype=FALSE)+
#   ggtitle('Cumulative Intake Curves from Simulated Bite Data: Avg Bite Size vs differing levels of variability') +
#   scale_y_continuous(name='Estimated Intake (E(t))') +
#   scale_x_continuous(name='Time (min)') +
#   labs(colour = "BiteSizeEstimate",linetype = "AverageYN") + 
#   theme(strip.text = element_text(size = 6, margin = margin(.05, 0, .05, 0, "cm")), legend.position="bottom", 
#         panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
#         panel.background = element_blank()) + facet_wrap(~id)
# 
# SimBites_Fogel2017_CumulativeIntake_BiteSizeVar_byIDfast = ggplot(SimBites_Fogel2017_bitesize_long[SimBites_Fogel2017_bitesize_long$EatRate_group == 'Fast', ], 
#                                                                   aes(y = EstimatedIntake, x = Time,
#                                                                       color = BiteSizeEstimate, linetype = AverageYN)) +
#   geom_smooth(method = 'loess', formula = y~x, se = F) +
#   scale_color_manual(values=c(brewer.pal(8, "Dark2"), 'black'))+
#   guides(linetype=FALSE)+
#   ggtitle('Cumulative Intake Curves from Simulated Bite Data: Avg Bite Size vs differing levels of variability') +
#   scale_y_continuous(name='Estimated Intake (E(t))') +
#   scale_x_continuous(name='Time (min)') +
#   labs(colour = "BiteSizeEstimate",linetype = "AverageYN") + 
#   theme(strip.text = element_text(size = 6, margin = margin(.05, 0, .05, 0, "cm")), legend.position="bottom", 
#         panel.border = element_blank(), panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
#         panel.background = element_blank()) + facet_wrap(~id)
# 
# 
