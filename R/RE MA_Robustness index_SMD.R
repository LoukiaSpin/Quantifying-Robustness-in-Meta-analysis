##########################################################################################
#                                                                                        #
#       Robustness Index of SMD comparing each informative IMDoM scenario with MAR       #
#                                                                                        #
##########################################################################################



## Load libraries
list.of.packages <- c("ggplot2", "ggthemes", "ggpubr", "reshape2")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("Quantifying-Robustness-in-Meta-analysis/R/KLD & Robustness Index.R")     # Estimate the Robustness Index for each scenario and comparison
source("Quantifying-Robustness-in-Meta-analysis/R/KLD plots_All scenarios.R")    # Panel of KLD.yx plots for all scenarios of one comparison (Makes sense for MA or for one reference-comparison in NMA)



## Load SMD and sd.SMD (Bayesian RE-MA results)
(SMD <- read.table("Quantifying-Robustness-in-Meta-analysis/data/MA-MOD_SMD.txt", head = T)[, 2:3])
drug.names <- c("inositol", "placebo")
nt <- 2



## Index of Robustness  
(RI <- RobustnessIndex(SMD, primary.scenar = 13, nt = nt)$RI)








