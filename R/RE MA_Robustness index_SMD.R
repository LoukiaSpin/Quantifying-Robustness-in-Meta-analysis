##########################################################################################
#                                                                                        #
#       Robustness Index of SMD comparing each informative IMDoM scenario with MAR       #
#                                                                                        #
##########################################################################################



## Load libraries
list.of.packages <- c("ggplot2", "ggthemes", "ggpubr", "reshape2")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("31_Functions/KLD & Robustness Index.R")     # Estimate the Robustness Index for each scenario and comparison
source("31_Functions/KLD plots_All scenarios.R")    # Panel of KLD.yx plots for all scenarios of one comparison (Makes sense for MA or for one reference-comparison in NMA)



## Load SMD and sd.SMD (Bayesian RE-MA results)
(SMD <- read.table("./30_Analysis & Results/MA - Continuous outcome/MA-MOD_SMD.txt", head = T)[, 2:3])
drug.names <- c("inositol", "placebo")
nt <- 2



## Index of Robustness  
(RI <- RobustnessIndex(SMD, primary.scenar = 13, nt = nt)$RI)
(KLD <- RobustnessIndex(SMD, primary.scenar = 13, nt = nt)$kldxy)



## Panel of KLD.yx plots for all scenarios in SMD
tiff("./30_Analysis & Results/MA - Continuous outcome/KLD Panel.tiff", height = 30, width = 38, units = "cm", compression = "lzw", res = 300)
KLD.plots(SMD, primary.scenar = 13, compar = 1, "continuous", drug.names)
dev.off()





