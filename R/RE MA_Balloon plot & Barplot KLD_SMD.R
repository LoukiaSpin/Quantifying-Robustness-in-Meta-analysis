##############################################################################################################
#                                                                                                            #
#                                                                                                            #
#              Enhanced balloon plot on summary treatment effects for all missingness scenarios              #
#                   Random-effects network meta-analysis on odds ratio (logarithmic scale)                   # 
#                                                                                                            #
#                                                                                                            #
##############################################################################################################



## Load functions
source("31_Functions/BalloonPlot_Sensitivity analysis_ES.R")
source("31_Functions/BalloonPlot_Sensitivity analysis_tau2.R")
source("31_Functions/KLD & Robustness Index.R")               # Estimate the Robustness Index for each scenario and comparison
source("31_Functions/Barplot KLD_All scenarios.R")          



## Load necessary estimated model parameters (Bayesian RE-MA results)
(SMD <- read.table("./30_Analysis & Results/MA - Continuous outcome/MA-MOD_SMD.txt", head = T)[, 2:5])
(tau2 <- read.table("./30_Analysis & Results/MA - Continuous outcome/MA-MOD_tausq.txt", head = T)[, 2:5])
drug.names <- c("placebo", "inositol") 



## Calculate the KLD measure 
(KLD <- RobustnessIndex(SMD, primary.scenar = 13, nt = length(drug.names))$kldxy)
summary(unlist(KLD))



## Balloonplot of SMD for all possible scenarios about IMDoM in both arms
(p1 <- BalloonPlot.Sensitivity.ES(SMD, compar = 1, "continuous", "negative", drug.names))



## Barplot of KLD measure for all scenarios
(p2 <- Barplot.KLD(unlist(KLD), outcome = "continuous", title = "") + theme(plot.title = element_blank()))



## Bring both plots together
tiff("./30_Analysis & Results/MA - Continuous outcome/BalloonPlot & Barplot KLD.tiff", height = 20, width = 35, units = "cm", compression = "lzw", res = 300)
ggarrange(p1, p2, ncol = 2, labels = c("A)", "B)"))
dev.off()



# Balloonplot of tau2
extent <- exp(0.049)   # Rhodes et al. 2015 - PMID: 25304503 (Table 3)

tiff("./30_Analysis & Results/MA - Continuous outcome/BalloonPlot_tau2.tiff", height = 30, width = 33, units = "cm", compression = "lzw", res = 300)
(p3 <- BalloonPlot.Sensitivity.tau2(tau2, extent, "continuous", drug.names))
dev.off()


