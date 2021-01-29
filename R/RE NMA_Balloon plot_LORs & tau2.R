##############################################################################################################
#                                                                                                            #
#                                                                                                            #
#           Enhanced balloon plot on summary treatment effects for all missingness scenarios & tau2          #
#                        Random-effects meta-analysis on standardised mean difference                        # 
#                                                                                                            #
#                                                                                                            #
##############################################################################################################



## Load necessary libraries
library("ggpubr")



## Load functions
source("31_Functions/BalloonPlot_Sensitivity analysis_ES.R")
source("31_Functions/BalloonPlot_Sensitivity analysis_tau2.R")



## Load database as dataframe (each row is a trial)
(LOR <- read.table("./30_Analysis & Results/NMA - Binary outcome/NMA-MOD_LOR.txt", head = T)[, 2:5])
(tau2 <- read.table("./30_Analysis & Results/NMA - Binary outcome/NMA-MOD_tausq.txt", head = T)[, 2:5])
drug.names <- c("placebo", "budesodine", "budesodine plus formoterol", "fluticasone", 
                "fluticasone plus salmeterol", "formoterol", "salmeterol", "tiotropium") # PBO is the reference
drug.names.tau2 <- c("placebo", "active")



# Balloonplot of LOR for the four comparisons with NO robustness
(p1 <- BalloonPlot.Sensitivity.ES(LOR, compar = 3, "binary", "negative", drug.names))
(p2 <- BalloonPlot.Sensitivity.ES(LOR, compar = 4, "binary", "negative", drug.names))
(p3 <- BalloonPlot.Sensitivity.ES(LOR, compar = 5, "binary", "negative", drug.names))
(p4 <- BalloonPlot.Sensitivity.ES(LOR, compar = 7, "binary", "negative", drug.names))


tiff("./30_Analysis & Results/NMA - Binary outcome/BalloonPlot_LOR.tiff", height = 30, width = 33, units = "cm", compression = "lzw", res = 300)
ggarrange(p1, p4, p3, p2, ncol = 2, nrow = 2, labels = c("A)", "B)", "C)", "D)"))
dev.off()



# Balloonplot of tau2
extent <- summary(rlnorm(1000, -2.06, 1.51))[3] # median

tiff("./30_Analysis & Results/NMA - Binary outcome/BalloonPlot_tau2.tiff", height = 30, width = 33, units = "cm", compression = "lzw", res = 300)
(p5 <- BalloonPlot.Sensitivity.tau2(tau2, extent, "binary", drug.names.tau2))
dev.off()


