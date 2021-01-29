###############################################################################################
#                                                                                             #
#                                                                                             #
#       Robustness Index of log OR comparing each informative log IMOR scenario with MAR      #
#                                                                                             #
#                                                                                             #
###############################################################################################



## Load libraries
list.of.packages <- c("ggpubr", "reshape2")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("31_Functions/KLD & Robustness Index.R")     # Estimate the Robustness Index for each scenario and comparison
source("31_Functions/Heatmap D_All comparisons.R")  # Lower triangular heatmap matrix (Makes sense for NMA)
source("31_Functions/KLD plots_All scenarios.R")    # Panel of KLD.xy plots for all scenarios of one comparison (Makes sense for MA or for one reference-comparison in NMA)



## Load database as dataframe (each row is a trial)
(LOR <- read.table("./30_Analysis & Results/NMA - Binary outcome/NMA-MOD_LOR.txt", head = T)[, 2:3])
drug.names <- c("placebo", "budesodine", "budesodine plus formoterol", "fluticasone", 
                "fluticasone plus salmeterol", "formoterol", "salmeterol", "tiotropium") # PBO is the reference
nt <- length(drug.names)



## Lower triangular heatmap matrix - Comparisons are read from the left to the right 
(RI <- RobustnessIndex(LOR, primary.scenar = 13, nt = nt)$RI)
(KLD <- RobustnessIndex(LOR, primary.scenar = 13, nt = nt)$kldxy)

tiff("./30_Analysis & Results/NMA - Binary outcome/Heatmap_D.tiff", height = 30, width = 38, units = "cm", compression = "lzw", res = 300)
HeatMap.AllComparisons.RI(RI, drug.names, 0.28)  
dev.off()



## Panel of KLD.yx plots for all scenarios in comparison-specific LOR 
tiff("./30_Analysis & Results/NMA - Binary outcome/KLD Panel.tiff", height = 20, width = 35, units = "cm", compression = "lzw", res = 300)
KLD.plots(LOR, primary.scenar = 13, compar = 3, "binary", drug.names)
dev.off()






