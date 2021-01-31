#----------------------------------------------------------------------------------------------------------------------------
#     R code to illustrate the proposed functions and reproduce the main and supplementary figures of the article
#     Author: Loukia Spineli
#     Date: January 2021
#----------------------------------------------------------------------------------------------------------------------------



## Load necessary libraries
list.of.packages <- c("ggplot2", "ggthemes", "reshape2", "ggpubr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("./R/enhanced.balloon.plot.ES.R")
source("./R/enhanced.balloon.plot.tau2.R")
source("./R/kld.robustness.index.R")               
source("./R/Barplot.Kullback.Leibler.divergence.R")  
source("./R/heatmap.nma.R") 



##################        PAIRWISE META-ANALYSIS EXAMPLE (Figure 2)        ##################

# Load necessary estimated model parameters (Bayesian 'random-effects meta-analysis' results)
(SMD <- read.table("./data/MA-MOD_SMD.txt", head = T)[, 2:5])     # Posterior mean, SD, 2.5% and 97.5% percentile (SMD; Standardised Mean Difference)
(tau2 <- read.table("./data/MA-MOD_tausq.txt", head = T)[, 2:5])  # Posterior median, SD, 2.5% and 97.5% percentile
drug.names <- c("placebo", "inositol") 



# Calculate the Kullback-Leibler Divergence (KLD) measure 
(KLD <- RobustnessIndex(SMD, primary.scenar = 13, nt = length(drug.names))$kldxy)  # MAR scenario is number 13



# Calculate the Robustness Index 
(RI <- RobustnessIndex(SMD, primary.scenar = 13, nt = length(drug.names))$RI)      # MAR scenario is number 13



# Figure 2A: Enhanced balloon-plot for SMD 
(p1 <- BalloonPlot.Sensitivity.ES(SMD, compar = 1, "continuous", "negative", drug.names))



# Figure 2B: Bar-plot of KLD measure for all scenarios
(p2 <- Barplot.KLD(unlist(KLD), outcome = "continuous", title = "", ylimit = 0.30) + theme(plot.title = element_blank()))



# Bring both plots together
ggarrange(p1, p2, ncol = 2, labels = c("A)", "B)"))



# Supplementary Figure S2: Enhanced balloon-plot for tau2
extent <- exp(0.049)   # Median of empirically-based prior for 'mental health indicators' and 'pharma vs placebo' [Rhodes et al. 2015 - PMID: 25304503 (Table 3)]
(p3 <- BalloonPlot.Sensitivity.tau2(tau2, extent, "continuous", drug.names))




##################        NETWORK META-ANALYSIS EXAMPLE (Figures 3 - 5)        ##################

# Load necessary estimated model parameters (Bayesian 'random-effects NMA' results)
(LOR <- read.table("./data/NMA-MOD_LOR.txt", head = T)[, 2:5])
(tau2.common <- read.table("./data/NMA-MOD_tausq.txt", head = T)[, 2:5])
drug.names.NMA <- c("placebo", "budesodine", "budesodine plus formoterol", "fluticasone", 
                    "fluticasone plus salmeterol", "formoterol", "salmeterol", "tiotropium")   # placebo is the reference
drug.names.tau2 <- c("placebo", "active")



# Calculate the Kullback-Leibler Divergence (KLD) measure 
(KLD.NMA <- RobustnessIndex(LOR, primary.scenar = 13, nt = length(drug.names.NMA))$kldxy)



# Calculate the Robustness Index for all possible comparisons
(RI.NMA <- RobustnessIndex(LOR, primary.scenar = 13, nt = length(drug.names.NMA))$RI)



# Figure 3: Heatmap on robustness index per possible comparison
HeatMap.AllComparisons.RI(RI.NMA, drug.names.NMA, 0.28)    # Robustness threshold equals 0.28



# Figure 4: Bar-plot of KLD measure for the four comparisons with lack of robustness 
(p4 <- Barplot.KLD(KLD.NMA[[3]], outcome = "binary", title = "Fluticasone vs placebo", ylimit = 0.35) + theme(axis.title.x = element_blank()))
(p5 <- Barplot.KLD(KLD.NMA[[4]], outcome = "binary", title = "Fluticasone plus salmeterol vs placebo", ylimit = 0.35) + theme(axis.title = element_blank()))
(p6 <- Barplot.KLD(KLD.NMA[[5]], outcome = "binary", title = "Formoterol vs placebo", ylimit = 0.35))
(p7 <- Barplot.KLD(KLD.NMA[[7]], outcome = "binary", title = "Tiotropium vs placebo", ylimit = 0.35) + theme(axis.title.y = element_blank()))



# Bring all plots together (4 by 2 matrix)
ggarrange(p4, p5, p6, p7, ncol = 2, nrow = 2, common.legend = T, legend = "bottom")



# Enhanced balloon-plots for log odds ratio for the four comparisons with lack of robustness 
(p8 <- BalloonPlot.Sensitivity.ES(LOR, compar = 3, "binary", "negative", drug.names.NMA)) 
(p9 <- BalloonPlot.Sensitivity.ES(LOR, compar = 4, "binary", "negative", drug.names.NMA)) 
(p10 <- BalloonPlot.Sensitivity.ES(LOR, compar = 5, "binary", "negative", drug.names.NMA)) 
(p11 <- BalloonPlot.Sensitivity.ES(LOR, compar = 7, "binary", "negative", drug.names.NMA)) 



# Figure 5: Bring all plots together 
ggarrange(p8, p11, p10, p9, ncol = 2, nrow = 2, labels = c("A)", "B)", "C)", "D)"))



# Supplementary Figure S3: Enhanced balloon-plot for tau2
extent.NMA <- summary(rlnorm(1000, -2.06, 1.51))[3]  # Median of empirically-based prior for 'symptoms reflecting the continuation of condition' and 'pharma vs placebo' [Turner et al. 2015 - PMID: 25475839]
(p12 <- BalloonPlot.Sensitivity.tau2(tau2.common, extent.NMA, "binary", drug.names.tau2))


