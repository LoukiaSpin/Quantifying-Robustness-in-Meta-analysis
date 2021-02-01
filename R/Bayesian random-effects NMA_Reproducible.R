#----------------------------------------------------------------------------------------------------------------------------
#     R code 1) to perform random-effects Bayesian network meta-analysis with consistency for aggregate binary outcomes 
#                  <Binomial likelihood, logit link, Random Effects> in Dias et al., 2013 in Appendix (PMID: 23104435)  
#                  One-stage pattern-mixture model with Informative Missingness Odds Ratio under several scenarios 
#                  <Hierarchical, intervention-specific prior log IMOR> in Turner et al., 2015 (PMID: 25809313)
#
#     R code 2) to produce the balloon-plots, heatmap of robustness, and bar-plots of Kullback-Leibler divergence per scenario
#
#     Author: Loukia Spineli
#     Date: October 2020
#----------------------------------------------------------------------------------------------------------------------------



## Load necessary libraries
list.of.packages <- c("dplyr", "R2jags", "mcmcplots", "ggplot2", "ggthemes", "reshape2", "ggpubr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load functions
source("./R/kld.robustness.index.R")
source("./R/heatmap.nma.R") 
source("./R/Barplot.Kullback.Leibler.divergence.R") 
source("./R/enhanced.balloon.plot.ES.R")
source("./R/enhanced.balloon.plot.tau2.R")
               
 

## Load database as dataframe (each row is a trial)
(data <- read.table("./data/19637942_Baker(2009).txt", head = T)) 
(r <- data %>% dplyr::select(starts_with("r")))                                            # Observed events per trial-arm
(m <- data %>% dplyr::select(starts_with("m")))                                            # Number of missing participant outcome data per trial-arm
(n <- data %>% dplyr::select(starts_with("n") & !ends_with("a")))                          # Number randomised per trial-arm
(t <- data %>% dplyr::select(starts_with("t")))                                            # Assigned intervention per trial-arm
(nt <- length(table(as.matrix(t))))                                                        # Number Of interventions in the network
(ns <- length(r[, 1]))                                                                     # Number of trials
(na <- data[, "na"])                                                                       # Number of interventions per trial
(ref <- which.max(table(as.matrix(t))))                                                    # Reference intervention is the most frequently investigated (here, placebo)
direction <- 0                                                                             # 0 negative outcome; 1 positive outcome
drug.names <- c("placebo", "budesodine", "budesodine plus formoterol", "fluticasone", 
                "fluticasone plus salmeterol", "formoterol", "salmeterol", "tiotropium")   
drug.names.tau2 <- c("placebo", "active")



## Selected predictive prior for tau2 - 'Sign reflecting end of condition' outcome type with 'pharma vs PBO' intervention-comparison type
mean.tausq <- -2.06 
sd.tausq <- 1.51
prec.tausq <- 1/sd.tausq^2
prec.logimor <- psi.logimor <- 1



## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
(scenarios <- c(-log(3), -log(2), log(0.9999), log(2), log(3)))
(logimor <- as.matrix(cbind(rep(scenarios, each = 5), rep(scenarios, 5)))) # 2nd column refers to the reference intervention



## Prepare parameters for JAGS
jagsfit <- data.jag <- list()



## Calculate time needed for all models 
start.time <- Sys.time()

set.seed(123)
memory.limit(size = 40000)

for(i in 1:length(logimor[, 1])){ 
  
  data.jag[[i]] <- list("r" = as.matrix(r), "mod" = as.matrix(m), "n" = as.matrix(n), "t" = as.matrix(t), "na" = na, "nt" = nt, "ns" = ns, "ref" = ref, "mean.tausq" = mean.tausq, "prec.tausq" = prec.tausq,
                        "logimor" = logimor[i, ], "prec.logimor" = prec.logimor, "psi.logimor" = psi.logimor, "D" = direction)
  
  param.jags <- c("LOR", "tausq")
  
  jagsfit[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = "./models/Full RE-NMA IMOR Pattern-mixture HIE-Arm.txt",
                       n.chains = 3, n.iter = 100000, n.burnin = 10000, n.thin = 10, DIC = T)
}
print(jagsfit)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken     



## Check autocorrelation  and traceplot in each scenario
jagsfit.mcmc <- as.mcmc(jagsfit[[1]])                          # For instance, for the first scenario
autplot1(jagsfit.mcmc[, c("LOR[2,1]", "tausq")], chain = 3)    # Another instance
traplot(jagsfit.mcmc, c("LOR[2,1]", "tausq"))



## Results on model parameters of interest
# Rhat is useful to infer whether the node converged!
(LOR <- do.call(rbind,lapply(1:length(logimor[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary[1:(nt*(nt - 1)/2), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
(tausq <- do.call(rbind,lapply(1:length(logimor[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tausq", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))

                               
                               
## Optional: Save the results as txt using write.table()
write.table(round(LOR, 4), file = "./NMA-MOD_LOR.txt", sep = "\t", quote = F)
write.table(round(tausq, 4), file = "./NMA-MOD_tausq.txt", sep = "\t", quote = F)



## Calculate the Kullback-Leibler Divergence (KLD) measure 
(KLD <- RobustnessIndex(ES.mat = LOR, primary.scenar = 13, nt = length(drug.names))$kldxy)  # primary analysis (here, MAR) is number 13



## Calculate the Robustness Index for all possible comparisons
(RI <- RobustnessIndex(ES.mat = LOR, primary.scenar = 13, nt = length(drug.names))$RI)      # primary analysis (here, MAR) is number 13



## Heatmap on robustness index per possible comparison
HeatMap.AllComparisons.RI(RI = RI, drug.names = drug.names, threshold = 0.28)               # Robustness threshold equals 0.28



## Bar-plot of KLD measure for the four comparisons with lack of robustness 
# 'KLD[[x]]' is a list of KLD measures for the x-th comparison under all re-analyses.
(p1 <- Barplot.KLD(KLD[[3]], outcome = "binary", title = "Fluticasone vs placebo", ylimit = 0.35) + theme(axis.title.x = element_blank()))
(p2 <- Barplot.KLD(KLD[[4]], outcome = "binary", title = "Fluticasone plus salmeterol vs placebo", ylimit = 0.35) + theme(axis.title = element_blank()))
(p3 <- Barplot.KLD(KLD[[5]], outcome = "binary", title = "Formoterol vs placebo", ylimit = 0.35))
(p4 <- Barplot.KLD(KLD[[7]], outcome = "binary", title = "Tiotropium vs placebo", ylimit = 0.35) + theme(axis.title.y = element_blank()))



## Bring all bar-plots together 
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = T, legend = "bottom")



## Enhanced balloon-plots for log odds ratio for the four comparisons with lack of robustness 
(p5 <- BalloonPlot.Sensitivity.ES(ES.mat = LOR, compar = 3, outcome = "binary", direction = "negative", drug.names = drug.names))  
(p6 <- BalloonPlot.Sensitivity.ES(ES.mat = LOR, compar = 4, outcome = "binary", direction = "negative", drug.names = drug.names))  
(p7 <- BalloonPlot.Sensitivity.ES(ES.mat = LOR, compar = 5, outcome = "binary", direction = "negative", drug.names = drug.names))   
(p8 <- BalloonPlot.Sensitivity.ES(ES.mat = LOR, compar = 7, outcome = "binary", direction = "negative", drug.names = drug.names)) 



## Bring all balloon-plots together 
ggarrange(p5, p8, p7, p6, ncol = 2, nrow = 2, labels = c("A)", "B)", "C)", "D)"))



## Enhanced balloon-plot for tau2 (assumed common in the network)
extent.NMA <- summary(rlnorm(1000, -2.06, 1.51))[3]  # Median of empirically-based prior for 'symptoms reflecting the continuation of condition' and 'pharma vs placebo' [Turner et al. 2015 - PMID: 25475839]
(p12 <- BalloonPlot.Sensitivity.tau2(tau2.mat = tausq, extent = extent, outcome = "binary", drug.names = drug.names.tau2))


