#----------------------------------------------------------------------------------------------------------------------------
#     R code 1) to perform random-effects Bayesian pairwise meta-analysis for aggregate continuous outcomes 
#                  <Normal likelihood, identity link, Random Effects> in Dias et al., 2013 in Appendix (PMID: 23104435)  
#                  Standardised Mean Difference after extending the aforementioned model
#                  One-stage pattern-mixture model with Informative Missingness Difference of Means under several scenarios 
#                  <Hierarchical, intervention-specific prior IMDoM> in Spineli et al., 2021 (PMID: 33406990)
#
#     R code 2) to calculate the robustness index, and produce the balloon-plots, and bar-plots of KLD measure per scenario
#
#     Author: Loukia Spineli
#     Date: October 2020
#----------------------------------------------------------------------------------------------------------------------------

## Load libraries
list.of.packages <- c("R2jags", "mcmcplots", "ggplot2", "ggthemes", "reshape2", "ggpubr")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)

## Load functions
source("./R/kld.robustness.index.R")               
source("./R/Barplot.Kullback.Leibler.divergence.R") 
source("./R/enhanced.balloon.plot.ES.R")
source("./R/enhanced.balloon.plot.tau2.R")
 

## Load database as dataframe (each row is a trial)
(data  <- read.table("./data/15106232_Taylor(2009).txt", head = T)) 
(y0    <- data[, c("y1", "y2")])                                                                          # Observed mean outcome per trial-arm
(sd0   <- data[, c("sd1", "sd2")])                                                                        # Observed standard deviation per trial-arm
(m     <- data[, c("m1", "m2")])                                                                          # Number of missing participant outcome data per trial-arm
(c     <- data[, c("c1", "c2")])                                                                          # Number of completers per trial-arm
(se0   <- round(sd0/sqrt(c), 2))                                                                          # Observed standard error per trial-arm
(N     <- c + m)                                                                                          # Number randomised per trial-arm
(ns    <- length(y0[, 1]))                                                                                # Number of included trials
(t     <- matrix(rep(c(1, 2), each = ns), nrow = ns, ncol = 2))                                           # Assigned intervention per trial-arm
(sigma <- sqrt(apply((sd0^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - unlist(2))))      # Trial-specific observed pooled standard deviation
drug.names <- c("placebo", "inositol") 


## Selected predictive prior for tau2 - 'Mental health outcome' outcome type with 'pharma vs PBO' intervention-comparison type (Table 3 in PMID: 25304503)
mean.tausq <- -2.99 
sd.tausq   <- 2.16
prec.tausq <- 1/sd.tausq^2
psi.imdom  <- 1


## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
(scenarios <- c(-2, -1, 0, 1, 2))
(imdom     <- as.matrix(cbind(rep(scenarios, each = 5), rep(scenarios, 5)))) # 2nd column refers to the reference intervention (control in MA)

## Prepare parameters for JAGS
jagsfit <- data.jag <- list()


## Calculate time needed for all models
start.time <- Sys.time()

set.seed(123)
memory.limit(size = 40000)

param.jags <- c("SMD", "tausq")

for(i in 1:length(imdom[, 1])){ 
  
  data.jag[[i]] <- list("y.o" = y0, "se.o" = se0, "m" = m, "t" = t, "N" = N, "ns" = ns, "mean.tausq" = mean.tausq, "prec.tausq" = prec.tausq,
                        "sigma" = sigma, "imdom" = imdom[i, ], "psi.imdom" = psi.imdom)
  
  jagsfit[[i]] <- jags(data = data.jag[[i]], parameters.to.save = param.jags, model.file = "./models/RE-MA SMD IMDoM Pattern-mixture HIE-Arm.txt",
                       n.chains = 3, n.iter = 100000, n.burnin = 10000, n.thin = 5, DIC = T, working.directory = getwd())
}
print(jagsfit)

end.time   <- Sys.time()
time.taken <- end.time - start.time
time.taken     


## Check autocorrelation and traceplot in each scenario
jagsfit.mcmc <- as.mcmc(jagsfit[[1]])                   # For instance, for the first scenario
autplot1(jagsfit.mcmc[, c("SMD", "tausq")], chain = 3)  # Another instance
traplot(jagsfit.mcmc, c("SMD", "tausq"))


## Results on model parameters of interest
# 'Rhat' is useful to infer whether the node has converged!
(SMD   <- do.call(rbind,lapply(1:length(imdom[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["SMD", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
(tausq <- do.call(rbind,lapply(1:length(imdom[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tausq", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))

                                                              
## Optional: save the results as txt 
write.table(round(SMD, 4), file = "./MA-MOD_SMD.txt", sep = "\t", quote = F)
write.table(round(tausq, 4), file = "./MA-MOD_tausq.txt", sep = "\t", quote = F)


## Calculate the Robustness Index 
(RI  <- RobustnessIndex(ES.mat = SMD, primary.scenar = 13, nt = length(drug.names))$RI)      # primary analysis (here, MAR) is number 13

## Enhanced balloon-plot for SMD 
(p1  <- BalloonPlot.Sensitivity.ES(ES.mat = SMD, compar = 1, outcome = "continuous", direction = "negative", drug.names = drug.names))

## Calculate the Kullback-Leibler Divergence (KLD) measure 
(KLD <- RobustnessIndex(ES.mat = SMD, primary.scenar = 13, nt = length(drug.names))$kldxy)  # primary analysis (here, MAR) is number 13

## Bar-plot of KLD measure for all scenarios
(p2  <- Barplot.KLD(unlist(KLD), outcome = "continuous", title = "", ylimit = 0.30) + theme(plot.title = element_blank()))

## Bring both plots together
ggarrange(p1, p2, ncol = 2, labels = c("A)", "B)"))

## Enhanced balloon-plot for tau2
extent <- exp(0.049)   # Median of empirically-based prior for 'mental health indicators' and 'pharma vs placebo' [Rhodes et al. 2015 - PMID: 25304503 (Table 3)]
(p3    <- BalloonPlot.Sensitivity.tau2(tau2.mat = tausq, extent = extent, outcome = "continuous", drug.names = drug.names))
