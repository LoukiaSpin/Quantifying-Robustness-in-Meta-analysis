######################################################################################################################
#                                                                                                                    #
#                  Perform random-effects Bayesian meta-analysis for aggregate continuous outcomes                   # 
#                                 <Normal likelihood, identity link, Random Effects>                                 # 
#                                 (Dias et al., 2013 in Appendix  - PMID: 23104435)                                  #    
#      One-stage pattern-mixture model with Informative Missingness Difference of Means under several scenarios      #
#                                       (Spineli et al., 2020 - PMID: 33406990)                                      #
#                                 <Hierarchical, intervention-specific prior IMDoM>                                  #
#                                                                                                                    #
######################################################################################################################



## Load libraries
list.of.packages <- c("R2jags", "mcmcplots")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load database as dataframe (each row is a trial)
(data <- read.table("./data/15106232_Taylor(2009).txt", head = T)) 
(y0 <- data[, c("y1", "y2")])                                          # Observed mean outcome per trial-arm
(sd0 <- data[, c("sd1", "sd2")])                                       # Observed standard deviation per trial-arm
(m <- data[, c("m1", "m2")])                                           # Number of missing participant outcome data per trial-arm
(c <- data[, c("c1", "c2")])                                           # Number of completers per trial-arm
(se0 <- round(sd0/sqrt(c), 2))                                         # Observed standard error per trial-arm
(N <- c + m)                                                           # Number randomised per trial-arm
(ns <- length(y0[, 1]))                                                # Number of included trials
(t <- matrix(rep(c(1, 2), each = ns), nrow = ns, ncol = 2))            # Assigned intervention per trial-arm
(sigma <- sqrt(apply((sd0^2)*(c - 1), 1, sum, na.rm = T)/(apply(c, 1, sum, na.rm = T) - unlist(2))))  # Trial-specific observed pooled standard deviation



## Selected predictive prior for tau2 - 'Mental health outcome' outcome type with 'pharma vs PBO' intervention-comparison type (Table 3 in PMID: 25304503)
mean.tausq <- -2.99 
sd.tausq <- 2.16
prec.tausq <- 1/sd.tausq^2
psi.imdom <- 1



## A 2x2 matrix of 25 reference-specific scenarios (PMID: 30223064)
(scenarios <- c(-2, -1, 0, 1, 2))
(imdom <- as.matrix(cbind(rep(scenarios, each = 5), rep(scenarios, 5)))) # 2nd column refers to the reference intervention (control in MA)



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

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken     



## Check autocorrelation and traceplot in each scenario
jagsfit.mcmc <- as.mcmc(jagsfit[[1]])                   # For instance, for the first scenario
autplot1(jagsfit.mcmc[, c("SMD", "tausq")], chain = 3)  # Another instance
traplot(jagsfit.mcmc, c("SMD", "tausq"))



## Results on model parameters of interest
(SMD <- do.call(rbind,lapply(1:length(imdom[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["SMD", c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
(tausq <- do.call(rbind,lapply(1:length(imdom[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tausq", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))

## Save the results as txt using write.table()!


