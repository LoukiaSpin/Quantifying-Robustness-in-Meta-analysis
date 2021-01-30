#####################################################################################################################
#                                                                                                                   #
#                Perform random-effects Bayesian network meta-analysis for aggregate binary outcomes                # 
#                                 <Binomial likelihood, logit link, Random Effects>                                 # 
#                                 (Dias et al., 2013 in Appendix  - PMID: 23104435)                                 #    
#          One-stage pattern-mixture model with Informative Missingness Odds Ratio under several scenarios          #
#                                       (Turner et al., 2015 - PMID: 25809313)                                      #
#                                <Hierarchical, intervention-specific prior log IMOR>                               #
#                                                                                                                   #
#####################################################################################################################



## Load libraries
list.of.packages <- c("dplyr", "R2jags", "mcmcplots")
lapply(list.of.packages, require, character.only = TRUE); rm(list.of.packages)



## Load database as dataframe (each row is a trial)
(data <- read.table("./data/19637942_Baker(2009).txt", head = T)) 
(r <- data %>% dplyr::select(starts_with("r")))                           # Observed events per trial-arm
(m <- data %>% dplyr::select(starts_with("m")))                           # Number of missing participant outcome data per trial-arm
(n <- data %>% dplyr::select(starts_with("n") & !ends_with("a")))         # Number randomised per trial-arm
(t <- data %>% dplyr::select(starts_with("t")))                           # Assigned intervention per trial-arm
(nt <- length(table(as.matrix(t))))                                       # Number Of interventions in the network
(ns <- length(r[, 1]))                                                    # Number of trials
(na <- data[, "na"])                                                      # Number of interventions per trial
(ref <- which.max(table(as.matrix(t))))                                   # Reference intervention is the most frequently investigated (here, placebo)
direction <- 0                                                            # 0 negative outcome; 1 positive outcome



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
(LOR <- do.call(rbind,lapply(1:length(logimor[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary[1:(nt*(nt - 1)/2), c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))
(tausq <- do.call(rbind,lapply(1:length(logimor[, 1]), function(i) jagsfit[[i]]$BUGSoutput$summary["tausq", c("50%", "sd", "2.5%", "97.5%", "Rhat", "n.eff")])))

## Save the results as txt using write.table()!

