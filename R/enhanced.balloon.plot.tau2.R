#----------------------------------------------------------------------------------------------------------------------------
#     R code to create the 'enhanced balloon plot' for the between-trial variance under all re-analyses
#         <For pairwise meta-analysis & network meta-analysis (assuming common between-trial variance)>
#     Author: Loukia Spineli
#     Date: December 2020
#----------------------------------------------------------------------------------------------------------------------------


BalloonPlot.Sensitivity.tau2 <- function(tau2.mat, extent, outcome, drug.names){  
  
    
  ## Define the number of the scenarios 
  scenarios <- c(1, 2, 3, 4, 5) 


  ## Each parameter is a matrix with rows referring to the scenarios 
  extent.tau2 <- upper.tau2 <- lower.tau2 <- sd.tau2 <- tau2 <- rep(NA, length(scenarios)^2)  
  # Effect estimate (e.g. posterior median of tau2)
  tau2 <- tau2.mat[, 1]
    
  # Uncertainty around tau2 (e.g. posterior standard deviation of tau2)
  sd.tau2 <- tau2.mat[, 2]
    
  # Lower bound of the 95% (credible or confidence) intrval of tau2
  lower.tau2 <- tau2.mat[, 3]
    
  # Upper bound of the 95% (credible or confidence) intrval of tau2
  upper.tau2 <- tau2.mat[, 4]
    
  # Dummy variable to indicate the extent of tau2
  extent.tau2 <- ifelse(tau2 < extent, "low", "considerable")
    
  
  ## Normalise tau2. We need this to weight the bubbles in the balloon plot (see, geom_point below).
  tau2.normalised <- (tau2 - min(tau2))/(max(tau2) - min(tau2)) 

  
  ## Indicate all combinations of scenarios for the 'active vs control' comparison
  missp <- data.frame(rep(scenarios, each = 5), rep(scenarios, 5));colnames(missp) <- c("active", "control")
  
  
  ## Now, bring all necessary input data in a dataframe to proceed with the creation of the balloon plot (via ggplot2)
  mat <- data.frame(missp, tau2, sd.tau2, extent.tau2); colnames(mat) <- c("active", "control", "value", "sd.value", "extent")

  
  ## Create the proposed balloon plot (separately, for binary and continuous outcomes)
  if(outcome == "binary"){
    
    bubble <- ggplot(mat, aes(x = active, y = control, color = sd.value, label = sprintf("%#.3f", value))) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 1, xmax = 5), ymin = 1, ymax = 5, color = "grey93", fill = "grey93", alpha = 0.1) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 2, xmax = 4), ymin = 2, ymax = 4, color = "grey100", fill = "grey100", alpha = 0.1) +
                geom_point(aes(size = tau2.normalised), stroke = 2, shape = ifelse(extent.tau2 == "low", "circle", "circle plus")) +  
                scale_size(range = c(0, 30)) +
                geom_text(colour = "black", fontface = "bold", size = 5) +
                geom_label(aes(3, 3, label = round(mat[13, 3], 3)), colour = "black", fontface = "bold",  size = 5) +
                scale_color_gradient(low = "deepskyblue", high = "brown1") +
                scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("1/3", "1/2", "1", "2", "3"), position = "bottom", expand = c(0.2, 0)) +
                scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("1/3", "1/2", "1", "2", "3"), expand = c(0.2, 0)) +
                coord_cartesian(ylim = c(1, 5), clip = 'off') +
                labs(x = paste("IMOR scenario in", drug.names[2]), y = paste("IMOR scenario in", drug.names[1]), color = "") +
                guides(shape = F, size = F) + 
                theme_bw() +
                theme(axis.text.x = element_text(size = 14, angle = 360, vjust = 0.8, hjust = 0.5), axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 1),
                      axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, angle = 90, face = "bold"), 
                      legend.position = "bottom", legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"),
                      legend.title = element_text(size = 14, face = "bold"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "grey86"))     
    
  } else {

    bubble <- ggplot(mat, aes(x = active, y = control, color = sd.value, label = sprintf("%#.3f", value))) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 1, xmax = 5), ymin = 1, ymax = 5, color = "grey93", fill = "grey93", alpha = 0.1) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 2, xmax = 4), ymin = 2, ymax = 4, color = "grey100", fill = "grey100", alpha = 0.1) +
                geom_point(aes(size = tau2.normalised), stroke = 2, shape = ifelse(extent.tau2 == "low", "circle", "circle plus")) +
                scale_size(range = c(0, 30)) +
                geom_text(colour = "black", fontface = "bold", size = 5) +
                geom_label(aes(3, 3, label = round(mat[13, 3], 3)), colour = "black", fontface = "bold",  size = 5) +
                scale_color_gradient(low = "deepskyblue", high = "brown1") +
                scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("-2", "-1", "0", "1", "2"), position = "bottom", expand = c(0.2, 0)) +
                scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("-2", "-1", "0", "1", "2"), expand = c(0.2, 0)) +
                coord_cartesian(ylim = c(1, 5), clip = 'off') +
                labs(x = paste("IMDoM scenario in", drug.names[2]), y = paste("IMDoM scenario in", drug.names[1]), color = "") +
                guides(shape = F, size = F) + 
                theme_bw() +
                theme(axis.text.x = element_text(size = 14, angle = 360, vjust = 0.8, hjust = 0.5), axis.text.y = element_text(size = 14, vjust = 0.5, hjust = 1),
                      axis.title.x = element_text(size = 14, face = "bold"), axis.title.y = element_text(size = 14, angle = 90, face = "bold"), 
                      legend.position = "bottom", legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"),
                      legend.title = element_text(size = 14, face = "bold"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "grey86")) 
    
  }
  
  return(bubble)
}
