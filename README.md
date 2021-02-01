***

# A novel decision framework for robustness of primary analysis' results

## Description of the repository

The repository offers the typical structure of separate folders for data, models and R (code/scripts), respectively.
* The _data_ folder includes two input text files: 15106232_Taylor(2009).txt and 19637942_Baker(2009).txt;
* The _models_ folder includes two text files, one for the Bayesian random-effects meta-analysis model of continuous outcome and one for the Bayesian random-effects network meta-analysis model of a binary outcome;
* The _R_ folder includes two analysis scripts (Bayesian random-effects MA_Reproducible.R and Bayesian random-effects NMA_Reproducible.R), which source the model and data scripts and perform all analyses, and the five scripts with functions to produce the relevant output; namely, the enhanced balloon-plots, the heatmap of robustness, and the bar-plots with the Kullback-Leibler divergence measure for each scenario analysis.

[JAGS](http://mcmc-jags.sourceforge.net/) must be installed to employ the [R2jags](https://github.com/suyusung/R2jags/issues/) package. After downloading/cloning the repo, the user can use the .Rproj file to source all code.

The next sections briefly illustrate the functions of our novel decision framework for robustness of the primary analysis results with emphasis on the summary treatment effects.

## Output 
Prerequisite R packages: [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggpubr](https://cran.r-project.org/web/packages/ggpubr/) and [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

### Robustness Index 

To calculate the __robustness index__ for the summary treatment effects, we have developed the function `RobustnessIndex()` which has the following syntax:

```r
RobustnessIndex(ES.mat, primary.scenar, nt)
```

#### Explaining the arguments

* ES.mat: The input is a _data.frame_ of the estimated summary treatment effect and its standard error under the primary analysis and alternative analyses__. Under a Bayesian framework, this input corresponds to the output of the R2jags package for that parameter. For instance, if we are interested in K alternative scenarios about the missing participant outcome data (MOD), in addition to primary analysis under the missing at random assumption (recommended starting point), the dataframe will have K+1 rows. In the case of network meta-analysis, the corresponding results refer to all possible pairwise comparisons of interventions in the investigated network under the primary analysis and alternative analyses. For instance, in the case of a triangle (i.e. three possible pairwise comparisons), the dataframe will have a total of (K+1)x3 rows, that is, K+1 rows for each possible comparison.
* primary.scenar: A number to indicate which one of the K+1 total analyses is the __primary analysis__. 
* nt: The __number of investigated interventions__. In the case of pairwise meta-analysis, we have `nt = 2`. In the case of network meta-analysis, `nt` equals the number of interventions in the investigated network.

#### Output of the function

The function `RobustnessIndex()` returns a list of two items:

1. a vector with the __robustness index__ calculated for all pairwise comparisons, and
2. a list of the __Kullback-Leibler Divergence measure__ for all alternative scenarios for every pairwise comparison.

#### Important notes

The robustness index is specific to the pairwise comparison. Therefore, we can calculate only __one robustness index__ for a pairwise meta-analysis, but __as many as the number of possible comparisons__ in the network meta-analysis. For instance, in a network of four interventions, we have six possible comparisons, and hence, we can calculate a total of six robustness indices. 

### Heatmap on Robustness Index 

To create the __heatmap on robustness index__, we have developed the function `HeatMap.AllComparisons.RI()` which has the following syntax:

```r
HeatMap.AllComparisons.RI(RI, drug.names, threshold)
```

#### Explaining the arguments 

* RI: A vector with the __robustness index__ (RI) calculated for all pairwise comparisons. You need to use, first, the `RobustnessIndex()` function. 
* drug.names: A vector with __the names of the interventions__ compared. The interventions should be in the same order that you considered to run pairwise or network meta-analysis. This is important particularly in the case of network meta-analysis so that you do not match the robustness index with the wrong comparisons.
* threshold: This refers to __the threshold of robustness__. Use __0.17 for a continuous outcome__, or __0.28 for a binary outcome__. You may consider more stringent thresholds, if they are clinically plausible.

#### Output of the function

The function `HeatMap.AllComparisons.RI()` returns a lower triangular heatmap matrix which should be read from left to right. Each cell illustrates the robustness index for the corresponding pairwise comparison. A robustness index below the threshold implies present robustness for the corresponding comparisons (__green__ cells), whereas a robustness index equal or above the threshold implies lack of robustness (__red__ cells). 

### Enhanced Balloon Plot

The function `BalloonPlot.Sensitivity.ES()` creates the __enhanced balloon plot__ for the summary treatment effects and has the following syntax:

```r
BalloonPlot.Sensitivity.ES(ES.mat, compar, outcome, direction, drug.names)
```
#### Explaining the arguments

* ES.mat: the same with the function `RobustnessIndex()`.
* compar: A number that indicates the __comparison of interest__. Evidently, for a pairwise meta-analysis, `compar=1`. In the case of a network meta-analysis, the user should specify the number that corresponds to the comparison of interest. To find the number that corresponds to the comparison of interest, you should use the vector with the names of the interventions (i.e. the argument __drug.names__) to create a matrix of the unique possible pairwise comparisons of T interventions in the investigated network:

```r
comparison <- matrix(combn(drug.names, 2), nrow = length(combn(drug.names, 2))/2, ncol = 2, byrow = T)
```

* direction: A character to specify the __type of outcome__. Use `direction="positive"` for a beneficial outcome, and`direction="negative"` for a harmful outcome.
* drug.names: the same with the function `HeatMap.AllComparisons.RI()`.

#### Important notes 

We provide a separate function to create the __enhanced ballon plot for the between-trial variance__, `BalloonPlot.Sensitivity.tau2()`, and it has the same arguments with the `BalloonPlot.Sensitivity.ES()` function, except for dropping the argument __compar__ and replacing the argument __direction__ with the argument __extent__ to indicate whether there is low statistical heterogeneity (i.e. posterior median of the between-trial variance < median of the selected empirically-based prior distribution for that parameter) or considerable.

The __enhanced balloon plot__ is relevant for every important parameter of the model. For a pairwise meta-analysis, the summary treatment effect and the between-trial variance (under a random-effects model) are the parameters of interest. For a network meta-analysis, the inconsistency factor from the node-splitting approach, and the surface under the cumulative ranking curve are additinal parameters of interest. However, the current functions are only applicable for the summary treatment effects and between-trial variance (assumed common in network meta-analysis). In the next version of the function for the enhanced balloon plot, we will consider other model parameters, as well.

<!--- Last but not least, we are working on bringing all the functions into an R package that will also allow for the joint synthesis of observed and missing outcomes in each trial's arm. -->

