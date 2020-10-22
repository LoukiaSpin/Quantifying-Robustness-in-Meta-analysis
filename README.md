***

# A novel decision framework for robustness 

## Description

In our article __"Quantifying the robustness of primary analysis results: a case study on missing outcome data in pairwise and network meta-analysis"__ we propose the robustness index alongside a novel decision framework to infer the robustness of the primary analysis results. 

The __robustness index__ quantifies the overall deviation of the primary analysis from alternative analyses regarding the summary treatment effect of a pairwise comparison of interventions. Being a function of the __Kullback-Leibler Divergence__ measure, the robustness index ranges from zero to infinity, with zero indicating a perfect match of the summary treatment effect under the primary analysis and _all_ alternative analyses. In a sense, the further away the robustness index is from zero, the more difficult it is to defend the robustness of the primary analysis results to alternative analyses. 

To conclude on the presence or lack of robustness of the primary analysis results, we compare the robustness index with the __threshold of robustness__ that we have specified using an intuitive and objective rule. 

## Installation 

You can download the __package__ with the proposed functions directly from [GitHub](https://github.com/LoukiaSpin/Quantifying-Robustness-in-Meta-analysis/) or use the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```r
install.packages("devtools")
devtools::install_github("LoukiaSpin/Quantifying-Robustness-in-Meta-analysis", build_vignettes = T)
```
The next sections briefly illustrate the functions of our novel decision framework for robustness of the primary analysis results with emphasis on the summary treatment effects: the `RobustnessIndex()`, and the `HeatMap.AllComparisons.RI()`.

## Robustness Index 

To calculate the __robustness index__, we have developed the function `RobustnessIndex()` which has the following syntax:

```r
RobustnessIndex(ES.mat, primary.scenar, nt)
```

### Explaining the arguments

* <span style="color: green;">ES.mat</span>: It is a dataframe with the <span style="color: blue;">pairwise meta-analysis</span> results on the __summary treatment effect__ (_first column_) and __its standard error__ (_second column_) __under the primary analysis and alternative analyses__. Under a Bayesian framework, the corresponding results should be the posterior mean and posterior standard deviation of the summary treatment effect. For instance, if we are interested in <font size="1">$K$</font>  alternative scenarios about the missing participant outcome data (MOD), in addition to primary analysis under the missing at random assumption (recommended starting point), the dataframe will have <font size="1">$K+1$</font> rows. In the case of <span style="color: blue;">network meta-analysis</span>, the corresponding results refer to all possible pairwise comparisons of interventions in the investigated network under the primary analysis and alternative analyses. For instance, in the case of a triangle (i.e. three possible pairwise comparisons), the dataframe will have a total of <font size="1">$(K+1)\times3$</font> rows, that is, <font size="1">$K+1$</font> rows for each possible comparison.
* <span style="color: green;">primary.scenar</span>: A number to indicate which one of the <font size="1">$K+1$</font> total analyses is the __primary analysis__. 
* <span style="color: green;">nt</span>: The __number of investigated interventions__. In the case of <span style="color: blue;">pairwise meta-analysis</span>, we have `nt = 2`. In the case of <span style="color: blue;">network meta-analysis</span>, `nt` equals the number of interventions in the investigated network.

### Output of the function

The function  `RobustnessIndex()` returns a list of two items:

1. a vector with the __robustness index__ calculated for all pairwise comparisons, and
2. a list of the __Kullback-Leibler Divergence measure__ for all alternative scenarios for every pairwise comparison.

### Important notes

The robustness index is specific to the pairwise comparison. Therefore, we can calculate only __one robustness index__ for a <span style="color: blue;">pairwise meta-analysis</span>, but __as many as the number of possible comparisons__ in the <span style="color: blue;">network</span>. For instance, in a network of four interventions, we have six possible comparisons, and hence, we can calculate a total of six robustness indices. 

## Heatmap on Robustness Index 

We propose the heatmap as a useful tool to illustrate the robustness index for all possible comparisons of interventions in the investigated network. We use green and red colours to indicate the comparisons with presence and lack of robustness, respectively. Make sure that you have already installed the R package [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) to be able to use the proposed function.
To create the __heatmap on robustness index__, we have developed the function `HeatMap.AllComparisons.RI()` which has the following syntax:

```r
HeatMap.AllComparisons.RI(RI, drug.names, threshold)
```

### Explaining the arguments 

* <span style="color: green;">RI</span>: A vector with the __robustness index__ (RI) calculated for all pairwise comparisons. You need to use, first, the `RobustnessIndex()` function. 
* <span style="color: green;">drug.names</span>: A vector with __the names of the interventions__ compared. The interventions should be in the same order that you considered to run pairwise or network meta-analysis. This is important particularly in the case of network meta-analysis so that you do not match the robustness index with the wrong comparisons.
* <span style="color: green;">threshold</span>: This refers to __the threshold of robustness__. Use __0.17 for a continuous outcome__, or __0.28 for a binary outcome__. You may consider more stringent thresholds, if you wish; however, <span style="color: red;">do not exceed the aforementioned 'default' values</span>, as they represent substantial deviations from the primary analysis, and <span style="color: red;">they may signal spurious robustness</span>.

### Output of the function

The function  `HeatMap.AllComparisons.RI()` returns a lower triangular heatmap matrix which should be read from left to right. Each cell illustrates the robustness index for the corresponding pairwise comparison. A robustness index below the threshold implies present robustness for the corresponding comparisons (<span style="background: green;">__green__</span> cells), whereas a robustness index equal or above the threshold implies lack of robustness (<span style="background: red;">__red__</span> cells).

### Important notes

The outcome-specific thresholds mentioned above refer to the corresponding square root of __the median of the empirically-based priors for the between-trial variance__ [1, 2]. We consider this median to reflect __low statistical heterogeneity__, and hence, an _acceptable_ overall deviation of the alternative analyses from the primary analysis. 

#### References> 

1. Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive distributions for between-study heterogeneity and simple methods for their application in Bayesian meta-analysis. Stat Med. 2015; 34(6):984–98. 
2. Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed for the extent of heterogeneity in meta-analyses of continuous outcome data. J Clin Epidemiol. 2015; 68(1):52–60. 

# Comprehensive sensitivity analysis for missing outcome data 

## Description

We propose the __enhanced balloon plot__ as a means to illustrate the parameter of interest (e.g. summary treatment effect), and uncertainty thereof, under progressively stringent yet clinically plausible scenarios for the missingness mechanisms. 

In our article, we have applied the __one-stage pattern-mixture model using Bayesian methods__ [1, 2] to obtain the results for _all scenarios_ in a single analysis. We opted for Bayesian methods for offering flexibility in the analysis of aggregate missing outcome data and for being popular in network meta-analysis. 

The proposed plot can be also obtained straightforwardly when a __two-stage pattern-mixture model__ has been applied to address the missing outcome data in a pairwise or network meta-analysis [3, 4].

#### References

1. Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for uncertainty due to missing binary outcome data in pairwise meta-analysis. Stat Med. 2015; 34(12):2062–80.
2. Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. BMC Med Res Methodol. 2019; 19(1):86.
3. White I, Higgins JPT, Wood AM. Allowing for uncertainty due to missing data in meta-analysis—Part 1: two-stage methods. Stat Med. 2008; 27(5):711–727.
4. Mavridis D, White I, Higgins J, Cipriani A, Salanti G. Allowing for uncertainty due to missing continuous outcome data in pairwise and network meta-analysis. Stat Med. 2015; 34(5):721–41.

## Enhanced balloon plot

The function to create the enhanced balloon plot is `BalloonPlot.Sensitivity.ES()`, and it is also found in the __package__ [Quantifying-Robustness-in-Meta-analysis](https://github.com/LoukiaSpin/Quantifying-Robustness-in-Meta-analysis/). 

Make sure that you have already installed the R package [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) to be able to use the proposed function. The function  `BalloonPlot.Sensitivity.ES()` has the following syntax:

```r
BalloonPlot.Sensitivity.ES(ES.mat, compar, outcome, direction, drug.names)
```

### Explaining the arguments

* <span style="color: green;">ES.mat</span>: It is a dataframe with the <span style="color: blue;">pairwise meta-analysis</span> results on the __summary treatment effect__ (_first column_), __its standard error__ (_second column_), the __lower bound__ (_third column_), and the __upper bound__ (_fourth column_) of the 95% confidence interval __under the primary analysis and alternative analyses__. Under a Bayesian framework, the corresponding results should be the posterior mean, posterior standard deviation and the bounds of the 95% credible interval of the summary treatment effect. For instance, if we are interested in <font size="1">$K$</font>  alternative scenarios about the missing participant outcome data (MOD), in addition to primary analysis under the missing at random assumption (recommended starting point), the dataframe will have <font size="1">$K+1$</font> rows. In the case of <span style="color: blue;">network meta-analysis</span>, the corresponding results refer to all possible pairwise comparisons of interventions in the investigated network under the primary analysis and alternative analyses. For instance, in the case of a triangle (i.e. three possible pairwise comparisons), the dataframe will have a total of <font size="1">$(K+1)\times3$</font> rows, that is, <font size="1">$K+1$</font> rows for each possible comparison.
* <span style="color: green;">compar</span>: A number that indicates the __comparison of interest__. Evidently, for a <span style="color: blue;">pairwise meta-analysis</span>, `compar=1`. In the case of a <span style="color: blue;">network meta-analysis</span>, the user should specify the number that corresponds to the comparison of interest. To find the number that corresponds to the comparison of interest, you should use the vector with the names of the interventions (i.e. the argument <span style="color: green;">drug.names</span>) to create a <font size="1">${T \choose 2}\times2$</font> matrix of the unique possible pairwise comparisons of <font size="1">$T$</font> interventions in the investigated network:

```r
comparison <- matrix(combn(drug.names, 2), nrow = length(combn(drug.names, 2))/2, ncol = 2, byrow = T)
```
* <span style="color: green;">direction</span>: A character to specify the __type of outcome__. Use `direction="positive"` for a beneficial outcome, and`direction="negative"` for a harmful outcome.
* <span style="color: green;">drug.names</span>: A vector with __the names of the interventions__ compared. The interventions should be in the same order that you considered to run pairwise or network meta-analysis. This is important particularly in the case of network meta-analysis so that you do not match the balloons with the wrong comparisons.

### Important notes 

The __enhanced balloon plot__ is relevant for every important parameter of the model. For a <span style="color: blue;">pairwise meta-analysis</span>, the summary treatment effect and the between-trial variance (under a random-effects model) are the parameters of interest. For a <span style="color: blue;">network meta-analysis</span>, the inconsistency factor from the node-splitting approach [1], and the surface under the cumulative ranking curve [2] are additinal parameters of interest. However, <span style="color: red;">the current function is only applicable for the summary treatment effects</span>. In the next version of the `BalloonPlot.Sensitivity.ES()` we will consider other model parameters, as well.

#### References

1. Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. Stat Med. 2010; 29(7-8):932-44. 
2. Salanti G, Ades AE, Ioannidis JPA. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. J Clin Epidemiol. 2011; 64(2):163–71
