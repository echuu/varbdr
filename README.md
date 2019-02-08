# Variational Bayes for Conditional Density Estimation

This project aims to use a Variational Bayes algorithm for fast inference in the context of Bayesian density regression (conditional density estimation). We use a mixture of experts with covariate dependent weights that enter the mixing weights through a logit link.


### Getting Started and Installation

In order to run all the simulations/diagonistics/analysis, you will need to to install the following libraries:

```
install.packages("ggplot2")
install.packages("reshape2")
install.packages("matrixcalc")
```

In order to make sure everything runs as expected, you need to set the `HOME_DIR` variable in `globals.R` to the directory `varbdr/`. Once `globals.R` is sourced, then all the simulation files should run without problems. 


