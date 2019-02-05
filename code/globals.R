
## globals.R

# libraries used
library(ggplot2)
library(reshape2)
library(matrixcalc)

# home directory of all the code
HOME_DIR = "~/varbdr/code"                # linux
# HOME_DIR = "C:/Users/chuu/varbdr/code"  # windows

# directories of each of the different algorithms
COV_INDEP = paste(HOME_DIR, "/cov_indep", sep = '')
COV_DEP   = paste(HOME_DIR, "/cov_dep", sep = '')

# file names (mostly common to all algorithms)
INIT_PRIORS      =  "initPriors.R"
INIT_VAR_PARAMS  =  "initVarParams.R"
VARBDR           =  "varbdr.R"
E_STEP           =  "eStep.R"
M_STEP           =  "mStep.R"
ELBO             =  "elbo.R"
MISC_FUNCS       =  "misc.R"
DENSITY          =  "density.R"

# data-generating schemes
BETA_1 = "beta1.R"



