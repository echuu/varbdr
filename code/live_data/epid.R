
library(dplyr)
setwd("C:/Users/chuu/varbdr/code")
source("globals.R")
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation

setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")                              # load debug functions
sourceCpp("getVarParams.cpp")                        # load C++ implementation

setwd("C:/Users/chuu/varbdr/code/live_data")
epi = read.csv("dde_data.csv", stringsAsFactors = FALSE) # 2380 x 18

# exclude observations with age > 45 weeks
epi_sub = epi %>% filter(GEST_DEL <= 45) # 2313 x 18

# Gestational age vs. DDE_A plot
ggplot(epi_sub, aes(x = DDE_A, y = GEST_DEL)) + geom_point() + theme_bw()



normalize <- function(x) { return ((x - min(x)) / (max(x) - min(x))) }


y = epi_sub$GEST_DEL
X = as.matrix(epi_sub$DDE_A)
X_norm = X - mean(X) # center the covariates

# store mean, sd of covariates to scale input x 
shift = attr(X_norm, 'scaled:center')  # mean of covariates
scale = attr(X_norm, 'scaled:scale')   # sd of covariates

# transform the input to be on same scale as design matrix
x_in = unname(quantile(X, c(0.1, 0.6, 0.9, 0.99))) - mean(X)
x_in_0 = (x_in - shift)

N = nrow(X)
D = ncol(X)
K = 2
intercept = FALSE
max_iter = 5e3 # make bigger after it converges

y_grid = seq(min(y) - 2.5, max(y) + 2.5, length.out = 2000)
theta_gest = varbdr_cpp(y, X_norm, N, D, K, intercept, max_iter)
theta_gest$curr

gest_cd = plotCD(theta_gest, K, x_in, y_grid, true_den = NULL, k_den = NULL)
multiplot(plotlist = gest_cd$cd_plots, cols = 2)


