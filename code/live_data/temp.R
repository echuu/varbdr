

library(dplyr)
setwd("C:/Users/chuu/varbdr/code")
source("plotCD.R")
source("globals.R")
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation

setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")                              # load debug functions
sourceCpp("getVarParams.cpp")                        # load C++ implementation

library(hdrcde)

X0 = maxtemp[1:3649] # 3649 x 1
y = maxtemp[2:3650] # 3649 x 1


xy_data = data.frame(y = y, x = X0)

ggplot(xy_data, aes(x, y)) + geom_point() + theme_bw() + 
    labs(x = "yesterday's temp", y = "today's temp")


X = cbind(1, as.matrix(X0) / 10)
N = nrow(X)
D = ncol(X)
K = 2
intercept = TRUE
max_iter = 1e4 # make bigger after it converges


theta_wthr = varbdr_cpp(y, X, N, D, K, intercept, max_iter)
theta_wthr$curr # 2302 iterations


x_in = quantile(X[,2], c(0.20, 0.5, 0.75, 0.99))
y_grid = seq(min(y) - 4, max(y) + 10, length.out = 2000)
wthr_cd = plotCD(theta_wthr, K, x_in, y_grid, true_den = NULL, k_den = NULL)
multiplot(plotlist = wthr_cd$cd_plots, cols = 2)



