
# test the np package for conditional density estimation
library("datasets")
library("np")

data("faithful")
attach(faithful)

bw = npcdensbw(formula = eruptions ~ waiting)
summary(bw)

fhat = npcdens(bws = bw)

plot(fhat)


# ------------------------------------------------------------------------------


source("globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)

## generate data from true density ---------------------------------------------
N = 500
K = 4
dp1_synth = r_dpmix1(N)
y0        = dp1_synth$y
X0        = dp1_synth$X
mu1       = dp1_synth$param_vec
y_grid    = seq(-1.5, 1.5, length.out = 500)


X = data.frame(X0[,2])
y = data.frame(y0)
bw = npcdensbw(xdat = X, ydat = y)


data.frame()

# compute condensity object
fhat = npcdens(bws = bw) 

# contains the estimated conditional density function, fhat$condens

summary(fhat)

plot(bw)

