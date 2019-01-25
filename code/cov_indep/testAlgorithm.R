
## testAlgorithm.R -- run vb algorithm for covariate-independent case

source("varbdr.R")


## load/generate the data

# response
y = 0

# design matrix
X = 0

# run algorithm
theta = varbdr(y = y, X = X)





# evalulate performance (?); this part still a little unsure of what to do
# do i need to have a predictive density ready? doesn't really make sense to 
# look at coefficients because we introduce those artificially as part of
# the model.


