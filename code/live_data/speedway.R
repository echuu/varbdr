


library(hdrcde)
library(np)

# 
# lane2 : (1318 x 2) lane 2 of the 4-lane California freeway I-880
# flow  : traffic flow in vehicles per lane per hour
# speed : speed in miles per hour

plot(lane2) # plot of traffic speed vs. traffic flow

X = cbind(as.matrix(lane2$flow))
y = lane2$speed
N = nrow(X)
D = ncol(X)
K = 2
intercept = TRUE
max_iter  = 1e3

y_grid = seq(min(y), max(y), length.out = 2000)
x_in = c(1400)

meanParams = generateParams(y, X, N, D, K, intercept = FALSE, max_iter)
(m_k  = meanParams$m_k)
(mu_k = meanParams$mu_k)
theta_R  = varbdr(y = y, X = X, K, intercept = FALSE, m_k = m_k, mu_k = mu_k)    

theta_sf = varbdr_cpp(y, X, N, D, K, intercept, max_iter)
