
setwd(getwd())

source("globals.R")

setwd(HOME_DIR)
source(UNI_CDE) 
source(DENSITY)


N = 1e4
K = 2

synth_data_bimodal = r_binorm(N)
y      = synth_data_bimodal$y
X      = synth_data_bimodal$X
params = synth_data_bimodal$param_mat
y_grid = seq(-5, 5, length.out = 1000)

x         = c(-0.6745, 0, 0.6745)
f_yx      = matrix(0, nrow = length(y_grid), ncol = length(x))

for(i in 1:length(x)) {
    f_yx[,i] = d_binorm(y_grid = y_grid, x = x[i])
}

df_y      = data.frame(x = y_grid, f_yx)
df_y_long = melt(df_y, id.vars = "x")
levels(df_y_long$variable) = c("x = -0.6745", "x = 0", "x = 0.6745")
ggplot(df_y_long, aes(x, y = value)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(min(y_grid), max(y_grid), by = 1)) + 
    facet_wrap(~variable)


# run cavi ---------------------------------------------------------------------

# (1.1) covariate-INDEPENENT vb
# source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
# theta1_1 = varbdr(y = y, X = X, K = 3)               # store CAVI results

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))               # load cov-dep varbdr.R
theta0 = varbdr(y = y, X = X, K, intercept = TRUE)    # store CAVI results

saveRDS(theta0, file = "bimodal_results/theta_N_1e4_K_2.RDS")


# read in variational objects for varying N ------------------------------------

theta1_500_2 = readRDS(file = "bimodal_results/theta_N_500_K_2.RDS") 
theta1_1e3_2 = readRDS(file = "bimodal_results/theta_N_1e3_K_2.RDS") 
theta1_1e4_2 = readRDS(file = "bimodal_results/theta_N_1e4_K_2.RDS")



# overlay approximations for varying values of N -------------------------------
library(dplyr)
overlayPlots = function(theta, K, den_label,
                        x = c(-0.6745, 0, 0.6745)) {
    
    source(DENSITY)
    
    n = length(theta)
    approx_d  = rep(list(py_bouch), n)
    # approx_d = list(py_bouch, py_bouch, py_bouch, py_bouch)
    # den_label = c("cov-dep-1", "cov-dep-2", "cov-dep-3", "cov-dep-4")
    p_list2 = xQuantileDensity(x, params, d_binorm,
                               approx_d, den_label, theta, K,
                               y_grid)
    
    p = multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), 
                  cols = 3)
    return(p)
}

overlayPlots(theta = list(theta1_500_2, theta1_1e3_2, theta1_1e4_2), 
             K = 2, den_label = c("N=500", "N=1e3", "N=1e4"))


# generate plots ---------------------------------------------------------------

# generate plot overlays of the percentiles
theta2_list      = list(theta1_1e4_2)      # list of var. params for each alg
source(DENSITY)
approx_d  = list(py_bouch)            # list of approx density functions
den_label = c("cov-dep")   # labels for each approx density
# params argument isn't used, function will redefine what it needs
# params as is corresponds to the parameters for the randomly generated x
# here, we care about the parameters associated with the quantiles
x         = c(-0.6745, 0, 0.6745)
y_grid = seq(-5, 5, length.out = 1000)
p_list3 = xQuantileDensity(x, params, d_binorm,
                           approx_d, den_label, theta2_list, K,
                           y_grid)

multiplot(plotlist = as.list(sapply(p_list3, function(x) x[1])), cols = 3)



# compare with np results ------------------------------------------------------
library("np")
X = data.frame(X[,2])
y = data.frame(y)
xy_df = data.frame(y = y, x = X)
names(xy_df) = c("x", "y")

fy.x = npcdens(y ~ x, xy_df)

overlays = vector("list", length(x))
y_grid = seq(-5, 5, length.out = 1000)

for (i in 1:length(x)) {
    y_eval = data.frame(y = y_grid, x = x[i])
    fy_eval = predict(fy.x, newdata = y_eval)
    
    fy_df = data.frame(y = y_grid, np = fy_eval)
    
    # ggplot(fy_df, aes(x = y, y = np)) + geom_line()
    
    # overlay np estimate with the other approximations ------------------------
    
    # x = 0.15
    full_approx_df = merge(p_list2[[i]]$approx_df, fy_df, by = "y")
    approx_df_long = melt(full_approx_df, 
                          measure.vars = c(den_label, "true", "np"))
    
    overlays[[i]] = ggplot(approx_df_long, 
                           aes(x = y, y = value, colour = variable)) + 
        geom_line(size = 0.8) + labs(x = "y", y = "p(y)") + theme_bw() +
        theme(legend.position = "none")
}

multiplot(plotlist = overlays, cols = 3)



