


% read in old faithful data

X = csvread("faithful.csv")';

n = size(X)(2); % num of observations
d = size(X)(1); % num of covariates
k = 25;         % num of clusters

% run variational gmm algorithm
[labels, model, elbo] = gm_vb(X, k);