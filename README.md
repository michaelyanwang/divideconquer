# divideconquer
Package for fast sparse Cox proportional hazards model used in Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "A Fast Divide-and-Conquer Sparse Cox Regression." arXiv preprint arXiv:1804.00735 (2018).

The dcalasso package is an updated package for the same method, with a more user-friendly wrapper.

## Simulated dataset (Time-independent survival data)
set.seed(1)

\# Sample size = 1M
N = 1e6

\# Number of covariate = 50
p.x = 50; 

\# Divide the data in 100 chunks for computation
K = 100; n = N/K;  

\# Correlation between covariates
cor = 0.5;  

\# Grid of tuning parameter for adaptive LASSO
lambda.grid = 10^seq(-10,3,0.01); 

\# True parameters of covariates
bb = c(rep(0.8, 3), rep(0.4, 3), rep(0.2, 3))
beta0 = c(1, bb, rep(0, p.x - length(bb)))



\# Generate Data
dat.mat0 = SIM.FUN(N, p.x = p.x, cor = cor, family='Cox',beta0 = beta0)

\# Divide data
dat.list = lapply(1:K,function(kk){dat.mat0[1:n+(kk-1)*n,]}); 
gc()
