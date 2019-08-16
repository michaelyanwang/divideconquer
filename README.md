# divideconquer
Package for fast sparse Cox proportional hazards model used in Wang, Yan, Chuan Hong, Nathan Palmer, Qian Di, Joel Schwartz, Isaac Kohane, and Tianxi Cai. "A Fast Divide-and-Conquer Sparse Cox Regression." arXiv preprint arXiv:1804.00735 (2018).

Please use the dcalasso package (michaelyanwang/dcalasso), which is an updated package for the same method, with a more user-friendly wrapper.
# Installation

```
require(devtools)
install_github("michaelyanwang/divideconquer")
require(divideconquer)
```

## Simulated dataset: Time-independent survival data
```
require(divideconquer)
set.seed(1)

# Sample size = 1M
N = 1e6

# Number of covariate = 50
p.x = 50; 

# Divide the data in 100 chunks for computation
K = 100; n = N/K;  

# Correlation between covariates
cor = 0.5;  

# Grid of tuning parameter for adaptive LASSO
lambda.grid = 10^seq(-10,3,0.01); 

# True parameters of covariates
bb = c(rep(0.8, 3), rep(0.4, 3), rep(0.2, 3))
beta0 = c(1, bb, rep(0, p.x - length(bb)))

# Generate Data
dat.mat0 = SIM.FUN(N, p.x = p.x, cor = cor, family='Cox',beta0 = beta0)

# Divide data
dat.list = lapply(1:K,function(kk){dat.mat0[1:n+(kk-1)*n,]}); 
```

## Analysis of simulated time-independent survival data
```
### Round 1: divide-and-conquer initial estimator
  ptm = proc.time()
  bini = coxph(Surv(time = dat.list[[1]][,1],event = dat.list[[1]][,2])~dat.list[[1]][,-c(1,2)])$coef
  update1 = iteration.fun.cox2(dat.list=dat.list,bini=bini,kk.list=1:K);
  bnew1 = apply(update1$b.k,1,mean);
  time1.0 = proc.time() - ptm;gc();print(time1.0[3])
  
### Round 2: one-step divide-and-conquer update 
  ptm = proc.time()
  update2 = iteration.fun.cox2(dat.list=dat.list,bini=bnew1,kk.list=1:K);
  bnew2 = apply(update2$b.k,1,mean)
  time2.0 = proc.time() - ptm;gc();print(time2.0[3])
  
### Adaptive LASSO estimation
  ptm = proc.time() 
  info2  = -update2$Ahat
  Ahalf2 = svd(-update2$Ahat); Ahalf2 = Ahalf2$u%*%diag(sqrt(Ahalf2$d))%*%t(Ahalf2$v);
  bhat2.DCOS = Est.ALASSO.Approx.GLMNET(Ahalf2%*%bnew2,Ahalf2,bnew2,N.adj=sum(dat.mat0[,2]),lambda.grid=lambda.grid,N.inflate = dim(dat.mat0)[1])
  time2.1 = proc.time() - ptm;gc();print(time2.1[3])
```

## View results
```
# Point estimate
bhat2.DCOS$bhat.BIC # Using BIC to tune adaptive LASSO parameter
bhat2.DCOS$lambda.BIC # Optimal lambda

# Variance-covariance matrix for active set
solve(info2[bhat2.DCOS$bhat.BIC!=0,bhat2.DCOS$bhat.BIC!=0]*n)

# Lower bound of 95% CI for active set
(bhat2.DCOS$bhat.BIC - qnorm(0.975) * sqrt(diag(solve(info2*n))))[bhat2.DCOS$bhat.BIC!=0]
# Upper bound of 95% CI for active set
(bhat2.DCOS$bhat.BIC + qnorm(0.975) * sqrt(diag(solve(info2*n))))[bhat2.DCOS$bhat.BIC!=0]
```

## Simulated dataset: Time-dependent survival data
```
require(divideconquer)
set.seed(1)


# Number of patients
n.subject = 1e6

# Number of time-independent covariates
p.ti = 50

# Number of time-dependent covariates
p.tv = 50

# Divide the data in 100 chunks for computation
K = 100; n = n.subject/K

# Correlation between covariates
cor = 0.5

# Grid of tuning parameter for adaptive LASSO
lambda.grid = 10^seq(-10,3,0.01)


# Generate Data
dat.mat0 = SIM.FUN.TVC(p.ti, p.tv, n.subject, cor, beta0.ti = NULL, beta0.tv = NULL) # First 50 covariates are time-independent, followed by 50 time-dependent covariates: beta0.ti = c(c(rep(0.08,3),rep(0.04,3),rep(0.02,3)),rep(0,p.ti-9)); beta0.tv = c(c(rep(0.08,3),rep(0.04,3),rep(0.02,3)),rep(0,p.tv-9))

# Divide data
dat.list = lapply(1:K,function(kk){dat.mat0[dat.mat0[,dim(dat.mat0)[2]] %in% ((1:n)+(kk-1)*n),-dim(dat.mat0)[2]]})
```

## Analysis of simulated time-dependent survival data
```
### Round 1
ptm = proc.time()
bini = coxph(Surv(dat.list[[1]][,1],dat.list[[1]][,2],dat.list[[1]][,3])~dat.list[[1]][,-c(1,2,3)])$coef
update1 = iteration.fun.cox3(dat.list=dat.list,bini=bini,kk.list=1:K);
bnew1 = apply(update1$b.k,1,mean);
time1.0 = proc.time() - ptm;gc();print(time1.0[3])


### Round 2
ptm = proc.time()
update2 = iteration.fun.cox3(dat.list=dat.list,bini=bnew1,kk.list=1:K);
bnew2 = apply(update2$b.k,1,mean)
time2.0 = proc.time() - ptm;gc();print(time2.0[3])


### Adaptive LASSO estimation
ptm = proc.time()
info2  = -update2$Ahat
Ahalf2 = svd(-update2$Ahat); Ahalf2 = Ahalf2$u%*%diag(sqrt(Ahalf2$d))%*%t(Ahalf2$v);
bhat2.DCOS = Est.ALASSO.Approx.GLMNET(Ahalf2%*%bnew2,Ahalf2,bnew2,N.adj=sum(dat.mat0[,2]),lambda.grid=lambda.grid,N.inflate = dim(dat.mat0)[1])
time2.1 = proc.time() - ptm;gc();print(time2.1[3])
```

## View results
```
# Point estimate
bhat2.DCOS$bhat.BIC # Using BIC to tune adaptive LASSO parameter
bhat2.DCOS$lambda.BIC # Optimal lambda

# Variance-covariance matrix for active set
solve(info2[bhat2.DCOS$bhat.BIC!=0,bhat2.DCOS$bhat.BIC!=0]*n)

# Lower bound of 95% CI for active set
(bhat2.DCOS$bhat.BIC - qnorm(0.975) * sqrt(diag(solve(info2*n))))[bhat2.DCOS$bhat.BIC!=0]
# Upper bound of 95% CI for active set
(bhat2.DCOS$bhat.BIC + qnorm(0.975) * sqrt(diag(solve(info2*n))))[bhat2.DCOS$bhat.BIC!=0]
```
