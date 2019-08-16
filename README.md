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
  # V1        V2        V3        V4        V5        V6        V7        V8        V9       V10       V11       V12       V13       V14       V15       V16       V17       V18
  # 0.8014559 0.8007799 0.7995421 0.3993714 0.4030890 0.3972277 0.2036702 0.1979250 0.1967845 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
  # V19       V20       V21       V22       V23       V24       V25       V26       V27       V28       V29       V30       V31       V32       V33       V34       V35       V36
  # 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
  # V37       V38       V39       V40       V41       V42       V43       V44       V45       V46       V47       V48       V49       V50
  # 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
bhat2.DCOS$lambda.BIC # Optimal lambda
  # [1] 0.004073803

# Variance-covariance matrix for active set
solve(info2[bhat2.DCOS$bhat.BIC!=0,bhat2.DCOS$bhat.BIC!=0]*n)
  # [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]          [,9]
  # [1,]  6.913335e-04  2.338834e-05  2.594814e-05 -2.330726e-05 -2.004602e-05 -2.184883e-05 -4.379367e-05 -4.372946e-05 -4.434652e-05
  # [2,]  2.338834e-05  6.922461e-04  2.316011e-05 -2.021649e-05 -2.129422e-05 -2.340024e-05 -4.253062e-05 -4.546958e-05 -4.237310e-05
  # [3,]  2.594814e-05  2.316011e-05  6.914615e-04 -2.071537e-05 -2.116795e-05 -2.167582e-05 -4.448787e-05 -4.536370e-05 -4.571958e-05
  # [4,] -2.330726e-05 -2.021649e-05 -2.071537e-05  6.215142e-04 -4.352129e-05 -4.619849e-05 -5.406780e-05 -5.604920e-05 -5.498740e-05
  # [5,] -2.004602e-05 -2.129422e-05 -2.116795e-05 -4.352129e-05  6.239145e-04 -4.319653e-05 -5.421495e-05 -5.513638e-05 -5.554325e-05
  # [6,] -2.184883e-05 -2.340024e-05 -2.167582e-05 -4.619849e-05 -4.319653e-05  6.205542e-04 -5.484947e-05 -5.344318e-05 -5.569948e-05
  # [7,] -4.379367e-05 -4.253062e-05 -4.448787e-05 -5.406780e-05 -5.421495e-05 -5.484947e-05  6.089977e-04 -6.442752e-05 -6.279552e-05
  # [8,] -4.372946e-05 -4.546958e-05 -4.536370e-05 -5.604920e-05 -5.513638e-05 -5.344318e-05 -6.442752e-05  6.054489e-04 -5.912090e-05
  # [9,] -4.434652e-05 -4.237310e-05 -4.571958e-05 -5.498740e-05 -5.554325e-05 -5.569948e-05 -6.279552e-05 -5.912090e-05  6.064989e-04

# Lower bound of 95% CI for active set
(bhat2.DCOS$bhat.BIC - qnorm(0.975) * sqrt(diag(solve(info2*n))))[bhat2.DCOS$bhat.BIC!=0]
  # V1        V2        V3        V4        V5        V6        V7        V8        V9
  # 0.7480016 0.7472489 0.7460384 0.3484353 0.3521123 0.3463074 0.1531643 0.1476124 0.1464508
# Upper bound of 95% CI for active set
(bhat2.DCOS$bhat.BIC + qnorm(0.975) * sqrt(diag(solve(info2*n))))[bhat2.DCOS$bhat.BIC!=0]
  # V1        V2        V3        V4        V5        V6        V7        V8        V9
  # 0.8549101 0.8543109 0.8530459 0.4503074 0.4540658 0.4481480 0.2541761 0.2482375 0.2471181

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
        V1         V2         V3         V4         V5         V6         V7         V8         V9        V10        V11        V12        V13        V14        V15        V16        V17 
0.07937987 0.08296920 0.08114094 0.03930308 0.03633204 0.04046889 0.01942764 0.01789658 0.02157509 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
       V18        V19        V20        V21        V22        V23        V24        V25        V26        V27        V28        V29        V30        V31        V32        V33        V34 
0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
       V35        V36        V37        V38        V39        V40        V41        V42        V43        V44        V45        V46        V47        V48        V49        V50        V51 
0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.07797483 
       V52        V53        V54        V55        V56        V57        V58        V59        V60        V61        V62        V63        V64        V65        V66        V67        V68 
0.08080146 0.08066730 0.04172541 0.04309695 0.04258372 0.01837601 0.01799859 0.01856933 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
       V69        V70        V71        V72        V73        V74        V75        V76        V77        V78        V79        V80        V81        V82        V83        V84        V85 
0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 
       V86        V87        V88        V89        V90        V91        V92        V93        V94        V95        V96        V97        V98        V99       V100 
0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000

```
