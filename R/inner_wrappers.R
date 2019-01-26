#' Least square approximated adaptive lasso
#'
#' \code{Est.ALASSO.Approx.GLMNET} fits adaptive lasso based on least square approximation (LSA)
#'
#' @param ynew a p by 1 vector of LSA-based outcome
#' @param xnew a p by p matrix of LSA-based covariates
#' @param bini 1/bini is used as the penalty for adaptive lasso
#' @param N.adj log(N.adj) is used to penalize BIC, where BIC = -2loglik + df*log(N.adj), and modified BIC, where modBIC = -2loglik + df*N.adj^0.1
#' @param lambda.grid the grid of lambda to put into glmnet
#' @param modBIC if the program also finds the lamda that minimizes the modified BIC. The default is TRUE.
#' @return bhat.BIC adaptive lasso estimator using BIC and bhat.modBIC adaptive lasso estimator using modified BIC. lambda.BIC and lambda.modBIC return the optimal lambda chosen by BIC and modified BIC.
#' @param N.inflate the factor to inflate sum residual squared
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples Est.ALASSO.Approx.GLMNET(ynew,xnew,bini,N.adj)
Est.ALASSO.Approx.GLMNET = function(ynew,xnew,bini,N.adj,lambda.grid, modBIC = T,N.inflate=N.adj){
  w.b = 1/abs(bini)
  tmpfit = glmnet(x=xnew, y=ynew, family='gaussian',penalty.factor = w.b,
                  alpha=1, lambda = lambda.grid, intercept=F)
  LL = apply((c(ynew) - predict(tmpfit,xnew,type = 'response'))^2, 2, sum)*N.inflate
  ### Modified BIC
  if (modBIC){
    BIC.lam = LL+min(N.adj^0.1,log(N.adj))*tmpfit$df
    m.opt = which.min(BIC.lam); bhat.modBIC = tmpfit$beta[,m.opt]; lamhat.modBIC = tmpfit$lambda[m.opt]
  }else{
    bhat.modBIC = lamhat.modBIC = NA
  }
  ### BIC
  BIC.lam = LL+log(N.adj)*tmpfit$df
  m.opt = which.min(BIC.lam); bhat.BIC = tmpfit$beta[,m.opt]; lamhat.BIC = tmpfit$lambda[m.opt]
  return (list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC,
               lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC))
}








#' A wrapper for adaptive lasso
#'
#' \code{Est.ALASSO.GLMNET} fits adaptive lasso based on glmnet. The best lambda (penalizing factor) is chosen by BIC or modified BIC.
#'
#' @param data matrix or data frame. For non-survival outcome, first column = response, rest = design matrix.
#' For continuous-time survival outcome, first column = time, second column = delta (0/1), rest = design matrix.
#' @param BIC.factor factor in modified BIC, BIC = -2 loglikelihood + df * N^BIC.factor
#' @param fam0 family of the response, taking "binomial", "Poisson", "Cox"
#' @param w.b w.b used to penalize adaptive lasso. If null, a glm/Cox model will be fitted and 1/abs(coefficients) will be used as w.b
#' @param lambda.grid the grid of lambda to put into glmnet
#' @param modBIC if the program also finds the lamda that minimizes the modified BIC. The default is TRUE.
#' @return bhat.BIC adaptive lasso estimator using BIC and bhat.modBIC adaptive lasso estimator using modified BIC. lambda.BIC and lambda.modBIC return the optimal lambda chosen by BIC and modified BIC.
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples Est.ALASSO.Approx.GLMNET(ynew,xnew,bini,N.adj)
Est.ALASSO.GLMNET = function(data,BIC.factor=0.1,fam0="binomial", w.b = NULL, lambda.grid, modBIC = T){
  if (fam0 != 'Cox'){
    data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x);gc()
    if (is.null(w.b)){
      bini = glm(y~x,family=fam0)$coef; w.b = 1/abs(bini[-1]); gc()
    }
    tmpfit = glmnet(x=x, y=y, family=fam0,penalty.factor = w.b,
                    alpha=1, lambda = lambda.grid, intercept=T)
    gc()

    N.adj = dim(x)[1]
  }else{
    data = as.matrix(data); y = cbind(time=data[,1],status=data[,2]); x = data[,-c(1,2),drop=F]; nn=length(y); pp = ncol(x);gc()
    if (is.null(w.b)){
      bini = coxph(Surv(time = y[,1],event = y[,2])~x)$coef; w.b = 1/abs(bini); gc()
    }
    tmpfit = glmnet(x,y, family='cox',penalty.factor = w.b,
                    alpha=1, lambda = lambda.grid)
    gc()

    N.adj = sum(y[,2])
  }
  dev = deviance(tmpfit)
  ### Modified BIC
  if (modBIC){
    BIC.lam = dev+min(N.adj^0.1,log(N.adj))*tmpfit$df
    m.opt = which.min(BIC.lam); bhat.modBIC = c(tmpfit$a0[m.opt],tmpfit$beta[,m.opt]); lamhat.modBIC = tmpfit$lambda[m.opt]
  }else{
    bhat.modBIC = lamhat.modBIC = NA
  }
  ### BIC
  BIC.lam = dev+log(N.adj)*tmpfit$df
  m.opt = which.min(BIC.lam); bhat.BIC = c(tmpfit$a0[m.opt],tmpfit$beta[,m.opt]); lamhat.BIC = tmpfit$lambda[m.opt]
  return (list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC,
               lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC))
}








#' A modified wrapper for adaptive lasso with the optimal penalty chosen by cross-validation
#'
#' \code{Est.ALASSO.GLMNET.CV} fits adaptive lasso based on cv.glmnet. The best lambda (penalizing factor) is chosen by 10-fold cross-validation.
#'
#' @param data matrix or data frame. For non-survival outcome, first column = response, rest = design matrix.
#' For continuous-time survival outcome, first column = time, second column = delta (0/1), rest = design matrix.
#' @param fam0 family of the response, taking "binomial", "Poisson", "Cox"
#' @param w.b w.b used to penalize adaptive lasso. If null, a glm/Cox model will be fitted and 1/abs(coefficients) will be used as w.b
#' @param lambda.grid the grid of lambda to put into glmnet
#' @param chunksize The prediction step in cv.glmnet could take a large amount of memory. chunksize specifies how many chunks you would like to split the prediction step. The predition step will be run in a loop if chunksize>1.
#' @return a list containing two arguments: bhat.cv adaptive lasso estimator using 10-fold cross-validation; lambda.cv is the optimal lambda chosen by cross-validation.
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples Est.ALASSO.GLMNET.CV(data, fam0="binomial", w.b = NULL, lambda.grid, chunksize = 50)
Est.ALASSO.GLMNET.CV = function(data, fam0="binomial", w.b = NULL, lambda.grid, chunksize = 50){

  if (fam0 != 'Cox'){
    data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x);gc()
    if (is.null(w.b)){
      bini = glm(y~x,family=fam0)$coef; w.b = 1/abs(bini[-1]); gc()
    }
    tmpfit = mycv.glmnet(x=x, y=y, family=fam0,lambda = lambda.grid,
                         penalty.factor = w.b,nfolds=10,alpha=1, intercept=T)
    gc()
    chunks=split(tmpfit$lambda,factor(1:ceiling(length(tmpfit$lambda)/chunksize)))
    lambdanew = cvmnew = cvsdnew = c()
    for (cc in 1:length(chunks)){
      tmpcv = mycv.lognet(tmpfit$outlist,chunks[[cc]],tmpfit$x,tmpfit$y,tmpfit$weights,tmpfit$offset,
                          tmpfit$foldid,tmpfit$type.measure,tmpfit$grouped,tmpfit$keep)
      lambdanew = c(lambdanew, chunks[[cc]])
      cvmnew = c(cvmnew, tmpcv$cvm)
      cvsdnew = c(cvsdnew, tmpcv$cvsd)
      print(cc)
    }
    nas = is.na(cvsdnew)
    if (any(nas)) {
      lambdanew = lambdanew[!nas]
      cvmnew = cvmnew[!nas]
      cvsdnew = cvsdnew[!nas]
    }
    od = order(lambdanew,decreasing=T)
    lambdanew = lambdanew[od]
    cvmnew = cvmnew[od]
    cvsdnew = cvsdnew[od]
    head(lambdanew,20)

    cvname = tmpcv$name
    out = list(lambda = lambdanew, cvm = cvmnew, cvsd = cvsdnew, cvup = cvmnew +
                 cvsdnew, cvlo = cvmnew - cvsdnew, nzero = tmpfit$nz, name = cvname,
               glmnet.fit = tmpfit$glmnet.object)
    if (tmpfit$keep)
      out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin = if (cvname == "AUC") {
      getmin(lambdanew, -cvmnew, cvsdnew)
    }  else {
      getmin(lambdanew, cvmnew, cvsdnew)
    }
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    bhat = coef(obj, s='lambda.min')
    return(list(bhat.cv = bhat, lambda.cv = lamin))
  }else{
    ##### Cox Part
    data = as.matrix(data); y = cbind(time=data[,1],status=data[,2]); x = data[,-c(1,2),drop=F]; nn=length(y); pp = ncol(x);gc()
    if (is.null(w.b)){
      bini = coxph(Surv(time = y[,1],event = y[,2])~x)$coef; w.b = 1/abs(bini); gc()
    }
    tmpfit = mycv.glmnet(x,Surv(time = y[,1],event = y[,2]), family='cox',penalty.factor = w.b,nfolds=10,
                         alpha=1, lambda = lambda.grid)
    chunks=split(tmpfit$lambda,factor(1:ceiling(length(tmpfit$lambda)/chunksize)))
    lambdanew = cvmnew = cvsdnew = c()
    for (cc in 1:length(chunks)){
      tmpcv = mycv.coxnet(tmpfit$outlist,chunks[[cc]],tmpfit$x,tmpfit$y,tmpfit$weights,tmpfit$offset,
                          tmpfit$foldid,tmpfit$type.measure,tmpfit$grouped,tmpfit$keep)
      lambdanew = c(lambdanew, chunks[[cc]])
      cvmnew = c(cvmnew, tmpcv$cvm)
      cvsdnew = c(cvsdnew, tmpcv$cvsd)
      print(cc)
    }
    nas = is.na(cvsdnew)
    if (any(nas)) {
      lambdanew = lambdanew[!nas]
      cvmnew = cvmnew[!nas]
      cvsdnew = cvsdnew[!nas]
    }
    od = order(lambdanew,decreasing=T)
    lambdanew = lambdanew[od]
    cvmnew = cvmnew[od]
    cvsdnew = cvsdnew[od]
    head(lambdanew,20)

    cvname = tmpcv$name
    out = list(lambda = lambdanew, cvm = cvmnew, cvsd = cvsdnew, cvup = cvmnew +
                 cvsdnew, cvlo = cvmnew - cvsdnew, nzero = tmpfit$nz, name = cvname,
               glmnet.fit = tmpfit$glmnet.object)
    if (tmpfit$keep)
      out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin = if (cvname == "AUC") {
      getmin(lambdanew, -cvmnew, cvsdnew)
    }  else {
      getmin(lambdanew, cvmnew, cvsdnew)
    }
    obj = c(out, as.list(lamin))
    class(obj) = "cv.glmnet"
    bhat = coef(obj, s='lambda.min')
    return(list(bhat.cv = bhat, lambda.cv = lamin))
  }
}
















#' A wrapper for divide-and-conquer adaptive lasso proposed by Chen and Xie (2014) and Tang et al. (2016).
#'
#' \code{Est.ALASSO.GLMNET.TANGXIE} fits adaptive lasso based on cv.glmnet. The best lambda (penalizing factor) is chosen by 10-fold cross-validation.
#'
#' @param dat.list a list of matrices. Each element of the list is a sub-dataset. In each sub-dataset, for non-survival outcome, first column = response, rest = design matrix.
#' In each sub-dataset, for continuous-time survival outcome, first column = time, second column = delta (0/1), rest = design matrix.
#' @param K Number of sub-datasets in dat.list
#' @param w.b w.b used to penalize adaptive lasso. If null, a glm/Cox model will be fitted and 1/abs(coefficients) will be used as w.b
#' @param BIC.factor factor in modified BIC, BIC = -2 loglikelihood + df * N^BIC.factor
#' @param fam0 family of the response, taking "binomial", "Poisson", "Cox"
#' @param lambda.grid the grid of lambda to put into glmnet
#' @param mvpct majority voting percentage used for Chen and Xie (2014)
#' @param modBIC if the program also finds the lamda that minimizes the modified BIC. The default is TRUE.
#' @param adaptive if adaptive lasso is used. The default is TRUE
#' @param onestep if one-step estimator should be used as the initial estimator for Cox fit
#' @return a list containing two arguments: bhat.cv adaptive lasso estimator using 10-fold cross-validation; lambda.cv is the optimal lambda chosen by cross-validation.
#' @references Chen, Xueying, and Min-ge Xie. "A split-and-conquer approach for analysis of extraordinarily large data." Statistica Sinica (2014): 1655-1684.
#' @references Tang, Lu, Ling Zhou, and Peter X-K. Song. "Method of Divide-and-Combine in Regularised Generalised Linear Models for Big Data." arXiv preprint arXiv:1611.06208 (2016).
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples Est.ALASSO.GLMNET.TANGXIE(dat.list,K,BIC.factor=0.1,fam0="binomial",lambda.grid,mvpct = 0.5)
Est.ALASSO.GLMNET.TANGXIE = function(dat.list,K,BIC.factor=0.1,fam0="binomial",lambda.grid,mvpct = 0.5, modBIC = T,
                                     adaptive = T, onestep = F){
  # Adaptive?
  if (adaptive){
    # Logistic
    if (fam0 != 'Cox'){
      bvec.DC = sapply(1:K,function(kk){
        glm(dat.list[[kk]][,1]~dat.list[[kk]][,-1],family=binomial)$coef
      });
      bini.DC = apply(bvec.DC,1,mean);  w.b = 1/abs(bini.DC[-1])
    }else{
      # Cox without one-step estimator
      if (!onestep){
        bvec.DC = sapply(1:K,function(kk){
          coxph(Surv(time = dat.list[[kk]][,1],event = dat.list[[kk]][,2])~dat.list[[kk]][,-c(1,2)])$coef
        });
        bini.DC = apply(bvec.DC,1,mean);  w.b = 1/abs(bini.DC)
      }else{
        # Cox with one-step estimator
        ### Round 1
        bini = coxph(Surv(time = dat.list[[1]][,1],event = dat.list[[1]][,2])~dat.list[[1]][,-c(1,2)])$coef
        update1 = iteration.fun.cox2(dat.list=dat.list,bini=bini,kk.list=2:K);
        bnew1 = apply(cbind(bini,update1$b.k),1,mean);

        ### Round 2
        update2 = iteration.fun.cox2(dat.list=dat.list,bini=bnew1,kk.list=2:K);
        bnew2 = apply(cbind(bini,update2$b.k),1,mean)

        ### Round 3
        update3 = iteration.fun.cox2(dat.list=dat.list,bini=bnew2,kk.list=2:K);
        bnew3 = apply(cbind(bini,update3$b.k),1,mean)

        ### Round 4
        update4 = iteration.fun.cox2(dat.list=dat.list,bini=bnew3,kk.list=2:K);
        bnew4 = apply(cbind(bini,update4$b.k),1,mean)

        ### Round 5
        update5 = iteration.fun.cox2(dat.list=dat.list,bini=bnew4,kk.list=2:K);
        bnew5 = apply(cbind(bini,update5$b.k),1,mean)

        w.b = 1/abs(bnew5)
      }
    }
  }else{
    nvars = dim(dat.list[[1]][,-c(1,2)])[2]
    w.b = rep(1, nvars)
  }


  if (fam0 != 'Cox'){
    bhat.BIC.all = bhat.modBIC.all = bhat_c.modBIC.all = bhat_c.BIC.all = matrix(NA, nrow = K, ncol = length(w.b)+1)
  }else{
    bhat.BIC.all = bhat.modBIC.all = bhat_c.modBIC.all = bhat_c.BIC.all = matrix(NA, nrow = K, ncol = length(w.b))
  }
  lambda.BIC = lambda.modBIC = rep(NA,K)


  for (kk in 1:K){
    if (fam0 != 'Cox'){
      data = dat.list[[kk]]; y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x);gc()
      tmpfit = glmnet(x=x, y=y, family=fam0,penalty.factor = w.b,
                      alpha=1, lambda = lambda.grid, intercept=T)
      gc()

      N.adj = dim(x)[1]
    }else{
      data = dat.list[[kk]]; y = cbind(time=data[,1],status=data[,2]); x = data[,-c(1,2),drop=F]; nn=length(y); pp = ncol(x);gc()
      tmpfit = glmnet(x,y, family='cox',penalty.factor = w.b,
                      alpha=1, lambda = lambda.grid)
      gc()

      N.adj = sum(y[,2])
    }


    dev = deviance(tmpfit)
    ### Modified BIC
    if (modBIC){
      BIC.lam = dev+min(N.adj^0.1,log(N.adj))*tmpfit$df
      m.opt = which.min(BIC.lam); bhat.modBIC = tmpfit$beta[,m.opt]; lambda.modBIC[kk] = lamhat = tmpfit$lambda[m.opt]
      if (fam0 != 'Cox'){
        a0.modBIC = tmpfit$a0[m.opt]
        xi = cbind(rep(1,length(y)),x)
        pred = as.vector(g.logit(xi %*% c(a0.modBIC,bhat.modBIC)))
        XPX.modBIC = t(xi) %*% ((pred*(1-pred)) * xi)
        bhat_c.modBIC = c(a0.modBIC,bhat.modBIC) + solve(XPX.modBIC) %*% (t(xi) %*% (y-pred))
        if (kk == 1){
          J.modBIC = XPX.modBIC
          I.modBIC = XPX.modBIC %*% bhat_c.modBIC
        }else{
          J.modBIC = J.modBIC + XPX.modBIC
          I.modBIC = I.modBIC + XPX.modBIC %*% bhat_c.modBIC
        }
      }
    }


    ### BIC
    BIC.lam = dev+log(N.adj)*tmpfit$df
    m.opt = which.min(BIC.lam); bhat.BIC = tmpfit$beta[,m.opt]; lambda.BIC[kk] = lamhat = tmpfit$lambda[m.opt]

    if (fam0!='Cox'){
      a0.BIC = tmpfit$a0[m.opt]
      xi = cbind(rep(1,length(y)),x)
      pred = as.vector(g.logit(xi %*% c(a0.BIC,bhat.BIC)))
      XPX.BIC = t(xi) %*% ((pred*(1-pred)) * xi)
      bhat_c.BIC = c(a0.BIC,bhat.BIC) + solve(XPX.BIC) %*% (t(xi) %*% (y-pred))

      if (kk == 1){
        J.BIC = XPX.BIC
        I.BIC = XPX.BIC %*% bhat_c.BIC
      }else{
        J.BIC = J.BIC + XPX.BIC
        I.BIC = I.BIC + XPX.BIC %*% bhat_c.BIC
      }
    }



    if (fam0 != 'Cox'){
      bhat.BIC.all[kk,] = c(a0.BIC,bhat.BIC)
      if (modBIC){
        bhat.modBIC.all[kk,] = c(a0.modBIC,bhat.modBIC)
      }
      bhat_c.BIC.all[kk,] = bhat_c.BIC
      if (modBIC){
        bhat_c.modBIC.all[kk,] = bhat_c.modBIC
      }
    }else{
      bhat.BIC.all[kk,] = bhat.BIC
      if (modBIC){
        bhat.modBIC.all[kk,] = bhat.modBIC
      }
    }

    print(kk);gc()
  }

  if (fam0 !='Cox'){
    beta_tang.BIC = as.vector(solve(J.BIC) %*% I.BIC)
    if (modBIC){
      beta_tang.modBIC = as.vector(solve(J.modBIC) %*% I.modBIC)
    }else{
      beta_tang.modBIC = rep(NA, length(w.b))
    }

    ind = colMeans(bhat.BIC.all!=0)>mvpct
    beta_xie.BIC = rep(0, length(w.b)+1)
    beta_xie.BIC[ind] = as.vector(solve(J.BIC[ind,ind]) %*% I.BIC[ind,1])
    if (modBIC){
      ind = colMeans(bhat.modBIC.all!=0)>mvpct
      beta_xie.modBIC = rep(0, length(w.b)+1)
      beta_xie.modBIC[ind] = as.vector(solve(J.modBIC[ind,ind]) %*% I.modBIC[ind,1])
    }else{
      beta_xie.modBIC = rep(NA, length(w.b))
    }
  }else{
    ind = colMeans(bhat.BIC.all!=0)>mvpct
    beta_xie.BIC = rep(0, length(w.b))
    beta_xie.BIC[ind] = colMeans(bhat.BIC.all[,ind])
    if (modBIC){
      ind = colMeans(bhat.modBIC.all!=0)>mvpct
      beta_xie.modBIC = rep(0, length(w.b))
      beta_xie.modBIC[ind] = colMeans(bhat.BIC.all[,ind])
    }else{
      beta_xie.modBIC = rep(NA, length(w.b))
    }
  }


  if (fam0 !='Cox'){
    return(list(beta_xie.modBIC=beta_xie.modBIC,beta_xie.BIC=beta_xie.BIC,beta_tang.modBIC=beta_tang.modBIC,beta_tang.BIC=beta_tang.BIC,lambda.BIC = lambda.BIC, lambda.modBIC =lambda.modBIC))
  }else{
    return(list(beta_xie.modBIC=beta_xie.modBIC,beta_xie.BIC=beta_xie.BIC,lambda.BIC = lambda.BIC, lambda.modBIC =lambda.modBIC))
  }
}

#' A wrapper for divide-and-conquer MCP Cox model proposed by Chen and Xie (2014) .
#'
#' \code{MCP.Xie} fits Cox model with an MCP penalty. The best lambda (penalizing factor) is chosen by BIC.
#'
#' @param dat.list a list of matrices. Each element of the list is a sub-dataset. In each sub-dataset, for non-survival outcome, first column = response, rest = design matrix.
#' In each sub-dataset, for continuous-time survival outcome, first column = time, second column = delta (0/1), rest = design matrix.
#' @param K Number of sub-datasets in dat.list
#' @param BIC.factor factor in modified BIC, BIC = -2 loglikelihood + df * N^BIC.factor
#' @param lambda.grid the grid of lambda to put into glmnet
#' @param mvpct majority voting percentage used for Chen and Xie (2014)
#' @return a list containing the vector of estimates.
#' @references Chen, Xueying, and Min-ge Xie. "A split-and-conquer approach for analysis of extraordinarily large data." Statistica Sinica (2014): 1655-1684.
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples MCP.Xie(dat.list,K,BIC.factor=0.1,lambda.grid,mvpct = 0.5)
MCP.Xie = function(dat.list,K,lambda.grid,mvpct = 0.5){
  bhat.BIC.all = matrix(NA, nrow = K, ncol = dim(dat.list[[1]])[2]-2)
  lambda.BIC = rep(NA,K)

  for (kk in 1:K){
    tmpfit = ncvsurv2(X = dat.list[[kk]][,-c(1,2),drop=F], y = cbind(dat.list[[kk]][,1],dat.list[[kk]][,2]),
                     penalty = 'MCP', lambda = lambda.grid)
    gc()

    N.adj = sum(dat.list[[kk]][,2])
    dev = 2*tmpfit$loss

    ### BIC
    BIC.lam = dev+log(N.adj)*colSums(tmpfit$beta!=0)
    m.opt = which.min(BIC.lam); bhat.BIC = tmpfit$beta[,m.opt]; lambda.BIC[kk] = lamhat = tmpfit$lambda[m.opt]

    bhat.BIC.all[kk,] = bhat.BIC
    print(kk);gc()
  }

  ind = colMeans(bhat.BIC.all!=0)>mvpct
  beta_xiemcp.BIC = rep(0, dim(dat.list[[1]])[2]-2)
  beta_xiemcp.BIC[ind] = colMeans(bhat.BIC.all[,ind])

  return(list(beta_xie.BIC=beta_xiemcp.BIC))

}

