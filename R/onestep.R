######## One-step estimator functions for Cox ########
#' Vector to matrix
#'
#' \code{VTM} Replicate vector vc by dm times and create a dm by length(vc) matrix
#' @export

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}
#' Cumulative sum
#'
#' \code{cumsum2} calculates cumulative sum
#' @export
cumsum2 <- function(mydat)     #cumsum by row, col remains the same
{
  if(is.null(dim(mydat))) return(cumsum(mydat))
  else{
    out <- matrix(cumsum(mydat), nrow=nrow(mydat))
    out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
    return(out)
  }
}
#' Inner function for one-step estimator
#'
#' \code{MySum} is an inner function to calculate score and negative information.
#' @export
MySum <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
{
  if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
  if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
  pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)

  if(is.null(Vi)){return(pos)}else{
    Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
    out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
    out[pos!=0,] <- Vi[pos,]
    if(is.null(dim(Vi))) out <- c(out)
    return(out) ## n.y x p
  }
}
#' Inner function for one-step estimator
#'
#' \code{PI.k.FUN} is an inner function to calculate score and negative information.
#' @export
PI.k.FUN <- function(tt,ebzi,xi,zi,k0=0,vi=NULL)
{
  out = ebzi; pz = ncol(zi); nn=length(out)
  if(k0==1){out=out*zi}
  if(k0==2){out = out*zi[,rep(1:pz,pz)]*zi[,rep(1:pz,rep(pz,pz))]}
  as.matrix(MySum(tt,"<=",xi,Vi=out)/nn)
}


#' Calculate score and negative information (A) matrix for a Cox model
#'
#' \code{Score.A.FUN} evaluates score and negative information (A)
#'
#' @param data a matrix with first column = U, second column = delta (0/1), rest = design
#' @param betahat at which coefficient level to evaluate score and negative information
#' @param rtn 'Score' returns only the score function, 'Score+A' returns score and negative information (second derivative of PL), 'Score+A+Approx' returns score and negative information (approximated by -SS')
#' @export
Score.A.FUN <- function(data,betahat,rtn="score")
{
  xi = data[,1]; di = data[,2]; zi = data[,-(1:2),drop=F]; ebzi = c(exp(zi%*%betahat)); pz = ncol(zi); nn = length(xi);
  tmpind = di==1; tj = xi[tmpind]; zj = zi[tmpind,,drop=F]
  pi0.tj   = c(PI.k.FUN(tj,ebzi,xi,zi,k0=0)); pi1.tj   = PI.k.FUN(tj,ebzi,xi,zi,k0=1)
  smat = zj - pi1.tj/pi0.tj
  score = apply(smat,2,sum)/dim(data)[1]
  if(rtn=="score"){
    return(score)
  }else if (rtn=='Score+A'){
    Ahat = PI.k.FUN(tj,ebzi,xi,zi,k0=2)*pi0.tj - pi1.tj[,rep(1:pz,pz)]*pi1.tj[,rep(1:pz,rep(pz,pz))]
    Ahat = -matrix(apply(Ahat/pi0.tj^2,2,sum),ncol=pz)/dim(data)[1]
    return(list("score"=score,"neg_info_mat"=Ahat))
  }else if (rtn=='Score+A+Approx'){
    Ahat = -t(smat) %*% smat/dim(data)[1]
    return(list("score"=score,"neg_info_mat"=Ahat))
  }
}


#' One-step iteration function for divide-and-conquer Cox model
#'
#' \code{iteration.fun.cox} estimates a one-step for a Cox regression for each subset.
#'
#' @param dat.list list of subsets (after dividing). In each subset, first column = U, second column = delta (0/1), rest = design matrix
#' @param bini initial estimator as starting point
#' @param kk.list which subsets of dat.list get one-step update
#' @param rtn return score, score+negative information, score+negative information using -SS' approximation
#' @return a list with b.k a matrix of one-step estimators and with Ahat the negative information matrix
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples iteration.fun.cox(dat.list=dat.list,bini=bini,kk.list=2:K,rtn='Score+A+Approx')
iteration.fun.cox = function(dat.list,bini,kk.list,rtn='Score+A+Approx'){
  K = length(dat.list)
  ScoreA.list = lapply(1:K,function(kk){Score.A.FUN(dat.list[[kk]],bini,rtn)})
  Uini.list = sapply(kk.list,function(kk){ScoreA.list[[kk]]$score})
  Aini.list = lapply(1:K,function(kk){ScoreA.list[[kk]]$neg_info_mat});Ahat.ini = Reduce("+",Aini.list)/K
  bhat.list = -solve(Ahat.ini)%*%Uini.list + bini;
  list("b.k"=bhat.list,"Ahat"=Ahat.ini)
}

#' One-step iteration function for divide-and-conquer Cox model
#'
#' \code{iteration.fun.cox2} estimates a one-step for a Cox regression for each subset. The core used a C function.
#'
#' @param dat.list list of subsets (after dividing). In each subset, first column = U, second column = delta (0/1), rest = design matrix
#' @param bini initial estimator as starting point
#' @param kk.list which subsets of dat.list get one-step update
#' @return a list with b.k a matrix of one-step estimators and with Ahat the negative information matrix
#' @author Yan Wang, Tianxi Cai
#' @importFrom Rcpp evalCpp
#' @useDynLib divideconquer ScoreAFUNC
#' @export
#' @examples iteration.fun.cox(dat.list=dat.list,bini=bini,kk.list=2:K,rtn='Score+A+Approx')
iteration.fun.cox2 = function(dat.list,bini,kk.list){
  K = length(dat.list)
  ScoreA.list = lapply(1:K,function(kk){
    sorted = order(dat.list[[kk]][,1])
    tmp = .Call("ScoreAFUNC",as.double(dat.list[[kk]][sorted,1]),
                  as.integer(dat.list[[kk]][sorted,2]),
                  dat.list[[kk]][sorted,-c(1,2)],
                  as.double(bini))
    list(score=tmp[[1]],neg_info_mat=matrix(tmp[[2]],length(bini),length(bini)))
  })
  Uini.list = sapply(kk.list,function(kk){ScoreA.list[[kk]]$score})
  Aini.list = lapply(1:K,function(kk){ScoreA.list[[kk]]$neg_info_mat});Ahat.ini = Reduce("+",Aini.list)/K
  bhat.list = -solve(Ahat.ini)%*%Uini.list + bini;
  list("b.k"=bhat.list,"Ahat"=Ahat.ini)
}



#' One-step iteration function for divide-and-conquer time-varying Cox model
#'
#' \code{iteration.fun.cox3} estimates a one-step for a time-varying Cox regression for each subset. The core used a C function.
#'
#' @param dat.list list of subsets (after dividing). In each subset, first column = start time, second column = stop time, third time = event (0/1), rest = design matrix
#' @param bini initial estimator as starting point
#' @param kk.list which subsets of dat.list get one-step update
#' @return a list with b.k a matrix of one-step estimators and with Ahat the negative information matrix
#' @author Yan Wang, Tianxi Cai
#' @importFrom Rcpp evalCpp
#' @useDynLib divideconquer ScoreAFUNC
#' @export
#' @examples iteration.fun.cox(dat.list=dat.list,bini=bini,kk.list=2:K,rtn='Score+A+Approx')
iteration.fun.cox3 = function(dat.list,bini,kk.list){
  K = length(dat.list)
  ScoreA.list = lapply(1:K,function(kk){
    sort.end <- order(-dat.list[[kk]][,2]) - 1L
    sort.start <- order(-dat.list[[kk]][,1]) - 1L
    y = Surv(dat.list[[kk]][,1],dat.list[[kk]][,2],dat.list[[kk]][,3])
    x = dat.list[[kk]][,-c(1:3)]
    storage.mode(y) <- storage.mode(x) <- "double"
    tmp = .Call("ScoreAFUNC_AG", y, x,
                as.integer(dim(dat.list[[kk]])[1]),
                sort.start,sort.end,as.double(bini))
    list(score=tmp[[1]],neg_info_mat=matrix(tmp[[2]],length(bini),length(bini)))
  })
  Uini.list = sapply(kk.list,function(kk){ScoreA.list[[kk]]$score})
  Aini.list = lapply(1:K,function(kk){ScoreA.list[[kk]]$neg_info_mat});Ahat.ini = Reduce("+",Aini.list)/K
  bhat.list = -solve(Ahat.ini)%*%Uini.list + bini;
  list("b.k"=bhat.list,"Ahat"=Ahat.ini)
}
######## One-step estimator functions for Logistic regression ########
#' Inner function for one-step estimator
#'
#' \code{g.logit} evaluates logit function.
#' @export
g.logit = function(xx){
  exp(xx)/(1+exp(xx))
}
#' Inner function for one-step estimator
#'
#' \code{g.logit} evaluates the derivative of logit function.
#' @export
dg.logit = function(xx){
  g.logit(xx)*(1-g.logit(xx))
}

#' Calculate score for logistic regression.
#'
#' \code{U.fun} evaluates scores for logistic regression.
#' @param bet is the beta for evaluation
#' @param dat is the dataset
#' @export
U.fun = function(bet,dat){
  yy = dat[,1]; xx.vec = cbind(1,dat[,-1])
  c(t(c(yy - g.logit(xx.vec%*%bet)))%*%xx.vec)/length(yy)
}

#' Calculate negative information (A) for logistic regression.
#'
#' \code{A.fun} evaluates negative information (A) for logistic regression.
#' @param bet is the beta for evaluation
#' @param dat is the dataset
#' @export
A.fun = function(bet,dat){
  yy = dat[,1]; xx.vec = cbind(1,dat[,-1])
  -t(c(dg.logit(xx.vec%*%bet))*xx.vec)%*%xx.vec/length(yy)
}


#' Calculate log-likelihood for logistic regression.
#'
#' \code{logitlik.fun} evaluates log-likelihood for logistic regression.
#' @param bet.mat is the beta matrix for evaluation
#' @param dat is the dataset
#' @export
logitlik.fun = function(bet.mat,dat){
  yi = dat[,1]; xi = dat[,-1]; pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
  apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
}


#' One-step iteration function for divide-and-conquer logistic model
#'
#' \code{iteration.fun} estimates a one-step for a logistic regression for each subset.
#'
#' @param dat.list list of subsets (after dividing). In each subset, first column = outcome, rest = design matrix
#' @param bini initial estimator as starting point
#' @param kk.list which subsets of dat.list get one-step update
#' @return a list with b.k a matrix of one step estimator and Ahat the negative information matrix
#' @author Yan Wang, Tianxi Cai
#' @export
#' @examples iteration.fun(dat.list=dat.list,bini=bini,kk.list=2:K)
iteration.fun = function(dat.list,bini,kk.list){
  K = length(dat.list)
  Uini.list = sapply(kk.list,function(kk){U.fun(bini,dat=dat.list[[kk]])}) ## p x (K-1) matrix
  Aini.list = lapply(1:K,function(kk){A.fun(bini,dat=dat.list[[kk]])}); Ahat.ini = Reduce("+",Aini.list)/K
  bhat.list = -solve(Ahat.ini)%*%Uini.list + bini;
  list("b.k"=bhat.list,"Ahat"=Ahat.ini)
}
