############# Cross-validation functions #################
#' Inner function for Est.ALASSO.GLMNET.CV
#'
#' \code{mycv.glmnet} is a modified function of cv.glmnet, such that the prediction step can be split and run.
#' @export
mycv.glmnet= function (x, y, weights, offset = NULL, lambda = NULL, type.measure = c("mse",
                                                                                     "deviance", "class", "auc", "mae"), nfolds = 10, foldid,
                       grouped = TRUE, keep = FALSE, parallel = FALSE, ...)
{
  if (missing(type.measure))
    type.measure = "default"
  else type.measure = match.arg(type.measure)
  if (!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv.glmnet")
  N = nrow(x)
  if (missing(weights))
    weights = rep(1, N)
  else weights = as.double(weights)
  y = drop(y)
  glmnet.call = match.call(expand.dots = TRUE)
  which = match(c("type.measure", "nfolds", "foldid", "grouped",
                  "keep"), names(glmnet.call), F)
  if (any(which))
    glmnet.call = glmnet.call[-which]
  glmnet.call[[1]] = as.name("glmnet")
  glmnet.object = glmnet(x, y, weights = weights, offset = offset,
                         lambda = lambda, ...)
  gc()
  glmnet.object$call = glmnet.call
  is.offset = glmnet.object$offset
  if (inherits(glmnet.object, "multnet") && !glmnet.object$grouped) {
    nz = predict(glmnet.object, type = "nonzero")
    nz = sapply(nz, function(x) sapply(x, length))
    nz = ceiling(apply(nz, 1, median))
  }
  else nz = sapply(predict(glmnet.object, type = "nonzero"),
                   length)
  if (missing(foldid))
    foldid = sample(rep(seq(nfolds), length = N))
  else nfolds = max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
    {
      which = foldid == i
      if (is.matrix(y))
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset)
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      glmnet(x[!which, , drop = FALSE], y_sub, lambda = lambda,
             offset = offset_sub, weights = weights[!which],
             ...)
    }
  }
  else {
    for (i in seq(nfolds)) {
      which = foldid == i
      if (is.matrix(y))
        y_sub = y[!which, ]
      else y_sub = y[!which]
      if (is.offset)
        offset_sub = as.matrix(offset)[!which, ]
      else offset_sub = NULL
      outlist[[i]] = glmnet(x[!which, , drop = FALSE],
                            y_sub, lambda = lambda, offset = offset_sub,
                            weights = weights[!which], ...)
      gc()
    }
  }
  return(list(outlist=outlist, lambda=lambda, x=x, y=y, weights=weights,
              offset=offset, foldid=foldid, type.measure=type.measure, grouped=grouped, keep=keep,
              glmnet.object=glmnet.object,nz=nz))
}
#' Inner function for Est.ALASSO.GLMNET.CV
#'
#' \code{mycv.lognet} is a modified function of cv.lognet, such that the prediction step can be split and run.
#' @export
mycv.lognet <- function (outlist, lambda, x, y, weights, offset, foldid, type.measure,
                         grouped, keep = FALSE)
{
  typenames = c(mse = "Mean-Squared Error", mae = "Mean Absolute Error",
                deviance = "Binomial Deviance", auc = "AUC", class = "Misclassification Error")
  if (type.measure == "default")
    type.measure = "deviance"
  if (!match(type.measure, c("mse", "mae", "deviance", "auc",
                             "class"), FALSE)) {
    warning("Only 'deviance', 'class', 'auc', 'mse' or 'mae'  available for binomial models; 'deviance' used")
    type.measure = "deviance"
  }
  prob_min = 1e-05
  prob_max = 1 - prob_min
  nc = dim(y)
  if (is.null(nc)) {
    y = as.factor(y)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ]
  }
  N = nrow(y)
  nfolds = max(foldid)
  if ((N/nfolds < 10) && type.measure == "auc") {
    warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",
            call. = FALSE)
    type.measure = "deviance"
  }
  if ((N/nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold",
            call. = FALSE)
    grouped = FALSE
  }
  if (!is.null(offset)) {
    is.offset = TRUE
    offset = drop(offset)
  }
  else is.offset = FALSE
  mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
  which_lam = lambda >= mlami
  predmat = matrix(NA, nrow(y), length(lambda))
  nlams = double(nfolds)
  for (i in seq(nfolds)) {
    which = foldid == i
    fitobj = outlist[[i]]
    if (is.offset)
      off_sub = offset[which]
    preds = predict(fitobj, x[which, , drop = FALSE], s = lambda[which_lam],
                    offset = off_sub, type = "response")
    gc()
    nlami = sum(which_lam)
    predmat[which, seq(nlami)] = preds
    rm(preds);gc()
    nlams[i] = nlami
  }
  if (type.measure == "auc") {
    cvraw = matrix(NA, nfolds, length(lambda))
    good = matrix(0, nfolds, length(lambda))
    for (i in seq(nfolds)) {
      good[i, seq(nlams[i])] = 1
      which = foldid == i
      for (j in seq(nlams[i])) {
        cvraw[i, j] = auc.mat(y[which, ], predmat[which,
                                                  j], weights[which])
      }
    }
    N = apply(good, 2, sum)
    weights = tapply(weights, foldid, sum)
  }
  else {
    ywt = apply(y, 1, sum)
    y = y/ywt
    weights = weights * ywt
    N = nrow(y) - apply(is.na(predmat), 2, sum)
    cvraw = switch(type.measure, mse = (y[, 1] - (1 - predmat))^2 +
                     (y[, 2] - predmat)^2, mae = abs(y[, 1] - (1 - predmat)) +
                     abs(y[, 2] - predmat), deviance = {
                       predmat = pmin(pmax(predmat, prob_min), prob_max)
                       lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                       ly = log(y)
                       ly[y == 0] = 0
                       ly = drop((y * ly) %*% c(1, 1))
                       2 * (ly - lp)
                     }, class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <=
                                                                       0.5))
    if (grouped) {
      cvob = cvcompute(cvraw, weights, foldid, nlams)
      cvraw = cvob$cvraw
      weights = cvob$weights
      N = cvob$N
    }
  }
  cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                    w = weights, na.rm = TRUE)/(N - 1))
  out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
  if (keep)
    out$fit.preval = predmat
  out
}
#' Inner function for Est.ALASSO.GLMNET.CV
#'
#' \code{mycv.coxnet} is a modified function of cv.coxnet, such that the prediction step can be split and run.
#' @export
mycv.coxnet = function (outlist, lambda, x, y, weights, offset, foldid, type.measure,
                        grouped, keep = FALSE)
{
  typenames = c(deviance = "Partial Likelihood Deviance")
  if (type.measure == "default")
    type.measure = "deviance"
  if (!match(type.measure, c("deviance"), FALSE)) {
    warning("Only 'deviance'  available for Cox models; changed to type.measure='deviance'")
    type.measure = "deviance"
  }
  if (!is.null(offset)) {
    is.offset = TRUE
    offset = drop(offset)
  }
  else is.offset = FALSE
  nfolds = max(foldid)
  if ((length(weights)/nfolds < 10) && !grouped) {
    warning("Option grouped=TRUE enforced for cv.coxnet, since < 3 observations per fold",
            call. = FALSE)
    grouped = TRUE
  }
  cvraw = matrix(NA, nfolds, length(lambda))
  mlami = max(sapply(outlist, function(obj) min(obj$lambda)))
  which_lam = lambda >= mlami
  for (i in seq(nfolds)) {
    which = foldid == i
    which_lam = lambda >= mlami
    fitobj = outlist[[i]]
    coefmat = predict(fitobj, type = "coeff", s = lambda[which_lam])
    which_lam = !is.na(colMeans(as.matrix(coefmat)))
    if (grouped) {
      plfull = coxnet.deviance(x = x, y = y, offset = offset,
                               weights = weights, beta = coefmat[,which_lam])
      plminusk = coxnet.deviance(x = x[!which, ], y = y[!which, ],
                                 offset = offset[!which], weights = weights[!which],
                                 beta = coefmat[,which_lam])
      cvraw[i, (1:length(lambda))[which_lam][seq(along = plfull)]] = plfull - plminusk
    }
    else {
      plk = coxnet.deviance(x = x[which, ], y = y[which,
                                                  ], offset = offset[which], weights = weights[which],
                            beta = coefmat)
      cvraw[i, seq(along = plk)] = plk
    }
  }
  status = y[, "status"]
  N = nfolds - apply(is.na(cvraw), 2, sum)
  weights = as.vector(tapply(weights * status, foldid, sum))
  cvraw = cvraw/weights
  cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
  cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                    w = weights, na.rm = TRUE)/(N - 1))
  out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
  if (keep)
    warning("keep=TRUE not implemented for coxnet")
  out
}
