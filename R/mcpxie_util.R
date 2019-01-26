############# MCP Xie modified function #################
#' Inner function for MCP.Xie
#'
#' \code{ncvsurv2} is a modified function of ncvsurv (to resolve the bug from lamNames)
#' @export
ncvsurv2 = function (X, y, penalty = c("MCP", "SCAD", "lasso"), gamma = switch(penalty,
                                                                    SCAD = 3.7, 3), alpha = 1, lambda.min = ifelse(n > p, 0.001,
                                                                                                                   0.05), nlambda = 100, lambda, eps = 1e-04, max.iter = 10000,
          convex = TRUE, dfmax = p, penalty.factor = rep(1, ncol(X)),
          warn = TRUE, returnX = FALSE, ...)
{
  penalty <- match.arg(penalty)
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0 + ., data = X), silent = TRUE)
    if (class(tmp)[1] == "try-error")
      stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X) == "integer")
    storage.mode(X) <- "double"
  if (class(y) != "matrix") {
    tmp <- try(y <- as.matrix(y), silent = TRUE)
    if (class(tmp)[1] == "try-error")
      stop("y must be a matrix or able to be coerced to a matrix")
    if (ncol(y) != 2)
      stop("y must have two columns for survival data: time-on-study and a censoring indicator")
  }
  if (storage.mode(y) == "integer")
    storage.mode(y) <- "double"
  if (storage.mode(penalty.factor) != "double")
    storage.mode(penalty.factor) <- "double"
  if (gamma <= 1 & penalty == "MCP")
    stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty == "SCAD")
    stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2)
    stop("nlambda must be at least 2")
  if (alpha <= 0)
    stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor) != ncol(X))
    stop("penalty.factor does not match up with X")
  if (any(is.na(y)) | any(is.na(X)))
    stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")
  tOrder <- order(y[, 1])
  yy <- as.numeric(y[tOrder, 1])
  Delta <- y[tOrder, 2]
  n <- length(yy)
  XX <- std(X[tOrder, , drop = FALSE])
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min,
                             nlambda, penalty.factor)
    user.lambda <- FALSE
  }
  else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  res <- .Call(ncvreg:::"cdfit_cox_dh", XX, Delta, penalty, lambda,
               eps, as.integer(max.iter), as.double(gamma), penalty.factor,
               alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor ==
                                                                        0)), as.integer(warn))
  b <- matrix(res[[1]], p, nlambda)
  loss <- -1 * res[[2]]
  iter <- res[[3]]
  Eta <- matrix(res[[4]], n, nlambda)
  ind <- !is.na(iter)
  b <- b[, ind, drop = FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind, drop = FALSE]
  if (warn & sum(iter) == max.iter)
    warning("Algorithm failed to converge for some values of lambda")
  convex.min <- if (convex)
    ncvreg:::convexMin(b, XX, penalty, gamma, lambda * (1 - alpha),
              "cox", penalty.factor, Delta = Delta)
  else NULL
  beta <- matrix(0, nrow = ncol(X), ncol = length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns, ] <- bb
  offset <- -crossprod(attr(XX, "center")[ns], bb)
  varnames <- if (is.null(colnames(X)))
    paste("V", 1:ncol(X), sep = "")
  else colnames(X)
  dimnames(beta) <- list(varnames, as.character(lambda))
  val <- structure(list(beta = beta, iter = iter, lambda = lambda,
                        penalty = penalty, gamma = gamma, alpha = alpha, convex.min = convex.min,
                        loss = loss, penalty.factor = penalty.factor, n = n,
                        time = yy, fail = Delta, order = tOrder), class = c("ncvsurv",
                                                                            "ncvreg"))
  val$Eta <- sweep(Eta, 2, offset, "-")
  if (returnX) {
    val$X <- XX
    val$y <- yy
  }
  val
}
