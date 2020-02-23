ComprMSE_fixedtau2 <- function(par, Y, X, sig.sq, tau.sq0, tau.sq1) {
  #############################################################
  ## Function to compute the sum of unit-specific MSEs for the compromise estimator between MLE and BPE
  ##
  ## Input:
  ##   par - vector of length 2 - (alpha, tau.sq)
  ##   Y - the vector (of length K) of responses
  ##   X - the K x p design matrix
  ##   W0 - diag(w0) of the w(0) weighting scheme
  ##   W1 - diag(w1) of the w(1) weighting scheme
  ##   sigma.sq - var(Y_k|\theta_k)
  ##  Output:
  ##     The estimated MSE for this choice of (alpha, tau.sq)
  ###################################################################
  alpha <- par[1]

  ## Define W0 and W1
  ww0 <- 1/(sig.sq + tau.sq0)
  ww1 <- (sig.sq/(sig.sq + tau.sq1))^2
  ww0 <- ww0/sum(ww0)
  ww1 <- ww1/sum(ww1)
  
  W0 <- diag(ww0)
  W1 <- diag(ww1)
  
  K <- nrow(X)
  W <- alpha*W0 + (1 - alpha)*W1
  XWX <- crossprod(X, W%*%X)
  XtW <- crossprod(X, W)
  tau.sq.comp <- alpha*tau.sq0 + (1 - alpha)*tau.sq1
  Bvec <- sig.sq/(sig.sq + tau.sq.comp)
  B <- diag(Bvec)
  
  UU <- B%*%X%*%solve(XWX, XtW) - B
  Uy <- UU%*%Y
  VY <- diag(sig.sq) ## Var(Y|\theta)
  TU <- (UU + t(UU) + diag(rep(1,K)))%*%VY
  
  MSE.hat <- sum(Uy*Uy) + sum(diag(TU))
  return(MSE.hat)
}
