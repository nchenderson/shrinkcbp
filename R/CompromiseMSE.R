
ComprMSE <- function(par, Y, X, sig.sq) {
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
  tau.sq <- par[2]

  ## Define W0 and W1
  ww0 <- 1/(sig.sq + tau.sq)
  ww1 <- (sig.sq/(sig.sq + tau.sq))^2
  ww0 <- ww0/sum(ww0)
  ww1 <- ww1/sum(ww1)

  #W0 <- diag(ww0)
  #W1 <- diag(ww1)

  #K <- nrow(X)
  #W <- alpha*W0 + (1 - alpha)*W1
  #XWX <- crossprod(X, W%*%X)
  #XtW <- crossprod(X, W)
  Bvec <- sig.sq/(sig.sq + tau.sq)
  B <- diag(Bvec)

  ww.compr <- alpha*ww0 + (1 - alpha)*ww1
  XWX <- crossprod(X, X*ww.compr)
  XtW <- t(X*ww.compr)
  PP <- X%*%solve(XWX, XtW)
  PPy <- PP%*%Y
  Uy <- Bvec*PPy - Bvec*Y

 
  #UU <- B%*%PP - B
  #Uy <- drop(UU%*%Y)
  #VY <- diag(sig.sq) ## Var(Y|\theta)
  #TU <- (UU + t(UU) + diag(rep(1,K)))%*%VY

  #TU <- UU%*%diag(sig.sq)
  #MSE.hat <- sum(Uy*Uy) + 2*sum(diag(TU))
  MSE.hat <- sum(Uy*Uy) + 2*sum(Bvec*sig.sq*diag(PP)) + sum(sig.sq) - 2*sum(Bvec*sig.sq)

  return(MSE.hat)
}



