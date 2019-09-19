
Q_URE <- function(tau.sq, Y, X, sig.sq) {

  ######################################################################
  ## Function that computes the objective function Q(A) (A = tau.sq) for the unbiased Risk estimate
  ## argmax Q(A) gives the value of the mle
  ##
  ## Input:
  ##    A - the value of A (or, equivalenty A = tau.sq = var(theta))
  ##    Y - the vector (of length K) of responses
  ##    X - the K x p design matrix
  ##    nn - the vector (of length K) of sample sizes
  ##    sigma.sq - var(Y_k) = sig.sq/nn
  ##
  ##  Output:
  ##    the value of Q(A)
  ########################################################################

  ## Define W0 and W1
  ww0 <- 1/(sig.sq + tau.sq)
  ww0 <- ww0/sum(ww0)

  #W <- diag(ww0)
  #K <- nrow(X)

  XWX <- crossprod(X, X*ww0)
  XtW <- t(X*ww0)
  #XtWy <- crossprod(X, ww0*Y)
  Bvec <- sig.sq/(sig.sq + tau.sq)
  #B <- diag(Bvec)

  #UU <- B%*%X%*%solve(XWX, XtW) - B
  #PPy <- X%*%solve(XWX, XtWy)
  PP <- X%*%solve(XWX, XtW)
  PPy <- PP%*%Y
  Uy <- Bvec*PPy - Bvec*Y
  #VY <- diag(sig.sq) ## Var(Y|\theta)
  #TU <- (UU + t(UU))%*%VY

  MSE.hat <- sum(Uy*Uy) + 2*sum(Bvec*sig.sq*diag(PP)) + sum(sig.sq) - 2*sum(Bvec*sig.sq)
  #MSE.hat <- sum(Uy*Uy) + sum(diag(TU)) + sum(sig.sq)
  return(MSE.hat)
}



