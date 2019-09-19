
Q_MLE <- function(A, Y, X, sig.sq) {
   ######################################################################
   ## Function that computes the objective function Q(A) (A = tau.sq) for the MLE
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

   #K <- length(Y)
   mvar <- sig.sq + A  ## marginal variances
   inv_mvar <- 1/mvar
  # DD <- diag(inv_mvar)
   XtDy <- crossprod(X, Y*inv_mvar)
   XtDX <- crossprod(X, X*inv_mvar)
   #PP <- diag(rep(1,K)) - X%*%solve(XtD%*%X, XtD)
   Py <- Y - X%*%solve(XtDX, XtDy)

   ans <- sum(log(mvar)) + sum(Py*Py*inv_mvar)
   return(ans)
}

