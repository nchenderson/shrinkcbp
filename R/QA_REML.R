

Q_REML <- function(A, Y, X, sig.sq) {
   ######################################################################
   ## Function that computes \tilde{Q}(A) (A = tau.sq) as defined in the BPE paper
   ## argmax \tilde{Q}(A) gives the value of A used in the BPE estimates
   ##
   ## Input:
   ##    A - the value of A (or, equivalenty A = tau.sq = var(theta))
   ##    Y - the vector (of length K) of responses
   ##    X - the K x p design matrix
   ##    sigma.sq - var(Y_k|\theta_k) = sig.sq
   ##
   ##  Output:
   ##    the value of \tilde{Q}(A)
   ########################################################################


   #K <- length(Y)
   mvar <- sig.sq + A  ## marginal variances
   inv_mvar <- 1/mvar
   #DD <- diag(inv_mvar)
   #XtD <- crossprod(X, DD)
   #XtDX <- XtD%*%X
   XtDX <- crossprod(X, X*inv_mvar)
   XtDy <- crossprod(X, Y*inv_mvar)
   #PP <- diag(rep(1,K)) - X%*%solve(XtD%*%X, XtD)
   #Py <- DD%*%PP%*%Y
   Py <- Y - X%*%solve(XtDX, XtDy)

   t1 <- sum(log(mvar))
   t2 <- determinant(XtDX, logarithm=TRUE)$modulus
   t3 <- sum(Y*inv_mvar*Py)
   ans <- t1 + t2 + t3
   return(ans)
}
