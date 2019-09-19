

Q_BPE <- function(A, Y, X, sig.sq) {
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


   gamma.vec <- sig.sq/(A + sig.sq)
   gamma.vec.sq <- gamma.vec*gamma.vec
   #GammaMat <- diag(gamma.vec)
   #GammaMat2 <- diag(gamma.vec*gamma.vec)
   #XGX <- crossprod(X, GammaMat2%*%X)
   #XGY <- crossprod(X, GammaMat2%*%Y)
   XGX <- crossprod(X, X*gamma.vec.sq)
   XGY <- crossprod(X, gamma.vec.sq*Y)
   bpe <- solve(XGX, XGY)
   #fitted_vals <- X%*%bpe
   #resids <- Y - fitted_vals
   resids <- Y - X%*%bpe

   ans <- sum(resids*resids*gamma.vec.sq) + 2*A*sum(gamma.vec)
   return(ans)
}
