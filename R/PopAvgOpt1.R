# Weighted average of estimates approach and regression adjustment approach
# Numerical evaluation



PopAvgOpt1 <- function(K, sig.sq=25, barn=5, sigma.n=4, rho=.25, xi=1) {
  ### sis is sigma_sq

  ff <- function(aa) {
      tmpU <- exp(seq(-aa,  aa, length = K))
      tmpn <- (K*barn*tmpU)/sum(tmpU)
      ans <- sd(tmpn) - sigma.n
      return(ans)
  }
  
  tmp <- uniroot(ff, lower=0.5, upper=20)
  U <- exp(seq(-tmp$root, tmp$root, length = K))
  n <- (K*barn*U)/sum(U)
  hh <- sd(log(n))/2
  
  ff <- pnorm(log(n), mean=mean(log(n)), sd=hh) - 1/2
  Vf <- var(ff)
  kappa <- mean(ff*n)
  rho.sq <- rho*rho

  c1 <- (xi*rho*sd(n))/kappa
  c2 <- sqrt(xi*xi - Vf*c1*c1)
  theta <- c1*ff + c2*rnorm(K)
  YY <- theta + (sig.sq/sqrt(n))*rnorm(K)   ## generate "observations"
  
  ans <- list(Y=YY, n=n, true.theta=mean(theta), thetas=theta)
  return(ans)
}
 






