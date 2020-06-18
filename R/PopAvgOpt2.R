PopAvgOpt2 <- function(K, sig.sq=25, barn=5, sigma.n=4, gamma=0, xi=1, p=3) {
  
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

    #gamma <- 0.2
    if(p==1) {
        etas <- 5*(n - barn)
    }
    else if(p==2) {
        etas <- 5*(n - barn) - 0.5*(n - barn)^2 
    } else if(p==3) {
        etas <- 5*(n - barn) - 0.5*(n - barn)^2 + 0.1*(n - barn)^3
    }
    if(gamma > 0) {
        z <- gamma*(log(n) + rnorm(K, sd=5))
        rho.ez <- cor(etas, z)
        thetas <- xi*(((etas - mean(etas))/sd(etas) + (z - mean(z))/sd(z))/sqrt(2*(1 + rho.ez)))
    } else if(gamma==0) {
         thetas <- xi*(etas - mean(etas))/sd(etas)
    }
    YY <- thetas + (sig.sq/sqrt(n))*rnorm(K)   ## generate "observations"

    ans <- list(Y=YY, n=n, true.theta=mean(thetas), thetas=thetas)
}



