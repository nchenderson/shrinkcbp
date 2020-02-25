
GenSimDatISS <- function(X, bbeta, tau.sq, rho, n.seq, sd.eta = 1, re.dist="normal", sig.sq=1) {
    ### Function to generate data Y under the informative sample size simulation setting
    ### Function to generate data Y and nn under the geometric series pattern for n...

    ## need to figure out u.min and u.max?

    K <- nrow(X)
    #nn <- runif(K, min=u.min, max=u.max)
    #if(ncase=="exp.unif") {
        #nn <- rgamma(K, shape=1,rate=1/2)
    #    nn <- exp(runif(K, min=u.min, max=u.max))
    #    sigma.nn <- sd(nn)  ## change this later
    #}
    #else {
    #    nn <- sample(c(u.min, u.max), size=K, replace=TRUE, prob=c(.5, .5))
    #    sigma.nn <- (u.max - u.min)/2
    #}
    nn <- n.seq
    sigma.nn <- sd(nn)
    bb <- (rho*sqrt(tau.sq))/sigma.nn
    aa <- sqrt(tau.sq)*sqrt(1 - rho*rho)
    lin.pred <- X%*%bbeta
    if(re.dist=="normal") {
         eta <- rnorm(K, sd=sd.eta)
         theta.vec <- lin.pred + aa*eta + bb*nn
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    } else if(re.dist=="normal.mix") {
         zz <- sample(0:1, size=K, replace=TRUE)
         eta <- rep(0.0, K)
         eta[zz==0] <- rnorm(sum(zz==0), mean=-4/sqrt(2), sd=1/sqrt(2))
         eta[zz==1] <- rnorm(sum(zz==1), mean=4/sqrt(2), sd=1/sqrt(2))
         eta <- eta*sd.eta
         theta.vec <- lin.pred + aa*eta + bb*nn
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    } else if(re.dist=="tdist") {
         eta <- rt(K, df=3)*(sd.eta/sqrt(3))
         theta.vec <- lin.pred + aa*eta + bb*nn
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    } else if(re.dist=="uniform") {
         eta <- runif(K, min=-sqrt(3), max=sqrt(3))
         theta.vec <- lin.pred + aa*eta + bb*nn
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    }
    rho.hat <- cor(theta.vec, nn)
    
    ans <- list()
    ans$n <- nn
    ans$Y <- Y
    ans$theta.vec <- theta.vec
    return(ans)
}

