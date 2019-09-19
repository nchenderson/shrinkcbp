
GenSimDatLC <- function(X, bbeta, tau.sq, aa,  u.min=1, u.max=5, re.dist="normal", ncase="exp.unif",sig.sq=1) {
    ### Function to generate data Y under the informative sample size simulation setting
    ### Function to generate data Y and nn under the geometric series pattern for n...

    ## need to figure out u.min and u.max?

    K <- nrow(X)
    tau <- sqrt(tau.sq)
    #nn <- runif(K, min=u.min, max=u.max)
    if(ncase=="exp.unif") {
        #nn <- rgamma(K, shape=1,rate=1/2)
        nn <- exp(runif(K, min=u.min, max=u.max))
        #nn <- rgamma(K, shape=1, rate=.25)
    }
    else {
        nn <- sample(c(u.min, u.max), size=K, replace=TRUE, prob=c(.5, .5))
    }
    z.latent <- rep(0, K) 
    for(k in 1:K) {
        ptmp <- c(1, nn[k], 1/nn[k])
        z.latent[k] <- sample(c(0,1,2), size=1, replace=TRUE, prob=ptmp)
    }
    lin.pred <- X%*%bbeta
    if(re.dist=="normal") {
         eta <- rnorm(K, sd=tau)
         theta.vec <- lin.pred + aa*as.numeric(z.latent==1) - aa*as.numeric(z.latent==2) + eta
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    } else if(re.dist=="normal.mix") {
         zz <- sample(0:1, size=K, replace=TRUE)
         eta <- rep(0.0, K)
         eta[zz==0] <- rnorm(sum(zz==0), mean=-4/sqrt(2), sd=1/sqrt(2))
         eta[zz==1] <- rnorm(sum(zz==1), mean=4/sqrt(2), sd=1/sqrt(2))
         eta <- eta*tau
         theta.vec <- lin.pred + aa*as.numeric(z.latent==1) - aa*as.numeric(z.latent==2) + eta 
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    } else if(re.dist=="tdist") {
         eta <- rt(K, df=3)*(tau/sqrt(3))
         theta.vec <- lin.pred + aa*as.numeric(z.latent==1) - aa*as.numeric(z.latent==2) + eta 
         Y <- rnorm(K, mean=theta.vec, sd=sqrt(sig.sq/nn))
    }
    ans <- list()
    ans$n <- nn
    ans$Y <- Y
    ans$theta.vec <- theta.vec
    return(ans)
}

