ShrinkageEstimates <- function(Y, X, sd.resid, weights=c("compromise","compromise.ftau","mle","bpe"), 
                               vc.est=c("mle", "reml","ure")) {
    ### Function to compute estimates of small area means using different weighting schemes
    ###
    ### Input:
    ###    Y - a K x 1 vector of responses
    ###    X - a K x p design matrix
    ###    sd.resid - a scalar or a K x 1 vector of sampling standard deviations
    ###                (i.e., Y_k|\theta_k ~ N(\theta_k, sd.resid_k^2)
    ###    weights - the method for finding weights ....
    ###
    ### Output:
    ###

    weights <- match.arg(weights)

    K <- nrow(X)
    VY <- var(Y)
    sis <- sd.resid*sd.resid

    ests <- switch(weights,
                  compromise = {
                       compromise_pars <- optim(par=c(1/2, VY/4), fn=ComprMSE, lower=c(0,0), upper=c(1, 4*VY), Y=Y, X=X, sig.sq=sis, method="L-BFGS-B")$par
                       alpha.opt <- compromise_pars[1]
                       tausq.opt <- compromise_pars[2]

                       ww.mle <- MLEWeights(K, tau.sq=tausq.opt, sigma.sq=sis)
                       ww.bpe <- BPEWeights(K, tau.sq=tausq.opt, sigma.sq=sis)
                       ww.compromise <- alpha.opt*ww.mle + (1 - alpha.opt)*ww.bpe
                     #  W <- diag(ww.compromise)

                       XWX <- crossprod(X, X*ww.compromise)
                       XtWy <- crossprod(X, ww.compromise*Y)
                       Bvec <- sis/(sis + tausq.opt)

                      # XWX <- crossprod(X, W%*%X)
                      # XtW <- crossprod(X, W)
                      # Bvec <- sis/(sis + tausq.opt)
                      # B <- diag(Bvec)

                       #Amat <- B%*%X%*%solve(XWX, XtW) + diag(rep(1,K)) - B
                       #theta.hat <- Amat%*%Y
                       PPy <- X%*%solve(XWX, XtWy)
                       theta.hat <- Bvec*PPy + (1 - Bvec)*Y
                       ests <- list(shrinkage.estimate=theta.hat, weights=ww.compromise, alpha=alpha.opt, tau.sq=tausq.opt)
                  }, compromise.ftau = {
                       tausq.reml <- optimize(Q_REML, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                       tausq.bpe <- optimize(Q_BPE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum

                       compromise_pars <- optim(par=c(1/2, VY/4), fn=ComprMSE_fixedtau, lower=c(0,0), upper=c(1, 4*VY), Y=Y, X=X, tau.sq0=tausq.reml,
                                                tau.sq1=tausq.bpe, sig.sq=sis, method="L-BFGS-B")$par
                       alpha.opt <- compromise_pars[1]
                       tausq.opt <- compromise_pars[2]

                       ww.mle <- MLEWeights(K, tau.sq=tausq.reml, sigma.sq=sis)
                       ww.bpe <- BPEWeights(K, tau.sq=tausq.bpe, sigma.sq=sis)
                       ww.compromise <- alpha.opt*ww.mle + (1 - alpha.opt)*ww.bpe
                       #W <- diag(ww.compromise)

                       XWX <- crossprod(X, X*ww.compromise)
                       XtWy <- crossprod(X, ww.compromise*Y)
                       Bvec <- sis/(sis + tausq.opt)
                       PPy <- X%*%solve(XWX, XtWy)
                       theta.hat <- Bvec*PPy + (1 - Bvec)*Y

                       #XWX <- crossprod(X, W%*%X)
                       #XtW <- crossprod(X, W)
                       #Bvec <- sis/(sis + tausq.opt)
                       #B <- diag(Bvec)

                       #Amat <- B%*%X%*%solve(XWX, XtW) + diag(rep(1,K)) - B
                       #theta.hat <- Amat%*%Y
                       ests <- list(shrinkage.estimate=theta.hat, weights=ww.compromise, alpha=alpha.opt, tau.sq=tausq.opt)
                   }, mle = {
                       if(vc.est=="mle") {
                           mle.tausq <- optimize(Q_MLE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                       } else if(vc.est=="reml") {
                           mle.tausq <- optimize(Q_REML, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                       } else if(vc.est=="ure") {
                           mle.tausq <- optimize(Q_URE, interval=c(0, 4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                       }
                       else if(vc.est=="fay.herriot") {
                           mle.tausq <- 1
                       }
                       ## else if (vc.est=="sure")
                       ww.mle <- MLEWeights(K, tau.sq=mle.tausq, sigma.sq=sis)
                       #W <- diag(ww.mle)

                      # XWX <- crossprod(X, W%*%X)
                      # XtW <- crossprod(X, W)
                       XWX <- crossprod(X, X*ww.mle)
                       XtWy <- crossprod(X, ww.mle*Y)
                       Bvec <- sis/(sis + mle.tausq)
                      # B <- diag(Bvec)

                       PPy <- X%*%solve(XWX, XtWy)

                       #Amat <- B%*%X%*%solve(XWX, XtW) + diag(rep(1,K)) - B
                       #theta.hat <- Amat%*%Y
                       theta.hat <- Bvec*PPy + (1 - Bvec)*Y
                       ests <- list(shrinkage.estimate=theta.hat, weights=ww.mle, alpha=0, tau.sq=mle.tausq)
                  }, bpe = {
                       bpe.tausq <- optimize(Q_BPE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                       ww.bpe <- BPEWeights(K, tau.sq=bpe.tausq, sigma.sq=sis)
                       #W <- diag(ww.bpe)

                       XWX <- crossprod(X, X*ww.bpe)
                       XtWy <- crossprod(X, ww.bpe*Y)
                       Bvec <- sis/(sis + bpe.tausq)

                       #XWX <- crossprod(X, W%*%X)
                       #XtW <- crossprod(X, W)
                       #Bvec <- sis/(sis + bpe.tausq)
                       #B <- diag(Bvec)

                       PPy <- X%*%solve(XWX, XtWy)

                       #Amat <- B%*%X%*%solve(XWX, XtW) + diag(rep(1,K)) - B
                       #theta.hat <- Amat%*%Y
                       theta.hat <- Bvec*PPy + (1 - Bvec)*Y
                       ests <- list(shrinkage.estimate=theta.hat, weights=ww.bpe, alpha=1, tau.sq=bpe.tausq)
                  })
    return(ests)
}
