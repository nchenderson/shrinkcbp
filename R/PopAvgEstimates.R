PopAvgEstimates <- function(Y, sigma, nn, weights=c("compromise","compromise.ftau", "compromise.ftau2", "compromise.difftau", "mle","bpe"), 
                            vc.est=c("mle", "reml","ure")) {
  
   weights <- match.arg(weights)
  
   K <- length(Y)
   X <- matrix(1, nrow=K, ncol=1)
   VY <- var(Y)
   sis <- (sigma*sigma)/nn
   sig.sq <- sigma*sigma
  
   ests <- switch(weights,
                  compromise = {
                     best.tausq <- optimize(Mhat.naught, interval=c(0, 4*VY), nn=nn, yy=Y, sig.sq=sig.sq)$minimum 
                     
                     K <- length(Y)
                     BB <- sig.sq/(sig.sq + nn*best.tausq)
                     Bdot <- sum(BB)
                     ss <- sig.sq/nn
                     ww0 <- (ss/(ss + best.tausq))^2
                     ww1 <- 1/(ss + best.tausq)
                     ww0 <- ww0/sum(ww0)
                     ww1 <- ww1/sum(ww1)
                     
                     Y.hat0 <- sum(ww0*Y)
                     C1 <- ((Bdot/K)*(sum((ww1 - ww0)*Y)))^2
                     C2 <- Bdot*mean((ww1 - ww0)*Y)*mean(BB*(Y - Y.hat0)) - ((Bdot*sig.sq)/(K^2))*sum((ww1 - ww0)/nn)
                     
                     alpha.prem <- C2/C1
                     alpha.opt <- min(max(0, alpha.prem), 1)
                     
                     ww.compromise <- alpha.opt*ww1 + (1 - alpha.opt)*ww0
                     Y.hat <- sum(ww.compromise*Y)
                     theta.hat <- (Bdot*Y.hat)/K + mean((1 - BB)*Y)
                     
                     #compromise_pars <- optim(par=c(1/2, VY/4), fn=ComprMSE, lower=c(0,0), upper=c(1, 4*VY), Y=Y, X=X, sig.sq=sis, method="L-BFGS-B")$par
                     #alpha.opt <- compromise_pars[1]
                     #tausq.opt <- compromise_pars[2]
                     
                     #ww.bpe <- BPEWeights(K, tau.sq=tausq.opt, sigma.sq=sis)
                     #ww.compromise <- alpha.opt*ww.mle + (1 - alpha.opt)*ww.bpe

                     #mu.hat <- sum(ww.compromise*Y)/sum(ww.compromise)
                     #Bvec <- sis/(sis + tausq.opt)
                     #theta.hats <- Bvec*mu.hat + (1 - Bvec)*Y
                     
                     #theta.hat <- mean(theta.hats)
                     ests <- list(shrinkage.estimate=theta.hat, weights=ww.compromise, alpha=alpha.opt, tau.sq=best.tausq)
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
                      
                      mu.hat <- sum(ww.compromise*Y)/sum(ww.compromise)
                      Bvec <- sis/(sis + tausq.opt)
                      theta.hats <- Bvec*mu.hat + (1 - Bvec)*Y
                     
                      theta.hat <- mean(theta.hats)
                      ests <- list(shrinkage.estimate=theta.hat, weights=ww.compromise, alpha=alpha.opt, tau.sq=tausq.opt)
                  }, compromise.ftau2 = {
                      tausq.reml <- optimize(Q_REML, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      tausq.bpe <- optimize(Q_BPE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                     
                      
                      #Mhat.naught.ftau <- function(alpha, tausq0, tausq1, nn, yy, sig.sq) {
                      compromise_pars <- optimize(f=Mhat.naught.ftau, interval=c(0,1), yy=Y, nn=nn, sig.sq=sig.sq, tausq0=tausq.bpe,
                                                  tausq1=tausq.reml)
                      alpha.opt <- compromise_pars$minimum
                      tausq.opt <- alpha.opt*tausq.reml + (1 - alpha.opt)*tausq.bpe
                     
                      ss <- sig.sq/nn
                      ww0 <- (ss/(ss + tausq.opt))^2
                      ww1 <- 1/(ss + tausq.opt)
                      ww.bpe <- ww0/sum(ww0)
                      ww.mle <- ww1/sum(ww1)
                      BB <- sig.sq/(sig.sq + nn*tausq.opt)
                      Bdot <- sum(BB)
                      
                      ww.compromise <- alpha.opt*ww.mle + (1 - alpha.opt)*ww.bpe
                      Y.hat <- sum(ww.compromise*Y)
                      theta.hat <- (Bdot*Y.hat)/K + mean((1 - BB)*Y)
                      
                      ests <- list(shrinkage.estimate=theta.hat, weights=ww.compromise, alpha=alpha.opt, tau.sq=tausq.opt)
                  }, compromise.difftau = {
                      tausq.reml <- optimize(Q_REML, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      tausq.bpe <- optimize(Q_BPE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      compromise_pars <- optim(par=c(1/2, tausq.bpe, tausq.reml), fn=ComprDiffMSE, lower=c(0,0), upper=c(1, 4*VY), Y=Y, X=X, sig.sq=sis, method="L-BFGS-B")$par
                      alpha.opt <- compromise_pars[1]
                      tausq.opt0 <- compromise_pars[2]
                      tausq.opt1 <- compromise_pars[3]
                      #tausq.shrink <- compromise_pars[4]
                      tausq.shrink <- alpha.opt*tausq.opt1 + (1 - alpha.opt)*tausq.opt0
                     
                      ww.mle <- MLEWeights(K, tau.sq=tausq.opt1, sigma.sq=sis)
                      ww.bpe <- BPEWeights(K, tau.sq=tausq.opt0, sigma.sq=sis)
                      ww.compromise <- alpha.opt*ww.mle + (1 - alpha.opt)*ww.bpe
                      
                      mu.hat <- sum(ww.compromise*Y)/sum(ww.compromise)
                      Bvec <- sis/(sis + tausq.shrink)
                      theta.hats <- Bvec*mu.hat + (1 - Bvec)*Y
                      theta.hat <- mean(theta.hats)
                      ests <- list(shrinkage.estimate=theta.hat, weights=ww.compromise, alpha=alpha.opt, tau.sq=tausq.shrink)
                  }, mle = {
                      if(vc.est=="mle") {
                         mle.tausq <- optimize(Q_MLE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      } else if(vc.est=="reml") {
                         mle.tausq <- optimize(Q_REML, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      } else if(vc.est=="ure") {
                        mle.tausq <- optimize(Q_URE, interval=c(0, 4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      }
                      ww.mle <- MLEWeights(K, tau.sq=mle.tausq, sigma.sq=sis)
                     
                      mu.hat <- sum(ww.mle*Y)/sum(ww.mle)
                      Bvec <- sis/(sis + mle.tausq)
                      theta.hats <- Bvec*mu.hat + (1 - Bvec)*Y
                      
                      theta.hat <- mean(theta.hats)
                      ests <- list(shrinkage.estimate=theta.hat, weights=ww.mle, alpha=0, tau.sq=mle.tausq)
                  }, bpe = {
                      bpe.tausq <- optimize(Q_BPE, interval=c(0,4*VY), Y=Y, X=X, sig.sq=sis)$minimum
                      ww.bpe <- BPEWeights(K, tau.sq=bpe.tausq, sigma.sq=sis)
                     
                      mu.hat <- sum(ww.bpe*Y)/sum(ww.bpe)
                      Bvec <- sis/(sis + bpe.tausq)
                      theta.hats <- Bvec*mu.hat + (1 - Bvec)*Y
                      theta.hat <- mean(theta.hats)
                      ests <- list(shrinkage.estimate=theta.hat, weights=ww.bpe, alpha=1, tau.sq=bpe.tausq)
                  })
   return(ests)  
}