Mhat.naught <- function(tau.sq, nn, yy, sig.sq) {
    K <- length(yy)
    BB <- sig.sq/(sig.sq + nn*tau.sq)
    Bdot <- sum(BB)
    ss <- sig.sq/nn
    ww0 <- (ss/(ss + tau.sq))^2
    ww1 <- 1/(ss + tau.sq)
    ww0 <- ww0/sum(ww0)
    ww1 <- ww1/sum(ww1)
  
    Y.hat0 <- sum(ww0*yy)
    C1 <- ((Bdot/K)*(sum((ww1 - ww0)*yy)))^2
    C2 <- Bdot*mean((ww1 - ww0)*yy)*mean(BB*(yy - Y.hat0)) - ((Bdot*sig.sq)/(K^2))*sum((ww1 - ww0)/nn)
    
    alpha.prem <- C2/C1
    alpha.opt <- min(max(0, alpha.prem), 1)
    ww <- alpha.opt*ww1 + (1 - alpha.opt)*ww0
      
    Y.hat <- sum(ww*yy)
    term1 <- (mean(BB*(yy - Y.hat)))^2
    term2 <- (Bdot*sig.sq*mean(ww/nn))/K
    term3 <- sum((sig.sq*sig.sq)/(nn*(sig.sq + nn*tau.sq)))/(K*K)
    
    ans <- term1 + 2*term2 - 2*term3
    return(ans)
}

Mhat.naught.ftau <- function(alpha, tausq0, tausq1, nn, yy, sig.sq) {
  K <- length(yy)
  tau.sq <- alpha*tausq1 + (1 - alpha)*tausq0
  BB <- sig.sq/(sig.sq + nn*tau.sq)
  Bdot <- sum(BB)
  ss <- sig.sq/nn
  ww0 <- (ss/(ss + tau.sq))^2
  ww1 <- 1/(ss + tau.sq)
  ww0 <- ww0/sum(ww0)
  ww1 <- ww1/sum(ww1)
  
  ww <- alpha*ww1 + (1 - alpha)*ww0
  
  Y.hat <- sum(ww*yy)
  term1 <- (mean(BB*(yy - Y.hat)))^2
  term2 <- (Bdot*sig.sq*mean(ww/nn))/K
  term3 <- sum((sig.sq*sig.sq)/(nn*(sig.sq + nn*tau.sq)))/(K*K)
  
  ans <- term1 + 2*term2 - 2*term3
  return(ans)
}

#tau.seq <- seq(.01, 10, length.out=100)
#ff <- rep(0, 100)
#for(k in 1:100) {
#   ff[k] <- Mhat.naught(tau.sq=tau.seq[k], nn=tmp$n, yy=tmp$Y, sig.sq=sig.sq)
#}
#plot(tau.seq, ff)



