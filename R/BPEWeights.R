
BPEWeights <- function(K, tau.sq, sigma.sq) {
    dd <- sigma.sq/(sigma.sq + tau.sq)
    tmp <- dd*dd
    ans <- tmp/sum(tmp)
    return(ans)
}


