
MLEWeights <- function(K, tau.sq, sigma.sq) {
    tmp <- 1/(sigma.sq + tau.sq)
    ans <- tmp/sum(tmp)
    return(ans)
}


