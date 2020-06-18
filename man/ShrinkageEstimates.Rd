\name{ShrinkageEstimates}
\alias{ShrinkageEstimates}
\title{
Shrinkage Estimation of Small Domain Parameters
}
\description{
Computes shrinkage estimates for small domain parameters using different choices for the regression weights.
}
\usage{
ShrinkageEstimates(Y, X, sd.resid, weights=c("cure", "cure.plugin", 
                  "cure.multitau", "mle","bpe"), vc.est=c("mle", "reml","ure"))
}
\arguments{
  \item{Y}{
The K x 1 vector of responses.
}
  \item{X}{
The K x p design matrix.
}
  \item{sd.resid}{
The K x 1 vector sampling standard deviations.
}
  \item{weights}{
The method by which to choose regression weights. 
}
  \item{vc.est}{
When the mle method is used for regression weights, this determines which method is used to estimate the variance component. This can be either maximum likelihood estimation, restricted maximum likelihood estimation, or minimization of an unbiased risk estimate.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following components
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
K <- 20
X <- cbind(rep(1, K), rnorm(K))
beta.coef <- rnorm(2)
nn <- 1 + rpois(K, lambda=4)
sig.sq <- 1/nn
ran.ef <- rnorm(K)

Y <- X\%*\%beta.coef + ran.ef + sqrt(sig.sq)*rnorm(K)
sd.resid <- sqrt(sig.sq)


cbp.est <- ShrinkageEstimates(Y, X, sd.resid, weights="cure")
eblup.est <- ShrinkageEstimates(Y, X, sd.resid, weights="mle")
obp.est <- ShrinkageEstimates(Y, X, sd.resid, weights="bpe")

par(mfrow=c(1,2))
plot(cbp.est$shrinkage.estimate, eblup.est$shrinkage.estimate)
plot(cbp.est$shrinkage.estimate, obp.est$shrinkage.estimate)
abline(0, 1)

}
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ robust }% __ONLY ONE__ keyword per line
