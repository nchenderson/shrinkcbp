\name{ShrinkageEstimates}
\alias{ShrinkageEstimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shrinkage Estimation
}
\description{
Computes shrinkage estimates for different weighting schemes
}
\usage{
ShrinkageEstimates(Y, X, sd.resid, weights = c("compromise", "mle", "bpe"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
the vector of responses
}
  \item{X}{
the design matrix
}
  \item{sd.resid}{
the sampling standard deviations
}
  \item{weights}{
the method by which to choose weights
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
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


a <- ShrinkageEstimates(Y, X, sd.resid)
b <- ShrinkageEstimates(Y, X, sd.resid, weights="mle")
d <- ShrinkageEstimates(Y, X, sd.resid, weights="bpe")

plot(a$shrinkage.estimate, b$shrinkage.estimate)
abline(0, 1)





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
