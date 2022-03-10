% File src/library/stats/man/Normal.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2014 R Core Team
% Distributed under GPL 2 or later

\name{xkcd}
\title{xkcd}
\alias{xkcd}
\alias{dxkcd}
\alias{pxkcd}
\alias{qxkcd}
\alias{rxkcd}
% These concepts are for the last example
%\concept{error function}
%\concept{erf}
%\concept{erfc}
%\concept{erfinv}
%\concept{erfcinv}
%\title{The Normal Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the xkcd distribution with standard deviation equal to \code{sd}.
}
\usage{
dxkcd(x, sd = 1, log = FALSE,swap.end.points = FALSE)
pxkcd(q, sd = 1, log.p = FALSE,swap.end.points = FALSE)
qxkcd(p, sd = 1, log.p = FALSE,swap.end.points = FALSE)
rxkcd(n, sd = 1 )
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{sd}{vector of standard deviations.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{swap.end.points}{logical; if TRUE, all computations are about the random variable \eqn{W = 1/[(2\pi)^{1/2}\sigma] - Y rather than Y. if FALSE, we simply use Y without transform.}.}
}
\value{
  \code{dxkcd} gives the density,
  \code{pxkcd} gives the distribution function,
  \code{qxkcd} gives the quantile function, and
  \code{rxkcd} generates random deviates.
  
  The length of the result is determined by \code{n} for
  \code{rxkcd}, \code{x,p,q} for corresponding functions.  
  
  The numerical arguments other than \code{n} are recycled to the
  length of the result.  Only the first elements of the logical
  arguments are used.
  
  \code{sd <= 0} or \code{sd = Inf} return an error without a value. 
  \code{NA} input will also return an error.
}
\details{
  \code{sd} has default value of 1, which is open to user specification.
  
  Suppose \code{f} is the probability density function (PDF) of a normal distribution, and we have a random vector \code{(X, Y)} uniformly distributed on the region bounded above by the graph of f and bounded below by the horizontal axis. Then the marginal distribution of \code{X} is this normal distribution, the conditional distribution of \code{Y} given \code{X} is uniform on the interval from zero to \code{f(X)}. The marginal distribution of \code{Y} is the xkcd distribution that we want R functions for.
  
  Let \code{h(y)} be the distance from the means of be the distance from the mean of \code{X} to either of the points where \code{f(x) = y}. Then the distribution function (DF) of \code{Y} is
  \deqn{
    G(y) =1-F(\mu +h(y))+F(\mu-h(y))+2yh(y) ,0<y\le f(\mu),
    }
and Mathematica says the probability density function (PDF) of \code{Y} simplifies to
  \deqn{
    g(y) =2h(y), 0\le y\le f(\mu).
    }
}
%\seealso{
%  \link{Distributions} for other standard distributions, including
%  \code{\link{dlnorm}} for the \emph{Log}normal distribution.
%}


\references{
  The xkcd is distribution introduced by https://xkcd.com/2118/.
}


\examples{
require(graphics)

dxkcd(0) == Inf
dxkcd(1/sqrt(2*pi)) == 0


## Using "log = TRUE" for an extended range :
par(mfrow = c(2,1))
plot(function(x) dxkcd(x, log = TRUE), 0.001, 1/sqrt(2*pi),
     main = "log { xkcd density }")
curve(log(dxkcd(x)), add = TRUE, col = "red", lwd = 2)
mtext("dxkcd(x, log=TRUE)", adj = 0)
mtext("log(dxkcd(x))", col = "red", adj = 1)


plot(function(x) pxkcd(x, log.p = TRUE), 0.001, 1/sqrt(2*pi),
     main = "log { xkcd Cumulative }")
curve(log(pxkcd(x)), add = TRUE, col = "red", lwd = 2)
mtext("pxkcd(x, log=TRUE)", adj = 0)
mtext("log(pxkcd(x))", col = "red", adj = 1)

}
