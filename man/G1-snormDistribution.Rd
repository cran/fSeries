\name{snormDistribution}

\alias{snormDistribution}

\alias{dsnorm}
\alias{psnorm}
\alias{qsnorm}
\alias{rsnorm}


\title{Skew Normal Distribution}


\description{
    
    A collection of functions to compute density, distribution function, 
    quantile function and to generate random variates for the skew
    normal distribution. 
    \cr
    
    The functions are:
    
    \tabular{rll}{
        1   \tab    \code{[dpqr]norm}   
            \tab    Normal distribution from R's base package, \cr
        2   \tab    \code{[dpqr]snorm}  
            \tab    Skew Normal distribution.}
}


\usage{
dsnorm(x, mean = 0, sd = 1, xi = 1.5)
psnorm(q, mean = 0, sd = 1, xi = 1.5)
qsnorm(p, mean = 0, sd = 1, xi = 1.5)
rsnorm(n, mean = 0, sd = 1, xi = 1.5)
}


\arguments{

    \item{mean, sd, xi}{
        location parameter \code{mean},
        scale parameter \code{sd},
        skewness parameter \code{xi}.
        }
    \item{n}{
        number of observations.
        } 
    \item{p}{
        a numeric vector of probabilities.
        }
    \item{x, q}{
        a numeric vector of quantiles.
        }
    
}


\value{

    All values are numeric vectors:
    \code{d*} returns the density,
    \code{p*} returns the distribution function,
    \code{q*} returns the quantile function, and
    \code{r*} generates random deviates.
    
}

\details{
  
    \bold{Symmetric Normal Distibution:}
    \cr\cr
    The functions for the normal distribution are part of R's
    base package.
    \cr
    
    \bold{Skew Normal Distribution:}
    \cr\cr
    The skew normal distribution functions are defined as described
    by Fernandez and Steel (2000).
    
}


\examples{
## snorm -
   xmpSeries("\nStart: Skew Normal Distribuion:  > ")
   par(mfrow = c(2, 2), cex = 0.75)
   set.seed(1953)
   r = rsnorm(n = 1000, mean = 1, sd = 0.5, xi = 1.5)
   plot(r, type = "l", main = "snorm: xi = 1.5")
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue4")
   x = seq(-4, 6, 0.1)
   lines(x, dsnorm(x = x, mean = 1, sd = 0.5, xi = 1.5))
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue4")
   lines(x, psnorm(x, mean = 1, sd = 0.5, xi = 1.5))
   # Compute quantiles:
   qsnorm(psnorm(q = -4:6, mean = 1, sd = 0.5, xi = 1.5), 
     mean = 1, sd = 0.5, xi = 1.5) 
}


\author{

    Diethelm Wuertz for this \R-port.
    
}


\references{
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages. 
}


\seealso{
  
    \code{sstdDistribution},
    \code{sgedDistribution}.
    
}


\keyword{distribution}

