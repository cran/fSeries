\name{sstdDistribution}

\alias{sstdDistribution}

\alias{dstd}
\alias{pstd}
\alias{qstd}
\alias{rstd}

\alias{dsstd}
\alias{psstd}
\alias{qsstd}
\alias{rsstd}


\title{Symmetric and Skew Student-t Distribution}


\description{
    
    A collection of functions to compute density, distribution function, 
    quantile function and to generate random variates for the symmetric
    and skew Sudent-t distribution with unit variance. 
    \cr
    
    The functions are:
    
    \tabular{rll}{
        1   \tab    \code{[dpqr]std}   
            \tab    Symmetric Student-t Distribution, \cr
        2   \tab    \code{[dpqr]sstd}  
            \tab    Skew Student-t Distribution.}
}


\usage{
dstd(x, mean = 0, sd = 1, nu = 5)
pstd(q, mean = 0, sd = 1, nu = 5)
qstd(p, mean = 0, sd = 1, nu = 5)
rstd(n, mean = 0, sd = 1, nu = 5)

dsstd(x, mean = 0, sd = 1, nu = 5, xi = 1.5)
psstd(q, mean = 0, sd = 1, nu = 5, xi = 1.5)
qsstd(p, mean = 0, sd = 1, nu = 5, xi = 1.5)
rsstd(n, mean = 0, sd = 1, nu = 5, xi = 1.5)
}


\arguments{

    \item{mean, sd, nu, xi}{
        location parameter \code{mean},
        scale parameter \code{sd},
        shape parameter \code{nu},
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
  
    \bold{Symmetric Student-t Distibution:}
    \cr\cr
    The functions for the symmetric Student-t distribution are 
    rescaled in such a way that they have unit variance in 
    contrast to the Student-t family \code{dt}, \code{pt}, 
    \code{qt} and \code{rt} which are part of R's base package.
    \cr
    
    \bold{Skew Student-t Distribution:}
    \cr\cr
    The skew Student-t distribution functions are defined as described
    by Fernandez and Steel (2000). Note that the function have unit
    variance.
    
}


\examples{
## sstd -
   xmpSeries("\nStart: Skew Student-t Distribuion:  > ")
   par(mfrow = c(2, 2), cex = 0.75)
   set.seed(1953)
   r = rsstd(n = 1000, nu =4, xi = 1.5)
   # Print Variance:
   var(r)
   plot(r, type = "l", main = "sstd: xi = 1.5")
   # Plot empirical density and compare with true density:
   hist(r, n = 30, xlim = c(-5, 5), probability = TRUE, 
     border = "white", col = "steelblue4")
   x = seq(-5, 5, 0.1)
   lines(x, dsnorm(x = x, xi = 1.5))
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue4")
   lines(x, psstd(x, xi = 1.5))
   # Compute quantiles:
   qsstd(psstd(q = -5:5, xi = 1.5), xi = 1.5) 
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
  
    \code{snormDistribution},
    \code{sgedDistribution}.
    
}


\keyword{distribution}

