\name{SkewGedDistribution}

\alias{SkewGedDistribution}

\alias{dged}
\alias{pged}
\alias{qged}
\alias{rged}

\alias{dsged}
\alias{psged}
\alias{qsged}
\alias{rsged}

\title{Symmetric and Skew Generalized Error Distribution}


\description{
    
    A collection and description of functions to compute 
    density, distribution function, quantile function and 
    to generate random variates for the symmetric and skew 
    generalized error distribution. 
    \cr
    
    The functions are:
    
    \tabular{ll}{
 	\code{[dpqr]ged} \tab Symmetric GED distribution, \cr
	\code{[dpqr]sged} \tab Skew GED distribution. }
}


\usage{
dged(x, mean = 0, sd = 1, nu = 2)
pged(q, mean = 0, sd = 1, nu = 2)
qged(p, mean = 0, sd = 1, nu = 2)
rged(n, mean = 0, sd = 1, nu = 2)

dsged(x, mean = 0, sd = 1, nu = 2, xi = 1.5)
psged(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
qsged(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
rsged(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
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
  
    \bold{Symmetric GED Distibution:}
    \cr\cr
    The generalized error distribution functions are defined as 
    described by Nelson (1991).
    \cr
    
    \bold{Skewed GED Distribution:}
    \cr\cr
    The skew generalized error distribution functions are defined as 
    described by Fernandez and Steel (2000).
    
}


\examples{
## sged -
   xmpSeries("\nStart: Skew Generalized Error Distribuion:  > ")
   par(mfrow = c(2, 2), cex = 0.75)
   set.seed(1953)
   r = rsged(n = 1000, mean = 1, sd = 0.5, xi = 1.5)
   plot(r, type = "l", main = "sged: xi = 1.5")
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue4")
   x = seq(-1, 5, 0.1)
   lines(x, dsged(x = x, mean = 1, sd = 0.5, xi = 1.5))
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue4")
   lines(x, psged(x, mean = 1, sd = 0.5, xi = 1.5))
   # Compute quantiles:
   qsged(psged(q = -1:5, mean = 1, sd = 0.5, xi = 1.5), 
     mean = 1, sd = 0.5, xi = 1.5) 
}


\author{

    Diethelm Wuertz for the Rmetrics \R-port.
    
}


\references{

Nelson D.B. (1991);
    \emph{Conditional Heteroscedasticity in Asset Returns: A New Approach},
    Econometrica, 59, 347--370.
    
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages. 
}


\seealso{
  
    \code{sstdDisytribution},
    \code{snormDistribution}.
    
}


\keyword{distribution}

