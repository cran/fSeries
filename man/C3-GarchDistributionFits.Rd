\name{GarchDistributionFits}

\alias{GarchDistributionFits}

\alias{normFit}
\alias{snormFit}
\alias{stdFit}
\alias{sstdFit}
\alias{gedFit}
\alias{sgedFit}

\title{Parameter Fit of a Distribution}


\description{

	A collection and description of maximum likelihood 
	estimators to fit the parameters of a distribution 
	and to compute basic statistical properties. Included 
	are estimators for the symmetric and skew normal, 
	the Student-t, and the generalized error distribution. 
	\cr

	The functions are:
	
	\tabular{ll}{
	\code{normFit} \tab MLE parameter fit for a Normal distribution, \cr
	\code{snormFit} \tab MLE parameter fit for a skew Normal distribution, \cr
	\code{stdFit} \tab MLE parameter fit for a Sudent-t distribution, \cr
	\code{stdFit} \tab MLE parameter fit for a skew Sudent-t distribution, \cr
	\code{gedFit} \tab MLE parameter fit for a generalized error distribution, \cr
	\code{nigFit} \tab MLE parameter fit for a skew generalized error distribution. }

}


\usage{
normFit(x, \dots)
snormFit(x, \dots)
stdFit(x, \dots)
sstdFit(x, \dots)
gedFit(x, \dots)
sgedFit(x, \dots)
}


\arguments{
  
	\item{x}{
		a numeric vector.
		}
	\item{\dots}{
		parameters parsed to the optimization function \code{nlm}.
		}

}


\value{
  
  	The functions \code{*Fit} return a list with the following components:
  
  	\item{estimate}{
		the point at which the maximum value of the log liklihood 
		function is obtained.
		}
  	\item{objective}{
		the value of the estimated maximum, i.e. the value of the
		log liklihood function.
		}
  	\item{message}{
		an integer indicating why the optimization process terminated.
		}
  	\item{code}{
		an integer indicating why the optimization process terminated.\cr
		1: relative gradient is close to zero, current iterate is probably 
			solution; \cr
		2: successive iterates within tolerance, current iterate is probably 
			solution; \cr
		3: last global step failed to locate a point lower than \code{estimate}. 
			Either \code{estimate} is an approximate local minimum of the 
			function or \code{steptol} is too small; \cr
 		4: iteration limit exceeded; \cr
 		5: maximum step size \code{stepmax} exceeded five consecutive times. 
 			Either the function is unbounded below, becomes asymptotic to a 
 			finite value from above in some direction or \code{stepmax} 
			is too small.
			}
  	\item{gradient}{
		the gradient at the estimated maximum.
		}
  	\item{steps}{
		number of function calls.
		}
  
}


\details{

	The function \code{\link{nlm}} is used to minimize the "negative" maximum 
	log-likelihood function. \code{nlm} carries out a minimization using 
	a Newton-type algorithm.
	
}


\examples{
## SOURCE("fBasics.A0-SPlusCompatibility")
## SOURCE("fSeries.D1-HeavisideFunction")
## SOURCE("fSeries.D2-SkewNormalDistribution")
## SOURCE("fSeries.D3-SkewStudentDistribution")
## SOURCE("fSeries.D4-SkewGedDistribution")
## SOURCE("fSeries.D5-GarchDistributionFits")

## snormFit -
   xmpSeries("\nStart: MLE Fit to skew Normal Density > ")
   options(warn = -1) # suppress negative logs from nlm
   normFit(rnorm(1000))
   snormFit(rnorm(1000))
   
## sstdFit -
   xmpSeries("\nNext: MLE Fit to skew Student-t Density > ")
   sstdFit(rsstd(1000, mean = 1, sd = 1.5, nu = 7, xi = 2/3))
   
## sgedFit -
   xmpSeries("\nNext: MLE Fit to skew Generalized Error Density > ")
   sgedFit(rsged(1000, mean = -1, sd = 0.5, nu = 3, xi = 3/2),
   	 print.level = 2)
}


\seealso{
  
	\code{\link{dsnorm}},
  	\code{\link{dsstd}}, 
  	\code{\link{dsged}}, 
  	\code{\link{nlm}}.
  	
}


\author{

	Diethelm Wuertz for the Rmetrics \R-port.
	
}


\keyword{distribution}

