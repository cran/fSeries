\name{armaStatistics}

\alias{armaStatistics}

\alias{armaRoots}
\alias{armaTrueacf}

\title{Statisitcs of the True ARMA Process}


\description{

    A collection of functions for the analysis of the true ARMA
    process. The first investigates the characteristic polynomial
    and the second the true autocorrelation function.
    \cr
    
    The functions are:
		
	\tabular{rll}{
		1	\tab 	\code{armaRoots}	 
			\tab  	Roots of the characteristic ARMA polynomial, \cr
    	2	\tab 	\code{armaTrueacf}
    		\tab  	True autocorrelation function of an ARMA process. }

}


\usage{
armaRoots(coefficients, n.plot = 400, digits = 4, \dots) 
armaTrueacf(model, lag.max = 20, type = "correlation", doplot = TRUE) 
}


\arguments{

  	\item{coefficients}{	
		a numeric vector with the coefficients of the 
		characterisitic polynomial.
		}
  	\item{digits}{	
		output precision, an integer value.
		}
	\item{doplot}{
		a logical. Should a plot be displayed?
		}
	\item{lag.max}{
		 maximum number of lags at which to calculate the acf or 
		 pacf, an integer value by default 20.
		 }
	\item{model}{
		a specification of the ARMA model with two elements:
		\code{model$ar} is the vector of the AR coefficients, and 
		\code{model$ma} is the vector of the MA coefficients.
		}
	\item{n.plot}{	
		the number of data point to plot the unit circle;
		an integer value.
		}
  	\item{type}{
		a character string, "correlation" to compute the true 
		autocorrelation function, "partial" to compute the true
		partial autocorrelation function, or "both" if both
		functions are desired. The start of one of the strings 
		will suffice.
		} 
	\item{\dots}{
		arguments past to the plot function.
		} 		
	
}


\value{
	
	\code{armaRoots} 
	\cr
	returns a three column data frame with the real, the imaginary part 
	and the radius of the roots. The number of rows corresponds
	to the coefficients.
	
	\code{armaTrueacf}
	\cr
	returns a two column data frame with the lag and the correlation 
	function. 
	
}


\author{

	B.D. Ripley for the underlying \code{ARMAacf} function, \cr
	Diethelm Wuertz for this \R-port.
	
}


\seealso{
  
	\code{\link{armaModelling}}.
	
}


\examples{
## armaRoots -
   xmpSeries("\nStart: Roots of the Characteristic Polynomial > ")
   # Calculate and plot the  of an ARMA process:
   par(mfrow = c(2, 2), cex = 0.7)
   coefficients = c(-0.5, 0.9, -0.1, -0.5)
   armaRoots(coefficients)
   
## armaTrueacf -
   xmpSeries("\nNext: True ACF of an ARMA Process > ")
   model = list(ar = c(0.3, +0.3), ma = 0.1)
   armaTrueacf(model)
   model = list(ar = c(0.3, -0.3), ma = 0.1)
   armaTrueacf(model)
}


\keyword{models}

