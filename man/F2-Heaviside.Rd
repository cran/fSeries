\name{Heaviside}

\alias{Heaviside}

\alias{H}
\alias{Sign}
\alias{Delta}
\alias{Boxcar}
\alias{Ramp}


\title{Haviside and Related Functions}


\description{

	This is a collection of functions which compute the Haeviside
	and related functions.
	
	The functions are:
	
	\tabular{rll}{
		1 	\tab 	\code{H}  	
			\tab  	Computes Heaviside unit step function, \cr
		2	\tab 	\code{Sign} 
			\tab  	Just another signum function \cr
		3 	\tab 	\code{Delta}  	
			\tab  	Computes delta function, \cr
		4	\tab 	\code{Boxcar} 
			\tab  	Computes boxcar function, \cr
		5	\tab 	\code{Ramp} 
			\tab  	Computes ramp function. }
			
}


\usage{
H(x, a = 0)
Sign(x, a = 0)
Delta(x, a = 0)
Boxcar(x, a = 0.5)
Ramp(x, a = 0)
}


\arguments{
 	 
  	\item{a}{
		a numeric value, the location of the break.\cr
		}
  	\item{x}{
		a numeric vector.
		}
}


\details{
  
	The Heaviside step function \code{H} is 1 for \code{x>a}, 
	\code{1/2} for \code{x=a}, and \code{0} for \code{x<a}.
	\cr
	The Sign function \code{Sign} is \code{1} for \code{x>a}, 
	\code{0} for \code{x=a}, and \code{-1} for \code{x<a}.
	\cr
	The delta function \code{Delta} is defined as: 
	\code{Delta(x) = d/dx H(x-a)}.
	\cr
	The boxcar function \code{boxcar} is defined as: 
	\code{Boxcar(x) = H(x+a) - H(x-a)}.
	\cr
	The ramp function is defined as: 
	\code{Ramp(x)= (x-a)*H(x-a)}.
	
}


\value{

	returns the function values of the selected function.
	
}


\references{

Weisstein W. (2004);
	\emph{http://mathworld.wolfram.com/HeavisideStepFunction.html}.
	Mathworld.
	
}


\author{

    Diethelm Wuertz for this \R-port.
    
}

\examples{
## H -
   x = sort(round(c(-1, -0.5, 0, 0.5, 1, 5*rnorm(5)), 2))
   h = H(x)
	
## Sign -
   s = Sign(x)
   
## Delta -
   d = Delta(x)
   
## Boxcar -
   Pi = Boxcar(x)

## ramp - 
   r = Ramp(x)
   cbind(x = x, Step = h, Signum = s, Delta = d, Pi = Pi, R = r)		
}


\keyword{programming}

