\name{fSeriesTools}


\alias{fSeriesTools}

\alias{xmpSeries}
\alias{xmpfSeries}


\title{fSeries Tools}


\description{

	Popups the example menu and starts or continues the example 
	and demo programs.
	
}


\usage{
xmpfSeries()
xmpSeries(prompt = "")
}


\arguments{

  	\item{prompt}{
		the string printed when prompting the user for input. 
		}	
}


\value{
  
	xmpfSeries
	\cr
	Popups the example menu.
	\cr
	
	xmpSeries
	\cr
	Nothing, the default, or the the prompt if you have set
	\code{xmpSeries = readline} on the command prompt. 
	
}


\details{

	The example in the manual pages may be interactive and ask for
	input from the user. To achieve this you have to type on the
	command line: \code{xmpSeries = readline}

}


\author{

	Diethelm Wuertz for this R-Port.

}


\examples{
## currently none ...
}


\keyword{programming}

