\name{garchOxModelling}

\alias{garchOxFit}
\alias{print.garchOx}
\alias{summary.garchOx}
\alias{plot.garchOx}


\title{R Interface for Garch Ox}


\description{

	A collection and description of functions to 
	to fit the parameters of an univariate time 
	series to GARCH models interfacing the GARCH
	OX Package. 
	\cr
	
	The family of GARCH time series models includes the following 
	processes:
	
	\tabular{rll}{
	1	\tab garch 	\tab  generalized AR conditional heteroskedastic models, \cr
	2	\tab egarch \tab  exponential GARCH models, \cr
	3	\tab aparch \tab  asymmetretic power ARCH models. }	
	
}


\usage{
garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ garch(1, 1), 
	series = x, cond.dist = c("gaussian", "t", "ged", "skewed-t"), 
	include.mean = TRUE, truncation = 100, trace = TRUE, title = NULL,
	description = NULL))
	
\method{print}{garchOx}(x, digits, \dots)
\method{summary}{garchOx}(object, \dots)
\method{plot}{garchOx}(x, \dots)
}


\arguments{

  	\item{cond.dist}{
	    a character string describing the distribution of innovations. 
	    By default the optimization is based on gaussian log likelihood 
	    parameter optimization denoted by "gaussian". Alternatively, a 
	    Student-t "t", a generalized error "sged", or a skewed
	    Student-t "skewed-t" can be chosen.
	    }	
	\item{description}{
        a character string which allows for a brief description.
        }
    \item{digits}{
  		the number of digits to be printed.
  		}
  	\item{formula.mean}{
  		formula object which specifies the mean. AR and MA models can
  		be specified as ARMA(0,q) and ARMA(p,0) respectively.
  		}
  	\item{formula.var}{
  		formula object which specifies the variance. Use \code{~garch(p,q)}
  		to specify GARCH(p,q) models, ...
  		}
  	\item{include.mean}{
  		should the mean be included? By default TRUE.
  		}	
  	\item{object}{
  		an object of class \code{garchOx} as returned from the function
  		code{garchOxFit}.
  		}
  	\item{series}{
  		the series to be modeled.
  		}
  	\item{title}{
        a character string which allows for a project title.
        }
    \item{trace}{
		a logical. Trace optimizer output? By default TRUE.
		}
  	\item{truncation}{
  		number of truncation points.
  		}	
  	\item{x}{
  		an object of class \code{garchOx} as returned from the function
  		code{garchOxFit}.
  		}
	\item{\dots}{
        additional arguments to be passed to the \code{print}, 
        \code{summary}, and \code{plot} methods.
        }
  	
}	
	

\details{

	\bold{Ox Interface:}
	\cr	
	the function \code{garchOxFit} interfaces a subset of the functionality 
	of the G@ARCH 3.0 Package written in Ox. 
	G@RCH 3.0 is one of the most sophisticated packages for modelling 
	univariate GARCH processes including GARCH, EGARCH, GJR, APARCH, 
	IGARCH, FIGARCH, FIEGARCH, FIAPARCH and HYGARCH models. Parameters
	can be estimated by approximate (Quasi-) maximum likelihood methods
	under four assumptions: normal, Student-t, GED or skewed Student-t 
	errors.
	\cr
	
	\bold{About Ox:}
	\cr
	"Ox" is an object-oriented matrix language with a comprehensive 
	mathematical and statistical function library. Many packages were 
	written for Ox including software mainly for econometric modelling. 
	The Ox packages for time series analysis and forecasting, Arfima,
	Garch and State Space Modelling are especially worth to note. 
	Since most of the R-users wan't to change to another Statistical 
	Computing environment, we made selected parts of the G@RCH OX software 
	available for them through an R-Interface. What you have to do, is 
	to check the "Ox citation and copyright" rules and if you agree and 
	fullfill the conditions, then download the OxConsole Software 
	together with the "OxGarch" Package, currently G@RCH 3.0. If you are 
	not qualified for a free license, order your copy from Timberlake 
	Consultants. We recommend to install the "Setup.exe" under the path 
	"C:\\Ox\\" and to unzip the OxGarch Package in the directory 
	"C:\\Ox\\Packages".
	\cr
	
	\bold{Ox Citation and Copyright Rules:}
	\cr
	Ox and all its components are copyright of Jurgen A. Doornik. The 
	Console (command line) versions may be used freely for academic 
	research and teaching purposes only. Commercial users and others 
	who do not qualify for the free version must purchase the Windows 
	version of Ox and GiveWin with documentation, regardless of which 
	version they use (so even when only using Ox on Linux or Unix). 
	Ox should be cited whenever it is used. Refer to the two references 
	given below. Note, failure to cite the use of Ox in published work 
	may result in loss of the right to use the free version, and an 
	invoice at the full commercial price. Ox is available from Timberlake 
	Consultants. The Ox syntax is public, and you may do with your own 
	Ox code whatever you wish.
	\cr
	
	\bold{Work to do:}
	\cr
	Note, only a small part of the functionalities are interfaced until
	now to R. But, principally it would be possible to interface also other
	functionalities offered by the Ox Garch Package.
	
}


\references{

Doornik, J.A. (2002), 
	Object-Oriented Matrix Programming Using Ox, 
	London, 3rd ed.: Timberlake Consultants Press and Oxford: 
	www.nuff.ox.ac.uk/Users/Doornik. 

Doornik, J.A. and Ooms, M. (1999), 
	A Package for Estimating, Forecasting and Simulating Arfima Models, 
	Oxford: www.nuff.ox.ac.uk/Users/Doornik. 
	
}


\author{

	Jurgen Doormik for the Ox Environment, \cr 
	NN for the Ox Garch package, and \cr
    Diethelm Wuertz for \R's Ox Garch interface. 
    
}


\seealso{
   
	\code{garchModelling}, and
   	\code{garch} from R's tseries package.
   	
}


\keyword{models}
