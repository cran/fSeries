\name{UnitrootTests}

\alias{UnitrootTests}

\alias{fURTEST}
\alias{fURTEST-class}

\alias{urTest}

\alias{adfTest}
\alias{unitrootTest}

\alias{tsadfTest}
\alias{tskpssTest}
\alias{tspoTest}
\alias{tsppTest}

\alias{urersTest}
\alias{urkpssTest}
\alias{urppTest}
\alias{urspTest}
\alias{urzaTest}


\alias{print.fURTEST}
\alias{summary.fURTEST}


\title{Unit Root and Cointegration Time Series Tests}


\description{

	A collection and description of functions for unit root testing.
	The family of tests includes ADF tests based on Banerjee's et al.
	tables and on J.G. McKinnons' numerical distribution functions.
	In addition we have included functions from the contributed
	R package \code{tseries} written by Adrian Trapletti and from
	the package \code{urca} written by Bernhard Pfaff.
	\cr
	
	The functions are:
	
	\tabular{ll}{
	\code{adfTest} \tab Augmented Dickey--Fuller test for unit roots, \cr
	\code{unitrootTest} \tab the same based on McKinnons's test statistics. }
	
	Added functions based on the 'tseries' package include:
	
	\tabular{ll}{
	\code{tsadfTest} \tab Augmented Dickey--Fuller test for unit roots, \cr
	\code{tskpssTest} \tab KPSS test for level or trend stationarity, \cr
	\code{tsppTest} \tab Phillips--Perron test for unit roots, \cr	
	\code{tspoTest} \tab Phillips--Ouliaris test for cointegration. }	
	
	Added functions based on the 'urca' package include:
	
	\tabular{ll}{
	\code{urersTest} \tab Elliott--Rothenberg--Stock test for unit roots, \cr
	\code{urkpssTest} \tab KPSS unit root test for stationarity, \cr
	\code{urppTest} \tab Phillips--Perron test for unit roots, \cr	
	\code{urspTest} \tab Schmidt--Phillips test for unit roots, \cr	
	\code{urzaTest} \tab Zivot--Andrews test for unit roots. }
		
	
	NOTE: The contributed R packages \code{tseries} and \code{urca} are
	not required, the functions are builtin. 
	
}


\note{

	The output of the various hypothesis tests is an object of class
	\code{fURTEST}. The associated \code{print} method gives an unique 
	report about the test results for all tests.
	
}


\usage{
urTest(x, method = c("unitroot", "adf", "tsadf", "tskpss", "tspp", 
	"urers", "urkpss", "urpp", "ursp", "urza"), title = NULL, 
	description = NULL, \dots)

adfTest(x, type = c("nc", "c", "ct"), lags = 1)
unitrootTest(x, trend = c("nc", "c", "ct"), statistic = c("t", "n"), 
	method = "adf", lags = 1)

tsadfTest(x, alternative = c("stationary", "explosive"), 
	k = trunc((length(x)-1)^(1/3)))
tskpssTest(x, nullhyp = c("level", "trend"), lshort = TRUE) 
tsppTest(x, alternative = c("stationary", "explosive"), 
	type = c("Z(alpha)", "Z(t_alpha)"), lshort = TRUE)
tspoTest(x, demean = TRUE, lshort = TRUE)

urersTest(x, type = c("DF-GLS", "P-test"), model = c("constant", "trend"),
	lag.max = 4)
urkpssTest(x, type = c("mu", "tau"), lags = c("short", "long", "nil"),
    use.lag = NULL)
urppTest(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"),
	lags = c("short", "long"))
urspTest(x, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4),
	signif = c(0.01, 0.05, 0.1))
urzaTest(x, model = c("intercept", "trend", "both"), lag)

\method{print}{fURTEST}(x, \dots)
\method{summary}{fURTEST}(object, \dots)
}


\arguments{
  
  	\item{alternative}{
  		[tsadfTest] - \cr
		the argument indicates an alternative hypothesis. The choice 
		must be one of the elements specified in the input vector, by 
		default the first entry is selected. \cr
		[tsppTest] - \cr
		indicates the alternative hypothesis and must be one of 
		\code{"stationary"}, the default, or \code{"explosive"}. 
		One can specify just the initial letter.
		}
  	\item{demean}{
  	    [tspoTest] - \cr
		a logical indicating whether an intercept is included in the 
		cointegration regression or not.
		}
  	\item{description}{
  		a character string which allows for a brief description.
  		}
	\item{k}{
  	    [tsadfTest] - \cr
		the lag order to calculate the test statistic.
		}
  	\item{lag.max}{
		[urersTest] - \cr
		the maximum numbers of lags used for testing of a decent lag 
		truncation for the \code{"P-test"}, BIC used, or the maximum 
		number of lagged differences to be included in the test 
		regression for \code{"DF-GLS"}.
  		}
  	\item{lag}{
		[urzaTest] - \cr
		the highest number of lagged endogenous differenced variables 
		to be included in the test regression.
  		}
  	\item{lags}{
  		[urkpssTest][urppTest] - \cr
  		the maximum number of lags used for error term correction.
  		}
  	\item{lshort}{
  	    [tsppTest] - \cr
		a logical indicating whether the short or long version of the 
		truncation lag parameter is used. \cr
		[tspoTest] - \cr
		a logical indicating whether the short or long version of the 
		truncation lag parameter is used.
		}
	\item{method}{
		[urTest] - \cr
  		a character string describing the desired method, one of:
  		\code{"unitroot"}, \code{"adf"}, \code{"tsadf"}, \code{"tskpss"}, 
  		\code{"tspp"}, \code{"urers"}, \code{"urkpss"}, \code{"urpp"}, 
  		\code{"ursp"}, \code{"urza"}, or \code{"W3SLS"}. \cr
  		[unitrootTest] - \cr
  		a character string specifying the method of the test. Valid 
  		choices are \code{"adf"} for the augmented Dickey-Fuller test, 
  		and \code{"pp"} for Phillips-Perron test. The default is 
  		\code{"adf"}. 
  		}
	\item{model}{
		[urersTest] - \cr
		a character string dennoting the deterministic model used for 
		detrending, either \code{"constant"}, the default, or 
		\code{"trend"}. \cr
		[urppTest] - \cr
		a character string which determines the deterministic part in 
		the test regression, either \code{"constant"}, the default, or 
		\code{"trend"}. \cr
		[urzaTest] - \cr
		a character string specifying if the potential break occured 
		in either the \code{"intercept"}, the linear \code{"trend"} or 
		in \code{"both"}.
  		}
  	\item{nullhyp}{
	    [tskpssTest] - \cr
		a character string which indicates the null hypothesis and must 
		be one of  \code{"level"} (default) or \code{"trend"}. One can 
		specify just the initial letter.
		}
  	\item{object}{
		an object of class \code{"fURTEST"}.
  		}
  	\item{pol.deg}{
		[urspTest] - \cr
		the polynomial degree in the test regression.
  		}
  	\item{signif}{
		[urspTest] - \cr
		the significance level for the critical value of the test 
		statistic.
  		}
  	\item{statistic}{
		[unitrootTest] - \cr
		a character string specifying the test statistic. Valid choices 
		are  \code{"t"} for the t-statistic, the default, and \code{"n"} 
		for the normalized statistic, sometimes also referred to as 
		the rho-statistic. 
  		}
  	\item{title}{
  		a character string which allows for a project title.
  		}
	\item{trend}{
		[unitrootTest] - \cr
		a character string describing the type of the unit root 
		regression. Valid choices are \code{"nc"} for a regression 
		with no intercept (constant) nor time trend, and \code{"c"} 
		for a regression with an intercept (constant) but no time 
		trend, \code{"ct"} for a regression with an intercept 
		(constant) and a time trend. The default is \code{"c"}. 
  		}
  	\item{type}{
  	    [tsppTest] - \cr
    	indicates which variant of a test will be computed. The choice 
    	must be one of the elements specified in the input vector, by 
    	default the first entry is selected. \cr
    	[urkpssTest] - \cr
    	a character string which denotes the type of deterministic part,
    	either \code{"mu"}, the default, or \code{"tau"}.
    	[urppTest] - \cr
    	a character string which specifies the test type, either 
    	\code{"Z-alpha"}, the default, or \code{"Z-tau"}.
    	[urspTest] - \cr
    	a character string which specifies the test type, either 
    	\code{"tau"}, the default, or \code{"rho"}.
    	} 
    \item{use.lag}{
		[urkpssTest] - \cr
		a character string specifying the number of lags. Allowed
		arguments are \code{lags=c("short", "long", "nil")}, for more 
		information see the details section.
  		}
  	\item{x}{
		a numeric vector or time series object.
		}
	\item{\dots}{
		[urTest] - \cr
		optional arguments passed to the underlying test functions.
		} 	
		
}


\details{

	\bold{ADF Tests:}
	
	The \code{adftest} computes test statistics and p values using
	the implementation from Trapletti's augmented Dickey--Fuller test 
	for unit roots, see below. In contrast to the function 
	\code{tsadfTest} three kind of test types can be selected.
	\cr
	
	
	\bold{Unit Root Tests from Adrian Trapletti's "tseries" Package:}
	\cr
	
	\emph{Augmented Dickey--Fuller Test for Unit Roots:} \cr
	The general regression equation which incorporates a constant and a
	linear trend is used and the t-statistic for a first order
	autoregressive coefficient equals one is computed. The number of lags
	used in the regression is \code{k}. The default value of
	\code{trunc((length(x)-1)^(1/3))} corresponds to the suggested upper
	bound on the rate at which the number of lags, \code{k}, should be
	made to grow with the sample size for the general \code{ARMA(p,q)}
	setup. Note that for \code{k} equals zero the standard Dickey-Fuller
	test is computed. The p-values are interpolated from Table 4.2, p. 103
	of Banerjee et al. (1993). If the computed statistic is outside the
	table of critical values, then a warning message is generated. 
	Note, that missing values are not allowed. \cr
	\code{[tseries:adf.test]} \cr
		
	\emph{Kwiatkowski--Phillips--Schmidt--Shin Test:} 
	\cr
	To estimate \code{sigma^2} the Newey-West estimator is used.
  	If \code{lshort} is \code{TRUE}, then the truncation lag parameter is
  	set to \code{trunc(3*sqrt(n)/13)}, otherwise
  	\code{trunc(10*sqrt(n)/14)} is used.  The p-values are interpolated
  	from Table 1 of Kwiatkowski et al. (1992). If the computed statistic
  	is outside the table of critical values, then a warning message is
  	generated. \cr
  	\code{[tseries:kpss.test]} \cr
		
	\emph{Phillips--Ouliaris Test:} 
	\cr
	The \code{poTest} Z(alpha) statistic for a unit root in the
  	residuals of the cointegration regression is computed, see also
  	\code{ppTest}. The unit root is estimated from a regression of
  	the first variable (column) of \code{x} on the remaining variables of
  	\code{x} without a constant and a linear trend.  To estimate
  	\code{sigma^2} the Newey-West estimator is used.  If \code{lshort} is
  	\code{TRUE}, then the truncation lag parameter is set to
  	\code{trunc(n/100)}, otherwise \code{trunc(n/30)} is used.  The
  	p-values are interpolated from Table Ia and Ib, p. 189 of Phillips and
  	Ouliaris (1990).  If the computed statistic is outside the table of
  	critical values, then a warning message is generated. The dimension 
  	of \code{x} is restricted to six variables. \cr
  	\code{[tseries:po.test]} \cr  
  	
	\emph{Phillips--Perron Test:} 
	\cr	
	The general regression equation which incorporates a constant and a
  	linear trend is used and the \code{Z(alpha)} or \code{Z(t_alpha)}
  	statistic for a first order autoregressive coefficient equals one are
  	computed.  To estimate \code{sigma^2} the Newey-West estimator is
  	used.  If \code{lshort} is \code{TRUE}, then the truncation lag
  	parameter is set to \code{trunc(4*(n/100)^0.25)}, otherwise
  	\code{trunc(12*(n/100)^0.25)} is used.  The p-values are interpolated
  	from Table 4.1 and 4.2, p. 103 of Banerjee et al. (1993). If the
  	computed statistic is outside the table of critical values, then a
  	warning message is generated. \cr
  	\code{[tseries:pp.test]} \cr
  	

	\bold{Unit Root Tests from Berhard Pfaff's "urca" Package:}
	\cr

	\emph{Elliott--Rothenberg--Stock Test for Unit Roots:} 
	\cr
	To improve the power of the unit root test, Elliot, Rothenberg and 
	Stock proposed a local to unity detrending of the time series. ERS 
	developed a feasible point optimal test, \code{"P-test"}, which 
	takes serial correlation of the error term into account. The second 
	test type is the \code{"DF-GLS"} test, which is an ADF-type test 
	applied to the detrended data without intercept. Critical values 
	for this test are taken from MacKinnon in case of \code{model="constant"}
	and else from Table 1 of Elliot, Rothenberg and Stock. \cr
	\code{[urca:ur.ers]} \cr
		
	\emph{KPSS Test for Unit Roots:} 
	\cr
	Performs the KPSS unit root test, where the Null hypothesis is 
	stationarity. The test types specify as deterministic component 
	either a constant \code{"mu"} or a constant with linear trend 
	\code{"tau"}. \code{lags="short"} sets the number of lags to 
	\emph{root 4 of [4 times (n/100)}, whereas \code{lags="long"} 
	sets the number of lags to \emph{root 4 of [12 times (n/100)]}. 
	If \code{lags="nil"} is choosen, then no error correction is made. 
	Furthermore, one can specify a different number of maximum lags 
	by setting use.lag accordingly. \cr
	\code{[urca:ur.kpss]} \cr
	
	\emph{Phillips--Perron Test for Unit Roots:} 
	\cr
	Performs the Phillips and Perron unit root test. Beside the 
	Z statistics Z-alpha and Z-tau, the Z statistics for the 
	deterministic part of the test regression are computed, too. 
	For correction of the error term a Bartlett window is used. \cr
	\code{[urca:ur.pp]} \cr
		
	\emph{Schmidt--Phillips Test for Unit Roots:} 
	\cr
	Performs the Schmidt and Phillips unit root test, where under 
	the Null and Alternative Hypothesis the coefficients of the 
	deterministic variables are included. Two test types are available: 
	the \code{"rho-test"} and the \code{"tau-test"}. Both tests are 
	extracted from the LM principle. \cr
	\code{[urca:ur.sp]} \cr
		
	\emph{Zivot--Andrews Test for Unit Roots:} 
	\cr
	Performs the Zivot and Andrews unit root test, which allows a 
	break at an unknown point in either the intercept, the linear 
	trend or in both. This test is based upon the recursive estimation 
	of a test regression. The test statistic is defined as the 
	minimum t-statistic of the coeffcient of the lagged endogenous 
	variable. \cr
	\code{[urca:ur.za]} 
 
}


\value{

	All tests return an object of class \code{"fURTEST"} with the
	following slots:


	\item{@call}{
		the function call.		
		}
	\item{@data}{
		a data frame with the input data.
		}
	\item{@data.name}{
		a character string giving the name of the data.
		}
	\item{@test}{
		a list object which holds the output of the underlying
		test function, usually the same entries as an object of
		class \code{"htest"}.
		}
	\item{@title}{
		a character string with the name of the test.
		}
	\item{@description}{
		a character string with a brief description of the
		test.
		}
		
	The entries of the \code{@test} slot include the following components:

  	\item{$statistic}{
		the value of the test statistic.
		}
  	\item{$parameter}{
		the lag order.
		}
  	\item{$p.value}{
		the p-value of the test.
		}
  	\item{$method}{
		a character string indicating what type of test was
		performed.
		} 
  	\item{$data.name}{
		a character string giving the name of the data.
		}
  	\item{$alternative}{
		a character string describing the alternative
		hypothesis.
		}
		
}


\references{

Banerjee A., Dolado J.J., Galbraith J.W., Hendry D.F. (1993);
	\emph{Cointegration, Error Correction, and the Econometric 
		Analysis of Non-Stationary Data},
	Oxford University Press, Oxford. 
	
Dickey, D.A., Fuller, W.A. (1979);
	\emph{Distribution of the estimators for autoregressive time 
		series with a unit root}, 
	Journal of the American Statistical Association 74, 427--431. 
	
Kwiatkowski D., Phillips P.C.B, Schmidt P., Shin Y. (1992);
	\emph{Testing the Null Hypothesis of Stationarity against 
	the Alternative of a Unit Root},
	Journal of Econometrics 54, 159--178.
	
MacKinnon, J.G. (1996);
	\emph{Numerical distribution functions for unit root and 
		cointegration tests},
	Journal of Applied Econometrics 11, 601--618.
	
Perron P. (1988);
	\emph{Trends and Random Walks in Macroeconomic Time Series},
	Journal of Economic Dynamics and Control 12, 297--332.

Phillips P.C.B., Ouliaris S. (1990);
	\emph{Asymptotic Properties of Residual Based Tests for Cointegration},
	Econometrica 58, 165--193.
	
Phillips P.C.B., Perron P. (1988);
	\emph{Testing for a unit root in time series regression}, 
	Biometrika 75, 335--346.
  	
Said S.E., Dickey D.A. (1984);
	\emph{Testing for Unit Roots in Autoregressive-Moving Average 
		Models of Unknown Order},
	Biometrika 71, 599--607.
	
Schwert G.W. (1989);
  	\emph{Tests for Unit Roots: A Monte Carlo Investigation},
    Journal of Business and Economic Statistics 2, 147--159.

	
}


\author{

	Adrian Trapletti for the tests from R's "tseries" package,\cr
	Bernhard Pfaff for the tests from R's "urca" package,\cr
	Diethelm Wuertz for the Rmetrics \R-port.
}


\seealso{

	Fur further details we refer to the manual pages of the 
	\code{"urca"} package.
	
	\code{classicalTests},
	\code{tseriesTests},
	\code{lmTests}.
	
}


\examples{
## Unit Root Tests build on Adrian Trapletti's Implementation:

## tsadfTest - 
   xmpSeries("\nStart: Augmented Dickey-Fuller Test for Unit Roots >")
   # A time series which contains no unit-root:
   x = rnorm(1000)  
   tsadfTest(x)
   # A time series which contains a unit-root:
   y = diffinv(x)
   tsadfTest(y)

## tskpssTest -
   xmpSeries("\nNext: KPSS test for Level or Trend Stationary >")
   # Time Series is Level Stationary:
   x = rnorm(1000)
   tskpssTest(x)
   # Time Series has Unit Root:
   tskpssTest(cumsum(x))
   # Time Series is Trend Stationary:
   x = 0.3*(1:1000) + rnorm(1000)
   tskpssTest(x, nullhyp = "trend")
   
## tspoTest -
   xmpSeries("\nNext: Phillips-Ouliaris Test for Cointegration >")
   # Non-Cointegrated Case:
   x = ts(diffinv(matrix(rnorm(2000), 1000, 2)))
   tspoTest(x)
   # Cointegrated Case:
   x = diffinv(rnorm(1000))
   y = 2.0 - 3.0 * x + rnorm(x, sd = 5)
   z = ts(cbind(x, y))
   tspoTest(z)
   
## tsppTest -
   xmpSeries("\nNext: Phillips-Perron Test for Unit Roots >")
   # The Time Series has no Unit Root:
   x = rnorm(1000)
   tsppTest(x)
   # The Time Series has Unit Root:
   y = cumsum(x)  
   tsppTest(y)     
 
## Unit Root Tests build on Bernhard Pfaff's Implementation:
    
## ur*Test - 
   # Examples can be found in the demo file "xmpTestUnitRoots".  
}


\keyword{htest}

