
#
# Example: 
#	Unit Root and Stationarit Tests
#
# Description:
#   These examples show how to perform the Augmented Dickey-Fuller unit 
#	root test, the Phillips-Perron unit root test and the KPSS test for 
#	stationarity using wrapper functions based on functions from Adrian 
#   Trapletti's 'tseries package.
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# Tests:

	
	# Augmented Dickey - Fuller Test:		
	# A time series which contains no unit root:
	x = rnorm(1000)  
	tsadfTest(x)
	# A time series which contains a unit-root:
	y = diffinv(x)
	tsadfTest(y)
	###
	
	
	# Phillips - Perron Unit Root Test:
	# A time series which contains no unit root:
	x = rnorm(1000)
	tsppTest(x)
	# A time series which contains a unit-root:
	y = cumsum(x)  
	tsppTest(y)
	###
	
		
	# KPSS Test for Stationarity:
	# A time series which is level stationary:
	x = rnorm(1000)
	tskpssTest(x)
	# A time series which contains a unit-root:
	y = cumsum(x)
	tskpssTest(y)
	# A time series which is trend stationary:
	x = 0.3*(1:1000)+rnorm(1000)
	tskpssTest(x, null = "Trend")
	###


################################################################################

	 