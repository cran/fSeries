
#
# Example: 
#	Perorm a Augmented Dickey - Fuller Unit Root Test, a 
#	Phillips - Perron Unit Root Test and a KPSS Test for 
#	Stationarity
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------

	
# Augmented Dickey - Fuller Test:	
	
	# A time series which contains no unit root:
	x = rnorm(1000)  
	adfTest(x)
	# A time series which contains a unit-root:
	y = diffinv(x)
	adfTest(y)
	Continue = readline("Press any key > ")
	
# Phillips - Perron Unit Root Test:

	# A time series which contains no unit root:
	x = rnorm(1000)
	ppTest(x)
	# A time series which contains a unit-root:
	y = cumsum(x)  
	ppTest(y)
	Continue = readline("Press any key > ")
	
# KPSS Test for Stationarity:

	# A time series which is level stationary:
	x = rnorm(1000)
	kpssTest(x)
	# A time series which contains a unit-root:
	y = cumsum(x)
	kpssTest(y)
	# A time series which is trend stationary:
	x = 0.3*(1:1000)+rnorm(1000)
	kpssTest(x, null = "Trend")

	 