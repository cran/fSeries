
#
# Example: 
#	BDS Time Series Test
#
# Description:
#	This example shows how to perform the BDS time series test.
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# Test:	
	

	# BDS Test: iid Time Series:	
	par(mfrow = c(3, 1))
	x = rnorm(100)
	ts.plot(x, type = "l", main = "iid Time Series")
	bdsTest(x, m = 3)
	###
	
	
	# BDS Test: Non Identically Distributed Time Series:
	x = c(rnorm(50), runif(50))
	ts.plot(x, type = "l", main = "Non-iid Time Series")
	bdsTest(x, m = 3)  
	###
	
	
	# BDS Test: Non Independent Innovations from Quadratic Map:
	x = rep(0.2, 100)
	for (i in 2:100) x[i] = 4*(1-x[i-1])*x[i-1]
	ts.plot(x, type="l", main="Quadratic Map")
	bdsTest(x, m = 3)
	###


################################################################################

