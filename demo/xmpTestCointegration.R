
#
# Example: 
#	Perform a Phillips - Ouliaris Cointegration Test
#	
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
	

	# Phillips-Ouliaris Test: Non-Cointegrated Case:
	x = ts(diffinv(matrix(rnorm(2000), 1000, 2)))
	tspoTest(x)
	###
	
	
	# Phillips - Ouliaris Test: Cointegrated Case:
	x = diffinv(rnorm(1000))
	y = 2.0 - 3.0*x + rnorm(x, sd = 5)
	z = ts(cbind(x, y))
	tspoTest(z)
	###
	
	
################################################################################

