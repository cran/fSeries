
#
# Example: 
# 	Perform a Jarque - Bera Normality Test
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------

	
# Jarque-Bera Normality Test: Null
	
	x = rnorm(100)
	jbTest(x)
	
# Jarque-Bera Normality Test: Alternative
	
	x = runif(100)-0.5
	jbTest(x)

