
#
# Example: 
# 	Normality Tests
#
# Description:
#	This example shows ho to perform the Jarque-Bera test for
#	normality.
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# Tests:
	

	# Jarque-Bera Normality Test: Null
	x = rnorm(100)
	jbTest(x)
	###
	
	
	# Jarque-Bera Normality Test: Alternative	
	x = runif(100)-0.5
	jbTest(x)
	###


################################################################################

