
#
# Example:
#
#	Fitting ARCH Models
#


# 	SETTINGS *******************************************************************
	# library(finGARCH)
	# require(tseries) # for comparison with A. Trapletti's garch() 
	DEBUG = FALSE
	data(dem2gbp)
	x <- dem2gbp[, 1]
	
	
# 	ARCH(1) - EXTRACT MEAN *****************************************************


# 	ARCH(1) - INCLUDE MEAN *****************************************************	
	f1 = garchFit(x ~ arch(1), trace=TRUE)
		
	# Final Estimate:
    #	          omega       alpha1           mu 
 	# f1:	0.146553862  0.371448240 -0.001590542 


# 	ARCH(3[1,3]) - SUBSET MODEL ************************************************
	g1 = garchFit(x ~ arch(3), fixed=c(F, F, T, F, F), trace=TRUE)	
	
 	# Final Estimate:
 	# Coefficient(s):
    #         Estimate   Std.Error  t.value   Pr(>|t|)    
	# omega   0.123791    0.006364   19.451    < 2e-16 ***
	# alpha1  0.329306    0.038357    8.585    < 2e-16 ***
	# alpha2  0.000000    0.035094    0.000       1.00    
	# alpha3  0.146979    0.027406    5.363   8.19e-08 ***
	# mu     -0.001908    0.008918   -0.214       0.83    

 
# ------------------------------------------------------------------------------

