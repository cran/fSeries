
# 	SETTINGS *******************************************************************
	library(garch)
	require(tseries)
	DEBUG = FALSE
	data(dem2gbp)
	x <- dem2gbp[,1]
		
		
# 	GARCH(2,1) *****************************************************************

	f1 = garchFit(x ~ garch(2, 1), trace=TRUE)
	f2 = garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ garch(2, 1), 
  	  	 series = x, trace = TRUE)
  	f3 = garch(x+0.005, c(2,1))
	
	# f1 - Coefficient(s):
    #               Estimate    Std.Error   t value   Pr(>|t|)    
	# omega          0.01125     0.003184     3.533   0.000412 ***
	# alpha1         0.1686      0.028664     5.882   4.05e-09 ***
	# beta1          0.4899      0.131274     3.732   0.000190 ***
	# beta2          0.2973      0.126126     2.357   0.018417 *  
	# mu            -0.00506     0.008524    -0.594   0.552638    

 	# f2 - Coefficient(s):	
  	#  	          Coefficient   Std.Error   t-value   t-prob
	# Cst(V)         0.01125     0.001543     7.293   0.0000
	# ARCH(Alpha1)   0.1686      0.016622    10.14    0.0000
	# GARCH(Beta1)   0.4899      0.11170      4.385   0.0000
	# GARCH(Beta2)   0.2973      0.10227      2.907   0.0037
	# Cst(M)        -0.00504    0.008485     -0.593   0.5530

	# f3 - Coefficient(s):
    # 		        Estimate   Std.Error    t value   Pr(>|t|)    
	# a0             0.01115    0.001487      7.498   6.46e-14 ***
	# a1             0.1681     0.016413     10.242    < 2e-16 ***
	# b1             0.4913     0.111972      4.388   1.15e-05 ***
	# b2             0.2969     0.102242      2.904   0.00369  ** 

	
# 	GARCH(1,2) *****************************************************************

	f1 = garchFit(x ~ garch(1, 2), trace=TRUE)
	f2 = garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ garch(1, 2), 
  	  	 series = x, trace = TRUE)
  	f3 = garch(x+0.0065, c(1,2)) 
	
	# f1 - Coefficient(s):
    #               Estimate    Std.Error   t value    Pr(>|t|)    
	# omega   	    0.01078     0.002190     4.925    8.44e-07 ***
	# alpha1        0.1534      0.03454      4.441    8.97e-06 ***
	# alpha2        0.0000      0.04829      0.000    1.00   
	# beta1         0.8058      0.03610     22.323     < 2e-16 ***
	# mu           -0.00628     0.008614    -0.729    0.466    
  

 	# f2 - Coefficient(s):	
  	#  	          Coefficient   Std.Error   t-value    t-prob
	# Cst(V)        0.01076     0.001546      6.957    0.0000
	# ARCH(1)       0.1534      0.01740       8.819    0.0000
	# GARCH(1)      0.0000      0.02429       0.000    1.0000
	# GARCH(2)      0.8059      0.02240      35.977    0.0000
	# Cst(M)       -0.00622     0.008436     -0.737    0.4613


	# f3 - Coefficient(s):
    # 		        Estimate   Std.Error    t value   Pr(>|t|)    
	# a0 			0.01435     0.001939      7.403   1.33e-13 ***
	# a1 			0.1774      0.01935       9.167    < 2e-16 ***
	# a2 			0.0000      0.02697       0.000   1.00    
	# b1 			0.7658      0.02642      28.990    < 2e-16 ***


# ******************************************************************************

