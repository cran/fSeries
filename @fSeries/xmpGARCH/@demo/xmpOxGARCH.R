

# garchOxFit Examples

# ------------------------------------------------------------------------------

# 	Load Benchmark Data Set:
  	data(dem2gbp)
   	x = dem2gbp[,1]

   	   	
# ------------------------------------------------------------------------------
   	
   	
# 	Example: ARMA/GARCH Models -  Gaussian Distribution
	
   	arch2 <- garchOxFit(
   		formula.mean = ~arma(0,0), formula.var = ~garch(0,2))
	 	
	garch20 <- garchOxFit(
   		formula.mean = ~arma(0,0), formula.var = ~garch(2,0))
	
	garch11 <- garchOxFit(
   		formula.mean = ~arma(0,0), formula.var = ~garch(1,1))
	 	
	ar1.garch11 <- garchOxFit(
   		formula.mean = ~arma(1,0), formula.var = ~garch(1,1))	 	
   		
   	ma1.garch11 <- garchOxFit(
   		formula.mean = ~arma(0,1), formula.var = ~garch(1,1))	 
   		
   	arma11.garch11 <- garchOxFit(
   		formula.mean = ~arma(1,1), formula.var = ~garch(1,1))
  
   		
# ------------------------------------------------------------------------------


# 	Example 2:	Other than Gaussian Distributions:
   	
	garch11.t <- garchOxFit(formula.var = ~garch(1,1), cond.dist = "t")  	
	  	
	garch11.ged <- garchOxFit(formula.var = ~garch(1,1), cond.dist = "ged")  			
   			
   	garch11 <- garchOxFit(formula.var = ~garch(1,1), cond.dist = "skewed-t")  

     
# ------------------------------------------------------------------------------


#  	Example 2: GARCH Models:

	egarch11 <- garchOxFit(formula.var = ~egarch(1,1))  	
   	
   	gjr11 <- garchOxFit(formula.var = ~gjr(1,1))  	
   	
   	aparch11 <- garchOxFit(formula.var = ~aparch(1,1))  
   	
   	igarch11 <- garchOxFit(formula.var = ~igarch(1,1)) # fails
	
# ------------------------------------------------------------------------------

	