
#
# Example: Aparch Model Simulations and Parameter Estimations
#

# Settings:
		par(mfrow=c(4,2))


# Bollerslev's t-GARCH(1,1):
# --------------------------------------------------------------------

	# Simulate:
  		set.seed(5066)
  		x <- aparchSim(model=list(omega=1e-6, 
  			alpha=0.1, gamma=0, alpha.lags=1, 
  			beta=0.8, beta.lags=1, delta=2), 
  	    	innov=rt(5000, df=10), start.innov=rt(1000, df=10))
  	    ts.plot(x, main="GARCH(1,1) - Student-t")
  	    acf(abs(x), lag.max=20)
  	    	
	# Estimate:
  	    aparchFit(x=x, 
  	    	order=list(alpha.lags = 1, beta.lags = 1, delta = 2),
  	    	opt=list(gamma=FALSE, delta=FALSE, disparm=TRUE),
			distribution="t", disparm=4,
    		doprint=FALSE)
	
    Continue <- readline("Press any key > ") 
    

# Taylor-Schwert TS-GARCH(1,5[1,5]) subset model, delta=1:
# --------------------------------------------------------------------
  		
	# Simulate:
  		set.seed(5066)
  		x <- aparchSim(model=list(omega=1e-6, 
  	    	alpha=0.1, gamma=0, alpha.lags=1, 
  	    	beta=c(0.5, 0.3), beta.lags=c(1, 5), delta=1),
  	    	innov=rnorm(5000), start.innov=rnorm(1000))
  	    ts.plot(x, main="Subset TS GARCH(1,5) - Mormal")
  	    acf(abs(x), lag.max=20)
  	    	
	# Estimate:
  	    aparchFit(x, 
  	    	order=list(alpha.lags = 1, beta.lags = c(1, 5), delta = 1),
			opt=list(gamma=FALSE, delta=FALSE, disparm=FALSE),
    		doprint=FALSE)
    		
    Continue <- readline("Press any key > ") 
 
    		
# asymmetric GJR(1,1):	    	
# --------------------------------------------------------------------
  	    	
	# Simulate:
  		set.seed(5066)
  		x <- aparchSim(model=list(omega=1e-6, 
  	    	alpha=0.1, gamma=0.2, alpha.lags=1, 
  	    	beta=0.6, beta.lags=1, delta=2), 
  	    	innov=rnorm(1000), start.innov=rnorm(100))
  	    ts.plot(x, main="GJR(1,1) - Normal")
  	    acf(abs(x), lag.max=20)
  	     	
	# Estimate:    		
	    aparchFit(x=x, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=TRUE, delta=FALSE, disparm=FALSE),
    		doprint=FALSE)

			
