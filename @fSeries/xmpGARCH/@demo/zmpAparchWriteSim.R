
#
# Example: GARCH Modeling - xmpAparchWriteSim
#

# Write a R-function which simulates an APARCH Process

.aparchSim <-
function(model = list(omega=1e-6, alpha=0.1, gamma=0, alpha.lags=1,
beta=0.8, beta.lags=1, delta=1), n = 1000, innov = NULL, n.start = 100, 
start.innov = NULL, rand.gen = rnorm, ...)
{	
    # A function written by D. Wuertz	

    # This is an example for a function without the needs of
    # a compiled external function.
    
  	# Innovations:		
	if (is.null(innov)) innov <- rand.gen(n, ...)
	if (is.null(start.innov)) start.innov <- rand.gen(n.start, ...)	
  	x <- z <- c(start.innov, innov)
  	h <- 0 * x 
	 
    # Orer maxpq:
    deltainv <- 1/model$delta	
    maxpq <- max( 
    	model$alpha.lags[length(model$alpha.lags)], 
    	model$beta.lags[length(model$beta.lags)] )
    
    # Simulate:
    for (n in (maxpq+1):length(z)) {     
    	h[n] <- model$omega + 	
			sum(model$alpha*(abs(abs(x[n-model$alpha.lags]) - 
				model$gamma*x[n-model$alpha.lags])^model$delta)) +
			sum(model$beta*h[n-model$beta.lags]) 
    	x[n] <- h[n]**deltainv * z[n] }
		
  	# Return Value:
	as.ts(x[(length(start.innov)+1):length(z)])
}

