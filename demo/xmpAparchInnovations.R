
#
# Example: 
#	Garch Model Simulations with several kind of innovations
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


# Settings:
  	
	par(mfrow = c(3, 2), cex = 0.6)
  
# GARCH(1,1) Normal Innovations

	# Simulate:
  	set.seed(4711)
  	x = aparchSim(model = list(omega = 1.02e-6, 
  		alpha = 0.1, gamma = 0, alpha.lags = 1, beta = 0.6, 
  		beta.lags = 1, delta = 2), innov = rnorm(2000), 
  		start.innov = rnorm(100))	
	# Plot:
	ts.plot(x, ylim = c(-0.015, 0.015), 
		main = "GARCH(1,1) Normal Innovations")
	acf(abs(x), lag.max = 12)
	var(x)		
	# Start Next:
	Continue = readline("Press any key > ") 

# GARCH(1,1) Student t(4) Innovations

	# Simulate:
  	set.seed(1858)
  	x = aparchSim(model = list(omega = 0.45e-6, 
  		alpha = 0.1, gamma = 0, alpha.lags = 1, beta = 0.6, 
  		beta.lags = 1, delta = 2), innov = rt(2000, df = 4), 
  		start.innov = rt(100, df = 4)) 		
	# Plot:
	ts.plot(x, ylim = c(-0.015, 0.015), 
		main = "GARCH(1,1) - Student t(4) Innovations")
	acf(abs(x), lag.max = 12)
	var(x)		
	# Start Next:
	Continue = readline("Press any key > ") 	
	
# GARCH(1,1) 1.9 Stable Innovations

	# Simulate:
  	set.seed(1245)
  	x = aparchSim(model = list(omega = 0.3e-6, 
  		alpha = 0.1, gamma = 0, alpha.lags = 1, beta = 0.6, 
  		beta.lags = 1, delta = 2), innov = rsymstb(2000, 
  		alpha = 1.95), start.innov = rsymstb(100, alpha = 1.95)) 		
	# Plot:
	ts.plot(x, ylim = c(-0.015, 0.015), 
		main = "GARCH(1,1) - 1.95 Stable Innovations")
	acf(abs(x), lag.max = 12)
	var(x)
	
