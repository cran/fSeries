
#
# Title: 
#	Simulate Artificial GARCH Processes
#   
# Description: 
#	This example simulates GARCH processes of the following type: 
#


# Settings:

	par(mfrow = c(2, 1))
	
# GARCH: GARCH(1,1):

	# Bollerslev's GARCH(1,1) Process:
   	model = list(omega = 1e-6, alpha = 0.1, beta = 0.8)
   	spec = garchSpec(model)
   	spec@title = "normal-GARCH(1,1)"
   	spec
   	ts = 100 * garchSim(n = 1000, spec = spec)[,"x"]
   	plot(x = ts, type = "l", main = spec@title, col = "steelblue3",)
   	abline(h = 0, col = "grey")
   	grid()
   	
    Continue <- readline("Press any key > ") 

# APARCH: TS-GARCH(1,5[1,5])

	# Taylor-Schwert's subset-normal-TS-GARCH(1,5[1,5]), delta=1:
   	model = list(omega = 1e-6, alpha = 0.1, gamma = 0, alpha.lags = 1, 
   		beta = c(0.5, 0.3), beta.lags = c(1, 5), mu = 0, delta=1)	
   	ts <- aparchSim(, innov = rnorm(5000), start.innov = rnorm(1000))
    tsPlot(ts, main = "Subset normal-TS-GARCH(1,5[1,5])")
    
    Continue <- readline("Press any key > ") 
    
    
    ################
    
    
    # Simulate ARCH(2):
	model = list(omega = 1e-6, alpha = c(0.1, 0.3))
	model
	x = garchSim(model, n = 1000)
	class(x)
	tsPlot(x, main = "ARCH(2) Model")

# Simulate GARCH(1,1):
	model = list(omega = 1e-6, alpha = 0.1, beta = 0.8, mu = 0.01)
	model
	x <- garchSim(model, n=1000)
	tsPlot(x, main = "GARCH(1,1) Model")
	
	
	
	
	#####################
	
	
	
	
	

# Settings:
  	par(mfrow=c(4,2))


# Bollerslev's GARCH(1,2) model
# h(t)^2 = omega + alpha(i)*eps(t-i)^2 + beta(i)*h(t-i)^2 
# --------------------------------------------------------------------

    # Simulate:
  	set.seed(4661)
  	x <- aparchSim(model=list(omega=0.8e-5, 
  	    alpha=0.1, gamma=0, alpha.lags=1, 
  	    beta=c(0.5, 0.2), beta.lags=c(1, 2), delta=2), 
  	    innov=rnorm(5000), start.innov=rnorm(1000))
	# Plot:
	ts.plot(x, ylim=c(-0.04, 0.04), main="GARCH(1,2) Model")
	acf(abs(x), lag.max=20)
	var(x)


# Taylor-Schwert GARCH(p,q) model
# h(t) = omega + alpha(i)*|eps(t-i)| + beta(j)*h(t-j) 
# --------------------------------------------------------------------

    # Simulate:
  	set.seed(5066)
  	x <- aparchSim(model=list(omega=1.4e-3, 
  	    alpha=0.1, gamma=0, alpha.lags=1, 
  	    beta=c(0.5, 0.2), beta.lags=c(1, 2), delta=1), 
  	    innov=rnorm(5000), start.innov=rnorm(1000))
	# Plot:
	ts.plot(x, ylim=c(-0.04, 0.04), main="TS GARCH(1,2) Model")
	acf(abs(x), lag.max=20)
	var(x)
	

# Glosten's et al. GJR(p,q) model
# h(t)^2 = omega + alpha'(i)*eps(t-i)^2 + beta(j)*h(t-j)^2 + gamma'(i)*S(i)*eps(t-i)
#   S(i) = 1 if eps(t-i)>0 ; S(i) = 0 otherwise
#   0 <= gamma(i) < 1
# --------------------------------------------------------------------

    # Simulate:
  	set.seed(8531)
  	x <- aparchSim(model=list(omega=0.8e-5, 
  		alpha=0.1, gamma=0.2, alpha.lags=1,
  		beta=c(0.5, 0.2), beta.lags=c(1, 2), delta=2), 
  		innov=rnorm(5000), start.innov=rnorm(1000))
	# Plot:
	ts.plot(x, ylim=c(-0.04, 0.04), main="GJR(1,2) Model")
	acf(abs(x), lag.max=20)
	var(x)


	
	
	
	##############################
	
	
	
	
	 

# GARCH(1,1) Normal Innovations
	# Simulate:
  	set.seed(4711)
  	x <- aparchSim(model=list(omega=1.02e-6, 
  		alpha=0.1, gamma=0, alpha.lags=1, beta=0.6, beta.lags=1, delta=2), 
  		innov=rnorm(2000), start.innov=rnorm(100))
	# Plot:
	ts.plot(x, ylim=c(-0.015, 0.015), main="GARCH(1,1) Normal Innovations")
	acf(abs(x), lag.max = 12)
	var(x)

# GARCH(1,1) Student t(4) Innovations
	# Simulate:
  	set.seed(1858)
  	x <- aparchSim(model=list(omega = 0.45e-6, alpha = 0.1, gamma = 0, 
  		alpha.lags = 1, beta = 0.6, beta.lags = 1, delta = 2), 
  		innov=rt(2000, df=4), start.innov=rt(100, df=4))
	# Plot:
	ts.plot(x, ylim = c(-0.015, 0.015), 
		main = "GARCH(1,1) - Student t(4) Innovations")
	acf(abs(x), lag.max = 12)
	var(x)
	
# GARCH(1,1) 1.9 Stable Innovations
	# Simulate:
  	set.seed(2245)
  	x <- aparchSim(model = list(omega = 0.3e-6, alpha = 0.1, gamma = 0, 
  		alpha.lags = 1, beta = 0.6, beta.lags = 1, delta = 2), innov = 
  		rsymstb(2000, alpha=1.95), start.innov=rsymstb(100, alpha=1.95))
	# Plot:
	ts.plot(x, ylim = c(-0.015, 0.015), 
		main = "GARCH(1,1) - 1.95 Stable Innovations")
	acf(abs(x), lag.max = 12)
	var(x)
	


    
    