
#
# Example: 
#	Aparch Model Simulations
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


# Settings:
  	par(mfrow = c(3, 2), cex = 0.6)

# Bollerslev's GARCH(1,2) model
# h(t)^2 = omega + alpha(i)*eps(t-i)^2 + beta(i)*h(t-i)^2 

    # Simulate:
  	set.seed(4661)
  	x = aparchSim(model = list(omega=0.8e-5, 
  	    alpha = 0.1, gamma = 0, alpha.lags=1, 
  	    beta = c(0.5, 0.2), beta.lags = c(1, 2), delta = 2), 
  	    innov = rnorm(5000), start.innov = rnorm(1000))
	# Plot:
	ts.plot(x, ylim = c(-0.04, 0.04), main = "GARCH(1,2) Model")
	acf(abs(x), lag.max = 20)
	var(x)

# Taylor-Schwert GARCH(p,q) model
# h(t) = omega + alpha(i)*|eps(t-i)| + beta(j)*h(t-j) 

    # Simulate:
  	set.seed(5066)
  	x = aparchSim(model=list(omega = 1.4e-3, 
  	    alpha = 0.1, gamma = 0, alpha.lags = 1, 
  	    beta=c(0.5, 0.2), beta.lags = c(1, 2), delta = 1), 
  	    innov = rnorm(5000), start.innov = rnorm(1000))
	# Plot:
	ts.plot(x, ylim = c(-0.04, 0.04), main = "TS GARCH(1,2) Model")
	acf(abs(x), lag.max = 20)
	var(x)
	
# Glosten's et al. GJR(p,q) model
# h(t)^2 = omega + alpha'(i)*eps(t-i)^2 + beta(j)*h(t-j)^2 + gamma'(i)*S(i)*eps(t-i)
#   S(i) = 1 if eps(t-i)>0 ; S(i) = 0 otherwise
#   0 <= gamma(i) < 1

    # Simulate:
  	set.seed(8531)
  	x = aparchSim(model=list(omega=0.8e-5, 
  		alpha = 0.1, gamma = 0.2, alpha.lags = 1,
  		beta = c(0.5, 0.2), beta.lags = c(1, 2), delta = 2), 
  		innov=rnorm(5000), start.innov=rnorm(1000))
	# Plot:
	ts.plot(x, ylim = c(-0.04, 0.04), main = "GJR(1,2) Model")
	acf(abs(x), lag.max = 20)
	var(x)

