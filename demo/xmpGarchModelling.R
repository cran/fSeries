
#
# Example: 
#	Model a GARCH time series process 
#
# Description:
#	PART I: Estimate GARCH models of the following type ARCH(2) 
#     and GARCH(1,1) with normal conditional distribution functions.
#   PART II: Simulate GARCH models of the following type, ARCH(2) 
#     and GARCH(1,1),
#	with normal conditional distribution functions.
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------

	
# PART I: Estimation:

	# Settings:
	set.seed(547)
    # Bollerslev's GARCH(1,1) with normal innovations:
	model = list(omega = 1e-6, alpha = 0.1, beta = 0.8, mu = 0)
	x = garchSim(model, n = 1000)
	fit = garchFit(as.numeric(x), order = c(1, 1))
	print(fit)
	# Summary and Diagnostic Analysis:
	summary(fit)
	# Plot Results:
	par(mfrow = c(2, 2))
	plot(fit)
	###


# PART II: Simulation

	# Settings:
  	par(mfrow = c(2, 1))
	set.seed(411)
	# Simulate ARCH(2):
	model = list(omega = 1e-6, alpha = c(0.1, 0.3), h0 = 1e-7)
	model
	x = garchSim(model)
	ts.plot(x, main = "ARCH(2) Model")
	# Simulate GARCH(1,1):
	model = list(omega = 1e-6, alpha = 0.1, beta = 0.8, h0 = 1e-7)
	model
	x = garchSim(model)
	ts.plot(x, main = "GARCH(1,1) Model")
	###
	
	