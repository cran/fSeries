
#
# Example: Garch Model Parameter Estimation
#	Estimate GARCH models of the ARCH(2) and GARCH(1,1) 
#	type with normal conditional distribution functions.
#


# Settings:
	require(fBasics)
	set.seed(547)
	
# Bollerslev's GARCH(1,1) with normal innovations:
    model = list(omega = 1e-6, alpha = 0.1, beta = 0.8)
	x = garchSim(model, n = 500)
	fit = garchFit(formula.var = ~garch(1, 1))
	print(fit)
	
# Summary and Diagnostic Analysis:	
	summary(fit)
	
# Plot Results:
	par(mfrow = c(4, 2))
	plot(fit)
	
	