

# 1-Step ahead Predictor AR(1):

	y0 = rnorm(1)
	a1hat = 0.4
	muhat = 0
	y1 = muhat + a1hat*(y0-muhat)	
	print(c(y0, y1))
	
	
# h-Step ahead Predictor AR(1):
	h = 1:5
	y1 = muhat + a1hat^h * (y0-muhat) 
	print(c(y0, y1))
	
	
# 1 Step ahead Predictor for GARCH(1,1)

	sigma1 = omegaHat + alpha1Hat * eps0^2 + beta1Hat * sigma0
	
	#  	eps2[t|t] = eps2[t] 
	#	