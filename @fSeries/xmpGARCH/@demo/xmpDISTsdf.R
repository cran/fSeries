
		
sdfMoments =
function(method = c("norm", "ged", "std"), mean=0, sd=1, nu=5, xi = 1)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	#	Compute different Moment statistics for standardized 
	#	distribution functions
	
	# Details:
	#	norm - Normal Distribution ("nu" unused)
	#	ged  - Generalized Error Distribution
	#   std  - Student-t Distribution
	
	# FUNCTION:
	
	# Select:
	method = method[1]
	
	# Normal Distribution:
	if (method == "norm") {
		# Standardized Centered Moments:
		cM = c(0, 1, 0, 3)
		
		# Standardized Absolute Moments:
		# NOTE: m[2n-1] = 2^n*Gamma(n)/sqrt(2*pi)
		aM = c(2/sqrt(2*pi), 1, 4/sqrt(2*pi), 3)
		
		# Skewed Standardized Centered Moments:
		# mu = aM[1]*(xi-1/xi)
		# sigma = sqrt( (aM[2]-aM[1]^2)*(xi^2+1/xi^2) + 2*aM[1]^2 - aM[2] )
		xiM = function(m, xi){
			xi.r = rep(NA, times = 4)
			for (r in 1:4){
				xi.r[r] = m[r]*( xi^(r+1) + (-1)^r / xi^(r+1) ) / (xi + 1/xi) }
			xi.r }
		MOM = xiM(m = aM, xi=xi)
		scM = (MOM[1] - mu) / sigma * sd
		scM = c(scM, (MOM[2] - 2*MOM[1]*mu + mu^2) / sigma^2 * sd^2 )
		scM = c(scM, (MOM[3] - 3*MOM[2]*mu + 3*MOM[1]*mu^2 -
			 mu^3) / sigma^3 * sd^3 )
	    scM = c(scM, (MOM[4] - 4*MOM[3]*mu + 6*MOM[2]*mu^2 - 
	    	4*MOM[1]*mu^3 + mu^4) / sigma^4 * sd^4	)
	    
	    # Result:	
	    result = cbind(cM, aM, scM)
	}
	
	# Generalized Error Distribution:
	if (method == "ged") {
		# Standardized Centered Moments:
		cM = c(0, 1, 0, 3)
		
		# Standardized Absolute Moments:
		# NOTE: m[2n-1] = 2^n*Gamma(n)/sqrt(2*pi)
		aM = c(2/sqrt(2*pi), 1, 4/sqrt(2*pi), 3)
		
		# Skewed Standardized Centered Moments:
		# mu = aM[1]*(xi-1/xi)
		# sigma = sqrt( (aM[2]-aM[1]^2)*(xi^2+1/xi^2) + 2*aM[1]^2 - aM[2] )
		xiM = function(m, xi){
			xi.r = rep(NA, times = 4)
			for (r in 1:4){
				xi.r[r] = m[r]*( xi^(r+1) + (-1)^r / xi^(r+1) ) / (xi + 1/xi) }
			xi.r }
		MOM = xiM(m = aM, xi=xi)
		scM = (MOM[1] - mu) / sigma * sd
		scM = c(scM, (MOM[2] - 2*MOM[1]*mu + mu^2) / sigma^2 * sd^2 )
		scM = c(scM, (MOM[3] - 3*MOM[2]*mu + 3*MOM[1]*mu^2 -
			 mu^3) / sigma^3 * sd^3 )
	    scM = c(scM, (MOM[4] - 4*MOM[3]*mu + 6*MOM[2]*mu^2 - 
	    	4*MOM[1]*mu^3 + mu^4) / sigma^4 * sd^4	)
	    
	    # Result:	
	    result = cbind(cM, aM, scM)
	}
	
	# Student-t Distribution
	if (method == "std") {
		# Standardized Centered Moments:
		cM = c(0, 1, 0, 3)
		
		# Standardized Absolute Moments:
		# NOTE: m[2n-1] = 2^n*Gamma(n)/sqrt(2*pi)
		aM = c(2/sqrt(2*pi), 1, 4/sqrt(2*pi), 3)
		
		# Skewed Standardized Centered Moments:
		# mu = aM[1]*(xi-1/xi)
		# sigma = sqrt( (aM[2]-aM[1]^2)*(xi^2+1/xi^2) + 2*aM[1]^2 - aM[2] )
		xiM = function(m, xi){
			xi.r = rep(NA, times = 4)
			for (r in 1:4){
				xi.r[r] = m[r]*( xi^(r+1) + (-1)^r / xi^(r+1) ) / (xi + 1/xi) }
			xi.r }
		MOM = xiM(m = aM, xi=xi)
		scM = (MOM[1] - mu) / sigma * sd
		scM = c(scM, (MOM[2] - 2*MOM[1]*mu + mu^2) / sigma^2 * sd^2 )
		scM = c(scM, (MOM[3] - 3*MOM[2]*mu + 3*MOM[1]*mu^2 -
			 mu^3) / sigma^3 * sd^3 )
	    scM = c(scM, (MOM[4] - 4*MOM[3]*mu + 6*MOM[2]*mu^2 - 
	    	4*MOM[1]*mu^3 + mu^4) / sigma^4 * sd^4	)
	    
	    # Result:	
	    result = cbind(cM, aM, scM)
	}
	
	# Return Value:
	result
}
