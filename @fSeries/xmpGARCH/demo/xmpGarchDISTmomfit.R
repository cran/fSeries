
#
# Title: 
#	Parameter Fit of a Density Function
#   
# Description: 
#	Fits the parameters of a distribution function by the
#	maximum log-likelihood approach using the "nlm"
#	optimizer from R's base package.

	
# ------------------------------------------------------------------------------


# Load Data Set:

	data(demgbp)
	demgbp = demgbp[,1]/100
	
	
# ------------------------------------------------------------------------------


# Basic Statistics:

	x = demgbp

	skewness = function(x) { sum((x - mean(x))^3/sqrt(var(x))^3)/length(x) }
	kurtosis = function(x) { sum((x - mean(x))^4/var(x)^2)/length(x) - 3 }
	cl.vals = function(x, ci) {
        t.val = qt((1 - ci)/2, n - 1)
        lcl = mean(x) + sqrt(var(x)/length(x)) * t.val
        ucl = mean(x) - sqrt(var(x)/length(x)) * t.val
        c(lcl, ucl) }

	z = c(x.length = length(x), min(x), max(x), as.numeric(quantile(x, 
        prob = 0.25, na.rm = TRUE)), as.numeric(quantile(x, prob = 0.75, 
        na.rm = TRUE)), mean(x), median(x), sum(x), sqrt(var(x)/length(x)), 
        cl.vals(x, 0.95)[1], cl.vals(x, 0.95)[2], var(x), sqrt(var(x)), 
        skewness(x), kurtosis(x))
    
    znames <- c(" nobs", " Minimum", " Maximum", " 1. Quartile", 
        " 3. Quartile", " Mean", " Median", " Sum", " SE Mean", 
        " LCL Mean", " UCL Mean", " Variance", " Stdev", 
        " Skewness", " Kurtosis")
        
    cat("\nBasic Statistics:\n")
    data.frame(z, row.names = znames)

    	
# ------------------------------------------------------------------------------

	
# Estimate the Parameters:

	fit = dFit(x = demgbp, density = dsstd, parm = c(mean = mean(x), 
		sd = mad(x), nu = 2, xi = 0), trace = TRUE, doplot = TRUE) 
		
	print(fit)
	
	
# ------------------------------------------------------------------------------


# Parameter Estimation by fitting the moments

	
gedFit = 
function(x, method = c("mom"," mle")) 
{	# A function implemented by Diethelm Wuertz

	# Compute nu:
	
	# Internal Functions:
	modelKurtosis = 
		function(nu) { gamma(1/nu)*gamma(5/nu)/gamma(3/nu)^2 - 3 }
	kurtosis = 
		function(x) { sum((x - mean(x))^4/var(x)^2)/length(x) - 3 }
	fkappa = 
		function(x, K) { modelKurtosis(x)  - K }
		
	# Estimate:	
	nu = c(1/8, 1/4, 1/2, 1, 2, 4, 8)	
	index = max(cumsum(sign(modelKurtosis(nu)-kurtosis(x))))	
	nu = c(nu, Inf)
	nu.estimated = uniroot(fkappa, interval = c(nu[index], nu[index+1]), 
		K = kurtosis(x))$root
	
	# Return value:
	fit = c(mean = mean(x), sd = sqrt(var(x)), nu = nu.estimated)
	
	if (method[1] == "mle") {
		fit = dFit(x = x, density = "dged", parm = fit, 
		trace = TRUE, doplot = TRUE) }	
		
	fit
	
	
}


	# Compute Excess Kurtosis:
	kappa = kurtosis(demgbp)
	kappa
	
	# Estimate Parameters:
	fit = gedFit(demgbp)
	fit
	
	# Test:
	gamma(1/nu)*gamma(5/nu)/gamma(3/nu)^2 - 3