

nleqnsFit = 
function(formulas, data = list(), 
method = c("OLS", "WLS", "SUR", "2SLS", "W2SLS", "3SLS", "W3SLS"), 
title = "", description = "", ...)
{	# A function implemented by Diethelm Wuertz

    # Description:
    #   Common function call for several system equation models.
    
    # Arguments:
    #	formulas - the list of formulas describing the system of 
    #		equations
    #	data - the input data set in form of a 'data.frame' or 
	#		'timeSeries' object
	#	method - a character string describing the desired method, 
	#		one of: "OLS", "WLS", "SUR", "2SLS", "W2SLS", "3SLS", 
	#		or "W3SLS".
	#	title - a character string which allows for a project title
	#	description - a character string which allows for a brief 
	#		description
	#	... - additional optional arguments to be passed to the 
	#		underlying function 'systemfit' 
	
	# Value:
	# The function 'eqnaFit' returns an object of class "fEQNS"
  	# with the following slots:
	#	@call - the matched function call
  	#	@data - the input data in form of a 'data.frame' or a 
  	#		'timeSeries' object
	#	@description - a character string which allows for a brief 
	#		project description
	#	@method - the character string describing the desired method
  	#	@formulas - the list of formulas describing the system of 
  	#		equations
	#	@title - a character string which allows for a project title
	# 	@fit - a summary of the  results as a list returned from the 
	#		underlying functions from the 'systemfit' package.	

    # FUNCTION:
    
    # Fit:
    fit = nlsystemfit(method = method[1], eqns = formulas, 
    	data = as.data.frame(data), ...)	
    class(fit) = c("list", "nlsystemfit.system")
    
    # Add:
    fit$coef = unclass(.coef.systemfit(fit))
    fit$confint = unclass(.confint.systemfit(fit))
    fit$fitted = unclass(.fitted.systemfit(fit))
    fit$residuals = unclass(.residuals.systemfit(fit))
    fit$vcov = unclass(.vcov.systemfit(fit))
    
    # Default Title:
    if (title == "") title = paste(method[1], "Estimation")
    
	# Return Value:
    new("fEQNS",     
        call = as.call(match.call()),
        formulas = formulas, 
        data = as.data.frame(data),
        method = as.character(method[1]), 
        fit = fit,
        title = as.character(title), 
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


NLSUR = 
function(formulas, data = list(), ...)
{
	# Fit:
	ans = nleqnsFit(formulas = formulas, data = data, method = "SUR", ...)
		
	# Return Value:
	ans
}


# ------------------------------------------------------------------------------