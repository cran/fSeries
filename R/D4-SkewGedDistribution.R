
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file

################################################################################
# FUNCTION:             GED Distribution:
#  dged                   Density for the Generalized Error Distribution
#  pged                   Probability function for the GED
#  qged                   Quantile function for the GED
#  rged                   Random Number Generator for the GED
# FUNCTION:             DESCRIPTION:
#  dsged                  Density for the skewed GED
#  psged                  Probability function for the skewed GED
#  qsged                  Quantile function for the skewed GED
#  rsged                  Random Number Generator for the skewed GED
################################################################################


################################################################################
# GED - Generalized Error Distribution:


dged =
function(x, mean = 0, sd = 1, nu = 2)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Compute the density for the 
    #   generalized error distribution.
    
    # FUNCTION:
    
    # Compute Density:
    z = (x - mean ) / sd
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    result = g * exp (-0.5*(abs(z/lambda))^nu) / sd
    
    # Return Value
    result
}

        
# ------------------------------------------------------------------------------


pged = 
function(q, mean = 0, sd = 1, nu = 2)
{   # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Compute the probability for the  
    #   generalized error distribution.
    
    # FUNCTION:
        
    # Compute Probability:
    q = (q - mean ) / sd
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    h = 2^(1/nu) * lambda * g * gamma(1/nu) / nu
    s = 0.5 * ( abs(q) / lambda )^nu
    result = 0.5 + sign(q) * h * pgamma(s, 1/nu)
    
    # Return Value:
    result
}
        

# ------------------------------------------------------------------------------


qged =
function(p, mean = 0, sd = 1, nu = 2)
{   # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Compute the quantiles for the  
    #   generalized error distribution.
    
    # FUNCTION:
    
    # Compute Quantiles:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    q = lambda * (2*qgamma((abs(2*p-1)), 1/nu))^(1/nu)
    result = q*sign(2*p-1) * sd + mean
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------

    
rged = 
function(n, mean = 0, sd = 1, nu = 2)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Generate GED random deviates. The function uses the 
    #   method based on the transformation of a Gamma random 
    #   variable.
    
    # FUNCTION:
    
    # Generate Random Deviates:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    # print(lambda)
    r = rgamma(n, 1/nu)
    z =  lambda * (2*r)^(1/nu) * sign(runif(n)-1/2)
    result = z * sd + mean

    
    # Return Value:
    result
}


################################################################################
# Skewed Generalized Error Distribution


.dsged = 
function(x, nu, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Standardize:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = x*sigma + mu  
	
	# Compute:
	Xi = xi^sign(z)
	g = 2  / (xi + 1/xi)	
	Density = g * dged(x = z/Xi, nu=nu)  
	
	# Return Value:
	Density * sigma 
}

	  
dsged =
function(x, mean = 0, sd = 1, nu = 2, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Compute the density function of the 
	#	skewed generalized error distribution
	
	# FUNCTION:
		
    # Shift and Scale:
	result = .dsged(x = (x-mean)/sd, nu = nu, xi = xi) / sd
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------


.psged =
 function(q, nu, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Standardize:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = q*sigma + mu
	
	# Compute:	
	Xi = xi^sign(z)
	g = 2  / (xi + 1/xi)	
	Probability = H(z) - sign(z) * g * Xi * pged(q = -abs(z)/Xi, nu=nu)
	
	# Return Value:
	Probability 
}

	  
psged =
function(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Compute the distribution function of the 
	#	skewed generalized error distribution
	
	# FUNCTION:
			  
	# Shift and Scale:
	result = .psged(q = (q-mean)/sd, nu = nu, xi = xi)
		  
	# Return Value:
	result
}


# ------------------------------------------------------------------------------	


.qsged =
function(p, nu, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Standardize:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	
	# Compute:	
	g = 2  / (xi + 1/xi)
	sig = sign(p-1/2) 
	Xi = xi^sig		  
	p = (H(p-1/2)-sig*p) / (g*Xi)
	Quantile = (-sig*qged(p=p, sd=Xi, nu=nu) - mu ) / sigma
	
	# Return Value:
	Quantile 
}

	  	
qsged =
function(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Compute the quantile function of the 
	#	skewed generalized error distribution
	
	# FUNCTION:
		
	# Shift and Scale:
	result = .qsged(p = p, nu = nu, xi = xi) * sd + mean
	
	# Return Value:
	result
}

	
# ------------------------------------------------------------------------------
	

.rsged =
function(n, nu, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Generate Random Deviates:
	weight = xi / (xi + 1/xi)
	z = runif(n, -weight, 1-weight)
	Xi = xi^sign(z)
	Random = -abs(rged(n, nu=nu))/Xi * sign(z)	
	
	# Scale:
	lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
	g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
	m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
	mu = m1*(xi-1/xi)
	sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	Random = (Random - mu ) / sigma	
	
	# Return value:
	Random 
}


rsged =
function(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Generate random deviates from the 
	#	skewed generalized error distribution
	
	# FUNCTION:
		
	# Shift and Scale:
	result = .rsged(n = n, nu = nu, xi = xi) * sd + mean
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------

