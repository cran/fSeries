
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
# FUNCTION:             NORMAL Distribution - part of R's base Package:
#  * dnorm                R: Density for the Normal Distribution
#  * pnorm                R: Probability function for the Normal Distribution
#  * qnorm                R: Quantile function for the Normal Distribution
#  * rnorm                R: Random Number Generator for the Normal Distribution  
# FUNCTION:             SKEW-DISTRIBUTIONS DESCRIPTION:
#  dsnorm                 Density for the skew normal Distribution
#  psnorm                 Probability function for the skew NORM
#  qsnorm                 Quantile function for the skew NORM
#  rsnorm                 Random Number Generator for the skew NORM
################################################################################


################################################################################
# Skewed Normal Distribution



.dsnorm = 
function(x, xi) 
{   # A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:

	# Standardize:
	m1 = 2/sqrt(2*pi)
	mu = m1 * (xi - 1/xi)
	sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	z = x*sigma + mu  
	# Compute:
	Xi = xi^sign(z)
	g = 2 / (xi + 1/xi)	
	Density = g * dnorm(x = z/Xi)  
	# Return Value:
	Density * sigma 
}
		

dsnorm =
function(x, mean = 0, sd = 1, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Compute the density function of the 
	#	skew normal distribution
	
	# Arguments:
	# 	x - a numeric vector of quantiles.
	#	mean, sd, xi - location parameter, scale parameter, and 
	#		skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .dsnorm(x = (x-mean)/sd, xi = xi) / sd
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------


.psnorm =
function(q, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Standardize:
	  m1 = 2/sqrt(2*pi)
	  mu = m1 * (xi - 1/xi)
	  sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	  z = q*sigma + mu
	# Compute:	
	  Xi = xi^sign(z)
	  g = 2  / (xi + 1/xi)	
	  Probability = H(z) - sign(z) * g * Xi * pnorm(q = -abs(z)/Xi)
	# Return Value:
	  Probability }
	

psnorm =
function(q, mean = 0, sd = 1, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Compute the distribution function of the 
	#	skew normal distribution
	
	# Arguments:
	# 	q - a numeric vector of quantiles.
	#	mean, sd, xi - location parameter, scale parameter, and 
	#		skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .psnorm(q = (q-mean)/sd, xi = xi)
		  
	# Return Value:
	result
}


# ------------------------------------------------------------------------------	


.qsnorm =
function(p, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Standardize:
	  m1 = 2/sqrt(2*pi)
	  mu = m1 * (xi - 1/xi)
	  sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	# Compute:	
	  g = 2  / (xi + 1/xi)
	  sig = sign(p-1/2) 
	  Xi = xi^sig		  
	  p = (H(p-1/2)-sig*p) / (g*Xi)
	  Quantile = (-sig*qnorm(p=p, sd=Xi) - mu ) / sigma
	# Return Value:
	  Quantile 
}
	
	
qsnorm =
function(p, mean = 0, sd = 1, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Compute the quantile function of the 
	#	skew normal distribution
	
	# Arguments:
	#   p - a numeric vector of probabilities.
	#	mean, sd, xi - location parameter, scale parameter, and 
	#		skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .qsnorm(p = p, xi = xi) * sd + mean
	
	# Return Value:
	result
}

	
# ------------------------------------------------------------------------------


.rsnorm =
function(n, xi) 
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Internal Function
	
	# FUNCTION:
	
	# Generate Random Deviates:
	  weight = xi / (xi + 1/xi)
	  z = runif(n, -weight, 1-weight)
	  Xi = xi^sign(z)
	  Random = -abs(rnorm(n))/Xi * sign(z)	
	# Scale:
	  m1 = 2/sqrt(2*pi)
	  mu = m1 * (xi - 1/xi)
	  sigma = sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
	  Random = (Random - mu ) / sigma	
	# Return value:
	  Random 
}
		

rsnorm =
function(n, mean = 0, sd = 1, xi = 1.5)
{	# A function implemented by Diethelm Wuertz 

	# Description:
	# 	Generate random deviates from the 
	#	skew normal distribution
	
	# Arguments:
	#	n - an integer value giving the number of observation.
	#	mean, sd, xi - location parameter, scale parameter, and 
	#		skewness parameter.
	
	# FUNCTION:
	
	# Shift and Scale:
	result = .rsnorm(n = n, xi = xi) * sd + mean
	
	# Return Value:
	result
}


# ------------------------------------------------------------------------------

