
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
# FUNCTION:		    		DESCRIPTION:
#  rollFun					 Compute Rolling Function Value
#   rollMean				  Compute Rolling Mean
#	rollVar				  	  Compute Rolling Variance
#	rollMin				      Compute Rolling Minimum
#	rollMax				      Compute Rolling Maximum
################################################################################


rollFun = 
function(x, n, FUN, ...) 
{	# A function imlemented by Diethelm Wuertz

	# Description:
	#	Compute rolling function value

	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Roll FUN:
	start = 1
	end = length(x)-n+1
	m = x[start:end]
	for (i in 2:n) {
		start = start + 1
		end = end + 1
		m = cbind(m, x[start:end])}
	
	# Result:
	ans = apply(m, MARGIN = 1, FUN = FUN, ...)
	
	# Return value:
	ans
}


# ------------------------------------------------------------------------------
	

rollMean = 
function(x, n = 9, trim = TRUE, na.rm = FALSE) 
{	# A function imlemented by Diethelm Wuertz

	# Description:
	#	Compute rolling mean

	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Roll Mean:
	if (na.rm) x = as.vector(na.omit(x))
	rmean = rollFun(x = x, n = n, FUN = mean) 
	if (!trim) rmean = c(rep(NA, (n-1)), rmean)
		
	# Return Value:
	rmean
}
	
	
# ------------------------------------------------------------------------------


rollVar  = 
function(x, n = 9, trim = TRUE, unbiased = TRUE, na.rm = FALSE) 
{ 	# A function imlemented by Diethelm Wuertz

	# Description:
	#	Compute rolling variance

	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Roll Var:
	if (na.rm) x = as.vector(na.omit(x))
	rvar = rollFun(x = x, n = n, FUN = var) 
	if (!unbiased) rvar = (rvar * (n-1))/n
	if (!trim) rvar = c(rep(NA, (n-1)), rvar)
	
	# Return Value:
	rvar 
}
	

# ------------------------------------------------------------------------------


rollMax  = 
function(x, n = 9, trim = TRUE, na.rm = FALSE) 
{	# A function imlemented by Diethelm Wuertz

	# Description:
	#	Compute rolling maximum

	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Roll Max:
	if (na.rm) x = as.vector(na.omit(x))
	rmax = rollFun(x = x, n = n, FUN = max)
	if (!trim) rmax = c(rep(NA, (n-1)), rmax)
	
	# Return Value:
	rmax 
}
	

# ------------------------------------------------------------------------------

rollMin  = 
function(x, n = 9, trim = TRUE, na.rm = FALSE) 
{ 	# A function imlemented by Diethelm Wuertz

	# Description:
	#	Compute rolling function minimum

	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Roll Min:
	if (na.rm) x = as.vector(na.omit(x))
	rmin = rollFun(x = x, n = n, FUN = min) 
	if (!trim) rmin = c(rep(NA, (n-1)), rmin)
	
	# Return Value:
	rmin 
}


# ******************************************************************************

