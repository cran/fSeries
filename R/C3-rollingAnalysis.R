
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
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
# for the code accessed (or partly included) from other R-ports:
#   R: see R's copyright and license file
#   ts: collected by Brian Ripley. See SOURCES
#   tseries: Compiled by Adrian Trapletti <a.trapletti@bluewin.ch>
#   fracdiff: S original by Chris Fraley <fraley@stat.washington.edu>
#     R-port: by Fritz Leisch <leisch@ci.tu-wien.ac.at>
#     since 2003-12: Martin Maechler
#   lmtest: Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>
#     Achim Zeileis <zeileis@ci.tuwien.ac.at>
#     David Mitchell
#   mda: S original by Trevor Hastie & Robert Tibshirani
#     R port by Friedrich Leisch, Kurt Hornik and Brian D. Ripley
#   mgcv: Simon Wood <simon@stats.gla.ac.uk>
#   modreg: Brian Ripley and the R Core Team
#   polspline: Charles Kooperberg <clk@fhcrc.org>
#   nnet: S original by Venables & Ripley. 
#     R port by Brian Ripley <ripley@stats.ox.ac.uk>
#       following earlier work by Kurt Hornik and Albrecht Gebhardt


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

	# FUNCTION:
	
	# Roll FUN:
	start = 1
	end = length(x)-n+1
	m = x[start:end]
	for (i in 2:n) {
		start = start + 1
		end = end + 1
		m = cbind(m, x[start:end])}
	
	# Return Value:
	apply(m, MARGIN = 1, FUN = FUN, ...)
}


# ------------------------------------------------------------------------------
	

rollMean = 
function(x, n = 9, trim = TRUE, na.rm = FALSE) 
{	# A function imlemented by Diethelm Wuertz

	# FUNCTION:
	
	# Roll Mean:
	if(na.rm) x = as.vector(na.omit(x))
	rvar = rollFun(x = x, n = n, FUN = mean) 
	if(!trim) rmean = c(rep(NA, (n-1)), rmean)
		
	# Return Value:
	rollFun(x = x, n = n, FUN = mean) 
}
	
	
# ------------------------------------------------------------------------------


rollVar  = 
function(x, n = 9, trim = TRUE, unbiased = TRUE, na.rm = FALSE) 
{ 	# A function imlemented by Diethelm Wuertz

	# FUNCTION:
	
	# Roll Var:
	if(na.rm) x = as.vector(na.omit(x))
	rvar = rollFun(x = x, n = n, FUN = var) 
	if(!unbiased) rvar = (rvar * (n-1))/n
	if(!trim) rvar = c(rep(NA, (n-1)), rvar)
	
	# Return Value:
	rvar 
}
	

# ------------------------------------------------------------------------------


rollMax  = 
function(x, n = 9, trim = TRUE, na.rm = FALSE) 
{	# A function imlemented by Diethelm Wuertz

	# FUNCTION:
	
	# Roll Max:
	if(na.rm) x = as.vector(na.omit(x))
	rmax = rollFun(x = x, n = n, FUN = max)
	if(!trim) rmax = c(rep(NA, (n-1)), rmax)
	
	# Return Value:
	rmax 
}
	

# ------------------------------------------------------------------------------

rollMin  = 
function(x, n = 9, trim = TRUE, na.rm = FALSE) 
{ 	# A function imlemented by Diethelm Wuertz

	# FUNCTION:
	
	# Roll Min:
	if(na.rm) x = as.vector(na.omit(x))
	rmin = rollFun(x = x, n = n, FUN = min) 
	if(!trim) rmin = c(rep(NA, (n-1)), rmin)
	
	# Return Value:
	rmin 
}


# ------------------------------------------------------------------------------

