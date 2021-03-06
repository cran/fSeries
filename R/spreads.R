
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
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 FINANCIAL TIME SERIES:
#  midquotes                 Computes mid quotes from a 'timeSeries' object
#  spreads                   Computes spreads from a 'timeSeries' object
# OLD FUNCTIONS:
#  midquoteSeries <- midquotes
#  spreadSeries <- spreads
################################################################################


midquotes = 
function(x, which = c("Bid", "Ask"))
{   
    # Compute Mid Quotes:
    midQuotes = 0.5 * ( x[, which[1]] + x[, which[2]] ) 
    
    # Return Value:
    midQuotes
}


# ------------------------------------------------------------------------------


midquoteSeries =
function(...)
{
    midquotes(...)
}


# ------------------------------------------------------------------------------


spreads = 
function(x, which = c("Bid", "Ask"), tickSize = NULL)
{   
    # Compute Spread:
    Spread = x[, which[2]] - x[, which[1]] 
    if (!is.null(tickSize)) Spread@Data = round(Spread@Data/tickSize)
    
    # Return Value:
    Spread
}


# ------------------------------------------------------------------------------


spreadSeries =
function(...)
{
    spreads(...)
}


################################################################################

