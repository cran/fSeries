
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
# Copyright (C)				Adrian Trapletti
# Copyright (C) 1998-2003  	 Diethelm Wuertz for this R-port
# Topic						 fSeries - TRADING INDICATORS AND TECHNICAL ANALYSIS
# ------------------------------------------------------------------------------
# FUNCTION:		    		BENCHMARK ANALYSIS FUNCTIONS:
#  maxDrawDown				 Computes the maximum drawdown
#  sharpeRatio				 Calculates the Sharpe Ratio
#  sterling Ratio			 Calculates the Sterling Ratio
#  ohlcPlot					 Creates a Open-High-Low-Close plot
################################################################################


getReturns = 
function(x, type = c("continuous", "discrete"), percentage = FALSE, 
trim = TRUE)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	Computes return series given a financial price series.
	
	# Note:
	#	Function for S-Plus Compatibility
	
	# Arguments:	Description:
	#	x			numeric vector
	#	type		not used [continuous], for Splus compatibility
	#	percentage 	FALSE, if TRUE the series will be expressed in %
	#	trim		not used [FALSE], for Splus compatibility	
	
	# FUNCTION:
	
	# Settings:
	type = type[1]
	
	# Check object:
	series = class(x)
	
	# timeSeries Object:
	if (series == "timeSeries") {
		pos = positions(x)
		if (trim) pos = pos[-1]
		x = seriesData(x)
		if (trim) x = x[-1] }

	# Continuous: Calculate Log Returns:
	if (type == "continuous") x = c(NA, diff(log(x)))
	
	# Discrete: Calculate Returns:
	if (type == "discrete") x = c(NA, diff(x))/x
	
	# Percentage Return ?
	if(percentage) x = x*100
	
	# Return as Time Series ?
	if (series == "ts") {
		s = start(x)
		f = frequency(x)
		x = ts(x, start=s, frequency=f) }
		
	# Return an "its" Irregular Time Series ?
	if (series == "its")
    	x = its(x, dates=attributes(x)$dates)
    	
    # Return a "timeSeries" object ?
	if (series == "timeSeries")
    	x = timeSeries(data = x, positions = pos)
	
	# Return Value:
	x
}


# ******************************************************************************


maxDrawDown = 
function(x)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    cmaxx = cummax(x)-x
    mdd = max(cmaxx)
    to = which(mdd == cmaxx)
    from = double(NROW(to))
    for (i in 1:NROW(to))
        from[i] = max(which(cmaxx[1:to[i]] == 0))
    return(list(maxdrawdown = mdd, from = from, to = to))
}


# ------------------------------------------------------------------------------


sharpeRatio = 
function(x, r = 0, scale = sqrt(250))
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(NROW(x) == 1)
        return(NA)
    else {
        y = diff(x)
        return(scale * (mean(y)-r)/sd(y))
    }
}


# ------------------------------------------------------------------------------


sterlingRatio = 
function(x)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(NROW(x) == 1)
        return(NA)
    else {
        return((x[NROW(x)]-x[1]) / maxdrawdown(x)$maxdrawdown)
    }
}


# ------------------------------------------------------------------------------


ohlcPlot = 
function(x, xlim = NULL, ylim = NULL, xlab = "Time", ylab, col = par("col"), 
bg = par("bg"), axes = TRUE, frame.plot = axes, ann = par("ann"), main = NULL,
date = c("calendar", "julian"), format = "%Y-%m-%d",
origin = "1899-12-30", ...)
{
	if ((!is.mts(x)) ||
	  (colnames(x)[1] != "Open") ||
	  (colnames(x)[2] != "High") ||
	  (colnames(x)[3] != "Low") ||
	  (colnames(x)[4] != "Close"))
	  stop("x is not a open/high/low/close time series")
	xlabel = if (!missing(x)) 
	  deparse(substitute(x))
	else NULL
	if (missing(ylab)) 
	  ylab = xlabel
	date = match.arg(date)
	time.x = time(x)
	dt = min(lag(time.x)-time.x)/3
	ylim = c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
	if (is.null(xlim)) 
	  xlim = range(time.x)
	if (is.null(ylim)) 
	  ylim = range(x[is.finite(x)])
	plot.new()
	plot.window(xlim, ylim, ...)
	for (i in 1:NROW(x)) {
	  segments(time.x[i], x[i,"High"], time.x[i], x[i,"Low"],
	           col = col[1], bg = bg)
	  segments(time.x[i] - dt, x[i,"Open"], time.x[i], x[i,"Open"],
	           col = col[1], bg = bg)
	  segments(time.x[i], x[i,"Close"], time.x[i] + dt, x[i,"Close"],
	           col = col[1], bg = bg)
	}
	if (ann) 
	  title(main = main, xlab = xlab, ylab = ylab, ...)  
	if (axes) {
	  if (date == "julian") {
	      axis(1, ...)
	      axis(2, ...)
	  }
	  else {
	      n = NROW(x)
	      lab.ind = round(seq(1, n, length=5))
	      D = as.vector(time.x[lab.ind]*86400) + as.POSIXct(origin, tz = "GMT")
	      DD = format.POSIXct(D, format = format, tz ="GMT")
	      axis(1, at=time.x[lab.ind], lab=DD, ...)
	      axis(2, ...)
	  }
	}
	if (frame.plot) 
	  box(...)
}


# ------------------------------------------------------------------------------

