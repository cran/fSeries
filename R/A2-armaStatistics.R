
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
# FUNCTION:					DESCRIPTION:
#  armaTrueacf				 True ARMA Autocorrelation Function
#  armaRoots				 Roots of the ARMA Characteristic Polynomial
################################################################################


armaTrueacf = 
function(model, lag.max = 20, type = "correlation", doplot = TRUE)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	A synonyme to ARMAacf

	# Notes:
	#	A synonyme for arma.tacf under R. See R's .First.lib.
	#	Implemented from ARMAacf
	
	# FUNCTION:
	
	# Settings:
	lag = 0:lag.max
	result = NA
	if (type=="partial" || type=="p" || type=="both" || type=="b") {
		main = ylab = "True PACF"
		lag = 1:lag.max
		pacf = ARMAacf(model$ar, model$ma, lag.max=lag.max, 
			pacf=TRUE)
		result = data.frame(cbind(lag, pacf))
		if (doplot) {
			plot(x=lag, y=pacf, type = "n", xlab = "Lag", 
				ylab = ylab, main = main, 
				ylim = c(min(c(pacf, 0)), 1) )
			lines(x = lag, y = pacf, type = "h")
			abline(h = 0)}}
	if (type == "correlation" || type == "c" || type == "both" || type=="b") {
		main = ylab = "True ACF"
		lag = 0:lag.max
		acf = ARMAacf(model$ar, model$ma, lag.max=lag.max, 
			pacf=FALSE)
		result = data.frame(cbind(lag, acf))
		if (doplot) {
			plot(x=lag, y = acf, type = "n", xlab = "Lag", 
				ylab = ylab, main = main, 
				ylim = c(min(c(acf, 0)), 1) )
			lines(x=lag, y=acf, type = "h")
			abline(h = 0) } }	
			
	# Return Value:
	result
}


# ------------------------------------------------------------------------------


armaRoots = 
function(coefficients, n.plot = 400, digits = 4, ...)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	Calculates the roots of a characteristc polynomial

	# FUNCTION:
	
	# Algorithm:
	root = polyroot(c(1, -coefficients))
	real.root = Re(root)
	im.root = Im(root)
	xrange = range(real.root)
	xrange = c(xrange[1] - 1.2*abs(xrange[1]), 
		xrange[2]+1.2 * abs(xrange[2]))
	xplot = seq(xrange[1], xrange[2], length = n.plot)
	fpoly = 1
	for(i in 1:length(coefficients)) {
		fpoly = fpoly - xplot^i * coefficients[i] }
	plot(xplot, fpoly, type = "l", xlab = "B", ylab = "Function", ...)
	title(main = "Polynomial Function vs. B")
	abline(h = 0)
	distance = sqrt(real.root^2 + im.root^2)
	root.mat = cbind(round(real.root, digits = digits),
		round(im.root, digits = digits), 
		round(distance, digits = digits))
	dimnames(root.mat) = list(1:nrow(root.mat), 
		c("re", "im", "dist"))
	size.limit = max(abs(real.root), 1.5, abs(im.root))
	plot(root, xlim = c( - size.limit, size.limit),
		ylim = c( - size.limit, size.limit), 
		xlab = "", ylab = "", ...)
	x = (2*pi/360)*(0:360)
	# symbols(0, 0, circles = 1, add = TRUE, inches = FALSE, col = 6)
	lines(sin(x), cos(x))
	abline(h = 0)
	abline(v = 0)
	title("Roots and Unit Circle", 
		xlab = "Real Part", ylab = "Imaginary Part")
	result = root.mat
		
	# Return Value:
	data.frame(result)
}


# -----------------------------------------------------------------------------

