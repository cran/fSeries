# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA. 

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
# FUNCTION:       DESCRIPTION:
#  henonSim        Simulate data from Henon Map 
#  ikedaSim        Simulate data from Ikeda Map
#  logisticSim     Simulate data from Logistic Map
#  lorentzSim      Simulate data from Lorentz Map
#  roesslerSim     Simulate data from Roessler Map
#  .rk4            Internal Funtion - Runge-Kutta Solver
################################################################################


henonSim = 
function(n = 1000, n.skip = 100, parms = c(a = 1.4, b = 0.3), 
start = runif(2), doplot = FALSE)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Simulate Data from Henon Map
    
    # Details:
    #   Creates iterates of the Henon map:
    #   *   x(n+1)  =  1 - a*x(n)^2 + b*y(n)
    #   *   y(n+1)  =  x(n)

    # Arguments:
    #   n - number of points x, y
    #   n.skip - number of transients discarded
    #   a - parameter a
    #   b - parameter b
    #   start[1] - initial x
    #   start[2] - initial y
    
    # FUNCTION:
    
    # Simulate Map: 
    a = parms[1]
    b = parms[2]
    x = rep(0, times = (n+n.skip))
    y = rep(0, times = (n+n.skip))
    x[1] = start[1]
    y[1] = start[2]
    for ( i in 2:(n+n.skip) ) {
        x[i]  =  1 - a*x[i-1]^2 + b*y[i-1]
        y[i]  =  x[i-1] }
    x = x[(n.skip+1):(n.skip+n)] 
    y = y[(n.skip+1):(n.skip+n)] 

    # Plot Map:
    if (doplot) {
        plot(x = x, y = y, type = "n", xlab = "x[n]", ylab = "y[n]", 
        main = "Henon Map")
        points(x = x, y = y, col = "red", cex = 0.25) }
    
    # Return Value:
    data.frame(x, y)    
}


# ------------------------------------------------------------------------------


ikedaSim = 
function(n = 1000, n.skip = 100, parms = c(a = 0.4, b = 6.0, c = 0.9), 
start = runif(2), doplot = FALSE)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Simulate Ikeda Map Data
    
    # Details:
    #   Prints iterates of the Ikeda map (Re(z) and Im(z)): 
    #                                        i*b
    #   z(n+1)  =  1 + c*z(n)* exp( i*a - ------------ )
    #                                     1 + |z(n)|^2

    # Arguments:
    #   n - number of points z
    #   n.skip - number of transients discarded
    #   a - parameter a
    #   b - parameter b; 6.0 
    #   c - parameter c; 0.9 
    #   start[1] - initial Re(z)
    #   start[2] - initial Im(z)
    
    # FUNCTION:
    
    # Simulate Map: 
    a = parms[1]
    b = parms[2]
    c = parms[3]
    a = complex(real = 0, imag = a)
    b = complex(real = 0, imag = b)
    z = rep(complex(real = start[1], imag = start[2]), times = (n+n.skip))
    for ( i in 2:(n+n.skip) ) {
        z[i] = 1 + c*z[i-1] * exp(a-b/(1+abs(z[i-1])^2)) }
    z = z[(n.skip+1):(n.skip+n)] 
    
    # Plot Map:
    if (doplot) {
        x = Re(z); y = Im(z)
        plot(x, y, type = "n", xlab = "x[n]", ylab = "y[n]", 
            main = "Ikeda Map")
        points(x, y, col = "red", cex = 0.25)
        x = Re(z)[1:(length(z)-1)]; y = Re(z)[2:length(z)]
        plot(x, y, type = "n", xlab = "x[n]", ylab = "x[n+1]", 
            main = "Ikeda Map")
        points(x, y, col = "red", cex = 0.25) }
    
    # Return Value:
    data.frame(Re(z), Im(z))   
}
 

# ------------------------------------------------------------------------------


logisticSim = 
function(n = 1000, n.skip = 100, parms = c(r = 4), start = runif(1), 
doplot = FALSE)
{   # A function written by Diethelm Wuertz
    
    # Description:
    #   Simulate Data from Logistic Map
    
    # Details:
    #   Creates iterates of the Logistic Map:
    #   *   x(n+1)  =  r * x[n] * ( 1 - x[n] )

    # Arguments:
    #   n - number of points x, y
    #   n.skip - number of transients discarded
    #   r - parameter r
    #   start - initial x
    
    # FUNCTION:
    
    # Simulate Map: 
    r = parms[1]
    x = rep(0, times = (n+n.skip))
    x[1] = start
    for ( i in 2:(n+n.skip) ) {
        x[i]  =  r * x[i-1] * ( 1 - x[i-1] ) }
    x = x[(n.skip+1):(n.skip+n)] 

    # Plot Map:
    if (doplot) {
        plot(x = x[1:(n-1)], y = x[2:n], type = "n", xlab = "x[n-1]", 
            ylab = "x[n]", main = "Logistic Map")
        points(x = x[1:(n-1)], y = x[2:n], col = "red", cex = 0.25) }
    
    # Return Value:
    data.frame(x)   
}
    
                       
# ------------------------------------------------------------------------------


lorentzSim = 
function(times = seq(0, 40, by = 0.01), parms = c(sigma = 16, r = 45.92, b = 4),
start = c(-14, -13, 47), doplot = TRUE, ...)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Simulates a Lorentz Map
    
    # Notes:
    #   Requires rk4 from R package "odesolve"
   
    # FUNCTION:
    
    # Requirements:
    # BUILTIN - require(odesolve)
    
    # Attractor:
    lorentz = function(t, x, parms) {
        X = x[1]; Y = x[2]; Z = x[3] 
        with(as.list(parms), {
            dX = sigma * ( Y - X )
            dY = -X*Z + r*X - Y
            dZ = X*Y - b*Z
            list(c(dX, dY, dZ))}) }

    # Classical RK4 with fixed time step:
    s = .rk4(start, times, lorentz, parms)

    # Display:
    if (doplot) {
        xylab = c("x", "y", "z", "x")
        for (i in 2:4) plot(s[, 1], s[, i], type = "l", 
            xlab = "t", ylab = xylab[i-1], main =  "Lorentz", ...)
        k = c(3, 4, 2)
        for (i in 2:4) plot(s[, i], s[, k[i-1]], type = "l", 
            xlab = xylab[i-1], ylab = xylab[i], main =  "Lorentz", ...) }
        
    # Return Value:
    s
}


# ------------------------------------------------------------------------------


roesslerSim = 
function(times = seq(0, 100, by = 0.01), parms = c(a = 0.2, b = 0.2, c = 8.0),
start = c(-1.894, -9.920, 0.0250), doplot = TRUE, ...)
{   # A function written by Diethelm Wuertz
    
    # Description:
    #   Simulates a Lorentz Map
    
    # Notes:
    #   Requires contributed R package "odesolve"
   
    # FUNCTION:
    
    # Attractor:
    roessler = function(t, x, parms) {
        X = x[1]; Y = x[2]; Z = x[3] 
        with(as.list(parms), {
            dX = -(Y+Z)
            dY = X + a*Y
            dZ = b + X*Z -c*Z
            list(c(dX, dY, dZ))}) }

    # Classical RK4 with fixed time step:
    s = .rk4(start, times, roessler, parms)

    # Display:
    if (doplot) {
        xylab = c("x", "y", "z", "x")
        for (i in 2:4) plot(s[, 1], s[, i], type = "l", 
            xlab = "t", ylab = xylab[i-1], main =  "Roessler", ...)
        k = c(3, 4, 2)
        for (i in 2:4) plot(s[, i], s[, k[i-1]], type = "l", 
            xlab = xylab[i-1], ylab = xylab[i], main =  "Roessler", ...) }
        
    # Return Value:
    s
}


# ------------------------------------------------------------------------------


.rk4 =  
function(y, times, func, parms) 
{
    # Description:
    #	Classical Runge-Kutta-fixed-step-integration
	
    # Autrhor:
    # 	R-Implementation by Th. Petzoldt,
	
	# Notes:
	#	From Package: odesolve
	#	Version: 0.5-12
	#	Date: 2004/10/25
	#	Title: Solvers for Ordinary Differential Equations
	#	Author: R. Woodrow Setzer <setzer.woodrow@epa.gov>
	#	Maintainer: R. Woodrow Setzer <setzer.woodrow@epa.gov>
	#	Depends: R (>= 1.4.0)   
	#	License: GPL version 2 
	#	Packaged: Mon Oct 25 14:59:00 2004
	
	# FUNCTION:

	# Checks:
	if (!is.numeric(y)) stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")
    if (!is.numeric(parms)) stop("`parms' must be numeric")

    # Dimension:
    n = length(y)

    # Call func once to figure out whether and how many "global"
    # results it wants to return and some other safety checks
    rho = environment(func)
    tmp = eval(func(times[1],y,parms),rho)
    if (!is.list(tmp)) stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
		stop(paste("The number of derivatives returned by func() (",
			length(tmp[[1]]),
             "must equal the length of the initial conditions vector (",
             length(y),")",sep=""))
	Nglobal = if (length(tmp) > 1) length(tmp[[2]]) else 0
     
    y0 = y
    out = c(times[1], y0)
    for (i in 1:(length(times)-1)) {
        t  = times[i]
        dt = times[i+1] - times[i]
        F1 = dt * func(t,      y0,            parms)[[1]]
        F2 = dt * func(t+dt/2, y0 + 0.5 * F1, parms)[[1]]
        F3 = dt * func(t+dt/2, y0 + 0.5 * F2, parms)[[1]]
        F4 = dt * func(t+dt  , y0 + F3,       parms)[[1]]
        dy = (F1 + 2 * F2 + 2 * F3 + F4)/6
        y1 = y0 + dy
        out<- rbind(out, c(times[i+1], y1))
        y0 = y1
    }

    nm = c("time", 
    	if (!is.null(attr(y, "names"))) names(y) 
    	else as.character(1:n))
    if (Nglobal > 0) {
        out2 = matrix(nrow=nrow(out), ncol=Nglobal)
        for (i in 1:nrow(out2))
            out2[i,] = func(out[i,1], out[i,-1], parms)[[2]]
        out = cbind(out, out2)
        nm = c(nm,
			if (!is.null(attr(tmp[[2]],"names"))) names(tmp[[2]])
			else as.character((n+1) : (n + Nglobal)))
    }
    dimnames(out) = list(NULL, nm)
    
    # Return Value:
    out
}


################################################################################

