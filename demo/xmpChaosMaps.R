
#
# Example:
#   Simulate chaotic time series models
#
# Description:
#   Write R functions which generate the following chaotic time
#   series models:
#
#       henonSim        Simulates data from Henon Map 
#       ikedaSim        Simulates data from Ikeda Map 
#       logisticSim     Simulates data from Logistic Map
#       lorentzSim      Simulates data from Lorentz Map
#       roesslerSim     Simulates data from Roessler Map
#   
# Notes:
#   To genererate the series for the Lorentz and Roessler maps
#   use function "rk4" R's contributed  package "odesolve" 
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


#
# Package Information:
#
#   Package: odesolve
#   Version: 0.5-10
#   Date: 2004/04/05
#   Title: Solvers for Ordinary Differential Equations
#   Author: R. Woodrow Setzer <setzer.woodrow@epa.gov>
#   Maintainer: R. Woodrow Setzer <setzer.woodrow@epa.gov>
#   Depends: R (>= 1.4.0)
#   Description: This package provides an interface for the ODE solver lsoda.  
#       ODEs are expressed as R functions or as compiled code.    
#   License: GPL version 2 
#


# ------------------------------------------------------------------------------


# Examples:

	# The functions can be found in "fSeries.R"

	# Plot Henon Map:
	x = henonMap(5000)
	plot(x, col = "steelblue4", main = "Henon Map")


# ------------------------------------------------------------------------------

