
#
# Example: 
#   Create Trading Patterns
#
# Description:
# 	This example shows how to calculate, to diplay and to write data 
#   to a file with selected indicators including responses and predictors
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# Create Pattern:

	# Settings:
	data(spc1970)
    dosave = FALSE
    file = "spc1970-pattern.csv"
    ###

    
	# Work with Logarithmic Data:
  	H = log(spc1970[, 3])  	# log High
  	L = log(spc1970[, 4])  	# log Low
  	C = log(spc1970[, 5])  	# log Close
  	R = c(0, diff(C, 1))    # log Return
  	###
  	

	# Plot Closing Prices and Returns:
  	par(mfcol = c(2, 2))
  	plot(C, main = "Log(Close)", type = "l")
  	plot(R, main = "Returns", type = "l")
  	###
  	
 
	# Select a Piece of the Time Series to Display:
  	n1 = 5900; n2 = 6025
  	plot(x = n1:n2, y = C[n1:n2], main = "Window - Log(Close)", type = "l")
  	plot(x = n1:n2, y = R[n1:n2], main = "Window - Returns", type = "l")
	###
	
	
	# Responses:
  	# Tomorrows returns - Shift One Day Back
  	response = c(diff(C), 0) 
 	###
 	
 	
	# Predictors: 
  	# Calculate and Plot(Window) some Selected Indicators:
  	par(mfrow = c(3, 2))
  	p01 = fpkTA(C, H, L, 12)
	  plot(p01[n1:n2], main = "Predictor: %K[12]", type = "l")
  	p02 = fpkTA(C, H, L, 12) - fpdTA(C, H, L, 12, 3)
	  plot(p02[n1:n2], main = "Predictor: %K[12]-%D[12, 3]", type = "l")
  	p03 = rsiTA(C, 6)-rsiTA(C, 12)
	  plot(p03[n1:n2], main = "Predictor: RSI[6]-RSI[12]", type = "l")
  	p04 = oscTA(C,3,10)
	  plot(p04[n1:n2], main = "Predictor: OSC[C,3,6]", type = "l")
  	p05 = cdoTA(C,11, 26, 9)
	  plot(p05[n1:n2], main = "Predictor: CDO[C,12,26,9]", type = "l")
  	p06 = wprTA(C, H, L, 5)
	  plot(p06[n1:n2], main = "Predictor: WPR[C,H,L,5]", type = "l")
	###
	
	
	# Save Responses and Predictors Pattern:
	if (dosave) {
	  z = cbind.data.frame(response, p01, p02, p03, p04, p05, p06)
      names(z) = c("R[NYSE|-1]", "%K[12]", "%K[12]-%D[12|3]",
        "RSI[6]-RSI[12]", "OSC[C|3|6]", "CDO[C|12|26|9]", "WPR[C|H|L|5]")
	  write.table(z, file, sep = ",", dimnames.write = "colnames") 
	}
    ###	

    
################################################################################

		