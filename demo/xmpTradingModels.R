
#
# Example: 
#	Trading Models:
#
# Description:
#	Compare two simple trading strategies based on trades with 
#   the MACD Oscillator and on trades in the direction of the trend.
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# 1. Trade with the MACD Oscillator:


	# Read Data:
	data(spc1970)
	data = spc1970
	index = "spc1970"
	###
	
	
	# Work with Logarithmic Closings:
   	O = log(data[,2]) # log Opening Price
	H = log(data[,3]) # log High Price
	L = log(data[,4]) # log Low Price
	C = log(data[,5]) # log Closing Price
	time = 1:length(C)
	date = 1970 + time/252 # Approximate decimal date
	par(mfrow = c(2, 2))
	plot(x = date, y = 100*(C-C[1]), main = paste("Index:", index), type = "l")
	###
	
	
	# Trading Positions and Average Trade Lengths:
 	# position = sign(cdoTA(C, 5, 34, 7))       # MACD Oscillator
	position = (c(0, diff(cdoTA(C, 5, 34, 7)))) # MACD Oscillator
	position = sign(position)
	signals = abs(c(0,diff(position),0))
	tl = emaTA(diff((1:length(signals))[signals>0.5]),2/1261)
	plot(tl, type = "l", xlab = "Number of Trades", 
	  ylab = "Length in Days", main = "Averaged Trade Length")   
	###
	
	
	# Cumulated Return:
	returns = sign(position)*c(diff(C), 0)
   	cumret = 100*cumsum(returns)
	plot(x = date, y = cumret, main = "Cumulative Return", type = "l")
	###
	
	
	# Annualized Returns:
	annualized = 100*252*emaTA(returns, 2/1261)
	plot(x = date, y = annualized, main = "Annualized Returns", type = "l")
	###
	
		
################################################################################
# 2. Trade in the Direction of the Trend


	# Read Data:
	data(spc1970)
	data = spc1970
	index = "spc1970"	
	###
	
	
	# Work with Logarithmic Closings:
 	C = log(data[,5])      # log Opening Price
  	time = 1:length(C)
	date = 1970 + time/252 # Approximate decimal date
	par(mfrow = c(2, 2))
	plot(x = date, y = 100*(C-C[1]), main = paste("Index:", index), type = "l")
	###
	
	
	# Averaged Trade Lengths:
	# Trade for tomorrow in today's direction - go with the trend
 	position = sign(c(0, diff(C)))
	signals = abs(diff(position))
	tl = emaTA(diff((1:length(signals))[signals > 0.5]), 2/1261)
	plot(tl, type = "l", xlab = "Number of Trades", ylab = "Length in Days",
	  main = "Averaged Trade Length")
	###
	
	
	# Cumulated Return:
	returns = sign(position)*c(diff(C), 0)
   	cumret = 100*cumsum(returns)
	plot(x = date, y = cumret, main = "Cumulative Return", type = "l")
	###
	
	
	# Annualized Returns:
	annualized = 100*252*emaTA(returns, 2/1261)
	plot(x = date, y = annualized, 
	  main = "Annualized Returns", type = "l")
	###
	
################################################################################

	