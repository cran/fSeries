
#
# Example: 
#	Perform correlation tests on windowed indicators 
#	to identify good indicators
#  
# Notes: 
#	Results from xmpIndicators1.ssc are needed
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# pattern Tests:


	# Settings:
	# Read file from example xmpIndicatos.ssc
	data(spcindis)
	data = spcindis
	###
	
	
	# Length of Window in Days:
 	win.length = 5*252  # 5 years windows
	win.shift =   2*21  # bi-monthly shifted
	###
	
	
	# Number of Window Cycles:
 	iw = floor((length(data[,1])-win.length)/win.shift) - 1
	statistics = p.values = 
	  matrix(rep(0, times = 6*iw), byrow = TRUE, ncol = 6)
	###
	
	
	# Loop over all Windows:			
	n1 = 1 - win.shift
	n2 = win.length - win.shift
   	for ( i in 1:iw ) {
	  n1 = n1 + win.shift  # Start Window
	  n2 = n2 + win.shift  # End Window
	  for (j in 1:6) {
		result = cor.test(data[,j+1][n1:n2], data[,1][n1:n2], 
		  method = "pearson")
		statistic = as.numeric(result$statistic)
		p.value = as.numeric(result$p.value)
 		statistics[i,j] = statistic
		p.values[i,j] = p.value 
	  }
	cat(i, "out of", iw, ":", n1,n2,"\n") 
	}
	###


	# Plot Result:
	par(mfrow = c(4, 3), cex = 0.5)
	for (i in 1:6) 
	  plot(statistics[,i], type = "l", main = "Statistics")
	for (i in 1:6) 
	  plot(p.values[,i], type = "l", main = "p.value")
	###
	
	
################################################################################

