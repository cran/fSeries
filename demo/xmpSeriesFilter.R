 
# 
# Example:
#   Time Series Filter available under R
#
# Description:
#   Discuss what kind of time series filter are availalbe under R. In
#   addition implement the Hodrick-Prescott filter which decomposes a 
#   macroeconomic time series into a smooth trend component and a 
#   cyclical component. 
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# 1 Filter available in the "base" package:


    # filter - 
    #	Applies linear filtering to a univariate time series or to 
    #	each series separately of a multivariate time series.
    if (FALSE) args(filter)
    ###
    
    
################################################################################
# 2 Filter available in the "pastecs" package:


    #
    # Package: pastecs
    # Title: Package for Analysis of Space-Time Ecological Series
    # Version: 1.1-1
    # Date: 2002-11-12
    # Author: Frederic Ibanez <ibanez@obs-vlfr.fr>, Philippe Grosjean
    #        <phgrosjean@sciviews.org> & Michele Etienne
    #        <etienne@obs-vlfr.fr>
    # Description: Regulation, decomposition and analysis of space-time
    #        series. The pastecs library is a PNEC-Art4 and IFREMER (Benoit
    #        Beliaeff <Benoit.Beliaeff@ifremer.fr>) initiative to translate
    #        PASSTEC 2000
    #        (http://www.obs-vlfr.fr/~enseigne/anado/passtec/passtec.htm).
    # URL: http://www.sciviews.org/pastecs
    # Maintainer: Philippe Grosjean <phgrosjean@sciviews.org>
    # License: GNU Public Licence 2.0 or above at your convenience
    # Depends: boot, ctest, mva, ts
    # Built: R 1.9.0; ; 2004-03-28 14:51:33
    #
    
    #   FILTER:
    #   decaverage    Decomposition using a moving average  
    #   deccensus     Decomposition using the CENSUS II method  
    #   decdiff       Decomposition using differences (trend elimination)  
    #   decevf        Decomposition using eigenvector filtering (EVF)  
    #   decloess      Decomposition by the LOESS method  
    #   decmedian     Decomposition using a running median  
    #   decreg        Decomposition using a regression model  
    
    
################################################################################
# 3 Hodrick-Prescott Filter


	# R function available in "funSeries.R"
	# Compute the cyclical component in the US CPI data for monthly data
    data(CPI.dat)
    logCPI = log(as.vector(CPI.dat[,2])) 
    cycCPI = logCPI - hpFilter(logCPI, 14400)
    par(mfrow = c(2, 1), cex = 0.5)
    ts.plot(logCPI, xlab = "Index", col ="steelblue4"); grid()
    ts.plot(cycCPI, xlab = "Index", col ="steelblue4"); grid()


################################################################################

