
#
# Example:
#   Write an interface to the Garch Ox Software Package.
#
# Details:
#   Currently, as of May 2004, I am writing a new R package for 
#   the analysis of GARCH processes. As tests I use the benchmarks
#   described and discussed by Brooks, Burke and Persand, and by
#   McCullough and Renfro. As a reference for my R implementation 
#   I used among other software products mainly the GARCH Ox Package 
#   of ...
#   To make the comparisons more easier for me I have implementd
#   an Interface for R to Garch Ox. Here you will find it. The 
#   concept is quite simple, the Garch function writes all 
#   parameters to an ASCII file, then the Ox engine is started, 
#   reads the parameters from the file and returns the results 
#   to other ASCII files. These again are read by R and can be 
#   used for the diagnostic and graphical analysis. Only a print 
#   and a plot method are available.
#   
# Notes:
#   The function "garchOxFit" interfaces a subset of the functionality 
#   of the G@ARCH 3.0 Package written in Ox. G@RCH 3.0 is to my 
#   opinion one of the most sophisticated packages for modelling 
#   univariate GARCH processes including GARCH, EGARCH, GJR, APARCH, 
#   IGARCH, FIGARCH, FIEGARCH, FIAPARCH and HYGARCH models. Parameters
#   can be estimated by Approximate (Quasi-) Maximum Likelihood Methods
#   under four assumptions: normal, Student-t, GED or skewed Student-t .
#
#   "Ox" is an object-oriented matrix language with a comprehensive 
#   mathematical and statistical function library. Many packages were 
#   written for Ox including software mainly for econometric modelling. 
#   The Ox packages for time series analysis and forecasting, Arfima,
#   Garch and State Space Modelling are especially worth to note. 
#
#   Before you can use Ox, you have to check the "Ox citation and 
#   "copyright" rules and if you agree and fullfill the conditions, 
#   then download the OxConsole Software together with the "OxGarch" 
#   Package, currently G@RCH 3.0. If you are not qualified for a free 
#   license, order your copy from Timberlake Consultants. 
#
#   Windows: I recommend to install the "Setup.exe" under the path 
#   "C:\\Ox\\" and to unzip the OxGarch Package in the directory 
#   "C:\\Ox\\Packages".
#   Linux: not used so far ...
#
#   Ox Citation and Copyright Rules: "Ox" and all its components 
#   are copyright of Jurgen A. Doornik. The Console (command line) 
#   versions may be used freely for academic research and teaching 
#   purposes only. Commercial users and others who do not qualify 
#   for the free version must purchase the Windows version of Ox 
#   and GiveWin with documentation, regardless of which version 
#   they use (so even when only using "Ox" on Linux or Unix). 
#   Ox must be cited whenever it is used. Refer to the references 
#   given below. Note, failure to cite the use of "Ox" in published 
#   work may result in loss of the right to use the free version, 
#   and an invoice at the full commercial price. The "Ox" syntax is 
#   public, and you may do with your own "Ox" code whatever you wish.
#   
#   Remember, only a small part of the functionalities are interfaced 
#   to R. But, principally it would be possible to interface also other
#   functionalities offered by the "Ox" Garch Package. The "Ox" library 
#   file "libsOxGarch.ox". Feel free to modify and to add additional 
#   functionalities.
#
# Notes:
#   OX and GARCH@OX ARE NOT PART OF THIS DISTRIBUTION!
#
# References:
#   Brooks C., Burke S.P, Persand G. (2001);
#       \emph{Benchmarks and the Accuracy of GARCH Model Estimation},
#       International Journal of Forecasting 17, 45--56.    
#   Doornik, J.A. (2002), 
#       Object-Oriented Matrix Programming Using Ox, 
#       London, 3rd ed.: Timberlake Consultants Press and Oxford: 
#       www.nuff.ox.ac.uk/Users/Doornik.    
#   Doornik, J.A. and Ooms, M. (1999), 
#       A Package for Estimating, Forecasting and Simulating Arfima Models, 
#       Oxford: www.nuff.ox.ac.uk/Users/Doornik.        
#   Laurent S.,PETERS J.P. (2002);
#       \emph{G@RCH 2.2: an Ox Package for Estimating and Forecasting 
#       Various ARCH Models}, 
#       Journal of Economic Surveys, 16, 447--485.      
#   McCullough B.D., Renfro C.G. (1998);
#       \emph{Benchmarks and Software Standards: A Case Study of GARCH 
#       Procedures},
#       Journal of Economic and Social Measurement 25, 59--71. 
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


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
# MA 02111-1307 USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  garchOxFit                Fits parameters of a garch model           
#   print.garchOx            S3 Print Method
#   plot.garchOx             S3 Plot Method
################################################################################


################################################################################
# garchOx
# The functions can be found in "funSeries.R"


# Windows Only, adapt it under Linux ...

    OXBIN <<- "C:/Ox/bin"
    OXCMD <<- "C:\\Ox\\bin\\oxl.exe library\\fSeries\\libs\\GarchOx.ox"
    
    # Creates Files: "OxParameter.txt", "OxSeries.csv", "OxSeries.csv"


################################################################################


## Examples
    
    # Description:
    #   The file "dem2gbp" contains daily observations of the 
    #   Deutschmark / British Pound foreign exchange log returns. 
    #   This data set has been promoted as an informal benchmark 
    #   for GARCH time-series software validation. See McCullough and 
    #   Renfro [1991], and Brooks, Burke, and Persand [2001] for details.
    #   The nominal returns are expressed in percent, as published in 
    #   Bollerslev and Ghysels [2001]. The data set is available from 
    #   the \emph{Journal of Business and Economic Statistics}, (JBES), 
    #   \emph{ftp://www.amstat.org}. A text file has one column of 
    #   data listing the percentual log-returns of the DEM/GBP exchange 
    #   rates. The sample period is from January 3, 1984, to December 
    #   31, 1991, for a total of 1975 daily observations of FX exchange 
    #   rates.

    # Load Benchmark Data Set:
    data(dem2gbp)
    x = dem2gbp[, 1]

    
## Example 1: ARMA/GARCH Models -  Gaussian Distribution
    
    arch2 = garchOxFit(
        formula.mean = ~ arma(0,0), formula.var = ~ garch(0,2))
        arch2
        
    garch20 = garchOxFit(
        formula.mean = ~ arma(0,0), formula.var = ~garch(2,0))
        garch20
    
    garch11 = garchOxFit(
        formula.mean = ~arma(0,0), formula.var = ~garch(1,1))
        garch11
        
    ar1.garch11 = garchOxFit(
        formula.mean = ~arma(1,0), formula.var = ~garch(1,1))   
        ar1.garch11     
        
    ma1.garch11 = garchOxFit(
        formula.mean = ~arma(0,1), formula.var = ~garch(1,1))    
        ma1.garch11
        
    arma11.garch11 = garchOxFit(
        formula.mean = ~arma(1,1), formula.var = ~garch(1,1))
        arma11.garch11
  
        
## Example 2: Other than Gaussian Distributions:
  
    
    garch11.t = garchOxFit(
        formula.var = ~garch(1,1), cond.dist = "t")  
        garch11.t   
        
    garch11.ged = garchOxFit(
        formula.var = ~garch(1,1), cond.dist = "ged")   
        garch11.ged     
            
    garch11.st = garchOxFit(
        formula.var = ~garch(1,1), cond.dist = "skewed-t")  
        garch11.st

     
## Example 3: Extended GARCH Models:
    
    # Sorry, print method is not yet working for the following models:
    
    egarch11 = garchOxFit(
        formula.var = ~ egarch(1,1))    
    
    gjr11 = garchOxFit(
        formula.var = ~ gjr(1,1))   
    
    aparch11 = garchOxFit(
        formula.var = ~ aparch(1,1))  
    
    # Fails ...
    # igarch11 = garchOxFit(
    #   formula.var = ~ igarch(1,1)) 
    
        
# ------------------------------------------------------------------------------

    