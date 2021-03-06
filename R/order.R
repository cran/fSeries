
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
# FUNCTION:                REARRANGEMENT:
#  sortColnames             Returns sorted column names of a time Series
#  sampleColnames           Returns sampled column names of a time Series
#  orderColnames            Returns ordered column names of a time Series
#  statsColnames            Returns statistically rearranged column names
#  pcaColnames              Returns PCA correlation ordered column names
#  hclustColnames           Returns hierarchical clustered column names
################################################################################


orderColnames =
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns ordered column names of a time Series
    
    # Arguments:
    
    # FUNCTION:
    
    # Order:
    ans = order(colnames(as.matrix(x)), ...)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


sortColnames = 
function(x, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns sorted column names of a time Series
    
    # Arguments:
    
    # FUNCTION:
    
    # Sort:
    ans = sort(colnames(as.matrix(x)), ...)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


sampleColnames = 
function(x, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns sampled column names of a time Series
    
    # Arguments:
    
    # FUNCTION:
    
    # Sample:
    ans = sample(colnames(as.matrix(x)), ...)
        
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


statsColnames =
function(x, FUN = colMeans, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns statistically rearranged column names
    
    # Arguments:
    
    # Note:
    #   Example of function Candidates:
    #   colStats        calculates column statistics,  
    #   colSums         calculates column sums,  
    #   colMeans        calculates column means,  
    #   colSds          calculates column standard deviations,  
    #   colVars         calculates column variances,  
    #   colSkewness     calculates column skewness,  
    #   colKurtosis     calculates column kurtosis,  
    #   colMaxs         calculates maximum values in each column,  
    #   colMins         calculates minimum values in each column,  
    #   colProds        computes product of all values in each column,  
    #   colQuantiles    computes quantiles of each column.  

    # FUNCTION:
    
    # Apply colStats Function:
    fun = match.fun(FUN)
    Sort = sort(fun(x, ...))
    Order = names(Sort)
    ans = colnames(as.matrix(x)[, Order]) 
    attr(ans, "control") <- Sort    
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pcaColnames =
function(x, robust = FALSE, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns PCA correlation ordered column names
    
    # Arguments:
    
    # FUNCTION:

    # Order:
    if (robust) {
        x.cor = robustbase::covMcd(as.matrix(x), cor = TRUE, ...)$cor
    } else {
        x.cor = cor(as.matrix(x), ...) 
    }
    x.eigen = eigen(x.cor)$vectors[,1:2]
    e1 = x.eigen[, 1]
    e2 = x.eigen[, 2]
    Order = order(ifelse(e1 > 0, atan(e2/e1), atan(e2/e1)+pi)) 
    ans = colnames(as.matrix(x))[Order]
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


hclustColnames =
function(x, method = c("euclidean", "complete"), ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns hierarchical clustered column names
    
    # Arguments:
    #   ...
    #       method - the agglomeration method to be used. This should 
    #           be (an unambiguous abbreviation of) one of "ward", "single", 
    #           "complete", "average", "mcquitty", "median" or "centroid". 
       
    # FUNCTION:
    
    # Order:
    Order = hclust(dist(t(as.matrix(x)), 
        method = method[1]), method = method[2], ...)$order  
    ans = colnames(as.matrix(x))[Order]
       
    # Return Value:
    ans
}


################################################################################

