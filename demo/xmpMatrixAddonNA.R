
#
# Example:
#   Learn to manage missing values in a data matrix. 
#
# Description:
#   Write R functions which remove, substitute, interpolate and
#   impute missing values in a matrix object:
#
#       removeNA        removes NAs from a matrix object
#       subtituteNA     substitutes NAs by zeroes, the column 
#                       mean or column median
#       interpNA        interpolate NAs using R's "approx" 
#                       function
#       knnNA           imputes NAs by the knn-Algorithm using R's
#                       contributed function "knn" from the "EMV" 
#                       Package
# 
# Notes:
#   We didn't take care, that NAs at the border of a matrix or
#   properly considered. That has still to be done!
#   The R functions can be found in "funSeries.R"
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


# Examples:
    
    # DON'T RUN:
    if (FALSE) {


    # Generate a matrix with missing values:
    M = matrix(round(rnorm(40), 2), ncol = 5)
    M[2:3, 1] = M[4, 3:4] = c(NA, NA)
    M[7, 1] =  NA
    colnames(M) = c("a", "b", "c", "d", "e")
    rownames(M) = paste("R", as.character(1:8), sep = "")
    M
    ###
    
    
    # Remove rows with missing values from M:
    M
    removeNA(M)
    # Remove Columns with NAs from M:
    t(removeNA(t(M)))
    ###
    
    
    # Substitute missing values with Zeros, the Median and Mean:
    M
    substituteNA(M)
    substituteNA(M, "median")
    round(substituteNA(M, "mean"), 2)
    ###
    
    
    # Interpolate by Columns:
    M
    interpNA(M, "before")
    interpNA(M, "after")
    interpNA(M, "linear")
    # By Rows:
    M
    t(interpNA(t(M), "before"))
    t(interpNA(t(M), "after"))
    t(interpNA(t(M), "linear"))
    ###
    
    
    # Use "knn" Algorithm:
    knnNA(M, k=2)
    round(knnNA(M, k = 2, correlation = TRUE), 2)
    # Apply it to the transposed matrix:
    t(round(knnNA(t(M), k = 2), 2))
    t(round(knnNA(t(M), k = 2, correlation = TRUE), 2))
    ###
    
    
    } 
    ### END OF DON'T RUN
    

# ******************************************************************************

    