
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

# Copyrights (C)
# for this R-port: 
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
# for the code accessed (or partly included) from other R-ports:
#   R: see R's copyright and license file
#   date: Terry Therneau <therneau@mayo.edu>
#     R port by Th. Lumley <thomas@biostat.washington.edu>  K. Halvorsen 
#       <khal@alumni.uv.es>, and Kurt Hornik <Kurt.Hornik@R-project.org>
#   ts: Collected by Brian Ripley. See SOURCES
#   tseries: Compiled by Adrian Trapletti <a.trapletti@bluewin.ch>
# for ical:
#   libical: Libical is an Open Source implementation of the IETF's 
#     iCalendar Calendaring and Scheduling protocols. (RFC 2445, 2446, 
#     and 2447). It parses iCal components and provides a C API for 
#     manipulating the component properties, parameters, and subcomponents.
#   Olsen's VTIMEZONE: These data files are released under the GNU 
#     General Public License, in keeping with the license options of 
#     libical. 
# for the holiday database:
#   holiday information collected from the internet and governmental 
#   sources obtained from a few dozens of websites
    
    
################################################################################
# GENERATION:           DESCRIPTION:
#  matrix               R  creates a matrix from the given set of values
#   diag                R  creates a diagonal matrix or extracts diagonals
#   triang              S  extracs the lower tridiagonal part from a matrix
#   Triang              S  extracs the upper tridiagonal part from a matrix
#   pascal              S  creates a pascal matrix
#   colVec              S  creates a column vector from a vector
#   rowVec              S  creates a row vector from a vector
#  as.matrix            R  attempts to turn its argument into a matrix     
#  is.matrix            R  tests if its argument is a (strict) matrix
#  dimnames             R  retrieves or sets the dimnames of an object
#  colnames|rownames    R  retrieves or sets the row or column names 
#  colIds|rowIds        S  ... use alternatively
#  colIds<-|rowIds<-    S  ... for assignments
# SUBSETS:              DESCRIPTION:
#  dim                  R  returns the dimension of a matrix object
#  ncol|nrow            R  counts columns|rows of a matrix object
#  length               R  counts elements of a matrix object
#   "["|"[["            R  subsets a matrix object
#   (Arith)             R  Elementwise Arithmetic: + - * /
#   (Lops)              R  Elementwise logical Ops: > < >= <= == !=
#  cbind|rbind          R  augments a matrix object by columns|rows
#  na.omit              R  removes NA from a matrix object
# BASIC STATISTICS:     DESCRIPTION:
#  var                  R  returns the variance matrix
#  cov                  R  returns the covariance matrix
#  col|rowStats         B  calculates column|row statistics 
#   col|rowMeans        R  calculates column|row means
#   col|rowAvgs         B  calculates column|row averages
#   col|rowVars         B  calculates column|row variances
#   col|rowStdevs       B  calculates column|row standard deviations
#   col|rowSkewness     B  calculates column|row skewness 
#   col|rowKurtosis     B  calculates column|row kurtosis 
#   col|rowCumsums      B  calculates column|row cumulated sums 
# LINEAR ALGEBRA:       DESCRIPTION:
#  det                  R  returns the determinante of a matrix
#  inv|chol2inv       S|R  returns the inverse of a matrix
#  norm                 S  returns the norm of a matrix
#  rk                   S  returns the rank of a matrix
#  tr                   S  returns trace of a matrix
#  t                    R  returns the transposed matrix
#  %*%                  R  returns the product of two matrices
#  %x%|kron           R|S  returns the Kronecker product
# MORE LINEAR ALGEBRA:  DESCRIPTION:
#  chol                 R  returns the Cholesky factor matrix
#  eigen                R  returns eigenvalues and eigenvectors
#  svd                  R  returns the singular value decomposition
#  kappa                R  returns the condition number of a matrix
#  qr                   R  returns the QR decomposition of a matrix
#  solve                R  solves a system of linear equations
#  backsolve            R ... use when the matrix is upper triangular
#  forwardsolve         R ... use when the matrix is lower triangular
################################################################################
# NOTES:
#  WHERE YOU FIND THE FUCTIONS?
#	R  Basic R Package
#   B  Rmetrics fBasics Package
#   S  Rmetrics fSeries Package
#  REQUIREMENTS:
#	fBasics
#   EMV
################################################################################


triang = 
function(x) 
{   # A function Implemented by Diethelm Wuertz
        
    # Description:
    #   Returns lower triangle matrix
        
    # FUNCTION:
    
    # Triangulate:
    x[row(x) < col(x)] = 0 
        
    # Return Value:
    x 
}
    

# ------------------------------------------------------------------------------

            
Triang = 
function(x) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Returns upper triangle matrix
    
    # FUNCTION:
    
    # Triangulate
    x[row(x) > col(x)] = 0 
    
    # Return Value:
    x 
} 
        

# ------------------------------------------------------------------------------


pascal = 
function(n) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Creates a Pascal matrix
    
    # Arguments:
    #	n - the dimension of the square matrix
    
    # Details:
    #   http://mathworld.wolfram.com/PascalMatrix.html
    #   Pascal matrices are symmetric and positive definite. 
    #   The determinant of a Pascal matrix is 1. 
    #   The inverse of a Pascal matrix has integer entries. 
    #   If lambda is an eigenvalue of a Pascal matrix, 
    #       then 1/lambda is also an eigenvalue of the matrix.
    #   The Cholesky factor of a Pascal matrix consists of 
    #       the elements of Pascal’s triangle
        
    # FUNCTION:
    
    # Pascal:
    N = n-1
    n.over.r = function(n, r) { 
        prod(1:n) / (prod(1:(n-r)) * prod(1:r) ) }
    X = rep(1, N)
    for ( i in 1:N )
        for ( j in 1:N )
        X = c(X, n.over.r(i+j, j))
        X = cbind(rep(1, N+1), matrix(X, byrow = TRUE, ncol = N))
        
    # Return Value:
    X 
}
 

# ------------------------------------------------------------------------------


colVec = 
function(x) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Converts a vector to a column vector
    
    # Details:
    #   A column vector is a matrix with one column.
    
    # Return Value:
    
    # FUNCTION:
    
    # Double Transpose:
    ans = t(t(x)) 
    
    # Return Value:
    ans
}
    

# ------------------------------------------------------------------------------


rowVec = 
function(x) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Converts a vector to a row vector
    
    # Details:
    #   A row vector is a matrix with one row.
    
    # FUNCTION:
    
    # Transpose:
    ans = t(x) 
    
    # Return Value:
    ans
}
    

# ------------------------------------------------------------------------------

       
colIds = 
function(x, ...) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Retrieves row names of a matrix-like object
    
    # FUNCTION:
    
    # Convert to Matrix
    x = as.matrix(x)
    
    # Return Value:
    colnames(x, ...) 
}
        

# ------------------------------------------------------------------------------

        
rowIds = 
function(x, ...) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Retrieves row names of a matrix-like object
    
    # FUNCTION:
    
    # Convert to Matrix
    x = as.matrix(x)
    
    # Return Value:
    rownames(x, ...) }
        

# ------------------------------------------------------------------------------


"colIds<-" = 
function(x, value)
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Sets column names of a matrix-like object
    
    # FUNCTION:
    
    # Column Names:
    dn = dimnames(x)
    if(is.null(dn)) {
        if(is.null(value)) return(x)
        if((nd = length(dim(x))) < 2)
            stop("Object has less than two dimensions")
        dn = vector("list", nd)
    }
    if(length(dn) < 2)
        stop("Object has less than two dimensions")
    if(is.null(value)) dn[2] = list(NULL) else dn[[2]] = value
    dimnames(x) = dn
    
    # Return Value:
    x
}
 

# ------------------------------------------------------------------------------

       
"rowIds<-" = 
function(x, value) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Sets row names of a matrix-like object
    
    # FUNCTION:
    
    # Row names:
    dn = dimnames(x)
    if(is.null(dn)) {
        if(is.null(value)) return(x)
        if((nd = length(dim(x))) < 1)
            stop("attempt to set rownames on object with no dimensions")
        dn = vector("list", nd) }
    if(length(dn) < 1)
        stop("attempt to set rownames on object with no dimensions")
    if(is.null(value)) dn[1] = list(NULL) else dn[[1]] = value
    dimnames(x) = dn
    
    # Return Value:
    x
}

 
# ------------------------------------------------------------------------------

            
inv = 
function(x) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Returns the inverse of a matrix
    
    # FUNCTION:
    
    # Inverse:
    ans = chol2inv(chol(x))
    
    # Return Value:
    ans 
}
        

# ------------------------------------------------------------------------------


norm = 
function(x, p = 2) 
{	# A function Implemented by Diethelm Wuertz
	
	# Description:
	#	Returns the spectral norm of a matrix
	
	# Details:
	#	http://mathworld.wolfram.com/MatrixNorm.html:
	# 	For p = 1
	#		The maximum absolute column sum norm |A|_1 is defined 
	#		as the maximum of the sum of the absolute valued elements
	#		of columns of the matrix.
	# 	For p = 2:
	#		The spectral |A|_2 norm is "the" of a matrix. This value
	#   	is computed as the square root of the maximum eigenvalue   
	#   	of A^H A where A^H is the conjugate transpose.
	#   For p = Inf:
	#		The maximum absolute row sum norm |A|_inf is defined 
	#		as the maximum of the sum of the absolute valued elements
	#		of rows of the matrix.

	# FUNCTION:
	
	# Compute Norm:
	ans = NA
	if (p == 1) {
		x = abs(x)
		ans = max(apply(x, 2, sum)) }
	if (p == 2) {
		ans = sqrt(max(eigen(t(x) %*% x)$values))}
	if (p == Inf) {
		x = abs(x)
		ans = max(apply(x, 1, sum)) }
	if (is.na(ans)) stop("Invalid value for p")
		
	# Return value:
	ans
}


# ------------------------------------------------------------------------------

        
rk = 
function(x, method = c("qr", "chol")) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Returns the rank of a matrix
    
    # FUNCTION:
    
    # Rank:
    method = method[1]
    if (method == "chol") {
        ans = attr(chol(x, pivot = TRUE), "rank") }
    else {
        ans = qr(x)$rank }
    
    # Return Value:
    ans 
}
        

# ------------------------------------------------------------------------------

    
tr = 
function(x) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Returns the trace of a matrix
    
    # FUNCTION:
    
    # Trace:
    if (dim(x)[1] != dim(x)[2] ) {
        return(NA) }
    else {
        return(sum(diag(x))) } }
        
    # Return Value:
    invisible()
            

# ------------------------------------------------------------------------------


kron = 
function(x, y) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Returns Kronecker product
    
    # FUNCTION:
    
    # Kronecker Product:
    ans = x %*% y 
    
    # Return Value:
    ans
}


################################################################################

