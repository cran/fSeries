\name{matrixAddon}

\alias{matrixAddon}

\alias{triang}
\alias{Triang}
\alias{pascal}
\alias{colVec}
\alias{rowVec}
\alias{colIds}
\alias{rowIds}
\alias{colIds<-}
\alias{rowIds<-}

\alias{inv}
\alias{norm}
\alias{rk}
\alias{tr}
\alias{kron}

\title{Matrix Arithmetics and Linear Algebra}


\description{
    
    This help page summarizes selected functions available under R for
    matrix arithmetics and linear algebra. Addititionally functions
    are documented which have been added by Rmetrics.
    
    The origin of the functions is marked in the following way:
    R: part of R's base packages, 
    B: part of Rmetrics' fBasics package,
    S: part of this Rmetrics package, fSeries.
    
    The functions are listed by topic.
    \cr
 
    
    \tabular{rll}{  
    
    		\tab	
			\tab	\bold{General Matrix Functions:} \cr
			
     
    	1   \tab    \code{matrix}   
            \tab    Creates a matrix from the given set of values, \cr
        2   \tab    \code{diag}  
            \tab    R: Creates a diagonal matrix or extracts diagonals, \cr
        3   \tab    \code{triang} 
            \tab    S: Extracs the lower tridiagonal part from a matrix, \cr
        4   \tab    \code{Triang}    
            \tab    S: Extracs the upper tridiagonal part from a matrix, \cr
        5   \tab    \code{pascal}  
            \tab    S: Creates a Pascal matrix, \cr
        6   \tab    \code{colVec}   
            \tab    S: Creates a column vector from a vector, \cr
        7   \tab    \code{rowVec}  
            \tab    S: Creates a row vector from a vector, \cr
        8   \tab    \code{as.matrix}   
            \tab    R: Attempts to turn its argument into a matrix, \cr
        9   \tab    \code{is.matrix}   
            \tab    R: Tests if its argument is a (strict) matrix, \cr
        10  \tab    \code{dimnames}   
            \tab    R: Retrieves or sets the dimnames of an object, \cr
        11  \tab    \code{colnames}, \code{rownames}   
            \tab    R: Retrieves or sets the row or column names, \cr   
        12  \tab    \code{colIds}, \code{rowIds}   
            \tab    R: ... use alternatively, \cr
    
    		\tab	
			\tab	\bold{Simple Matrix Operations:} \cr
			
    	1   \tab    \code{dim}   
            \tab    R: Returns the dimension of a matrix object, \cr
        2   \tab    \code{ncol}, \code{nrow}   
            \tab    R: Counts columns|rows of a matrix object, \cr 
        3   \tab    \code{length}   
            \tab    R: Counts elements of a matrix object, \cr 
        4   \tab    \code{"["}, \code{"[["}   
            \tab    R: Subsets a matrix object, \cr
        5   \tab    \code{cbind}, \code{rbind}   
            \tab    R: Augments a matrix object by columns|rows, \cr   

    		\tab	
			\tab	\bold{Basic Statistical Functions:} \cr
			
    	1   \tab    \code{var}   
            \tab    R: Returns the variance matrix, \cr
        2   \tab    \code{cov}  
            \tab    R: Returns the covariance matrix, \cr 
        3   \tab    \code{colStats}, \code{rowStats} 
            \tab    B: Calculates column|row statistics, \cr
        4   \tab    \code{colMeans}, \code{rowMeans}   
            \tab    B: Calculates column|row means, \cr 
        5   \tab    \code{colAvgs}, \code{rowAvgs}   
            \tab    B: Calculates column|row averages, \cr
        6   \tab    \code{colVars}, \code{rowVars}   
            \tab    B: Calculates column|row variances, \cr 
        7   \tab    \code{colStdevs}, \code{rowStdevs}   
            \tab    B: Calculates column|row standard deviations, \cr
        8   \tab    \code{colSkewness}, \code{rowSkewness}  
            \tab    B: Calculates column|row skewness, \cr
        9   \tab    \code{colKurtosis}, \code{rowKurtosis}   
            \tab    B: Calculates column|row kurtosis, \cr 
        10  \tab    \code{colCumsums}, \code{rowCumsums} 
            \tab    B: Calculates column|row cumulated sums, \cr  
            
			\tab	
			\tab	\bold{Linear Algebra:} \cr
			
    	1   \tab    \code{\%*\%}
            \tab    R: Returns the product of two matrices, \cr
        2   \tab    \code{\%x\%}, \code{kron}
            \tab    R: Returns the Kronecker product, \cr 
        3   \tab    \code{det}   
            \tab    R: Returns the determinante of a matrix, \cr
        4   \tab    \code{inv} 
            \tab    S: Returns the inverse of a matrix, \cr
        5   \tab    \code{norm} 
            \tab    S: Returns the norm of a matrix, \cr 
        6   \tab    \code{rk} 
            \tab    S: Returns the rank of a matrix, \cr
        7   \tab    \code{tr} 
            \tab    S: Returns trace of a matrix, \cr 
        8   \tab    \code{t} 
            \tab    R: Returns the transposed matrix, \cr
	
			\tab	
			\tab	\bold{More Linear Algebra:} \cr
	
    	1   \tab    \code{chol}   
            \tab    R: Returns the Cholesky factor matrix, \cr
        2   \tab    \code{eigen}  
            \tab    R: Returns eigenvalues and eigenvectors, \cr 
        3   \tab    \code{svd}
            \tab    R: Returns the singular value decomposition, \cr
        4   \tab    \code{kappa} 
            \tab    R: Returns the condition number of a matrix, \cr 
        5   \tab    \code{qr}
            \tab    R: Returns the QR decomposition of a matrix, \cr
        6   \tab    \code{solve} 
            \tab    R: Solves a system of linear equations, \cr
        7   \tab    \code{backsolve} 
            \tab    R: ... use when the matrix is upper triangular, \cr
        8   \tab    \code{forwardsolve}
            \tab    R: ... use when the matrix is lower triangular. }  

}


\usage{
triang(x)
Triang(x)
pascal(n)

colVec(x)
rowVec(x)
colIds(x, \dots)
rowIds(x, \dots)

inv(x)
norm(x, p = 2)
rk(x, method = c("qr", "chol"))
tr(x)
kron(x, y)
}


\arguments{

	\item{method}{
	 	[rk] - \cr
	 	a character value, the dimension of the square matrix. One can
	 	choose from two methods: For \code{method = "qr"} the rank is
	 	computed as \code{qr(x)$rank}, or alternatively for 
	 	\code{method = "chol"} the rank is
	 	computed as \code{attr(chol(x, pivot = TRUE), "rank")}.
	 	}
	 \item{n}{
	 	[pascal] - \cr
	 	an integer value, the dimension of the square matrix
	 	}
	\item{p}{
	 	[norm] - \cr
	 	an integer value, \code{1}, \code{2} or \code{Inf}.
	 	\code{p=1} - The maximum absolute column sum norm which is defined 
		as the maximum of the sum of the absolute valued elements of columns 
		of the matrix.
		\code{p=2} - The spectral norm is "the norm" of a matrix \code{X}. 
		This value is computed as the square root of the maximum eigenvalue 
		of \code{CX} where \code{C} is the conjugate transpose.
		\code{p=Inf} - The maximum absolute row sum norm is defined 
		as the maximum of the sum of the absolute valued elements
		of rows of the matrix.
	 	}
	\item{x, y}{
		a numeric matrix.
		}
	\item{\dots}{
		arguments to be passed.
		}
}


\details{

	\bold{Function from R's Base Package:}
	\cr\cr
	Most of the functions are described in their R help pages which we
	recommend to consult for further information. For the additiotnal 
	functions added by Rmetrics we give a brief introduction.
	\cr
	
	
	\bold{General Functions:}
	\cr
	
	Functions to generate matrices and related functions are described 
	in the help page \code{\link{matrix}}. To "decorate" these objects 
	several naming functions are available, a description can be found 
	on the help pages \code{\link{dimnames}} and \code{\link{rownames}}.
	\cr
	
	The function \code{pascal} generates a Pascal matrix of order \code{n} 
	which is a symmetric, positive, definite matrix with integer entries 
	made up from Pascal's triangle. The determinant of a Pascal matrix is 1. 
	The inverse of a Pascal matrix has integer entries. If \code{lambda} 
	is an eigenvalue of a Pascal matrix, then \code{1/lambda} is also an 
	eigenvalue of the matrix. Pascal matrices are ill-conditioned. 
	\cr
	
	The functions \code{triang} and \code{Triang} allow to transform a
	square matrix to a lower or upper triangular form. 
	A triangular matrix is either an upper triangular matrix or lower 
	triangular matrix. For the first case all matrix elements \code{a[i,j]}
	of matrix \code{A} are zero for \code{i>j}, whereas in the second case
	we have just the opposite situation. A lower triangular matrix is 
	sometimes also called left triangular. In fact, triangular matrices 
	are so useful that much computational linear algebra begins with 
	factoring or decomposing a general matrix or matrices into triangular 
	form. Some matrix factorization methods are the Cholesky factorization 
	and the LU-factorization. Even including the factorization step, 
	enough later operations are typically avoided to yield an overall 
	time savings. Triangular matrices have the following properties: the 
	inverse of a triangular matrix is a triangular matrix, the product of 
	two triangular matrices is a triangular matrix, the determinant of a 
	triangular matrix is the product of the diagonal elements, the 
	eigenvalues of a triangular matrix are the diagonal elements.
	\cr
	
	The functions \code{colVec} and \code{rowVec} transform a vector into 
	a column and row vector, respectively. A column vector is a matrix 
	object with one column, and a row vector is a matrix object with one row.
	\cr
	
	The functions \code{dimnames}, \code{colname}, and \code{rowname}
	can be used to retrieve and set the names of matrices. The functions 
	\code{rowIds}, \code{colIds}, are S-Plus like synonyms.
	\cr
	
	
	\bold{Simple Matrix Operations:}
	\cr
	
	The functions \code{\link{dim}}, \code{\link{nrow}} and 
	\code{\link{ncol}} are functions to extract the dimension and  
	the number of rows and columns of a matrix. 
	\cr
	
	The usual arithmetic operators, logical operators and mathematical 
	functions like \code{sqrt} or \code{exp} and \code{log}
	operate on matrices element by element. Note, that \code{"*"}
	is not matrix multiplication, instead we have to use \code{"\%*\%"}.
	\cr
	
	The methods \code{"["} and \code{"[["} are suited to extract
	subsets from a matrix, to delete rows and columns, or to permute
	rows and columns. 
	\cr
	
	
	\bold{Basic Statistical Functions:}
	\cr
	
	The functions \code{var} and \code{cov} compute the variance and
	covariance of a matrix.
	\cr
	
	Beside these functions \code{\link{colMeans}} and \code{\link{rowMeans}} 
	are R functions which compute the mean of columns and rows of a matrix. 
	Rmetrics has added further functions to compute column- or rowwise 
	variances, standard deviations, skewness, kurtosis and cumulated sums. 
	Two general  functions named \code{rowStats} and 
	\code{colStats} allow to apply through the argument list any 
	function to compute row and column statistics from matrices.
	\cr
	
	
	\code{Linear Algebra:}
	\cr
	
	Matrix multiplication is done using the operator \code{\%*\%}. 
	\cr
	
	The \emph{Kronecker product} can be computed using the operator 
	\code{\%x\%} or alternatively using the function \code{kron}.
	\cr

	The function \code{det} computes the determinant of a matrix.
	\cr
	
	The inverse of a square matrix \code{inv(X)} of dimension
	\code{n} is defined so that 
	\code{X \%*\% inv(X) = inv(X) \%*\% X = diag(n)} where
	the matrix \code{diag(n)} is the \code{n}-dimensional identity matrix.
	A precondition for the existence of the matrix inverse is that the
	determinant of the matrix \code{det(X)} is nonzero. 
	\cr

	The function \code{t} computes the transposed of a square matrix.
	\cr
	
	The function \code{tr} computes the trace of a square matrix which
	is the sum of the diagonal elements of the matrix under consideration.
	\cr
	
	The function \code{rk} computes the rank of a matrix which is 
	the dimension of the range of the matrix corresponding to the number 
	of linearly independent rows or columns of the matrix, or to the 
	number of nonzero singular values. The rank of a matrix is also
	named inear map. 
	\cr
	
	The function \code{norm} computes the norm of a matrix. Three choices 
	are possible: 
	\code{p=1} - The maximum absolute column sum norm which is defined 
	as the maximum of the sum of the absolute valued elements of columns 
	of the matrix.
	\code{p=2} - The spectral norm is "the norm" of a matrix \code{X}. 
	This value is computed as the square root of the maximum eigenvalue 
	of \code{CX} where \code{C} is the conjugate transpose.
	\code{p=Inf} - The maximum absolute row sum norm is defined 
	as the maximum of the sum of the absolute valued elements
	of rows of the matrix.
	\cr
	

	\code{More Linear Algebra:}
	\cr
	
	The function \code{chol} returns the Cholesky factor matrix,
	\code{eigen} returns eigenvalues and eigenvectors,
	\code{svd} returns the singular value decomposition,
	\code{kappa} estimate the condition number of a matrix,
	\code{qr} returns the QR decomposition of a matrix,
	\code{ginv} returns the Moore-Penrose generalized inverse,
	\code{solve} solves a system of linear equations, use
	\code{backsolve} when the matrix is upper triangular, and
	\code{forwardsolve} when the matrix is lower triangular.
	 
}


\references{

Higham N.J., (2002);
	\emph{Accuracy and Stability of Numerical Algorithms}, 
	2nd ed., SIAM.

Golub, van Loan, (1996);
	\emph{Matrix Computations}, 
	3rd edition. Johns Hopkins University Press. 
}


\examples{
## Examples:

## Create Pascal Matrix:
   xmpSeries("\nStart: Pascal Matrix > ")
   P = pascal(3)
   P
   # Create lower triangle matrix
   L = triang(P)
   L
   # Extract diagonal part
   diag(P)
   
## Add/Subtract/Multiply/Divide:  
   xmpSeries("\nNext: Math Operations > ")
   X = P
   # Multiply matrix with a constant
   3 * X
   # Multiply two matrices elementwise
   X * P                     
   # Multiplies rows/columns of a matrix by a vector
   X \%*\% diag(P)            
   diag(P) \%*\% X           
       
## Operate on Subsets of a Matrix:
   xmpSeries("\nNext: Matrix Subsets > ")
   n = 3; i = 2; j = 3
   D = diag(1:3)
   # Return the dimension of a matrix
   dim(P)                         
   # Get the last colum of a matrix
   P[, ncol(P)]                   
   # Delete a column of a matrix
   P[, -i]                      
   # Permute the columns of a matrix
   P[c(3, 1, 2), ]              
   # Augments matrix horizontally 
   cbind(P, D)                           
      
## Apply a function to all Elements of a Matrix: 
   xmpSeries("\nNext: Operate Element by Element > ")
   # Return square root for each element
   sqrt(P)
   # Exponentiate the matrix elementwise
   exp(P)
   # Compute the median of each column
   apply(P, 2, "median") 
   # Test on all elements of a matrix       
   all( P > 2 )   
   # test on any element in a matrix                
   any( P > 2 )                  
     
## More Matrix Operations:
   xmpSeries("\nNext: More Matrix Operations > ")
   # Return the product of two matrices
   P \%*\% D   
   # Return the Kronecker Product                     
   P \%x\% D                        
   # Return the transposed matrix
   t(P)                           
   # Return the inverse of a matrix
   inv(P)  
   # Return the norm of a matrix                      
   norm(P)                        
   # Return the determinante of a matrix
   det(P)                         
   # Return the rank of a matrix
   rk(P)                            
   # Return trace of a matrix
   tr(P)                          
   # Return the variance matrix
   var(P)     
   # Return the covariance matrix                   
   cov(P)                         
    
## More Linear Algebra:
   xmpSeries("\nNext: Linear Algebra > ")
   X = P; b = c(1, 2, 3)
   # Return the Cholesky factor matrix
   chol(X)                        
   # Return eigenvalues and eigenvectors
   eigen(X)                       
   # Return the singular value decomposition
   svd(X)                         
   # Return the condition number of a matrix
   kappa(X)                       
   # Return the QR decomposition of a matrix
   qr(X)                          
   # Solve a system of linear equations
   # ... use backsolve when the matrix is upper triangular
   # ... use forwardsolve when the matrix is lower triangular
   solve(X, b)  
   backsolve(Triang(X), b)
   solve(Triang(X), b)                 
   forwardsolve(triang(X), b) 
   solve(triang(X), b)                                        
}


\author{

    Diethelm Wuertz for the addons of this R-port.
    
}


\keyword{math}

