
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

# Copyright (C) 1998-2003 by Diethelm Wuertz


################################################################################
# FUNCTION:
#  vech				Stack the lower triange of a matrix as vector
#  vec				Stack a matrix as a column vector
#  %x%              Kronecker Product, is part of R's base package
################################################################################


vech = 
function(X)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	vech is the operator that stacks the lower triangle
	#	of a NxN matrix as an N(N+1)/2x1 vector:
	#	vech(X) =(X11, X21, X22, X31, ..., XNN)'
	
	# Note:
	# 	Example for a 3x3 Matrix:
	#	X11, X21, X22, X31, X32, X33
	
	# FUNCTION:
	
	# Return Value:
	t(X[!upper.tri(X)])

}


# ------------------------------------------------------------------------------


vec = 
function(X)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	vec is the operator that stacks a matrix
	#	as a column vector:
	#	vec(X) = (X11, X21, ..., XN1, X12, X22, ..., XNN)'

	# Note:
	# 	Example for a 3x3 Matrix:
	#	X11, X21, X22, X31, X32, X33
	
	# FUNCTION:
	
	# Return Value:
	t(t(as.vector(X)))
}


################################################################################

# 
# Example:
#


# Show, that the following relation holds:
	# vec ( A %*% B %*% C ) = ( t(C) %x% A ) %*% vec(B)
		
# Three random Matrixes:
	A = matrix(rnorm(9), 3)
	B = matrix(rnorm(9), 3)
	C = matrix(rnorm(9), 3)
	
# Print vech:
	A; vech(A)
	
# Print vec:
	B; vec(B)
	
# Vector Product:
	A %*% B

# Kronecker Product:
	A %x% B

# Left Hand Side:
	LHS = vec ( A %*% B %*% C )

# Right Hand Side:
	# %x# denotes the Kronecker Product
	RHS = ( t(C) %x% A ) %*% vec(B)
	
# Test the Difference:
	data.frame(LHS, RHS, LHS-RHS)
