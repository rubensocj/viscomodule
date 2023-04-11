import numpy as np
import math as mt
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve

"""
 Builds the matrix for the linear equation system
 kz - deformation constant rate
 time - time array
 pp - relaxation time array
 num - number of terms in Prony series
 
 Returns - a 2-dimension numpy float array as the matrix from system
"""
def matrixA(kz, time, pp, num):
	nn = num+1
	ma = np.empty(shape=(nn, nn), dtype=float)
	
	# formation law for the first index in the matrix
	# A[1,1]
	ss = 0
	for i in np.arange(0, len(time), 1):
		ss = ss+time[i]**2

	ma[0][0] = kz*ss
	
	# formation law for line 1 and column 1
	# line 1: A[1,n], n = 2, ... , n+1
	# column 1:	A[n,1], n = 2, ... , n+1
	for n in np.arange(0, len(pp), 1): # 1 to 11 in Park
		ss = 0
		for j in np.arange(0, len(time), 1):
			ss = ss+time[j]*(1-mt.exp(-time[j]/pp[n]))
		
		ai = kz*pp[n]*ss
		# column 1:	line n, n = 2, ... , n+1
		ma[n+1][0] = ai
		# line 1:	column n, n = 2, ... , n+1
		ma[0][n+1] = ai
		
	# formation law for n x m sub-matrix where n = 2, ... , n+1 and m = 2, ... , m+1
	# A[n,n], n = 2, ... , n+1
	for n in np.arange(0, len(pp), 1): # row loop
		for m in np.arange(0, len(pp), 1): # column loop
			ss = 0
			for j in np.arange(0, len(time), 1): # summatory loop
				ss = ss+(1-mt.exp(-time[j]/pp[n]))*(1-mt.exp(-time[j]/pp[m]))
			
			# line n, n = 2, ... , n+1
			# column m, m = 2, ... , m+1
			ma[n+1][m+1] = kz*pp[n]*pp[m]*ss
		
	return ma

"""
 Builds the independent vector for the linear equation system
 tension - output tension array
 time - time array
 pp - relaxation time array
 num - number of terms in Prony series
 
 Returns - a 1-dimension numpy float array as the vector from system
"""
def vectorB(tension, time, pp, num):
	nn = num+1
	bb = np.empty(shape=(nn), dtype=float)
	
	# formation law for the first index in the vector
	# b[1]
	ss = 0
	for j in np.arange(0, len(time), 1):
		ss = ss+tension[j]*time[j]
	
	bb[0] = ss
	
	# formation law for n-th indexes, n = 2, ... , n+1
	# b[n], n = 2, ... , n+1
	for n in np.arange(1, len(pp)+1, 1): # 1 to 11 in Park
		ss = 0
		for j in np.arange(0, len(time), 1):
			ss = ss+tension[j]*(1-mt.exp(-time[j]/pp[n-1]))
			
		bb[n] = pp[n-1]*ss
		
	return bb

"""
 Builds the reduced matrix for the linear equation system
 This approach does not determinate the E_inf module
 
 kz - deformation constant rate
 time - time array
 pp - relaxation time array
 num - number of terms in Prony series
 
 Returns - a 2-dimension numpy float array as the matrix from system
"""
def matrixAred(kz, time, pp, num):
	nn = num
	ma = np.empty(shape=(nn, nn), dtype=float)
	
	# formation law for n x n matrix
	# A[i,k]
	for i in np.arange(0, len(pp), 1): # row loop
		for k in np.arange(0, len(pp), 1): # column loop
			ss = 0
			for j in np.arange(0, len(time), 1): # summatory loop
				ss = ss+(1-mt.exp(-time[j]/pp[i]))*(1-mt.exp(-time[j]/pp[k]))
				
			ma[i][k] = kz*pp[i]*pp[k]*ss
		
	return ma

"""
 Builds the reduced independent vector for the linear equation system
 This approach does not determinate the E_inf module
 
 tension - output tension array
 time - time array
 pp - relaxation time array
 num - number of terms in Prony series
 
 Returns - a 1-dimension numpy float array as the vector from system
"""
def vectorBred(tension, time, pp, Einf, kz, num):
	nn = num
	bb = np.empty(shape=(nn), dtype=float)
	
	# formation law for i-th indexes
	# b[i]
	for i in np.arange(0, len(pp), 1): # 1 to 11 in Park
		ss1 = 0
		ss2 = 0
		for j in np.arange(0, len(time), 1):
			ss1 = ss1+tension[j]*(1-mt.exp(-time[j]/pp[i]))
			ss2 = ss2+time[j]*(1-mt.exp(-time[j]/pp[i]))
			
		bb[i] = pp[i]*ss1 - kz*Einf*pp[i]*ss2
		
	return bb

"""
 Applies partial pivoting on the A matrix and the corresponding elimination
 a - input matrix, dtype=float
 b - input independent vector, dtype=float
 
 Returns - a tuple containing: scheduled A matrix aka U matrix, L matrix, permuted b vector and permutation vector
"""
def pivot(a, b):
	aa = a
	bb = b
	nn = aa.shape[0]
	ll = np.zeros(shape=(nn,nn))
	
	# permutation vector:	1, 2, ... , n + 1
	per = list(range(0, nn))
	
	# pivot vector
	# each index corresponds to the maximum row element in modulus
	# mdmax = np.empty(shape=nn, dtype=float)
	
	# fills the pivot vector by finding the maximum modulus values for each row
	# np.arange(ini, end, step)  generates a sequence from ini to end - 1 by step
	# for i in np.arange(0, nn, 1):
		# ss = 0
		# for j in np.arange(0, nn, 1):
			# ss = max(ss, abs(aa[i][j]))
		
		# mdmax[i] = ss

	# i loop throw the columns of matrix aa
	for i in np.arange(0, nn-1, 1):
		id = 0
		rzmax = 0
		for j in np.arange(i, nn, 1):
			# obtain the ratio between the j-th element of column i
			# and maximum value in modulus corresponding to that j-th line
			rz = abs(aa[per.index(j)][i])
			
			# get the smaller index corresponding to bigger ratio 
			# tries if the i-th ratio is bigger than i-1-th
			# if it is, then get the index j
			# the j-th element of perm vector will be changed by the first one
			if rz>rzmax:
				rzmax = rz
				id = j # the index of line where the pivot is
				
		# permutes the vector only if the pivot is in another line
		# permutes the independent vector too
		if id != i:
			per = permut(per, per.index(id), i)
			bb = permut(bb, per.index(id), i)

		lr = 0
		# loop throw the i-th line indexing by the permutation vector
		for k in np.arange(i+1, nn, 1):
			# obtain the element that will fill the L matrix
			lr = aa[per.index(k)][i]/aa[per.index(i)][i]
			
			# make the other elements in that line be zeros
			# loop throw the columns
			for j in np.arange(i, nn, 1):
				aa[per.index(k)][j] = aa[per.index(k)][j] - lr*aa[per.index(i)][j]
			bb[per.index(k)] = bb[per.index(k)] - lr*bb[per.index(i)]
			ll[per.index(k)][i] = lr
		
	# fills the diagonal of L matrix with 1
	# corresponds to sum the (L-I)+I
	for i in np.arange(0, nn, 1):
		for j in np.arange(i, nn, 1):
			ll[per.index(i)][i] = 1
			
	return (aa, ll, bb, per)
		
"""
 Permutes two elements on a vector
 
 per - the vector to be permuted
 j - first position
 i - second position
 
 Returns - the permuted vector
"""
def permut(p, j, i):
	vv = p[j]
	ww = p[i]
	
	p[j] = ww
	p[i] = vv
	
	return p

"""
 Solves a linear equation system by progressive substitution
 In LU decomposing, it is used to perform the Ly = b solving
 
 l - the L matrix obtained from LU decomposing
 b - the permuted independent vector of the Ax = b system
 p - the permutation vector associated to the LU decomposing
 
 Returns - a y permuted vector to be used in Ux = y system
"""
def progressub(l, b, p):
	# Ly = b	: progressive substitution
	ll = l
	bb = b
	per = p
	
	# vector dimension, n
	# initializes empty vector xx to keep the results
	nn = bb.shape[0]
	yy = np.empty(shape=(nn), dtype=float)
	
	yy[0] = bb[per.index(0)]/ll[per.index(0)][0]
	for i in np.arange(1, nn, 1):
		ss = bb[per.index(i)]
		for j in np.arange(0, i, 1):
			ss = ss - ll[per.index(i)][j]*yy[j]
		yy[i] = ss
	
	# reorder the vector by the permutation indexes
	yy = yy[per]
	
	return yy

"""
 Solves a linear equation system by regressive substitution
 In LU decomposing, it is used to perform the Ux = y solving
 
 u - the U matrix obtained from LU decomposing
 y - the resulting vector from progressive substitution
 p - the permutation vector associated to the LU decomposing
 
 Returns - a non permuted x vector as solving of the original Ax = b system
"""
def retrosubst(u, y, p):
	# Ux = y	: regressive substitution
	uu = u
	yy = y
	per = p
	
	# vector dimension, n
	# initializes empty vector xx to keep the results
	nn = yy.shape[0]
	xx = np.empty(shape=(nn), dtype=float)
	
	# last element of vector
	xx[nn-1] = yy[per.index(nn-1)]/uu[per.index(nn-1)][nn-1]
	
	ss = 0
	# np.arange(i, e, step) with step = -1 starts in i and goes to e+1. So this loop must start at nn-1. Then uses nn-2
	for i in np.arange(nn-2, -1, -1):
		ss = yy[per.index(i)]
		for j in np.arange(i+1, nn, 1):
			ss = ss - uu[per.index(i)][j]*xx[j]
		xx[i] = ss/uu[per.index(i)][i]
		
	return xx

"""
 Solves the system. Tests if the matrix A is singular
 aa - input matrix
 bb - input independent vector
 
 Returns - a 1-dimension numpy float array as the result vector from system
"""
def solve(aa, bb):
        if isInvertible(aa):
                return np.linalg.solve(aa, bb)
        else:
                return np.ones((len(bb),), dtype=int)

"""
 Solves the system by LU factorization and obtain the Prony series constants
 aa - input matrix
 bb - input independent vector
 
 Returns - a 1-dimension numpy float array as the result vector from system
"""
def solveSystem(aa, bb):
	# LU method
	# Gauss scheduling with or no pivoting
	# to get the U matrix resulting from scheduling
	# and the scheduling coefficients to builds the L matrix
	# then if Ax = b and A = LU then (LU)x = b
	# doing Ux = y then Ly = b
	# so solving the system:
	# 
	# Ly = b	: progressive substitution
	# Ux = y	: regressive substitution
	
	piv = pivot(aa, bb)
	pro = progressub(piv[1], piv[2], piv[3])
	ret = retrosubst(piv[0], pro, piv[3])
	
	return ret

"""
 Solves the system by LU factorization implemented by scipy
"""
def solveLU(aa, bb):
        lu, piv = lu_factor(aa)
        x = lu_solve((lu, piv), bb)

        return x

"""
 Check if matrix A is invertible
 If a nxn matrix has rank (posto) equals to n, wich means, has no linear dependent lines

 return - boolean
"""
def isInvertible(ma):
        return np.linalg.matrix_rank(ma) > 2 and abs(np.linalg.det(ma)) > 10e-20
    
def jacobi(A, b, x, N):
	D = np.diag(A)
	R = A - np.diagflat(D)

	for i in range(N):
		x = (b - np.dot(R,x))/D
	return x
