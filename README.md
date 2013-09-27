Polynomial Numerical Linear Algebra package for Matlab/Octave
-------------------------------------------------------------

Author: Kim Batselier
------

Multi-functional package for solving all kinds of problems with multivariate polynomials in double precision. 


1. Basics
---------
Only 1 monomial ordering is supported: the degree negative lex ordering. This graded ordering is defined for two n-variate monomials x^a and x^b as

a > b

if 

|a| = \sum_i^n a_i > |b| = \sum_i^n b_i or

|a|=|b| and the leftmost nonzero entry of a-b is negative.

In this definition x^a stands for

x^a = x_1^a_1 * x_2^a_2 * ... * x_n^a_n,

where a is an n-tuple of nonnegative integers.

A system of s multivariate polynomials is represented by an s-by-2 Any Array. The first column elements are vectors containing the coefficients. The second column elements are matrices containing the corresponding exponents. For example, the system

f_1 = 5.3 x_1^2 + 9 x_2 x_3 -1
f_2 = 2 x_1^3 + .5 x_2^2 - 7.89 x_3 - 94
f_3 = x_1 - 2.13

is represented by

polysys{1,1}=[5.3,9,-1]
polysys{1,2}=[2 0 0;0 1 1;0 0 0]
polysys{2,1}=[2,.5,-7.89,-94]
polysys{2,2}=[3 0 0;0 2 0;0 0 1;0 0 0]
polysys{3,1}=[1,-2.13]
polysys{3,2}=[1 0 0;0 0 0]

2. Available functions
----------------------

index = feti(exponent)

Converts an exponent to its corresponding linear index with respect to the degree negative lex ordering. If the exponent argument is a matrix of exponents, then feti() is applied to each row of the matrix.

exponent = fite(index,n)

Converts a linear index with respect to the degree negative lex ordering to the corresponding n-variate exponent.

d0 = getD0(polysys)

Returns the maximal total degree of the given polynomial system.

dorig = getDorig(polysys)

dorig[i] (i=1:s) contains the maximal total degree of the multivariate polynomial corresponding with polysys[i,1],polysys[i,2].

(p,q) = getMDim(polysys,d)

Returns the number of columns p and the number of rows q of the Macaulay matrix of degree d.




