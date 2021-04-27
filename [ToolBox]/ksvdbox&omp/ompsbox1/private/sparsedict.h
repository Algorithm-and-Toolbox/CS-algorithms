/**************************************************************************
 *
 * File name: sparsedict.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 18.8.2009
 *
 * Implementation of sparse dictionary operators.
 *
 *************************************************************************/


#ifndef __SPARSEDICT_H__
#define __SPARSEDICT_H__

#include "mex.h"



/**************************************************************************
 * Dictionary operator for sparse-separable dictionary.
 *
 * Compute D*x where D=B*A is a sparse dictionary with separable base
 * dictionary, and x is a sparse vector.
 * 
 * Parameters:
 *   B[]   - array of either 2 or 3 matrices. B[i] is the base dictionary for
 *           the i-th dimension, of size n[i] X k[i]
 *   Apr,Air,Ajc - sparse representation of the matrix A, of size prod{k[i]} X m
 *   xpr,xir,xjc - sparse representation of the vector x, of length m
 *   Dx    - output signal of length prod{n[i]}
 *   ndims - length of the array B (2 or 3), also the # of dimensions of D*x
 *   n, k  - sizes of the matrices B[i]
 *   m     - number of columns in A
 *
 * Note: This function re-writes the contents of Dx.
 *
 **************************************************************************/
void dict_vec_sep(double *B[], double Apr[], mwIndex Air[], mwIndex Ajc[], 
        double xpr[], mwIndex xir[], mwIndex xjc[], double Dx[], int ndims, 
        mwSize n[], mwSize k[], mwSize m);



/**************************************************************************
 * Dictionary-transpose operator for sparse-separable dictionary.
 *
 * Compute D'*x where D=B*A is a sparse dictionary with separable base
 * dictionary, and x is a 2 or 3 dimensional signal.
 * 
 * Parameters:
 *   B[]   - array of either 2 or 3 matrices. B[i] is the base dictionary for
 *           the i-th dimension, of size n[i] X k[i]
 *   Apr,Air,Ajc - sparse representation of the matrix A, of size prod{k[i]} X m
 *   x     - signal of length prod{n[i]}
 *   DtX   - output signal of length m
 *   ndims - length of the array B (2 or 3), also the # of dimensions of x
 *   n, k  - sizes of the matrices B[i]
 *   m     - number of columns in A
 *
 * Note: This function re-writes the contents of DtX.
 *
 **************************************************************************/
void dictT_vec_sep(double *B[], double Apr[], mwIndex Air[], mwIndex Ajc[],
        double x[], double DtX[], int ndims, mwSize n[], mwSize k[], mwSize m);


#endif

