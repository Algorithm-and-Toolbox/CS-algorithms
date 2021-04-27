/**************************************************************************
 *
 * File name: ompscore.h
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 25.8.2009
 *
 * Contains the core implementation of sparse Batch-OMP / OMP-Cholesky.
 *
 *************************************************************************/


#ifndef __OMPS_CORE_H__
#define __OMPS_CORE_H__


#include "mex.h"



/**************************************************************************
 * Perform sparse-separable Batch-OMP or OMP-Cholesky on the specified set
 * of signals, using either a fixed number of atoms or an error bound. 
 * The function supports 2D and 3D signals.
 *
 * Parameters:
 *
 *   B - the separable base dictionary, represented as an array of of 1-D dictionaries.
 *     B[i] is the dictionary for the i-th dimension, and should be of size n[i] X k[i].
 *     The array B itself should be of length 2 or 3
 *   Apr, Air, Ajc - the sparse representation of the matrix A, of size prod{k[i]} X m
 *   x - the column signals, arranged in column-major order, of size prod{n[i]} X L
 *   G - D'*D, of size m X m
 *   T - target sparsity, or maximal number of atoms for error-based OMP
 *   eps - target residual norm for error-based OMP
 *   gamma_mode - one of the constants FULL_GAMMA or SPARSE_GAMMA
 *   profile - if non-zero, profiling info is printed
 *   msg_delta - positive: the # of seconds between status prints, otherwise: nothing is printed
 *   erroromp - if nonzero indicates error-based OMP, otherwise fixed sparsity OMP
 *
 * Usage:
 *
 *   This function can be called either with or without the parameter G. 
 *   When G is specified, Batch-OMP is performed. When G is null,
 *   OMP-Cholesky is performed.
 *
 *   Fixed-sparsity usage:
 *   ---------------------
 *   The number of atoms must be specified in T. The value of eps is ignored.
 *   Set erroromp to 0.
 *
 *   Error-OMP usage:
 *   ----------------
 *   The target error must be specified in eps. A hard limit on the number
 *   of atoms can also be specified via the parameter T. Otherwise, T should 
 *   be negative. Finally, set erroromp to nonzero.
 *
 *
 * Returns: 
 *   An mxArray containing the sparse representations of the signals in x
 *   (allocated using the appropriate mxCreateXXX() function).
 *   The array is either full or sparse, depending on gamma_mode.
 *
 **************************************************************************/
mxArray* ompscore_sep(double *B[], double Apr[], mwIndex Air[], mwIndex Ajc[], double x[], double G[], int ndims,
             mwSize n[], mwSize k[], mwSize m, mwSize L, int T, double eps, int gamma_mode, int profile, double msg_delta, int erroromp);


#endif

