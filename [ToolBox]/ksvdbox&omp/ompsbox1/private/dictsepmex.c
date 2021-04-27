/**************************************************************************
 *
 * File name: dictsepmex.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 31.8.2009
 *
 *************************************************************************/


#include "sparsedict.h"
#include "omputils.h"


/* Input Arguments */

#define IN_B          prhs[0]
#define IN_A          prhs[1]
#define IN_X          prhs[2]


/* Output Arguments */

#define	Y_OUT         plhs[0]


/***************************************************************************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

{
  char Bi_str[8];
  mxArray *Bi;
  double *B[3], *y, *xpr;
  mwIndex *xir, *xjc, xjc_current[2];
  mwSize n[3], k[3], m, L, N, K, ndims, i;   /* B{i} is n[i] x k[i], B is N x K, A is K x m, X is N x L */
  
  
  /* check parameters */
  
  checkcell_1d(IN_B, "DICTSEP", "B");
  checkmatrix(IN_A, "DICTSEP", "A");
  checksparse(IN_A, "DICTSEP", "A");
  checkmatrix(IN_X, "DICTSEP", "X");
  checksparse(IN_X, "DICTSEP", "X");
  
  ndims = mxGetM(IN_B)*mxGetN(IN_B);
  if (ndims<2 || ndims>3) {
    mexErrMsgTxt("Signal dimension must be 2 or 3.");
  }
  
  for (i=0; i<ndims; ++i) {
    Bi = mxGetCell(IN_B,i);
    sprintf(Bi_str, "B{%d}", i+1);
    checkmatrix(Bi, "DICTSEP", Bi_str);
    
    B[i] = mxGetPr(Bi);
    n[i] = mxGetM(Bi);
    k[i] = mxGetN(Bi);
  }
  
  
  
  /* check sizes */
  
  if (ndims==2) {
    N = n[0]*n[1];
    K = k[0]*k[1];
  }
  else {
    N = n[0]*n[1]*n[2];
    K = k[0]*k[1]*k[2];
  }
  
  m = mxGetN(IN_A);
  L = mxGetN(IN_X);
  
  if (K != mxGetM(IN_A)) {
    mexErrMsgTxt("Base dictionary and matrix A have incompatible sizes.");
  }
  
  if (m != mxGetM(IN_X)) {
    mexErrMsgTxt("A and X have incompatible sizes.");
  }
  
  
  /* allocate output matrix */
  
  Y_OUT = mxCreateDoubleMatrix(N, L, mxREAL);
  y = mxGetPr(Y_OUT);
  
  
  /* Compute */
  
  xpr = mxGetPr(IN_X);
  xir = mxGetIr(IN_X);
  xjc = mxGetJc(IN_X);
  
  xjc_current[0] = 0;
  
  for (i=0; i<L; ++i) {
    xjc_current[1] = xjc[i+1]-xjc[i];
    if (xjc_current[1]==0) {
      continue;
    }
    dict_vec_sep(B, mxGetPr(IN_A), mxGetIr(IN_A), mxGetJc(IN_A),
        xpr+xjc[i], xir+xjc[i], xjc_current, y+N*i, ndims, n, k, m);
  }
  
  return;
}

