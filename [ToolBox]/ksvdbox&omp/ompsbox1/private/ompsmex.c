/**************************************************************************
 *
 * File name: ompsmex.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last Updated: 25.8.2009
 *
 *************************************************************************/

#include "ompscore.h"
#include "omputils.h"
#include "mexutils.h"


/* Input Arguments */

#define IN_B          prhs[0]
#define IN_A          prhs[1]
#define IN_X          prhs[2]
#define IN_G          prhs[3]
#define IN_T          prhs[4]
#define IN_SPARSE_G   prhs[5]
#define IN_MSGDELTA   prhs[6]
#define IN_PROFILE    prhs[7]


/* Output Arguments */

#define	GAMMA_OUT     plhs[0]


/***************************************************************************************/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])

{
  char Bi_str[12];
  mxArray *Bi;
  double *B[3], msgdelta;
  int gmode, profile, g_specified, T;
  mwSize n[3], k[3], m, L, N, K, ndims, i;   /* B{i} is n[i] x k[i], B is N x K, A is K x m, X is N x L */
  
  
  /* check parameters */
  
  checkcell_1d(IN_B, "OMPS", "Bsep");
  checkmatrix(IN_A, "OMPS", "A");
  checksparse(IN_A, "OMPS", "A");
  checkmatrix(IN_X, "OMPS", "X");
  checkmatrix(IN_G, "OMPS", "G");
  
  checkscalar(IN_T, "OMPS", "T");
  checkscalar(IN_SPARSE_G, "OMPS", "sparse_g");
  checkscalar(IN_MSGDELTA, "OMPS", "msgdelta");
  checkscalar(IN_PROFILE, "OMPS", "profile");
  
  g_specified = !mxIsEmpty(IN_G);
  
  ndims = mxGetM(IN_B)*mxGetN(IN_B);
  if (ndims<2 || ndims>3) {
    mexErrMsgTxt("Signal dimension must be 2 or 3.");
  }
  
  for (i=0; i<ndims; ++i) {
    Bi = mxGetCell(IN_B,i);
    sprintf(Bi_str, "Bsep{%d}", i+1);
    checkmatrix(Bi, "OMPS", Bi_str);
    
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
  
  if (N != mxGetM(IN_X)) {
    mexErrMsgTxt("B and X have incompatible sizes.");
  }
  
  if (g_specified) {
    if (mxGetN(IN_G)!= mxGetM(IN_G)) {
      mexErrMsgTxt("G must be a square matrix.");
    }
    if (mxGetN(IN_G) != m) {
      mexErrMsgTxt("A and G have incompatible sizes.");
    }
  }
  
  
  
  /* get parameters */
  
  T = (int)(mxGetScalar(IN_T)+1e-2);
  if ((int)(mxGetScalar(IN_SPARSE_G)+1e-2)) {
    gmode = SPARSE_GAMMA;
  }
  else {
    gmode = FULL_GAMMA;
  }
  msgdelta = mxGetScalar(IN_MSGDELTA);
  profile = (int)(mxGetScalar(IN_PROFILE)+1e-2);
  
  
  
  /* Do OMP! */
  
  GAMMA_OUT = ompscore_sep(B, mxGetPr(IN_A), mxGetIr(IN_A), mxGetJc(IN_A), mxGetPr(IN_X), g_specified ? mxGetPr(IN_G) : 0,
                ndims, n, k, m, L, T, 0, gmode, profile, msgdelta, 0);
  
  return;
}

