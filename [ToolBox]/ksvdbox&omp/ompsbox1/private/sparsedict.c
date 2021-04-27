/**************************************************************************
 *
 * File name: sparsedict.c
 *
 * Ron Rubinstein
 * Computer Science Department
 * Technion, Haifa 32000 Israel
 * ronrubin@cs
 *
 * Last updated: 14.8.2009
 *
 *************************************************************************/


#include "sparsedict.h"
#include "myblas.h"



/* compute D*x for a sparse dictionary and sparse x */

void dict_vec_sep(double *B[], double Apr[], mwIndex Air[], mwIndex Ajc[], double xpr[], mwIndex xir[], mwIndex xjc[], double Dx[], int ndims, mwSize n[], mwSize k[], mwSize m)
{
  mwIndex i;
  double *Ax, *tmp, *tmp2, *Bt;
  
  tmp = mxMalloc(n[0]*k[1]*sizeof(double));
  Bt = mxMalloc(n[1]*k[1]*sizeof(double));
  transpose(B[1], Bt, n[1], k[1]);
  
  if (ndims==2) {
    Ax = mxMalloc(k[0]*k[1]*sizeof(double));
    mat_sp_vec_sp(1, Apr, Air, Ajc, xpr, xir, xjc, Ax, k[0]*k[1], m);
    mat_mat(1, B[0], Ax, tmp, n[0], k[0], k[1]);
    mat_mat(1, tmp, Bt, Dx, n[0], k[1], n[1]);
  }
  else {   /* ndims==3 */
    
    Ax = mxMalloc(k[0]*k[1]*k[2]*sizeof(double));
    tmp2 = mxMalloc(n[0]*n[1]*k[2]*sizeof(double));
    
    mat_sp_vec_sp(1, Apr, Air, Ajc, xpr, xir, xjc, Ax, k[0]*k[1]*k[2], m);
    
    for (i=0; i<k[2]; ++i) {
      mat_mat(1, B[0], Ax+i*k[0]*k[1], tmp, n[0], k[0], k[1]);
      mat_mat(1, tmp, Bt, tmp2+i*n[0]*n[1], n[0], k[1], n[1]);
    }
    
    tens_mat(1, tmp2, B[2], Dx, n[0], n[1], k[2], n[2]);
    
    mxFree(tmp2);
    
  }
  
  mxFree(Ax);
  mxFree(Bt);
  mxFree(tmp);
}



/* compute D'*x for a sparse dictionary */

void dictT_vec_sep(double *B[], double Apr[], mwIndex Air[], mwIndex Ajc[], double x[], double DtX[], int ndims, mwSize n[], mwSize k[], mwSize m)
{
  mwIndex i; 
  double *tmp, *tmp2, *BtX;
 
  tmp = mxMalloc(k[0]*n[1]*sizeof(double));
  
  if (ndims==2) {   
    BtX = mxMalloc(k[0]*k[1]*sizeof(double));
    matT_mat(1, B[0], x, tmp, n[0], k[0], n[1]);
    mat_mat(1, tmp, B[1], BtX, k[0], n[1], k[1]);
    matT_sp_vec(1, Apr, Air, Ajc, BtX, DtX, k[0]*k[1], m);
  }
  else {   /* ndims==3 */

    BtX = mxMalloc(k[0]*k[1]*k[2]*sizeof(double));
    tmp2 = mxMalloc(k[0]*k[1]*n[2]*sizeof(double));
    
    for (i=0; i<n[2]; ++i) {
      matT_mat(1, B[0], x+i*n[0]*n[1], tmp, n[0], k[0], n[1]);
      mat_mat(1, tmp, B[1], tmp2+i*k[0]*k[1], k[0], n[1], k[1]);
    }
    
    tens_matT(1, tmp2, B[2], BtX, k[0], k[1], n[2], k[2]);
    
    matT_sp_vec(1, Apr, Air, Ajc, BtX, DtX, k[0]*k[1]*k[2], m);
    
    mxFree(tmp2);
    
  }
  
  mxFree(BtX);
  mxFree(tmp);
}



