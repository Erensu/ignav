/*---------------------------------------------------------------------------
 * null.cc : compute null-space of matrix
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/11 1.0 new
-----------------------------------------------------------------------------*/
#include "carvig.h"

#define EPS  1E-10

/* null-space of matrix------------------------------------------------------
 * args:    double *A  I  input matrix (mxn)
 *          int m,n    I  size of matrix
 *          double *N  O  null-space matrix of `A'
 *          int *p,*q  I  size of N (pxq)
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int null(const double *A,int m,int n,double *N,int *p,int *q)
{
    int i,j,k,*index;
    double *U,*V,*W;

    trace(3,"null: m=%d n=%d\n",m,n);

    if (m>n) return 0;
    U=mat(m,m); V=mat(n,n); W=mat(m,n); index=imat(1,MAX(m,n));
    svd(A,m,n,U,W,V);

    for (i=0,k=0;i<m;i++) {
        if (fabs(W[i])<EPS) index[k++]=i;
    }
    for (i=0;i<n-m;i++) index[k++]=m+i;
    for (i=0;i<k;i++) {
        for (j=0;j<n;j++) N[i*n+j]=V[index[i]*n+j];
    }
    *p=n; *q=k;
    free(U); free(V); free(W); free(index);
    return 1;
}
