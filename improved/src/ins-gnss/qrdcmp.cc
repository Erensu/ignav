/*---------------------------------------------------------------------------
 * QR decomposition program, svdcmp, from "Numerical Recipes in C"
 * (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
 * and B.P. Flannery
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/01/04 1.0 new
-----------------------------------------------------------------------------*/
#include "carvig.h"

/* helper function of QR decomposition----------------------------------------*/
static int qrdcmp(double **A, int m, int n,double **Qo,double **Ro)
{
    int i,j,k,nn,jj;
    double u,alpha,w,t,**Q=NULL;
    if (m<n) return 0;

    Q=dmat(m,m);
    for (i=0;i<=m-1;i++) for (j=0;j<=m-1;j++) {
            Q[i][j]=0.0;
            if (i==j) Q[i][j]=1.0;
        }
    nn=n;
    if (m==n) nn=m-1;
    for (k=0;k<=nn-1;k++) {
        u=0.0;
        for (i=k;i<=m-1;i++) {
            w=fabs(A[i][k]); if (w>u) u=w;
        }
        alpha=0.0;
        if (fabs(u)<1E-20) continue;
        for (i=k;i<=m-1;i++) {
            t=A[i][k]/u;
            alpha=alpha+t*t;
        }
        if (A[k][k]>0.0) u=-u; alpha=u*sqrt(alpha);
        if (fabs(alpha)+1.0==1.0) {
            free(Q[0]);
            free(Q);
            return 0;
        }
        u=sqrt(2.0*alpha*(alpha-A[k][k]));

        if ((u+1.0)!=1.0) {
            A[k][k]=(A[k][k]-alpha)/u;

            for (i=k+1;i<=m-1;i++) A[i][k]=A[i][k]/u;
            for (j=0;j<=m-1;j++) {
                t=0.0;
                for (jj=k;jj<=m-1;jj++) t=t+A[jj][k]*Q[jj][j];
                for (i=k;i<=m-1;i++) {
                    Q[i][j]=Q[i][j]-2.0*t*A[i][k];
                }
            }
            for (j=k+1;j<=n-1;j++) {
                t=0.0;
                for (jj=k;jj<=m-1;jj++) t=t+A[jj][k]*A[jj][j];
                for (i=k;i<=m-1;i++) A[i][j]=A[i][j]-2.0*t*A[i][k];
            }
            A[k][k]=alpha;
            for (i=k+1;i<=m-1;i++) {
                A[i][k]=0.0;
            }
        }
    }
    for (i=0;i<=m-2;i++)
        for (j=i+1;j<=m-1;j++) {
            t=Q[i][j];
            Q[i][j]=Q[j][i];
            Q[j][i]=t;
        }
    for (i=0;i<m;i++) {
        for (j=0;j<m;j++) Qo[i][j]=Q[i][j];
    }
    for (i=0;i<m;i++) {
        for (j=0;j<n;j++) {
            Ro[i][j]=A[i][j];
        }
    }
    free(Q[0]);
    free(Q);
    return 1;
}
/* QR decomposition function--------------------------------------------------
 * args:    double *A     I  input need to decomposition matrix
 *          int m,n       I  size of A (mxn)
 *          double *Q,*R  O  output matrix: Q (mxm), R (mxn)
 *          int flag      I  flag=0: col-major,
 *                           flag=1: row-major
 *
 * return :status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int qr(const double *A,int m,int n,double *Q,double *R,int flag)
{
    double **Ao,**Qo,**Ro;
    int i,j;

    Ao=dmat(m,n);
    Qo=dmat(m,m);
    Ro=dmat(m,n);

    for (i=0;i<m;i++) {
        for (j=0;j<n;j++) Ao[i][j]=(flag==0?A[i+j*m]:A[j+i*n]);
    }
    if (!qrdcmp(Ao,m,n,Qo,Ro)) {
        free(Qo[0]); free(Ro[0]);
        free(Ao[0]);
        free(Qo);
        free(Ao);
        free(Ro); return 0;
    }
    for (i=0;i<m;i++) {
        for (j=0;j<m&&Q;j++) Q[i+j*m]=Qo[i][j];
        for (j=0;j<n&&R;j++) R[i+j*m]=Ro[i][j];
    }
    free(Qo[0]); free(Ro[0]);
    free(Ao[0]);
    free(Qo);
    free(Ao);
    free(Ro); return 1;
}

