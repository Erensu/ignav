/*---------------------------------------------------------------------------
 * Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
 * (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
 * and B.P. Flannery
 * 
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2017/11/30 1.0 new
-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants and macros -----------------------------------------------------*/
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1)>(dmaxarg2)?(dmaxarg1):(dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1)<(iminarg2)?(iminarg1):(iminarg2))
#define MAX(x,y)  ((x)>=(y)?(x):(y))
#define MIN(x,y)  ((x)<=(y)?(x):(y))

#ifdef LAPACK
extern "C"
{
extern void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a,
                   int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
                   double* work, int* lwork, int* info);
};
#endif

/* compute (a2 + b2)^1/2 without destructive underflow or overflow-----------*/
static double pythag(double a,double b)
{
    double absa,absb;
    absa=fabs(a);
    absb=fabs(b);
    if (absa>absb) {
        return absa*SQRT(1.0+SQR(absb/absa));
    }
    else return (absb==0.0?0.0:absb*SQRT(1.0+SQR(absa/absb)));
}
/*---------------------------------------------------------------------------
 * Given a matrix a[1..m][1..n], this routine computes its singular value
 * decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
 * matrix of singular values W is output as a vector w[1..n].  The matrix V (not
 * the transpose VT) is output as v[1..n][1..n].
 * --------------------------------------------------------------------------*/
static void svdcmp(double **a, int m, int n,double *W, double **V)
{
    int flag,i,its,j,jj,k,l,nm,s2,inc=1;
    double anorm,c,f,g,h,s,scale,x,y,z,*rv1,sw,*su,*sv,*w;

    rv1=(double*)malloc(n*sizeof(double));
    w  =(double*)malloc(n*sizeof(double));
    su =(double*)malloc(m*sizeof(double));
    sv =(double*)malloc(n*sizeof(double));

    g=scale=anorm=0.0; /* householder reduction to bidiagonal form */
    for (i=0;i<n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i<m) {
            for (k=i;k<m;k++) scale+=fabs(a[k][i]);
            if (scale) {
                for (k=i;k<m;k++) {
                    a[k][i]/=scale;
                    s+=a[k][i]*a[k][i];
                }
                f=a[i][i];
                g=-SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=i;k<m;k++) s+=a[k][i]*a[k][j];
                    f=s/h;
                    for (k=i;k<m;k++) a[k][j]+=f*a[k][i];
                }
                for (k=i;k<m;k++) a[k][i]*=scale;
            }
        }
        w[i]=scale*g;
        g=s=scale=0.0;
        if (i<m&&i!=n-1) {
            for (k=l;k<n;k++) scale+=fabs(a[i][k]);
            if (scale) {
                for (k=l;k<n;k++) {
                    a[i][k]/=scale;
                    s+=a[i][k]*a[i][k];
                }
                f=a[i][l];
                g=-SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
                for (j=l;j<m;j++) {
                    for (s=0.0,k=l;k<n;k++) s+=a[j][k]*a[i][k];
                    for (k=l;k<n;k++) a[j][k]+=s*rv1[k];
                }
                for (k=l;k<n;k++) a[i][k]*=scale;
            }
        }
        anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) { /* accumulation of right-hand transformations. */
        if (i<n-1) {
            if (g) {
                for (j=l;j<n;j++) /* double division to avoid possible underflow. */
                    V[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s+=a[i][k]*V[k][j];
                    for (k=l;k<n;k++) V[k][j]+=s*V[k][i];
                }
            }
            for (j=l;j<n;j++) V[i][j]=V[j][i]=0.0;
        }
        V[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=IMIN(m,n)-1;i>=0;i--) { /* accumulation of left-hand transformations. */
        l=i+1;
        g=w[i];
        for (j=l;j<n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            for (j=l;j<n;j++) {
                for (s=0.0,k=l;k<m;k++) s+=a[k][i]*a[k][j];
                f=(s/a[i][i])*g;
                for (k=i;k<m;k++) a[k][j]+=f*a[k][i];
            }
            for (j=i;j<m;j++) a[j][i]*=g;
        }
        else for (j=i;j<m;j++) a[j][i]=0.0;
        ++a[i][i];
    }
    for (k=n-1;k>=0;k--) { /* diagonalization of the bidiagonal form. */
        for (its=0;its<30;its++) {
            flag=1;
            for (l=k;l>=0;l--) { /* Test for splitting. */
                nm=l-1; /* note that rv1[1] is always zero. */
                if ((fabs(rv1[l])+anorm)==anorm) {
                    flag=0;
                    break;
                }
                if ((fabs(w[nm])+anorm)==anorm) break;
            }
            if (flag) {
                c=0.0; /* cancellation of rv1[l], if l > 1. */
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    rv1[i]=c*rv1[i];
                    if ((fabs(f)+anorm)==anorm) break;
                    g=w[i];
                    h=pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s=-f*h;
                    for (j=0;j<m;j++) {
                        y=a[j][nm];
                        z=a[j][i];
                        a[j][nm]=y*c+z*s;
                        a[j][i]=z*c-y*s;
                    }
                }
            }
            z=w[k];
            if (l==k) { /* convergence. */
                if (z<0.0) { /* singular value is made nonnegative. */
                    w[k]=-z;
                    for (j=0;j<n;j++) V[j][k]=-V[j][k];
                }
                break;
            }
            if (its==29) fprintf(stderr,"no convergence in 30 svdcmp iterations\n");
            x=w[l]; /* shift from bottom 2-by-2 minor. */
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0; /* next QR transformation: */
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=pythag(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y*=c;
                for (jj=0;jj<n;jj++) {
                    x=V[jj][j];
                    z=V[jj][i];
                    V[jj][j]=x*c+z*s;
                    V[jj][i]=z*c-x*s;
                }
                z=pythag(f,h);
                w[j]=z; /* rotation can be arbitrary if z=0. */
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;
                for (jj=0;jj<m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
    do {inc*=3; inc++;}while (inc<=n);
    do {
        inc/=3;
        for (i=inc;i<n;i++) {
            sw=w[i];
            for (k=0;k<m;k++) su[k]=a[k][i];
            for (k=0;k<n;k++) sv[k]=V[k][i];
            j=i;
            while (w[j-inc]<sw) {
                w[j]=w[j-inc];
                for (k=0;k<m;k++) a[k][j]=a[k][j-inc];
                for (k=0;k<n;k++) V[k][j]=V[k][j-inc];
                j-=inc;
                if (j<inc) break;
            }
            w[j]=sw;
            for (k=0;k<m;k++) a[k][j]=su[k];
            for (k=0;k<n;k++) V[k][j]=sv[k];
        }
    } while (inc>1);
    for (k=0;k<n;k++) {
        s2=0;
        for (i=0;i<m;i++) if (a[i][k]<0.0) s2++;
        for (j=0;j<n;j++) if (V[j][k]<0.0) s2++;
        if (s2>(m+n)/2) {
            for (i=0;i<m;i++) a[i][k]=-a[i][k];
            for (j=0;j<n;j++) V[j][k]=-V[j][k];
        }
    }
    for (i=0;i<MIN(m,n);i++) W[i]=w[i];

    free(rv1);
    free(w);
    free(su);
    free(sv);
}
/* allocate memory for matrix------------------------------------------------*/
extern double** dmat(int m,int n)
{
    int i; double **A;
    if (m<=0||n<=0) return NULL;
    A=(double**)malloc(m*sizeof(double*));
    A[0]=(double*)calloc(m*n,sizeof(double));

    for (i=1;i<m;i++) A[i]=A[i-1]+n;
    return A;
}
/* svd decomposition---------------------------------------------------------*/
extern int svd(const double *A,int m,int n,double *U,double *W,double *V)
{
#ifdef LAPACK
    register int info=-1,lwork=-1,lda=m,ldu=m,ldvt=n,i,j;
    double *s,*u,*vt,*a;
    double wkopt;
    double* work;

    s=mat(1,MAX(m,n)); u=mat(1,m*m); vt=mat(1,n*n);
    a=mat(m,n);

    matcpy(a,A,m,n);

    dgesvd_("All","All",&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,&wkopt,&lwork,&info);
    lwork=(int)wkopt;
    work =(double*)malloc(lwork*sizeof(double));
    dgesvd_("All","All",&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);
    if (info>0) {
        trace(3,"the algorithm computing SVD failed to converge.\n");

        free((void*)work);
        free(a); free(u); free(vt); free(s);
        return 0;
    }
    for (i=0;i<m&&U;i++) {
        for (j=0;j<MIN(m,n);j++) U[i+j*m]=u[i+j*ldu];
    }
    for (i=0;i<n&&V;i++) for (j=0;j<n;j++) V[i+j*n]=vt[j+i*ldvt];
    for (i=0;i<MIN(m,n)&&W;i++) W[i]=s[i];

    free(work);
    free(a); free(u); free(vt); free(s);
    return 1;
#else
    register int i,j;
    double **a,**v,**u,*w;

    w=mat(1,MAX(m,n));
    a=dmat(m,n); v=dmat(n,n); u=dmat(m,m);

    mat2dmat(A,a,m,n); svdcmp(a,m,n,w,v);
    if (U) {
        for (i=0;i<m;i++) {
            for (j=0;j<MIN(m,n);j++) U[i+j*m]=a[i][j];
        }
    }
    if (W) {
        for (i=0;i<MIN(m,n);i++) W[i]=w[i];
    }
    if (V) {
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) V[i+j*n]=v[i][j];
        }
    }
    free(w);
    free(a[0]); free(a);
    free(v[0]); free(v);
    free(u[0]); free(u); return 1;
#endif
}
