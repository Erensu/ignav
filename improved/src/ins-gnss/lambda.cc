/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*     [3] Paul de Jonge, Christian Tiberius. The LAMBDA method for integer
*         ambiguity estimation: implementation aspects
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants/macros ----------------------------------------------------------*/
#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {
#if 0
            info=-1; break;
#else
            D[i]=fabs(A[i+i*n]);
#endif
        }
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info<0) {
        trace(2,"\"%s : LD factorization error, trying UD from Gibbs"
                "(col major UD=LD)\n",__FILE__);
        double Qc[n*n];
        memcpy(Qc,Q,n*n*sizeof(double));
        matrix_udu(n,Qc,L,D);
        info=0; 
    }
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/(fabs(D[k])==0.0?SGN(D[k])*1E-20:D[k]);
        if (newdist<maxdist) {
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        trace(2,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s)
{
    int info;
    double *L,*D,*Z,*z,*E;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction */
        reduction(n,L,D,Z);
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
        
        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}
/*---------------------------------------------------------------------------*/
static int exctract_L(const double *L,int k,double *LL,int m,int n)
{
    int i;
    for (i=0;i<n;i++) LL[i]=L[m-1+k*i]; return n;
}
/* compute the bootstrapped success rate ------------------------------------*/
static double amb_bs_success(const double *D,int n)
{
    double s=1.0;
    int i; for (i=0;i<n;i++) s*=(2.0*norm_distri(0.5/SQRT(D[i]))-1.0);
    return s;
}
/* bootstrap search ambiguity ------------------------------------------------
 * args   : int    n      I  number of float parameters
 *          double *a     I  float parameters (n x 1)
 *          double *Q     I  covariance matrix of float parameters (n x n)
 *          double *F     O  fixed solutions (n x m)
 *          double *s     O  sum of squared residulas of fixed solutions (1 x m)
 *          double *Ps    O  success ratio
 * return : status (0:ok,other:error)
 * ---------------------------------------------------------------------------*/
extern int bootstrap(int n,const double *a, const double *Q, double *F,double *Ps)
{
    double *L,*D,*Z,*an,*ZT;
    double *af,*afc,*S,*zhat,*LL;
    int i,j,k,info;

    if (n<=0) return -1;

    L=zeros(n,n); D=mat(n,1); Z=eye(n);

    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D); free(Z); return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);

    if (Ps) *Ps=amb_bs_success(D,n);

    zhat=mat(n,1);
    an=mat(n,1);
    ZT=mat(n,n);

    for (i=0;i<n;i++) an[i]=a[i]-int(a[i]);

    matmul("TN",n,1,n,1.0,Z,an,0.0,zhat); /* z=Z'*a */

    af=zeros(n,1); afc=zeros(n,1);
    S=zeros(n,1); LL=mat(n,1);

    afc[n-1]=zhat[n-1]; af[n-1]=ROUND(afc[n-1]);

    for (i=n-1;i>0;i--) {
        j=exctract_L(L,n,LL,i+1,i);

        for (k=0;k<j;k++) S[k]=S[k]+(af[j]-afc[j])*LL[k];
        afc[j-1]=zhat[j-1]+S[j-1];
        af[j-1]=ROUND(afc[j-1]);
    }

    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) ZT[i+j*n]=Z[j+i*n];
    }
    if (!(info=matinv(ZT,n))) {
        matmul("NN",n,1,n,1.0,ZT,af,0.0,F);
        for (i=0;i<n;i++) F[i]+=int(a[i]);
    }
    free(L); free(D); free(Z);
    free(af); free(afc);
    free(S); free(LL);
    return info;
}
/*---------------------------------------------------------------------------*/
static int exctract_vec(const double *D,int n,int k,double *Do)
{
    int i,j; for (i=k-1,j=0;i<n;i++) Do[j++]=D[i]; return j;
}
/*---------------------------------------------------------------------------*/
static void exctract_mat(const double *A,int r,int c,int a1,int b1,
                         int a2,int b2,double *Ao)
{
    int i,j,rs,re,cs,ce,nr,nc;

    rs=MIN(a1,a2); re=MAX(a1,a2);
    cs=MIN(b1,b2); ce=MAX(b1,b2);
    nr=re-rs+1; nc=ce-cs+1;

    for (i=0;i<nc;i++) {
        for (j=0;j<nr;j++) Ao[i*nr+j]=A[(cs-1+i)*r+rs-1+j];
    }
}
/*---------------------------------------------------------------------------*/
static int parsearch(const double *zhat,const double *Qzhat,const double *Z,
                     const double *L,const double *D,double P0,int n,int m,
                     double *F,double*s)
{
    int i,j,k=1,info=0,l;
    double ps,*DD,*zz,*LL,*Qzpar,*QP,*zpar;

    ps=amb_bs_success(D,n);
    DD=mat(1,n); zz=mat(1,n); LL=mat(n,n);

    while(ps<P0&&k<n) {
        k=k+1;
        j=exctract_vec(D,n,k,DD);
        ps=amb_bs_success(DD,j);
    }
    Qzpar=mat(n,n); QP=mat(n,n);
    zpar =mat(m,n);

    if (ps>P0) {
        if (k==1) {
            search(n,m,L,D,zhat,F,s); 
        }
        else {
            exctract_vec(zhat,n,k,zz); exctract_vec(D,n,k,DD);

            exctract_mat(L,n,n,k,k,n,n,LL);

            search(n-k+1,m,LL,DD,zz,zpar,s);

            exctract_mat(Qzhat,n,n,k,k,n,n,Qzpar);
            exctract_mat(Qzhat,n,n,k-1,k,1,n,LL);

            if (!(info=matinv(Qzpar,n-k+1))) {
                matmul("NN",k-1,n-k+1,n-k+1,1.0,LL,Qzpar,0.0,QP);

                for (i=0;i<m;i++) {
                    for (j=0;j<n-k+1;j++) {
                        DD[j]=zz[j]-zpar[i*(n-k+1)+j];
                    }
                    matmul("NN",k-1,1,n-k+1,1.0,QP,DD,0.0,LL);
                    for (j=0;j<k-1;j++) {
                        F[i*n+j]=ROUND(zhat[j]-LL[j]);
                    }
                    for (j=k-1,l=0;j<n;j++) {
                        F[i*n+j]=zpar[i*(n-k+1)+l++];
                    }
                }
            }
        }
    }
    free(DD); free(LL); free(zz);
    free(Qzpar); free(QP);
    free(zpar);
    return info;
}
/* partial ambiguity search---------------------------------------------------
 * args   : int    n      I  number of float parameters
 *          double *a     I  float parameters (n x 1)
 *          double *Q     I  covariance matrix of float parameters (n x n)
 *          double *F     O  fixed solutions (n x m)
 *          double *s     O  sum of squared residuals of fixed solutions (1 x m)
 *          double *p0    O  success ratio
 * return : status (0:ok,other:error)
 * --------------------------------------------------------------------------*/
extern int plambda(const double *a,const double *Qa,int n,int m,double *F,
                   double *s,double p0)
{
    double *L,*D,*Z,*zhat,*Qz,*Q,*E;
    int info;

    trace(3,"par_lambda:\n");

    if (n<=0) return -1;

    L=zeros(n,n); D=mat(n,1); Z=eye(n); E=mat(m,n);

    /* LD factorization */
    if ((info=LD(n,Qa,L,D))) {
        free(L); free(D); free(Z); return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);

    zhat=mat(n,1);
    Qz=mat(n,n); Q=mat(n,n);

    matmul("TN",n,1,n,1.0,Z,a,0.0,zhat); /* z=Z'*a */
    matmul("TN",n,n,n,1.0,Z,Qa,0.0,Q);
    matmul("NN",n,n,n,1.0,Q,Z,0.0,Qz); /* Qz=Z'*Qa*Z */

    parsearch(zhat,Qz,Z,L,D,p0,n,m,E,s);

    info=solve("T",Z,E,n,m,F); /* F=Z'\E */

    free(L); free(D); free(Z);
    free(zhat); free(Qz); free(Q);
    free(E);
    return info;
}
