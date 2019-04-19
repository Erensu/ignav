/*---------------------------------------------------------------------------
* ins-vo-mono.cc : monocular camera estimate motion functions
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/05/10 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

#define MAXVOGDOP     1.5       /* reject threshold of gdop */
#define REFINE_F      1         /* refine fundamental matrix */

/* struct type---------------------------------------------------------------*/
typedef struct {                /* distance data type for 3-D points */
    double d;                   /* distance */
    int id;                     /* id of point */
} dist_t;

typedef struct {                /* structure for storing matches */
    double u1p,v1p;             /* u,v-coordinates in previous image */
    double u1c,v1c;             /* u,v-coordinates in current image */
} fea_t;

typedef struct {                /* structure for storing a frame image */
    gtime_t gtime;              /* time stamp */
    int n,nmax;                 /* number of feature data/allocated */
    fea_t *data;                /* feature data records */
} frame_t;

/* convert (**)matrix to (*)matrix-------------------------------------------
 * args   :  double **a  I  (**) matrix
 *           double *b   O  (*) matrix
 *           int m       I  rows of matrix
 *           int n       I  cols of matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void dmat2mat(double **a,double *b,int m,int n)
{
    register int i,j;
    for (i=0;i<m;i++) for (j=0;j<n;j++) b[i+j*m]=a[i][j];
}
/* convert (*) matrix to (**) matrix-----------------------------------------
 * args   :  double *a   I  (*) matrix
 *           double **b  O  (**) matrix
 *           int m,n     I  rows and cols of matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void mat2dmat(const double *a,double **b,int m,int n)
{
    register int i,j;
    for (i=0;i<m;i++) for (j=0;j<n;j++) b[i][j]=a[i+j*m];
}
/* construct dialog matrix from given vector--------------------------------*/
extern void dialog(const double *v,int n,double *D)
{
    register int i,j;
    for (i=0;i<n;i++) for (j=0;j<n;j++) D[i+j*n]=i==j?v[i]:0.0;
}
/* compute fundamental matrix------------------------------------------------
 * args   :  frame_t *frame  I  input feature points
 *           int *active     I  active index of feature points
 *           int n           I  number of active feature points
 *           double *F       O  fundamental matrix
 * return : none
 * --------------------------------------------------------------------------*/
static void fundamental(const frame_t *frame,const int *active,int n,double *F)
{
    register int i;
    double *A,*U,*W,*V,F_[9];
    fea_t p;

    A=mat(n,9); U=mat(n,n);
    W=mat(1,9); V=mat(9,9);

    /* create constraint matrix A */
    for (i=0;i<n;i++) {
        p=frame->data[active[i]];

        A[i+0*n]=p.u1c*p.u1p; A[i+1*n]=p.u1c*p.v1p;
        A[i+2*n]=p.u1c;       A[i+3*n]=p.v1c*p.u1p;
        A[i+4*n]=p.v1c*p.v1p; A[i+5*n]=p.v1c;
        A[i+6*n]=p.u1p;       A[i+7*n]=p.v1p;
        A[i+8*n]=1.0;
    }
    /* compute singular value decomposition of A */
    svd(A,n,9,U,W,V);

    /* extract fundamental matrix from the column of
     * V corresponding to the smallest singular value
     * */
    F_[0]=V[0+8*9]; F_[3]=V[1+8*9]; F_[6]=V[2+8*9];
    F_[1]=V[3+8*9]; F_[4]=V[4+8*9]; F_[7]=V[5+8*9];
    F_[2]=V[6+8*9]; F_[5]=V[7+8*9]; F_[8]=V[8+8*9];

    svd(F_,3,3,U,W,V);

    /* enforce rank 2 */
    W[2]=0.0;

    dialog(W,3,F_);
    matmul33("NNT",U,F_,V,3,3,3,3,F);

    free(A); free(U);
    free(W); free(V);
}
/* get camera calibration matrix-------------------------------------------*/
static void getK(const voopt_t *opt,double *K)
{
    K[0]=opt->calib.fu; K[1]=0.0;           K[2]=0.0;
    K[3]=0.0;           K[4]=opt->calib.fv; K[5]=0.0;
    K[6]=opt->calib.cu; K[7]=opt->calib.cv; K[8]=1.0;
}
/* compute essential matrix -------------------------------------------------*/
static void essential(const double *Tc,const double *Tp,const double *F,const double *K,
                      double *E)
{
    double *E1,*E2,*U,*V,*W,*D;

    E1=mat(3,3); E2=mat(3,3);
    U=mat(3,3); V=mat(3,3); W=mat(1,3);
    D=zeros(3,3);

    /* de-normalise and extract essential matrix */
    matmul33("TNN",Tc,F,Tp,3,3,3,3,E1);
    matmul33("TNN",K,E1,K,3,3,3,3,E2);

    /* re-enforce rank 2 constraint on essential matrix */
    svd(E2,3,3,U,W,V);

    W[2]=0.0; dialog(W,3,D);

    /* essential matrix */
    matmul33("NNT",U,D,V,3,3,3,3,E);

    free( W); free( U); free(V); free(D);
    free(E1); free(E2);
}
/* triangulate by chieral-method--------------------------------------------
 * args   :  frame_t *frame  I  input feature points
 *           double *K       I  camera calibration matrix
 *           double *R       I  rotation matrix
 *           double *t       I  translation matrix
 *           double *X       O  output triangulate points
 * return : number of inliers of triangular points
 * --------------------------------------------------------------------------*/
static int triangulate(const frame_t *frame,const double *K,const double *R,
                       const double *t,double *X)
{
    register int i,j,n;
    double *P1,*P2,*P3,*U,*W,*V,*J,*AX1,*BX1;

    /* number of matched feature points */
    n=frame->n;

    P1=zeros(3,4);
    P2=zeros(3,4); P3=zeros(3,4);
    U=mat(4,4); W=mat(4,4);
    V=mat(4,4); J=mat(4,4);

    AX1=mat(3,n); BX1=mat(3,n);

    asi_blk_mat(P1,3,4,K,3,3,0,0);
    asi_blk_mat(P2,3,4,R,3,3,0,0);
    asi_blk_mat(P2,3,4,t,3,1,0,3);

    matmul("NN",3,4,3,1.0,K,P2,0.0,P3);

    /* triangulation via orthogonal regression */
    for (i=0;i<frame->n;i++) {
        for (j=0;j<4;j++) {
            J[0+j*4]=P1[2+3*j]*frame->data[i].u1p-P1[0+3*j];
            J[1+j*4]=P1[2+3*j]*frame->data[i].v1p-P1[1+3*j];
            J[2+j*4]=P3[2+3*j]*frame->data[i].u1c-P3[0+3*j];
            J[3+j*4]=P3[2+3*j]*frame->data[i].v1c-P3[1+3*j];
        }
        /* singular value decomposition */
        svd(J,4,4,U,W,V);

        /* triangulation point */
        asi_blk_mat(X,4,n,V+12,4,1,0,i);
    }
    /* compute inliers */
    matmul("NN",3,n,4,1.0,P1,X,0.0,AX1);
    matmul("NN",3,n,4,1.0,P3,X,0.0,BX1);

    for (j=0,i=0;i<n;i++) {
        if (AX1[2+i*3]*X[3+i*4]>0&&BX1[2+i*3]*X[3+i*4]>0) j++;
    }
    free(P1); free(P2); free(P3);
    free(U); free(W); free(V);
    free(J);
    return j; /* return number of inliers */
}
/* matrix copy--------------------------------------------------------------*/
extern void matcop(const double *A,int m,int n,const double s,double *B)
{
    int i,j; for (i=0;i<m;i++) for (j=0;j<n;j++) B[i+j*m]=s*A[i+j*m];
}
/* get R|t------------------------------------------------------------------*/
static void getRt(const double *Ra,const double *Rb,const double *t,int i,
                  double *Ro,double *to)
{
    if (i==0) {matcpy(Ro,Ra,3,3); matcpy(to,t,1,3);}
    if (i==1) {matcpy(Ro,Ra,3,3); matcop(t,1,3,-1.0,to);}
    if (i==2) {matcpy(Ro,Rb,3,3); matcpy(to,t,1,3);}
    if (i==3) {matcpy(Ro,Rb,3,3); matcop(t,1,3,-1.0,to);}
}
/* compute 3d points X and R|t up to scale----------------------------------*/
static int e2rt(const frame_t *frame,const double *E,const double *K,double *X,
                double *R,double *t)
{
    register int i,n,m,max=0;
    double W[9]={0,1,0,-1,0,0,0,0,1},Z[9]={0,-1,0,1,0,0,0,0,0};
    double *U,*S,*V,*T,*Ra,*Rb,*Xc,*Ro,*to;

    n=frame->n;

    Xc=mat(4,n); Ro=mat(3,3); to=mat(1,3);
    U =mat(3,3); S =mat(1,3); V =mat(3,3);
    T =mat(3,3); Ra=mat(3,3); Rb=mat(3,3);

    /* extract T,R1,R2 (8 solutions) */
    svd(E,3,3,U,S,V);

    matmul33("NNT",U,Z,U,3,3,3,3,T);
    matmul33("NNT",U,W,V,3,3,3,3,Ra);
    matmul33("NTT",U,W,V,3,3,3,3,Rb);

    /* convert T to t */
    t[0]=T[5]; t[1]=T[6]; t[2]=T[1];

    /* assure determinant to be positive */
    if (det(Ra,3)<0) for (i=0;i<9;i++) Ra[i]=-Ra[i];
    if (det(Rb,3)<0) {
        for (i=0;i<9;i++) Rb[i]=-Rb[i];
    }
    /* create vector containing all 4 solutions */
    for (i=0;i<4;i++) {
        getRt(Ra,Rb,t,i,Ro,to);
        m=triangulate(frame,K,Ro,to,Xc);
        if (m>max) {
            max=m;
            matcpy(X,Xc,4,n); matcpy(R,Ro,3,3); matcpy(t,to,1,3);
        }
    }
    free(Ra); free(Rb); free(U);
    free(Ro); free(to);
    free(S ); free(V );
    free(T ); free(Xc);

    /* return number of inliers */
    return max;
}
/* compare distance ---------------------------------------------------------*/
static int cmpdist(const void *p1, const void *p2)
{
    dist_t *q1=(dist_t *)p1,*q2=(dist_t *)p2;
    double ds=q1->d-q2->d;
    if (fabs(ds)>DTTOL) return ds<0?-1:1;
}
/* smaller than median-------------------------------------------------------*/
static int smallerthanmedian(const double *X,int n,double *mid,double *Xs)
{
    register int i; dist_t *dist;

    dist=(dist_t*)malloc(sizeof(dist_t)*n);

    for (i=0;i<n;i++) {
        dist[i].d=fabs(X[0+4*i])+fabs(X[1+4*i])+fabs(X[2+4*i]);
        dist[i].id=i;
    }
    /* sort elements */
    qsort(dist,n,sizeof(dist_t),cmpdist);

    /* get median */
    *mid=dist[n/2].d;

    /* create matrix containing elements closer than median */
    if (Xs) {
        for (i=0;i<=n/2;i++) matcpy(Xs+4*i,X+4*dist[i].id,1,4);
    }
    free(dist);
    return (n/2+1);
}
/* get inlier feature points from frame--------------------------------------
 * args   :  frame_t *frame  I  input image frame
 *           voaid_t *opt    I  visual odometry options
 *           double *F       I  fundamental matrix
 *           int *index      O  index of inlier feature points
 * return : number of inlier feature points
 * --------------------------------------------------------------------------*/
static int getinlier(const frame_t *frame,const voopt_t *opt,const double *F,
                     int *index)
{
    register int i,ni;
    double Fx[3],Ftx[3],x1[3],x2[3],x2tFx1,d;

    for (ni=0,i=0;i<frame->n;i++) {

        /* F*x1 */
        x1[0]=frame->data[i].u1p;
        x1[1]=frame->data[i].v1p; x1[2]=1.0;

        matmul("NN",3,1,3,1.0,F,x1,0.0,Fx);

        /* F'*x2 */
        x2[0]=frame->data[i].u1c;
        x2[1]=frame->data[i].v1c; x2[2]=1.0;

        matmul("TN",3,1,3,1.0,F,x2,0.0,Ftx);

        /* x2'*F*x1 */
        matmul33("TNN",x2,F,x1,1,3,3,1,&x2tFx1);

        /* sampson distance */
        d=SQR(x2tFx1)/(SQR(Fx[0])+SQR(Fx[1])+SQR(Ftx[0])+SQR(Ftx[1]));

        /* check threshold */
        if (fabs(d)<opt->inlier_thres) {
            index[ni++]=i;
        }
    }
    return ni;
}
/* get random and unique sample of num numbers from 1:N ---------------------
 * args   :  int n        I  number of feature points
 *           int num      I  number of sample feature sample
 *           int *sample  O  index list of sample feature points
 * return : number of sampled feature points
 * --------------------------------------------------------------------------*/
static int getrandsample(int n,int num,int *sample)
{
    int i,j,ns,*list;
    list=imat(n,1);

    /* create vector containing all indices */
    for (i=0;i<n;i++) list[i]=i;

    /* add num indices to current sample */
    for (ns=0,i=0;ns<num&&i<n;i++) {
        j=(int)random()%n;
        if (list[j]<0) continue;
        sample[ns++]=j; /* add sample index */
        list[j]=-1; 
    }
    free(list);
    return ns;
}
/* add feature point to frame------------------------------------------------*/
static int addfeature(frame_t *frame,const fea_t *p)
{
    fea_t *fea_data;
    if (frame->nmax<=frame->n) {

        if (frame->nmax<=0) frame->nmax=64;
        else frame->nmax+=256;

        if (!(fea_data=(fea_t *)realloc(frame->data,sizeof(fea_t)*frame->nmax))) {
            trace(1,"addfeature: add feature fail\n");
            free(frame->data);
            frame->data=NULL; frame->n=frame->nmax=0;
            return -1;
        }
        frame->data=fea_data;
    }
    frame->data[frame->n++]=*p;
    return 1;
}
/* free frame struct---------------------------------------------------------*/
static void freeframe(frame_t *frame)
{
    if (frame->data) free(frame->data);
    frame->n=frame->nmax=0;
}
/* normalize feature points--------------------------------------------------
 * args   :  frame_t *frame  I  input image feature points
 *           voaid_t *opt    I  visual odometry options
 *           double *Tp      O  precious image transformation matrices
 *           double *Tc      O  current image transformation matrices
 *           frame_t *nframe O  output normalize feature points
 * retutn : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
static int normfeature(const frame_t *frame,const voopt_t *opt,
                       double *Tp,double *Tc,frame_t *nframe)
{
    register int i;
    double cpu=0,cpv=0,ccu=0,ccv=0;
    double sp=0,sc=0;
    fea_t p;

    trace(3,"normfeature:\n");

    if (frame->n<=0) {
        trace(2,"no feature points to normalize\n");
        return 0;
    }
    /* shift origins to centroids */
    for (i=0;i<frame->n;i++) {
        cpu+=frame->data[i].u1p; cpv+=frame->data[i].v1p;
        ccu+=frame->data[i].u1c; ccv+=frame->data[i].v1c;
    }
    cpu/=frame->n;
    cpv/=frame->n;
    ccu/=frame->n;
    ccv/=frame->n;

    for (i=0;i<frame->n;i++) {
        p.u1p=frame->data[i].u1p-cpu; p.v1p=frame->data[i].v1p-cpv;
        p.u1c=frame->data[i].u1c-ccu; p.v1c=frame->data[i].v1c-ccv;
        addfeature(nframe,&p);
    }
    /* scale features such that mean distance from origin is sqrt(2) */
    for (i=0;i<frame->n;i++) {
        sp+=SQRT(SQR(nframe->data[i].u1p)+SQR(nframe->data[i].v1p));
        sc+=SQRT(SQR(nframe->data[i].u1c)+SQR(nframe->data[i].v1c));
    }
    if (fabs(sp)<1E-10||fabs(sc)<1E-10) return 0;

    sp=SQRT(2.0)*frame->n/sp;
    sc=SQRT(2.0)*frame->n/sc;
    for (i=0;i<frame->n;i++) {
        nframe->data[i].u1p*=sp; nframe->data[i].v1p*=sp;
        nframe->data[i].u1c*=sc; nframe->data[i].v1c*=sc;
    }
    /* transform matrix */
    if (Tp) {
        Tp[0]=sp;      Tp[1]=0.0;     Tp[2]=0.0;
        Tp[3]=0.0;     Tp[4]=sp;      Tp[5]=0.0;
        Tp[6]=-sp*cpu; Tp[7]=-sp*cpv; Tp[8]=1.0;
    }
    if (Tc) {
        Tc[0]=sc;      Tc[1]=0.0;     Tc[2]=0.0;
        Tc[3]=0.0;     Tc[4]=sc;      Tc[5]=0.0;
        Tc[6]=-sc*ccu; Tc[7]=-sc*ccv; Tc[8]=1.0;
    }
    return 1;
}
/* check match feature geometric distribution---------------------------------*/
static int chkgeodistribt(const voopt_t *opt,const frame_t *frame,const int *active,int n)
{
    double xyz[3],*Q,*H,gdop=1E9;
    int nv=0,i,j;

    H=zeros(3,n); Q=zeros(3,3);
    for (i=0;i<n;i++) {
        xyz[0]=frame->data[active[i]].u1c;
        xyz[1]=frame->data[active[i]].v1c;
        xyz[2]=1.0;

        for (j=0;j<3;j++) H[3*nv+j]=xyz[j]/norm(xyz,3);
        nv++;
    }
    matmul("NT",3,3,n,1.0,H,H,0.0,Q);
    if (!matinv(Q,3)) {
        gdop=SQRT(Q[0]+Q[4]+Q[8]); /* gdop */
    }
    free(H); free(Q);
    return gdop<MAXVOGDOP;
}
/* estimate motion-----------------------------------------------------------
 * args   :  voopt_t *opt    I  visual odometry options
 *           frame_t *frame  I  feature points frame
 *           double *Tr      O  rotation and translation
 *           double *ratio   O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
static int estmono(const voopt_t *opt,const frame_t *frame,double *Tr,double *ratio)
{
    int i,j,n,*idx,*in,*inl,ni,m;
    double K[9],Tp[9],Tc[9],E[9],F[9],U[9],W[3],V[9],D[9];
    double R[9],t[3],*X,*x,*d,mid;
    frame_t nf={0};

    trace(3,"estmono:\n");

    if ((n=frame->n)<8) {
        trace(2,"no enough feature points\n");
        return 0;
    }
    idx=imat(8,1); in=imat(1,n); inl=imat(1,n);
    X=mat(4,n);
    x=mat(2,n);
    d=mat(1,n);

    /* create calibration matrix */
    getK(opt,K);

    /* normalize feature points and return on errors */
    if (!normfeature(frame,opt,Tp,Tc,&nf)) {

        trace(2,"normalize fail\n");
        free(idx); free(in); free(X); free(x);
        free(inl); free(d );
        freeframe(&nf);
        return 0;
    }
    /* initial RANSAC estimate of F */
    for (ni=0,i=0;i<opt->ransac_iters;i++) {

        /* draw random sample set */
        getrandsample(n,8,idx);

        if (!chkgeodistribt(opt,&nf,idx,8)) continue;

        /* estimate fundamental matrix and get inliers */
        fundamental(&nf,idx,8,F); m=getinlier(&nf,opt,F,inl);

        /* update model if we are better */
        if (m>ni) {
            ni=m;
            imatcpy(in,inl,1,m);
        }
    }
    if (ni<10) {
        trace(2,"no enough inliers\n");
        free(idx); free(in); free(X); free(x);
        free(inl); free(d );
        freeframe(&nf);
        return 0;
    }
    /* ratio of inliers */
    if (ratio) *ratio=(double)ni/nf.n;

#if REFINE_F
    /* refine F using all inliers */
    fundamental(&nf,in,ni,F);
#endif
    /* de-normalise and extract essential matrix */
    matmul33("TNN",Tc,F,Tp,3,3,3,3,F);
    matmul33("TNN",K,F,K,3,3,3,3,E);

    /* re-enforce rank 2 constraint */
    svd(E,3,3,U,W,V); W[2]=0.0;
    dialog(W,3,D);
    matmul33("NNT",U,D,V,3,3,3,3,E);

    /* compute 3d points X and R|t up to scale */
    e2rt(frame,E,K,X,R,t);

    /* normalize 3D points and remove points behind image plane */
    for (i=0;i<n;i++) {
        for (j=0;j<4;j++) X[4*i+j]/=X[4*i+3];
    }
    for (i=0,j=0;i<n;i++) {
        if (X[2+i*4]>0) matcpy(X+4*j++,X+i*4,1,4);
    }
    if (j<10) {
        trace(2,"no enough feature points to process\n");

        free(idx); free(in); free(X); free(x);
        free(inl); free(d );
        freeframe(&nf);
        return 0;
    }
#if 0
    smallerthanmedian(X,j,&mid,NULL);

    if (mid>opt->motion_thres) {
        trace(2,"fail due to little motion\n");
        free(idx); free(in); free(X); free(x);
        free(inl); free(d );
        freeframe(&nf);
        return 0;
    }
#endif
    /* compute rotation angles */
    U[0]=asin( R[0+3*2]);
    U[1]=asin(-R[1+3*2]/cos(U[0]));
    U[2]=asin(-R[0+3*1]/cos(U[0]));

    if (Tr) {
        Tr[0]=U[0]; Tr[1]=U[1]; Tr[2]=U[2];
        Tr[3]=t[0]; Tr[4]=t[1]; Tr[5]=t[2];
    }
    free(idx); free(in); free(X); free(x);
    free(inl); free(d );
    freeframe(&nf);
    return 1;
}
/* compute transformation matrix from transformation vector------------------*/
static void tfvec2mat(const double *tr,double *T)
{
    double rx,ry,rz,tx,ty,tz,sx,cx,sy,cy,sz,cz;

    rx=tr[0]; ry=tr[1]; rz=tr[2];
    tx=tr[3]; ty=tr[4]; tz=tr[5];

    /* precompute sine/cosine */
    sx=sin(rx),cx=cos(rx),sy=sin(ry);
    cy=cos(ry),sz=sin(rz),cz=cos(rz);

    /* compute transformation */
    T[0]=+cy*cz;          T[4]=-cy*sz;          T[8]=+sy;     T[12]=tx;
    T[1]=+sx*sy*cz+cx*sz; T[5]=-sx*sy*sz+cx*cz; T[9]=-sx*cy;  T[13]=ty;
    T[2]=-cx*sy*cz+sx*sz; T[6]=+cx*sy*sz+sx*cz; T[10]=+cx*cy; T[14]=tz;
    T[3]=0.0;             T[7]=0.0;             T[11]=0.0;    T[15]=1.0;
}
/* estimate mono visual odometry motion(R|t)---------------------------------
 * args   :  voopt_t *opt      I  visual odometry options
 *           match_set_t *mf   I  feature list
 *           double *Tr        O  rotation and translation
 *           double *ratio     O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int estmonort(const voopt_t *opt,const match_set_t *mf,double *Tr,double *ratio)
{
    frame_t frame={0};
    fea_t f={0};
    double tr[6];
    int i,flag;

    trace(3,"estmonort:\n");

    /* do for each match feature point */
    for (i=0;i<mf->n;i++) {
        f.u1c=mf->data[i].uc; f.v1c=mf->data[i].vc;
        f.u1p=mf->data[i].up; f.v1p=mf->data[i].vp;

        /* add new feat. */
        if (!addfeature(&frame,&f)) continue;
    }
    if (!(flag=estmono(opt,&frame,tr,ratio))) {
        trace(3,"mono visual odometry motion estimate fail\n");
        freeframe(&frame);
        return 0;
    }
    tfvec2mat(tr,Tr);
    trace(3,"Tr=\n"); tracemat(3,Tr,4,4,12,6);

    freeframe(&frame);
    return flag;
}
/* ---------------------------------------------------------------------------
 * triangulate Triangulates 3D points from two sets of feature vectors and a
 * a frame-to-frame transformation
 * args:    double *C21    I  rotation matrix from 2-frame to 1-frame
 *          double *t12_1  I  translation from 1-frame to 2-frame expressed
 *                            in 1-frame
 *          double *obs1   I  (u,v) observation in 1-frame
 *          double *obs2   I  (u,v) observation in 2-frame
 *          double *K      I  camera calibration matrix
 *          double *X      O  output feature 3D-point
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int triangulate3D(const double *C21,const double *t12_1,
                         const double *obs1, const double *obs2,
                         const double *K,double *X)
{
    int flag=0,i; double Xp[4];
    frame_t frame={0};
    fea_t feat;

    trace(3,"triangulate3D:\n");

    /* feature point measurement data */
    feat.u1p=obs1[0]; feat.v1p=obs1[1];
    feat.u1c=obs2[0];
    feat.v1c=obs2[1];

    addfeature(&frame,&feat);
    flag=triangulate(&frame,K,C21,t12_1,Xp);

    /* normalize 3D points */
    for (i=0;i<3;i++) X[i]=Xp[i]/Xp[3];

    freeframe(&frame);
    return flag&&X[2]>0.0;
}


