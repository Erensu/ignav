/*---------------------------------------------------------------------------
* vo-stereo.cc : stereo camera estimate motion functions
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/12/27 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

#define UPDATED                  1
#define CONVERGED                2
#define FAILED                   3

typedef struct {                 /* bucketing parameters */
    int max_features;            /* maximal number of features per bucket */
    float  bucket_width;         /* width of bucket */
    float  bucket_height;        /* height of bucket */
} bucket_t;

typedef struct {                 /* visual odometry options */
    match_t match;               /* feature match options */
    bucket_t bucket;             /* bucketing parameters */
    calib_t calib;               /* camera parameters */
    double height;               /* camera height above ground (meters) */
    double pitch;                /* camera pitch (rad, negative=pointing down) */
    double base;                 /* baseline (meters) for stereo camera */
    int    ransac_iters;         /* number of RANSAC iterations */
    double inlier_thres;         /* fundamental matrix inlier threshold */
    double motion_thres;         /* directly return false on small motions */
    bool   reweighting;          /* lower border weights (more robust to calibration errors) */
} vopt_t;

typedef struct {                 /* structure for storing matches */
    double u1p,v1p;              /* u,v-coordinates in previous left  image */
    double u2p,v2p;              /* u,v-coordinates in previous right image */
    double u1c,v1c;              /* u,v-coordinates in current  left  image */
    double u2c,v2c;              /* u,v-coordinates in current  right image */
    int i1p;                     /* feature index (for tracking) */
    int i2p;                     /* feature index (for tracking) */
    int i1c;                     /* feature index (for tracking) */
    int i2c;                     /* feature index (for tracking) */
} fea_t;

typedef struct {                 /* structure for storing a frame image */
    gtime_t gtime;               /* time stamp */
    int n,nmax;                  /* number of feature data/allocated */
    fea_t *data;                 /* feature data records */
} frame_t;

/* get stereo camera observation data----------------------------------------
 * args   :  frame_t *frame  I  frame observation data
 *           int *active     I  active index of feature points
 *           int n           I  number of active feature points
 *           double *pobs    O  output feature points observation
 * return : none
 * --------------------------------------------------------------------------*/
static void getobs(const frame_t *frame,const int *active,int n,double *pobs)
{
    int i;

    for (i=0;i<n;i++) {
        pobs[4*i+0]=frame->data[active[i]].u1c; /* u1 */
        pobs[4*i+1]=frame->data[active[i]].v1c; /* v1 */
        pobs[4*i+2]=frame->data[active[i]].u2c; /* u2 */
        pobs[4*i+3]=frame->data[active[i]].v2c; /* v2 */
    }
}
/* get residuals and jacobian of feature point observation data--------------
 * args   :  double *tr   I  motion parameters
 *           int *active  I  active index of feature points
 *           int n        I  number of active feature points
 *           vopt_t *vopt I  visual odometry option
 *           double *XYZ  I  3D points
 *           double *pobs I  feature point observation data
 *           double *J    O  jacobian of feature points observation
 *           double *pres O  predict residuals
 *           double *prex O  predict feature points
 * return : none
 * --------------------------------------------------------------------------*/
static void res_jacobian(const double *tr,const int *active,int n,
                         const vopt_t *vopt,const double *XYZ,const double *pobs,
                         double *J,double *prex,double *pres)
{
    int i,j;
    double rx=tr[0],ry=tr[1],rz=tr[2];
    double tx=tr[3],ty=tr[4],tz=tr[5];

    double sx=sin(rx),cx=cos(rx),sy=sin(ry);
    double cy=cos(ry),sz=sin(rz),cz=cos(rz);

    /* compute rotation matrix and derivatives : R=Rx*Ry*Rz */
    double r00=+cy*cz;          double r01=-cy*sz;          double r02= +sy;
    double r10=+sx*sy*cz+cx*sz; double r11=-sx*sy*sz+cx*cz; double r12=-sx*cy;
    double r20=-cx*sy*cz+sx*sz; double r21=+cx*sy*sz+sx*cz; double r22=+cx*cy;

    double rdrx10=+cx*sy*cz-sx*sz; double rdrx11=-cx*sy*sz-sx*cz; double rdrx12=-cx*cy;
    double rdrx20=+sx*sy*cz+cx*sz; double rdrx21=-sx*sy*sz+cx*cz; double rdrx22=-sx*cy;
    double rdry00=-sy*cz;          double rdry01=+sy*sz;          double rdry02=+cy;
    double rdry10=+sx*cy*cz;       double rdry11=-sx*cy*sz;       double rdry12=+sx*sy;
    double rdry20=-cx*cy*cz;       double rdry21=+cx*cy*sz;       double rdry22=-cx*sy;
    double rdrz00=-cy*sz;          double rdrz01=-cy*cz;
    double rdrz10=-sx*sy*sz+cx*cz; double rdrz11=-sx*sy*cz-cx*sz;
    double rdrz20=+cx*sy*sz+sx*cz; double rdrz21=+cx*sy*cz-sx*sz;

    /* loop variables */
    double X1p,Y1p,Z1p;
    double X1c,Y1c,Z1c,X2c;
    double X1cd,Y1cd,Z1cd;
    double w=1.0;

    /* for all observations do */
    for (i=0;i<n;i++) {

        /* get 3d point in previous coordinate system */
        X1p=XYZ[3*active[i]+0];
        Y1p=XYZ[3*active[i]+1];
        Z1p=XYZ[3*active[i]+2];

        /* compute 3d point in current left coordinate system */
        X1c=r00*X1p+r01*Y1p+r02*Z1p+tx;
        Y1c=r10*X1p+r11*Y1p+r12*Z1p+ty;
        Z1c=r20*X1p+r21*Y1p+r22*Z1p+tz;

        /* weighting */
        if (vopt->reweighting) {
            w=1.0/(fabs(pobs[4*i+0]-vopt->calib.cu)/fabs(vopt->calib.cu)+0.05);
        }
        /* compute 3d point in current right coordinate system */
        X2c=X1c-vopt->base;

        /* for all paramters do (rx,ry,rz,tx,ty,tz)*/
        for (j=0;j<6;j++) {

            /* derivatives of 3d pt. in curr. left coordinates wrt. param j */
            switch (j) {
                case 0:
                    X1cd=0;
                    Y1cd=rdrx10*X1p+rdrx11*Y1p+rdrx12*Z1p;
                    Z1cd=rdrx20*X1p+rdrx21*Y1p+rdrx22*Z1p;
                    break;
                case 1:
                    X1cd=rdry00*X1p+rdry01*Y1p+rdry02*Z1p;
                    Y1cd=rdry10*X1p+rdry11*Y1p+rdry12*Z1p;
                    Z1cd=rdry20*X1p+rdry21*Y1p+rdry22*Z1p;
                    break;
                case 2:
                    X1cd=rdrz00*X1p+rdrz01*Y1p;
                    Y1cd=rdrz10*X1p+rdrz11*Y1p;
                    Z1cd=rdrz20*X1p+rdrz21*Y1p;
                    break;
                case 3: X1cd=1; Y1cd=0; Z1cd=0; break;
                case 4: X1cd=0; Y1cd=1; Z1cd=0; break;
                case 5: X1cd=0; Y1cd=0; Z1cd=1; break;
            }
            if (J) {
                /* set jacobian entries (project via K) */
                J[(4*i+0)*6+j]=w*vopt->calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c); /* left u' */
                J[(4*i+1)*6+j]=w*vopt->calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c); /* left v' */
                J[(4*i+2)*6+j]=w*vopt->calib.f*(X1cd*Z1c-X2c*Z1cd)/(Z1c*Z1c); /* right u' */
                J[(4*i+3)*6+j]=w*vopt->calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c); /* right v' */
            }
        }
        /* set prediction (project via K) */
        if (prex) {
            prex[4*i+0]=vopt->calib.f*X1c/Z1c+vopt->calib.cu; /* left u */
            prex[4*i+1]=vopt->calib.f*Y1c/Z1c+vopt->calib.cv; /* left v */
            prex[4*i+2]=vopt->calib.f*X2c/Z1c+vopt->calib.cu; /* right u */
            prex[4*i+3]=vopt->calib.f*Y1c/Z1c+vopt->calib.cv; /* right v */
        }
        if (pres&&prex) {
            /* set residuals */
            pres[4*i+0]=w*(pobs[4*i+0]-prex[4*i+0]);
            pres[4*i+1]=w*(pobs[4*i+1]-prex[4*i+1]);
            pres[4*i+2]=w*(pobs[4*i+2]-prex[4*i+2]);
            pres[4*i+3]=w*(pobs[4*i+3]-prex[4*i+3]);
        }
    }
}
/* update stereo motion parameters-------------------------------------------
 * args   :  frame_t *frame  I   frame observation data
 *           vopt_t *vopt    I   visual odometry option
 *           int *active     I   active index of feature points
 *           int n           I   number of active feature points
 *           double *tr      IO  motion parameters
 *           double step     I   estimate parameters step
 *           double eps      I   thresholds for estimating parameters
 * return : 0: fail, 1: converged, 2: update
 * --------------------------------------------------------------------------*/
static int updparam(const frame_t *frame,const vopt_t *vopt,const int *active,
                    const double *XYZ,int na,double *tr,double step,double eps)
{
    int i,j,k,info=0;
    double A[6*6],B[6],C[6],*J,*px,*pr,*obs,a,b;

    /* we need at least 3 observations */
    if (na<3) return 0;

    obs=mat(4,na); px=mat(4,na); pr=mat(4,na);
    J=mat(4,6*na);

    /* extract observations and compute predictions */
    getobs(frame,active,na,obs);
    res_jacobian(tr,active,na,vopt,XYZ,obs,J,px,pr);

    /* fill matrices A and B */
    for (i=0;i<6;i++) {
        for (j=0;j<6;j++) {
            a=0.0;
            for (k=0;k<4*na;k++) a+=J[k*6+i]*J[k*6+j];
            A[6*i+j]=a;
        }
        b=0.0;
        for (k=0;k<4*na;k++) {
            b+=J[k*6+i]*(pr[k]);
        }
        B[i]=b;
    }
    /* perform elimination */
    if (!solve("T",A,B,6,1,C)) {
        info=1;
        for (i=0;i<6;i++) {
            tr[i]+=step*C[i];
            if (fabs(C[i])>eps) info=0;
        }
        /* free local variables */
        free(obs);
        free(px); free(pr); free(J);

        if (info) return CONVERGED; /* converged */
        else      return UPDATED; /* update */
    }
    /* free local variables */
    free(obs); free(px);
    free(pr); free(J);
    return FAILED;
}
/* get inliers of feature points --------------------------------------------*/
static int getInlier(const frame_t *frame,const double *XYZ,const vopt_t *vopt,
                     const double *tr,int *inliers)
{
    int i,j,n,*active;
    double *obs,*J,*px,*pr;

    if ((n=frame->n)<=0) {
        trace(2,"no feature points\n");
        return 0;
    }
    active=imat(n,1);
    obs=mat(4,n); J=mat(4,6*n);
    px=mat(4,n); pr=mat(4,n);

    /* index of active feature points */
    for (i=0;i<n;i++) active[i]=i;

    /* extract observations and compute predictions */
    getobs(frame,active,n,obs);
    res_jacobian(tr,active,n,vopt,XYZ,obs,J,px,pr);

    /* compute inliers */
    for (j=0,i=0;i<n;i++) {
        if (pow(obs[4*i+0]-px[4*i+0],2)+pow(obs[4*i+1]-px[4*i+1],2)+
            pow(obs[4*i+2]-px[4*i+2],2)+pow(obs[4*i+3]-px[4*i+3],2)
            <vopt->inlier_thres*vopt->inlier_thres) {
            inliers[j++]=i;
        }
    }
    free(obs); free(J);
    free(px); free(pr);
    return j;
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

    trace(3,"getrandsample:\n");

    list=imat(n,1);

    /* create vector containing all indices */
    for (i=0;i<n;i++) list[i]=i;

    /* add num indices to current sample */
    for (ns=0,i=0;ns<num&&i<n;i++) {
        j=rand()%n;
        if (list[j]<0) continue;
        sample[ns++]=j; /* add sample index */
        list[j]=-1; /* disable this feature point */
    }
    free(list);
    return ns;
}
/* stereo camera estimate motion---------------------------------------------
 * args   :  vopt_t *opt     I  visual odometry options
 *           frame_t *frame  I  feature points frame
 *           double *Tr      O  rotation and translation
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
static int eststereo(const vopt_t *opt,const frame_t *frame,double *Tr)
{
    double d,*XYZ,Trc[6];
    int info=1,i,n,iter,rslt,*idx,m=0,mp=0,*inliers,*inlierc;

    trace(3,"eststereo:\n");

    if ((n=frame->n)<=0) {
        trace(2,"no feature points\n");
        return 0;
    }
    XYZ=mat(3,n);
    idx=imat(1,n);
    inliers=imat(1,n); inlierc=imat(1,n);

    /* project matches of previous image into 3d */
    for (i=0;i<frame->n;i++) {

        d=MAX(frame->data[i].u1p-frame->data[i].u2p,0.0001);

        XYZ[3*i+0]=(frame->data[i].u1p-opt->calib.cu)*opt->base/d;
        XYZ[3*i+1]=(frame->data[i].v1p-opt->calib.cv)*opt->base/d;
        XYZ[3*i+2]=opt->calib.f*opt->base/d;
    }
    /* initial RANSAC estimate */
    for (i=0;i<opt->ransac_iters;i++) {

        /* draw random sample set */
        getrandsample(n,3,idx);

        /* minimize reprojection errors */
        iter=0;
        rslt=UPDATED;
        while (rslt==UPDATED) {
            rslt=updparam(frame,opt,idx,XYZ,3,Trc,1,1E-6);
            if (iter++>20||rslt==CONVERGED) break;
        }
        /* overwrite best parameters if we have more inliers */
        if (rslt!=FAILED) {
            m=getInlier(frame,XYZ,opt,Trc,inlierc);
            if (m>mp) {
                imatcpy(inliers,inlierc,1,m);
                matcpy(Tr,Trc,1,6);
                mp=m;
            }
        }
    }
    /* final optimization (refinement) */
    if (mp>6) {
        iter=0;
        rslt=UPDATED;
        while (rslt==UPDATED) {
            rslt=updparam(frame,opt,idx,XYZ,3,Tr,1,1E-6);
            if (iter++>20||rslt==CONVERGED) break;
        }
        /* not converged */
        if (rslt!=CONVERGED) info=0;
    }
    else {
        /* no enough inliers */
        info=0;
    }
    free(XYZ);
    free(idx);
    free(inliers); free(inlierc);
    return info;
}
