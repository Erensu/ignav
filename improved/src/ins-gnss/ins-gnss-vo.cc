/*--------------------------------------------------------------------------------
 * ins-gnss-vo.cc : ins-gnss-vo coupled common functions
 *
 * reference :
 *    [01] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *         Navigation System, Artech House, 2008
 *    [02] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *         for IMU calibration without external equipments,2014.
 *    [03] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *         INS 2008.
 *    [04] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [05] Li M, Mourikis A I. High-precision, consistent EKF-based visual–inertial
 *         odometry[J].International Journal of Robotics Research,2013,32(6):690-711.
 *    [06] Monocular Visual Inertial Odometry on a Mobile Device.
 *    [07] Mourikis A / , Roumeliotis S / . A Multi-State Constraint Kalman Filter
 *         for Vision-aided Inertial Navigation[C]// IEEE International Conference
 *         on Robotics & Automation. IEEE, 2007.
 *    [08] Forster C , Carlone L , Dellaert F , et al. On-Manifold Preintegration
 *         for Real-Time Visual-Inertial Odometry[J]. IEEE Transactions on Robotics,
 *         2015, 33(1):1-21.
 *    [09] Observability-constrained vision-aided inertial navigation.
 *    [10] Li M ,Mourikis A I. Improving the accuracy of EKF-based visual-inertial
 *         odometry[C]// IEEE International Conference on Robotics & Automation.
 *         IEEE, 2012.
 *    [11] Li M ,Mourikis A I. High-precision, consistent EKF-based visual-inertial
 *         odometry[M].Sage Publications, Inc. 2013.
 *    [12] Pizzoli M , Forster C , Scaramuzza D . REMODE: Probabilistic, Monocular
 *         Dense Reconstruction in Real Time[C]// IEEE International Conference on
 *         Robotics and Automation (ICRA), Hong Kong, 2014. IEEE, 2014.
 *    [13] Li M , Kim B H , Mourikis A I . Real-time motion tracking on a cellphone
 *         using inertial sensing and a rolling-shutter camera[C]// IEEE International
 *         Conference on Robotics & Automation. IEEE, 2013.
 *    [14] Monocular Visual-Inertial SLAM and Self Calibration for Long Term Autonomy
 *    [15] Lupton T , Sukkarieh S . Visual-Inertial-Aided Navigation for High-Dynamic
 *         Motion in Built Environments Without Initial Conditions[J]. IEEE Transactions
 *         on Robotics, 2012, 28(1):61-76.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define NUMPOSE              10   /* numbers of camera transform matrix for detect motion*/
#define THRES_ANG     (3.0*D2R)   /* thresholds of angular rate for detect motion */
#define THRES_TRANS       (0.1)   /* thresholds of translation for detect motion */
#define MAXVEL              0.1   /* max velocity for using non-holonomic constraint */
#define MINVEL              5.0   /* min velocity for initial camera pose */
#define MAXGYRO      (10.0*D2R)   /* max rotation speed for using non-holonomic constraint */
#define VARVEL        SQR(0.05)   /* initial variance of receiver vel ((m/s)^2) */

#define MIN_LEN        3          /* min length of tracking data */
#define MAX_LEN        20         /* max length of tracking data */
#define MAXRES_POSE    30.0*D2R   /* max residual for pose filter (rad) */
#define VARPOSE        1.5*D2R    /* variance of pose measurement (rad) */
#define VARPOS         3.0        /* variance of position measurement (m) */
#define THRES_RATIO    0.6        /* threshold of ratio inliers to update camera motion */
#define THRES_TRACK_RATIO  0.4    /* threshold of ratio for check camera track status */
#define USE_BUCKET_FEATS 1        /* use bucket matched feature points data */

#define OPENCV           0
#define LIBVISO          0
#define GMS              0

/* type definitions-----------------------------------------------------------*/
typedef struct {
    gtime_t gtime;                      /* pose measurement time */
    int status;                         /* visual odometry estimate status */
    double rpy[3];                      /* roll/pitch/yaw of precious to current frame */
    double tra[3];                      /* translation of precious to current frame */
} pose_t;

typedef struct {                        /* mono visual odometry data */
    gtime_t ts,te;                      /* start/end time of mono visual odometry data */
    double dT[16];                      /* transformation matrix (4x4) */
    double rc[2][3],Cce[2][9];          /* camera position/attitude in ecef */
} voss_t;

typedef struct {                        /* pre-integration IMU measurement data type */
    gtime_t ts,te;                      /* start/end time of mono visual odometry data */
    double dp[3],dC[9];                 /* pre-integration IMU measurement data */
    double Rp[9],Ra[9];                 /* measurement variance */
} preint_t;

typedef struct {                        /* filter workspace */
    double Cbe[2][9],re[2][3],ve[2][3]; /* ins states in precious and current epoch */
    double Cce[2][9],rc[2][3];          /* camera states in precious and current epoch */
    insstate_t *insdata;                /* ins states data */
    voss_t *vodata;                     /* mono visual odometry data */
    preint_t *pimudata;                 /* pre-integration IMU measurement data */
    int ni,nimax;                       /* number of ins states in precious */
    int nv,nvmax;                       /* number of mono visual odometry */
    int np,npmax;                       /* number of pre-integration IMU measurement data */
    int flag;                           /* flag of filter */
} filt_t;

/* global variables------------------------------------------------------------*/
static track_t  tracks={0};             /* all tracking feature points data */
static match_t  matchs={0};             /* match feature points data */
static filt_t   filts ={0};             /* vo aid filter workspace */

/* predict image point using homography---------------------------------------
 * args:    insopt_t *opt   I   ins options
 *          double *R,*t    I   rotation and translation between frames
 *          double *K       I   camera calibration parameters
 *          double *uv      I   image coordinate of precious frame
 *          double *pre_uv  O   predict image coordinate in current frame
 * return: status (1: ok, 0: fail)
 * ---------------------------------------------------------------------------*/
extern int predictfeat(const double *R,const double *t,const double *K,
                       const double *uv,double *uvp)
{
    double H[9],Ki[9],p[3],pp[3],T[16];
    double C[9],r[9];

    trace(3,"predictfeatH:\n");

    /* pose of current relative to precious */
    rt2tf(R,t,T); if (matinv(T,4)) return 0;
    tf2rt(T,C,r);

    /* homography by rotation */
    matcpy(Ki,K,3,3); if (matinv(Ki,3)) return 0;
    matmul33("NNN",K,C,Ki,3,3,3,3,H);

    trace(3,"H=\n"); tracemat(3,H,3,3,12,6);

    /* predict image coordinate */
    p[0]=uv[0];
    p[1]=uv[1];
    p[2]=1.0;
    matmul("NN",3,1,3,1.0,H,p,0.0,pp);

    uvp[0]=pp[0]/pp[2]; uvp[1]=pp[1]/pp[2];
    return 1;
}
/* is in ROI------------------------------------------------------------------*/
extern int inroi(const float u,const float v,const matchopt_t *opt)
{
    return u>opt->roi[0][0]&&u<opt->roi[1][0]&&
           v>opt->roi[0][1]&&v<opt->roi[1][1];
}
/* update camera pose by using visual odometry-------------------------------
 * args  :  insstate_t *ins  IO  ins states
 *          insopt_t *opt    I   ins options
 *          double *dT       I   transform of precious to current frame
 *          gtime_t time     I   time of transform measurement
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int updatecamera(insstate_t *ins,const insopt_t *opt,const double *dT,gtime_t time)
{
    double Tp[16],T[16],TT[16],dt;

    trace(3,"updatecamera:\n");

    if ((dt=timediff(time,ins->vo.time))>10.0) {
        trace(2,"update camera pose failed\n");
        return 0;
    }
    if (fabs(dt)<1E-5) {
        trace(2,"duplicate update\n");
        return 0;
    }
    rt2tf(ins->vo.Cce,ins->vo.rc,Tp);
    matcpy(TT,dT,4,4);
    matinv(TT,4);
    matmul("NN",4,4,4,1.0,Tp,TT,0.0,T);

    tf2rt(T,ins->vo.Cce,ins->vo.rc);
    ins->vo.time=time;
    return 1;
}
/* ins states convert to camera states-----------------------------------------*/
extern void ins2camera(const insstate_t *ins,double *rc,double *Cce)
{
    double T[3];
    int i;
    matmul("NT",3,3,3,1.0,ins->Cbe,ins->Cbc,0.0,Cce);
    matmul("NN",3,1,3,1.0,ins->Cbe,ins->lbc,0.0,T);
    for (i=0;i<3;i++) {
        rc[i]=ins->re[i]+T[i];
    }
}
/* ins states convert to camera states-----------------------------------------*/
extern void camera2ins(const double *Cce,const double *rc,const double *Cbc,
                       const double *lbc,double *Cbe,double *re)
{
    double T[3];
    int i;
    matmul("NN",3,3,3,1.0,Cce,Cbc,0.0,Cbe);
    matmul("NN",3,1,3,1.0,Cbe,lbc,0.0,T);
    for (i=0;i<3&&re;i++) re[i]=rc[i]-T[i];
}
/* initial camera attitude (relative to n-frame) based on ins states----------
 * args  :  insstate_t *ins  IO  ins states
 *          insopt_t *opt    I   ins options
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int initcampose(insstate_t *ins,const insopt_t *opt)
{
    trace(3,"initcampose:\n");
    int i;
    double T[3];

    matmul("NT",3,3,3,1.0,ins->Cbe,ins->vo.Cbc,0.0,ins->vo.Cce);
    matmul("NN",3,1,3,1.0,ins->Cbe,ins->vo.lbc,0.0,T);
    for (i=0;i<3;i++) {
        ins->vo.rc[i]=ins->re[i]+T[i];
    }
    ins->vo.time=ins->time;
    return 1;
}
/* transform matrix convert to pose-------------------------------------------*/
static int convpose(const double *T,pose_t *pose)
{
    double R[9];
    tf2rt(T,R,pose->tra);
    dcm2rpy(R,pose->rpy);
}
/* detect motion by using camera motions--------------------------------------*/
static int chkstatic(const double *dT,int status,gtime_t time)
{
    static pose_t T[NUMPOSE]={0};
    int i;

    for (i=0;i<NUMPOSE-1;i++) T[i]=T[i+1]; convpose(dT,&T[i]);
    T[i].status=status;
    T[i].gtime =time;

    for (i=0;i<NUMPOSE;i++) {
        if (!T[i].gtime.time||T[i].status<=0) return 0;
    }
    for (i=0;i<NUMPOSE;i++) {
        if (norm(T[i].rpy,3)>THRES_ANG||
            norm(T[i].tra,3)>THRES_TRANS) {

            /* detect motion */
            return 0;
        }
    }
    /* detect static */
    return 1;
}
/* zero velocity update based on camera motion detect-------------------------
 * args  :  insopt_t *opt   I  ins options
 *          double *dT      I  transform matrix of precious to current frame
 *          int status      I  flag of visual odometry estimate status
 *          imud_t *imu     I  imu measurement data
 *          gtime_t tc      I  camera transform matrix time
 *          insstate_t *ins IO ins states
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int zvucam(const insopt_t *opt,const double *dT,int status,const imud_t *imu,
                  gtime_t tc,insstate_t *ins)
{
    double *x,*H,*R,*v,I[9]={-1,0,0,0,-1,0,0,0,-1};
    int nx=ins->nx,info=0;

    trace(3,"zvucam: time=%s",time_str(tc,3));

    if (!chkstatic(dT,status,tc)) {
        trace(2,"no static motion detected\n");
        return 0;
    }
    x=zeros(1,nx); H=zeros(3,nx);
    R=zeros(3,3); v=zeros(3,1);

    /* sensitive matrix */
    asi_blk_mat(H,3,nx,I,3,3,0,3);

    /* variance matrix */
    R[0]=R[4]=R[8]=VARVEL;

    v[0]=ins->ve[0];
    v[1]=ins->ve[1];
    v[2]=ins->ve[2]; /* residual vector */

    if (norm(v,3)<MAXVEL&&norm(imu->gyro,3)<MAXGYRO) {

        /* ekf filter */
        info=filter(x,ins->P,H,v,R,nx,3);

        /* solution fail */
        if (info) {
            trace(2,"zero velocity update filter error\n");
            info=0;
        }
        else {
            /* solution ok */
            ins->stat=INSS_ZVU;
            info=1;
            clp(ins,opt,x);
            trace(3,"zero velocity update ok\n");
        }
    }
    free(x); free(H);
    free(R); free(v);
    return info;
}
/* add ins states-------------------------------------------------------------*/
static int addinstate(const insstate_t *ins,filt_t *filt)
{
    insstate_t inss,*p;

    /* copy ins states to buffer data */
    inss=*ins;

    inss.P=zeros(ins->nx,ins->nx);
    matcpy(inss.P,ins->P,ins->nx,ins->nx);

    /* expend buffer data */
    if (filt->nimax<=filt->ni) {
        if (filt->nimax<=0) filt->nimax=16*2; else filt->nimax*=2;
        if (!(p=(insstate_t*)realloc(filt->insdata,sizeof(insstate_t)*filt->nimax))) {

            /* add ins fail */
            filt->ni=filt->nimax=0;
            free(filt->insdata);
            filt->insdata=NULL;
            return -1;
        }
        filt->insdata=p;
    }
    /* add ins states */
    filt->insdata[filt->ni++]=inss;
    return 1;
}
/* add mono visual odometry data----------------------------------------------*/
static int addvodata(const vostate_t *vo,gtime_t ts,gtime_t te,filt_t *filt)
{
    voss_t voss,*p;
    matcpy(voss.dT,vo->dT,4,4);
    voss.ts=ts;
    voss.te=te;

    /* expend buffer data */
    if (filt->nvmax<=filt->nv) {
        if (filt->nvmax<=0) filt->nvmax=16*2; else filt->nvmax*=2;
        if (!(p=(voss_t*)realloc(filt->vodata,sizeof(voss_t)*filt->nvmax))) {

            /* add ins fail */
            filt->nv=filt->nvmax=0;
            free(filt->vodata);
            filt->vodata=NULL;
            return -1;
        }
        filt->vodata=p;
    }
    /* add ins states */
    filt->vodata[filt->nv++]=voss;
    return 1;
}
/* initial vo-aid-------------------------------------------------------------*/
extern void initvoaidlc(insopt_t *opt)
{
    trace(3,"initvoaidlc:\n");

    opt->voopt.match.f =opt->voopt.calib.f;
    opt->voopt.match.fu=opt->voopt.calib.fu;
    opt->voopt.match.fv=opt->voopt.calib.fv;

    opt->voopt.match.cu=opt->voopt.calib.cu;
    opt->voopt.match.cv=opt->voopt.calib.cv;

    init_match(&matchs,&opt->voopt.match);
    inittrackimgbuf(&opt->voopt);
    return;
}
/* free vo-aid----------------------------------------------------------------*/
extern void freevoaidlc()
{
    int i;
    trace(3,"freevoaidlc:\n");

    freetrackset(&tracks);
    free_match(&matchs);
    freetrackimgbuf();
}
/* initial filter workspace---------------------------------------------------*/
static void initfiltws(insstate_t *ins)
{
    filts.flag=1;
    filts.ni=0;
    filts.nv=0;
    filts.np=0;
    ins2camera(ins,ins->vo.rc ,ins->vo.Cce );
    ins2camera(ins,ins->vo.rc0,ins->vo.Cce0);
    matcpy(ins->vo.Cbc,ins->Cbc,3,3);
    matcpy(ins->vo.lbc,ins->lbc,3,1);
    ins->vo.time=ins->time;

    addinstate(ins,&filts);
}
/* get camera calibration matrix----------------------------------------------*/
static void getcameraK(const insstate_t *ins,double *K)
{
    /* camera calibration matrix */
    K[0]=ins->fx; K[4]=ins->fy;
    K[6]=ins->ox; K[7]=ins->oy;
    K[8]=1.0;
}
/* get camera calibration parameters-------------------------------------------*/
static void getcamcalibp(const insstate_t *ins,cam_t *cam)
{
    cam->k1=ins->k1;
    cam->k2=ins->k2;
    cam->p1=ins->p1;
    cam->p2=ins->p2;
}
/* undistort match feature points----------------------------------------------*/
static void undisfeatpoint(const double *uvd,const double *K,const cam_t *cam,
                           double *uvu)
{
    double uv[3],pf[3],pfn[3];
    double Ki[9];

    matcpy(Ki,K,3,3); if (matinv(Ki,3)) return;
    uv[0]=uvd[0];
    uv[1]=uvd[1]; uv[2]=1.0;

    matmul3v("N",Ki,uv,pf);
    undistortradtan(cam,pf,pfn,NULL); pfn[2]=1.0;
    matmul3v("N",K,pfn,uvu);
}
static void undisdortmatch(match_point *mp,const cam_t *cam,const double *K)
{
    double uv[3],uvu[3];;

    /* current feature point */
    uv[0]=mp->uc; uv[1]=mp->vc;
    undisfeatpoint(uv,K,cam,uvu);
    mp->uc=(float)uvu[0];
    mp->vc=(float)uvu[1];

    /* precious feature point */
    uv[0]=mp->up;
    uv[1]=mp->vp;
    undisfeatpoint(uv,K,cam,uvu);
    mp->up=(float)uvu[0];
    mp->vp=(float)uvu[1];
}
/* undistort feature point-----------------------------------------------------*/
static void undistortfeats(match_set *mset,const voopt_t *opt,const insstate_t *ins)
{
    cam_t cam; double K[9]={0};
    int i;

    /* camera calibration parameters */
    getcamcalibp(ins,&cam);
    getcameraK(ins,K);

    /* undistort feature point */
    for (i=0;i<mset->n;i++) {
        undisdortmatch(&mset->data[i],&cam,K);
    }
}
/* update all track data------------------------------------------------------*/
static int updatetrack(const insopt_t *opt,const img_t *img,const insstate_t *ins)
{
    match_set *pset=&matchs.mp_dense;

#if USE_BUCKET_FEATS
    pset=&matchs.mp_bucket;
#endif
    /* update track data */
    if (!match2track(pset,matchs.pt,matchs.time,img->id,&matchs.Ip,&matchs.Ic,
                     &opt->voopt,&tracks)) {

        trace(2,"update track fail\n");
        return 0;
    }
#if GMS
    /* get matches from GMS-matcher */
    getgmsmatches(&opt->voopt,&matchs.mp_bucket,&matchs.Ip,&matchs.Ic);
#endif
#if 0
    drawalltrack(&tracks);
#endif
    /* undistort feature point */
    undistortfeats(&matchs.mp_bucket,&opt->voopt,ins);
    return 1;
}
/* estimate mono-camera motion-------------------------------------------------*/
static int estmotion(const match_set *pset,const voopt_t *opt,const insstate_t *ins,
                     double *dT,double *ratio)
{
#if OPENCV
    return solveRt(opt,ins,pset,dT,ratio);
#else
    return estmonort(opt,pset,dT,ratio);
#endif
}
static int estmonocamera(const voopt_t *opt,const img_t *img,double *dT,double *ratio)
{
    return estvo(opt,img,dT,ratio);
}
static void camdp2insdp(const double *dpc,const double *Cbc,const double *lbc,const double *dC,
                        double *dpi)
{
    double dp1[3],dp2[3],dp[3];
    int i;

    matmul("NN",3,1,3,-1.0,Cbc,lbc,0.0,dp1);
    matmul33("NNN",dC,Cbc,lbc,3,3,3,1,dp2);
    for (i=0;i<3;i++) dp[i]=dpc[i]-dp1[i]-dp2[i];
    matmul("TN",3,1,3,1.0,Cbc,dp,0.0,dpi);
}
/* update track data-----------------------------------------------------------*/
static int updatecameratrack(insstate_t *ins,const insopt_t *opt,const img_t *img,
                             vostate_t *vo)
{
    /* match feature points */
    if (img==NULL&&img->feat==NULL) {
#if !LIBVISO
        trace(2,"no feature points measurement data\n");
        return 0;
#endif
    }
    if (matchfeats(&matchs,img)<=0) {
        trace(2,"match feature points fail\n");
        return 0;
    }
    /* update track data */
    if (!updatetrack(opt,img,ins)) return 0;
#if LIBVISO
    /* using libviso to solve R|t */
    vo->status=estmonocamera(&opt->voopt,img,vo->dT,&vo->ratio);
#else
    /* estimate mono-camera motion */
    vo->status=estmotion(&matchs.mp_bucket,&opt->voopt,ins,vo->dT,&vo->ratio);
#endif
    if (!vo->status||vo->ratio<THRES_RATIO) {
        trace(2,"estimate motion fail\n");

        vo->status=0;
        return 0;
    }
    /* update camera pose */
    updatecamera(ins,opt,vo->dT,img->time);
    return 1;
}
/* jacobians of bg------------------------------------------------------------*/
static void jacob_bgk(const double *Rk_1,const double *Rj,const double *omgb,
                      const double dt,double *Jbgk)
{
    double dR[9],phi[3],Jr[9];

    matmul("TN",3,3,3,1.0,Rk_1,Rj,0.0,dR);
    phi[0]=omgb[0]*dt;
    phi[1]=omgb[1]*dt;
    phi[2]=omgb[2]*dt;
#if 1
    so3jac(phi,Jr);
#else
    so3_jac(phi,Jr,NULL);
#endif
    matmul("TN",3,3,3,dt,dR,Jr,0.0,Jbgk);
}
static void jacob_bg(const double *z,const double dt,double *Jbg)
{
    insstate_t *pins=filts.insdata;
    double Jbgk[9],Jbgi[9]={0},Jr[9],*Rk_1,*Rj,Rz[9];
    int i,j;

    for (Rj=pins[filts.ni-1].Cbe,i=1;i<filts.ni;i++) {
        Rk_1=pins[i].Cbe;
        jacob_bgk(Rk_1,Rj,pins[i-1].omgb,dt,Jbgk);
        for (j=0;j<9;j++) {
            Jbgi[j]+=Jbgk[j];
        }
    }
#if 1
    so3jac(z,Jr);
#else
    so3_jac(z,Jr,NULL);
#endif
    so3_exp(z,Rz);
    matmul33("NTN",Jr,Rz,Jbgi,3,3,3,3,Jbg);
}
/* check estimated states----------------------------------------------------*/
static int chkest_state(const double *dx,const double *P,const insopt_t *opt,const insstate_t *ins)
{
    static int iba=xiBa(opt),nba=xnBa(opt);
    static int ibg=xiBg(opt),nbg=xnBg(opt);
    static double factor=3.0;
    static double pba[3]={0},pbg[3]={0};
    int i,nx=xnX(opt),flag=0;

    /* check estimated states */
    if (     dx[  0]==DISFLAG&&norm(dx+  0,3)>15.0*D2R) flag|=1;
    if (nba&&dx[iba]==DISFLAG&&norm(dx+iba,3)>1E5*Mg2M) flag|=1;
    if (nbg&&dx[ibg]==DISFLAG&&norm(dx+ibg,3)>15.0*D2R) flag|=1;
    for (i=0;i<3;i++) {
        if (SQR(dx[i+xiBg(opt)])<=SQR(factor)*P[i+xiBg(opt)+(i+xiBg(opt))*nx]) continue;
        flag|=1;
    }
    for (i=0;i<3;i++) {
        if (SQR(dx[i+xiBa(opt)])<=SQR(factor)*P[i+xiBa(opt)+(i+xiBa(opt))*nx]) continue;
        flag|=1;
    }
#if 0
    /* check accl bias and gyro bias */
    if (norm(pba,3)>0.0) {
        for (i=0;i<3;i++) if (fabs(dx[iba+i])/fabs(pba[i])>0.3) flag|=1;
    }
    if (norm(pbg,3)>0.0) {
        for (i=0;i<3;i++) if (fabs(dx[ibg+i])/fabs(pbg[i])>0.3) flag|=1;
    }
    matcpy(pba,ins->ba,1,3);
    matcpy(pbg,ins->bg,1,3);
#endif
    if (flag) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    return 1;
}
/* trace filter workspace-----------------------------------------------------*/
static void trace_fltworkspace(const insstate_t *ins,const insopt_t *opt)
{
    double dp[3],da[3],dC[9];
    int i;

    trace(3,"imu data : n=%d\n",filts.ni);
    for (i=0;i<filts.ni;i++) {

        fprintf(stderr,"%6.3lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf\n",
                time2gpst(filts.insdata[i].time,NULL),
                filts.insdata[i].fb[0],
                filts.insdata[i].fb[1],
                filts.insdata[i].fb[2],
                filts.insdata[i].omgb[0],
                filts.insdata[i].omgb[1],
                filts.insdata[i].omgb[2]);
    }
    trace(3,"vo data : n=%d\n",filts.nv);
    for (i=0;i<filts.nv;i++) {

        tf2rt(filts.vodata[i].dT,dC,dp);
        so3_log(dC,da,NULL);
        fprintf(stderr,"%.3lf-%.3lf: %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf  %8.4lf\n",
                time2gpst(filts.vodata[i].ts,NULL),
                time2gpst(filts.vodata[i].te,NULL),
                da[0],da[1],da[2],
                dp[0],dp[1],dp[2]);
    }
}
/* get current max track length-----------------------------------------------*/
static int gettrackstatus(const track_t *track)
{
    double ratio;
    int i,j,n,nlost=0;

    for (j=0,n=0,i=0;i<track->n;i++) {
        if (track->data[i].flag!=TRACK_UPDATED) continue;
        if (track->data[i].flag==TRACK_LOST) nlost++;
        if (track->data[i].n>=MIN_LEN) j++;
        n++;
    }
    trace(3,"nlost=%d  j=%d  n=%d\n",nlost,j,n);

    ratio=double(j)/double(n);
    return ratio>THRES_TRACK_RATIO;
}
/* construct measurement data-------------------------------------------------*/
static int bldmeasuredata(gtime_t ts,gtime_t te,int *start,int *end,int *index,int n,double *dCi,double *dCc,double *dpc)
{
    double T0[16],dT[16],dTi[16],*pCe,*cCe;
    int i;

    for (*start=0,i=0;i<filts.ni;i++) {
        if (fabs(timediff(filts.insdata[i].time,ts))<=DTTOL) {
            *start=i;
            break;
        }
    }
    for (*end=0,i=0;i<filts.ni;i++) {
        if (fabs(timediff(filts.insdata[i].time,te))<=DTTOL) {
            *end=i;
            break;
        }
    }
    pCe=filts.insdata[*start].Cbe;
    cCe=filts.insdata[  *end].Cbe;
    matmul("TN",3,3,3,1.0,pCe,cCe,0.0,dCi);
    for (seteye(T0,4),i=0;i<n;i++) {
        matcpy(dTi,filts.vodata[index[i]].dT,4,4);
        matinv(dTi,4);

        matmul("NN",4,4,4,1.0,T0,dTi,0.0,dT);
        matcpy(T0,dT,4,4);
    }
    tf2rt(dT,dCc,dpc);
    return *end>*start;
}
/* vo aid filter--------------------------------------------------------------*/
static int voflt_att(insstate_t *ins,const insopt_t *opt,const vostate_t *vo,
                     int start,int end,
                     const double *dCc,const double *dCi)
{
    double dCz[9],dC[9],r[3],z[3],Jbg[9];
    double phiz[3],phi[3],phib[3];
    double *H,*v,*R,*x,*P;
    int i,j,nv,nx=ins->nx;
    int ibg,nbg;

    trace(3,"voflt_att:\n");

    matmul33("TNN",ins->Cbc,dCc,ins->Cbc,3,3,3,3,dCz);

    so3_log(dCi,phi ,NULL);
    so3_log(dCz,phiz,NULL);

    trace(3,"phi =\n"); tracemat(3,phi ,1,3,12,5);
    trace(3,"phiz=\n"); tracemat(3,phiz,1,3,12,5);
#if OPENCV
    matcpy(phib,phiz,1,3);
#else
    matmul("TN",3,1,3,1.0,ins->Cbc,phiz,0.0,phib);
    trace(3,"phib=\n"); tracemat(3,phib,1,3,12,5);
#endif
    so3_exp(phib,dCz);
    matmul("TN",3,3,3,1.0,dCi,dCz,0.0,dC);
    so3_log(dC,z,NULL);

    H=zeros(3,nx); v=zeros(3, 1);
    R=zeros(3, 3); x=zeros(nx,1);
    P=zeros(nx,nx);

    ibg=xiBg(opt);
    nbg=xnBg(opt);

    jacob_bg(z,ins->dt,Jbg);

    for (nv=0,i=0;i<3;i++) {
        if (fabs(v[nv]=z[i])>MAXRES_POSE) continue;

        for (j=ibg;j<ibg+nbg;j++) H[j+nv*nx]=-Jbg[i+(j-ibg)*3];
        r[nv++]=SQR(VARPOSE);
    }
    if (v&&nv) {
        trace(3,"v=\n"); tracemat(3,v,3,1,12,6);
    }
    if (H&&nv) {
        trace(3,"H=\n"); tracemat(3,H,nx,nv,12,6);
    }
    if (R&&nv) {
        for (i=0;i<nv;i++) {
            R[i+i*nv]=r[i];
        }
        trace(3,"R=\n"); tracemat(3,R,nv,nv,12,6);
    }
    if (nv<=0) {
        trace(2,"pose fusion filter fail\n");
        free(x); free(v); free(P);
        free(H); free(R);
        return 0;
    }
    matcpy(P,ins->P,nx,nx);
    if (!filter(x,P,H,v,R,nx,nv)) {
        if (!chkest_state(x,P,opt,ins)) goto exit;

        /* corrections */
        clp(ins,opt,x);

        /* update covariance matrix */
        matcpy(ins->P,P,nx,nx);

        trace(3,"dx=\n");
        tracemat(3,x,nx,1,12,6);

        trace(3,"Px=\n");
        tracemat(3,P,nx,nx,12,6);
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(v); free(P);
    free(H); free(R);
    return 1;
}
/* jacobians of velocity wrt ba-----------------------------------------------*/
static void jacob_dvdba(const double dt,int i,int j,double *Jdvdba)
{
    double *Rk,*Ri,dR[9];
    int k,p;

    setzero(Jdvdba,3,3);

    for (Ri=filts.insdata[i].Cbe,k=i;k<j;k++) {
        Rk=filts.insdata[k].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);

        for (p=0;p<9;p++) {
            Jdvdba[p]+=dR[p]*dt;
        }
    }
}
/* jacobians of position wrt ba-----------------------------------------------*/
static void jacob_dpdba(const double dt,int start,int end,double *Jdpdba)
{
    double Jdvdba[9],dR[9],*Rk,*Ri;
    int i,j;

    setzero(Jdpdba,3,3);

    for (Ri=filts.insdata[0].Cbe,i=start;i<=end;i++) {
        Rk=filts.insdata[i].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);

        jacob_dvdba(dt,start,i,Jdvdba);
        for (j=0;j<9;j++) {
            Jdpdba[j]+=Jdvdba[j]*dt-0.5*dR[j]*dt*dt;
        }
    }
}
/* jacobians of attitude wrt bg-----------------------------------------------*/
static void jacob_dadbg(int p,int q,double *Jbg)
{
    insstate_t *pins=filts.insdata;
    double Jbgk[9],*Rk_1,*Rj;
    int i,j;

    setzero(Jbg,3,3);

    for (Rj=pins[q].Cbe,i=p+1;i<q;i++) {
        Rk_1=pins[i].Cbe;
        jacob_bgk(Rk_1,Rj,pins[i-1].omgb,pins[i-1].dt,Jbgk);
        for (j=0;j<9;j++) Jbg[j]+=Jbgk[j];
    }
}
/* jacobians of velocity wrt bg-----------------------------------------------*/
static void jacob_dvdbg(const double dt,int i,int j,double *Jdvdbg)
{
    double *Rk,*Ri,dR[9],W[9],Jdadbg[9];
    double T[9];
    int k,p;

    setzero(Jdvdbg,3,3);

    for (Ri=filts.insdata[i].Cbe,k=i;k<j;k++) {
        Rk=filts.insdata[k].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);
        skewsym3(filts.insdata[k].fb,W);

        jacob_dadbg(i,k,Jdadbg);

        matmul33("NNN",dR,W,Jdadbg,3,3,3,3,T);
        for (p=0;p<9;p++) {
            Jdvdbg[p]+=T[p]*dt;
        }
    }
}
/* jacobians of position wrt bg-----------------------------------------------*/
static void jacob_dpdbg(const double dt,int start,int end,double *Jdpdbg)
{
    double Jdvdbg[9],Jdadbg[9],dR[9],*Rk,*Ri;
    double W[9],T[9];
    int k,j;

    setzero(Jdpdbg,3,3);

    for (Ri=filts.insdata[start].Cbe,k=start;k<=end;k++) {
        Rk=filts.insdata[k].Cbe;
        matmul("TN",3,3,3,1.0,Ri,Rk,0.0,dR);
        jacob_dvdbg(dt,start,k,Jdvdbg);

        jacob_dadbg(start,k,Jdadbg);

        skewsym3(filts.insdata[k].fb,W);
        matmul33("NNN",dR,W,Jdadbg,3,3,3,3,T);
        for (j=0;j<9;j++) {
            Jdpdbg[j]+=Jdvdbg[j]*dt-0.5*T[j]*dt*dt;
        }
    }
}
/* delta velocity due to accl measurement--------------------------------------*/
static void deltavel(int i,int j,double *dv)
{
    double dRk[9],dvk[3];
    int k;

    setzero(dv,1,3);
    for (k=i;k<j;k++) {
        matmul("TN",3,3,3,1.0,filts.insdata[i].Cbe,filts.insdata[k].Cbe,0.0,dRk);
        matmul("NN",3,1,3,filts.insdata[k].dt,dRk,filts.insdata[k].fb,0.0,dvk);
        dv[0]+=dvk[0];
        dv[1]+=dvk[1];
        dv[2]+=dvk[2];
    }
}
static void deltapos(int start,int end,double *dp)
{
    double dRk[9],dvk[3],dpk[3];
    int k;

    setzero(dp,1,3);
    for (k=start;k<end;k++) {
        matmul("TN",3,3,3,1.0,filts.insdata[start].Cbe,filts.insdata[k].Cbe,0.0,dRk);
        matmul("NN",3,1,3,0.5*filts.insdata[k].dt*filts.insdata[k].dt,dRk,filts.insdata[k].fb,0.0,dpk);
        deltavel(start,k,dvk);
        dp[0]+=dpk[0]+dvk[0]*filts.insdata[k].dt;
        dp[1]+=dpk[1]+dvk[1]*filts.insdata[k].dt;
        dp[2]+=dpk[2]+dvk[2]*filts.insdata[k].dt;
    }
}
/* filter for delta position--------------------------------------------------*/
static int voflt_pos(insstate_t *ins,const insopt_t *opt,const vostate_t *vo,
                     int start,int end,
                     const double *dpc,const double *dC)
{
    double dpi[3],ge[3],dpk[3],dpz[3],Jdpdba[9],Jdpdbg[9],Jdpda[9];
    double dgr[3],Wg[9],dt;
    double dvp[3],dr[3],T[3];
    double *v,*H,*R,*x,*P,r[3];
    int ibg,nbg,iba,nba,iax,nax,ivx,nvx,nv=0;
    int i,j,nx=ins->nx;

    trace(3,"voflt_pos:\n");

    camdp2insdp(dpc,ins->Cbc,ins->lbc,dC,dpi);

    trace(3,"dpi=\n");
    tracemat(3,dpi,1,3,12,5);

    for (dt=0.0,i=start;i<end;i++) dt+=filts.insdata[i].dt;
    deltapos(start,end,dpk);
    pregrav(filts.insdata[start].re,ge);

    T[0]=filts.insdata[end].re[0]-filts.insdata[start].re[0];
    T[1]=filts.insdata[end].re[1]-filts.insdata[start].re[1];
    T[2]=filts.insdata[end].re[2]-filts.insdata[start].re[2];
    matmul("TN",3,1,3,1.0,filts.insdata[start].Cbe,T,0.0,dr);

    trace(3,"dr=\n");
    tracemat(3,dr,1,3,12,5);

    matmul("TN",3,1,3,dt,filts.insdata[start].Cbe,filts.insdata[start].ve,0.0,dvp);
    matmul("TN",3,1,3,0.5*dt*dt,filts.insdata[start].Cbe,ge,0.0,dgr);
    for (i=0;i<3;i++) {
        dpz[i]=dpk[i]+dgr[i]+dvp[i];
    }
    trace(3,"dpz=\n");
    tracemat(3,dpz,1,3,12,5);

    jacob_dpdba(ins->dt,start,end,Jdpdba);
    jacob_dpdbg(ins->dt,start,end,Jdpdbg);

    skewsym3(dgr,Wg);
    matmul("TN",3,3,3,1.0,filts.insdata[start].Cbe,Wg,0.0,Jdpda);

    ibg=xiBg(opt); nbg=xnBg(opt);
    iba=xiBa(opt); nba=xnBa(opt);
    iax= xiA(opt); nax= xnA(opt);
    ivx= xiV(opt); nvx= xnV(opt);

    H=zeros(3,nx); v=zeros(3, 1);
    R=zeros(3, 3); x=zeros(nx,1);
    P=zeros(nx,nx);

    for (i=0;i<3;i++) {
        if (fabs(v[nv]=dpz[i]-dpi[i])>10.0) continue;
        for (j=iba;j<iba+nba;j++) H[j+nv*nx]=Jdpdba[i+(j-iba)*3];
        for (j=ibg;j<ibg+nbg;j++) H[j+nv*nx]=Jdpdbg[i+(j-ibg)*3];

        for (j=ivx;j<ivx+nvx;j++) H[j+nv*nx]=filts.insdata[start].Cbe[j-ivx+i*3]*dt;
        for (j=iax;j<iax+nax;j++) H[j+nv*nx]=Jdpda[i+(j-iax)*3];
        r[nv++]=SQR(VARPOS);
    }
    if (v&&nv) {
        trace(3,"v=\n"); tracemat(3,v,3,1,12,6);
    }
    if (H&&nv) {
        trace(3,"H=\n");
        tracemat(3,H,nx,nv,12,6);
    }
    if (R&&nv) {
        for (i=0;i<nv;i++) R[i+i*nv]=r[i];
        trace(3,"R=\n");
        tracemat(3,R,nv,nv,12,6);
    }
    if (nv<=0) {
        trace(2,"position fusion filter fail\n");
        free(x); free(v); free(P);
        free(H); free(R);
        return 0;
    }
    matcpy(P,ins->P,nx,nx);
    if (!filter(x,P,H,v,R,nx,nv)) {
        if (!chkest_state(x,P,opt,ins)) goto exit;
#if 1
        /* partial feedback */
        for (i=0;i<nbg;i++) {
            x[ibg+i]*=1E-3;
        }
        for (i=0;i<nax;i++) {
            x[iax+i]*=1E-3;
        }
#endif
        /* corrections */
        clp(ins,opt,x);
#if 1
        /* update covariance matrix */
        matcpy(ins->P,P,nx,nx);
#endif
        trace(3,"dx=\n");
        tracemat(3,x,nx,1,12,6);

        trace(3,"Px=\n");
        tracemat(3,P,nx,nx,12,6);
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(v); free(P);
    free(H); free(R);
    return 1;
}
/* filter update--------------------------------------------------------------*/
static int voflt(insstate_t *ins,const insopt_t *opt,const vostate_t *vo)
{
    int i,j=0,index[MAX_LEN][MAX_LEN]={0},n[MAX_LEN]={0},info=0,start,end;
    double dCc[9],dCi[9];
    double dpc[3];
    gtime_t tt[MAX_LEN][2]={0},ts,te;

    trace(3,"voflt:\n");

    trace_fltworkspace(ins,opt);

    ts=filts.vodata[0].ts;
    te=filts.vodata[0].te;
    tt[j][0]=ts;
    tt[j][1]=te;
    index[j][n[j]++]=0;

    for (i=1;i<filts.nv;i++) {
        if (fabs(timediff(filts.vodata[i].ts,te))<1E-5) {
            te=filts.vodata[i].te;
            tt[j][1]=te;
            index[j][n[j]++]=i;
            continue;
        }
        else {
            j++;
            ts=filts.vodata[i].ts;
            te=filts.vodata[i].te;
            tt[j][0]=ts;
            tt[j][1]=te;
            index[j][n[j]++]=i;
        }
    }
    j++;
    for (i=0;i<j;i++) {
        if (n[i]<2) continue;
        if (!bldmeasuredata(tt[i][0],tt[i][1],&start,&end,index[i],n[i],dCi,dCc,dpc)) continue;
        info|=voflt_att(ins,opt,vo,start,end,dCc,dCi);
#if 0
        info|=voflt_pos(ins,opt,vo,start,end,dpc,dCc);
#endif
    }
    if (info) ins->stat=INSS_VO;
    return info;
}
/* update.--------------------------------------------------------------------*/
static int updateall(insstate_t *ins,const insopt_t *opt,const img_t *img)
{
    int info=0,trackstat,reset=0;

    trace(3,"updateall:\n");

    /* initial filter workspace */
    if (filts.flag==0) initfiltws(ins);
    
    info=updatecameratrack(ins,opt,img,&ins->vo);
    if (info) {
        /* add vo data */
        addvodata(&ins->vo,matchs.pt,matchs.time,&filts);
#if 1
        /* vo filter */
        trackstat=gettrackstatus(&tracks);

        if (trackstat&&filts.nv>=MIN_LEN) {
            info=voflt(ins,opt,&ins->vo);
            reset=1;
        }
#endif
    }
    else reset=1;
    if (reset) {

        /* reset filter */
        filts.flag=0;
    }
    /* re-initial filter */
    if (filts.flag==0) initfiltws(ins);
    return info;
}
/* propagate covariance matrix of ins-camera pose in sliding windows----------*/
static int propagate(insstate_t *ins,const insopt_t *opt)
{
    /* add a ins states */
    if (!addinstate(ins,&filts)) return 0;
    return 1;
}
/* using visual odometry to aid ins/gnss pose estimating-----------------------
 * args:    insopt_t *opt   I   ins options
 *          insstate_t *ins IO  ins states
 *          imud_t *imu     I   imu measurement data
 *          img_t *img      I   image measurement data
 *          int flag        I   update flag
 * return: status (1: ok, 0: fail)
 * ----------------------------------------------------------------------------*/
extern int voigposlc(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                     const img_t *img,int flag)
{
    trace(3,"voigposlc: time=%s\n",time_str(imu->time,4));

    switch (flag) {
        case 0: return propagate(ins,opt);
        case 1: return propagate(ins,opt)&&updateall(ins,opt,img);
        default: {
            trace(2,"not support mode\n");
        }
    }
    return 0;
}

