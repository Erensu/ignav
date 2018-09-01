/*------------------------------------------------------------------------------
* ins-vo-aid.cc : visual odometry aid ins-gnss coupled functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*    [5] Weiss,Real-Time Metric State Estimation for Modular Vision-Inertial
*        Systems
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/03/15 1.0 new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants -----------------------------------------------------------------*/
#define NUMPOSE              10   /* numbers of camera transform matrix for detect motion*/
#define THRES_ANG     (3.0*D2R)   /* thresholds of angular rate for detect motion */
#define THRES_TRANS       (0.1)   /* thresholds of translation for detect motion */
#define MAXVEL              0.1   /* max velocity for using non-holonomic constraint */
#define MINVEL              3.0   /* min velocity for initial camera pose */
#define MAXGYRO      (10.0*D2R)   /* max rotation speed for using non-holonomic constraint */
#define VARVEL        SQR(0.05)   /* initial variance of receiver vel ((m/s)^2) */

/* type definitions-----------------------------------------------------------*/
typedef struct {
    gtime_t gtime;                /* pose measurement time */
    int status;                   /* visual odometry estimate status */
    double rpy[3];                /* roll/pitch/yaw of precious to current frame */
    double tra[3];                /* translation of precious to current frame */
} pose_t;

/* transform matrix convert to pose-------------------------------------------*/
static int convpose(const double *T,pose_t *pose)
{
    double R[9]; tf2rt(T,R,pose->tra); dcm2rpy(R,pose->rpy);
}
/* detect motion by using camera motions--------------------------------------*/
static int chk_static(const double *dT,int status,gtime_t time)
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
        if (norm(T[i].rpy,3)>THRES_ANG||norm(T[i].tra,3)>THRES_TRANS) {
            return 0;
        }
    }
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
extern int zvu_cam(const insopt_t *opt,const double *dT,int status,
                   const imud_t *imu,gtime_t tc,insstate_t *ins)
{
    double *x,*H,*R,*v,I[9]={-1,0,0,0,-1,0,0,0,-1};
    int nx=ins->nx,info=0;

    trace(3,"zvu_cam: time=%s",time_str(tc,3));

    if (!chk_static(dT,status,tc)) {
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
/* initial camera attitude (relative to n-frame) based on ins states----------
 * args  :  insstate_t *ins  IO  ins states
 *          insopt_t *opt    I   ins options
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int initcampose(insstate_t *ins,const insopt_t *opt)
{
    double T[3],W[3],w[9];
    int i;

    trace(3,"initcamatt:\n");

    seteye(ins->vo.T,4);
    matmul("NN",3,3,3,1.0,ins->Cbe,ins->vo.Ccb,0.0,ins->vo.Cce);
#if 1
    /* expressed in tangent space */
    so3_log(ins->vo.Cce,ins->vo.phi,NULL);
#endif

    matmul("NN",3,1,3,1.0,ins->Cbe,ins->vo.lever,0.0,T);
    for (i=0;i<3;i++) {
        ins->vo.rc[i]=ins->re[i]+T[i];
    }
#if 1
    matcpy(ins->vo.vc,ins->ve,3,1);
#else
    matmul33("NNN",Omge,ins->Cbe,ins->vo.lever,3,3,3,1,T);
    skewsym3(ins->omgb,w);
    matmul33("NNN",ins->Cbe,w,ins->vo.lever,3,3,3,1,W);
    for (i=0;i<3;i++) {
        ins->vo.vc[i]=ins->ve[i]+W[i]+T[i];
    }
#endif
    ins->vo.time=ins->time;
    return 1;
}
/* initial visual odometry options-------------------------------------------*/
extern void initvo(insstate_t *ins,const insopt_t *opt)
{
    matcpy(ins->vo.lever,opt->voopt.lever,1,3);
    rpy2dcm(opt->voopt.rpy,ins->vo.Ccb);
    seteye(ins->vo.T,4);
    setzero(ins->vo.phi,1,3);
}
/* initial camera pose using gps solutions-----------------------------------
 * args   :  gtime_t t       I  time of gps solution data
 *           double *r       I  gps position solution
 *           insopt_t *opt   I  ins options
 *           insstate_t *ins I  ins states
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int initcamgps(gtime_t time,const double *r,const insopt_t *opt,
                      insstate_t *ins)
{
    static gtime_t ts[2]={0};
    static double rr[3*2]={0};
    double ve[3]={0},vn[3]={0},pos[3],C[9],rpy[3]={0};
    double Cnb[9],Cbe[9],T[3];
    double W[3],w[9];
    int i;

    if (!time.time||norm(r,3)<=0.0) return 0;

    matcpy(&rr[0],&rr[3],1,3);
    matcpy(&rr[3],r,1,3);

    ts[0]=ts[1]; ts[1]=time;

    if (fabs(timediff(ts[1],ts[0]))<1E-5) return 0;
    for (i=0;i<2;i++) {
        if (norm(rr+3*i,3)<=0.0) return 0;
    }
    for (i=0;i<3;i++) {
        ve[i]=(rr[3+i]-rr[i])/timediff(ts[1],ts[0]);
    }
    if (norm(ve,3)<MINVEL) return 0;
    ecef2pos(r,pos);
    ned2xyz(pos,C);

    /* velocity in navigation frame */
    matmul("TN",3,1,3,1.0,C,ve,0.0,vn);

    rpy[2]=vel2head(vn); /* yaw */

    /* initial camera attitude */
    rpy2dcm(rpy,Cnb);
    matmul("NT",3,3,3,1.0,C,Cnb,0.0,Cbe);

    matmul("NN",3,3,3,1.0,Cbe,ins->vo.Ccb,0.0,ins->vo.Cce);
#if 1
    /* expressed in tangent space */
    so3_log(ins->vo.Cce,ins->vo.phi,NULL);
#endif

    /* initial camera position */
    matmul("NN",3,1,3,1.0,Cbe,ins->vo.lever,0.0,T);
    for (i=0;i<3;i++) {
        ins->vo.rc[i]=r[i]+T[i];
    }
#if 1
    /* initial camera velocity */
    matcpy(ins->vo.vc,ins->ve,3,1);
#else
    /* initial camera velocity */
    matmul33("NNN",Omge,ins->Cbe,ins->vo.lever,3,3,3,1,T);
    skewsym3(ins->omgb,w);
    matmul33("NNN",ins->Cbe,w,ins->vo.lever,3,3,3,1,W);
    for (i=0;i<3;i++) {
        ins->vo.vc[i]=ins->ve[i]+W[i]+T[i];
    }
#endif
    ins->vo.time=time;

    setzero(rr,1,6); return 1;
}
/* update camera pose by using visual odometry-------------------------------
 * args  :  insstate_t *ins  IO  ins states
 *          insopt_t *opt    I   ins options
 *          double *dT       I   transform of precious to current frame
 *          gtime_t time     I   time of transform measurement
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int updatecam(insstate_t *ins,const insopt_t *opt,
                     const double *dT,gtime_t time)
{
    double Tp[16],T[16],TT[16],rp[3],dt;
    int i;
    trace(3,"updatecam:\n");

    if ((dt=timediff(time,ins->vo.time))>10.0) {
        trace(2,"update camera pose failed\n");
        return 0;
    }
    if (fabs(dt)<1E-5) {
        trace(2,"duplicate update\n");
        return 0;
    }
    matcpy(rp,ins->vo.rc,1,3);
    rt2tf(ins->vo.Cce,ins->vo.rc,Tp);
    matcpy(TT,dT,4,4);
    matinv(TT,4);
    matmul("NN",4,4,4,1.0,Tp,TT,0.0,T);

    tf2rt(T,ins->vo.Cce,ins->vo.rc);
#if 1
    /* expressed in tangent space */
    so3_log(ins->vo.Cce,ins->vo.phi,NULL);
#endif
    for (i=0;i<3;i++) {
        if (dt>0.0) ins->vo.vc[i]=(ins->vo.rc[i]-rp[i])/dt;
        else {
            ins->vo.vc[i]=0.0;
        }
    }
    matmul("NN",4,4,4,1.0,ins->vo.T,TT,0.0,T);
    matcpy(ins->vo.T,T,4,4);

    ins->vo.time=time;
    return 1;
}




