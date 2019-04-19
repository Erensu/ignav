/*------------------------------------------------------------------------------
* ins-odo.cc : odometry aid for ins navigation functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/11/13 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define VARVEL      SQR(0.1)     /* variance of odometry velocity measurement */
#define MAXINOV     30.0         /* max innovations for updates ins state (m/s) */

/* global states index -------------------------------------------------------*/
static int ibg=0,nbg=0;          /* index and number of gyro bias states */
static int ios=0,nos=0;          /* index and number of scale factor for odometry */
static int iol=0,nol=0;          /* index and number of lever arm for odometry */
static int ioa=0,noa=0;          /* index and number of misalignment between odometry and imu body frame */

/* convert velocity in rear-wheel frame to body frame--------------------------
 * args    :  double *vr     I  velocity in rear-wheel frame
 *            double *lever  I  lever arm from odometry to body
 *            double *Cbe    I  transform matrix of body frame to ecef frame
 *            imud_t *imu    I  imu measurement data
 *            double s       I  odometry scale factor
 *            double *vb     O  velocity in body frame
 * return  : none
 * --------------------------------------------------------------------------*/
extern void vr2vb(const double *vr,const double *lever,const double *Cbe,
                  const imud_t *imu,const double s,double *vb)
{
    double T[3],S[9];

    matcpy(T,imu->gyro,1,3);
    matmul("TN",3,1,3,-1.0,Cbe,Omge,1.0,T);

    skewsym3(T,S);
    matmul("NN",3,1,3,1.0,S,lever,0.0,T);

    vb[0]=(1.0+s)*(vr[0]-T[0]);
    vb[1]=(1.0+s)*(vr[1]-T[1]);
    vb[2]=(1.0+s)*(vr[2]-T[2]);
}
/* convert velocity in ecef frame to rear frame------------------------------
 * args    :  double *ve     I  velocity in ecef frame
 *            double *lever  I  lever arm from odometry to body
 *            double *Cbe    I  transform matrix of body frame to ecef frame
 *            double *Cbr    I  transform matrix of body frame to rear frame
 *            imud_t *imu    I  imu measurement data
 *            double s       I  odometry scale factor
 *            double *vr     O  velocity in rear-wheel frame
 * return  : none
 * --------------------------------------------------------------------------*/
extern int ve2vr(const double *ve,const double *Cbe,const double *Cbr,
                 const double *lever,const imud_t *imu,const double s,
                 double *vr)
{
    double T[3],S[9],omg[3]={0,0,-OMGE},vb[3];

    matcpy(T,imu->gyro,1,3);
    matmul("TN",3,1,3,1.0,Cbe,omg,1.0,T);

    skewsym3(T,S);
    matmul("TN",3,1,3,1.0,Cbe,ve,0.0,T);
    matcpy(vb,T,1,3);
    matmul("NN",3,1,3,1.0,S,lever,1.0,vb);
    matmul("NN",3,1,3,1.0+s,Cbr,vb,0.0,vr);
}
/* initial odometry estimated states-----------------------------------------
 * args   : odopt_t *opt     I  odometry options
 *          insstate_t *ins  IO ins states
 * return : none
 * --------------------------------------------------------------------------*/
extern void initodo(const odopt_t *opt,insstate_t *ins)
{
    double I[9]={1,0,0,0,1,0,0,0,1};

    trace(3,"initodo:\n");

    matcpy(ins->rbl,opt->lever,1,3);
    matcpy(ins->Cbr,I,3,3);

    ins->os=opt->s;
}
/* jacobian of odometry velocity measurement by odometry scale factor--------*/
static void jacobian_ov_ds(const double *vr,double *dvdds)
{
    int i; for (i=0;i<3;i++) dvdds[i]=-vr[i];
}
/* jacobian of odometry velocity measurement by gyro bias--------------------*/
static void jacobian_ov_bg(const double *Cbe,const double *Cbr,const double *lever,
                           double *dvdbg)
{
    double sl[9];
    skewsym3(lever,sl);
    matmul("NN",3,3,3,1.0,Cbr,sl,0.0,dvdbg);
}
/* jacobian of odometry velocity measurement by attitude error---------------*/
static void jacobian_ov_da(const double *ve,const double *Cbe,const double *Cbr,
                           const double *lever,double *dvdda)
{
    double T[9],sl[9];

    skewsym3(lever,T);
    matmul33("NTN",T,Cbe,Omge,3,3,3,3,sl);

    matcpy(dvdda,sl,3,3);
    skewsym3(ve,T);
    matmul("TN",3,3,3,-1.0,Cbe,T,-1.0,sl);
    matmul("NN",3,3,3,1.0,Cbr,sl,0.0,dvdda);
}
/* jacobian of odometry velocity measurement by ins velocity error-----------*/
static void jacobian_ov_dv(const double *Cbe,const double *Cbr,double *dvddv)
{
    int i; double T[9];
    matt(Cbe,3,3,T);
    for (i=0;i<9;i++) T[i]=-T[i];
    matmul("NN",3,3,3,1.0,Cbr,T,0.0,dvddv);
}                                                      
/* jacobian of odometry velocity measurement by misalignment for odometry frame
 * and body frame
 * args    :  double *Cbe    I  transform matrix from body frame to ecef frame
 *            double *Cbr    I  transform matrix form body frame to rear frame
 *            double *vb     I  odometry velocity in body frame
 *            double *dvdma  O  output jacobian matrix
 * return  : none
 * --------------------------------------------------------------------------*/
static void jacobian_ov_dma(const double *vb,const double *Cbe,const double *Cbr,
                            double *dvdma)
{
    double sl[3];
    matmul("NN",3,1,3,1.0,Cbr,vb,0.0,sl);
    skewsym3(sl,dvdma);
}
/* jacobian of odometry velocity measurement by lever arm--------------------*/
static void jacobian_ov_lever(const double *Cbe,const double *Cbr,
                              const double *gyro,double *dvdl)
{
    double T[3],sl[9];

    matcpy(T,gyro,1,3);
    matmul("TN",3,1,3,1.0,Cbe,Omge,-1.0,T);

    skewsym3(T,sl);
    matmul("NN",3,3,3,1.0,Cbr,sl,0.0,dvdl);
}
/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jacobian_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
}
/* sensitive matrix for odometry measurement---------------------------------*/
static int odoHVR(const insopt_t *opt,const double *Cbe,const double *lever,
                  const double *Cbr,const double *ve,const double s,
                  const imud_t *imu,const odod_t *odo,int nx,double *v,
                  double *H,double *R)
{
    int i,j,IV,IA,NV,NA,nv=0;
    double vr[3],r[3],S[9];
    double ds[3],dbg[9],da[9],dma[9],dl[9],ddv[9],dap[3];

    trace(3,"odoHVR:\n");

    ibg=xiBg(opt); nbg=xnBg(opt);
    ios=xiOs(opt); nos=xnOs(opt);
    iol=xiOl(opt); nol=xnOl(opt);
    ioa=xiOa(opt); noa=xnOa(opt);

    IV=xiV(opt); NV=xnV(opt);
    IA=xiA(opt); NA=xnA(opt);

    ve2vr(ve,Cbe,Cbr,lever,imu,s,vr);

    jacobian_ov_ds (vr,ds);
    jacobian_ov_da (ve,Cbe,Cbr,lever,da);
    jacobian_ov_dma(vr,Cbe,Cbr,dma);
    jacobian_ov_bg (Cbe,Cbr,lever,dbg);
    jacobian_ov_dv (Cbe,Cbr,ddv);
    jacobian_ov_lever(Cbe,Cbr,imu->gyro,dl);

#if UPD_IN_EULER
    jacobian_prot_pang(Cbe,S);
    matcpy(dap,da,1,3);
    matmul("NN",1,3,3,1.0,dap,S,0.0,da);
#endif

    for (i=0;i<(opt->odopt.all?3:1);i++) {

        if (fabs(v[nv]=(odo->vr[i]-vr[i]))>MAXINOV) {
            trace(2,"too large innovations for velocity\n");
        }
        if (H) {
            for (j=ios;j<ios+nos;j++) H[j+nv*nx]=ds [i];
            for (j=ioa;j<ioa+noa;j++) H[j+nv*nx]=dma[i+(j-ioa)*3];
            for (j=iol;j<iol+nol;j++) H[j+nv*nx]=dl [i+(j-iol)*3];
            for (j=ibg;j<ibg+nbg;j++) H[j+nv*nx]=dbg[i+(j-ibg)*3];
            for (j=IA ;j<IA+NA  ;j++) H[j+nv*nx]=da [i+(j-IA )*3];
            for (j=IV ;j<IV+NV  ;j++) H[j+nv*nx]=ddv[i+(j-IV )*3];
        }
        r[nv]=VARVEL;
        nv++;
    }
    if (nv) {
        for (i=0;i<nv;i++) R[i+i*nv]=r[i];
        trace(3,"R=\n"); tracemat(5,R,nv,nv,12,6);
        trace(3,"v=\n"); tracemat(5,v,nv,1,12,6);
        trace(3,"H=\n"); tracemat(5,H,nx,nv,15,6);
    }
    return nv;
}
/* close loop for ins states-------------------------------------------------*/
static void odoclp(const double *x,const insopt_t *opt,insstate_t *ins)
{
    int i;

    if (nos) ins->os-=x[ios];
    if (noa) {
        corratt(x+ioa,ins->Cbr);
    }
    for (i=iol;i<iol+nol;i++) ins->rbl[i-iol]+=x[i];
    clp(ins,opt,x);
}
/* odometry velocity measurement update ins states---------------------------*/
static int odofilt(const insopt_t *opt,const imud_t *imu,const odod_t *odo,
                   insstate_t *ins)
{
    int nx=ins->nx,info=0,nv;
    double *v,*H,*R,*x;

    trace(3,"odofilt:\n");

    if (norm(imu->gyro,3)>30.0*D2R) {
        trace(2,"filter fail due to vehicle have large turn\n");
        return 0;
    }
    v=zeros(3,1); H=zeros(3,nx);
    R=zeros(3,3); x=zeros(1,nx);

    matcpy(ins->vr,odo->vr,1,3);

    if ((nv=odoHVR(opt,ins->Cbe,ins->rbl,ins->Cbr,ins->ve,ins->os,
                   imu,odo,nx,v,H,R))) {

        /* ekf filter */
        info=filter(x,ins->P,H,v,R,nx,nv);

        /* solution fail */
        if (info) {
            trace(2,"odometry aid ins fail\n");
            info=0;
        }
        else {
            trace(3,"odometry aid ins ok\n");

            /* filter ok */
            ins->stat=INSS_ODO;
            info=1;
            odoclp(x,opt,ins);
        }
    }
    free(v); free(H); free(R); free(x);
    return info;
}
/* adjust odometry measurement data----------------------------------------*/
static int odomeas(const insopt_t *opt,const odod_t *odo,odod_t *odom)
{
    static double dt=0.0,dr=0.0;

    dt+=odo->dt; dr+=odo->dr;

    if (dt>(opt->odopt.odt==0.0?
            0.5:opt->odopt.odt)) {
        odom->dt=dt;
        odom->dr=dr;
        odom->vr[0]=dr/dt;
        dt=dr=0.0; return 1;
    }
    return 0;
}
/* odometry velocity measurement aid ins updates states----------------------
 * args   :  insopt_t *opt    I   ins options
 *           imud_t *imu      I   imu measurement data
 *           odod_t *odo      I   odometry measurement data
 *           insstate_t *ins  IO  ins states
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int odo(const insopt_t *opt,const imud_t *imu,const odod_t *odo,
               insstate_t *ins)
{
    odod_t odom={0};

    trace(3,"odo:\n");

    if (odomeas(opt,odo,&odom)) {
        return odofilt(opt,imu,&odom,ins);
    }
    return 0;
}