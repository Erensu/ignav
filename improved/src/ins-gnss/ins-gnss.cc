/*------------------------------------------------------------------------------
* ins-gnss.cc : ins-gnss loosely coupled common functions
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
* history : 2017/10/02 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants/macros ----------------------------------------------------------*/
#define MAXDT       3600.0                 /* max time difference for ins-gnss coupled */
#define UNC_ATT     (10.0*D2R)             /* default initial attitude variance */
#define UNC_VEL     (30.0)                 /* default initial velocity variance */
#define UNC_POS     (30.0)                 /* default initial position variance */
#define UNC_BA      (1000.0*Mg2M)          /* default initial accl bias variance */
#define UNC_BG      (10.0*D2R/3600.0)      /* default initial gyro bias uncertainty
                                            * per instrument (deg/hour, converted to rad/sec) */
#define UNC_DT      (0.5)                  /* default initial time synchronization error variance (s) */
#define UNC_SG      (1E-4)                 /* default initial residual scale factors of gyroscopes variance */
#define UNC_SA      (1E-4)                 /* default initial residual scale factors of gyroscopes variance */
#define UNC_RG      (1E-4)                 /* default initial non-orthogonal between sensor axes for gyro */
#define UNC_RA      (1E-4)                 /* default initial non-orthogonal between sensor axes for accl */
#define UNC_LEVER   (1.0)                  /* default initial lever arm for body to ant. uncertainty (m) */
#define UNC_OS      (0.01)                 /* default initial odometry scale factor uncertainty */
#define UNC_OA      (5.0*D2R)              /* default initial odometry misalignment uncertainty */
#define UNC_CLK     (100.0)                /* default initial receiver clock uncertainty (m) */
#define UNC_CLKR    (10.0)                 /* default initial receiver clock drift uncertainty (/s) */
#define UNC_CMA     (30.0*D2R)             /* default initial misalignment from camera to imu body */
#define UNC_VMA     (30.0*D2R)             /* default initial misalignment from v-frame to b-frame */
#define UNC_LCM     (0.5)                  /* default initial lever arm from camera to imu body uncertainty (m) */
#define UNC_CFO     (10.0)                 /* default initial camera calibration parameters: fx,fu,ox,oy uncertainty (pixel) */
#define UNC_CKP     (1.0)                  /* default initial camera calibration parameters: k1,k2,p1,p2 uncertainty */

#define NMP         3                      /* number of position measurements for ins-gnss coupled */
#define NMV         3                      /* number of velocity measurements for ins-gnss coupled */
#define NM          (NMP+NMV)              /* number of all measurements for ins-gnss coupled */
#define IMP         0                      /* index of position measurements */
#define IMV         NMP                    /* index of velocity measurements */

#define NNAC        3                      /* number of accl process noise */
#define NNGY        3                      /* number of gyro process noise */
#define NNBA(opt)   ((opt)->baopt==INS_BAEST?3:0)    /* number of accl bias process noise */
#define NNBG(opt)   ((opt)->bgopt==INS_BGEST?3:0)    /* number of gyro bias process noise */
#define NNPX(opt)   (NNBA(opt)+NNBG(opt)+NNAC+NNGY)
#define INAC        0                      /* index of accl process noise */
#define INGY        NNAC                   /* index of gyro process noise */
#define INBA        (NNAC+NNGY)            /* index of accl bias process noise */
#define INBG(opt)   (NNAC+NNGY+NNBA(opt))  /* index of gyro bias process noise */

#define WGS_E       0.0818191908425        /* WGS84 eccentricity */
#define MATLAB      0                      /* MATLAB generate jacobians for gnss-velocity measurements by attitude */
#define STD_POS     2.5                    /* pos-std-variance of gnss positioning (m)*/
#define STD_VEL     0.1                    /* vel-std-variance of gnss positioning (m/s) */

#define MAXINOP         1000.0             /* max innovations for updates ins states */
#define MAXINOV         100.0              /* max innovations for updates ins states */
#define MAXVEL          5.0                /* max velocity for refine attitude */
#define MINVEL          3.0                /* min velocity for refine attitude */
#define MAXGYRO         (5.0*D2R)          /* max gyro measurement for refine attitude */
#define MAXANG          3.0                /* max angle for refine attitude */
#define MAXUPDTIMEINT   60.0               /* max update time internal is same as correlation time */
#define MAXSYNDIFF      1.0                /* max time difference of ins and gnss time synchronization */
#define CORRETIME       360.0              /* correlation time for gauss-markov process */
#define MAXROT          (10.0*D2R)         /* max rotation of vehicle when velocity matching alignment */
#define MAXSOLS         5                  /* max number of solutions for reboot lc  */
#define MAXDIFF         10.0               /* max time difference between solution */
#define MAXVARDIS       (10.0)             /* max variance of disable estimated state */
#define REBOOT          0                  /* ins loosely coupled reboot if always update fail */
#define USE_MEAS_COV    0                  /* use measurement covariance {xx,yy,zz,xy,xz,yz} for update, not just {xx,yy,zz} */
#define CHKNUMERIC      1                  /* check numeric for given value */
#define NOINTERP        0                  /* no interpolate ins position/velocity when gnss measurement if need */
#define COR_IN_ROV      0                  /* correction attitude in rotation vector,otherwise in euler angles */
#define UPD_INS_E       0                  /* updates ins states in e-frame */

/* global states index -------------------------------------------------------*/
static int IA=0,NA=0;                      /* index and number of attitude states */
static int IV=0,NV=0;                      /* index and number of velocity states */
static int IP=0,NP=0;                      /* index and number of position states */
static int iba=0,nba=0;                    /* index and number of accl bias states */
static int ibg=0,nbg=0;                    /* index and number of gyro bias states */
static int idt=0,ndt=0;                    /* index and number of ins-gnss time synchronization error */
static int isa=0,nsa=0;                    /* index and number of residual scale factors of accl states */
static int isg=0,nsg=0;                    /* index and number of residual scale factors of gyroscopes states */
static int irg=0,nrg=0;                    /* index and number of non-orthogonal between sensor axes for gyro */
static int ira=0,nra=0;                    /* index and number of non-orthogonal between sensor axes for accl */
static int ila=0,nla=0;                    /* index and number of lever arm for body to ant. */
static int iso=0,nos=0;                    /* index and number of scale factor for odometry */
static int iol=0,nol=0;                    /* index and number of lever arm for odometry */
static int ioa=0,noa=0;                    /* index and number of misalignment between odometry and imu body frame */
static int irc=0,nrc=0;                    /* index and number of receiver clock state */
static int irr=0,nrr=0;                    /* index and number of receiver clock drift state */
static int icm=0,ncm=0;                    /* index and number of misalignment from camera to imu body */
static int ivm=0,nvm=0;                    /* index and number of misalignment from b-frame to v-frame */
static int icl=0,ncl=0;                    /* index and number of lever-arm from camera to imu body */
static int ifo=0,nfo=0;                    /* index and number of camera calibration parameters: fx,fy,ox,oy */
static int ikp=0,nkp=0;                    /* index and number of camera calibration parameters: k1,k2,p1,p2 */

/* do nothing----------------------------------------------------------------*/
extern void none(void) {}
/* system noise covariance matrix--------------------------------------------*/
static void sysQ(int is,int n,int nx,double v,double dt,double *Q)
{
    int i; for (i=is;i<is+n;i++) Q[i+i*nx]=v*fabs(dt);
}
/* determine approximate system noise covariance matrix ---------------------*/
static void getQ(const insopt_t *opt,double dt,double* Q)
{
    int nx=xnX(opt);

    trace(3,"getQ:\n");

    if (fabs(dt)>=MAXDT) {
        trace(3,"too large time difference of ins and gnss measurements\n");
    }
    setzero(Q,nx,nx);

    sysQ(IA,NA,nx,opt->psd.gyro,dt,Q);
    sysQ(IV,NV,nx,opt->psd.accl,dt,Q);
    sysQ(iba,nba,nx,opt->psd.ba,dt,Q);
    sysQ(ibg,nbg,nx,opt->psd.bg,dt,Q);
    sysQ(irc,nrc,nx,opt->psd.clk ,dt,Q);
    sysQ(irr,nrr,nx,opt->psd.clkr,dt,Q);
}
/* process noise gain matrix-------------------------------------------------*/
static void getGn(const insopt_t *opt,const insstate_t *ins,const double dt,
                  double *Gn)
{
    int nx=xnX(opt),nprn=NNPX(opt);
    double *I=eye(3);

    setzero(Gn,nx,nprn);
    asi_blk_mat(Gn,nx,nprn,ins->Cbe,3,3,0,9);
    asi_blk_mat(Gn,nx,nprn,ins->Cbe,3,3,3,0);

    opt->baopt==INS_BAEST?asi_blk_mat(Gn,nx,nprn,I,3,3,9 ,6):none();
    opt->bgopt==INS_BGEST?asi_blk_mat(Gn,nx,nprn,I,3,3,12,9):none();
    free(I);
}
/* process noise covariance matrix-------------------------------------------*/
static void getprn(const insstate_t *ins,const insopt_t *opt,double dt,double* Q)
{
    int nprn=NNPX(opt),nx=xnX(opt),i;
    double *Qn=zeros(nprn,nprn),*Gn=mat(nprn,nx);

    for (i=INAC;i<INGY+NNAC;i++) Qn[i+i*nprn]=opt->psd.gyro*fabs(dt);
    for (i=INGY;i<INGY+NNGY;i++) Qn[i+i*nprn]=opt->psd.accl*fabs(dt);
    for (i=INBA     ;i<INBA     +NNBA(opt);i++) Qn[i+i*nprn]=opt->psd.ba*fabs(dt);
    for (i=INBG(opt);i<INBG(opt)+NNBG(opt);i++) Qn[i+i*nprn]=opt->psd.bg*fabs(dt);

    getGn(opt,ins,fabs(dt),Gn);
    matmul33("NNT",Gn,Qn,Gn,nx,nprn,nprn,nx,Q);
    free(Qn); free(Gn);
}
/* functions for initial error covariance matrix ----------------------------*/
extern void initP(int is,int ni,int nx,double unc,double unc0,double *P0)
{
    int i,j;
    for (i=is;i<is+ni;i++) for (j=0;j<nx;j++) {
        if (j==i) P0[j+i*nx]=SQR(unc==0.0?unc0:unc);
        else      P0[j+i*nx]=P0[i+j*nx]=0.0;
    }
}
/* reset index of estimate states--------------------------------------------*/
static void resetindex(void)
{
    iba=nba=ibg=nbg=idt=ndt=isa=nsa=0;
    isg=nsg=irg=nrg=ira=nra=ila=nla=0;
    iso=nos=iol=nol=ioa=noa=irc=nrc=0;
    irr=nrr=0;
    icm=ncm=0;
    ivm=nvm=0;
    icl=ncl=0;
    ifo=nfo=0;
    ikp=nkp=0;
    IA=NA=0; IV=NV=0; IP=NP=0;
}
/* initial error covariance matrix-------------------------------------------*/
extern void getP0(const insopt_t *opt,double *P0)
{
    int nx=xnX(opt),i;

    trace(3,"getP0:\n");

    setzero(P0,nx,nx);

    initP(IA,NA,nx,opt->unc.att,UNC_ATT,P0);
    initP(IV,NV,nx,opt->unc.vel,UNC_VEL,P0);
    initP(IP,NP,nx,opt->unc.pos,UNC_POS,P0);
    initP(iba,nba,nx,opt->unc.ba,UNC_BA,P0);
    initP(ibg,nbg,nx,opt->unc.bg,UNC_BG,P0);
    initP(idt,ndt,nx,opt->unc.dt,UNC_DT,P0);
    initP(isg,nsg,nx,opt->unc.sg,UNC_SG,P0);
    initP(isa,nsa,nx,opt->unc.sa,UNC_SA,P0);
    initP(irg,nrg,nx,opt->unc.rg,UNC_RG,P0);
    initP(ira,nra,nx,opt->unc.ra,UNC_RA,P0);
    initP(ila,nla,nx,opt->unc.lever,UNC_LEVER,P0);
    initP(iso,nos,nx,opt->unc.os,UNC_OS,P0);
    initP(ioa,noa,nx,opt->unc.oa,UNC_OA,P0);
    initP(iol,nol,nx,opt->unc.ol,UNC_LEVER,P0);
    initP(irc,nrc,nx,opt->unc.rc,UNC_CLK,P0);
    initP(irr,nrr,nx,opt->unc.rr,UNC_CLKR,P0);
    initP(icm,ncm,nx,opt->unc.cma,UNC_CMA,P0);
    initP(ivm,nvm,nx,opt->unc.vma,UNC_VMA,P0);
    initP(icl,ncl,nx,opt->unc.lma,UNC_LCM,P0);

    for (i=0;i<4&&nfo;i++) {
        initP(ifo+i,1,nx,opt->unc.fo[i],UNC_CFO,P0);
    }
    for (i=0;i<4&&nkp;i++) {
        initP(ikp+i,1,nx,opt->unc.kp[i],UNC_CKP,P0);
    }
}
/* geocentric radius---------------------------------------------------------*/
extern double georadi(const double *pos)
{
    return RE_WGS84/sqrt(1.0-SQR(WGS_E*sin(pos[0])))*
           sqrt(SQR(cos(pos[0]))+SQR(1.0-SQR(WGS_E))*SQR(sin(pos[0])));
}
/* initial ins-gnss coupled ekf estimated states and it covariance-----------
 * args  :  insopt_t *opt    I  ins options
 *          insstate_t *ins  IO ins states
 * return : none
 * note   : it also can initial ins tightly coupled
 * --------------------------------------------------------------------------*/
extern void initlc(insopt_t *opt,insstate_t *ins)
{
    int i; gtime_t t0={0};

    trace(3,"initlc:\n");

    /* initial global states index and numbers */
    IA=xiA  (opt); NA=xnA  (opt);
    IV=xiV  (opt); NV=xnV  (opt);
    IP=xiP  (opt); NP=xnP  (opt);
    iba=xiBa(opt); nba=xnBa(opt);
    ibg=xiBg(opt); nbg=xnBg(opt);
    isg=xiSg(opt); nsg=xnSg(opt);
    isa=xiSa(opt); nsa=xnSa(opt);
    idt=xiDt(opt); ndt=xnDt(opt);
    ira=xiRa(opt); nra=xnRa(opt);
    irg=xiRg(opt); nrg=xnRg(opt);
    ila=xiLa(opt); nla=xnLa(opt);
    iso=xiOs(opt); nos=xnOs(opt);
    iol=xiOl(opt); nol=xnOl(opt);
    ioa=xiOa(opt); noa=xnOa(opt);
    irc=xiRc(opt); nrc=xnRc(opt);
    irr=xiRr(opt); nrr=xnRr(opt);
    icm=xiCm(opt); ncm=xnCm(opt);
    ivm=xiVm(opt); nvm=xnVm(opt);

    icl=xiCl (opt); ncl=xnCla(opt);
    ifo=xiCfo(opt); nfo=xnCfo(opt);
    ikp=xiCkp(opt); nkp=xnCkp(opt);

    ins->nx=xnX (opt);
    ins->nb=xnRx(opt);
    ins->x =mat(ins->nx,1); ins->P =mat(ins->nx,ins->nx);
    ins->xa=mat(ins->nx,1); ins->Pa=mat(ins->nx,ins->nx);
    ins->xb=mat(ins->nb,1); ins->Pb=mat(ins->nb,ins->nb);
    ins->F =eye  (ins->nx);
    ins->P0=zeros(ins->nx,ins->nx);

    ins->ptime=ins->ptct=ins->plct=t0;
    ins->dtrr=0.0;
    ins->gstat=ins->ns=0;
    ins->gmeas.data=NULL;
    ins->rtkp=NULL;
    ins->gmeas.n=ins->gmeas.nmax=0;

    for (i=0;i<ins->nx;i++) ins->x [i]=0.0;
    for (i=0;i<ins->nx;i++) ins->xa[i]=0.0;
    for (i=0;i<6;i++) ins->dtr[i]=0.0;
    getP0(opt,ins->P);
    getP0(opt,ins->Pa);

    matcpy(ins->lever,opt->lever,3,1);
    ins->len=opt->len;

    matcpy(ins->ba,opt->imuerr.ba,1,3);
    matcpy(ins->bg,opt->imuerr.bg,1,3);
    matcpy(ins->Ma,opt->imuerr.Ma,3,3);
    matcpy(ins->Mg,opt->imuerr.Mg,3,3);
    matcpy(ins->Gg,opt->imuerr.Gg,3,3);

    odores(opt->odopt.res);
    odod(opt->odopt.d);

    if (((prcopt_t*)opt->gopt)->mode>=PMODE_INS_LGNSS_VO) {

        rpy2dcm(opt->voopt.ebc,ins->Cbc);
        matcpy(ins->lbc,opt->voopt.lbc,1,3);

        ins->fx=opt->voopt.calib.fu; ins->fy=opt->voopt.calib.fv;
        ins->ox=opt->voopt.calib.cu; ins->oy=opt->voopt.calib.cv;

        ins->k1=opt->voopt.cam.k1;
        ins->k2=opt->voopt.cam.k2;
        ins->p1=opt->voopt.cam.p1;
        ins->p2=opt->voopt.cam.p2;
        initvoaid  (opt);
        initvoaidlc(opt);
    }
}
/* free ins-gnss coupled ekf estimated states and it covariance--------------
 * args   : insstate_t *ins  IO  ins states
 * return : none
 * note   : it also can initial ins tightly coupled
 * --------------------------------------------------------------------------*/
extern void freelc(insstate_t *ins)
{
    if (ins->x ) free(ins->x ); ins->x =NULL;
    if (ins->P ) free(ins->P ); ins->P =NULL;
    if (ins->Pa) free(ins->Pa); ins->Pa=NULL;
    if (ins->Pb) free(ins->Pb); ins->Pb=NULL;
    if (ins->xb) free(ins->xb); ins->xb=NULL;
    if (ins->xa) free(ins->xa); ins->xa=NULL;
    if (ins->gmeas.data) free(ins->gmeas.data); ins->gmeas.data=NULL;

    ins->nx=ins->nb=0;
    ins->gmeas.n=ins->gmeas.nmax=0;
    resetindex ();
    freevoaid  ();
    freevoaidlc();
}
/* propagate matrix for stochastic parameters--------------------------------*/
static void stochasticPhi(int opt,int ix,int nix,int nx,double dt,double *phi)
{
    int i; if (nix<=0) return;
    for (i=ix;i<ix+nix;i++) {
        if (opt==INS_RANDOM_CONS ) phi[i+i*nx]=1.0;
        if (opt==INS_RANDOM_WALK ) phi[i+i*nx]=1.0;
        if (opt==INS_GAUSS_MARKOV) phi[i+i*nx]=exp(-fabs(dt)/CORRETIME);
    }
}
/* determine transition matrix(Phi=I+F*dt) in euler angle space--------------*/
static void getphi_euler(const insopt_t *opt, double dt, const double *Cbe,
                         const double *pos, const double *omgb,
                         const double *fib, double *phi)
{
    double T1[9],T2[9],T3[18],T4[9],W[18]={0},S[9]={0},Sinv[9];
    double rpy[3],T[9],WC[18]={0},W2[9];
    double fs[9]={0},omgs[9]={0},Ts[9];
    double omega[3],re,rn[3],ge[3];
    int i,j,nx=xnX(opt);

    trace(3,"getphi_euler:\n");

    if (fabs(dt)>MAXDT) {
        trace(3,"large time difference between ins and gnss measurements\n");
    }
    /* initial phi to unit matrix */
    seteye(phi,nx);

    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
    matcpy(Sinv,S,3,3);
    if (matinv(Sinv,3)) {
        trace(2,"transition matrix error\n");
        return;
    }
    matmul("NN",3,3,3,1.0,Sinv,Cbe,0.0,T4);
    matmul("NN",3,3,3,1.0,Omge,Omge,0.0,W2);

    W[0 ]=omgb[1]; W[3] =omgb[2];
    W[7 ]=omgb[0]; W[10]=omgb[2];
    W[14]=omgb[0]; W[17]=omgb[1];
    matmul("NN",3,6,3,1.0,Cbe,W,0.0,WC);
    matmul("NN",3,6,3,1.0,Sinv,WC,0.0,T3);

    matmul("NN",3,3,3,1.0,Sinv,Cbe,0.0,T1);
    matmul33("NNN",Sinv,Omge,S,3,3,3,3,T2);

    omgs[0]=omgb[0];
    omgs[4]=omgb[1];
    omgs[8]=omgb[2];
    matmul("NN",3,3,3,1.0,Cbe,omgs,0.0,Ts);

    /* attitude transmit matrix */
    for (i=IA;i<IA+NA;i++) {
        for (j=IA;j<  IA+ NA;j++) phi[i+j*nx]-=T2[i-IA+(j-IA )*3]*dt;
        for (j=ibg;j<ibg+nbg;j++) phi[i+j*nx] =T1[i-IA+(j-ibg)*3]*dt;
        for (j=isg;j<isg+nsg;j++) phi[i+j*nx] =Ts[i-IA+(j-isg)*3]*dt;
        for (j=irg;j<irg+nrg;j++) phi[i+j*nx] =T3[i-IA+(j-irg)*3]*dt;
    }
    /* velocity transmit matrix */
    matmul3("NN",Cbe,fib,omega);
    skewsym3(omega,T); ecef2pos(pos,rn);
    pregrav(pos,ge); re=georadi(rn);

    matmul("NN",3,3,3,1.0,T,S,0.0,T1);

    W[0 ]=fib[1]; W[3] =fib[2];
    W[7 ]=fib[0]; W[10]=fib[2];
    W[14]=fib[0]; W[17]=fib[1];
    matmul("NN",3,6,3,1.0,Cbe,W,0.0,WC);
    fs[0]=fib[0];
    fs[4]=fib[1];
    fs[8]=fib[2];
    matmul("NN",3,3,3,1.0,Cbe,fs,0.0,Ts);

    for (i=IV;i<IV+NV;i++) {
        for (j=IV; j<IV+NV;  j++) phi[i+j*nx]-=2.0*Omge[i-IV+(j-IV)*3]*dt;
        for (j=IA; j<IA+NA;  j++) phi[i+j*nx] =-T1[i-IV+(j-IA)*3]*dt;
        for (j=IP; j<IP+NP;  j++) phi[i+j*nx] =-2.0*dt/(re*norm(pos,3))*ge[i-IV]*pos[j-IP]-W2[i-IV+(j-IP)*3]*dt;
        for (j=iba;j<iba+nba;j++) phi[i+j*nx] =Cbe[i-IV+(j-iba)*3]*dt;
        for (j=isa;j<isa+nsa;j++) phi[i+j*nx] =Ts [i-IV+(j-isa)*3]*dt;
        for (j=ira;j<ira+nra;j++) phi[i+j*nx] =WC [i-IV+(j-ira)*3]*dt;
    }
    /* position transmit matrix */
    for (i=IP;i<IP+NP;i++) {
        for (j=IV;j<IV+NV;j++) phi[i+j*nx]=(i-IP)==(j-IV)?dt:0.0;
    }
    /* propagate matrix for stochastic parameters */
    stochasticPhi(opt->baproopt,iba,nba,nx,dt,phi);
    stochasticPhi(opt->bgproopt,ibg,nbg,nx,dt,phi);
    stochasticPhi(opt->sgproopt,isg,nsg,nx,dt,phi);
    stochasticPhi(opt->saproopt,isa,nsa,nx,dt,phi);
    stochasticPhi(opt->dtproopt,idt,ndt,nx,dt,phi);
    stochasticPhi(opt->raproopt,ira,nra,nx,dt,phi);
    stochasticPhi(opt->rgproopt,irg,nrg,nx,dt,phi);
    stochasticPhi(opt->osproopt,iso,nos,nx,dt,phi);
    stochasticPhi(opt->olproopt,iol,nol,nx,dt,phi);
    stochasticPhi(opt->oaproopt,ioa,noa,nx,dt,phi);
    stochasticPhi(opt->cmaopt  ,icm,ncm,nx,dt,phi);
    stochasticPhi(opt->vmaopt  ,ivm,nvm,nx,dt,phi);
    stochasticPhi(opt->claopt  ,icl,ncl,nx,dt,phi);

    stochasticPhi(opt->cfoopt  ,ifo,nfo,nx,dt,phi);
    stochasticPhi(opt->ckpopt  ,ikp,nkp,nx,dt,phi);
}
/* determine transition matrix(first-order approx: Phi=I+F*dt)---------------*/
static void getPhi1(const insopt_t *opt, double dt, const double *Cbe,
                    const double *pos, const double *omgb, const double *fib,
                    double *phi)
{
    int i,j,nx=xnX(opt);
    double omega[3]={0},T[9],ge[3],re,rn[3],W[18]={0},WC[18]={0},Cbv[9];
    double fs[9]={0},omgs[9]={0},Ts[9];
    double W2[9];

    trace(3,"getPhi1:\n");

    if (fabs(dt)>MAXDT) {
        trace(3,"large time difference between ins and gnss measurements\n");
    }
    /* initial phi to unit matrix */
    seteye(phi,nx);
    seteye(Cbv,3 );

    /* attitude transmit matrix */
    W[0 ]=omgb[1]; W[3] =omgb[2];
    W[7 ]=omgb[0]; W[10]=omgb[2];
    W[14]=omgb[0]; W[17]=omgb[1];
    matmul("NN",3,6,3,1.0,Cbe,W,0.0,WC);

    omgs[0]=omgb[0];
    omgs[4]=omgb[1];
    omgs[8]=omgb[2];
    matmul("NN",3,3,3,1.0,Cbe,omgs,0.0,Ts);

    for (i=IA;i<IA+NA;i++) {
        for (j= IA;j<IA+NA  ;j++) phi[i+j*nx]-=Omge[i-IA+(j-IA )*3]*dt;
        for (j=ibg;j<ibg+nbg;j++) phi[i+j*nx] =Cbe [i-IA+(j-ibg)*3]*dt;
        for (j=isg;j<isg+nsg;j++) phi[i+j*nx] =Ts  [i-IA+(j-isg)*3]*dt;
        for (j=irg;j<irg+nrg;j++) phi[i+j*nx] =WC  [i-IA+(j-irg)*3]*dt;
    }
    /* velocity transmit matrix */
    matmul3("NN",Cbe,fib,omega);
    skewsym3(omega,T); ecef2pos(pos,rn);
    pregrav(pos,ge); re=georadi(rn);
    matmul("NN",3,3,3,1.0,Omge,Omge,0.0,W2);

    W[0 ]=fib[1]; W[3] =fib[2];
    W[7 ]=fib[0]; W[10]=fib[2];
    W[14]=fib[0]; W[17]=fib[1];
    matmul("NN",3,6,3,1.0,Cbe,W,0.0,WC);
    fs[0]=fib[0];
    fs[4]=fib[1];
    fs[8]=fib[2];
    matmul("NN",3,3,3,1.0,Cbe,fs,0.0,Ts);
    for (i=IV;i<IV+NV;i++) {
        for (j=IV; j<IV+NV;  j++) phi[i+j*nx]-=2.0*Omge[i-IV+(j-IV)*3]*dt;
        for (j=IA; j<IA+NA;  j++) phi[i+j*nx] =-T[i-IV+(j-IA)*3]*dt;
        for (j=IP; j<IP+NP;  j++) phi[i+j*nx] =-2.0*dt/(re*norm(pos,3))*ge[i-IV]*pos[j-IP]-W2[i-IV+(j-IP)*3]*dt;
        for (j=iba;j<iba+nba;j++) phi[i+j*nx] =Cbe[i-IV+(j-iba)*3]*dt;
        for (j=isa;j<isa+nsa;j++) phi[i+j*nx] =Ts [i-IV+(j-isa)*3]*dt;
        for (j=ira;j<ira+nra;j++) phi[i+j*nx] =WC [i-IV+(j-ira)*3]*dt;
    }
    /* position transmit matrix */
    for (i=IP;i<IP+NP;i++) {
        for (j=IV;j<IV+NV;j++) phi[i+j*nx]=(i-IP)==(j-IV)?dt:0.0;
    }
    /* propagate matrix for stochastic parameters */
    stochasticPhi(opt->baproopt,iba,nba,nx,dt,phi);
    stochasticPhi(opt->bgproopt,ibg,nbg,nx,dt,phi);
    stochasticPhi(opt->sgproopt,isg,nsg,nx,dt,phi);
    stochasticPhi(opt->saproopt,isa,nsa,nx,dt,phi);
    stochasticPhi(opt->dtproopt,idt,ndt,nx,dt,phi);
    stochasticPhi(opt->raproopt,ira,nra,nx,dt,phi);
    stochasticPhi(opt->rgproopt,irg,nrg,nx,dt,phi);
    stochasticPhi(opt->osproopt,iso,nos,nx,dt,phi);
    stochasticPhi(opt->olproopt,iol,nol,nx,dt,phi);
    stochasticPhi(opt->oaproopt,ioa,noa,nx,dt,phi);
    stochasticPhi(opt->cmaopt  ,icm,ncm,nx,dt,phi);
    stochasticPhi(opt->vmaopt  ,ivm,nvm,nx,dt,phi);

    stochasticPhi(opt->claopt  ,icl,ncl,nx,dt,phi);
    stochasticPhi(opt->cfoopt  ,ifo,nfo,nx,dt,phi);
    stochasticPhi(opt->ckpopt  ,ikp,nkp,nx,dt,phi);
}
/* system matrix for accl-bias,gyro-bias,sg,sa and so on --------------------*/
static void stochasticF(int opt,int ix,int nix,int nx,double *F)
{
    int i; if (nix<=0) return;
    for (i=ix;i<ix+nix;i++) {
        if (opt==INS_RANDOM_CONS ) F[i+i*nx]=1E-10;
        if (opt==INS_RANDOM_WALK ) F[i+i*nx]=1E-10;
        if (opt==INS_GAUSS_MARKOV) F[i+i*nx]=-1.0/CORRETIME;
    }
}
/* system matrix of ins states propagate in ecef-frame-----------------------*/
static void getF(const insopt_t *opt,const double *Cbe,const double *pos,
                 const double *omgb,const double *fib,double *F)
{
    int i,j,nx=xnX(opt);
    double F21[9],F23[9],I[9]={1,0,0,0,1,0,0,0,1},omega[3],rn[3],ge[3],re;
    double W[18]={0},WC[18]={0};
    double fs[9]={0},omgs[9]={0},Ts[9];
    double W2[9];

    trace(3,"getF:\n");

    setzero(F,nx,nx);

    matmul3("NN",Cbe,fib,omega);
    skewsym3(omega,F21);

    ecef2pos(pos,rn);
    pregrav(pos,ge); re=georadi(rn);
    matmul("NT",3,3,1,-2.0/(re*norm(pos,3)),ge,pos,0.0,F23);
    matmul("NN",3,3,3,1.0,Omge,Omge,0.0,W2);

    /* ins attitude system matrix */
    W[0 ]=omgb[1]; W[3 ]=omgb[2];
    W[7 ]=omgb[0]; W[10]=omgb[2];
    W[14]=omgb[0]; W[17]=omgb[1];
    matmul("NN",3,6,3,1.0,Cbe,W,0.0,WC);

    omgs[0]=omgb[0];
    omgs[4]=omgb[1];
    omgs[8]=omgb[2];
    matmul("NN",3,3,3,1.0,Cbe,omgs,0.0,Ts);

    for (i=IA;i<IA+NA;i++) {
        for (j=IA; j<IA+NA;  j++) F[i+j*nx]=-Omge[i-IA+(j-IA )*3];
        for (j=ibg;j<ibg+nbg;j++) F[i+j*nx]= Cbe [i-IA+(j-ibg)*3];
        for (j=isg;j<isg+nsg;j++) F[i+j*nx]= Ts  [i-IA+(j-isg)*3];
        for (j=irg;j<irg+nrg;j++) F[i+j*nx]= WC  [i-IA+(j-irg)*3];
    }
    /* ins velocity system matrix */
    W[0 ]=fib[1]; W[3] =fib[2];
    W[7 ]=fib[0]; W[10]=fib[2];
    W[14]=fib[0]; W[17]=fib[1];
    matmul("NN",3,6,3,1.0,Cbe,W,0.0,WC);

    fs[0]=fib[0];
    fs[4]=fib[1];
    fs[8]=fib[2];
    matmul("NN",3,3,3,1.0,Cbe,fs,0.0,Ts);
    for (i=IV;i<IV+NV;i++) {
        for (j=IA; j<IA+NA;  j++) F[i+j*nx]=-F21[i-IV+(j-IA )*3];
        for (j=IV; j<IV+NV;  j++) F[i+j*nx]=-2.0*Omge[i-IV+(j-IV )*3];
        for (j=IP; j<IP+NP;  j++) F[i+j*nx]= F23[i-IV+(j-IP )*3]-W2[i-IV+(j-IP)*3];
        for (j=iba;j<iba+nba;j++) F[i+j*nx]= Cbe[i-IV+(j-iba)*3];
        for (j=isa;j<isa+nsa;j++) F[i+j*nx]= Ts [i-IV+(j-isa)*3];
        for (j=ira;j<ira+nra;j++) F[i+j*nx]= WC [i-IV+(j-ira)*3];
    }
    /* ins position system matrix */
    for (i=IP;i<IP+NP;i++) {
        for (j=IV;j<IV+NV;j++) F[i+j*nx]=I[i-IP+(j-IV)*3];
    }
    /* stochastic parameters system matrix */
    stochasticF(opt->baproopt,iba,nba,nx,F);
    stochasticF(opt->bgproopt,ibg,nbg,nx,F);
    stochasticF(opt->sgproopt,isg,nsg,nx,F);
    stochasticF(opt->saproopt,isa,nsa,nx,F);
    stochasticF(opt->dtproopt,idt,ndt,nx,F);
    stochasticF(opt->raproopt,ira,nra,nx,F);
    stochasticF(opt->rgproopt,irg,nrg,nx,F);
    stochasticF(opt->osproopt,iso,nos,nx,F);
    stochasticF(opt->olproopt,iol,nol,nx,F);
    stochasticF(opt->oaproopt,ioa,noa,nx,F);
    stochasticF(opt->cmaopt  ,icm,ncm,nx,F);
    stochasticF(opt->vmaopt  ,ivm,nvm,nx,F);
}
/* propagate state estimation error covariance-------------------------------*/
static void propP(const insopt_t *opt,const double *Q,const double *phi,
                  const double *P0,double *P)
{
    int i,j,nx=xnX(opt);
    double *PQ=mat(nx,nx),*Phi2=mat(nx,nx);

    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) PQ[i+j*nx]=P0[i+j*nx]+0.5*Q[i+j*nx];
    }
    matmul33("NNT",phi,PQ,phi,nx,nx,nx,nx,Phi2);
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) P[i+j*nx]=Phi2[i+j*nx]+0.5*Q[i+j*nx];
    }
    /* initialize every epoch for clock (white noise) */
    initP(irc,nrc,nx,opt->unc.rc,UNC_CLK,P);
    free(PQ); free(Phi2);
}
/* propagate state estimates noting that all states are zero due to closed-loop
 * correction----------------------------------------------------------------*/
static void propx(const insopt_t *opt,const double *x0,double *x)
{
    int i; for (i=0;i<xnCl(opt);i++) x[i]=1E-20;
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
/* jacobian of position measurement by attitude error term ------------------*/
static void jacobian_p_att(const double *Cbe,const double *lever,double *dpdatt)
{
    double cl[3];
    matmul3v("N",Cbe,lever,cl); skewsym3(cl,dpdatt);
}
/* jacobian of position measurement by ins-gnss time synchronization error---*/
static void jacobian_p_dt(const double *omgb,const double *lever,const double *Cbe,
                          const double *ve,double *dpddt)
{
    int i;
    double wl[9],cl[9];
    skewsym3(omgb,wl);
    matmul3v("N",wl,lever,cl);
    matmul3v("N",Cbe,cl,dpddt);
    for (i=0;i<3;i++) dpddt[i]+=ve[i];
}
/* jacobian of velocity measurement by attitude error term ------------------*/
static void jacobian_v_att(const double *Cbe,const double *lever,
                           const double *omgb,double *dvdatt)
{
    int i;
    double cl[9],wl[9],omgie[3];
    
    skewsym3(omgb,cl);
    matmul3("NN",cl,lever,wl);
    matmul3("NN",Cbe,wl,cl);
    
    matmul33("NNN",Omge,Cbe,lever,3,3,3,1,wl);
    for (i=0;i<3;i++) omgie[i]=cl[i]-wl[i];
    skewsym3(omgie,dvdatt);

#if MATLAB
    dvdatt[0]+=(Cbe[2]*lever[0]+Cbe[5]*lever[1]+Cbe[8]*lever[2])*OMGE;
    dvdatt[2]-=(Cbe[0]*lever[0]+Cbe[3]*lever[1]+Cbe[6]*lever[2])*OMGE;
    dvdatt[4]+=(Cbe[2]*lever[0]+Cbe[5]*lever[1]+Cbe[8]*lever[2])*OMGE;
    dvdatt[5]-=(Cbe[1]*lever[0]+Cbe[4]*lever[1]+Cbe[7]*lever[2])*OMGE;
#endif
}
/* jacobian of velocity measurement by accl bias term------------------------*/
static void jacobian_v_bg(const double *Cbe,const double *lever,double *dvdbg)
{
    double cl[9];
    skewsym3(lever,cl);
    matmul3("NN",Cbe,cl,dvdbg);
}
/* jacobian of velocity measurement by ins-gnss time synchronization error term-*/
static void jacobian_v_dt(const double *omgb,const double *Cbe,const double *lever,
                          const double *ae,double *dvddt)
{
    int i;
    double wl[9],cl[9];
    skewsym3(omgb,wl);
    matmul33("NNN",Cbe,wl,wl,3,3,3,3,cl);
    matmul3v("N",cl,lever,dvddt);
    for (i=0;i<3;i++) dvddt[i]+=ae[i];
}
/* jacobian of velocity measurement by residual scale factors of gyroscopes----*/
static void jacobian_v_dsg(const double *Cbe,const double *lever,const double *omgb,
                           double *dvds)
{
    /* dwib=bg+diag(wib)*sg+T_g*r_g (ref[4]) */
    int i;
    double dvdbg[9],wib[9]={0};
    
    jacobian_v_bg(Cbe,lever,dvdbg);
    for (i=0;i<3;i++) wib[i+i*3]=omgb[i];
    matmul3("NN",dvdbg,wib,dvds);
}
/* jacobian of velecity measurement by non-orthogonal between sensor axes of gyro.*/
static void jacobian_v_drg(const double *Cbe,const double *lever,const double *omgb,
                           double *dvdrg)
{
    /* dwib=bg+diag(wib)*sg+T_g*r_g (ref[4]) */
    double T[18]={0},dvdbg[9];

    jacobian_v_bg(Cbe,lever,dvdbg);
    T[0 ]=omgb[1]; T[3] =omgb[2];
    T[7 ]=omgb[0]; T[10]=omgb[2];
    T[14]=omgb[0]; T[17]=omgb[1];
    matmul("NN",3,6,3,1.0,dvdbg,T,0.0,dvdrg);
}
/* jacobian of velocity measurement by lever arm--------------------------------*/
static void jacobian_v_dla(const double *Cbe,const double *omgb,double *dvdla)
{
    double T[9];
    matmul3("NN",Omge,Cbe,dvdla);
    skewsym3(omgb,T);
    matmul("NN",3,3,3,-1.0,Cbe,T,1.0,dvdla);
}
/* measurement sensitive-matrix----------------------------------------------
 * set-up measurement sensitive mtarix
 * args       :  double *pos    I  imu-position (ecef) (m)
 *               double *Cbe    I  bidy-frame to ecef-frame tansformation matrix
 *               double *lever  I  gnss-antenna to body-frame lever arm (m)
 *               double *omgb   I  imu gyro measurements (corrected)
 *               double *fib    I  imu accl measurements (corrected)
 *               double *meas   I  measurements from gnss positioning {(0,2):pos,
 *                                 (3,5):vel}
 *               double *re     I  position of body-frame
 *               double *vr     I  velocity of body-frame
 *               double *std    I  position/velocity std from gnss positioning
 *               double *cov    I  position/velocity covariance matrix 
 *               double *v      O  measurement innovations
 *               double *R      O  measurement variance matrix
 *               double *H      O  sensitive-matrix about estimated states
 * return     :  int            O  numbers of rows of H
 * -------------------------------------------------------------------------*/
static int build_HVR(const insopt_t *opt,const double *pos,const double *Cbe,
                     const double *lever,const double *omgb,const double *fib,
                     const double *meas,const double *re,const double *ve,
                     const double *ae,const double *std,const double *cov,
                     const double *P,double *H,double *v,double *R)
{
    int i,j,nm=0,nx=xnX(opt),ind[NM];
    double r1[9],v1[9],v5[9],*I=eye(3),R_[NM]={0};
    double re_[3],ve_[3];
    double dt1[3],dt2[3],ds[9],drg[18],dla[9];
    double S[9],r1p[9];

    trace(3,"build_HVR:\n");

    /* remove lever-arm effects */
    rmlever(NULL,re,ve,lever,Cbe,omgb,re_,ve_);

    /* for position measurement */
    jacobian_p_att(Cbe,lever,r1);
    jacobian_p_dt (omgb,lever,Cbe,ve,dt1);

#if UPD_IN_EULER
    /* jacobian of perturb rotation wrt. perturb euler angle */
    jacobian_prot_pang(Cbe,S);

    matcpy(r1p,r1,3,3);
    matmul("NN",3,3,3,1.0,r1p,S,0.0,r1);
#endif
    for (i=IMP;i<IMP+NMP;i++) {
        if (meas[i]!=0.0) {
            if (fabs(v[nm]=(meas[i]-re_[i-IMP]))>=MAXINOP) {
                trace(2,"too large innovations for position\n");
            }
            if (H) {
                for (j=IA; j<IA+NA;  j++) H[j+nm*nx]= r1 [i-IMP+(j-IA)*3];
                for (j=IP; j<IP+NP;  j++) H[j+nm*nx]=-I  [i-IMP+(j-IP)*3];
                for (j=idt;j<idt+ndt;j++) H[j+nm*nx]= dt1[i-IMP];
                for (j=ila;j<ila+nla;j++) H[j+nm*nx]=-Cbe[i-IMP+(j-ila)*3];
            }
            R_[nm]=std[i]==0.0?SQR(STD_POS):SQR(std[i]);
            ind[nm++]=i;
        }
    }
    /* for velocity measurement */
    jacobian_v_bg (Cbe,lever,v5);
    jacobian_v_att(Cbe,lever,omgb,v1);
    jacobian_v_dt (omgb,Cbe,lever,ae,dt2);
    jacobian_v_dsg(Cbe,lever,omgb,ds);
    jacobian_v_drg(Cbe,lever,omgb,drg);
    jacobian_v_dla(Cbe,omgb,dla);

#if UPD_IN_EULER
    /* jacobian of perturb rotation wrt. perturb euler angle */
    jacobian_prot_pang(Cbe,S);
    
    matcpy(r1p,v1,3,3);
    matmul("NN",3,3,3,1.0,r1p,S,0.0,v1);
#endif
    for (i=IMV;i<IMV+NMV;i++) {
        if (meas[i]!=0.0) {
            if (fabs(v[nm]=(meas[i]-ve_[i-IMV]))>=MAXINOV) {
                trace(2,"too large innovations for velocity\n");
            }
            if (H) {
                for (j=IA; j<IA+NA;  j++) H[j+nm*nx]=v1 [i-IMV+(j-IA )*3];
                for (j=IV; j<IV+NV;  j++) H[j+nm*nx]=-I [i-IMV+(j-IV )*3];
                for (j=ibg;j<ibg+nbg;j++) H[j+nm*nx]=v5 [i-IMV+(j-ibg)*3];
                for (j=idt;j<idt+ndt;j++) H[j+nm*nx]=dt2[i-IMV];
                for (j=isg;j<isg+nsg;j++) H[j+nm*nx]=ds [i-IMV+(j-isg)*3];
                for (j=irg;j<irg+nrg;j++) H[j+nm*nx]=drg[i-IMV+(j-irg)*3];
                for (j=ila;j<ila+nla;j++) H[j+nm*nx]=dla[i-IMV+(j-ila)*3];
            }
            R_[nm]=std[i]==0.0?SQR(STD_VEL):SQR(std[i]);
            ind[nm++]=i;
        }
    }
    if (nm&&R) {
        if (cov) {
            for (i=0;i<nm;i++) for (j=0;j<nm;j++) {
                R[i+nm*j]=cov[ind[i]+NM*ind[j]];
            }
        }
        else {
            for (i=0;i<nm;i++) R[i+i*nm]=R_[i];
        }
        trace(3,"R=\n");
        tracemat(5,R,nm,nm,15,8);
    }
    if (nm&&H) {
        trace(3,"H=\n");
        tracemat(5,H,nx,nm,15,8);
    }
    if (nm&&v) {
        trace(3,"v=\n");
        tracemat(5,v,nm,1,15,8);
    }
    free(I); return nm;
}
/* validation of ins-gnss coupled solution-----------------------------------*/
static int valsol(double *x,double *P,const double *R,const double *v,
                  int nv,double thres)
{
    int i=0,j;
    double fact=thres*thres;

    trace(3,"valsol:nv=%d\n",nv);

    /* check estimated states */
    if (     x[  0]==DISFLAG&&norm(x+  0,3)> 5.0*D2R) i|=1;
    if (nba&&x[iba]==DISFLAG&&norm(x+iba,3)>1E5*Mg2M) i|=1;
    if (nbg&&x[ibg]==DISFLAG&&norm(x+ibg,3)>10.0*D2R) i|=1;
    if (i) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    /* post-fit residual test */
    for (j=0,i=0;i<nv&&thres;i++) {
        if (v[i]*v[i]<fact*R[i+i*nv]) continue;
        trace(2,"large residual (v=%6.3f sig=%.3f)\n",v[i],SQRT(R[i+i*nv]));
        j++;
    }
    return j<nv;
}
/* updates ins states--------------------------------------------------------*/
static void updinss(insstate_t *ins,const double *re,const double *ve,
                    const double *Cbe,const double *fib,const double *omgb,
                    const double *ba,const double *bg,const double *Mgc,
                    const double *Mac,const double *leverc)
{
    trace(3,"updinss:\n");

    quat_t q;
    matcpy(ins->re,re,1,3); matcpy(ins->ve,ve,1,3);
    matcpy(ins->ba,ba,1,3); matcpy(ins->bg,bg,1,3);

    matcpy(ins->fb,fib,1,3);
    matcpy(ins->omgb,omgb,1,3);
    matcpy(ins->Cbe,Cbe,1,9);
#if 1
    /* normalization dcm */
    normdcm(ins->Cbe);
#else
    dcm2quat(ins->Cbe,&q);
    quat_normalize_self(&q);
    quat2dcm(&q,ins->Cbe);
#endif
    matcpy(ins->Ma,Mac,1,9);
    matcpy(ins->Mg,Mgc,1,9);
    matcpy(ins->lever,leverc,1,3);
}
/* interpolate ins position/velocity when gnss measurement time and ins's -----
 * is different
 * --------------------------------------------------------------------------*/
static void interppv(const insstate_t *ins,double dt,double *pe, double *ve,
                     double *ae,double *Cbe)
{
    quat_t q1,q2,q; int i;
    double tt=timediff(ins->time,ins->ptime),rpy[3]={0},t=(dt+tt)/tt;

    trace(3,"interppv:\n");
    
#if NOINTERP
    t=1.0; dt=0.0;
#endif
    if (fabs(tt)<1E-6) {
        matcpy(pe,ins->re,1,3); matcpy(ve ,ins->ve ,1,3);
        matcpy(ae,ins->ae,1,3); matcpy(Cbe,ins->Cbe,3,3);
        return;
    }
    for (i=0;i<3;i++) {
        pe[i]=ins->re[i]+(ins->re[i]-ins->pins[i  ])/tt*dt;
        ve[i]=ins->ve[i]+(ins->ve[i]-ins->pins[i+3])/tt*dt;
        ae[i]=ins->ae[i]+(ins->ae[i]-ins->pins[i+6])/tt*dt;
    }
    /* attitude-quaternion of current time */
    dcm2rpy(ins->pCbe,rpy);
    rpy2quat(rpy,&q1);
    
    /* attitude-quaternion of precious time */
    dcm2rpy(ins->Cbe,rpy);
    rpy2quat(rpy,&q2);

    /* interpolate ins attitude by quaternion */
    quat_slerp(&q,&q1,&q2,t);
    quat_to_rh_rot_matrix(&q,Cbe);
}
/* copy ins parameters for backup--------------------------------------------*/
static void prepara(const insstate_t *ins,double *fib,double *omgb,double *Mg,
                    double *Ma,double *Gg,double *ba,double *bg,double *lever)
{
    /* copy imu accl and gyro measurements data */
    matcpy(fib,ins->fb,1,3);
    matcpy(omgb,ins->omgb,1,3);

    matcpy(Mg,ins->Mg,3,3);
    matcpy(Ma,ins->Ma,3,3);
    matcpy(Gg,ins->Gg,3,3);

    /* copy accl and gyro bias */
    matcpy(ba,ins->ba,1,3);
    matcpy(bg,ins->bg,1,3);

    matcpy(lever,ins->lever,1,3);
}
/* get the acceleration of body in ecef-frame by input acceleration measurements
 * args     :  double *fib    I  input acceleration measurements
 *             double *Cbe    I  body-frame to ecef-frame convert matrix
 *             double *re     I  position of imu-body in ecef-frame
 *             double *ve     I  velocity of imu-body in ecef-frame
 *             double *ae     O  acceleration of imu-body in ecef-frame
 * return :none
 * --------------------------------------------------------------------------*/
extern void getaccl(const double *fib,const double *Cbe,const double *re,
                    const double *ve,double *ae)
{
    int i; double fe[3],ge[3],cori[3];

    matmul3v("N",Cbe,fib,fe);
    pregrav(re,ge);
    matmul3v("N",Omge,ve,cori);

    for (i=0;i<3;i++) {
        ae[i]=fe[i]+ge[i]-2.0*cori[i];
    }
}
/* close-loop for states ----------------------------------------------------*/
extern void lcclp(double *x,double *Cbe,double *re,double *ve,double *fib,
                  double *omgb,double *Gg, double *rec,double *vec,double *aec,
                  double *bac, double *bgc,double *Mac,double *Mgc,
                  double *leverc,double *Cbec,double *fibc,double *omgbc,
                  const insopt_t *opt)
{
    double T[9],*I=eye(3),ang[3];
    int i;

    trace(3,"lcclp:\n");

#if COR_IN_ROV
    /* attitude correction */
    if (x[0]!=DISFLAG) {
        skewsym3(x,T); for (i=0;i<9;i++) I[i]-=T[i];
        matmul3("NN",I,Cbe,Cbec);
    }
#else
    /* attitude correction */
    if (x[0]!=DISFLAG) {
        correctatt(x,Cbe,Cbec);
    }
#endif

#if UPD_IN_EULER
    c2rpy(Cbe,ang);
    for (i=0;i<3;i++) ang[i]=ang[i]+x[i];
    rpy2c(ang,Cbec);
#endif
    /* velocity and position correction */
    vec[0]=ve[0]-x[IV+0];
    vec[1]=ve[1]-x[IV+1];
    vec[2]=ve[2]-x[IV+2];

    rec[0]=re[0]-x[IP+0];
    rec[1]=re[1]-x[IP+1];
    rec[2]=re[2]-x[IP+2];

    /* accl and gyro bias correction */
    if (nba&&x[iba]!=DISFLAG) {
        bac[0]+=x[iba]; bac[1]+=x[iba+1]; bac[2]+=x[iba+2];
    }
    if (nbg&&x[ibg]!=DISFLAG) {
        bgc[0]+=x[ibg]; bgc[1]+=x[ibg+1]; bgc[2]+=x[ibg+2];
    }
    /* residual scale factors of gyroscopes and accl correction */
    if (nsg&&x[isg]!=DISFLAG) {
        for (i=isg;i<isg+nsg;i++) Mgc[i-isg+(i-isg)*3]+=x[i];
    }
    if (nsa&&x[isa]!=DISFLAG) {
        for (i=isa;i<isa+nsa;i++) Mac[i-isa+(i-isa)*3]+=x[i];
    }
    /* correction non-orthogonal between sensor axes */
    if (opt->estrg&&x[irg]!=DISFLAG) {
        Mgc[3]+=x[irg  ]; Mgc[6]+=x[irg+1]; Mgc[1]+=x[irg+2];
        Mgc[7]+=x[irg+3]; Mgc[2]+=x[irg+4]; Mgc[5]+=x[irg+5];
    }
    if (opt->estra&&x[ira]!=DISFLAG) {
        Mac[3]+=x[ira  ]; Mac[6]+=x[ira+1]; Mac[1]+=x[ira+2];
        Mac[7]+=x[ira+3]; Mac[2]+=x[ira+4]; Mac[5]+=x[ira+5];
    }
    /* correction for lever arm */
    if (nla&&x[ila]!=DISFLAG) {
        leverc[0]+=x[ila+0];
        leverc[1]+=x[ila+1];
        leverc[2]+=x[ila+2];
    }
    /* correction imu accl and gyro measurements */
    if (fib&&omgb&&fibc&&omgbc&&Gg) {
        ins_errmodel2(fib,omgb,Mac,Mgc,bac,bgc,Gg,fibc,omgbc);

        /* correction imu-body accelerometer */
        getaccl(fibc,Cbec,rec,vec,aec);
    }
    free(I);
}
/* loosely-coupled INS/GNSS integration--------------------------------------
 * args  : insopt_t *opt     I  ins-gnss coupled options
 *         imudata_t *data   I  measurements from imu (corrected)
 *         insstate_t *ins   IO ins states from measurements updates
 *         double *meas      I  measurements from gnss positioning
 *         double *std       I  measurements std from gnss positioning
 *         double *cov       I  measurements covariance matrix
 *         double dt         I  time difference of gnss and IMU measurement data
 *                              (dt=t_gnss-t_imu)
 *         double *x0        I  propagate state estimates
 *         double *P0        I  propagate state estimates covariance
 * return: int               O  1:ok,0:failed
 * note  : if meas[i] is 0.0,then means no measurements to update
 *         if std [i] is 0.0,then means no measurements to update
 * --------------------------------------------------------------------------*/
static int lcfilt(const insopt_t *opt, insstate_t *ins, const double *meas,const double *std,
                  const double *cov,const double dt,
                  double *x, double *P)
{
    int nm,info=0,stat,nx=ins->nx,i;
    double *H,*v,*R;
    double re[3],ve[3],ae[3],Cbe[9],fib[3],omgb[3],stde[NM]={0};
    
    /* close-loop correction states */
    double rec[3],vec[3],aec[3],Cbec[9],fibc[3],omgbc[3],
           bac[3],bgc[3],Mac[9],Mgc[9],Gg[9],leverc[3];

    trace(3,"lcfiltr: dt=%.3lf\n",dt);

    /* check gnss and ins time synchronization */
    if (fabs(ins->age=dt)>=MAXSYNDIFF) {
        trace(2,"gnss and ins time synchronization fail,dt=%5.3lf\n",dt);
        return 0;
    }
    /* interpolate ins states: attitude,position and velocity */
    interppv(ins,dt,re,ve,ae,Cbe);

    prepara(ins,fib,omgb,Mgc,Mac,Gg,bac,bgc,leverc);

    /* build H,v and R matrix from input measurements */
    H=zeros(NM,nx); v=zeros(NM,1);
    R=zeros(NM,NM);

    if ((nm=build_HVR(opt,NULL,Cbe,leverc,omgb,fib,meas,re,ve,ae,std,cov,P,H,v,R))) {
        if (cov) {
            for (i=0;i<NM;i++) stde[i]=cov[i+i*NM];
        }
        else {
            matcpy(stde,std,1,NM);
        }
        /* disable gyro/accl/att state if measurement is bad */
        if (norm(stde,NM)>=MAXVARDIS) {
            for (i=IA ;i<IA+NA  ;i++) unusex(opt,i,ins,x);
            for (i=iba;i<iba+nba;i++) unusex(opt,i,ins,x);
            for (i=ibg;i<ibg+nbg;i++) unusex(opt,i,ins,x);
        }
        /* ekf filter */
        if ((info=filter(x,P,H,v,R,nx,nm))) {
            trace(2,"filter error (info=%d)\n",info);
            free(H); free(v); free(R);
            return 0;
        }
        trace(3,"dx=\n");
        tracemat(3,x,1,nx,12,5);
    }
    else {
        trace(3,"no gnss position and velocity measurement\n");
        free(H); free(v); free(R);
        return 0;
    }
    /* close-loop for estimated states */
    lcclp(x,Cbe,re,ve,fib,omgb,Gg,rec,vec,aec,bac,bgc,Mac,Mgc,
          leverc,Cbec,fibc,omgbc,opt);

    /* post-fit residuals for ins-gnss coupled */
    if ((nm=build_HVR(opt,NULL,Cbec,leverc,omgbc,fibc,meas,rec,vec,aec,std,cov,P,H,v,R))) {

        /* validation of solutions */
        if ((stat=valsol(x,P,R,v,nm,10.0))) {
            
            /* updates ins states */
            updinss(ins,rec,vec,Cbec,fibc,omgbc,bac,bgc,
                    Mgc,Mac,leverc);
        }
    }
    else {
        stat=0; /* update fail */
    }
    free(H); free(v); free(R);
    return stat;
}
/* precise system propagate matrix-------------------------------------------*/
static void precPhi(const insopt_t *opt,double dt,const double *Cbe,
                    const double *pos,const double *omgb,const double *fib,
                    double *Phi)
{
    int i,nx=xnX(opt);
    double *F=zeros(nx,nx),*FF=zeros(nx,nx);
    double *FFF=zeros(nx,nx),*I=eye(nx);

    getF(opt,Cbe,pos,omgb,fib,F);
    for (i=0;i<nx*nx;i++) F[i]*=dt;
#if 0
    /* third-order approx */
    matmul("NN",nx,nx,nx,1.0,F,F,0.0,FF);
    matmul33("NNN",F,F,F,nx,nx,nx,nx,FFF);

    for (i=0;i<nx*nx;i++) {
        Phi[i]=I[i]+F[i]+0.5*FF[i]+1.0/6.0*FFF[i];
    }
#else
    expmat(F,nx,Phi);
#endif
    free(F); free(FF);
    free(FFF); free(I);
}
/* get system propagate matrix-----------------------------------------------
 * args:    insopt_t *opt   I  ins options
 *          insstate_t *ins I  ins states
 *          double *Phi     O  system propagate matrix
 * return: none
 * -------------------------------------------------------------------------*/
extern void getPhi(const insopt_t *opt,const insstate_t *ins,double *Phi)
{
    return precPhi(opt,ins->dt,ins->Cbe,ins->re,ins->omgb,ins->fb,Phi);    
}
/* get system noise matrix---------------------------------------------------
 * args:    insopt_t *opt   I  ins options
 *          insstate_t *ins I  ins states
 *          double *Qn      O  system noise matrix
 * return: none
 * -------------------------------------------------------------------------*/
extern void getQn(const insopt_t *opt,const insstate_t *ins,double *Qn)
{
    opt->exprn?getprn(ins,opt,ins->dt,Qn):getQ(opt,ins->dt,Qn);
}
/* updates phi,P,Q of ekf----------------------------------------------------*/
static void updstat(const insopt_t *opt,insstate_t *ins,const double dt,
                    const double *x0,const double *P0,double *phi,double *P,
                    double *x,double *Q)
{
    /* determine approximate system noise covariance matrix */
    opt->exprn?getprn(ins,opt,dt,Q):
               getQ(opt,dt,Q);

    /* determine transition matrix
     * using last epoch ins states (first-order approx) */
    opt->exphi?precPhi(opt,dt,ins->Cbe,ins->re,ins->omgb,ins->fb,phi):
               getPhi1(opt,dt,ins->Cbe,ins->re,ins->omgb,ins->fb,phi);

#if UPD_IN_EULER
    getphi_euler(opt,dt,ins->Cbe,ins->re,ins->omgb,ins->fb,phi);
#endif
    /* propagate state estimation error covariance */
    if (fabs(dt)>=MAXUPDTIMEINT) {
        getP0(opt,P);
    }
    else {
        propP(opt,Q,phi,P0,P);
    }
    /* propagate state estimates noting that
     * all states are zero due to close-loop correction */
    if (x) propx(opt,x0,x);

    /* predict info. */
    if (ins->P0) matcpy(ins->P0,P  ,ins->nx,ins->nx);
    if (ins->F ) matcpy(ins->F ,phi,ins->nx,ins->nx);
}
/* propagate ins states and its covariance matrix----------------------------
 * args  : insstate_t *ins  IO  ins states
 *         insopt_t *opt    I   ins options
 *         double dt        I   time difference between current and precious
 *         double *x        O   updates ins states
 *         double *P        O   upadtes ins states covariance matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void propinss(insstate_t *ins,const insopt_t *opt,double dt,
                     double *x,double *P)
{
    int nx=ins->nx; double *phi,*Q;

    Q=mat(nx,nx); phi=mat(nx,nx);

    updstat(opt,ins,dt,ins->x,ins->P,phi,P,x,Q);
    free(Q); free(phi);
}
/* swap gps position/velocity measurement----------------------------------*/
static void swap(gmea_t *p1,gmea_t *p2)
{
    gmea_t tmp; tmp=*p1; *p1=*p2; *p2=tmp;
}
/* save precious epoch gps position/velocity measurements--------------------
 * args  :  insstate_t *ins    IO  ins sates
 *          sol_t *sol         I   gps solution (NULL: no input)
 *          gnss_meas_t *gmea  I   gps measurement (NULL: no input)
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int savegmeas(insstate_t *ins,const sol_t *sol,const gmea_t *gmea)
{
    int i; gmea_t meas={0};

    trace(3,"savegmeas:\n");

    if (gmea) {
        if (addgmea(&ins->gmeas,gmea)<0) {
            trace(2,"add gps position/velocity measurement fail\n");
            return 0;
        }
    }
    if (sol) {
        matcpy(meas.pe,sol->rr+0,1,3);
        meas.t=sol->time;

        for (i=0;i<3;i++) meas.std[i]=SQRT(sol->qr[i  ]);
        for (i=3;i<6;i++) meas.std[i]=SQRT(sol->qv[i-3]);

        if (addgmea(&ins->gmeas,&meas)<0) {
            trace(2,"add gps position/velocity measurement fail\n");
            return 0;
        }
    }
    /* here we only save NPOS epochs of precious */
    if (ins->gmeas.n>NPOS) {
        for (i=0;i<NPOS;i++) {
            swap(&ins->gmeas.data[i],&ins->gmeas.data[i+1]);
        }
        /* reset number of gps position/velocity measurements */
        ins->gmeas.n=NPOS;
    }
#ifdef TRACE
    trace(3,"precious epochs position/velocity measurement:\n");
    for (i=0;i<ins->gmeas.n;i++) {
       
        trace(3,"%3d time=%s\n",i+1,
              time_str(ins->gmeas.data[i].t,4));
    }
#endif
    return 1;
}
/* re-check attitude---------------------------------------------------------
 * args   :  insstate_t *ins  IO  ins states
 *           imud_t    *imu   I   imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int rechkatt(insstate_t *ins,const imud_t *imu)
{
    int i,j;
    double dt,vel[3*NPOS]={0},llh[3];
    double C[9],yaw,vn[3],rpy[3]={0};
    double vb[3],pvb[3],Cbe[9];

    trace(3,"rechkatt:\n");

    /* check gps solution status */
    for (i=0;i<NPOS;i++) {
        if (ins->gmeas.data[i].stat<=SOLQ_DGPS) return 0;
    }
    /* recheck ins attitude if need */
    if (ins->gmeas.n==NPOS) {

        /* velocity for trajectory*/
        for (i=NPOS;i>=2;i--) {
            if ((dt=timediff(ins->gmeas.data[i-1].t,ins->gmeas.data[i-2].t))>3.0
                ||fabs(dt)<=1E-5) {
                continue;
            }
            for (j=0;j<3;j++) {
                vel[3*(NPOS-i)+j]=(ins->gmeas.data[i-1].pe[j]-ins->gmeas.data[i-2].pe[j])/dt;
            }
        }
        /* check velocity whether is straight driving  */
        if (!chksdri(vel,NPOS-1)) {
            trace(2,"no straight driving\n");
            return 0;
        }
        /* check velocity */
        if (norm(vel,3)>MAXVEL&&norm(imu->gyro,3)<MAXGYRO) {

            /* velocity convert to attitude */
            ecef2pos(ins->gmeas.data[NPOS-1].pe,llh);
            ned2xyz(llh,C);

            /* yaw */
            matmul("TN",3,1,3,1.0,C,vel,0.0,vn);
            yaw=NORMANG(vel2head(vn)*R2D);

            /* attitude for current ins states */
            getatt(ins,rpy);

            if (fabs(yaw-NORMANG(rpy[2]*R2D))<MAXANG) return 0;
#if 0
            rpy[2]=yaw*D2R; /* reset yaw */
#else
            rpy[2]=(NORMANG(rpy[2]*R2D)+yaw)/2.0*D2R;
#endif
            rpy2dcm(rpy,C);
            matt(C,3,3,ins->Cbn);

            ned2xyz(llh,C);
            matmul("NN",3,3,3,1.0,C,ins->Cbn,0.0,Cbe);
            matmul("TN",3,1,3,1.0,Cbe,ins->ve,0.0,vb);
            matmul("TN",3,1,3,1.0,ins->Cbe,ins->ve,0.0,pvb);

            /* check again */
            if (fabs(norm(vb,3)-norm(pvb,3))<MINVEL&&(fabs(vb[1])<fabs(pvb[1]))
                &&(fabs(vb[2])<fabs(pvb[2]))) {

                matcpy(ins->Cbe,Cbe,3,3);
                trace(3,"recheck attitude ok\n");
                return 1;
            }
        }
    }
    trace(3,"no recheck attitude\n");
    return 0;
}
/* gps position measurement convert to velocity------------------------------*/
static void gm2vel(const gmea_t *g1,const gmea_t *g2,double *v)
{
    if (norm(g1->ve,3)>0.0&&norm(g2->ve,3)>0.0) {
        matcpy(v,g1->ve,1,3); return;
    }
    v[0]=(g1->pe[0]-g2->pe[0])/timediff(g1->t,g2->t);
    v[1]=(g1->pe[1]-g2->pe[1])/timediff(g1->t,g2->t);
    v[2]=(g1->pe[2]-g2->pe[2])/timediff(g1->t,g2->t);
}
/* reboot ins states---------------------------------------------------------*/
static void rebootsta(prcopt_t *opt,insstate_t *ins)
{
    trace(3,"initinsrt:\n");
    freelc(ins);
    initlc (&opt->insopt      ,ins);
    initodo(&opt->insopt.odopt,ins);
}
/* check imu body velocity---------------------------------------------------*/
static int chkvb(const insstate_t *ins)
{
#if 1
    double vb[3];
    matmul("TN",3,1,3,1.0,ins->Cbe,ins->ve,0.0,vb);
    return fabs(vb[1])<MAXVEL&&fabs(vb[2])<MAXVEL;
#else
    return 0;
#endif
}
/* check gps measurement data valid------------------------------------------*/
static int chkgmea(const gmea_t *data)
{
    if (data->stat==SOLQ_NONE) return 0;
    if (norm(data->pe,3)==0.0&&norm(data->ve,3)==0.0) return 0;
    if (norm(data->std,6)==0.0) return 0;
    if (data->t.time==0) return 0; return 1;
}
/* reboot ins-gnss loosely coupled if no coupling for a long time------------
 * args   :  insopt_t *opt    I   ins update options
 *           gmea_t *data     I   gps measurement data
 *           imud_t *imu      I   imu measurement data
 *           insstate_t *ins  IO  ins states
 * return : 1: reboot,0: no reboot
 * --------------------------------------------------------------------------*/
static int rebootlc(const insopt_t *opt,const gmea_t *data,const imud_t *imu,
                    insstate_t *ins)
{
    static gmea_t sols[MAXSOLS]={0};
    int i; double v[3];

    trace(3,"rebootlc: time=%s\n",time_str(imu->time,4));

    if (data==NULL) return 0;
    if (chkvb(ins)) return 0;

    /* check gps measurement data valid */
    if (!chkgmea(data)) {
        trace(2,"invalid gnss measurement data\n");
        return 1;
    }
    /* save gps measurement data to buffer */
    for (i=0;i<MAXSOLS-1;i++) sols[i]=sols[i+1]; sols[i]=*data;

    /* check solution status */
    for (i=0;i<MAXSOLS;i++) {
        if (sols[i].stat>opt->iisu) {
            trace(2,"gnss measurement status check fail\n");
            return 1;
        }
    }
    /* check velocity solution ok? */
    gm2vel(&sols[MAXSOLS-1],&sols[MAXSOLS-2],v);
    if (norm(v,3)<MINVEL) {
        trace(2,"velocity is not enough to reboot\n");
        return 1;
    }
    if (norm(imu->gyro,3)>MAXGYRO) {
        trace(2,"reboot fail due to large rotation\n");
        return 1;
    }
#if 1
    /* check measurement timestamp */
    for (i=0;i<MAXSOLS-1;i++) {
        if (timediff(sols[i+1].t,sols[i].t)>MAXDIFF||fabs(timediff(sols[i+1].t,sols[i].t))<1E-5) {
            trace(2,"measurement time check fail \n");
            return 1;
        }
    }
#endif
    /* reboot ins states */
    rebootsta((prcopt_t*)opt->gopt,ins);
    if (!ant2inins(sols[MAXSOLS-1].t,sols[MAXSOLS-1].pe,v,opt,
                   NULL,ins,NULL)) {
        trace(2,"reboot ins loosely coupled fail\n");
        return 1;
    }
    ins->time=sols[MAXSOLS-1].t;
    return 2;
}
/* ins-gnss couple function -------------------------------------------------
 * args  : insopt_t *opt      I  ins-gnss loosely/tightly coupled options
 *         imud_t *data       I  measurements from imu (corrected)
 *         insstate_t *ins    IO ins states from measurements updates
 *         gmea_t *gnss       I  measurements from gnss positioning
 *         int upd            I  updates flag (INSUPD_???)
 * return: 1:ok,0:failed
 * note  : if gnss.meas[i] is 0.0,then means no measurements to update
 *         if gnss.std [i] is 0.0,then means no measurements to update
 *         here ins and gnss time have already been synchronous
 *         if gnss==NULL ,then no updates ins states by gnss measurements
 * --------------------------------------------------------------------------*/
extern int lcigpos(const insopt_t *opt, const imud_t *data, insstate_t *ins,
                   gmea_t *gnss, int upd)
{
    double *phi,*Q,*P,*x,*pcov;
    double meas[NM],std[NM],cov[NM*NM]={0};
    int stat=1,i,nx=ins->nx,flag;

    trace(3,"lcigpos: upd=%d,time=%s\n",upd,time_str(data->time,4));

    /* backup current ins state for RTS */
    bckupinsinfo(ins,opt,2);

#if CHKNUMERIC
    /* check numeric of estimate state */
    for (i=0;i<3;i++) {
        if (isnan(ins->re[i])||isnan(ins->ve[i])||isnan(ins->ae[i])||
            isinf(ins->re[i])||isinf(ins->ve[i])||isinf(ins->ae[i])) {
            fprintf(stderr,"check numeric error: nan or inf\n");
            return 0;
        }
    }
#endif
    P=mat(nx,nx); x  =mat(nx, 1);
    Q=mat(nx,nx); phi=mat(nx,nx);

    updstat(opt,ins,ins->dt,ins->x,ins->P,phi,P,x,Q);
    matcpy(ins->x,x,nx, 1);
    matcpy(ins->P,P,nx,nx);

    /* ins mechanization update */
    ins->stat=INSS_NONE;
    if ((opt->soltype==0||opt->soltype==3)?
#if UPD_INS_E
        /* update ins states based on llh position mechanization */
        !updateinsn(opt,ins,data):
#else
        /* update ins states in e-frame */
        !updateins(opt,ins,data):
#endif
        !updateinsb(opt,ins,data)) {
        trace(2,"ins mechanization update fail\n");
        return 0;
    }
    /* backup predict ins state for RTS */
    bckupinsinfo(ins,opt,1);

    /* ins mechanization */
    if (upd==INSUPD_INSS) return stat;
    
#if REBOOT
    /* reboot ins loosely coupled if need */
    if ((flag=rebootlc(opt,upd==INSUPD_TIME?NULL:gnss,data,ins))) {
        if (flag==1) {
            trace(2,"ins/gnss loosely coupled still reboot\n");
            stat=0;
            goto exit;
        }
        ins->stat=INSS_REBOOT;
        stat=1;
        trace(3,"ins/gnss loosely coupled reboot ok\n");
        goto exit;
    }
#endif
    if (upd==INSUPD_TIME) { /* only ins mechanization */
        ins->stat=INSS_TIME;
    }
    else if (upd==INSUPD_MEAS&&gnss!=NULL) {

        trace(3,"gnss measurement data:\n");
        trace(3,"position=(%10.4lf %10.4lf %10.4lf)-(%6.4lf %6.4lf %6.4lf)\n",
              gnss->pe[0],gnss->pe[1],gnss->pe[2],
              gnss->std[0],
              gnss->std[1],
              gnss->std[2]);
        trace(3,"velocity=(%10.4lf %10.4lf %10.4lf)-(%6.4lf %6.4lf %6.4lf)\n",
              gnss->ve[0],gnss->ve[1],gnss->ve[2],
              gnss->std[3],
              gnss->std[4],
              gnss->std[5]);

        /* position and velocity for antenna */
        for (i=0;i<3;i++) {
            meas[i+0]=gnss->pe[i]; meas[i+3]=gnss->ve[i];
        }
#if USE_MEAS_COV
        asi_blk_mat(cov,NM,NM,gnss->covp,3,3,0,0);
        asi_blk_mat(cov,NM,NM,gnss->covv,3,3,3,3);
        pcov=cov;
#else
        matcpy(std,gnss->std,1,NM);
        pcov=NULL;
#endif
        /* determine estimated and covariance matrix */
        if (opt->updint==UPDINT_GNSS) {
            updstat(opt,ins,timediff(gnss->t,ins->plct),ins->xa,ins->Pa,phi,P,x,Q);
        }
        /* ins-gnss loosely coupled */
        if ((stat=lcfilt(opt,ins,meas,std,pcov,timediff(gnss->t,ins->time),x,P))) {

            /* propagate by gps epoch time internal */
            if (opt->updint==UPDINT_GNSS) {
                matcpy(ins->xa,x,nx, 1);
                matcpy(ins->Pa,P,nx,nx);
            }
            else {
                matcpy(ins->x,x,nx,1);
                matcpy(ins->P,P,nx,nx);
            }
            ins->stat=INSS_LCUD;
        }
        ins->plct=ins->time;
        if (stat) {
            savegmeas(ins,NULL,gnss);
            rechkatt(ins,data); /* recheck attitude */
        }
    }
    else {
        trace(2,"no gnss measurement data\n");
    }
exit:
    free(x); free(P);
    free(Q); free(phi);
    return stat;
}
/* convert ins solution status to sol_t struct-------------------------------
 * args   :  insstate_t *ins  IO ins solution status
 *           insopt_t *opt    I  ins options
 *           sol_t *sol       O  solution status
 * return : none
 * --------------------------------------------------------------------------*/
extern void ins2sol(insstate_t *ins,const insopt_t *opt,sol_t *sol)
{
    double Pa[9],Pp[9],Pv[9];
    int i;

    trace(3,"ins2sol: time=%s\n",time_str(ins->time,4));

    sol->ista=(unsigned char)ins->stat;
    sol->stat=(unsigned char)ins->gstat;
    sol->ns  =(unsigned char)ins->ns;

    /* receiver position/velocity */
    matcpy(sol->rr  ,ins->re,1,3);
    matcpy(sol->rr+3,ins->ve,1,3);

    /* vehicle attitude */
    getatt(ins,sol->att);
    insstatocov(opt,ins,Pa,Pv,Pp);

    for (i=0;i<3;i++) sol->qa[i]=(float)Pa[i+i*3];
    sol->qa[3]=(float)Pa[3]; /* xy */
    sol->qa[4]=(float)Pa[7]; /* zy */
    sol->qa[5]=(float)Pa[2]; /* zx */

    for (i=0;i<3;i++) sol->qr[i]=(float)Pp[i+i*3];
    sol->qr[3]=(float)Pp[3]; /* xy */
    sol->qr[4]=(float)Pp[7]; /* zy */
    sol->qr[5]=(float)Pp[2]; /* zx */

    for (i=0;i<3;i++) sol->qv[i]=(float)Pv[i+i*3];
    sol->qv[3]=(float)Pv[3]; /* xy */
    sol->qv[4]=(float)Pv[7]; /* zy */
    sol->qv[5]=(float)Pv[2]; /* zx */

    /*  velocity/acceleration in body frame */
    matmul3v("T",ins->Cbe,ins->ve,sol->vb);
    matmul3v("T",ins->Cbe,ins->ae,sol->ab);

    /* imu gyro/accl bias */
    matcpy(sol->bg,ins->bg,1,3);
    matcpy(sol->ba,ins->ba,1,3);

    /* Ma and Mg */
    matcpy(sol->Ma,ins->Ma,3,3);
    matcpy(sol->Mg,ins->Mg,3,3);

    /* imu raw data */
    matcpy(sol->imu.gyro,ins->omgb,1,3);
    matcpy(sol->imu.accl,ins->fb,1,3);

    sol->dtr[0]=ins->dtr[0]; /* receiver clock bias (s) */
    sol->dtr[1]=ins->dtr[1]; /* glo-gps time offset (s) */
    sol->dtr[2]=ins->dtr[2]; /* gal-gps time offset (s) */
    sol->dtr[3]=ins->dtr[3]; /* bds-gps time offset (s) */

    sol->dtrr=ins->dtrr;
    sol->time=ins->time;
    sol->vo=ins->vo; 
}
/* given gps antenna position/velocity to initial ins states-----------------
 * args    :  gtime_t time     I  time of antenna position/velocity
 *            double *rr       I  gps antenna position (ecef)
 *            double *vr       I  gps antenna velocity (ecef)
 *            insopt_t *opt    I  ins options
 *            imu_t *imu       I  imu measurement data
 *            insstate_t *ins  O  initialed ins states
 *            int *iimu        O  index of align imu measurement data
 * return  : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int ant2inins(gtime_t time,const double *rr,const double *vr,
                     const insopt_t *opt,const imu_t *imu,insstate_t *ins,
                     int *iimu)
{
    double llh[3],vn[3],C[9],rpy[3]={0};
    imud_t imus={0};
    int i;

    trace(3,"ant2inins: time=%s\n",time_str(time,4));

    /* velocity in navigation frame */
    ecef2pos(rr,llh);
    ned2xyz(llh,C);
    matmul("TN",3,1,3,1.0,C,vr,0.0,vn);

    matcpy(ins->rn,llh,1,3);
    matcpy(ins->vn,vn ,1,3);

    rpy[2]=vel2head(vn); /* yaw */

    rpy2dcm(rpy,C);
    matt(C,3,3,ins->Cbn);

    ned2xyz(llh,C);
    matmul("NN",3,3,3,1.0,C,ins->Cbn,0.0,ins->Cbe);

    /* find closest imu measurement index */
    if (imu) {
        for (i=0;i<imu->n;i++)  {
            if (fabs(timediff(time,imu->data[i].time))<DTTOL) break;
        }
        if (i>=imu->n) {
            trace(2,"no common time of imu and observation\n");
            return 0;
        }
        /* align imu measurement data */
        imus=imu->data[i];
        if (norm(imus.gyro,3)>MAXROT) return 0;

        /* index of align imu data */
        if (iimu) *iimu=i;
    }
    /* initial ins position */
    gapv2ipv(rr,vr,ins->Cbe,opt->lever,
             &imus,
             ins->re,ins->ve);
    return 1;
}
