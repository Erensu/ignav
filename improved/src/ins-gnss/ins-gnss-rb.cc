/*--------------------------------------------------------------------------------
 * ins-gnss-rb.cc : ins-gnss loosely coupled with robocentric formulation
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
 *----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants/macros ----------------------------------------------------------*/
#define NNAC            3                  /* number of accl process noise */
#define NNGY            3                  /* number of gyro process noise */
#define NNPX(opt)       (NNAC+NNGY)
#define INAC            0                  /* index of accl process noise */
#define INGY            NNAC               /* index of gyro process noise */
#define CORRETIME       360.0              /* correlation time for gauss-markov process */
#define MAXUPDTIMEINT   60.0               /* max update time internal is same as correlation time */
#define MAXSYNDIFF      1.0                /* max time difference of ins and gnss time synchronization */
#define MAXVARDIS       10.0               /* max variance of disable estimated state */
#define MAXINOP         1000.0             /* max innovations for updates ins states */
#define MAXINOV         100.0              /* max innovations for updates ins states */
#define MAXVEL          5.0                /* max velocity for refine attitude */
#define MINVEL          3.0                /* min velocity for refine attitude */
#define MAXGYRO         (5.0*D2R)          /* max gyro measurement for refine attitude */

#define NMP             3                  /* number of position measurements for ins-gnss coupled */
#define NMV             3                  /* number of velocity measurements for ins-gnss coupled */
#define NM              (NMP+NMV)          /* number of all measurements for ins-gnss coupled */
#define IMP             0                  /* index of position measurements */
#define IMV             NMP                /* index of velocity measurements */

#define MAXSOLS         5                  /* max number of solutions for reboot lc  */
#define MAXDIFF         10.0               /* max time difference between solution */

/* propagate matrix for stochastic parameters---------------------------------*/
static void stochasticPhi(int opt,int ix,int nix,int nx,double dt,double *phi)
{
    int i; if (nix<=0) return;
    for (i=ix;i<ix+nix;i++) {
        if (opt==INS_RANDOM_CONS ) phi[i+i*nx]=1.0;
        if (opt==INS_RANDOM_WALK ) phi[i+i*nx]=1.0;
        if (opt==INS_GAUSS_MARKOV) phi[i+i*nx]=exp(-fabs(dt)/CORRETIME);
    }
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
/* the error state system dynamics matrix-------------------------------------*/
static void stateF(const insopt_t *opt,const double *Cbe,const double *re,const double *ve,
                   const double *omgb,const double *fib,
                   double *F)
{
    int i,j,nx=xnX(opt),iax,nax,ibg,iba,isg,irg,nbg,nba,nsg,nrg,ivx,nvx,ipx,npx,isa,nsa,ira,nra;
    double vb[3],ge[3],rn[3],res;
    double W[18]={0},WC[18]={0},Wb[9],Wg[9],W2[9],I[9]={1,0,0,0,1,0,0,0,1};
    double wieb[3],we[3]={0,0,OMGE};
    double T1[9],T2[9],T3[9],t[3],T[9]={0},fs[9]={0},omgs[9]={0};
    double Ts[9],Tr[18];

    setzero(F,nx,nx);
    iax=xiA(opt); ibg=xiBg(opt); iba=xiBa(opt); isg=xiSg(opt); irg=xiRg(opt);
    nax=xnA(opt); nbg=xnBg(opt); nba=xnBa(opt); nsg=xnSg(opt); nrg=xnRa(opt);
    ivx=xiV(opt); isa=xiSa(opt);
    nvx=xnV(opt); nsa=xnSa(opt);
    ipx=xiP(opt); ira=xiRa(opt);
    npx=xnP(opt); nra=xnRa(opt);
    matmul("NN",3,3,3,1.0,Omge,Omge,0.0,W2);
    T[0]=omgb[0];
    T[4]=omgb[1];
    T[8]=omgb[2];
    matmul("NN",3,3,3,1.0,Cbe,T,0.0,omgs);

    W[0 ]=omgb[1]; W[3 ]=omgb[2];
    W[7 ]=omgb[0]; W[10]=omgb[2];
    W[14]=omgb[0]; W[17]=omgb[1];
    for (i=iax;i<iax+nax;i++) {
        for (j=iax;j<iax+nax;j++) F[i+j*nx]=-Omge[i-iax+(j-iax)*3];
        for (j=ibg;j<ibg+nbg;j++) F[i+j*nx]= Cbe[i-iax+(j-ibg)*3];
        for (j=isg;j<isg+nsg;j++) F[i+j*nx]=omgs[i-iax+(j-isg)*3];
        for (j=irg;j<irg+nrg;j++) F[i+j*nx]=  WC[i-iax+(j-irg)*3];
    }
    matmul("TN",3,1,3,1.0,Cbe,ve,0.0,vb);
    skewsym3(vb,Wb);
    pregrav(re,ge);
    skewsym3(ge,Wg);
    matmul33("NTN",Wb,Cbe,Omge,3,3,3,3,T1);
    matmul("TN",3,3,3,1.0,Cbe,Wg,0.0,T2);

    matmul("NN",3,1,3,1.0,W2,re,0.0,t);
    skewsym3(t,Wb);
    matmul("TN",3,3,3,1.0,Cbe,Wb,0.0,T3);

    for (i=0;i<9;i++) T[i]=-3.0*T1[i]+T2[i]-T3[i];
    for (i=ivx;i<ivx+nvx;i++) {
        for (j=iax;j<iax+nax;j++) F[i+j*nx]=T[i-ivx+(j-iax)*3];
    }
    matmul("TN",3,1,3,1.0,Cbe,we,0.0,wieb);
    for (i=0;i<3;i++) t[i]=omgb[i]-3.0*wieb[i];
    skewsym3(t,Wg);
    for (i=ivx;i<ivx+nvx;i++) {
        for (j=ivx;j<ivx+nvx;j++) F[i+j*nx]=-Wg[i-ivx+(j-ivx)*3];
    }
    ecef2pos(re,rn);
    res=georadi(rn);
    matmul("NT",3,3,1,-2.0/(res*norm(re,3)),ge,re,0.0,T1);
    matmul("TN",3,3,3,1.0,Cbe,T1,0.0,T2);
    matmul("TN",3,3,3,1.0,Cbe,W2,0.0,T3);
    skewsym3(vb,Wb);

    fs[0]=fib[0];
    fs[4]=fib[1];
    fs[8]=fib[2];
    matmul("NN",3,3,3,1.0,Wb,omgs,0.0,Ts);
    matmul("NN",3,6,3,1.0,Wb,W,0.0,Tr);

    for (i=ivx;i<ivx+nvx;i++) {
        for (j=ipx;j<ipx+npx;j++) F[i+j*nx]=T2[i-ivx+(j-ipx)*3]-T3[i-ivx+(j-ipx)*3];
        for (j=ibg;j<ibg+nbg;j++) F[i+j*nx]=Wb[i-ivx+(j-ibg)*3];
        for (j=iba;j<iba+nba;j++) F[i+j*nx]= I[i-ivx+(j-iba)*3];
        for (j=isa;j<isa+nsa;j++) F[i+j*nx]=fs[i-ivx+(j-isa)*3];
        for (j=ira;j<ira+nra;j++) F[i+j*nx]= W[i-ivx+(j-ira)*3];
        for (j=isg;j<isg+nsg;j++) F[i+j*nx]=Ts[i-ivx+(j-isg)*3];
        for (j=irg;j<irg+nrg;j++) F[i+j*nx]=Tr[i-ivx+(j-irg)*3];
    }
    skewsym3(ve,Wb);
    for (i=ipx;i<ipx+npx;i++) {
        for (j=iax;j<iax+nax;j++) F[i+j*nx]=-Wb[i-ipx+(j-iax)*3];
        for (j=ivx;j<ivx+nvx;j++) F[i+j*nx]=Cbe[i-ipx+(j-ivx)*3];
    }
    stochasticF(opt->baproopt,iba,nba,nx,F);
    stochasticF(opt->bgproopt,ibg,nbg,nx,F);
    stochasticF(opt->sgproopt,isg,nsg,nx,F);
    stochasticF(opt->saproopt,isa,nsa,nx,F);
    stochasticF(opt->raproopt,ira,nra,nx,F);
    stochasticF(opt->rgproopt,irg,nrg,nx,F);
}
/* determine transition matrix(first-order approx: Phi=I+F*dt)---------------*/
static void statePhi(const insopt_t *opt,const double dt,const double *Cbe,const double *re,
                     const double *ve,
                     const double *omgb,const double *fib,
                     double *Phi)
{
    int i,nx=xnX(opt);
    double *F=zeros(nx,nx),*FF=zeros(nx,nx);
    double *FFF=zeros(nx,nx),*I=eye(nx);

    stateF(opt,Cbe,re,ve,omgb,fib,F);
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
    free(F); free(FF); free(FFF);
    free(I);
}
/* process noise gain matrix-------------------------------------------------*/
static void getGn(const insopt_t *opt,const double *Cbe,const double *ve,double *Gn)
{
    int i,j,nx=xnX(opt),nprn=NNPX(opt),ivx,nvx,iax,nax;
    double I[9]={1,0,0,0,1,0,0,0,1},Wb[9],vb[3];

    iax=xiA(opt); nax=xnA(opt);
    ivx=xiV(opt); nvx=xnV(opt);

    setzero(Gn,nx,nprn);
    matmul("TN",3,1,3,1.0,Cbe,ve,0.0,vb);
    skewsym3(vb,Wb);

    for (i=iax;i<iax+nax;i++) {
        for (j=INAC;j<INAC+NNAC;j++) Gn[i+j*nx]=Cbe[i-iax+(j-INAC)*3];
    }
    for (i=ivx;i<ivx+nvx;i++) {
        for (j=INAC;j<INAC+NNAC;j++) Gn[i+j*nx]= I[i-ivx+(j-INAC)*3];
        for (j=INGY;j<INGY+NNGY;j++) Gn[i+j*nx]=Wb[i-ivx+(j-INGY)*3];
    }
    trace(3,"Gn=\n");
    tracemat(3,Gn,nx,nprn,12,5);
}
/* process noise covariance matrix-------------------------------------------*/
static void stateQ(const insstate_t *ins,const insopt_t *opt,double dt,double* Q)
{
    int nprn=NNPX(opt),nx=xnX(opt),i;
    double *Qn=zeros(nprn,nprn),*Gn=mat(nprn,nx);

    for (i=INAC;i<INGY+NNAC;i++) Qn[i+i*nprn]=opt->psd.gyro*fabs(dt);
    for (i=INGY;i<INGY+NNGY;i++) Qn[i+i*nprn]=opt->psd.accl*fabs(dt);

    getGn(opt,ins->Cbe,ins->ve,Gn);
    matmul33("NNT",Gn,Qn,Gn,nx,nprn,nprn,nx,Q);
    free(Qn); free(Gn);
}
/* propagate state estimation error covariance-------------------------------*/
static void propagateP(const insopt_t *opt,const double *Q,const double *phi,
                       const double *P0,double *P)
{
    int i,j,nx=xnX(opt),irc,nrc;
    double *PQ=mat(nx,nx),*Phi2=mat(nx,nx);

    irc=xiRc(opt); nrc=xnRc(opt);
    
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) PQ[i+j*nx]=P0[i+j*nx]+0.5*Q[i+j*nx];
    }
    matmul33("NNT",phi,PQ,phi,nx,nx,nx,nx,Phi2);
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) P[i+j*nx]=Phi2[i+j*nx]+0.5*Q[i+j*nx];
    }
    /* initialize every epoch for clock (white noise) */
    initP(irc,nrc,nx,opt->unc.rc,30.0,P);
    free(PQ); free(Phi2);
}
/* propagate state estimates noting that all states are zero due to closed-loop
 * correction----------------------------------------------------------------*/
static void propagatex(const insopt_t *opt,const double *x0,double *x)
{
    int i; for (i=0;i<xnCl(opt);i++) x[i]=1E-20;
}
/* updates phi,P,Q of ekf----------------------------------------------------*/
static void updstat(const insopt_t *opt,insstate_t *ins,const double dt,const double *x0,
                    const double *P0,double *phi,
                    double *P,double *x,double *Q)
{
    trace(3,"updstat:\n");

    /* system transmit and noise matrix */
    statePhi(opt,dt,ins->Cbe,ins->re,ins->ve,ins->omgb,ins->fb,phi);
    stateQ(ins,opt,dt,Q);

    trace(3,"Phi=\n");
    tracemat(3,phi,ins->nx,ins->nx,12,5);
    
    /* propagate error state covariance */
    if (fabs(dt)>=MAXUPDTIMEINT) getP0(opt,P);
    else {
        propagateP(opt,Q,phi,P0,P);
    }
    /* update error states */
    if (x) propagatex(opt,x0,x);

    /* backup predict information */
    if (ins->P0) matcpy(ins->P0,P  ,ins->nx,ins->nx);
    if (ins->F ) matcpy(ins->F ,phi,ins->nx,ins->nx);

    trace(3,"P=\n");
    tracemat(3,P,ins->nx,ins->nx,12,5);

    trace(3,"Q=\n");
    tracemat(3,Q,ins->nx,ins->nx,12,5);
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
/* measurement sensitive-matrix----------------------------------------------*/
static int bldHvR(const insopt_t *opt,const double *pos,const double *Cbe,
                  const double *lever,const double *omgb,const double *fib,
                  const double *meas,const double *re,const double *ve,
                  const double *ae,const double *std,const double *cov,const double *P,
                  double *H,double *v,double *R)
{
    int i,j,nx=xnX(opt),nm=0,ind[NM];
    int iax,nax,ibg,isg,irg,nbg,nsg,nrg,ivx,nvx,ipx,npx,idt,ndt,ila,nla;
    double r1[9],v1[9],v5[9],I[9]={1,0,0,0,1,0,0,0,1},r[NM]={0};
    double zr[3],zv[3];
    double dt1[3],dt2[3],ds[9],drg[18],dla[9];
    double We[9],T[9];

    /* remove lever-arm effects */
    rmlever(NULL,re,ve,lever,Cbe,omgb,zr,zv);
    isg=xiSg(opt); irg=xiRg(opt); ibg=xiBg(opt);
    nsg=xnSg(opt); nrg=xnRa(opt); nbg=xnBg(opt);

    iax=xiA(opt); nax=xnA(opt);
    ivx=xiV(opt); nvx=xnV(opt);
    ipx=xiP(opt); npx=xnP(opt);
    idt=xiDt(opt); ndt=xnDt(opt);
    ila=xiLa(opt); nla=xnLa(opt);

    /* for position measurement */
    jacobian_p_att(Cbe,lever,r1);
    jacobian_p_dt (omgb,lever,Cbe,ve,dt1);

    for (i=IMP;i<IMP+NMP;i++) {
        if (meas[i]!=0.0) {
            if (fabs(v[nm]=(meas[i]-zr[i-IMP]))>=MAXINOP) {
                trace(2,"too large innovations for position\n");
            }
            if (H) {
                for (j=iax;j<iax+nax;j++) H[j+nm*nx]= r1[i-IMP+(j-iax)*3];
                for (j=ipx;j<ipx+npx;j++) H[j+nm*nx]=-I [i-IMP+(j-ipx)*3];
                for (j=idt;j<idt+ndt;j++) H[j+nm*nx]= dt1[i-IMP];
                for (j=ila;j<ila+nla;j++) H[j+nm*nx]=-Cbe[i-IMP+(j-ila)*3];
            }
            r[nm]=std[i]==0.0?SQR(30.0):SQR(std[i]);
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
    skewsym3(ve,We);
    
    for (i=IMV;i<IMV+NMV;i++) {
        if (meas[i]!=0.0) {
            if (fabs(v[nm]=(meas[i]-zv[i-IMV]))>=MAXINOV) {
                trace(2,"too large innovations for velocity\n");
            }
            if (H) {
                for (j=iax;j<iax+nax;j++) H[j+nm*nx]= v1 [i-IMV+(j-iax)*3]+We[i-IMV+(j-iax)*3];
                for (j=ivx;j<ivx+nvx;j++) H[j+nm*nx]=-Cbe[i-IMV+(j-ivx)*3];
                for (j=ibg;j<ibg+nbg;j++) H[j+nm*nx]= v5 [i-IMV+(j-ibg)*3];
                for (j=isg;j<isg+nsg;j++) H[j+nm*nx]= ds [i-IMV+(j-isg)*3];
                for (j=irg;j<irg+nrg;j++) H[j+nm*nx]= drg[i-IMV+(j-irg)*3];
                for (j=idt;j<idt+ndt;j++) H[j+nm*nx]= dt2[i-IMV];
                for (j=ila;j<ila+nla;j++) H[j+nm*nx]= dla[i-IMV+(j-ila)*3];
            }
            r[nm]=std[i]==0.0?SQR(30.0):SQR(std[i]);
            ind[nm++]=i;
        }
    }
    if (nm&&R) {
        if (cov) {
            for (i=0;i<nm;i++) {
                for (j=0;j<nm;j++) {
                    R[i+nm*j]=cov[ind[i]+NM*ind[j]];
                }
            }
        }
        else {
            for (i=0;i<nm;i++) R[i+i*nm]=r[i];
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
    return nm;
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
/* validation of ins-gnss coupled solution-----------------------------------*/
static int valsol(const insopt_t *opt,double *x,double *P,const double *R,const double *v,
                  int nv,double thres)
{
    int iba,ibg,nba,nbg,i=0,j;
    double fact=thres*thres;

    ibg=xiBg(opt); nbg=xnBg(opt);
    iba=xiBa(opt); nba=xnBa(opt);

    /* check estimated states */
    if (     x[  0]==DISFLAG&&norm(x+  0,3)>15.0*D2R) i|=1;
    if (nba&&x[iba]==DISFLAG&&norm(x+iba,3)>1E5*Mg2M) i|=1;
    if (nbg&&x[ibg]==DISFLAG&&norm(x+ibg,3)>15.0*D2R) i|=1;
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
/* close-loop correction-----------------------------------------------------*/
static void lcclprb(double *x,double *Cbe,double *re,double *vb,double *fib,double *omgb,
                    double *Gg, double *rec,double *vbc,double *aec,
                    double *bac, double *bgc,double *Mac,double *Mgc,
                    double *leverc,double *Cbec,double *fibc,
                    double *omgbc,double *vec,const insopt_t *opt)
{
    double T[9],I[9]={1,0,0,0,1,0,0,0,1},ang[3];
    int i,ibg,iba,nbg,nba,ipx,ivx;

    ibg=xiBg(opt); nbg=xnBg(opt); iba=xiBa(opt); nba=xnBa(opt);
    ivx=xiV (opt); ipx=xiP (opt);

#if COR_IN_ROV
    /* attitude correction */
    if (x[0]!=DISFLAG) {
        skewsym3(x,T); for (i=0;i<9;i++) I[i]-=T[i];
        matmul3("NN",I,Cbe,Cbec);
    }
#else
    /* attitude correction */
    if (x[0]!=DISFLAG) correctatt(x,Cbe,Cbec);
#endif
    /* accl and gyro bias correction */
    if (nba&&x[iba]!=DISFLAG) for (i=0;i<3;i++) bac[i]+=x[iba+i];
    if (nbg&&x[ibg]!=DISFLAG) for (i=0;i<3;i++) bgc[i]+=x[ibg+i];

    /* velocity and position correction */
    for (i=0;i<3;i++) {
        rec[i]=re[i]-x[ipx+i];
        vbc[i]=vb[i]-x[ivx+i];
    }
    matmul("NN",3,1,3,1.0,Cbec,vbc,0.0,vec);
    if (fib&&omgb&&fibc&&omgbc&&Gg) {

        /* correction imu accl and gyro measurements */
        ins_errmodel2(fib,omgb,Mac,Mgc,bac,bgc,Gg,fibc,omgbc);
        getaccl(fibc,Cbec,rec,vec,aec);
    }
}
/* position filter-----------------------------------------------------------*/
static int lcfilt(const insopt_t *opt, insstate_t *ins, const double *meas,
                  const double *std, const double *cov,const double dt,
                  double *x, double *P)
{
    int nm,info=0,stat,nx=ins->nx,i;
    int iax,iba,ibg,nax,nba,nbg;
    double re[3],ve[3],vb[3],ae[3],Cbe[9],fib[3],omgb[3],stde[6]={0};
    double *H,*v,*R;

    /* close-loop correction states */
    double rec[3],vec[3],aec[3],Cbec[9],fibc[3],omgbc[3],vbc[3];
    double bac[3],bgc[3],Mac[9],Mgc[9],Gg[9],leverc[3];

    ibg=xiBg(opt); nbg=xnBg(opt);
    iba=xiBa(opt); nba=xnBa(opt);
    iax=xiA (opt); nax=xnA (opt);

    trace(3,"lcfiltr: dt=%.3lf\n",dt);

    /* check gnss and ins time synchronization */
    if (fabs(ins->age=dt)>=MAXSYNDIFF) {
        trace(2,"gnss and ins time synchronization fail,dt=%5.3lf\n",dt);
        return 0;
    }
    /* interpolate ins states: attitude,position and velocity */
    interppv(ins,dt,re,ve,ae,Cbe);
    prepara(ins,fib,omgb,Mgc,Mac,Gg,bac,bgc,leverc);
    matmul("TN",3,1,3,1.0,Cbe,ve,0.0,vb);

    /* build H,v and R matrix from input measurements */
    H=zeros(NM,nx); v=zeros(NM,1);
    R=zeros(NM,NM);

    if ((nm=bldHvR(opt,NULL,Cbe,leverc,omgb,fib,meas,re,ve,ae,std,cov,P,H,v,R))) {
        if (cov) {
            for (i=0;i<6;i++) stde[i]=cov[i+i*6];
        }
        else {
            matcpy(stde,std,1,6);
        }
        /* disable gyro/accl/att state if measurement is bad */
        if (norm(stde,6)>=MAXVARDIS) {
            for (i=iax;i<iax+nax;i++) unusex(opt,i,ins,x);
            for (i=iba;i<iba+nba;i++) unusex(opt,i,ins,x);
            for (i=ibg;i<ibg+nbg;i++) unusex(opt,i,ins,x);
        }
        /* ekf filter */
        if ((info=filter(x,P,H,v,R,nx,nm))) {
            trace(2,"filter error (info=%d)\n",info);
            free(H); free(v); free(R);
            return 0;
        }
    }
    else {
        trace(3,"no gnss position and velocity measurement\n");
        free(H); free(v); free(R);
        return 0;
    }
    /* close-loop for estimated states */
    lcclprb(x,Cbe,re,vb,fib,omgb,Gg,rec,vbc,aec,bac,bgc,Mac,Mgc,leverc,Cbec,fibc,omgbc,vec,opt);

    /* post-fit residuals for ins-gnss coupled */
    if ((nm=bldHvR(opt,NULL,Cbec,leverc,omgbc,fibc,meas,rec,vec,aec,std,cov,P,H,v,R))) {

        /* validation of solutions */
        if ((stat=valsol(opt,x,P,R,v,nm,10.0))) {

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
/* ins-gnss couple function ---------------------------------------------------
 * args  : insopt_t *opt      I  ins-gnss loosely/tightly coupled options
 *         imud_t *data       I  measurements from imu (corrected)
 *         insstate_t *ins    IO ins states from measurements updates
 *         gmea_t *gnss       I  measurements from gnss positioning
 *         int upd            I  updates flag (INSUPD_???)
 * return: 1:ok,0:failed
 * --------------------------------------------------------------------------*/
extern int lcigposrb(const insopt_t *opt, const imud_t *data, insstate_t *ins,gmea_t *gnss,
                     int upd)
{
    double *Phi,*Q,*P,*x,*pcov,meas[NM],std[NM],cov[NM*NM]={0};
    int stat=1,i,nx=ins->nx,flag;

    trace(3,"lcigpos: upd=%d,time=%s\n",upd,time_str(data->time,4));

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
    P=mat(nx,nx); Phi=mat(nx,nx);
    Q=mat(nx,nx); x=mat(nx,1);

    updstat(opt,ins,ins->dt,ins->x,ins->P,Phi,P,x,Q);
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
    if (upd==INSUPD_TIME) {
        ins->stat=INSS_TIME; /* only ins mechanization */
    }
    else if (upd==INSUPD_MEAS&&gnss) {
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
            meas[i+0]=gnss->pe[i];
            meas[i+3]=gnss->ve[i];
        }
#if USE_MEAS_COV
        asi_blk_mat(cov,NM,NM,gnss->covp,3,3,0,0);
        asi_blk_mat(cov,NM,NM,gnss->covv,3,3,3,3);
        pcov=cov;
#else
        matcpy(std,gnss->std,1,NM);
        pcov=NULL;
#endif
        /* ins-gnss loosely coupled */
        if ((stat=lcfilt(opt,ins,meas,std,pcov,timediff(gnss->t,ins->time),x,P))) {

            /* update ins states covariance */
            matcpy(ins->x,x,nx,1);
            matcpy(ins->P,P,nx,nx);

            trace(3,"P=\n");
            tracemat(3,ins->P,nx,nx,12,5);

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
    free(x); free(P); free(Phi);
    free(Q);
    return stat;
}

