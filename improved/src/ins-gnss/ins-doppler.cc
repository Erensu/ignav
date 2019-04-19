/*-----------------------------------------------------------------------------
* ins-doppler.cc : ins-gnss doppler coupled common functions
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
* history : 2017/11/28 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define MATLAB      1          /* MATLAB generate jacobians for velocity measurements by attitude */
#define MAXIONV     1E3        /* max innovations for doppler measurement */
#define STDOPP      0.1        /* doppler measurement standard variance */
#define UNC_CLKR    (10.0)     /* default initial receiver clock drift uncertainty (s/s) */

/* jacobian of doppler measurement by ins position---------------------------
 * args   :  double *rs     I  satellite position
 *           double *rr     I  receiver position
 *           double *vs     I  satellite velocity
 *           double *dopdp  O  output jacobian matrix
 * return : none
 * --------------------------------------------------------------------------*/
static void jaco_dop_dp(const double *rs,const double *rr,const double *vs,
                        const double *vr,double *dopdp)
{
    int i;
    double dv[3],dp[3],pr,t;

#if MATLAB
    for (i=0;i<3;i++) dp[i]=rs[i]-rr[i]; pr=norm(dp,3);
    for (i=0;i<3;i++) dv[i]=vs[i]-vr[i];

    t=pow(SQR(pr),1.5);

    dopdp[0]=2.0*dv[2]*dp[2]*dp[0]/(2.0*t)-dv[0]/pr+2.0*dv[0]*dp[0]*dp[0]/(2.0*t)+
             2.0*dv[1]*dp[1]*dp[0]/(2.0*t)+OMGE/CLIGHT*vs[1];
    dopdp[1]=2.0*dv[1]*dp[1]*dp[0]/(2.0*t)-dv[1]/pr+2.0*dv[2]*dp[2]*dp[1]/(2.0*t)+
             2.0*dv[0]*dp[0]*dp[1]/(2.0*t)-OMGE/CLIGHT*vs[0];
    dopdp[2]=2.0*dv[0]*dp[0]*dp[2]/(2.0*t)-dv[1]/pr+2.0*dv[1]*dp[1]*dp[2]/(2.0*t)+
             2.0*dv[2]*dp[2]*dp[2]/(2.0*t);
    for (i=0;i<3;i++) dopdp[i]=-dopdp[i];
#else
    for (i=0;i<3;i++) dopdp[i]=0.0;
#endif
}
/* jacobian of doppler measurement by attitude error-------------------------
 * args   :  double *dr    I  jacobian matrix of ins position
 *           double *Cbe   I  transform matrix from b-frame to ecef
 *           double *l     I  lever arm of b-frame to antenna position
 *           double *dopdv I  jacobian matrix of velocity
 *           double *gyro  I  gyro measurement data
 *           double *dopda O  output jacobian matrix
 * return : none
 * --------------------------------------------------------------------------*/
static void jaco_dop_da(const double *dr,const double *Cbe,const double *l,
                        const double *dopdv,const double *gyro,double *dopda)
{
    int i;
    double T[9],sl[3],S[9];

    /* partial matrix by ins position */
    matmul("NN",3,1,3,-1.0,Cbe,l,0.0,dopda);
    skewsym3(dopda,T);
    matmul("NN",1,3,3,1.0,dr,T,0.0,dopda);

    /* partial matrix by ins velocity */
    skewsym3(gyro,T);
    matmul33("NNN",Cbe,T,l,3,3,3,1,sl);
    skewsym3(sl,T);
    matmul("NN",1,3,3,1.0,dopdv,T,0.0,S);
    for (i=0;i<3;i++) dopda[i]-=S[i];

    matmul("NN",3,1,3,1.0,Cbe,l,0.0,sl);
    skewsym3(sl,T);
    matmul("NN",3,3,3,1.0,Omge,T,0.0,S);
    matmul("NN",1,3,3,1.0,dopdv,S,0.0,T);
    for (i=0;i<3;i++) dopda[i]+=T[i];
}
/* jacobian of doppler measurement by ins velocity---------------------------*/
static void jaco_dop_dv(const double *e,const double *rs,double *dopdv)
{
    dopdv[0]=e[0]-OMGE/CLIGHT*rs[1];
    dopdv[1]=e[1]+OMGE/CLIGHT*rs[0];
    dopdv[2]=e[2];
}
/* jacobian of doppler measurement by gyro bias------------------------------*/
static void jaco_dop_bg(const double *dopdv,const double *Cbe,const double *l,
                        double *dopbg)
{
    double T[9],S[9];

    skewsym3(l,T);
    matmul("NN",3,3,3,-1.0,Cbe,T,0.0,S);
    matmul("NN",1,3,3, 1.0,dopdv,S,0.0,dopbg);
}
/* jacobian of doppler measurement by lever arm------------------------------*/
static void jaco_dop_dl(const double *dopdv,const double *dr,const double *Cbe,
                        const double *gyro,double *dopdl)
{
    int i;
    double T[9],S[9];
    matmul("NN",1,3,3,1.0,dr,Cbe,0.0,dopdl);

    skewsym3(gyro,T);
    matmul("NN",3,3,3,1.0,Cbe,T,0.0,S);
    matmul("NN",1,3,3,1.0,dopdv,S,0.0,T);
    for (i=0;i<3;i++) dopdl[i]+=T[i];

    matmul("NN",3,3,3,1.0,Omge,Cbe,0.0,S);
    matmul("NN",1,3,3,1.0,dopdv,S,0.0,T);
    for (i=0;i<3;i++) dopdl[i]+=T[i];
}
/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jaco_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
}
/* get receiver position/velocity in ecef------------------------------------*/
static void getrcvp(const insstate_t *ins,double *rr,double *vr)
{
    rmlever(NULL,ins->re,ins->ve,ins->lever,ins->Cbe,ins->omgb,rr,vr);
}
/* sensitive matrix for doppler measurement----------------------------------*/
static int doppHVR(const prcopt_t *opt,const insstate_t *ins,const double *rs,
                   const double *dts,const int *svh,int n,const obsd_t *obs,
                   const nav_t *nav,double *v,double *H,double *R)
{
    double rr[3],vr[3],vs[3],e[3],azel[2],pos[3],rate,lam,dtr,*r;
    double dopbg[3],dopda[3],dopdl[3],dopdv[3],dopdp[3];
    double S[9],dap[3];
    int i,j,nv=0,nx=ins->nx;
    int ibg,nbg,ila,nla,irr,IP,IA,IV,NP,NA,NV;
    const insopt_t *opt_=&opt->insopt;

    trace(3,"doppHVR:\n");

    r=mat(n,1);

    getrcvp(ins,rr,vr);
    ecef2pos(rr,pos); 

    dtr=ins->dtrr; /* receiver clock drift */

    ibg=xiBg(opt_);
    nbg=xnBg(opt_);
    ila=xiLa(opt_);
    nla=xnLa(opt_);
    irr=xiRr(opt_); /* estimated states index */

    IP=xiP(opt_); NP=xnP(opt_);
    IA=xiA(opt_); NA=xnA(opt_);
    IV=xiV(opt_); NV=xnV(opt_);

    /* sensitive matrix for doppler measurement */
    for (i=0;i<n&&i<MAXOBS;i++) {

        lam=nav->lam[obs[i].sat-1][0];

        if (obs[i].D[0]==0.0||lam==0.0||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        if (!satsys(obs[i].sat,NULL)) continue;

        /* geometric distance/azimuth/elevation angle */
        if (geodist(rs+i*6,rr,e)<=0.0||satazel(pos,e,azel)<opt->elmin) continue;

        /* excluded satellite */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* satellite velocity relative to receiver in ecef */
        vs[0]=rs[0+3+i*6]-vr[0];
        vs[1]=rs[1+3+i*6]-vr[1];
        vs[2]=rs[2+3+i*6]-vr[2];

        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*vr[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*vr[1]);

        /* doppler residual */
        v[nv]=-lam*obs[i].D[0]-(rate+dtr-CLIGHT*dts[1+i*2]);

        if (fabs(v[nv])<MAXIONV) {

            /* jacobian matrix */
            if (H) {

                jaco_dop_dp(rs+6*i,rr,rs+3+6*i,vr,dopdp);
                jaco_dop_dv(e,rs+6*i,dopdv);
                jaco_dop_da(dopdp,ins->Cbe,ins->lever,dopdv,ins->omgb,dopda);
                jaco_dop_bg(dopdv,ins->Cbe,ins->lever,dopbg);
                jaco_dop_dl(dopdv,dopdp,ins->Cbe,ins->omgb,dopdl);
#if UPD_IN_EULER
                jaco_prot_pang(ins->Cbe,S);
                matcpy(dap,dopda,1,3);
                matmul("NN",1,3,3,1.0,dap,S,0.0,dopda);
#endif
                for (j=IP;j<IP+NP;j++) H[j+nv*nx]=dopdp[j-IP];
                for (j=IA;j<IA+NA;j++) H[j+nv*nx]=dopda[j-IA];
                for (j=IV;j<IV+NV;j++) H[j+nv*nx]=dopdv[j-IV];

                H[irr+nv*nx]=1.0; /* receiver clock drift */

                if (nbg) {
                    H[ibg+0+nv*nx]=dopbg[0];
                    H[ibg+1+nv*nx]=dopbg[1];
                    H[ibg+2+nv*nx]=dopbg[2];
                }
                if (nla) {
                    H[ila+0+nv*nx]=dopdl[0];
                    H[ila+1+nv*nx]=dopdl[1];
                    H[ila+2+nv*nx]=dopdl[2];
                }
            }
            /* measurement variance */
            r[nv++]=STDOPP/sin(azel[1]);
        }
    }
    if (R) {
        for (i=0;i<nv;i++) R[i+i*nv]=r[i];
        trace(3,"R=\n");
        tracemat(5,R,nv,nv,15,6);
    }
    if (nv&&v) {
        trace(3,"v=\n");
        tracemat(5,v,nv,1 ,15,6);
    }
    if (nv&&H) {
        trace(3,"H=\n");
        tracemat(5,H,nx,nv,15,6);
    }
    free(r);
    return nv;
}
/* doppler updates solutions valid-------------------------------------------*/
static int valdopp(const insopt_t *opt,const double *x,const double *R,
                   const double *v,int nv,double thres)
{
    double fact=SQR(thres);
    int nba,nbg,iba,ibg,i;

    trace(3,"valdoop:\n");

    nba=xnBa(opt); iba=xiBa(opt);
    nbg=xnBg(opt); ibg=xiBg(opt);

    /* check estimated states */
    if (norm(x,3)>15.0*D2R||(nba?norm(x+iba,3)>1E5*Mg2M:false)
        ||(nbg?norm(x+ibg,3)>30.0*D2R:false)) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    /* post-fit residual test */
    for (i=0;i<nv;i++) {
        if (v[i]*v[i]<fact*R[i+i*nv]) continue;
        trace(2,"large residual (v=%6.3f sig=%.3f)\n",v[i],SQRT(R[i+i*nv]));
        return 0;
    }
    return 1;
}
/* doppler filter for ins states updates-------------------------------------*/
static int doppfilt(insstate_t *ins,const prcopt_t *opt,const obsd_t *obs,int n,
                    const nav_t *nav,const double *rs,const double *dts,
                    const int *svh)
{
    double *v,*H,*R,*x,*P;
    int nx=ins->nx,nv,info=0,irr=0,nrr=0,IP,IV;
    static insstate_t inss;

    trace(3,"doppfilt:\n");

    v=mat(n,1); R=zeros(n,n); H=zeros(n,nx);
    x=zeros(nx,1); P=mat(nx,nx);

    matcpy(P,ins->P,nx,nx);

    irr=xiRr(&opt->insopt);
    nrr=xnRr(&opt->insopt);
    IP=xiP(&opt->insopt); IV=xiV(&opt->insopt);

    /* sensitive matrix for doppler measurement */
    nv=doppHVR(opt,ins,rs,dts,svh,n,obs,nav,v,H,R);

    if (nv>3) {

        /* receiver clock drift and non-zero due to close-loop */
        x[irr]=ins->dtrr;

        /* initialize every epoch for clock drift (white noise) */
        initP(irr,nrr,nx,opt->insopt.unc.rr,UNC_CLKR,P);

        /* ekf filter */
        if (filter(x,P,H,v,R,nx,nv)) {
            trace(2,"filter error\n");
        }
        else {
            inss.re[0]=ins->re[0]-x[IP+0];
            inss.re[1]=ins->re[1]-x[IP+1];
            inss.re[2]=ins->re[2]-x[IP+2];

            inss.ve[0]=ins->ve[0]-x[IV+0];
            inss.ve[1]=ins->ve[1]-x[IV+1];
            inss.ve[2]=ins->ve[2]-x[IV+2];

            inss.dtrr=x[irr];

            /* postfit doppler residuals */
            if (doppHVR(opt,&inss,rs,dts,svh,n,obs,nav,v,NULL,NULL)) {

                /* valid solutions */
                info=valdopp(&opt->insopt,x,R,v,nv,4.0);

                /* close loop */
                if (info) {
                    matcpy(ins->P,P,nx,nx);

                    ins->gstat=SOLQ_DOP;
                    ins->dtrr=x[irr];
                    ins->ns=nv;
                    clp(ins,&opt->insopt,x);
                }
            }
        }
    }
    free(x); free(P);
    free(v); free(H); free(R);
    return info;
}
/* ins states updates combination with doppler measurement-------------------
 * args    :  obsd_t *obs      I   observation data
 *            nav_t *nav       I   navigation data
 *            prcopt_t *opt    I   options
 *            insstate_t *ins  IO  ins states
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int doppler(const obsd_t *obs,int n,const nav_t *nav,const prcopt_t *opt,
                   insstate_t *ins)
{
    double *rs,*dts,*var;
    int i,no,stat,svh[MAXOBS];
    static obsd_t obsd[MAXOBS];

    trace(3,"doppler:\n");

    trace(4,"obs=\n"); traceobs(4,obs,n);
    trace(5,"nav=\n"); tracenav(5,nav);

    /* collect rover station observation */
    for (i=0,no=0;i<n;i++) {
        if (obs[i].rcv==1) obsd[no++]=obs[i];
    }
    if (no<=3) {
        trace(2,"no observation data\n");
        return 0;
    }
    rs=mat(6,n); dts=mat(2,n);
    var=mat(1,n);

    /* satellite positons,velocities and clocks */
    satposs(obs[0].time,obsd,no,nav,opt->sateph,rs,dts,var,svh);

    /* doppler measurement aid ins states updates */
    stat=doppfilt(ins,opt,obsd,no,nav,rs,dts,svh);

    if (stat) {
        trace(3,"doppler aid ins ok\n");
        ins->stat=INSS_DOPP;
    }
    else {
        trace(2,"doppler aid ins fail\n");
        stat=0;
    }
    free(rs); free(dts); free(var);
    return stat;
}
