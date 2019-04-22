/*------------------------------------------------------------------------------
* ins-gnss-tc.c : ins-gnss tightly coupled common functions
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
* history : 2017/10/19 1.0 new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants -----------------------------------------------------------------*/
#define MAXVAR       1E10         /* max variance for reset covariance matrix */
#define MAXSOLR      2            /* max number of reboot solutions */
#define MINVEL       3.0          /* min velocity for initial ins states */
#define MAXGYRO      (30.0*D2R)   /* max rotation speed value for initial */
#define MAXDIFF      30.0         /* max time difference between solution */
#define REBOOT       1            /* ins tightly coupled reboot if always update fail */
#define CHKNUMERIC   1            /* check numeric for given value */

/* solution convert to velocity----------------------------------------------*/
static void sol2vel(const sol_t *sol1,const sol_t *sol2,double *v)
{
    v[0]=(sol1->rr[0]-sol2->rr[0])/timediff(sol1->time,sol2->time);
    v[1]=(sol1->rr[1]-sol2->rr[1])/timediff(sol1->time,sol2->time);
    v[2]=(sol1->rr[2]-sol2->rr[2])/timediff(sol1->time,sol2->time);
}
/* check ins states covariance matrix----------------------------------------*/
static int chkpcov(int nx,const insopt_t *opt,double *P)
{
    int i; double var=0.0;
    for (i=xiP(opt);i<xiP(opt)+3;i++) var+=SQRT(P[i+i*nx]);

    if ((var/3)>MAXVAR) if (P) getP0(opt,P);
}
/* check imu body velocity---------------------------------------------------*/
static int chkvb(const insstate_t *ins)
{
#if 1
    double vb[3];
    matmul("TN",3,1,3,1.0,ins->Cbe,ins->ve,0.0,vb);
    return fabs(vb[1])<MINVEL&&fabs(vb[2])<MINVEL;
#else
    return 0;
#endif
}
/* initialization position mode/ionosphere and troposphere option------------*/
static void initrtkpos(rtk_t *rtk,prcopt_t *prcopt)
{
    prcopt->mode   =PMODE_SINGLE;
    prcopt->ionoopt=IONOOPT_BRDC;
    prcopt->tropopt=TROPOPT_SAAS;
    rtkinit(rtk,prcopt);
}
/* reboot ins states---------------------------------------------------------*/
static void rebootsta(prcopt_t *opt,insstate_t *ins)
{
    trace(3,"initinsrt:\n");
    freelc(ins);
    initlc (&opt->insopt      ,ins);
    initodo(&opt->insopt.odopt,ins);
}
/* reboot tightly couple if no coupling for a long time----------------------
 * args   : insstate_t *ins  IO ins states
 *          prcopt_t *opt    I  options
 *          obsd_t *obs      I  observation data
 *          int n            I  number of observation data
 *          imud_t *imu      I  imu measurement data
 *          nav_t *nav       I  navigation data
 * return : 0: no reboot, 1: reboot ok
 * --------------------------------------------------------------------------*/
static int rebootc(insstate_t *ins,const prcopt_t *opt,const obsd_t *obs,int n,
                   const imud_t *imu,const nav_t *nav)
{
    int i; double vr[3]={0};

    static prcopt_t prcopt=*opt;
    static sol_t sols[MAXSOLR]={0},sol0={0};
    static int first=1;
    static rtk_t rtk={0};

    trace(3,"rebootc: n=%d\n",n);

    if (obs==NULL||n<=0||chkvb(ins)) return 0;

    /* initial rtk position when first time  */
    if (first) {
        initrtkpos(&rtk,&prcopt); first=0;
    }
    rtkpos(&rtk,obs,n,nav);

    /* save position solution to buffer */
    for (i=0;i<MAXSOLR-1;i++) sols[i]=sols[i+1]; sols[i]=rtk.sol;

    /* check solution status */
    for (i=0;i<MAXSOLR;i++) {
        if (sols[i].stat>SOLQ_PPP||sols[i].stat==SOLQ_NONE) {
            trace(2,"solution status check fail\n");
            return 1;
        }
    }
    /* check time-continuity of solutions */
    for (i=0;i<MAXSOLR-1;i++) {
        if (timediff(sols[i+1].time,sols[i].time)>MAXDIFF) {
            trace(2,"large time difference of solution\n");
            return 1;
        }
        if (fabs(timediff(sols[i+1].time,sols[i].time))<1E-5) {
            trace(2,"duplicate gps measurement\n");
            return 1;
        }
    }
    matcpy(vr,sols[MAXSOLR-1].rr+3,1,3);
    if (norm(vr,3)==0.0) {
        sol2vel(sols+MAXSOLR-1,sols+MAXSOLR-2,vr);
    }
    if (norm(vr,3)<MINVEL||norm(imu->gyro,3)>MAXGYRO) {
        trace(2,"reboot ins state fail\n");
        return 1;
    }
    /* reboot ins states */
    rebootsta(&prcopt,ins);
    if (!ant2inins(sols[MAXSOLR-1].time,sols[MAXSOLR-1].rr,
                   vr,&opt->insopt,NULL,ins,NULL)) {
        trace(2,"reboot ins state fail\n");
        return 1;
    }
    ins->time=sols[MAXSOLR-1].time;

    for (i=0;i<MAXSOLR;i++) sols[i]=sol0;
    rtkfree(&rtk); first=1;

    trace(3,"reboot ins tightly coupled ok\n");
    return 2;
}
/* ins and gnss tightly coupled-----------------------------------------------
 * args   :  prcopt_t *opt   I  processing options
 *           obsd_t *obs     I  observation data
 *           int n           I  number of observation data
 *           nav_t *nav      I  navigation data
 *           imud_t *imu     I  imu measurement data
 *           rtk_t *rtk      IO rtk struct data
 *           insstate_t *ins IO ins states
 *           int upd         I  updates flag (INSUPD_???)
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int tcigpos(const prcopt_t *opt,const obsd_t *obs,int n,const nav_t *nav,
                   const imud_t *imu,rtk_t *rtk,insstate_t *ins,int upd)
{
    int i,nx=ins->nx,info=1,nu,nr,flag;
    double *P,dt;
    const insopt_t* insopt=&opt->insopt;

    trace(3,"tcigpos: update=%d,time=%s\n",upd,time_str(imu->time,4));

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
    ins->stat=INSS_NONE; /* start ins mechanization */
    if (
#if 0
        /* update ins states based on llh position mechanization */
        !updateinsn(insopt,ins,imu);
#else
        /* update ins states in e-frame */
        !updateins(insopt,ins,imu)
        ) {
#endif
        trace(2,"ins mechanization update fail\n");
        return 0;
    }
    P=zeros(nx,nx);

    /* propagate ins states */
    propinss(ins,insopt,ins->dt,ins->x,ins->P);

    /* check variance of estimated position */
    chkpcov(nx,insopt,ins->P);

    /* updates ins states using imu or observation data */
    if (upd==INSUPD_TIME) {
        ins->stat=INSS_TIME;
        info=1;
    }
    else {
        for (i=0;i<6;i++) rtk->sol.pqr[i]=rtk->sol.qr[i];
        rtk->sol.pstat=rtk->sol.stat;
        ins->gstat=SOLQ_NONE;
        ins->ns=0;
#if REBOOT
        /* reboot tightly coupled if need */
        if ((flag=rebootc(ins,opt,obs,n,imu,nav))) {
            if (flag==1) {
                trace(2,"ins tightly coupled still reboot\n");
                info=0; goto EXIT;
            }
            trace(3,"ins tightly coupled reboot ok\n");
            ins->stat=INSS_REBOOT; info=1;
            goto EXIT;
        }
#endif
        /* updates by measurement data */
        if (obs&&imu&&n) {

            /* count rover/base station observations */
            for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;
            for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;

            dt=timediff(obs[0].time,ins->time);

            /* check synchronization */
            if (fabs(dt)>3.0) {
                trace(2,"observation and imu sync error\n");
                info=0;
            }
            /* tightly coupled */
            if (info) {
                info=rtkpos(rtk,obs,nu+nr,nav);
            }
        }
        else info=0;

        /* update coupled solution status */
        if (info) {
            ins->ptct=ins->time;
            ins->stat=ins->stat==INSS_REBOOT?INSS_REBOOT:INSS_TCUD;

            trace(3,"tightly couple ok\n");

            /* lack satellites but tightly-coupled run */
            if (ins->ns<4) {
                ins->stat=INSS_LACK;
            }
            /* save precious epoch gps measurement */
            savegmeas(ins,&rtk->sol,NULL);
#if 1
            /* recheck attitude */
            rechkatt(ins,imu);
#endif
            /* ins state in n-frame */
            update_ins_state_n(ins);
        }
        else {
            trace(2,"tightly coupled fail\n");
            info=0;
        }
    }
EXIT:
    free(P); return info;
}

