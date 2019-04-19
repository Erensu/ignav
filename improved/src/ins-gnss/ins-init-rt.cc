/*-----------------------------------------------------------------------------
* ins-init-rt.cc : real-time initialization for ins navigation
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
* history : 2017/01/21 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants ----------------------------------------------------------------*/
#define MINSOL        5                  /* max number of solution data */
#define MINVEL        5.0                /* min velocity for initial ins states */
#define MAXGYRO       (30.0*D2R)         /* max rotation speed value for initial */
#define MAXVAR_POSE   (5.0*D2R)          /* max variance of pose measurement */
#define MAXDIFF       10.0               /* max time difference between solution */
#define ADJOBS        1                  /* adjust observation data */
#define ONLY_INIT_YAW 0                  /* only initial yaw from dual ants. measurement data */

/* coordinate rotation matrix ------------------------------------------------*/
#define Rx(t,X) do {                                     \
    (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0;         \
    (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do {                                     \
    (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0;         \
    (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do {                                     \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0;         \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

/* solution convert to velocity----------------------------------------------*/
static void sol2vel(const sol_t *sol1,const sol_t *sol2,double *v)
{
    v[0]=(sol1->rr[0]-sol2->rr[0])/timediff(sol1->time,sol2->time);
    v[1]=(sol1->rr[1]-sol2->rr[1])/timediff(sol1->time,sol2->time);
    v[2]=(sol1->rr[2]-sol2->rr[2])/timediff(sol1->time,sol2->time);
}
/* initial ins states--------------------------------------------------------*/
static void initinsrt(rtksvr_t *svr)
{
    trace(3,"initinsrt:\n");
    freelc (&svr->rtk.ins);
    initlc (&svr->rtk.opt.insopt      ,&svr->rtk.ins);
    initodo(&svr->rtk.opt.insopt.odopt,&svr->rtk.ins);
}
/* check solution valid------------------------------------------------------*/
static int checksol(const sol_t *sols)
{
    if (sols->stat==SOLQ_NONE) return 0;
    if (sols->time.time==0) return 0;
    if (norm (sols->rr,6)==0.0) return 0;
    if (normf(sols->qr,3)==0.0) return 0;
    return 1;
}
/* initialization ins states for real-time navigation use PVT solution--------
 * args   :  rtksvr_t *svr    IO  rtk server
 *           sol_t *sol       I   PVT solution data
 *           imud_t *imu      I   imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int insinitrt(rtksvr_t *svr,const sol_t *sol,const imud_t *imu)
{
    insopt_t *iopt=&svr->rtk.opt.insopt;
    insstate_t *ins=&svr->rtk.ins;
    static sol_t sols[MINSOL]={0};
    double vr[3]={0};
    int i;

    trace(3,"insinitrt: time=%s\n",time_str(imu->time,4));

    /* check solution valid */
    if (!checksol(sol)) {
        trace(2,"invalid solution data\n"); return 0;
    }
    svr->rtk.ins.stat=INSS_INIT;
    
    /* save pvt solution buffer */
    for (i=0;i<MINSOL-1;i++) sols[i]=sols[i+1]; sols[i]=*sol;
    for (i=0;i<MINSOL;i++) {
        if (sols[i].stat>iopt->iisu||sols[i].stat==SOLQ_NONE) return 0;
    }
    /* compute velocity from solutions */
    matcpy(vr,sols[MINSOL-1].rr+3,1,3);
    if (norm(vr,3)==0.0) {
        sol2vel(sols+MINSOL-1,sols+MINSOL-2,vr);
    }
    /* check velocity ok? */
    if (norm(vr,3)<MINVEL||norm(imu->gyro,3)>MAXGYRO) return 0;
    
    for (i=0;i<MINSOL-1;i++) {
        if (timediff(sols[i+1].time,sols[i].time)>MAXDIFF) {
            trace(2,"large time difference of solution\n");
            return 0;
        }
        if (fabs(timediff(sols[i+1].time,sols[i].time))<1E-5) {
            trace(2,"duplicate gps measurement\n");
            return 0;
        }
    }
    /* initial ins states */
    initinsrt(svr);
    if (!ant2inins(sols[MINSOL-1].time,sols[MINSOL-1].rr,vr,iopt,NULL,ins,NULL)) {
        trace(2,"initial ins state fail\n");
        return 0;
    }
    ins->time=sols[MINSOL-1].time;

    /* update ins state in n-frame */
    updinsn(ins);
    matcpy(ins->oxyz,ins->re,3,1);

    trace(3,"initial ins state ok\n");
    return 1;
}
/* initialization position mode/ionosphere and troposphere option------------*/
static void initrtkpos(rtk_t *rtk,prcopt_t *prcopt)
{
    prcopt->mode   =PMODE_KINEMA;
    prcopt->ionoopt=IONOOPT_BRDC;
    prcopt->tropopt=TROPOPT_SAAS;
#if ADJOBS
    prcopt->adjobs=1;
#endif
    rtkinit(rtk,prcopt);
}
/* initialization ins states for real-time navigation use observation--------
 * args   :  rtksvr_t *svr  IO  rtk server
 *           obsd_t *obs    I   observation data
 *           int n          I   number of observation data
 *           imud_t *imu    I   imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int insinirtobs(rtksvr_t *svr,const obsd_t *obs,int n,const imud_t *imu)
{
    double vr[3]={0}; insstate_t *ins=&svr->rtk.ins;

    /* global variables for rtk positioning */
    static int first=1,i;
    static prcopt_t popt=svr->rtk.opt;
    static rtk_t rtk={0};
    static sol_t sols[MINSOL]={0};

    trace(3,"insinirtobs: n=%d\n",n);

    svr->rtk.ins.stat=INSS_INIT;
    if (n<=0) {
        trace(2,"no observation data to initial\n"); return 0;
    }
    /* initial gps position options */
    if (first) {
        initrtkpos(&rtk,&popt); first=0;
    }
    rtkpos(&rtk,obs,n,&svr->nav);

    /* save position solution to buffer */
    for (i=0;i<MINSOL-1;i++) sols[i]=sols[i+1]; sols[i]=rtk.sol;
    for (i=0;i<MINSOL;i++) {
        if (sols[i].stat>popt.insopt.iisu||sols[i].stat==SOLQ_NONE) {
            trace(2,"check solution status fail\n");
            return 0;
        }
    }
    for (i=0;i<MINSOL-1;i++) {
        if (timediff(sols[i+1].time,sols[i].time)>MAXDIFF) {
            return 0;
        }
    }
    /* compute velocity from solutions */
    matcpy(vr,sols[MINSOL-1].rr+3,1,3);
    if (norm(vr,3)==0.0) {
        sol2vel(sols+MINSOL-1,sols+MINSOL-2,vr);
    }
    if (norm(imu->gyro,3)>MAXGYRO||norm(vr,3)<MINVEL) {
        return 0;
    }
    /* initialize ins states */
    initinsrt(svr);
    if (!ant2inins(sols[MINSOL-1].time,sols[MINSOL-1].rr,vr,&popt.insopt,
                   NULL,ins,NULL)) {
        return 0;
    }
    ins->time=sols[MINSOL-1].time;

    /* update ins state in n-frame */
    updinsn(ins);
    matcpy(ins->oxyz,ins->re,3,1);

    /* reset rtk position options */
    rtkfree(&rtk); first=1;

    trace(3,"initial ins state ok\n");
    return 1;
}
/* initial ins states from dual antennas pose measurement--------------------
 * args:  rtksvr_t *svr     IO  rtk server
 *        pose_meas_t *pose I   pose measurement data
 *        sol_t *sol        I   solution data
 *        imud_t *imu       I   imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int insinitdualant(rtksvr_t *svr,const pose_meas_t *pose,const sol_t *sol,
                          const imud_t *imu)
{
    double Cne[9],vr[3],Cvn[9],Cbn[9],Ry[9],Rz[9],pos[3];
    static sol_t sols[MINSOL]={0};
    insopt_t *iopt=&svr->rtk.opt.insopt;
    insstate_t *ins=&svr->rtk.ins;
    pose_meas_t posem;
    int i;

    trace(3,"insinitdualant:\n");

    /* check solution valid */
    if (!checksol(sol)) {trace(2,"invalid solution data\n"); return 0;}
    svr->rtk.ins.stat=INSS_INIT;

    /* save pvt solution buffer */
    for (i=0;i<MINSOL-1;i++) sols[i]=sols[i+1];
    sols[i]=*sol;

    for (i=0;i<MINSOL;i++) {
        if (sols[i].stat>iopt->iisu||sols[i].stat==SOLQ_NONE) return 0;
    }
    /* check solution time continuity */
    for (i=0;i<MINSOL-1;i++) {
        if (timediff(sols[i+1].time,sols[i].time)>MAXDIFF) {
            return 0;
        }
    }
    /* velocity in ecef */
    matcpy(vr,sols[MINSOL-1].rr+3,1,3);
    if (norm(vr,3)==0.0) {
        sol2vel(sols+MINSOL-1,sols+MINSOL-2,vr);
    }
#if 0
    if (norm(imu->gyro,3)>MAXGYRO||norm(vr,3)<MINVEL) {
        return 0;
    }
#endif

#if ONLY_INIT_YAW
    if (SQRT(pose->var[2])>MAXVAR_POSE) {
        trace(2,"initial yaw fail due to large variance\n");
        return 0;
    }
    posem=*pose;
    posem.rpy[0]=posem.rpy[1]=0.0;
#else
    /* check pose measurement availability */
    if (SQRT(pose->var[0]+pose->var[1]+pose->var[2])>MAXVAR_POSE) {
        trace(2,"large pose variance\n");
        return 0;
    }
    posem=*pose;
#endif
    /* initial carvig-sever */
    initinsrt(svr);
    
    /* initial ins attitude from dual ant. */
    ecef2pos(sols[MINSOL-1].rr,pos);
    ned2xyz(pos,Cne);

    Ry(-(posem.rpy[1]-iopt->mis_euler[1]*D2R),Ry);
    Rz(-(posem.rpy[2]-iopt->mis_euler[2]*D2R),Rz);
    matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cbn);

    Ry(-posem.rpy[1],Ry);
    Rz(-posem.rpy[2],Rz);
    matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cvn);

    matmul("TN",3,3,3,1.0,Cbn,Cvn,0.0,ins->Cvb);

    matmul33("NNT",Cne,Cvn,ins->Cvb,3,3,3,3,ins->Cbe);
    gapv2ipv(sols[MINSOL-1].rr,vr,ins->Cbe,ins->lever,imu,
             ins->re,ins->ve);

    /* update ins state in n-frame */
    updinsn(ins);
    matcpy(ins->oxyz,ins->re,3,1);

    ins->time=imu->time;
    trace(3,"initial ins state ok\n");
    return 1;
}
/* initial ins states from user given pose------------------------------------
 * args:  rtksvr_t *svr     IO  rtk server
 *        imud_t *imu       I   imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int insinitgiven(rtksvr_t *svr,const imud_t *imu)
{
    insopt_t *iopt=&svr->rtk.opt.insopt;
    insstate_t *ins=&svr->rtk.ins;
    double Cnb[9],Cne[9],rpy[3],rn0[3];
    int i;

    trace(3,"insinitgiven:\n");

    if (fabs(iopt->t0-time2gst(imu->time,NULL))>DTTOL) return 0;
    ins->stat=INSS_INIT;

    /* initial carvig-sever */
    initinsrt(svr);

    for (i=0;i<3;i++) rpy[i]=iopt->att0[i]*D2R;
    rn0[0]=iopt->rn0[0]*D2R;
    rn0[1]=iopt->rn0[1]*D2R;
    rn0[2]=iopt->rn0[2];

    rpy2dcm(rpy,Cnb);
    ned2xyz(rn0,Cne);
    matmul("NT",3,3,3,1.0,Cne,Cnb,0.0,ins->Cbe);
    pos2ecef(rn0,ins->re);
    matmul("NN",3,1,3,1.0,Cne,iopt->vn0,0.0,ins->ve);

    /* update ins state in n-frame */
    updinsn(ins);
    matcpy(ins->oxyz,ins->re,3,1);

    ins->time=imu->time;
    trace(3,"initial ins state ok\n");
    return 1;
}

