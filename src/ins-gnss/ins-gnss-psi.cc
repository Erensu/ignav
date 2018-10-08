/*------------------------------------------------------------------------------
* ins-gnss-psi.cc : ins-gnss loosely coupled based on psi formulation
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
* history : 2017/10/06 1.0 new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants-----------------------------------------------------------------*/
#define E_SQR   0.00669437999014    /* sqr of linear eccentricity of the ellipsoid */
#define MAXDT   60.0                /* max interval to update imu (s) */
#define NORM_QUAT 1                 /* quaternion normalization due to numerical errors */
#define UPD_POS_IN_QUAT 0           /* update position in quaternion formation */

/* constants/global variables -----------------------------------------------*/
static double I[9]={1,0,0,0,1,0,0,0,1};

/* computes the curvature matrix and the gravity-----------------------------*/
static __inline void geoparam(const double *pos,const double *vn,double *wen_n,
                              double *wie_n,double *g,double *Reo,double *Rno)
{
    static double a1,a2,a3,a4,a5,a6;
    double sr,Re,Rn,sL2,sL4;

    /* local radii */
    sr=sqrt(1.0-E_SQR*SQR(sin(pos[0])));
    Re=(RE_WGS84/sr)+pos[2];
    Rn=(RE_WGS84*(1.0-E_SQR)/(sr*sr*sr))+pos[2];

    if (Reo) *Reo=Re;
    if (Rno) *Rno=Rn;

    /* local geodetic frame implementation */
    if (wen_n) {
        wen_n[0]= vn[1]/Re;
        wen_n[1]=-vn[0]/Rn;
        wen_n[2]=-vn[1]*tan(pos[0])/Re;
    }
    if (wie_n) {
        wie_n[0]= OMGE*cos(pos[0]);
        wie_n[1]= 0.0;
        wie_n[2]=-OMGE*sin(pos[0]);
    }
    /* plump bob gravity (see:Heiskanen and Moritz (1967))*/
    if (g) {
        sL2=SQR(sin(pos[0]));
        sL4=SQR(sL2);

        a1= 9.7803267715;
        a2= 0.0052790414;
        a3= 0.0000232718;
        a4=-0.0000030876910891;
        a5= 0.0000000043977311;
        a6= 0.0000000000007211;
        *g=a1*(1+a2*sL2+a3*sL4)+(a4+a5*sL2)*pos[2]+a6*SQR(pos[2]);
    }
}
/* convert quaternion to {lat,lon}-------------------------------------------*/
static __inline void quat2pos(double q[4],double *r)
{
    r[0]=-2*atan(q[2]/q[0])-PI/2;
    r[1]=2*atan2(q[3],q[0]);
}
/* update velocity-----------------------------------------------------------*/
static void updatevel(insstate_t *ins,const imud_t *data,double dt,double *rn2)
{
    double g,wen_n[3],wie_n[3],dap[3],dvk[3],dvp[3],dak[3];
    double dv1[3],dv2[3],dv3[3],w[3];
    double qne[4],dqn[4],dqe[4],Cne[9],dq1[4],vn2[3];
    double win_n[3],wie_e[3]={0,0,OMGE};
    double dC[9];
    int i;

    trace(3,"updatevel:\n");

    /* geo parameters */
    geoparam(ins->rn,ins->vn,wen_n,wie_n,&g,NULL,NULL);

    /* 1.velocity increment due to specific force  */
    for (i=0;i<3;i++) {
        dap[i]=ins->omgbp[i]*ins->dt; dvp[i]=ins->fbp[i]*ins->dt;
    }
    for (i=0;i<3;i++) {
        dak[i]=data->gyro[i]*dt; dvk[i]=data->accl[i]*dt;
    }
    cross3(dak,dvk,dv1);
    cross3(dap,dvk,dv2);
    cross3(dvp,dak,dv3);
    for (i=0;i<3;i++) {
        dvk[i]=dvk[i]+0.5*dv1[i]+1.0/12.0*(dv2[i]+dv3[i]);
    }
    /* midway velocity */
    for (i=0;i<3;i++) vn2[i]=ins->vn[i]+0.5*(ins->dvn[i]);

    /* extrapolated position in midway */
    ned2xyz(ins->rn,Cne); dcm2quatx(Cne,qne);
    for (i=0;i<3;i++) {
        win_n[i]=(wie_n[i]+wen_n[i])*dt/2.0; wie_e[i]*=-dt/2.0;
    }
    rv2quat(win_n,dqn); rv2quat(wie_e,dqe);

    quatmulx(qne,dqn,dq1);
    quatmulx(dqe,dq1,qne);

    /* get midway position */
    quat2pos(qne,rn2);
    rn2[2]=ins->rn[2]-0.5*ins->vn[2]*dt;

    /* geo parameters in midway */
    geoparam(rn2,vn2,wen_n,wie_n,&g,NULL,NULL);

    for (i=0;i<3;i++) w[i]=(wie_n[i]+wen_n[i])*dt/2.0;
    skewsym3(w,dC);
    for (i=0;i<9;i++) dC[i]=I[i]-dC[i];

    matmul33("NNN",dC,ins->Cbn,dvk,3,3,3,1,dv1);

    /* 2.velocity increment due to gravity and Coriolis force */
    for (i=0;i<3;i++) w[i]=2.0*wie_n[i]+wen_n[i];

    cross3(w,vn2,dv2);
    for (i=0;i<3;i++) {
        dv2[i]=(i<2?-dv2[i]:g-dv2[i])*dt; /* gn=[0,0,g] */
    }
    /* update velocity in n-frame */
    ins->vn[0]+=(ins->dvn[0]=dv1[0]+dv2[0]);
    ins->vn[1]+=(ins->dvn[1]=dv1[1]+dv2[1]);
    ins->vn[2]+=(ins->dvn[2]=dv1[2]+dv2[2]);

    /* update acceleration in n-frame */
    for (i=0;i<3;i++) ins->an[i]=ins->dvn[i]/dt;
}
/* update position-----------------------------------------------------------*/
static void updatepos(insstate_t *ins,const imud_t *data,double dt,double *rn2)
{
    double vn2[3],Cne[9],Re,Rn,qne[4],wie_n[3],wen_n[3];
    double win_n[3],wie_e[3]={0,0,OMGE};
    double dqn[4],dqe[4],q[4];
    int i;

    trace(3,"updatepos:\n");

    ned2xyz(ins->rn,Cne); dcm2quatx(Cne,qne);

    /* velocity in midway */
    for (i=0;i<3;i++) {
        vn2[i]=(ins->vn[i]*2.0-ins->dvn[i])/2.0;
    }
    /* geo parameters in midway */
    geoparam(rn2,vn2,wen_n,wie_n,NULL,&Re,&Rn);

#if UPD_POS_IN_QUAT
    for (i=0;i<3;i++) {
        win_n[i]=(wie_n[i]+wen_n[i])*dt; wie_e[i]*=-dt;
    }
    rv2quat(win_n,dqn);
    rv2quat(wie_e,dqe);

    quatmulx(qne,dqn,q); quatmulx(dqe,q,qne);

    /* update position */
    quat2pos(qne,ins->rn);
    ins->rn[2]-=vn2[2]*dt;
#else
    /* update position */
    ins->rn[0]+=1.0/Rn*vn2[0]*dt;
    ins->rn[1]+=1.0/(Re*cos(rn2[0]))*vn2[1]*dt;
    ins->rn[2]-=vn2[2]*dt;
#endif
}
/* quaternion normalization--------------------------------------------------*/
static __inline void norm_quat(double *qbn)
{
    int i; double eq;
    for (eq=0.5*(norm(qbn,4)-1.0),i=0;i<4;i++) qbn[i]*=(1.0-eq);
}
/* update attitude ----------------------------------------------------------*/
static void updateatt(insstate_t *ins,const imud_t *data,double dt)
{
    double dqb[4],dqn[4],dq[4],qbn[4],dap[3],dak[3];
    double win_n[3],wen_n[3],wie_n[3],da[3];
    int i;

    trace(3,"updateatt:\n");

    for (i=0;i<3;i++) {
        dap[i]=ins->omgbp[i]*ins->dt; dak[i]=data->gyro[i]*dt;
    }
    /* geo parameters */
    geoparam(ins->rn,ins->vn,wen_n,wie_n,NULL,NULL,NULL);
    for (i=0;i<3;i++) {
        win_n[i]=-(wie_n[i]+wen_n[i])*dt/2.0;
    }
    rv2quat(win_n,dqn);

    cross3(dap,dak,da);
    for (i=0;i<3;i++) da[i]=1.0/12.0*da[i]+dak[i];
    rv2quat(da,dqb);

    dcm2quatx(ins->Cbn,qbn);

    quatmulx(qbn,dqb,dq);
    quatmulx(dqn,dq,qbn);
#if NORM_QUAT
    /* quaternion normalization to eliminate numerical error */
    norm_quat(qbn);
#endif
    /* update ins attitude */
    quat2dcmx(qbn,ins->Cbn);
}
/* save precious epoch ins states--------------------------------------------*/
static void savepins(insstate_t *ins,const imud_t *data)
{
    matcpy(ins->omgbp,data->gyro,1,3);
    matcpy(ins->fbp  ,data->accl,1,3);
    matcpy(ins->pins  ,ins->re ,1,3);
    matcpy(ins->pins+3,ins->ve ,1,3);
    matcpy(ins->pCbe  ,ins->Cbe,3,3);
}
/* update ins states ----------------------------------------------------------
* update ins states based on Llh position mechanization
* args   : insopt   *insopt I   ins updates options
*          insstate_t *ins  IO  ins states
*          imudata_t *data  I   imu measurement data
* return : 0 (fail) or 1 (ok)
*----------------------------------------------------------------------------*/
extern int updateinsn(const insopt_t *insopt,insstate_t *ins,const imud_t *data)
{
    double dt,rmid[3];
    int i;

    trace(3,"updateins:\n");

    trace(5,"ins(-)=\n"); traceins(5,ins);

    if ((dt=timediff(data->time,ins->time))>MAXDT||fabs(dt)<1E-6) {
        trace(2,"time difference too large: %.0fs\n",dt);
        return 0;
    }
    for (i=0;i<3;i++) {
        ins->omgb0[i]=data->gyro[i];
        ins->fb0  [i]=data->accl[i];
        if (insopt->exinserr) {
            ins_errmodel(data->accl,data->gyro,ins->fb,ins->omgb,ins);
        }
        else {
            ins->omgb[i]=data->gyro[i]-ins->bg[i];
            ins->fb[i]  =data->accl[i]-ins->ba[i];
        }
    }
    /* save precious ins states */
    savepins(ins,data);

    /* update velocity/position/attitude in n-frame */
    updatevel(ins,data,dt,rmid);
    updatepos(ins,data,dt,rmid);
    updateatt(ins,data,dt);

    /* update ins states in e-frame */
    update_ins_state_e(ins);

    ins->dt=timediff(data->time,ins->time);
    ins->ptime=ins->time;
    ins->time =data->time;
    ins->stat =INSS_MECH;

    trace(5,"ins(+)=\n"); traceins(5,ins);
    return 1;
}