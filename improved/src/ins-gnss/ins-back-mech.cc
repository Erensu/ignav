/*----------------------------------------------------------------------------
* ins-back-mech.cc : backward ins mechanization functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*    [5] Shin E H. Estimation techniques for low-cost inertial navigation,2005
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/03/02 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

#define MAXDT      60.0               /* max interval to update imu (s) */
#define INSUPDPRE  1                  /* inertial navigation equations precision */
#define E_SQR      0.00669437999014   /* sqr of linear eccentricity of the ellipsoid */
#define UPD_INS_E  0                  /* updates ins states in e-frame */

/* save precious epoch ins states--------------------------------------------*/
static void savepins(insstate_t *ins,const imud_t *data)
{
    matcpy(ins->omgbp ,ins->omgb,1,3);
    matcpy(ins->fbp   ,ins->fb  ,1,3);
    matcpy(ins->pins  ,ins->re  ,1,3);
    matcpy(ins->pins+3,ins->ve  ,1,3);
    matcpy(ins->pCbe  ,ins->Cbe ,3,3);
}
/* update ins attitude --------------------------------------------------------
 * args    : double t     I    time interval between epochs (s)
 *           double *Cbe  I    previous body-to-ecef coordinate transformation matrix
 *                        O    uipdates body-to-ecef coordinate transformation matrix
 *           double *omgb I    angular rate of body frame (rad/s) w.r.t eci-frame
 *                             expressed in ecef-frame
 *           double *das  I    rotational and sculling motion correction
 * return  :none
 * ---------------------------------------------------------------------------*/
static void updateatt(double dt, double *Cbe, const double *omgb,const double *das)
{
    double alpha[3],a,a1,a2,Ca[9],Ca2[9],Cbep[9];
    double Cbb[9]={1,0,0,0,1,0,0,0,1},Cei[9]={0};
    int i;

    trace(3,"updateatt: dt=%.4lf\n",dt);

    for (i=0;i<3;i++) alpha[i]=omgb[i]*dt+das[i];
    skewsym3(alpha,Ca);
    matmul3("NN",Ca,Ca,Ca2);
    a=norm(alpha,3);
    if (a<1E-8) {
        a1=1.0-a*a/6.0; a2=0.5-a*a/24.0;
    }
    else {
        a1=sin(a)/a; a2=(1.0-cos(a))/(a*a);
    }
    for (i=0;i<9;i++) {
        Cbb[i]+=a1*Ca[i]+a2*Ca2[i];
    }
    Cei[0]= cos(OMGE*dt); Cei[3]=-sin(OMGE*dt);
    Cei[1]= sin(OMGE*dt); Cei[4]= cos(OMGE*dt); Cei[8]=1.0;
    matmul3("NN",Cei,Cbe,Cbep);
    matmul3("NT",Cbep,Cbb,Cbe);
}
/* backward update ins states ------------------------------------------------
* update ins states with imu measurement data in e-frame in backward
* args   : insopt   *insopt I   ins updates options
*          insstate_t *ins  IO  ins states
*          imud_t    *data  I   imu measurement data
* return : 0 (fail) or 1 (ok)
*----------------------------------------------------------------------------*/
extern int updateinsbe(const insopt_t *insopt,insstate_t *ins,const imud_t *data)
{
    double dt,Cbe[9],fe[3],ge[3],cori[3],Cbb[9]={1,0,0,0,1,0,0,0,1};
    double Ca[9],Ca2[9],a1,a2,a,alpha[3]={0},Omg[9]={0},ae[3]={0};
    double das[3]={0},dvs[3]={0},fb[3];
    int i;

    trace(3,"updateinsb:\n");

    trace(5,"ins(-)=\n"); traceins(5,ins);

    /* save precious ins states */
    savepins(ins,data);

    if ((dt=-timediff(data->time,ins->time))>MAXDT||fabs(dt)<1E-6) {

        /* update time information */
        ins->dt=timediff(data->time,ins->time);
        ins->ptime=ins->time;
        ins->time =data->time;

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
            ins->fb  [i]=data->accl[i]-ins->ba[i];
        }
    }
    matcpy(Cbe,ins->Cbe,3,3);
    ae[2]=OMGE*dt;

#if SCULL_CORR
    /* rotational and sculling motion correction */
    rotscull_corr(ins,insopt,dt,dvs,das);
#endif
    /* update attitude */
    updateatt(dt,ins->Cbe,ins->omgb,das);

#if INSUPDPRE
    for (i=0;i<3;i++) alpha[i]=ins->omgb[i]*dt+das[i];
    skewsym3(alpha,Ca);
    /* check if the value is too small to keep numerical robustness */
    if ((a=norm(alpha,3))>1E-8) {
        a1=(1.0-cos(a))/SQR(a); a2=1.0/SQR(a)*(1.0-sin(a)/a);
        matmul3("NN",Ca,Ca,Ca2);
        for (i=0;i<9;i++) Cbb[i]+=a1*Ca[i]+a2*Ca2[i]; 
        skewsym3(ae,Omg);
        matmul3("NN",Cbe,Cbb,Ca);
        matmul3("NN",Omg,Cbe,Ca2);
        for (i=0;i<9;i++) Cbe[i]=Ca[i]-0.5*Ca2[i];
    }
    else {
        skewsym3(ae,Omg);
        matmul3("NN",Omg,Cbe,Ca);
        for (i=0;i<9;i++) Cbe[i]-=0.5*Ca[i];
    }
#else
    for (i=0;i<9;i++) Cbe[i]=(Cbe[i]+ins->Cbe[i])/2.0;
#endif
    /* specific-force/gravity in e-frame */
    for (i=0;i<3;i++) fb[i]=ins->fb[i]+dvs[i]/dt;
    matmul3v("N",Cbe,fb,fe);
    if (insopt->gravityex) {
        pregrav(ins->re,ge); /* precious gravity model */
    }
    else gravity(ins->re,ge);

    /* update velocity/position */
    matmul3v("N",Omge,ins->ve,cori);
    for (i=0;i<3;i++) {
        ins->ae[i]=fe[i]+ge[i]-2.0*cori[i]+dvs[i]/dt;
        ins->ve[i]-=ins->ae[i]*dt;
        ins->re[i]-=ins->ve[i]*dt+ins->ae[i]/2.0*dt*dt;
    }
    matcpy(ins->omgbp,data->gyro,1,3);
    matcpy(ins->fbp  ,data->accl,1,3);

    /* update ins stats in n-frame */
    updinsn(ins);

    ins->dt=timediff(data->time,ins->time);
    ins->ptime=ins->time;
    ins->time =data->time;
    ins->stat =INSS_MECH;

    trace(5,"ins(+)=\n"); traceins(5,ins);
    return 1;
}
/* computes the curvature matrix and the gravity-----------------------------*/
static void geoparam(const double *pos,const double *vn,double *wen_n,
                     double *wie_n,double *g)
{
    static double a1,a2,a3,a4,a5,a6;
    double sr,Re,Rn,sL2,sL4;

    /* local radii */
    sr=sqrt(1.0-E_SQR*SQR(sin(pos[0])));
    Re=(RE_WGS84/sr)+pos[2];
    Rn=(RE_WGS84*(1-E_SQR)/(sr*sr*sr))+pos[2];

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
/* quaternion multiply-------------------------------------------------------*/
extern void quatmulx(const double *qab,const double *qca,double *qcb)
{
    setzero(qcb,4,1);
    qcb[0]=(qab[0])*qca[0]+(-qab[1])*qca[1]+(-qab[2])*qca[2]+(-qab[3])*qca[3];
    qcb[1]=(qab[1])*qca[0]+( qab[0])*qca[1]+(-qab[3])*qca[2]+( qab[2])*qca[3];
    qcb[2]=(qab[2])*qca[0]+( qab[3])*qca[1]+( qab[0])*qca[2]+(-qab[1])*qca[3];
    qcb[3]=(qab[3])*qca[0]+(-qab[2])*qca[1]+( qab[1])*qca[2]+( qab[0])*qca[3];
}
/* convert quaternion into rotation vector-----------------------------------
 * note: q[1:4]=(w,x,y,z)
 * --------------------------------------------------------------------------*/
extern void quat2rot(const double *q,double *v)
{
    double sa,theta;

    sa=sqrt(SQR(q[1])+SQR(q[2])+SQR(q[3]));
    if (fabs(sa)<1E-6) setzero(v,3,1);
    else {
        theta=2*atan2(sa,q[0]);
        v[0]=theta*q[1]/sin(theta/2.0);
        v[1]=theta*q[2]/sin(theta/2.0);
        v[2]=theta*q[3]/sin(theta/2.0);
    }
}
extern void quatrot(const double *qab,double *va,int dir,double *vb)
{
    double v[4]={0.0,va[0],va[1],va[2]},q[4],vb_[4];
    double vr_a[4];
    int i;
    if (dir==0) {
        for (i=0;i<4;i++) q[i]=qab[i];
    }
    else {
        q[0]=qab[0];
        for (i=1;i<4;i++) q[i]=-qab[i];
    }
    quatmulx(q,v,vr_a);
    for (i=1;i<4;i++) q[i]=-qab[i];
    quatmulx(vr_a,q,vb_);
    for (i=0;i<3;i++) vb[i]=vb_[i+1];
}
/* convert rotation vector to quaternion-------------------------------------*/
extern void rv2quat(const double *rv,double *q)
{
    double rot_ang=sqrt(SQR(rv[0])+SQR(rv[1])+SQR(rv[2]));
    double cR,sR,rx,ry,rz;
    if (fabs(rot_ang)<1E-8) {
        q[0]=1.0; q[1]=q[2]=q[3]=0.0;
    }
    else {
        cR=cos(rot_ang/2.0);
        sR=sin(rot_ang/2.0);
        rx=rv[0]/rot_ang;
        ry=rv[1]/rot_ang;
        rz=rv[2]/rot_ang;

        q[0]=cR;
        q[1]=rx*sR; q[2]=ry*sR; q[3]=rz*sR;
    }
}
/* direction cosine matrix of rotation vector--------------------------------*/
extern void rot2dcm(const double *rot,double *dcm)
{
    static double I[9]={1,0,0,0,1,0,0,0,1};
    double Sr[9],Sr2[9];
    double rot_norm=norm(rot,3),sr_a,sr_b;
    int i;

    skewsym3(rot,Sr);
    matmul("NN",3,3,3,1.0,Sr,Sr,0.0,Sr2);

    if (rot_norm>0) {
        sr_a=sin(rot_norm)/rot_norm;
        sr_b=(1.0-cos(rot_norm))/SQR(rot_norm);
        for (i=0;i<9;i++) {
            dcm[i]=I[i]+sr_a*Sr[i]+sr_b*Sr2[i];
        }
    }
    else {
        seteye(dcm,3);
    }
}
/* quaternion convert to dcm-------------------------------------------------*/
extern void quat2dcmx(const double *qba, double *Cba)
{
    double a=qba[0];
    double b=qba[1];
    double c=qba[2];
    double d=qba[3];

    seteye(Cba,3);
    Cba[0]=SQR(a)+SQR(b)-SQR(c)-SQR(d);
    Cba[1]=2*(b*c+a*d);
    Cba[2]=2*(b*d-a*c);

    Cba[3]=2*(b*c-a*d);
    Cba[4]=SQR(a)-SQR(b)+SQR(c)-SQR(d);
    Cba[5]=2*(c*d+a*b);

    Cba[6]=2*(b*d+a*c);
    Cba[7]=2*(c*d-a*b);
    Cba[8]=SQR(a)-SQR(b)-SQR(c)+SQR(d);
}
/* dcm convert to quaternion-------------------------------------------------*/
extern void dcm2quatx(const double *dcm, double* quat)
{
    double tr=dcm[0]+dcm[4]+dcm[8];
    double Pa=1.0+tr;
    double Pb=1.0+2.0*dcm[0]-tr;
    double Pc=1.0+2.0*dcm[4]-tr;
    double Pd=1.0+2.0*dcm[8]-tr;
    int i;

    setzero(quat,1,4);
    if (Pa>=Pb&&Pa>=Pc&&Pa>=Pd) {
        quat[0]=0.5*sqrt(Pa);
        quat[1]=(dcm[5]-dcm[7])/4.0/quat[0];
        quat[2]=(dcm[6]-dcm[2])/4.0/quat[0];
        quat[3]=(dcm[1]-dcm[3])/4.0/quat[0];
    }
    else if (Pb>=Pc&&Pb>=Pd) {
        quat[1]=0.5*sqrt(Pb);
        quat[2]=(dcm[1]+dcm[3])/4.0/quat[1];
        quat[3]=(dcm[6]+dcm[2])/4.0/quat[1];
        quat[0]=(dcm[5]-dcm[7])/4.0/quat[1];
    }
    else if (Pc>=Pd) {
        quat[2]=0.5*sqrt(Pc);
        quat[3]=(dcm[5]+dcm[7])/4.0/quat[2];
        quat[0]=(dcm[6]-dcm[2])/4.0/quat[2];
        quat[1]=(dcm[1]+dcm[3])/4.0/quat[2];
    }
    else {
        quat[3]=0.5*sqrt(Pd);
        quat[0]=(dcm[1]-dcm[3])/4.0/quat[3];
        quat[1]=(dcm[6]+dcm[2])/4.0/quat[3];
        quat[2]=(dcm[5]+dcm[7])/4.0/quat[3];
    }
    if (quat[0]<=0) {
        for (i=0;i<4;i++) quat[i]=-quat[i];
    }
}
/* get dcm matrix of b-frame to n-frame--------------------------------------*/
static void getqbn(const insstate_t *ins, double* qbn)
{
    double pos[3],Cne[9],Cbn[9];
    ecef2pos(ins->re,pos);
    ned2xyz(pos,Cne);
    matmul3("TN",Cne,ins->Cbe,Cbn);
    dcm2quatx(Cbn,qbn);
}
/* update ins states---------------------------------------------------------*/
static void update(insstate_t *ins,const double *qbn,const double *vn,
                   const double *Cen,const double h)
{
    double Cbn[9],rn[3];
    quat2dcmx(qbn,Cbn);

    matcpy(ins->Cbn,Cbn,3,3); matcpy(ins->vn,vn,1,3);

    matmul("TN",3,3,3,1.0,Cen,Cbn,0.0,ins->Cbe);
    matmul("TN",3,1,3,1.0,Cen,vn,0.0,ins->ve);

    rn[0]=acos(Cen[6]);
    rn[1]=acos(Cen[4]);
    rn[2]=h;
    matcpy(ins->rn,rn,1,3);

    pos2ecef(rn,ins->re);
#if 1
    getaccl(ins->fb,ins->Cbe,ins->re,ins->ve,ins->ae);
#else
    matmul("TN",3,1,3,1.0,Cen,ins->an,0.0,ins->ae);
#endif
    getvn(ins,ins->vn);
}
/* backward update ins states in n-frame-------------------------------------
 * update ins states with imu measurement data in n-frame in backward
 * args   : insopt   *insopt I   ins updates options
 *          insstate_t *ins  IO  ins states
 *          imudata_t *data  I   imu measurement data
 * return : 0 (fail) or 1 (ok)
 *----------------------------------------------------------------------------*/
extern int updateinsbn(const insopt_t *insopt,insstate_t *ins,const imud_t *data)
{
    double wen_n[3],wie_n[3],g,vn[3],dv[3],da[3],dt;
    double dv1[3],dv2[3],qbn[4],w[3],rn[3];
    double qb[4],qn[4],q[4],vmid[3];
    double Cen[9],Cne[9],Cn[9],h;
    double das[3]={0},dvs[3]={0},fb[3];
    int i;

    trace(3,"updateinsbn:\n");

    trace(5,"ins(-)=\n"); traceins(5,ins);

    updinsn(ins);

    /* save precious ins states */
    savepins(ins,data);

    if ((dt=-timediff(data->time,ins->time))>MAXDT||fabs(dt)<1E-6) {

        /* update information */
        ins->dt=timediff(data->time,ins->time);
        ins->ptime=ins->time;
        ins->time =data->time;

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
    /* position/velocity in n-frame */
    ecef2pos(ins->re,rn);
    getvn(ins,vn);

    geoparam(rn,vn,wen_n,wie_n,&g);
#if SCULL_CORR
    /* rotational and sculling motion correction */
    rotscull_corr(ins,insopt,dt,dvs,das);
#endif
    for (i=0;i<3;i++) {
        dv[i]=ins->fb[i]*dt+dvs[i]; da[i]=-ins->omgb[i]*dt+das[i];
    }
    /* update velocity */
    getqbn(ins,qbn);
    quatrot(qbn,dv,0,dv1);

    for (i=0;i<3;i++) w[i]=2.0*wie_n[i]+wen_n[i];
    cross3(vn,w,dv2);
    dv2[2]+=g;

    for (i=0;i<3;i++) dv2[i]=dv2[i]*dt;
    for (i=0;i<3;i++) {
        ins->vn[i]=vn[i]-dv1[i]-dv2[i]; ins->an[i]=(dv1[i]+dv2[i])/dt;
    }
    /* update the attitude */
    rvec2quat(da,qb);
    quatmulx(qbn,qb,q);

    for (i=0;i<3;i++) w[i]=(wen_n[i]+wie_n[i])*dt;
    rvec2quat(w,qn);
    quatmulx(qn,q,qbn);

    /* position updates */
    for (i=0;i<3;i++) {
        vmid[i]=(vn[i]+ins->vn[i])/2.0;
    }
    geoparam(rn,vmid,wen_n,NULL,NULL);
    for (i=0;i<3;i++) {
        w[i]=wen_n[i]*dt;
    }
    rot2dcm(w,Cn);
    ned2xyz(rn,Cne);
    matmul("NT",3,3,3,1.0,Cn,Cne,0.0,Cen); h=rn[2]+vmid[2]*dt;

    update(ins,qbn,ins->vn,Cen,h);

    ins->dt=timediff(data->time,ins->time);
    ins->ptime=ins->time;
    ins->time =data->time;
    ins->stat =INSS_MECH;

    trace(5,"ins(+)=\n"); traceins(5,ins);
    return 1;
}
/* backward update ins states in e/n-frame-------------------------------------
 * update ins states with imu measurement data in backward
 * args   : insopt   *insopt I   ins updates options
 *          insstate_t *ins  IO  ins states
 *          imudata_t *data  I   imu measurement data
 * return : 0 (fail) or 1 (ok)
 *----------------------------------------------------------------------------*/
extern int updateinsb(const insopt_t *insopt,insstate_t *ins,const imud_t *data)
{
    trace(3,"updateinsb:\n");

#if UPD_INS_E
    return updateinsbe(insopt,ins,data);
#else
    return updateinsbn(insopt,ins,data);
#endif
}

