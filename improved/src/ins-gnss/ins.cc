/*------------------------------------------------------------------------------
* ins.cc : ins common functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *       for IMU calibration without external equipments,2014.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/09/29 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define MAXDT   60.0                /* max interval to update imu (s) */
#define RE      6378137.0           /* earth semi-major axis (WGS84) (m) */
#define FE      (1.0/298.257223563) /* earth flattening (WGS84) */
#define INSUPDPRE  1                /* inertial navigation equations precision */
#define MU      3.986004418E14      /* Earth’s gravitational constant and its WGS-84 value */
#define J2      1.082627E-3         /* Earth’s second gravitational constant */
#define E       0.0818191908425     /* WGS84 eccentricity */
#define RP      6356752.31425       /* WGS84 Polar radius in meters */
#define FLAT    1.0/298.257223563   /* WGS84 flattening */
#define E_SQR   0.00669437999014    /* sqr of linear eccentricity of the ellipsoid */
#define SCULL_CORR 1                /* rotational and sculling motion correction */

/* global variable -----------------------------------------------------------*/
extern const double Omge[9]={0,OMGE,0,-OMGE,0,0,0,0,0}; /* (5.18) */

/* multiply 3d matries -------------------------------------------------------*/
extern void matmul3(const char *tr, const double *A, const double *B, double *C)
{
    matmul(tr,3,3,3,1.0,A,B,0.0,C);
}
/* multiply 3d matrix and vector ---------------------------------------------*/
extern void matmul3v(const char *tr, const double *A, const double *b, double *c)
{
    char t[]="NN";
    t[0]=tr[0];
    matmul(t,3,1,3,1.0,A,b,0.0,c);
}
/* 3d skew symmetric matrix --------------------------------------------------*/
extern void skewsym3(const double *ang, double *C)
{
    C[0]=0.0;     C[3]=-ang[2]; C[6]=ang[1];
    C[1]=ang[2];  C[4]=0.0;     C[7]=-ang[0];
    C[2]=-ang[1]; C[5]=ang[0];  C[8]=0.0;
}
/* 3d skew symmetric matrix --------------------------------------------------*/
extern void skewsym3x(double x,double y,double z,double *C)
{
    C[0]=0.0; C[3]=-z;  C[6]=y;
    C[1]=z;   C[4]=0.0; C[7]=-x;
    C[2]=-y;  C[5]=x;   C[8]=0.0;
}
/* set all matrix elements to zero ------------------------------------------*/
extern void setzero(double *A,int n,int m)
{
    int i,j;
    for (i=0;i<n;i++) for (j=0;j<m;j++) A[i*m+j]=0.0;
}
/* D=A*B*C---------------------------------------------------------------------*/
extern void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D)
{
    char tr_[8];
    double *T=mat(n,q);
    matmul(tr,n,q,p,1.0,A,B,0.0,T);
    sprintf(tr_,"N%c",tr[2]);
    matmul(tr_,n,m,q,1.0,T,C,0.0,D); free(T);
}
/* conversion matrix of ned frame to ecef frame--------------------------------
 * conversion matrix of ned to ecef
 * args   : double *pos     I   position {lat,lon,height} (rad/m)
 *          double *Cne     O   convertion matrix between frame
 * return : none
 * --------------------------------------------------------------------------*/
extern void ned2xyz(const double *pos,double *Cne)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);

    Cne[0]=-sinp*cosl; Cne[3]=-sinl; Cne[6]=-cosp*cosl;
    Cne[1]=-sinp*sinl; Cne[4]= cosl; Cne[7]=-cosp*sinl;
    Cne[2]=      cosp; Cne[5]= 0.0;  Cne[8]=-sinp;
}
/* ned frame to body frame---------------------------------------------------
 * args    : double rpy      I  attitude {roll,picth,yaw} (rad)
 *           double Cnb      O  matrix of ned to body frame
 * return  : none
 * -------------------------------------------------------------------------*/
extern void rpy2dcm(const double *rpy,double *Cnb)
{
    double sin_phi=sin(rpy[0]),cos_phi=cos(rpy[0]),
           sin_theta=sin(rpy[1]),cos_theta=cos(rpy[1]),
           sin_psi=sin(rpy[2]),cos_psi=cos(rpy[2]);
    Cnb[0]= cos_theta*cos_psi;
    Cnb[3]= cos_theta*sin_psi;
    Cnb[6]=-sin_theta;
    Cnb[1]=-cos_phi*sin_psi+sin_phi*sin_theta*cos_psi;
    Cnb[4]= cos_phi*cos_psi+sin_phi*sin_theta*sin_psi;
    Cnb[7]= sin_phi*cos_theta;
    Cnb[2]= sin_phi*sin_psi+cos_phi*sin_theta*cos_psi;
    Cnb[5]=-sin_phi*cos_psi+cos_phi*sin_theta*sin_psi;
    Cnb[8]= cos_phi*cos_theta;
}
/* convert transform matrix to quaternion----------------------------------
 * args    : double *C  I  transform matrix
 *           quat_t *q  O  quaternion
 * return  : none
 * -----------------------------------------------------------------------*/
extern void dcm2quat(const double *C,quat_t *q)
{
    double a[3][3]={0},s,tr;

    a[0][0]=C[0]; a[0][1]=C[3]; a[0][2]=C[6];
    a[1][0]=C[1]; a[1][1]=C[4]; a[1][2]=C[7];
    a[2][0]=C[2]; a[2][1]=C[5]; a[2][2]=C[8];

    tr=a[0][0]+a[1][1]+a[2][2];
    if (tr>0) {
        s=0.5/sqrt(tr+1.0);
        q->w=0.25/s;
        q->x=(a[2][1]-a[1][2])*s;
        q->y=(a[0][2]-a[2][0])*s;
        q->z=(a[1][0]-a[0][1])*s;
    }
    else {
        if (a[0][0]>a[1][1]&&a[0][0]>a[2][2]) {
            s=2.0*sqrt(1.0+a[0][0]-a[1][1]-a[2][2]);
            q->w=(a[2][1]-a[1][2])/s;
            q->x=0.25*s;
            q->y=(a[0][1]+a[1][0])/s;
            q->z=(a[0][2]+a[2][0])/s;
        }
        else if (a[1][1]>a[2][2]) {
            s=2.0*sqrt(1.0+a[1][1]-a[0][0]-a[2][2]);
            q->w=(a[0][2]-a[2][0])/s;
            q->x=(a[0][1]+a[1][0])/s;
            q->y=0.25*s;
            q->z=(a[1][2]+a[2][1])/s;
        }
        else {
            s=2.0*sqrt(1.0+a[2][2]-a[0][0]-a[1][1]);
            q->w=(a[1][0]-a[0][1])/s;
            q->x=(a[0][2]+a[2][0])/s;
            q->y=(a[1][2]+a[2][1])/s;
            q->z=0.25*s;
        }
    }
}
/* conversion quaternion to matrix ---------------------------------------
 * args    : quat_t *q    I  input quaternion
 *           double *C    O  output transform matrix
 * return  : none
 * ----------------------------------------------------------------------*/
extern void quat2dcm(const quat_t *q,double *C)
{
    double sqw=q->w*q->w,sqx=q->x*q->x,sqy=q->y*q->y,sqz=q->z*q->z;
    double invs,tmp1,tmp2;

    invs=1.0/(sqx+sqy+sqz+sqw);
    C[0]=( sqx-sqy-sqz+sqw)*invs;
    C[4]=(-sqx+sqy-sqz+sqw)*invs;
    C[8]=(-sqx-sqy+sqz+sqw)*invs;

    tmp1=q->x*q->y; tmp2=q->z*q->w;
    C[1]=2.0*(tmp1+tmp2)*invs;
    C[3]=2.0*(tmp1-tmp2)*invs;

    tmp1=q->x*q->z; tmp2=q->y*q->w;
    C[2]=2.0*(tmp1-tmp2)*invs;
    C[6]=2.0*(tmp1+tmp2)*invs;

    tmp1=q->y*q->z; tmp2=q->x*q->w;
    C[5]=2.0*(tmp1+tmp2)*invs;
    C[7]=2.0*(tmp1-tmp2)*invs;
}
/* right-forward-up convert to forward-right-down--------------------------
 * args    : double *rfu       I  right-forward-up coordinates
 *           double *frd       O  forward-right-down coordinates
 *           double *C         O  transform matrix (NULL: no output)
 * return : none
 * -----------------------------------------------------------------------*/
extern void rfu2frd(const double *rfu,double *frd,double *C)
{
    if (rfu&&frd) {
        frd[0]=rfu[1]; frd[1]=rfu[0]; frd[2]=-rfu[2];
    }
    if (C) setzero(C,3,3),C[1]=C[3]=1.0,C[8]=-1.0;
}
/* enu-frame convert to ned-frame ----------------------------------------
 * args    : double *enu       I  enu coordinates
 *           double *ned       O  ned coordinates
 *           double *C         O  transform matrix (NULL: no output)
 * return  : none
 * ----------------------------------------------------------------------*/
extern void enu2ned(const double *enu,double *ned,double *C)
{
    if (enu&&ned) {
        ned[0]=enu[1]; ned[1]=enu[0]; ned[2]=-enu[2];
    }
    if (C) setzero(C,3,3),C[1]=C[3]=1.0,C[8]=-1.0;
}
/* ins error model (ref [2])-----------------------------------------------
 * args:   : double *accl      I  force measurements in b-frame (m/s^2)
 *           double *cor_accl  O  corrected force measurements in b-frame
 *           double *gyro      I  angular rate measurements in b-frame (rad/s)
 *           double *cor_gyro  O  corrected angular rate measurements in
 *                                b-frame (rad/s)
 *           insstate_t *ins   I  ins updates states
 * return  :none
 * -----------------------------------------------------------------------*/
extern void ins_errmodel(const double *accl,const double *gyro,double *cor_accl,
                         double *cor_gyro,insstate_t* ins)
{
    int i,j;
    double Mai[9],Mgi[9],*I=eye(3),T[9]={0},Gf[3]={0};

    for (i=0;i<3;i++) for (j=0;j<3;j++) Mai[i+j*3]=I[i+j*3]+ins->Ma[i+j*3];
    for (i=0;i<3;i++) for (j=0;j<3;j++) Mgi[i+j*3]=I[i+j*3]+ins->Mg[i+j*3];

    if (!matinv(Mai,3)&&!matinv(Mgi,3)) {
        matmul("NN",3,1,3,1.0,Mai,accl,0.0,T);
        matmul("NN",3,1,3,1.0,Mgi,gyro,0.0,T+3);
    }
    if (cor_accl) {
        for (i=0;i<3;i++) cor_accl[i]=T[i]-ins->ba[i];
        matmul("NN",3,1,3,1.0,ins->Gg,accl,0.0,Gf);
    }
    if (cor_gyro) {
        for (i=0;i<3;i++) {
            cor_gyro[i]=T[3+i]-ins->bg[i]-Gf[i];
        }
    }
    free(I);
}
/* correction imu accl. and gyro. measurements-----------------------------
 * args     :    I  double *accl  imu accl. measurements
 *               I  double *gyro  imu gyro. measurements
 *               I  double *Ma    non-orthogonal between sensor axes and
 *                                body frame for accelerometers
 *               I  double *Mg    non-orthogonal between sensor axes and
 *                                body frame for gyroscopes
 *               I  double *ba    accelemeter-bias
 *               I  double *bg    gyro-bias
 *               I  double *Gg    g-dependent bias for a gyro triad
 *               O  double *cor_accl
 *                         *cor_gyro  corrected imu accl. and gyro. measurements
 * return : none
 * -----------------------------------------------------------------------*/
extern void ins_errmodel2(const double *accl,const double *gyro,const double *Ma,
                          const double *Mg,const double *ba,const double *bg,
                          const double *Gg,double *cor_accl,double *cor_gyro)
{
    int i,j;
    double Mai[9],Mgi[9],*I=eye(3),T[9]={0},Gf[3]={0};

    for (i=0;i<3;i++) for (j=0;j<3;j++) Mai[i+j*3]=I[i+j*3]+Ma[i+j*3];
    for (i=0;i<3;i++) for (j=0;j<3;j++) Mgi[i+j*3]=I[i+j*3]+Mg[i+j*3];

    if (!matinv(Mai,3)&&!matinv(Mgi,3)) {
        matmul("NN",3,1,3,1.0,Mai,accl,0.0,T);
        matmul("NN",3,1,3,1.0,Mgi,gyro,0.0,T+3);
    }
    if (cor_accl) {
        for (i=0;i<3;i++) cor_accl[i]=T[i]-ba[i];
        matmul("NN",3,1,3,1.0,Gg,accl,0.0,Gf);
    }
    if (cor_gyro) {
        for (i=0;i<3;i++) cor_gyro[i]=T[3+i]-bg[i]-Gf[i];
    }
    free(I);
}
/* converts a coordinate transformation matrix to the corresponding set of
 * Euler angles -------------------------------------------------------------
 * args    : double *Cnb      I matrix of ned-frame to body-frame
 *           double *rpy      O attitude {roll,pitch,yaw}
 * return  : none
 * ------------------------------------------------------------------------*/
extern void dcm2rpy(const double *Cnb,double *rpy)
{
    rpy[0]= atan2(Cnb[7],Cnb[8]);  /* roll */
    rpy[1]=-asin (Cnb[6]);         /* pitch */
    rpy[2]= atan2(Cnb[3],Cnb[0]);  /* yaw */
}
/* converts a rotate vector to transformation matrix-----------------------
 * args    : double *rv       I rotate vector
 *           double *C        O transformation matrix
 * return  : none
 * -----------------------------------------------------------------------*/
extern void rov2dcm(const double *rv,double *C)
{
    int i;
    double a,*I=eye(3),rs[9],rs2[9],n1,n2;

    if ((a=norm(rv,3))<1E-8) {
        n1=1.0-SQR(a)*(1.0/6.0 -SQR(a)/120.0);
        n2=0.5-SQR(a)*(1.0/24.0-SQR(a)/720.0);
    }
    else n1=sin(a)/a,n2=(1.0-cos(a))/SQR(a);

    skewsym3(rv,rs);
    matmul("NN",3,3,3,1.0,rs,rs,0.0,rs2);
    for (i=0;i<9;i++) C[i]=I[i]-n1*rs[i]+n2*rs2[i];
    free(I);
}
/* convert rotation vector to transformation quaternion-------------------
 * args  : double *rv       I rotate vector
 *         quat_t *q        O quaternion
 * return : none
 * note   : this function is the same as quat_init_axis(),see quaternion.c
 * -----------------------------------------------------------------------*/
extern void rov2qua(const double *rv,quat_t *q)
{
    int i;
    double a,s;

    if ((a=norm(rv,3))<1E-8) {
        q->w=1.0-SQR(a)*(1.0/8.0-SQR(a)/384.0);
        s=0.5-SQR(a)*(1.0/48.0-SQR(a)/3840.0);
    }
    else q->w=cos(a/2.0),s=sin(a/2.0)/a;
    for (i=1;i<4;i++) q->vec[i]=s*rv[i-1];
}
/* 3x1 vector coordinate transformation by quaternion------------------------
 * args  : double *vi       I  vector to be transformed
 *         quat_t *q        I  transformation quaternion
 *         double *vo       I  output vector
 * return ：none
 * --------------------------------------------------------------------------*/
extern void qmulv(const double *vi,const quat_t *quat,double *vo)
{
    double C[9];
    quat_to_rh_rot_matrix(quat,C);
    matmul("NN",3,1,3,1.0,C,vi,0.0,vo);
}
/* convert quaternion to euler -----------------------------------------------
 * args  : quat_t *quat    I  transformation quaternion
 *         double *rpy     O  output vector (roll,pitch,yaw)
 * return : none
 * --------------------------------------------------------------------------*/
extern void quat2rpy(const quat_t *quat,double *rpy)
{
    euler_t euler;
    quat_to_euler(&euler,quat);
    rpy[0]=euler.roll; rpy[1]=euler.pitch; rpy[2]=euler.yaw;
}
/* attitude quaternion updating using rotation vector -----------------------
 * args  :  double  *vi  I  vector to be transformed
 *          quat_t  *qo  IO input and output quaternion,
 * return : none
 * -------------------------------------------------------------------------*/
extern void quatupd(const double *vi,quat_t *qo)
{
    quat_t q1,q2; rov2qua(vi,&q1); quat_copy(&q2,qo); quat_mul(qo,&q1,&q2);
}
/* calculates  acceleration due to gravity resolved about north, east, and down
 * args  : double *pos      I  position in n-frame {lat,lon,h} (rad/m)
 *         double *gn       O  gravity in n-frame (m/s2)
 * return: none
 * --------------------------------------------------------------------------*/
extern void gravity_ned(const double *pos,double *gn)
{
    double g0;
    g0=9.7803253359*(1.0+0.001931853*SQR(sin(pos[0])))/sqrt(1.0-SQR(E*sin(pos[0]))); /* (2.85) */
    /* calculate north gravity */
    gn[0]=-8.08E-9*pos[2]*sin(2.0*pos[0]); /* (2.140)-v2.0 */
    /* east gravity is zero */
    gn[1]=0.0;
    /* calculate down gravity */
    gn[2]=g0*(1.0-2.0/RE*(1.0+FLAT*(1.0-2*SQR(sin(pos[0])))+
          SQR(OMGE*RE)*RP/MU)*pos[2]+3.0/SQR(RE)*SQR(pos[2])); /* (2.139)-v2.0 */
}
/* gravity at earth surface ----------------------------------------------------
* gravity at earth surface
* args   : double *pos      I   position {lat,lon,height} (rad/m)
* return : gravity acceleration at earth surface (m/s^2)
*-----------------------------------------------------------------------------*/
extern double gravity0(const double *pos)
{
    double sinp,e2=FE*(2.0-FE);
    sinp=sin(pos[0]); /* (2.85) */
    return 9.7803253359*(1.0+0.001931853*sinp*sinp)/sqrt(1-e2*sinp*sinp);
}
/* gravity model ---------------------------------------------------------------
* gravity model in e-frame
* args   : double *re       I   position (ecef) (m)
*          double *ge       O   gravity acceleration (ecef) (m/s^2)
* return : none
*-----------------------------------------------------------------------------*/
extern void gravity(const double *re, double *ge)
{
    double pos[3],gn[3]={0},Cne[9];
    ecef2pos(re,pos);
    gn[2]=gravity0(pos)*(1.0-2.0*pos[2]/RE);
    ned2xyz(pos,Cne); 
    matmul3v("N",Cne,gn,ge);
}
/* calculates  acceleration due to gravity resolved about ecef-frame ---------
 * args  : double *pos     I   cartesian position of body frame w.r.t. ecef frame,
 *                             resolved about ecef-frame axes (m)
 *         double *g       O   Acceleration due to gravity (m/s^2) in ecef-frame
 * return: none
 * --------------------------------------------------------------------------*/
extern void pregrav(const double *pos, double *g)
{
    double r,z,gamma;

    /* calculate distance from center of the earth */
    if ((r=norm(pos,3))<RE_WGS84/2.0) {
        g[0]=g[1]=0.0; g[2]=-9.81;
        return;
    }
    z=-MU/(r*r*r);
    gamma=1.5*J2*SQR(RE)/SQR(r);
    g[0]=z*(pos[0]+gamma*(1.0-5.0*SQR(pos[2]/r))*pos[0]);
    g[1]=z*(pos[1]+gamma*(1.0-5.0*SQR(pos[2]/r))*pos[1]);
    g[2]=z*(pos[2]+gamma*(3.0-5.0*SQR(pos[2]/r))*pos[2]); /* (2.91) */
    g[0]+=SQR(OMGE)*pos[0];
    g[1]+=SQR(OMGE)*pos[1]; /* (2.84) */
}
/* gives the heading Euler angle, in terms of the roll, pitch,
 * and gyro measurements ------------------------------------------------------
 * args  : double roll     I   roll of body-to-ned frame {rad}    I
 *         double pitch    I   pitch of body-to-ned frame {rad}
 *         double *gyro    I   angular rate measurements in b-frame (rad/s)
 *         double *head    O   head of body-to-ned frame {rad}
 * return: none
 * --------------------------------------------------------------------------*/
extern void rp2head(const double roll,const double pitch,const double *gyro,
                    double *head)
{
    double sinh=-gyro[1]*cos(roll )+gyro[2]*sin(roll);
    double cosh= gyro[0]*cos(pitch)+gyro[1]*sin(roll)*sin(pitch)+
                 gyro[2]*cos(roll )*sin(pitch);
    *head=atan2(sinh,cosh); /* (5.91) */
}
/* initialize ins states ------------------------------------------------------
* initialize ins states with stationaly imu measurement data
* args   : insstate_t *ins  IO  ins states
*          double *re       I   initial ins position (ecef) (m)
*          double angh      I   initial heading angle (rad)
*          imud_t *data     I   imu measurement data (static case)
*          int    n         I   number of data (0: no use of data)
* return : none
*----------------------------------------------------------------------------*/
extern void initins(insstate_t *ins, const double *re, double angh,
                    const imud_t *data, int n,const insopt_t *opt)
{
    gtime_t time={0};
    double rpy[3]={0},ab[3]={0},fb[3]={0},pos[3],Cnb[9],Cne[9],ge[3],gb[3];
    int i,j;

    trace(3,"initins: n=%d angh=%.1f\n",n,angh*R2D);
    ecef2pos(re,pos);
    rpy[2]=angh;

    /* initial ned position/velocity/acceleration */
    for (i=0;i<3;i++) {
        ins->rn[i]=pos[i];
        ins->vn[i]=ins->an[i]=ins->ba[i]=ins->bg[i]=ins->fb[i]=0.0;
    }
    /* initial ins position/velocity/acceleration */
    for (i=0;i<3;i++) {
        ins->re[i]=re[i];
        ins->ve[i]=ins->ae[i]=ins->ba[i]=ins->bg[i]=ins->fb[i]=0.0;
    }
    if (n>0) { /* gyro and accl bias */
        for (i=0;i<3;i++) {
            for (j=0;j<n;j++) {
                fb[i]+=data[j].accl[i];
                ab[i]+=data[j].gyro[i];
            }
            fb[i]/=n; ab[i]/=n;
        }
        rpy[0]=atan2(-fb[1],-fb[2]); /* (5.89) */
        rpy[1]=atan(fb[0]/norm(fb+1,2));
        
        if (angh==0.0) {
            rp2head(rpy[0],rpy[1],ab,rpy+3);
        }
        time=data[n-1].time;
    }
    /* initial body-frame to n-frame */
    rpy2dcm(rpy,Cnb); matt(Cnb,3,3,ins->Cbn);
    ned2xyz(pos,Cne);
    matmul3("NT",Cne,Cnb,ins->Cbe); /* conversion order */
    if (n>0) {
        gravity(re,ge);
        matmul3v("T",ins->Cbe,ge,gb);
        for (i=0;i<3;i++) {
            ins->ba[i]=fb[i]+gb[i];
            ins->bg[i]=ab[i];
        }
    }
    ins->time=time;
}
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
static void updateatt(double t, double *Cbe, const double *omgb,const double *das)
{
    double alpha[3],a,a1,a2,Ca[9],Ca2[9],Comg[9],Cbep[9];
    double Cbb[9]={1,0,0,0,1,0,0,0,1},Cei[9]={0};
    int i;
    trace(3,"updateatt:\n");
    for (i=0;i<3;i++) alpha[i]=omgb[i]*t+das[i];
    skewsym3(alpha,Ca);
    matmul3("NN",Ca,Ca,Ca2);
    a=norm(alpha,3);
    if (a<1E-8) {
        a1=1.0-a*a/6.0; a2=0.5-a*a/24.0; /* 5.62 */
    }
    else {
        a1=sin(a)/a; a2=(1.0-cos(a))/(a*a);
    }
    for (i=0;i<9;i++) Cbb[i]+=a1*Ca[i]+a2*Ca2[i]; /* (5.63) */

#if INSUPDPRE
    /* determine the earth rotation over the update interval */
    Cei[0]= cos(OMGE*t); Cei[3]=sin(OMGE*t);
    Cei[1]=-sin(OMGE*t); Cei[4]=cos(OMGE*t); Cei[8]=1.0;
    matmul3("NN",Cei,Cbe,Cbep);
    matmul3("NN",Cbep,Cbb,Cbe);
#else
    matmul3("NN",Cbe,Cbb,Cbep);
    matmul3("NN",Omge,Cbe,Comg); /* (5.65) */
    for (i=0;i<9;i++) Cbe[i]=Cbep[i]-Comg[i]*t; /* (5.20) */
#endif
}
/* update ins states ----------------------------------------------------------
* updata ins states with imu measurement data in e-frame
* args   : insopt   *insopt I   ins updates options
*          insstate_t *ins  IO  ins states
*          imudata_t *data  I   imu measurement data
* return : 0 (fail) or 1 (ok)
*----------------------------------------------------------------------------*/
extern int updateins(const insopt_t *insopt,insstate_t *ins,const imud_t *data)
{
    double dt,Cbe[9],fe[3],ge[3],cori[3],Cbb[9]={1,0,0,0,1,0,0,0,1};
    double Ca[9],Ca2[9],a1,a2,a,alpha[3]={0},Omg[9]={0},ae[3]={0};
    double das[3]={0},dvs[3]={0},fb[3];
    int i;

    trace(3,"updateins:\n");

    trace(5,"ins(-)=\n"); traceins(5,ins);

    /* save precious epoch ins states */
    savepins(ins,data);

    if ((dt=timediff(data->time,ins->time))>MAXDT||fabs(dt)<1E-6) {

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
    ae[2]=OMGE*dt;

#if SCULL_CORR
    /* rotational and sculling motion correction */
    rotscull_corr(ins,insopt,dt,dvs,das);
#endif
    /* update attitude */
    for (i=0;i<9;i++) Cbe[i]=ins->Cbe[i];
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

        /* precious gravity model */
        pregrav(ins->re,ge);
    }
    else gravity(ins->re,ge); /* common gravity model */

    /* update velocity/position */
    matmul3v("N",Omge,ins->ve,cori);
    for (i=0;i<3;i++) {
        ins->ae[i]=fe[i]+ge[i]-2.0*cori[i]; 
        ins->ve[i]+=ins->ae[i]*dt;
        ins->re[i]+=ins->ve[i]*dt+ins->ae[i]/2.0*dt*dt;
    }
    /* update ins state in n-frame */
    updinsn(ins);

    ins->dt=timediff(data->time,ins->time);
    ins->ptime=ins->time;
    ins->time =data->time;
    ins->stat =INSS_MECH;  

    trace(5,"ins(+)=\n"); traceins(5,ins);
    return 1;
}
/* Calculates the meridian and transverse radii of curvature------------------
 * args   : double *rn        I    position in n-frame (lat,lon,h) {rad/m}
 *          double *R_N       O    meridian radius of curvature (m)
 *          double *R_E       O    transverse radius of curvature (m)
 * return : none
 * --------------------------------------------------------------------------*/
extern void radii(const double *rn,double *R_N,double *R_E)
{
    *R_N=RE*(1.0-SQR(E))/pow(1.0-SQR(E)*SQR(sin(rn[0])),1.5);
    *R_E=RE/sqrt(1.0-SQR(E)*SQR(sin(rn[0])));
}
/* set matrix to eye-matrix--------------------------------------------------*/
extern void seteye(double* A,int n)
{
    int i,j; for (i=0;i<n;i++) for (j=0;j<n;j++) A[i+j*n]=(i==j?1.0:0.0);
}
/* runs precision local-navigation-frame inertial navigation equations---------
 * args  : insopt_t *insopt   I    ins updates options
 *         insstate_t *ins    I O  ins updates states
 *         imudata *data      I    imu measurement data
 * return : 0 (failed) or 1 (ok)
 * note   : only the attitude update and specific force
 *          frame transformation phases are precise
 * --------------------------------------------------------------------------*/
extern int insupdate_ned(const insopt_t *insopt,insstate_t *ins,const imud_t *data)
{
    return updateinsn(insopt,ins,data);
}
/* remove effects of lever-arm------------------------------------------------
 * args  :         double *pos      I  imu body position in ecef-frame
 *                 double *re       I  imu body position in ecef-frame
 *                 double *ve       I  imu body velecity in ecef-frame
 *                 double *lever    I  imu body-frame to gnss antenna lever-arm
 *                 double *omgb     I  imu gyro measurements data
 *                 double *rec      O  imu body position in ecef-frame after
 *                                     removing effects of lever-arm (NULL: no output)
 *                 double *vec      O  imu body velecity in ecef-frame after
 *                                     removing effects of lever-arm (NULL: no output)
 * return : none
 * note: this function can convert ins body position/velocity to
 *       gps antenna position/velocity
 * --------------------------------------------------------------------------*/
extern void rmlever(const double *pos, const double *re, const double *ve,
                    const double *lever, const double *Cbe, const double *omgb,
                    double *rec, double *vec)
{
    int i; double wl[3],T[3],Omg[9];

    trace(3,"rmlever :\n");

    /* correct position */
    matmul3v("N",Cbe,lever,T);
    if (rec) for (i=0;i<3;i++) rec[i]=re[i]+T[i];

    /* correct velecity */
    skewsym3(omgb,Omg);
    matmul3v("N",Omg,lever,wl);
    matmul3v("N",Cbe,wl,T);
    matmul33("NNN",Omge,Cbe,lever,3,3,3,1,wl);
    if (vec) for (i=0;i<3;i++) vec[i]=ve[i]+T[i]-wl[i];
}
/* coning and sculling compensation------------------------------------------
 * args    :  double *gyro      I  gyro angular increments
 *            double *accl      I  acc velocity increments (may not exist)
 *            double *phim      O  rotation vector after coning compensation
 *            double *dvbm      O  velocity increment after rotation and
 *                                 sculling compensation
 *            double *wm0       I  precious gyro angular increments
 *            double *vm0       I  precious acc velocity increments
 *            double *dphim     O  attitude coning error
 *            double *rotm      O  velocity rotation error
 *            double *scullm    O  velocity sculling error
 * return :none---------------------------------------------------------------*/
extern void cnscl(const double *gyro,const double *accl,double *phim,double *dvbm,
                  const double *wm0,const double *vm0,double *dphim,double *rotm,
                  double *scullm)
{
    int i; double tmp[3],tmp1[3];

    cross3(wm0,gyro,tmp);
    for (i=0;i<3;i++) dphim[i]=tmp[i]*1.0/12.0;

    for (i=0;i<3;i++) phim[i]=gyro[i]+dphim[i];

    cross3(gyro,accl,tmp); cross3(vm0,gyro,tmp1);
    for (i=0;i<3;i++) scullm[i]=(tmp[i]+tmp1[i])/12.0;

    cross3(gyro,accl,tmp);
    for (i=0;i<3;i++) rotm[i]=1.0/2.0*tmp[i];
    for (i=0;i<3;i++) dvbm[i]=accl[i]+rotm[i]+scullm[i];
}
/* trace ins states ----------------------------------------------------------*/
extern void traceins(int level, const insstate_t *ins)
{
#if TRACE_INS
    double pos[3],Cne[9],Cnb[9],rpy[3],vel[3],acc[3];
    char s[64];
    time2str(ins->time,s,3);
    trace(level,"time  =%s\n",s);
    ecef2pos(ins->re,pos);
    ned2xyz(pos,Cne);
    matmul3("TN",ins->Cbe,Cne,Cnb); 
    dcm2rpy(Cnb,rpy);
    trace(level,"attn =%8.5f %8.5f %8.5f\n",rpy[0]*R2D,rpy[1]*R2D,rpy[2]*R2D);

    matmul3v("T",Cne,ins->ve,vel);
    matmul3v("T",Cne,ins->ae,acc);
    trace(level,"veln =%8.5f %8.5f %8.5f accn =%8.5f %8.5f %8.5f\n",
          vel[0],vel[1],vel[2],acc[0],acc[1],acc[2]); /* n-frame */

    matmul3v("T",ins->Cbe,ins->ve,vel);
    matmul3v("T",ins->Cbe,ins->ae,acc);
    trace(level,"velb =%8.5f %8.5f %8.5f accb =%8.5f %8.5f %8.5f\n",
          vel[0],vel[1],vel[2],acc[0],acc[1],acc[2]); /* b-frame */

    trace(level,"gbias =%8.5f %8.5f %8.5f abias=%8.5f %8.5f %8.5f\n",
          ins->bg[0],ins->bg[1],ins->bg[2],ins->ba[0],ins->ba[1],ins->ba[2]);
#endif
}
/* read imu measurement log file -----------------------------------------------
 * read imu measurement log file
 * args   : char   *file     I   imu measurement log file
 *          imu_t  *imu      O   imu measurement data in frd-frame
 *          int    decfmt    I   imu measurement data decode format
 *          int    format    I   imu measurement data format
 *          int    coor      I   imu body coordinate frame
 *          int    valfmt    I  imu gyro measurement data format
 * return : status (1:ok,0:no data/error)
 *-----------------------------------------------------------------------------*/
extern int readimu(const char *file, imu_t *imu,int decfmt,int format,int coor,
                   int valfmt)
{
    FILE *fp;
    imud_t *p;
    double v[12],sow;
    size_t siz;
    int week;
    char buff[1024];
    gtime_t t0;

    trace(3,"readimulog:s=%s\n",file);
    imu->n=imu->nmax=0; imu->data=NULL;
    imu->format=format; imu->coor=coor;
    imu->decfmt=decfmt; imu->valfmt=valfmt;

    if (!(fp=fopen(file,"r"))) {
        fprintf(stderr,"file open error : %s\n",file);
        return 0;
    }
    t0=timeget(); sow=time2gpst(t0,&week);

    while (fgets(buff,sizeof(buff),fp)) {
        /* time gyrox gyroy gyroz accx accy accz odometry */
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf \n",
                   v,v+1,v+2,v+3,v+4,v+5,v+6,v+7)<8) {
            continue;
        }
        if (imu->n>=imu->nmax) {
            trace(5,"readimulog:imu->n=%d nmax=%d\n",imu->n,imu->nmax);
            siz=sizeof(imud_t)*(imu->nmax+=4096);
            if (!(imu->data=(imud_t *)realloc(imu->data,siz))) {
                free(imu->data); imu->n=imu->nmax=0;
                fprintf(stderr,"memory allocation error\n");
                break;
            }
        }
        p=imu->data+imu->n++;

        p->gyro[0]=v[1]; /* rad */
        p->gyro[1]=v[2];
        p->gyro[2]=v[3];
        p->accl[0]=v[4]; /* m/s */
        p->accl[1]=v[5];
        p->accl[2]=v[6];

        /* time record */
        sow=v[0]; p->time=gpst2time(week,sow);
    }
    fclose(fp);
    return imu->n<=0?0:1;
}
/* use imu stationaly imu measurement data to estimate attitude--------------
 * args  :  imud_t *data  I  stationaly imu measurement data
 *          int n         I  number of imu measurement data
 *          double *Cbn   O  transform matrix
 * return : none
 * -------------------------------------------------------------------------*/
extern void estatt(const imud_t *data,int n,double *Cbn)
{
    int i,j;
    double fb[3]={0},ab[3]={0},rpy[3]={0},Cnb[9];

    trace(3,"estatt :\n");

    if (n>0) { /* gyro and accl bias */
        for (i=0;i<3;i++) {
            for (j=0;j<n;j++) {
                fb[i]+=data[j].accl[i];
                ab[i]+=data[j].gyro[i];
            }
            fb[i]/=n; ab[i]/=n;
        }
        rpy[0]=atan2(-fb[1],-fb[2]); /* (5.89) */
        rpy[1]=atan(fb[0]/norm(fb+1,2));
        rp2head(rpy[0],rpy[1],ab,rpy+3);
    }
    rpy2dcm(rpy,Cnb); matt(Cnb,3,3,Cbn);
}
/* gnss antenna position/velecity transform to ins position/velecity---------
 * args  :  double *pos    I  position of gnss antenna (ecef)
 *          double *vel    I  velecity of gnss antenna (ecef)
 *          double *Cbe    I  transform matrix of body-frame to ecef-frame
 *          double *lever  I  arm lever
 *          imud_t *imu    I  imu measurement data
 *          double *posi   O  ins position (ecef) (NULL :no output)
 *          double *veli   O  ins velecity (ecef) (NULL :no output)
 * return : none
 * -------------------------------------------------------------------------*/
extern void gapv2ipv(const double *pos,const double *vel,const double *Cbe,
                     const double *lever,const imud_t *imu,double *posi,
                     double *veli)
{
    int i; double T[9],TT[9];

    if (posi) {
        matcpy(posi,pos,3,1);
        matmul("NN",3,1,3,-1.0,Cbe,lever,1.0,posi);
    }
    if (veli) {
        skewsym3(imu->gyro,T);
        matmul("NN",3,1,3,1.0,T,lever,0.0,TT);
        matmul("NN",3,1,3,1.0,Cbe,TT,0.0,T);

        matmul33("NNN",Omge,Cbe,lever,3,3,3,1,TT);
        for (i=0;i<3;i++) veli[i]=vel[i]-T[i]+TT[i];
    }
}
/* correction direction cosine matrix by attitude errors---------------------
 * args   :  double *dx  I   attitude errors
 *           double *C   IO  direction cosine matrix
 * return : none
 * --------------------------------------------------------------------------*/
extern void corratt(const double *dx,double *C)
{
    int i;
    double T[9],*I=eye(3);

    skewsym3(dx,T);
    for (i=0;i<9;i++) I[i]-=T[i];

    matcpy(T,C,3,3);
    matmul("NN",3,3,3,1.0,I,T,0.0,C);
    free(I);
}
/* get attitude from ins states----------------------------------------------
 * args   :  insstate_t *ins  I  ins states
 *           double *rpy      O  ins attitude
 * return : none
 * --------------------------------------------------------------------------*/
extern void getatt(const insstate_t *ins,double *rpy)
{
    double llh[3],C[9],Cnb[9];

    ecef2pos(ins->re,llh);
    ned2xyz(llh,C);
    matmul("TN",3,3,3,1.0,ins->Cbe,C,0.0,Cnb);
    dcm2rpy(Cnb,rpy);
}
/* adjust imu data to frd-ned frame and convert to angular rate/acceleration
 * args   :  prcopt_t *opt  I   options
 *           imud_t *imu    IO  imu measurement data
 * return : none
 * --------------------------------------------------------------------------*/
extern void adjustimu(const prcopt_t *opt,imud_t *imu)
{
    double gyro[3],accl[3],dt=1.0/opt->insopt.hz;
    int i;

    trace(3,"adjustimu:\n");

    if (opt->insopt.imucoors==IMUCOOR_RFU) { /* convert to frd-ned-frame */
        matcpy(gyro,imu->gyro,1,3);
        matcpy(accl,imu->accl,1,3);
        matmul("NN",3,1,3,1.0,Crf,gyro,0.0,imu->gyro);
        matmul("NN",3,1,3,1.0,Crf,accl,0.0,imu->accl);
    }
    if (opt->insopt.imudecfmt==IMUDECFMT_INCR) {
        for (i=0;i<3;i++) {
            imu->gyro[i]/=dt;
            imu->accl[i]/=dt; /* convert to rate/acceleration */
        }
    }
    if (opt->insopt.imuvalfmt==IMUVALFMT_DEG) {
        for (i=0;i<3;i++) imu->gyro[i]*=D2R; /* convert to rad */
    }
}
/* estimate heading using the navigation frame velocity-----------------------
 * args   :  double *vel  I  velocity in navigation frame (ned-frame)
 * return :  heading (rad)
 * --------------------------------------------------------------------------*/
extern double vel2head(const double *vel)
{
    return atan2(vel[1],fabs(vel[0])<1E-4?1E-4:vel[0]);
}
/* check vehicle whether is straight driving --------------------------------*/
extern int chksdri(const double *vel,int n)
{
    int i;
    double *head=mat(n,1),hstd=0.0;

    for (i=0;i<n;i++) {
        head[i]=vel2head(vel+3*i); hstd=stds(head,n);
    }
    free(head); return hstd<3.0*D2R;
}
/* rotate vector convert to quaternion---------------------------------------*/
extern void rvec2quat(const double *rvec,double *q)
{
    double rot_ang,cR,sR,rx,ry,rz;

    /* assume always positive */
    rot_ang=sqrt(SQR(rvec[0])+SQR(rvec[1])+SQR(rvec[2]));
    if (fabs(rot_ang)<1E-8) {
        q[0]=1.0; q[1]=q[2]=q[3]=0.0;
    }
    else {
        cR=cos(rot_ang/2.0);
        sR=sin(rot_ang/2.0);
        rx=rvec[0]/rot_ang;
        ry=rvec[1]/rot_ang;
        rz=rvec[2]/rot_ang;
        q[0]=cR; q[1]=rx*sR; q[2]=ry*sR; q[3]=rz*sR;
    }
}
/* sign of given value-------------------------------------------------------*/
static double sign(double val) {return val==0.0?0.0:val<0.0?-1.0:1.0;}
/* dcm (Cba) to rotation vector (r^a)----------------------------------------*/
extern void dcm2rot(const double *C,double *rv)
{
    double sinPHI,cosPHI,PHI,F,sr_a,mu1,mu2,mu3;
    double u1,u2,u3;

    sinPHI=0.5*sqrt(SQR(C[5]-C[7])+SQR(C[6]-C[2])+SQR(C[1]-C[3]));
    cosPHI=0.5*(C[0]+C[4]+C[8]-1.0);

    PHI=atan2(sinPHI,cosPHI);
    if (cosPHI>=0) {
        if (sinPHI==0) F=1;
        else {
            F=PHI/sinPHI;
        }
        rv[0]=0.5*F*(C[5]-C[7]);
        rv[1]=0.5*F*(C[6]-C[2]);
        rv[2]=0.5*F*(C[1]-C[3]);
    }
    else {
        sr_a=1-cosPHI;
        mu1=sqrt((C[0]-1)/sr_a+1);
        mu2=sqrt((C[4]-1)/sr_a+1);
        mu3=sqrt((C[8]-1)/sr_a+1);
        if (mu1>=mu2&&mu1>=mu3) {
            u1=mu1*sign(C[5]-C[7]);
            u2=1.0/2.0/u1*(C[3]+C[1])/sr_a;
            u3=1.0/2.0/u1*(C[6]+C[2])/sr_a;
        }
        else if (mu2>=mu3) {
            u2=mu2*sign(C[6]-C[2]);
            u3=1/2/u2*(C[7]+C[5])/sr_a;
            u1=1/2/u2*(C[3]+C[1])/sr_a;
        }
        else {
            u3=mu3*sign(C[1]-C[3]);
            u1=1/2/u3*(C[6]+C[2])/sr_a;
            u2=1/2/u3*(C[7]+C[5])/sr_a;
        }
        rv[0]=PHI*u1;
        rv[1]=PHI*u2;
        rv[2]=PHI*u3;
    }
}
/* get imu body velocity in ned-frame----------------------------------------
 * args  :  insstate_t *ins  I  ins states
 *          double *vn       O  velocity in ned-frame
 * return: none
 * --------------------------------------------------------------------------*/
extern void getvn(const insstate_t *ins,double *vn)                                  
{
    double pos[3],C[9];

    trace(3,"getvn:\n");

    ecef2pos(ins->re,pos);
    ned2xyz(pos,C);
    matmul("TN",3,1,3,1.0,C,ins->ve,0.0,vn); 
}
/* update ins states in n-frame----------------------------------------------*/
extern void updinsn(insstate_t *ins)
{
    double Cne[9];

    /* position */
    ecef2pos(ins->re,ins->rn);

    /* attitude/velocity */
    ned2xyz(ins->rn,Cne);
    matmul("TN",3,3,3,1.0,Cne,ins->Cbe,0.0,ins->Cbn);
    matmul("TN",3,1,3,1.0,Cne,ins->ve ,0.0,ins->vn );

    /* acceleration */
    matmul("TN",3,1,3,1.0,Cne,ins->ae,0.0,ins->an);
}
/* update ins states in e-frame----------------------------------------------*/
extern void updinse(insstate_t *ins)
{
    double Cne[9];

    /* attitude and velocity */
    ned2xyz(ins->rn,Cne);
    matmul("NN",3,3,3,1.0,Cne,ins->Cbn,0.0,ins->Cbe);
    matmul("NN",3,1,3,1.0,Cne,ins->vn,0.0,ins->ve);

    /* position */
    pos2ecef(ins->rn,ins->re);

#if 1 /* update acceleration in e-frame */
    getaccl(ins->fb,ins->Cbe,ins->re,ins->ve,ins->ae);
#else
    /* update acceleration in n-frame */
    matmul("NN",3,1,3,1.0,Cne,ins->an,0.0,ins->ae);
#endif
}
/* computes the curvature matrix and the gravity-----------------------------*/
static void geoparam(const double *pos,const double *vn,double *wen_n,
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
/* get dcm matrix of b-frame to n-frame--------------------------------------*/
static void getqbn(const insstate_t *ins, double* qbn)
{
    double pos[3],Cne[9],Cbn[9];
    ecef2pos(ins->re,pos);
    ned2xyz(pos,Cne);
    matmul3("TN",Cne,ins->Cbe,Cbn);
    dcm2quatx(Cbn,qbn);
}
/* rotational and sculling motion correction --------------------------------*/
extern void rotscull_corr(insstate_t *ins,const insopt_t *opt,double dt,
                          double *dv,double *da)
{
    double dap[3],dvp[3],dak[3],dvk[3],dv1[3],dv2[3],dv3[3];
    int i;
    for (i=0;i<3;i++) {
        dap[i]=ins->omgbp[i]*ins->dt;
        dvp[i]=ins->fbp  [i]*ins->dt;
    }
    for (i=0;i<3;i++) {
        dak[i]=ins->omgb[i]*dt; dvk[i]=ins->fb[i]*dt;
    }
    cross3(dak,dvk,dv1);
    cross3(dap,dvk,dv2);
    cross3(dvp,dak,dv3);
    for (i=0;i<3&&dv;i++) {
        dv[i]=0.5*dv1[i]+1.0/12.0*(dv2[i]+dv3[i]);
    }
    if (da) {
        cross3(dap,dak,da);
        for (i=0;i<3;i++) da[i]=1.0/12.0*da[i];
    }
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
    double dt,vmid[3],wen_n[3],wie_n[3],Rn,Re,g;
    double dv[3],dv1[3],dv2[3],qbn[4],w[3],rn[3];
    double da[3],qb[4],qn[4],q[4];
    double das[3]={0},dvs[3]={0};
    int i;

    trace(3,"updateinsn:\n");

    trace(5,"ins(-)=\n"); traceins(5,ins);

    /* update ins state in n-frame */
    updinsn(ins);

    /* save precious epoch ins states */
    savepins(ins,data);

    if ((dt=timediff(data->time,ins->time))>MAXDT||fabs(dt)<1E-6) {

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
    matcpy(vmid,ins->vn,3,1);
    matcpy(rn  ,ins->rn,3,1);

    /* geo parameters */
    geoparam(rn,vmid,wen_n,wie_n,&g,&Re,&Rn);

#if SCULL_CORR
    /* rotational and sculling motion correction */
    rotscull_corr(ins,insopt,dt,dvs,das);
#endif
    /* update ins velocity */
    for (i=0;i<3;i++) dv[i]=ins->fb[i]*dt+dvs[i];

    getqbn(ins,qbn);
    quatrot(qbn,dv,0,dv1);

    for (i=0;i<3;i++) w[i]=2.0*wie_n[i]+wen_n[i];
    cross3(ins->vn,w,dv2);
    dv2[2]=g+dv2[2];

    for (i=0;i<3;i++) dv2[i]=dv2[i]*dt;
    for (i=0;i<3;i++) {
        ins->vn[i]=vmid[i]+dv1[i]+dv2[i]; ins->an[i]=(dv1[i]+dv2[i])/dt;
    }
    /* update ins attitude */
    for (i=0;i<3;i++) vmid[i]=(vmid[i]+ins->vn[i])/2.0;
    for (i=0;i<3;i++) da[i]=ins->omgb[i]*dt+das[i]; rvec2quat(da,qb);
    quatmulx(qbn,qb,q);

    /* geo parameters */
    geoparam(rn,vmid,wen_n,wie_n,&g,&Re,&Rn);
    for (i=0;i<3;i++) {
        da[i]=-(wen_n[i]+wie_n[i])*dt;
    }
    rvec2quat(da,qn);
    quatmulx(qn,q,qbn); quat2dcmx(qbn,ins->Cbn);

    /*  update ins position */
    ins->rn[0]+=1.0/Rn*vmid[0]*dt;
    ins->rn[1]+=1.0/Re/cos(rn[0])*vmid[1]*dt;
    ins->rn[2]-=vmid[2]*dt;

    /* update ins states in e-frame */
    updinse(ins);

    ins->dt=dt;
    ins->ptime=ins->time;
    ins->time =data->time;
    ins->stat =INSS_MECH;

    trace(5,"ins(+)=\n"); traceins(5,ins);
    return 1;
}
/* calculates specific force and angular rate from input w.r.t and resolved----
 * along ecef-frame axes
 * args:    double dt     I  time interval between epochs (s)
 *          double *Cbe1  I  body-to-ecef-frame coordinate transformation matrix
 *          double *Cbe0  I  previous body-to-ecef-frame coordinate transformation
 *                           matrix
 *          double *ve1   I  velocity of body frame w.r.t. ecef frame, resolved along
 *                           ecef-frame axes (m/s)
 *          double *ve0   I  previous velocity of body frame w.r.t. ecef frame,
 *                           resolved along ECEF-frame axes (m/s)
 *          double *re    I  cartesian position of body frame w.r.t. ecef frame,
 *                           resolved along ecef-frame axes (m)
 *          double *fb    O  specific force of body frame
 *          double *omgb  O  angular rate of body frame
 * return : status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int kinematicsecef(const double dt,const double *Cbe1,const double *Cbe0,
                          const double *ve1,const double *ve0,
                          const double *re1,
                          double *fb,double *omgb)
{
    double Cbb[9]={1,0,0,0,1,0,0,0,1};
    double we=OMGE*dt,We[9],Con[9],omg[3],fe[3],wie[3]={0};
    double s,ge[3],W[9],T[9],a,a1,a2,Ab[9];
    double fibe[3],Ab2[9],ae[3]={0},Omg[9],Ca[9];
    int i;

    trace(3,"kinematicsecef:\n");

    if (dt<=0.0) {
        setzero(fb,1,3); setzero(omgb,1,3);
        return 1;
    }
    We[0]= cos(we); We[3]=sin(we); We[6]=0.0;
    We[1]=-sin(we); We[4]=cos(we); We[7]=0.0;
    We[2]=0.0;
    We[5]=0.0;
    We[8]=1.0;

    /* approximate angular rate w.r.t. inertial frame */
    matmul33("TNN",Cbe1,We,Cbe0,3,3,3,3,Con);
    omg[0]=0.5*(Con[7]-Con[5]);
    omg[1]=0.5*(Con[2]-Con[6]);
    omg[2]=0.5*(Con[3]-Con[1]);

    s=acos(0.5*(Con[0]+Con[4]+Con[8]-1.0)); /* calculate and apply the scaling factor */
    if (s>1E-8) {
        for (i=0;i<3;i++) omg[i]=omg[i]*s/sin(s);
    }
    if (omgb) {
        for (i=0;i<3;i++) omgb[i]=omg[i]/dt; /* calculate the angular rate */
    }
    for (i=0;i<3;i++) fe[i]=(ve1[i]-ve0[i])/dt;
    pregrav(re1,ge); wie[2]=OMGE;
    skewsym3(wie,W);

    matmul("NN",3,1,3,1.0,W,ve0,0.0,T);
    for (i=0;i<3;i++) {
        fibe[i]=fe[i]-ge[i]+2.0*T[i];
    }
    skewsym3(omg,Ab); ae[2]=OMGE*dt;

    /* average body-to-ecef coordinate transformation */
    if ((a=norm(omg,3))>1E-8) {
        a1=(1.0-cos(a))/SQR(a); a2=1.0/SQR(a)*(1.0-sin(a)/a);

        matmul3("NN",Ab,Ab,Ab2);
        for (i=0;i<9;i++) Cbb[i]+=a1*Ab[i]+a2*Ab2[i];
        skewsym3(ae,Omg);
        matmul3("NN",Cbe0,Cbb,Ab);
        matmul3("NN",Omg,Cbe0,Ab2);
        for (i=0;i<9;i++) {
            Ca[i]=Ab[i]-0.5*Ab2[i];
        }
    }
    else {
        skewsym3(ae,Omg); matmul3("NN",Omg,Cbe0,Ca);
        for (i=0;i<9;i++) {
            Ca[i]=Cbe0[i]-0.5*Ca[i];
        }
    }
    if (matinv(Ca,3)) {
        return 0;
    }
    if (fb) {
        matmul("NN",3,1,3,1.0,Ca,fibe,0.0,fb);
    }
    return 1;
}
/* generates a uniform distributed random number between min and max --------*/
static double getuniform(double min,double max)
{
    return 1.0*rand()/RAND_MAX*(max-min)+min;
}
/* creates gaussian distributed random numbers (Box-Müller method)-----------*/
static double getgaussian(double std)
{
    if (std<0.0) std=-std;
    double x1,x2,w,y1;
    do {
        x1=getuniform(-1.0,1.0); x2=getuniform(-1.0,1.0); w=x1*x1+x2*x2;
    } while (w>=1.0);
    w=sqrt((-2.0*log(w))/w); y1=x1*w; return std*y1;
}
/* simulates an inertial measurement unit-------------------------------------
 * args:    double *fb      IO  specific force of body frame (m/s^2)
 *          double *omgb    IO  angular rate of body frame (rad/s)
 *          imu_err_t *err  I   imu error model
 *          double dt       I   time interval between epochs (s)
 * return: status (1: ok,0: fail)
 * --------------------------------------------------------------------------*/
extern int simimumeas(double *fb,double *omgb,const imu_err_t *err,double dt)
{
    double ng[3],na[3],Sg[9],Sa[9],I[9]={1,0,0,0,1,0,0,0,1};
    double domg[3],fb0[3],omg0[3];
    int i;

    trace(3,"simimumeas:\n");

    if (fabs(dt)<1E-10) return 0;
    for (i=0;i<3;i++) {
        ng[i]=err->wgn[i]/sqrt(dt)*getgaussian(0.5);
        na[i]=err->wan[i]/sqrt(dt)*getgaussian(0.5);
    }
    for (i=0;i<9;i++) {
        Sg[i]=err->Mg[i]+I[i]; Sa[i]=err->Ma[i]+I[i];
    }
    matmul("NN",3,1,3,1.0,Sg,omgb,0.0,omg0);
    matmul("NN",3,1,3,1.0,Sa,fb,0.0,fb0);
    matmul("NN",3,1,3,1.0,err->Gg,fb,0.0,domg);

    for (i=0;i<3;i++) {
        omgb[i]=omg0[i]+ng[i]+domg[i]+err->bg[i];
        fb[i]=fb0[i]+na[i]+err->ba[i];
    }
    return 1;
}
/* trace ins states ----------------------------------------------------------*/
extern void traceinss(int level, const double *Cbe,const double *re,
                      const double *ve,gtime_t time)
{
    double pos[3],Cne[9],Cnb[9],rpy[3],vel[3];
    char s[64];
    time2str(time,s,3);
    trace(level,"time  =%s\n",s);
    ecef2pos(re,pos);
    ned2xyz(pos,Cne);
    matmul3("TN",Cbe,Cne,Cnb);
    dcm2rpy(Cnb,rpy);
    trace(level,"attn =%8.5f %8.5f %8.5f\n",rpy[0]*R2D,rpy[1]*R2D,rpy[2]*R2D);

    if (ve) {
        matmul3v("T",Cne,ve,vel);
        trace(level,"veln =%8.5f %8.5f %8.5f\n",
              vel[0],vel[1],vel[2]); /* n-frame */

        matmul3v("T",Cbe,ve,vel);
        trace(level,"velb =%8.5f %8.5f %8.5f\n",
              vel[0],vel[1],vel[2]); /* b-frame */
    }
}

