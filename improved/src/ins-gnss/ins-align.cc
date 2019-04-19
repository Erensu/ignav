/*------------------------------------------------------------------------------
* ins-align.cc : ins initial alignment functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *       for IMU calibration without external equipments,2014.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/10/12 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define MIN_ACCL    (0.05)           /* check whether it is static for accl imu data (m/s^2) */
#define MIN_GRYO    (0.05*D2R)       /* check whether it is static for gyro imu data (rad/s) */
#define VAR_VEL     SQR(1.0)         /* initial velocity variance (m/s)^2 */
#define VAR_MVEL    SQR(0.1)         /* velocity measurements variance (m/s)^2 */
#define EPS         1E-2             /* tolerance of check two quaternions whether is close */
#define EPSX        1E-3             /* tolerance of check estimated states */

/* global variable ----------------------------------------------------------*/
static double I[9]={1,0,0,0,1,0,0,0,1};
extern const double Cen[9]={0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0};
extern const double Crf[9]={0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,-1.0};

/* check imu data whether is stationaly--------------------------------------*/
static int chkstatic(const imud_t *data,const insopt_t *opt,const double gn)
{
    return (     data->accl[0]    <MIN_ACCL)&&(data->accl[1]<MIN_ACCL)&&
           (fabs(data->accl[2]+gn)<MIN_ACCL)&&(data->gyro[0]<MIN_GRYO)&&
           (     data->gyro[1]    <MIN_GRYO)&&(data->gyro[2]<MIN_GRYO);
}
/* normalization dcm transform matrix-----------------------------------------*/
extern void normdcm(double *C)
{
    int i; double Cp[9],*T=eye(3),*CC=eye(3);

    matcpy(Cp,C,3,3);
    matmul("TN",3,3,3,1.0,Cp,Cp,-1.0,CC);

    for (i=0;i<9;i++) T[i]=T[i]-0.5*CC[i];
    matmul("NN",3,3,3,1.0,T,Cp,0.0,C);
    free(T); free(CC);
}
/* coarse alignment for ins navigation initial with stationaly imu measurement data
 * args  : insstate_t *ins  IO  ins states
 *         imud_t *data     I   imu measurement data (static case) {m/s^2,rad/s}
 *         int    n         I   number of data (0: no use of data)
 *         insopt_t  *opt   I   ins options
 * return : 1:ok, 0:fail
 * note   : this model is suitable for small attitude error's
 * ---------------------------------------------------------------------------*/
extern int coarse_align(insstate_t *ins,const imud_t *data,int n,
                        const insopt_t *opt)
{
    int i,j,k;
    double pos[3],fb[3]={0},omg[3]={0},gn[3],g,A[9]={0},B[9],Cne[9];
    double fbc[3],omgc[3];

    trace(3,"coarse_asign: n=%d\n",n);

    if (n<=0) return 0;
    ecef2pos(ins->re,pos); ned2xyz(pos,Cne);
    gravity_ned(pos,gn); g=norm(gn,3);

    /* average imu measurement data for remove noise */
    for (k=0,i=0;i<n;i++) {

        /* check static imu data */
        if (opt->align.chkstatic) {
            if (!chkstatic(data+i,opt,g)) continue;
        }
        for (k++,j=0;j<3;j++) {
            fb[j]+=data[i].accl[j]; omg[j]+=data[i].gyro[j];
        }
    }
    if (k<=0) {
        trace(2,"no enough static imu data\n");
        return 0;
    }
    for (i=0;i<3;i++) fb[i]/=k,omg[i]/=k;
    
    /* correction for accl and gyro bias */
    ins_errmodel(fb,omg,fbc,omgc,ins);

    /* get the initial attitude of body-frame */
    A[0]=-tan(pos[0])/g;          A[2]=-1.0/g;
    A[3]= 1.0/(OMGE*cos(pos[0])); A[7]=-1.0/(g*OMGE*cos(pos[0]));

    cross3(fbc,omgc,gn);
    for (i=0;i<3;i++) B[0+i*3]=fbc[i],B[1+i*3]=omgc[i];
    for (i=0;i<3;i++) B[2+i*3]=gn [i];

    /* normalization */
    matmul("NN",3,3,3,1.0,A,B,0.0,ins->Cbn);
    normdcm(ins->Cbn);

    /* normalization */
    matmul("NN",3,3,3,1.0,Cne,ins->Cbn,1.0,ins->Cbe);
    normdcm(ins->Cbe);

    ins->time=data[k].time;
    return 1;
}
/* asign a block matrix to another matrix by giving index -----------------
 * args    :  double *A   IO  input matrix for asign in elements
 *            int m,n     I   rows and cols of input matrix
 *            double *B   I   asign matrix
 *            int p,q     I   rows and cols of asign matrix
 *            int isr,isc I   start row and col to asign matrix
 * return  : none
 * ----------------------------------------------------------------------*/
extern void asi_blk_mat(double *A,int m,int n,const double *B,int p ,int q,
                        int isr,int isc)
{
    int i,j; for (i=isr;i<isr+p;i++) {
        for (j=isc;j<isc+q;j++) A[i+j*m]=B[(i-isr)+(j-isc)*p];
    }
}
/* convert imu measurements to incruments format --------------------------*/
static void imu_con_inc(const imud_t * data,int n,const double dt,double *phim,
                        double *dvbm)
{
    int i,j;
    for (i=0;i<3;i++) phim[i]=dvbm[i]=0.0;
    for (i=0;i<n;i++) for (j=0;j<3;j++) {
        phim[j]+=data[i].gyro[j]*dt; dvbm[j]+=data[i].accl[j]*dt;
    }
}
/* initial align fn method-------------------------------------------------*/
static void afnkfinit(const double dt,const double *pos,const double *phi0,
                      const insopt_t* opt,double *phi,double *Q,double *R,
                      double *P0,double *H,double *x0)
{
    int i;
    double we[3],W[9],I[4]={-1,0,0,-1},*II=eye(5),gn[3];
    const ins_align_t *pa=&opt->align;

    trace(3,"afnkfinit:\n");

    /* initial covariance,measurement variance and system noise matrix */
    for (i=0;i<3;i++) Q [i+i*5]=SQR(pa->web[i])*dt;
    for (i=0;i<2;i++) R [i+i*2]=SQR(pa->wdb[i]/SQRT(dt));
    for (i=0;i<3;i++) P0[i+i*5]=SQR(phi0[i]);
    for (i=3;i<5;i++) P0[i+i*5]=SQR(pa->eb[i-3]);
    for (i=0;i<5;i++) x0[i]=1E-20;

    /* system propagate matrix */
    we[0]=OMGE*cos(pos[0]); we[1]=0.0; we[2]=-OMGE*sin(pos[0]);
    skewsym3(we,W);
    asi_blk_mat(phi,5,5,W,3,3,0,0);
    asi_blk_mat(phi,5,5,I,2,2,0,3);

    gravity_ned(pos,gn); H[5]=gn[2]; H[1]=-gn[2];

    for (i=0;i<5*5;i++) phi[i]=II[i]+phi[i]*dt;
    free(II);
}
/* align fiter--------------------------------------------------------------*/
static int afnfilter(double *x,double *P,const double *Q,const double *R,
                    const double *H,const double *fn,const double *Phi,
                    int nx,int nv)
{
    int i;
    double *x_=mat(nx,1),*P_=mat(nx,nx),*v=mat(nv,1);

    matmul("NN",nx,1,nx,1.0,Phi,x,0.0,x_);
    matmul33("NNT",Phi,P,Phi,nx,nx,nx,nx,P_);
    for (i=0;i<nx*nx;i++) P_[i]+=Q[i];

    matcpy(v,fn,nv,1);
    matmul("TN",nv,1,nx,-1.0,H,x_,1.0,v);

    if (filter(x_,P_,H,v,R,nx,nv)) {
        free(x_); free(P_); free(v);
        return 0;
    }
    matcpy(x,x_,nx,1); matcpy(P,P_,nx,nx);
    free(x_); free(P_); free(v);
    return 1;
}
/* qbn=qpb-phi--------------------------------------------------------------*/
static void qdelphi(quat_t *qpb,const double *phi)
{
    quat_t q,qpbc; quat_copy(&qpbc,qpb); rov2qua(phi,&q); quat_mul(qpb,&q,&qpbc);
}
/* ins initial align uses Kalman filter with fn as measurement --------------
 * args  :  insstate_t *ins  IO  ins states
 *          imud_t *data     I   imu measurement data (static case) {m/s^2,rad/s}
 *          int    n         I   number of data (0: no use of data)
 *          double *phi0     I   initial misalignment angles estimation
 *          insopt_t  *opt   I   ins options
 *          double dt        I   imu measurement data time internal
 *          double *att0     O   attitude asign result
 * return : 1 :ok, 0:fail
 * -------------------------------------------------------------------------*/
static int alignfn(insstate_t *ins,const imud_t *data,int n,const double *phi0,
                   const insopt_t *opt,const double dt,double *rpyo,quat_t *qo)
{
    int i,j,k,ns;
    double phim[3]={0},dvbm[3]={0},Cnn[9],gn[3],wie[3];
    double tv[3],fn[3],t[3],rpy[3],xb[5];
    double *P,*Q,*R,*H,*x,*Phi;
    quat_t q,qinv;

    trace(3,"asignfn: n=%d,dt=%lf\n",n,dt);

    gravity_ned(ins->rn,gn);
    ns=opt->align.ns==0?1:opt->align.ns;

    /* ins alignment initial */
    x=zeros(5,1); P=zeros(5,5); Q=zeros(5,5); R=zeros(2,2); H=zeros(5,2);
    Phi=zeros(5,5);

    afnkfinit(dt,ins->rn,phi0,opt,Phi,Q,R,P,H,x);

    wie[0]= OMGE*cos(ins->rn[0]);
    wie[1]= 0.0;
    wie[2]=-OMGE*sin(ins->rn[0]);

    for (i=0;i<3;i++) tv[i]=-wie[i]*ns*dt/2.0;
    rov2dcm(tv,Cnn);

    /* ins initial attitude */
    dcm2quat(ins->Cbn,&q);

    /* start alignment filter */
    for (k=0,i=0;i<n;i+=ns) {

        /* check imu data whether is static */
        if (opt->align.chkstatic) {
            if (!chkstatic(data+i,opt,gn[2])) continue;
        }
        /* get the accl and gyro incruments */
        imu_con_inc(data+i,ns,dt,phim,dvbm);

        /*  update attitude*/
        qinv=invquat(&q); qmulv(wie,&qinv,tv);
        for (j=0;j<3;j++) {
            phim[j]-=tv[j]*dt*ns;
        }
        quatupd(phim,&q);

        /* get the special force in body-frame */
        for (j=0;j<3;j++) tv[j]=dvbm[j]/(ns*dt);
        qmulv(tv,&q,t);
        matmul("NN",3,1,3,1.0,Cnn,t,0.0,fn);

        /* filter updates */
        if (!afnfilter(x,P,Q,R,H,fn,Phi,5,2)) continue;

        /* backup attitude error estimated state */
        matcpy(xb,x,1,3);

        /* feedback attitude error */
        qdelphi(&q,x); for (j=0;j<3;j++) x[j]=1E-15;
        
        k++; /* index of time */
        quat2rpy(&q,rpy);
        trace(3,"rpy=\n"); tracemat(5,rpy,1,3,15,6);
    }
    /* check estimated states */
    if (norm(xb,3)<EPSX) {

        /* save alignment results */
        quat_copy(qo,&q); matcpy(rpyo,rpy,3,1);
    }
    else{
        trace(3,"asignfn : ins initial align uses kalman filter with fn as"
                " measurement fail\n");
        free(x); free(P); free(Q);
        free(R); free(H); free(Phi);
        return 0;
    }
    free(x); free(P); free(Q);
    free(R); free(H); free(Phi);
    return k;
}
/* initial asign vn method --------------------------------------------------*/
static void avnkfinit(const double dt,const double *pos,const double *phi0,
                      const insopt_t* opt,const double *wvn,double *Phi,
                      double *Q,double *R,double *P0,double *H,double *x0)
{
    int i,j;
    double wie[3],W[9],ddt,*I=eye(12),*II=eye(3),*Ht=zeros(3,12);
    const ins_align_t *pa=&opt->align;

    trace(3,"avnkfinit: dt=%lf\n",dt);

    ddt=(pa->ns==0?1.0:pa->ns)*dt;

    for (i=0;i<3;i++) Q[i+i*12]=SQR(pa->web[i  ])*ddt;
    for (i=3;i<6;i++) Q[i+i*12]=SQR(pa->wdb[i-3])*ddt;

    for (i=0;i<3;i++) R[i+i*3]=SQR(wvn[i]);

    for (i=0;i<3; i++) P0[i+i*12]=SQR(phi0[i]);
    for (i=3;i<6; i++) P0[i+i*12]=VAR_VEL;
    for (i=6;i<9; i++) P0[i+i*12]=SQR(pa->eb[i-6]);
    for (i=9;i<12;i++) {
        P0[i+i*12]=SQR(pa->db[i-9]);
    }
    for (i=0;i<12;i++) x0[i]=1E-15;

    wie[0]=-OMGE*cos(pos[0]);
    wie[1]=0.0;
    wie[2]=OMGE*sin(pos[0]);
    skewsym3(wie,W);
    asi_blk_mat(Phi,12,12,W,3,3,0,0);

    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        Phi[i+j*12]=Phi[i+j*12]*ddt+I[i+j*12];
    }
    asi_blk_mat(Ht,3,12,II,3,3,0,3);
    matt(Ht,12,3,H);

    free(I); free(II); free(Ht);
}
/* ned-frame rotation rate in i-frame resolved in ned-frame------------------*/
static void wnin(const double *rn,const double *vn,double *w)
{
    double M,N;
    radii(rn,&M,&N);
    w[0]=OMGE*cos(rn[0])+vn[1]/(N+rn[2]);
    w[1]=-vn[0]/(M+rn[2]);
    w[2]=-OMGE*sin(rn[0])-vn[1]*tan(rn[0])/(N+rn[2]);
}
/* ins initial align uses Kalman filter with vn as measurement --------------*/
static int alignvn(insstate_t *ins,const imud_t *data,int n,const double *phi0,
                   const double *wvn,const insopt_t *opt,const double dt,
                   quat_t *qo,double *rpyo)
{
    int i,j,k,ns;
    double Cnn[9],tv[3],we[3],phim[3]={0},dvbm[3]={0},rpy[3],xb[12];
    double gn[3],Cbn[9],Ct[9],Ct1[9],Wv[9],dvn[3],vn[3],wn[3];
    double *x,*P,*H,*Q,*R,*Phi;
    quat_t q;

    trace(3,"alignvn: n=%d,dt=%lf\n",n,dt);

    gravity_ned(ins->rn,gn);
    ns=opt->align.ns==0?1:opt->align.ns;

    /* ins alignment initial */
    x=zeros(12,1); P=zeros(12,12); Q=zeros(12,12);
    R=zeros(3,3); H=zeros(12,3); Phi=eye(12);
    avnkfinit(dt,ins->rn,phi0,opt,wvn,Phi,Q,R,P,H,x);

    we[0]= OMGE*cos(ins->rn[0]);
    we[1]= 0.0;
    we[2]=-OMGE*sin(ins->rn[0]);
    for (i=0;i<3;i++) tv[i]=-we[i]*ns*dt/2.0;
    rov2dcm(tv,Cnn);

    wnin(ins->rn,ins->vn,wn);

    /* ins initial attitude,position/velocity in ned-frame */
    dcm2quat(ins->Cbn,&q);
    matcpy(Cbn,ins->Cbn,3,3);
    matcpy(vn,ins->vn,1,3);

    for (k=0,i=0;i<n;i++) {
        /* check imu data whether is static */
        if (opt->align.chkstatic) {
            if (!chkstatic(data+i,opt,gn[2])) continue;
        }
        /* get the accl and gyro incruments */
        imu_con_inc(data+i,ns,dt,phim,dvbm);

        /* updates velocity in ned-frame */
        quat2dcm(&q,Cbn);
        matmul33("NNN",Cnn,Cbn,dvbm,3,3,3,1,dvn);
        for (j=0;j<3;j++) {
            vn[j]+=(dvn[j]+gn[j]*ns*dt);
        }
        /* updates attitude */
        matmul("TN",3,1,3,-dt*ns,Cbn,wn,1.0,phim);
        quatupd(phim,&q);

        /* update system propagate matrix */
        for (j=0;j<9;j++) Ct[j]=-Cbn[j]*ns*dt,Ct1[j]=-Ct[j];
        skewsym3(dvn,Wv);
        asi_blk_mat(Phi,12,12,Wv, 3,3,3,0);
        asi_blk_mat(Phi,12,12,Ct, 3,3,0,6);
        asi_blk_mat(Phi,12,12,Ct1,3,3,3,9);

        /* filter updates */
        if (!afnfilter(x,P,Q,R,H,vn,Phi,12,3)) {
            continue;
        }
        /* backup estimated states */
        matcpy(xb,x,1,12);

        /* feedback attitude error and velocity */
        for (j=0;j<2;j++) x[j]=-x[j]; qdelphi(&q,x);
        for (j=0;j<3;j++) x[j]=1E-15;
        for (j=3;j<6;j++) vn[j-3]-=x[j],x[j]=1E-15;
#if 1
        /* feedback wnin */
        wnin(ins->rn,vn,wn);
#endif
        k++; /* index of time */
        quat2rpy(&q,rpy);
        trace(3,"rpy=\n"); tracemat(5,rpy,1,3,15,6);
    }
    /* check estimated states */
    if (norm(xb,6)<EPSX) {
        /* save alignment results */
        quat_copy(qo,&q); matcpy(rpyo,rpy,3,1);
    }
    else{
        trace(3,"ins initial align use kalman filter with fn as measurement fail\n");
        free(x); free(P); free(Q);
        free(R); free(H); free(Phi);
        return 0;
    }
    free(x); free(P); free(Q);
    free(R); free(H); free(Phi);
    return k;
}
/* returns true if the two input quaternions are close to each other --------*/
static int arequatclose(const quat_t *q1,const quat_t *q2)
{
    quat_t dq,qinv; quat_inv(q1,&qinv); quat_mul(&dq,&qinv,q2);
    if ((fabs(1.0-dq.w)<EPS)&&(fabs(norm(dq.vec+1,3))<EPS)) return 1;
    return 0;
}
/* get an average (mean) from more then two quaternions (with two, slerp would
 * be used),note: this only works if all the quaternions are relatively close
 * together -----------------------------------------------------------------*/
static void avgquat(quat_t *qavg,quat_t *q1,quat_t *q2)
{
    quat_normalize_self(q1);
    quat_normalize_self(q2);
    int i; for (i=0;i<4;i++) {
        qavg->vec[i]=0.5*(q1->vec[i]+q2->vec[i]);
    }
}
/* precise asignment for ins navigation initial with stationaly imu measurement
 * data
 * args  : insstate_t *ins  IO  ins states
 *         imud_t *data     I   imu measurement data (static case) {m/s^2,rad/s}
 *         int    n         I   number of data (0: no use of data)
 *         insopt_t  *opt   I   ins options
 * return : 1:ok, 0:fail
 * note   : if no estimate accl. bias and gyro. bias,then consider them as input
 * ---------------------------------------------------------------------------*/
extern int fine_align(insstate_t *ins,const imud_t *data,int n,const insopt_t *opt)
{
    int k1=1,k2=1;
    double rpyfn[3],rpyvn[3];
    const ins_align_t *pas=&opt->align;
    quat_t qfn,qvn,qavg;

    trace(3,"fine_align: n=%d\n",n);

    if (n<=0||pas->dt<=0.0) {
        trace(3,"ins initial align fail\n");
        return 0;
    }
    /* align uses kalman filter with fn as measurement */
    if (opt->align_fn) {
        k1=alignfn(ins,data,n,pas->phi0,opt,pas->dt,rpyfn,&qfn);
    }
    /* align uses kalman filter with vn as measurement */
    if (opt->align_vn) {
        k2=alignvn(ins,data,n,pas->phi0,pas->wvn,opt,pas->dt,
                   &qvn,rpyvn);
    }
    /* check fn and vn method whether is coincidence */
    if (k1&&k2) {
        if (opt->align_fn&&opt->align_vn) {
            if (arequatclose(&qvn,&qfn)) {

                /* average quaternions, now is a simple way*/
                avgquat(&qavg,&qfn,&qvn);
                quat2dcm(&qavg,ins->Cbn);

                /* ins initial time */
                ins->time=data[MIN(k1,k2)].time;
            }
            else {
                trace(3,"fn-alignment and vn-alignment is not coincidence\n");
                return 0;
            }
        }
        else if (opt->align_vn) quat2dcm(&qvn,ins->Cbn);
        else if (opt->align_fn) quat2dcm(&qfn,ins->Cbn);
        else {
            trace(3,"ins initial align failed,because ins initial options "
                    "don't set\n");
            return 0;
        }
    }
    else {
        trace(3,"ins initial align failed\n");
        return 0;
    }
    return 1;
}
/* convert euler angles to attitude quaternion in enu-rfu-frame --------------
 * args   :   double *att  I  att=[pitch; roll; yaw] in radians
 *            double *q    O  attitude quaternion
 * return : none
 * -------------------------------------------------------------------------*/
extern void a2qua(const double *att,double *q)
{
    int i;
    double att2[3],sp,sr,sy,cp,cr,cy;

    for (i=0;i<3;i++) att2[i]=att[i]/2.0;
    sp=sin(att2[0]); sr=sin(att2[1]); sy=sin(att2[2]);
    cp=cos(att2[0]); cr=cos(att2[1]); cy=cos(att2[2]);
    q[0]=cp*cr*cy-sp*sr*sy;
    q[1]=sp*cr*cy-cp*sr*sy;
    q[2]=cp*sr*cy+sp*cr*sy;
    q[3]=cp*cr*sy+sp*sr*cy;
}
/* convert rotation vector to transformation matrix in enu-rfu-frame-----------
 * args   :   double *rv  I  rotation vector
 *            double *m   O  corresponding DCM,such that
 *                           m=I+sin(|rv|)/|rv|*(rvx)+[1-cos(|rv|)]/|rv|^2*(rvx)^2
 *                           where rvx is the askew matrix or rv
 * return : none
 * ---------------------------------------------------------------------------*/
extern void rv2m(const double *rv,double *m)
{
    double xx,yy,zz,n2,a,b,n,arvx,arvy,arvz,bxx,bxy,bxz,byy,byz,bzz;

    xx=rv[0]*rv[0]; yy=rv[1]*rv[1]; zz=rv[2]*rv[2]; n2=xx+yy+zz;
    if (n2<1E-8) {
        a=1.0-n2*(1.0/6.0-n2/120.0); b=0.5-n2*(1/24.0-n2/720.0);
    }
    else {
        n=sqrt(n2); a=sin(n)/n; b=(1.0-cos(n))/n2;
    }
    arvx=a*rv[0]; arvy=a*rv[1];       arvz=a*rv[2];
    bxx =b*xx;    bxy =b*rv[0]*rv[1]; bxz =b*rv[0]*rv[2];
    byy =b*yy;    byz =b*rv[1]*rv[2]; bzz =b*zz;

    m[0]=1     -byy-bzz; m[3]= -arvz+bxy;     m[6]=  arvy+bxz;
    m[1]=  arvz+bxy;     m[4]=1     -bxx-bzz; m[7]= -arvx+byz;
    m[2]= -arvy+bxz;     m[5]=  arvx+byz;     m[8]=1     -bxx-byy;
}
/* quaternion conjugation ---------------------------------------------------
 * args   :   double *qi  I  input quaternion
 *            double *qo  O  output quaternion
 * return : none
 * -------------------------------------------------------------------------*/
extern void qconj(const double *qi,double *qo)
{
    qo[0]=qi[0]; qo[1]=-qi[1]; qo[2]=-qi[2]; qo[3]=-qi[3];
}
/* 3x1 vector coordinate transformation by quaternion in enu-rfu-frame-----
 * args   :   double *q  I  transformation quaternion
 *            double *vi I  vector to be transformed
 *            double *vo O  output vector, such that vo=q*vi*conjugation(q)
 * return : none
 * -----------------------------------------------------------------------*/
extern void qmulve(const double *q,const double *vi,double *vo)
{
    double qo1,qo2,qo3,qo4;
    qo1=-q[1]*vi[0]-q[2]*vi[1]-q[3]*vi[2];
    qo2= q[0]*vi[0]+q[2]*vi[2]-q[3]*vi[1];
    qo3= q[0]*vi[1]+q[3]*vi[0]-q[1]*vi[2];
    qo4= q[0]*vi[2]+q[1]*vi[1]-q[2]*vi[0];
    vo[0]=-qo1*q[1]+qo2*q[0]-qo3*q[3]+qo4*q[2];
    vo[1]=-qo1*q[2]+qo3*q[0]-qo4*q[1]+qo2*q[3];
    vo[2]=-qo1*q[3]+qo4*q[0]-qo2*q[2]+qo3*q[1];
}
/* quaternion multiplication: q=q1*q2 ned-enu-frame-----------------------*/
extern void qmul(const double *q1,const double *q2,double *q)
{
    q[0]=q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
    q[1]=q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2];
    q[2]=q1[0]*q2[2]+q1[2]*q2[0]+q1[3]*q2[1]-q1[1]*q2[3];
    q[3]=q1[0]*q2[3]+q1[3]*q2[0]+q1[1]*q2[2]-q1[2]*q2[1];
}
/* attitude quaternion updating using rotation vector enu-rfu-frame--------
 * args   :   double *rv  I  roation vector
 *            double *q   IO output quaternion, such that q=qmul(q1,rv2q(rv))
 * return : none
 * ----------------------------------------------------------------------*/
extern void qupdt(const double *rv,double *q)
{
    int i;
    double n2,c,s,n_2,n,q2[4],nq2,q1[4];

    matcpy(q1,q,1,4);
    matmul("NT",1,1,3,1.0,rv,rv,0.0,&n2);
    if (n2<1E-8) {
        c=1.0-n2*(1.0/8.0-n2/384.0); s=1.0/2.0-n2*(1.0/48.0-n2/3840.0);
    }
    else{
        n=sqrt(n2); n_2=n/2.0;
        c=cos(n_2); s=sin(n_2)/n;
    }
    q2[0]=c; for (i=0;i<3;i++) q2[i+1]=s*rv[i];
    q[0]=q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
    q[1]=q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2];
    q[2]=q1[0]*q2[2]+q1[2]*q2[0]+q1[3]*q2[1]-q1[1]*q2[3];
    q[3]=q1[0]*q2[3]+q1[3]*q2[0]+q1[1]*q2[2]-q1[2]*q2[1];
    matmul("NT",1,1,4,1.0,q,q,0.0,&nq2);
    if (nq2>1.000001||nq2<0.99999) {
        for (i=0;i<4;i++) q[i]=q[i]/sqrt(nq2);
    }
}
/* convert attitude quaternion to euler attitude angles in enu-rfu-frame ----
 * args   :   double *q   I  attitude quaternion
 *            double *att O  euler angles att=[pitch; roll; yaw] in radians
 * return : none
 * -------------------------------------------------------------------------*/
extern void q2att(const double *q,double *att)
{
    double q11,q22,q33,q44,q12,q23,q34,q13,q24,q14,C12,C22,C31,C32,C33;
    q11=q[0]*q[0]; q12=q[0]*q[1]; q13=q[0]*q[2]; q14=q[0]*q[3];
    q22=q[1]*q[1]; q23=q[1]*q[2]; q24=q[1]*q[3];
    q33=q[2]*q[2]; q34=q[2]*q[3];
    q44=q[3]*q[3];
    C12=2.0*(q23-q14); C22=q11-q22+q33-q44;
    C31=2.0*(q24-q13); C32=2.0*(q34+q12); C33=q11-q22-q33+q44;
    att[0]=asin(C32);
    att[1]=atan2(-C31,C33);
    att[2]=atan2(-C12,C22);
}
/* onvert rotation vector to transformation quaternion in enu-efu-frame------
 * args   :   double *rv  I  rotation vector
 *            double *q   O  corresponding transformation quaternion, such that
 *                           q=[cos(|rv|/2);sin(|rv|/2)/|rv|*rv]
 * return : none
 * -------------------------------------------------------------------------*/
extern void rv2q(const double *rv,double *q)
{
    double n2,n_2,s,n;
    n2=rv[0]*rv[0]+rv[1]*rv[1]+rv[2]*rv[2];
    if (n2<1E-8) {
        q[0]=1.0-n2*(1.0/8.0-n2/384.0); s=1.0/2.0-n2*(1.0/48.0-n2/3840.0);
    }
    else {
        n=sqrt(n2); n_2=n/2; q[0]=cos(n_2); s=sin(n_2)/n;
    }
    q[1]=s*rv[0]; q[2]=s*rv[1]; q[3]=s*rv[2];
}
/* get the accurate quaternion from calculated quaternion and misalignment--
 * angles. It can be Denoted as 'qnb = qpb - phi', where qpb is calculated
 * quaternion and phi is misalignment angles
 * args   :   double *qpb  I  attitude quaternion from body-frame to computed
 *                            nav-frame
 *            double *phi  I  platform misalignment angles from computed
 *                            nav-frame to ideal nav-frame
 *            double *qnb  O  attitude quaternion from body-frame to ideal
 *                            nav-frame
 * return : none
 * note : this function is in enu-rfu-frame
 * -----------------------------------------------------------------------*/
extern void qdelphi(const double *phi,double *qnb)
{
    double q[4],qt[4]; matcpy(qt,qnb,1,4); rv2q(phi,q); qmul(q,qt,qnb);
}
/* convert attitude quaternion to direction cosine matrix(DCM)------------
 * args   :   double *q  I  attitude quaternion
 *            double *C  O  DCM from body-frame to navigation-frame
 * return : none
 * ----------------------------------------------------------------------*/
extern void q2mat(const double *qnb,double *C)
{
    double q11,q12,q13,q14,q22,q33,q44,q23,q24,q34;
    q11=qnb[0]*qnb[0]; q12=qnb[0]*qnb[1]; q13=qnb[0]*qnb[2]; q14=qnb[0]*qnb[3];
    q22=qnb[1]*qnb[1]; q23=qnb[1]*qnb[2]; q24=qnb[1]*qnb[3];
    q33=qnb[2]*qnb[2]; q34=qnb[2]*qnb[3];
    q44=qnb[3]*qnb[3];
    C[0]=q11+q22-q33-q44; C[1]=2.0*(q23+q14);   C[2]=2.0*(q24-q13);
    C[3]=2.0*(q23-q14);   C[4]=q11-q22+q33-q44; C[5]=2.0*(q34+q12);
    C[6]=2.0*(q24+q13);   C[7]=2.0*(q34-q12);   C[8]=q11-q22-q33+q44;
}
/* ins initial align uses Kalman filter with fn as measurement -----------
 * kalman filter states: [phiE, phiN, phiU, eby, ebz]
 * args   :   double *imu       I  imu measurements data
 *            double *q         I  coarse attitude quaternion (or att)
 *            double *pos       I  ins position (lat,lon,h)
 *            double *phi0      I  initial misalignment angles estimation
 *            ins_align_t *pasi I  ins initial alignments options
 *            double dt         I  imu data sampling interval
 *            int n             I  numbers of imu measurements data
 *            double *att0      O  attitude align result
 *            double *qo        O  attitude align quaternion
 * return : >1:ok,0:fail
 * note : this function is implemented in enu-rfu-frame
 * ---------------------------------------------------------------------*/
extern int alignfnex(const double *imu,const double *q0,const double *pos,
                     const double *phi0,const ins_align_t *pasi,double dt,
                     int n,double *att0,double *qo)
{
    int i,j,k;
    double wnie[3],Cnn[9],phim[3],dvbm[3],tv[3],ndt,T[9],II[4]={-1,0,0,-1};
    double gn[3],q[4],qinv[4],fn[3],tt[3],rpy[3],xb[5];
    double *x,*P,*Q,*R,*H,*F,*I=eye(5);
    imud_t data;

    trace(3,"alignfnex: n=%d,dt=%lf\n",n,dt);

    ndt=(pasi->ns==0?1.0:pasi->ns)*dt;

    gravity_ned(pos,gn); /* gravity in ned-frame */
    wnie[0]=0.0;
    wnie[1]=OMGE*cos(pos[0]);
    wnie[2]=OMGE*sin(pos[0]);
    for (i=0;i<3;i++) tv[i]=-wnie[i]*ndt/2.0; rv2m(tv,Cnn);

    Q=zeros(5,5); P=zeros(5,5); F=zeros(5,5);
    x=zeros(5,1); R=zeros(2,2); H=zeros(2,5);

    for (i=0;i<3;i++) Q[i+i*5]=SQR(pasi->web[i])*ndt;
    for (i=0;i<2;i++) R[i+i*2]=SQR(pasi->wdb[i]/SQRT(ndt));
    for (i=0;i<3;i++) P[i+i*5]=SQR(phi0[i]);
    for (i=3;i<5;i++) P[i+i*5]=SQR(pasi->eb[i-3]);
    for (i=0;i<5;i++) x[i]=1E-15;

    for (i=0;i<3;i++) tv[i]=-wnie[i];
    skewsym3(tv,T); matcpy(q,q0,1,4);
    asi_blk_mat(F,5,5,T, 3,3,0,0);
    asi_blk_mat(F,5,5,II,2,2,1,3);

    for (i=0;i<5*5;i++) F[i]=I[i]+F[i]*ndt;
    H[1]=-gn[2]; H[5]=gn[2];

    /* start alignment filter */
    for (k=0,i=0;i<n;i++) {

        for (j=0;j<3;j++) {
            data.gyro[j]=imu[7*i+j]; data.accl[j]=imu[7*i+3+j];
        }
        /* check imu data whether is static */
        if (!chkstatic(&data,NULL,-gn[2])) continue;

        /* get the accl and gyro incruments */
        imu_con_inc(&data,pasi->ns,dt,phim,dvbm);
        
        /* attitude updates */
        qconj(q,qinv); qmulve(qinv,wnie,tv);
        for (j=0;j<3;j++) tv[j]*=ndt;
        for (j=0;j<3;j++) phim[j]-=tv[j];
        qupdt(phim,q);

        /* updates gravity in enu-frame */
        for (j=0;j<3;j++) tv[j]=dvbm[j]/ndt;
        qmulve(q,tv,tt); matmul("NN",3,1,3,1.0,Cnn,tt,0.0,fn);

        /* align filter */
        if (!afnfilter(x,P,Q,R,H,fn,F,5,2)) {
            continue;
        }
        /* backup attitude error estimated states */
        matcpy(xb,x,1,3);

        /* feedback attitude error */
        qdelphi(x,q); for (j=0;j<3;j++) x[j]=1E-15;

        k++; /* index of time */
        q2att(q,rpy);
        trace(3,"rpy=\n"); tracemat(5,rpy,1,3,15,6);
    }
    /* check estimated states */
    if (norm(xb,3)<EPSX) {
        /* save alignment results */
        matcpy(qo,q,1,4); matcpy(att0,rpy,1,3);
    }
    else{
        trace(3,"ins initial align use kalman filter with fn as measurement fail\n");
        free(x); free(P); free(Q);
        free(R); free(H); free(F); free(I);
        return 0;
    }
    free(x); free(P); free(Q);
    free(R); free(H); free(F); free(I);
    return k;
}
/* ins initial align uses kalman filter with vn as measurement---------------
 * kalman filter states:
 *     [phiE,phiN,phiU, dvE,dvN,dvU, ebx,eby,ebz, dbx,dby,dbz]
 * args   :   double *imu     I  imu measurements data
 *            int n           I  numbers of imu measurements data
 *            double *q0      I  coarse attitude quaternion (or att)
 *            double *pos     I  ins position (lat,lon,h)
 *            double *phi0    I  initial misalignment angles estimation
 *            double *wvn     I  velocity measurement noise (3x1 vector)
 *            ins_align_t pa  I  ins initial alignments options
 *            double *att0    O  attitude align result
 *            double *qo      O  attitude align quaternion
 * return : >1:ok,0:fail
 * note : this function is implemented in enu-rfu-frame
 * ------------------------------------------------------------------------*/
extern int alignvnex(const double *imu,int n,double dt,const double *q0,
                     const double *pos,const double *phi0,const double *wvn,
                     const ins_align_t *pas,double *att0,double *qo)
{
    int k,i,j;
    double Cnn[9],q[4],Cbn[9],gn[3],phim[3],dvbm[3],dvn[3],vn[3]={0};
    double we[3],tv[3],T[9],Ct[9],Ct1[9],Wv[9],rpy[3],xb[12];
    double *x,*P,*Q,*R,*H,*I=eye(3),*Phi=zeros(12,12),*II=eye(12);
    imud_t data;

    trace(3,"alignvnex: n=%d,dt=%lf\n",n,dt);

    gravity_ned(pos,gn); matcpy(q,q0,1,4);
    we[0]=0.0; we[1]=OMGE*cos(pos[0]); we[2]=OMGE*sin(pos[0]);
    for (i=0;i<3;i++) tv[i]=-we[i]*dt/2.0; rv2m(tv,Cnn);

    x=zeros(12,1); P=zeros(12,12); Q=zeros(12,12);
    R=zeros(3,3); H=zeros(3,12);

    for (i=0;i<12;i++) x[i]=1E-15;
    for (i=0;i<3; i++) P[i+i*12]=SQR(phi0[i]);
    for (i=3;i<6; i++) P[i+i*12]=VAR_VEL;
    for (i=6;i<9; i++) P[i+i*12]=SQR(pas->eb[i-6]);
    for (i=9;i<12;i++) P[i+i*12]=SQR(pas->db[i-9]);

    for (i=0;i<3;i++) R[i+i*3 ]=SQR(wvn[i]);
    for (i=0;i<3;i++) Q[i+i*12]=SQR(pas->web[i])*dt;
    for (i=3;i<6;i++) {
        Q[i+i*12]=SQR(pas->wdb[i-3])*dt;
    }
    for (i=0;i<3;i++) tv[i]=-we[i];
    skewsym3(tv,T);
    asi_blk_mat(Phi,12,12,T,3,3,0,0);
    for (i=0;i<12*12;i++) Phi[i]=Phi[i]*dt+II[i];

    H[3]=H[16]=H[29]=1.0;

    /* start alignment filter */
    for (k=0,i=0;i<n;i++) {

        for (j=0;j<3;j++) {
            data.gyro[j]=imu[7*i+j]; data.accl[j]=imu[7*i+3+j];
        }
        /* check imu data whether is static */
        if (!chkstatic(&data,NULL,-gn[2])) continue;

        /* get the accl and gyro incruments */
        imu_con_inc(&data,pas->ns,dt,phim,dvbm);

        /* velocity updates */
        q2mat(q,Cbn);
        matmul33("NNN",Cnn,Cbn,dvbm,3,3,3,1,dvn);
        for (j=0;j<3;j++) vn[j]+=(dvn[j]-gn[j]*dt);

        /* updates attitude */
        matmul("TN",3,1,3,-dt,Cbn,we,1.0,phim);
        qupdt(phim,q);

        /* update system propagate matrix */
        for (j=0;j<9;j++) Ct[j]=-Cbn[j]*dt,Ct1[j]=-Ct[j];
        skewsym3(dvn,Wv);
        asi_blk_mat(Phi,12,12,Wv, 3,3,3,0);
        asi_blk_mat(Phi,12,12,Ct, 3,3,0,6);
        asi_blk_mat(Phi,12,12,Ct1,3,3,3,9);

        /* filter updates */
        if (!afnfilter(x,P,Q,R,H,vn,Phi,12,3)) {
            continue;
        }
        /* backup estimated states */
        matcpy(xb,x,12,1);

        /* feedback attitude error and velocity */
        qdelphi(x,q);
        for (j=0;j<3;j++) x[j]=1E-15;
        for (j=3;j<6;j++) vn[j-3]-=x[j],x[j]=1E-15;

        k++; /* index of time */
        q2att(q,rpy);
        trace(3,"rpy=\n"); tracemat(5,rpy,1,3,15,6);
    }
    /* check estimated states */
    if (norm(xb,6)<EPSX) {

        /* save alignment results */
        matcpy(qo,q,1,4); matcpy(att0,rpy,1,3);
    }
    else{
        trace(2,"ins initial align use kalman filter with fn as measurement fail\n");
        free(x); free(P); free(Q); free(II);
        free(R); free(H); free(Phi); free(I);
        return 0;
    }
    free(x); free(P); free(Q); free(II);
    free(R); free(H); free(Phi); free(I);
    return k;
}
/* precise alignment for ins navigation initial with stationaly imu measurement
 * data
 * args  : insstate_t *ins  IO  ins states
 *         imudata_t *data  I   imu measurement data (static case) {m/s^2,rad/s}
 *         int    n         I   number of data (0: no use of data)
 *         insopt_t  *opt   I   ins options
 * return : 1:ok, 0:fail
 * note   : if no estimate accl bias and gyro bias,then consider them as input
 *          and this function is implemented in ned-frd-frame
 * ---------------------------------------------------------------------------*/
extern int fine_alignex(insstate_t *ins,const imud_t *data,int n,
                        const insopt_t *opt)
{
    int i,k1=1,k2=1;
    double *imu=mat(7,n),qvn[4],qfn[4],avn[4],afn[4],q0[4],qavg[4],Ct[9];
    const ins_align_t *pas=&opt->align;
    quat_t q;

    trace(3,"fine_asignex: n=%d\n",n);

    if (n<=0||pas->dt<=0.0) {
        trace(3,"ins initial align fail\n");
        free(imu);
        return 0;
    }
    /* get the imu measurements data */
    for (i=0;i<n;i++) {
        imu[7*i+6]=time2secs(data[i].time);
        matmul("TN",3,1,3,1.0,Crf,data[i].gyro,0.0,imu+7*i  );
        matmul("TN",3,1,3,1.0,Crf,data[i].accl,0.0,imu+7*i+3);
    }
    /* ins initial attitude */
    matmul33("TNN",Cen,ins->Cbn,Crf,3,3,3,3,Ct);
    dcm2quat(Ct,&q);
    matcpy(q0,q.vec,1,4);

    /* align uses kalman filter with fn as measurement */
    if (opt->align_fn) {
        k1=alignfnex(imu,q0,ins->rn,pas->phi0,pas,pas->dt,n,afn,qfn);
    }
    /* align uses kalman filter with vn as measurement */
    if (opt->align_vn) {
        k2=alignvnex(imu,n,pas->dt,q0,ins->rn,pas->phi0,
                     pas->wvn,pas,avn,qvn);
    }
    /* check fn and vn method whether is coincidence */
    if (k1&&k2) {
        if (opt->align_fn&&opt->align_vn) {

            if (fabs(qfn[0]-qvn[0])<EPS&&fabs(qfn[1]-qvn[1])<EPS&&
                fabs(qfn[2]-qvn[2])<EPS&&
                fabs(qfn[3]-qvn[3])<EPS) {

                /* average quaternions, now is a simple way*/
                for (i=0;i<4;i++) qavg[i]=(qfn[i]+qvn[i])/2.0;
                q2mat(qavg,Ct);
                matmul33("NNT",Cen,Ct,Crf,3,3,3,3,ins->Cbn);

                /* ins initial time */
                ins->time=data[MIN(k1,k2)].time;
            }
            else {
                trace(2,"fn-alignment and vn-alignment is not coincidence\n");
                free(imu); return 0;
            }
        }
        else if (opt->align_vn) {
            q2mat(qvn,Ct);
            matmul33("NNT",Cen,Ct,Crf,3,3,3,3,ins->Cbn);
        }
        else if (opt->align_fn) {
            q2mat(qfn,Ct);
            matmul33("NNT",Cen,Ct,Crf,3,3,3,3,ins->Cbn);
        }
        else {
            trace(2,"ins initial align failed,due to ins initial options don't set\n");
            free(imu); return 0;
        }
    }
    else {
        trace(2,"ins initial align fail\n");
        free(imu); return 0;
    }
    free(imu); return k1&&k2;
}
/* jacobians of large yaw misalignment angle ---------------------------------*/
static void jacobian_lym(double *x,const double *wnie,const double *fn,double ts,
                         double *F)
{
    int i;
    double ax=x[0],ay=x[1],az=x[2],dvx=x[3],dvy=x[4];
    double caz=cos(az),saz=sin(az),jF[5*5]={0};
    double wN=wnie[1],wU=wnie[2],fnx=fn[0],fny=fn[1],fnz=fn[2];
    double *I=eye(5);

    /* 1st-order state updating */
    x[0]=ax+(-saz*wN+ay*wU)*ts;
    x[1]=ay+((1.0-caz)*wN-ax*wU)*ts;
    x[2]=az+ax*caz*wN*ts;
    x[3]=dvx+((1.0-caz)*fnx+saz*fny-(ay*caz+ax*saz)*fnz)*ts;
    x[4]=dvy+(-saz*fnx+(1.0-caz)*fny-(ay*saz-ax*caz)*fnz)*ts;

    /* jacobian matrix */
    jF[5]= wU; jF[10]=-caz*wN;
    jF[1]=-wU; jF[11]= saz*wN;
    jF[2]= caz*wN;  jF[12]=-ax*saz*wN;
    jF[3]=-saz*fnz; jF[ 8]=-caz*fnz; jF[13]= saz*fnx+caz*fny-(-ay*saz+ax*caz)*fnz;
    jF[4]= caz*fnz; jF[ 9]=-saz*fnz; jF[14]=-caz*fnx+saz*fny-( ay*caz+ax*saz)*fnz;

    /* compute the system matrix */
    for (i=0;i<5*5;i++) F[i]=jF[i]*ts+I[i]; free(I);
}
/* extended kalman filter(EKF) simulation with large yaw misalignment angle-----
 * args  : insstate_t *ins  IO  ins states
 *         imud_t *data     I   imu measurement data (static case) {m/s^2,rad/s}
 *         int    n         I   number of data (0: no use of data)
 *         insopt_t  *opt   I   ins options
 * return : 1:ok, 0:fail
 * note : this function is implemented in ned-frd-frame
 * ----------------------------------------------------------------------------*/
extern int fine_align_lym(insstate_t *ins,const imud_t *data,int n,
                          const insopt_t *opt)
{
    int i,j,info=1;
    double Cnn[9],Cbn[9],q[4],vn[3]={0},we[3],gn[3],tv[3],dt,Ct[9];
    double phim[3],dvbm[3],*imu=mat(7,n),dvn[3],fn[3],rpy[3],xb[5];
    double *x,*P,*Q,*R,*H,*F,*v;
    const ins_align_t *pas=&opt->align;
    imud_t sd;
    quat_t q0;

    trace(3,"fine_align_lym: n=%d\n",n);

    if (n<=0||(dt=pas->dt)<=0.0) {
        trace(3,"ins initial align fail\n");
        free(imu);
        return 0;
    }
    gravity_ned(ins->rn,gn);
    we[0]=0.0; we[1]=OMGE*cos(ins->rn[0]); we[2]=OMGE*sin(ins->rn[0]);
    for (i=0;i<3;i++) tv[i]=-we[i]*dt/2.0; rv2m(tv,Cnn);

    /* convert the imu measurement data to enu-rfu-frame */
    for (i=0;i<n;i++) {
        imu[7*i+6]=time2secs(data[i].time);
        matmul("TN",3,1,3,1.0,Crf,data[i].gyro,0.0,imu+7*i  );
        matmul("TN",3,1,3,1.0,Crf,data[i].accl,0.0,imu+7*i+3);
    }
    /* ins initial attitude */
    matmul33("TNN",Cen,ins->Cbn,Crf,3,3,3,3,Ct);
    dcm2quat(Ct,&q0);
    matcpy(q,q0.vec,1,4);

    x=zeros(5,1); P=zeros(5,5); R=zeros(2,2);
    H=zeros(2,5); Q=zeros(5,5); F=zeros(5,5); v=zeros(2,1);

    for (i=0;i<5;i++) x[i]=1E-15;
    for (i=0;i<3;i++) P[i+i*5]=SQR(pas->phi0[i]);
    for (i=3;i<5;i++) P[i+i*5]=VAR_VEL;

    for (i=0;i<2;i++) R[i+i*2]=VAR_MVEL;
    for (i=0;i<3;i++) Q[i+i*5]=SQR(pas->web[i  ])*dt;
    for (i=3;i<5;i++) {
        Q[i+i*5]=SQR(pas->wdb[i-3])*dt;
    }
    H[3]=H[9]=1.0;

    /* start alignment filter */
    for (i=0;i<n;i++) {

        for (j=0;j<3;j++) {
            sd.gyro[j]=imu[7*i+j]; sd.accl[j]=imu[7*i+3+j];
        }
        /* check imu data whether is static */
        if (!chkstatic(&sd,NULL,-gn[2])) {
            continue;
        }
        /* get the accl and gyro incruments */
        imu_con_inc(&sd,pas->ns,dt,phim,dvbm);

        /* velocity updating */
        q2mat(q,Cbn); matmul33("NNN",Cnn,Cbn,dvbm,3,3,3,3,dvn);
        for (j=0;j<3;j++) vn[j]+=(dvn[j]-gn[j]*dt);

        /* attitude updating */
        matmul("TN",3,1,3,-dt,Cbn,we,1.0,phim);
        qupdt(phim,q);

        /* kalman filter */
        for (j=0;j<3;j++) fn[j]=dvn[j]/dt; jacobian_lym(x,we,fn,dt,F);

        matcpy(v,vn,1,2);
        matmul("TN",2,1,5,-1.0,H,x,1.0,v);
        matmul33("NNT",F,P,F,5,5,5,5,P);
        for (j=0;j<5*5;j++) P[j]+=Q[j];

        if (filter(x,P,H,v,R,5,2)) continue;

        /* backup attitude error estimated states */
        matcpy(xb,x,1,3);

        /* feedback attitude */
        qdelphi(x,q); for (j=0;j<3;j++) x[j]=1E-15;

        q2att(q,rpy);
        trace(3,"rpy=\n"); tracemat(5,rpy,1,3,15,6);

        ins->time=data[i].time;
    }
    /* save ins alignment result */
    if (norm(xb,3)<EPSX) {
        /* convert to ned-frd-frame */
        q2mat(q,Ct);
        matmul33("NNT",Cen,Ct,Crf,3,3,3,3,ins->Cbn);
    }
    else {
        trace(3,"too large yaw misalignment angle\n");
        info=0;
    }
    free(x); free(P); free(H);
    free(R); free(Q); free(F);
    free(imu);
    return info;
}