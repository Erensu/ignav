/*------------------------------------------------------------------------------
 * ins-pose.cc : pose measurement aid for ins navigation functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *        for IMU calibration without external equipments,2014.
 *    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *        INS 2008.
 *    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [5] Timothy D. Barfoot. State Estimation for Robotics.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2017/11/13 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants ----------------------------------------------------------------*/
#define ORDERS        10           /* orders of approximate exponential of matrix */
#define MAXRES_POSE        10.0*D2R     /* max residual for pose filter */
#define VARPOSE       30.0*D2R     /* variance of pose measurement */
#define MAXTIMEDIFF   1.0          /* max difference between ins and pose measurement time */
#define JACOB_LEFT    1            /* use left jacobians of misalignment of camera to imu */

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

/* factorial operation-------------------------------------------------------*/
static double factorial(int n)
{
    if (n==0) return 1; return n*factorial(n-1);
}
/* exponential of a matrix----------------------------------------------------
 * args  :  double *A  I  a input matrix (nxn)
 *          int n      I  rows and cols of matrix A
 *          double *E  O  exponential of a matrix
 * return: none
 * --------------------------------------------------------------------------*/
extern void expmat(const double *A,int n,double *E)
{
    double s,*B,*C;
    int i,j,k;

    trace(3,"expmat:\n");

    C=mat(n,n); B=mat(n,n); seteye(E,n);
    seteye(B,n);

    for (i=0;i<ORDERS;i++) {

        s=1.0/factorial(i+1);
        matmul("NN",n,n,n,s,B,A,0.0,C);

        for (j=0;j<n;j++) {
            for (k=0;k<n;k++) E[j+k*n]+=C[j+k*n];
        }
        matcpy(B,C,n,n);
    }
    free(B); free(C);
}
/* roll-pitch-yaw convention from 1-frame to 2-frame-------------------------
 * args:  double *rpy  I  roll,pitch,yaw {rad}
 *        double *C    O  direction cosine matrix
 * note: C=Rz*Ry*Rx
 * return: none
 * --------------------------------------------------------------------------*/
extern void rpy2c(const double *rpy,double *C)
{
    double c1,c2,c3,s1,s2,s3;

    c1=cos(rpy[0]); c2=cos(rpy[1]); c3=cos(rpy[2]);
    s1=sin(rpy[0]); s2=sin(rpy[1]); s3=sin(rpy[2]);

    C[0]= c2*c3; C[3]= c1*s3+s1*s2*c3; C[6]=s1*s3-c1*s2*c3;
    C[1]=-c2*s3; C[4]= c1*c3-s1*s2*s3; C[7]=s1*c3+c1*s2*s3;
    C[2]=    s2; C[5]=         -s1*c2; C[8]=         c1*c2;
}
/* convert dcm to roll-pitch-yaw---------------------------------------------
 * args  :  double *C   I  direction cosine matrix
 *          double *rpy O  roll,pitch,yaw {rad}
 * note: C=Rz*Ry*Rx
 * return: none
 * --------------------------------------------------------------------------*/
extern void c2rpy(const double *C,double *rpy)
{
    rpy[0]=atan2(-C[5],C[8]);
    rpy[1]=asin ( C[2]);
    rpy[2]=atan2(-C[1],C[0]);
}
/* correct attitude use euler angle-------------------------------------------
 * args  :  double *dphi  I  estimate attitude error (rotation vector)
 *          double *C     I  attitude direction cosine matrix
 *          double *Cc    O  corrected attitude direction cosine matrix
 * return: none
 * --------------------------------------------------------------------------*/
extern void correctatt(const double *dphi,const double *C,double *Cc)
{
    double da[3]={0},ang[3]={0},C1[9],C2[9],C3[9],S[9];
    static const double I1[3]={1,0,0},
                        I2[3]={0,1,0},
                        I3[3]={0,0,1};

    trace(3,"correctatt:\n");

    c2rpy(C,ang);

    Ry(ang[1],C2); Rz(ang[2],C3);

    matmul33("NNN",C3,C2,I1,3,3,3,1,&S[0]);
    matmul("NN",3,1,3,1.0,C3,I2,0.0,&S[3]);
    matcpy(&S[6],I3,1,3);

    if (!matinv(S,3)) {
        matmul("NN",3,1,3,1.0,S,dphi,0.0,da);
    }
#if 0
    /* just for debugs */
    Rx(da[0],C1);
    Ry(da[1],C2);
    Rz(da[2],C3);
    matmul33("NNN",C3,C2,C1,3,3,3,3,S);
    matmul("NN",3,3,3,1.0,S,C,0.0,Cc);
#else
    ang[0]=ang[0]+da[0];
    ang[1]=ang[1]+da[1];
    ang[2]=ang[2]+da[2];
    rpy2c(ang,Cc);
#endif
}
/* addition operation in lie algebras----------------------------------------
 * args:  double *a1  I  first element in lie algebras
 *        double *a2  I  second element in lie algebras
 *        double *a   O  a=a1+a2
 * note:  exp(a)=exp(a1)*exp(a2)
 * return: none
 * --------------------------------------------------------------------------*/
extern void addlie(const double *a1,const double *a2,double *a)
{
    double s1[3],s2[3],s3[3],T1[9],T2[9];
    int i;

    trace(3,"addlie:\n");

    skewsym3(a1,T1); skewsym3(a2,T2);
    matmul("NN",3,1,3,1.0,T1,a2,0.0,s1);
    matmul33("NNN",T1,T1,a2,3,3,3,1,s2);
    matmul33("NNN",T2,T2,a1,3,3,3,1,s3);

    /* approximately value of lie algebras */
    for (i=0;i<3;i++) {
        a[i]=a1[i]+a2[i]+1.0/2.0*s1[i]+1.0/12.0*s2[i]+1.0/12.0*s3[i];
    }
}
/* inverse right/left jacobians of SO3---------------------------------------
 * args:  double *phi  I  lie algebras of SO3
 *        double *Jr   O  right jacobians of lie algebras of SO3
 *        double *Jl   O  left  jacobians of lie algebras of SO3
 * return: none
 * --------------------------------------------------------------------------*/
extern void so3_jac(const double *phi,double *Jr,double *Jl)
{
    double I[9]={1,0,0,0,1,0,0,0,1},T1[9],T2[9];
    double a[3],b,b2;
    int i;

    trace(3,"so3jac:\n");

    b=norm(phi,3);
    if (b<=0.0) {
        if (Jr) {seteye(Jr,3); return;}
        if (Jl) {seteye(Jl,3); return;}
    }
    for (i=0;i<3;i++) a[i]=phi[i]/b;

    matmul("NT",3,3,1,1.0,a,a,0.0,T1);
    skewsym3(a,T2);

    for (b2=b/2.0,i=0;i<9;i++) {
        if (Jr) {
            Jr[i]=b2/tan(b2)*I[i]+(1.0-b2*tan(b2))*T1[i]+b2*T2[i];
        }
        if (Jl) {
            Jl[i]=b2/tan(b2)*I[i]+(1.0-b2*tan(b2))*T1[i]-b2*T2[i];
        }
    }
}
/* initial covariance of camera pose expressed in tangent space from ins states
 * args:  insopt_t *opt    I  ins options
 *        insstate_t *ins  I  ins states
 *        vostate_t *vo    O  visual odometry states
 * return: none
 * --------------------------------------------------------------------------*/
extern void initcamposevar(const insopt_t *opt,insstate_t *ins,vostate_t *vo)
{
    double Jl[9],Pins[9],Cce[9],phi[3];
    int i,j;

    trace(3,"initcamposevar:\n");

    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            Pins[i+3*j]=ins->P[(xiA(opt)+i)+3*(xiA(opt)+j)];
        }
    }
    matmul("NN",3,3,3,1.0,ins->Cbe,vo->Cbc,0.0,Cce);
    so3_log(Cce,phi,NULL);

    so3_jac(phi,NULL,Jl);
    matmul33("NNT",Jl,Pins,Jl,3,3,3,3,vo->P);
}
/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jac_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
}
/* update covariance of camera pose expressed in tangent space based on------
 * measurement variance of vo update from precious to current frame
 * args:  insopt_t *opt  I  ins options
 *        double *var    I  euler angle measurement variance from precious to
 *                          current
 *        vostate_t *vo  O  update visual odometry states
 * note: euler angle from precious to current consider as small
 * return: none
 * --------------------------------------------------------------------------*/
extern void updcamposevar(const insopt_t *opt,const double *var,vostate_t *vo)
{
    double S[9],T[9],Jr[9],phi[3],cov[9]={0};
    int i;

    trace(3,"updcamposevar:\n");

    so3_log(vo->Cce,phi,NULL);
    so3_jac(phi,Jr,NULL);

    jac_prot_pang(vo->Cce,S);
    matmul("NN",3,3,3,1.0,Jr,S,0.0,T);

    cov[0]=var[0];
    cov[4]=var[1];
    cov[8]=var[2];

    matmul33("NNT",T,cov,T,3,3,3,3,S);
    for (i=0;i<9;i++) {
        vo->P[i]=vo->P[i]+S[i];
    }
}
/* update camera pose use measurement data from camera visual odometry--------*/
static void updatecampose(const insopt_t *opt,const pose_meas_t *data,
                          insstate_t *ins)
{
    double C[9]={0};

    trace(3,"updatecampose:\n");

    matcpy(C,ins->vo.Cce,3,3);
    matmul("NT",3,3,3,1.0,C,data->C,0.0,ins->vo.Cce);

    updcamposevar(opt,data->var,&ins->vo);
}
/* update dual antennas pose use measurement data from dual ant.-------------*/
static void updatedualant(const insopt_t *opt,const pose_meas_t *data,
                          insstate_t *ins)
{
    trace(3,"updatedualant:\n");
}
/* jacobian of phi residual wrt. attitude error------------------------------*/
static void jac_dphi_da(const double *phi,double *dphida)
{
    so3_jac(phi,NULL,dphida);
}
/* jacobian of phi residual wrt. misalignment of camera to imu body----------*/
static void jac_dphi_dmar(const double *phi,double *dphidma)
{
    so3_jac(phi,dphidma,NULL);
}
/* jacobian of phi residual wrt. misalignment of camera to imu (left mul.)---*/
static void jac_dphi_dmal(const double *phi,const double *Cbe,double *dphidma)
{
    double J[9]; so3_jac(phi,NULL,J);
    matmul("NN",3,3,3,1.0,Cbe,J,0.0,dphidma);
}
/* interpolate ins pose------------------------------------------------------*/
static void interp_ins_pose(const insstate_t *ins,double dt,double *Cbe)
{
    double tt,rpy[3]={0},t;
    quat_t q1,q2,q;

    trace(3,"interp_ins_pose:\n");

    tt=timediff(ins->time,ins->ptime); t=(dt+tt)/tt;
    if (fabs(tt)<1E-6) {
        matcpy(Cbe,ins->Cbe,3,3); return;
    }
    dcm2rpy(ins->pCbe,rpy); rpy2quat(rpy,&q1);
    dcm2rpy(ins->Cbe ,rpy); rpy2quat(rpy,&q2);

    quat_slerp(&q,&q1,&q2,t);
    quat_to_rh_rot_matrix(&q,Cbe);
}
/* jacobians of phi residual wrt. ins attitude error-------------------------*/
static void jac_dphi_datt(const double *phi,const double *Cbe,const double *Cne,
                          double *dphidatt)
{
    double Jl[9]; so3_jac(phi,NULL,Jl);
    matmul("NT",3,3,3,-1.0,Jl,Cne,0.0,dphidatt);
}
/* jacobians of phi residual wrt. misalignment from b-frame to v-frame-------*/
static void jac_dphi_dvma(const double *phi,const double *Cbe,const double *Cne,
                          double *dphidvma)
{
    double Jl[9],T[9]; so3_jac(phi,NULL,Jl);
    matmul("TN",3,3,3, 1.0,Cne,Cbe,0.0,T);
    matmul("NN",3,3,3,-1.0,Jl,T,0.0,dphidvma);
}
/* convert variance of euler to tangent space--------------------------------*/
static void convar_tgt(const double *Cbe,const double *var,double *covtgt)
{
    double S[9],covelr[9]={0};

    jac_prot_pang(Cbe,S);

    covelr[0]=var[0]==0.0?VARPOSE:var[0];
    covelr[4]=var[1]==0.0?VARPOSE:var[1];
    covelr[8]=var[2]==0.0?VARPOSE:var[2];
    matmul33("NNT",S,covelr,S,3,3,3,3,covtgt);
}
/* check estimated states----------------------------------------------------*/
static int chkest_state(const double *x,const insopt_t *opt)
{
    static int iba=xiBa(opt),nba=xnBa(opt);
    static int ibg=xiBg(opt),nbg=xnBg(opt);
    static int ivm=xiVm(opt),nvm=xnVm(opt);
    int flag=0;

    /* check estimated states */
    if (     x[  0]==DISFLAG&&norm(x+  0,3)>15.0*D2R) flag|=1;
    if (nba&&x[iba]==DISFLAG&&norm(x+iba,3)>1E5*Mg2M) flag|=1;
    if (nbg&&x[ibg]==DISFLAG&&norm(x+ibg,3)>15.0*D2R) flag|=1;
    if (nvm&&x[ivm]==DISFLAG&&norm(x+ivm,3)>15.0*D2R) flag|=1;
    if (flag) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    return 1;
}
/* pose fusion filter based on dual ant. pose measurement--------------------*/
static int poseantfilt(const insopt_t *opt,const pose_meas_t *data,
                       insstate_t *ins,double dt)
{
    double Cbe[9],Cne[9],Cvn[9],Cvn0[9],Ry[9],Rz[9],pos[3];
    double phi[3],phi0[3],omg,omg0,cov[9];
    double dphidvma[9],dphidatt[9];
    double *v,*H,*R,*x,*P;
    int i,j,nx,im,nm,ia,na,nv,ind[3],info=0;

    trace(3,"poseantfilt:\n");

    nx=ins->nx;
    im=xiVm(opt); nm=xnVm(opt); ia=xiA(opt); na=xnA(opt);

    interp_ins_pose(ins,dt,Cbe);

    Ry(-data->rpy[1],Ry);
    Rz(-data->rpy[2],Rz);
    matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cvn);

    ecef2pos(ins->re,pos);
    ned2xyz(pos,Cne);

    matmul33("TNN",Cne,Cbe,ins->Cvb,3,3,3,3,Cvn0);
    so3_log(Cvn0,phi0,&omg0);
    so3_log(Cvn ,phi ,&omg );
    if (fabs(omg-omg0)>MAXRES_POSE) {
        trace(2,"large residual\n");
        return 0;
    }
    v=zeros(1,3); H=zeros(3,nx);
    R=zeros(3,3);
    x=zeros(1 ,nx);
    P=zeros(nx,nx);

    jac_dphi_datt(phi0,Cbe,Cne,dphidatt);
    jac_dphi_dvma(phi0,Cbe,Cne,dphidvma);
    convar_tgt(Cvn,data->var,cov);

    /* H,v and R matrix from measurements */
    for (i=0,nv=0;i<3;i++) {
        if (fabs(v[nv]=phi[i]-phi0[i])>MAXRES_POSE) {
            continue;
        }
        if (H) {
            for (j=ia;j<ia+na;j++) H[j+nv*nx]=dphidatt[i+(j-ia)*3];
            for (j=im;j<im+nm;j++) H[j+nv*nx]=dphidvma[i+(j-im)*3];
        }
        ind[nv]=i;
        nv++;
    }
    if (v&&nv) {
        trace(3,"v=\n"); tracemat(3,v,3,1,12,6);
    }
    if (H&&nv) {
        trace(3,"H=\n"); tracemat(3,H,nx,nv,12,6);
    }
    if (R&&nv) {
        for (i=0;i<nv;i++) {
            for (j=0;j<nv;j++) {
                R[i+j*nv]=cov[ind[i]+ind[j]*3];
            }
        }
        trace(3,"R=\n");
        tracemat(3,R,nv,nv,12,6);
    }
    if (nv<=0) {
        trace(2,"pose fusion filter fail\n");
        free(x); free(P); free(v);
        free(H); free(R);
        return 0;
    }
    matcpy(P,ins->P,nx,nx);

    /* ekf filter for pose fusion */
    if (!filter(x,P,H,v,R,nx,nv)) {

        if (!chkest_state(x,opt)) {
            goto exit;
        }
        if (nm&&x[im]!=DISFLAG) corratt(x+im,ins->Cvb);

        dcm2rpy(ins->Cvb,phi);
        phi[0]*=R2D;
        phi[1]*=R2D;
        phi[2]*=R2D;

        trace(3,"v-frame and b-frame misalignment rpy=\n");
        tracemat(3,phi,1,3,12,6);

        matcpy(ins->P,P,nx,nx);
        clp(ins,opt,x);

        /* pose fusion ok */
        ins->pose=1;
        info=1;
    }
    else {
        trace(2,"filter fail\n");
    }
exit:
    free(x); free(P); free(v);
    free(H); free(R);
    return info;
}
/* pose fusion filter based on camera pose measurement-----------------------*/
static int posecamfilt(const insopt_t *opt,const pose_meas_t *data,
                       insstate_t *ins,double dt)
{
#if 0
    double *v,*x,*H,*R,*P,phi0[3],phi[3],Cce[9],r[3],S[9],T[9];
    double dphida[9],dphidma[9],Cbe[9],phi1[3];
    int i,j,nv=0,nx=ins->nx,ia,na,im=0,nm=0;
    int info=0;

    trace(3,"posecamfilt:\n");

    interp_ins_pose(ins,dt,Cbe);

    matmul("NN",3,3,3,1.0,ins->Cbe,ins->vo.Ccb,0.0,Cce);
#if JACOB_LEFT
    so3_log(ins->vo.Ccb,phi0,NULL);
    matmul("NN",3,1,3,0.0,Cbe,phi0,0.0,phi1);
#else
    so3_log(Cce,phi0,NULL);
#endif
    so3_log(ins->vo.Cce,phi,NULL);

    v=zeros(1,3); x=zeros(nx,1);
    H=zeros(3,nx);
    R=zeros(3,3);
    P=zeros(nx,nx);

    ia=xiA (opt); na=xnA (opt);
    im=xiCm(opt); nm=xnCm(opt);

    jac_dphi_da  (phi0,dphida );
    jac_dphi_dmar(phi0,dphidma);

#if JACOB_LEFT
    jac_dphi_dmal(phi0,Cbe,dphidma);
#endif

#if UPD_IN_EULER
    jac_prot_pang(ins->Cbe,S);
    matcpy(T,dphida,3,3);
    matmul("NN",3,3,3,1.0,T,S,0.0,dphida);
#endif
    /* build H,v and R matrix from input measurements */
    for (i=0;i<3;i++) {

#if JACOB_LEFT
        if (fabs(v[nv]=phi[i]-phi1[i])>MAXRES_POSE) {
            continue;
        }
#else
        if (fabs(v[nv]=phi[i]-phi0[i])>MAXRES_POSE) {
            continue;
        }
#endif
        if (H) {
            for (j=ia;j<ia+na;j++) H[j+nv*nx]=-dphida [i+(j-ia)*3];
            for (j=im;j<im+nm;j++) H[j+nv*nx]=-dphidma[i+(j-im)*3];
        }
        r[nv]=data->var[i]==0?SQR(VARPOSE):data->var[i];
        nv++;
    }
    if (v&&nv) {
        trace(3,"v=\n"); tracemat(3,v,3,1,12,6);
    }
    if (H&&nv) {
        trace(3,"H=\n"); tracemat(3,H,nx,nv,12,6);
    }
    if (R&&nv) {
        for (i=0;i<nv;i++) R[i+i*nv]=r[i];
        trace(3,"R=\n");
        tracemat(3,R,nv,nv,12,6);
    }
    if (nv<=0) {
        trace(2,"pose fusion filter fail\n");
        free(x); free(P); free(v);
        free(H); free(R);
        return 0;
    }
    matcpy(P,ins->P,nx,nx);
    if (!filter(x,P,H,v,R,nv,nv)) {

        /* update estimated states */
        if (nm&&x[im]!=DISFLAG) {

#if JACOB_LEFT
            corratt(x+im,ins->vo.Ccb);
#else
            so3_exp(x+im,T);
            matcpy(S,ins->vo.Ccb,3,3);
            matmul("NN",3,3,3,1.0,S,T,0.0,ins->vo.Ccb);
#endif
        }
        matcpy(ins->P,P,nx,nx);
        clp(ins,opt,x);

        /* pose fusion ok */
        ins->pose=1;
        info=1;
    }
    else {
        trace(2,"filter fail\n");
    }
    free(x); free(P); free(v);
    free(H); free(R);
    return info;
#endif
}
/* external pose measurement to update ins navigation------------------------*/
extern int posefusion(const insopt_t *opt,const pose_meas_t *data,
                      insstate_t *ins,int flag)
{
    trace(3,"posefusion: time=%s\n",time_str(ins->time,3));
    double dt;

    /* reset pose fusion status flag */
    ins->pose=0;

    if (!data->time.time||norm(data->C,9)<=0.0) {
        trace(2,"pose fusion fail\n");
        return 0;
    }
    if (flag==INSUPD_TIME) {

        /* update camera or others pose use measurement */
        switch (data->type) {
            case POSE_CAM     : updatecampose(opt,data,ins); break;
            case POSE_DUAL_ANT: updatedualant(opt,data,ins); break;
        }
    }
    else {
        /* pose measurement fusion */
        if (fabs(dt=timediff(data->time,ins->time))>MAXTIMEDIFF) {
            trace(2,"time isn't synchronize\n");
            return 0;
        }
        switch (data->type) {
            case POSE_CAM     : {
                return posecamfilt(opt,data,ins,dt);
            }
            case POSE_DUAL_ANT: {
                return poseantfilt(opt,data,ins,dt);
            }
        }
    }
    return 1;
}
