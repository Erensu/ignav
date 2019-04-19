/*-----------------------------------------------------------------------------
* ins-dualant-calib.cc : dual antennas calibration functions
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
* history : 2018/05/29 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

#define MAXTIMEDIFF  1.0          /* max difference between ins and pose measurement time */
#define VARPOSE      30.0*D2R     /* variance of pose measurement */
#define MINVEL       5.0          /* min velocity for calibration dual antennas */
#define MINNUMMEAS   100          /* min number of pose measurement for initial */
#define ADD_NOISE    1            /* add system noise to covariance */

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

/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jac_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
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
/* covariance of attitude----------------------------------------------------*/
static void covatt(const sol_t *sol,double *cov)
{
    cov[0]=sol->qa[0];
    cov[4]=sol->qa[1];
    cov[8]=sol->qa[2];

    cov[1]=cov[3]=sol->qa[3];
    cov[2]=cov[6]=sol->qa[4];
    cov[5]=cov[7]=sol->qa[5];
}
/* pose measurement covariance-----------------------------------------------*/
static void meascov(const pose_meas_t *pose,const sol_t *sol,double *cov,
                    double *Cvn)
{
    double Cne[9],Cnb[9],Cbn[9],Jl[9],phi[3],T[9],att[3];
    double Ry[9],Rz[9],I[18],W[36]={0};
    double cov_be[9];
    double cov_bn[9],cov_vn[9];
    int i,j;

    trace(3,"meascov: time=%s\n",time_str(pose->time,4));

    ned2xyz(sol->rr,Cne);

    for (i=0;i<3;i++) {
        att[i]=sol->att[i]*D2R;
    }
    rpy2dcm(att,Cnb); matt(Cnb,3,3,Cbn);
    so3_log(Cbn,phi,NULL);
    so3_jac(phi,NULL,Jl);

    covatt(sol,cov_be);
    matmul("NT",3,3,3,1.0,Jl,Cne,0.0,T);
    matmul33("NNT",T,cov_be,T,3,3,3,3,cov_bn);

    Ry(-pose->rpy[1],Ry);
    Rz(-pose->rpy[2],Rz);
    matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cvn);

    so3_log(Cvn,phi,NULL);
    so3_jac(phi,NULL,Jl);
    matmul("NN",3,3,3,1.0,Jl,Cnb,0.0,T);

    convar_tgt(Cvn,pose->var,cov_vn);

    for (i=0;i<3;i++) for (j=0;j<3;j++) W[i+i*6]=cov_vn[i+j*3];
    for (i=3;i<6;i++) {
        for (j=3;j<6;j++) W[i+j*6]=cov_bn[(i-3)+(j-3)*3];
    }
    for (i=0;i<9;i++) {
        I[i]=T[i]; I[9+i]=-T[i];
    }
    matmul33("NNT",I,W,I,3,6,6,3,cov);
}
/* ekf filter matrix---------------------------------------------------------*/
static int build_ekf(const insopt_t *opt, const sol_t* sol,
                     const pose_meas_t *pose, const double *Cvb0,
                     double *v,double *H,double *R)
{
    double phi[3],phi0[3],Cnb[9],Cvn[9],Cvb[9],att[3];
    double Jl[9],R0[9];
    int i,j,ind[3],nv=0;

    trace(3,"build_ekf:\n");

    meascov(pose,sol,R0,Cvn);
    so3_log(Cvb0,phi0,NULL);
    so3_jac(phi0,NULL,Jl);

    for (i=0;i<3;i++) {
        att[i]=sol->att[i]*D2R; /* convert degree to rad */
    }
    rpy2dcm(att,Cnb);
    matmul("NN",3,3,3,1.0,Cnb,Cvn,0.0,Cvb);
    so3_log(Cvb,phi,NULL);

    for (i=0;i<3;i++) {
        if (norm(phi0,3)&&fabs(v[nv]=phi[i]-phi0[i])>VARPOSE) {
            continue;
        }
        if (H) {
            for (j=0;j<3;j++) H[j+i*3]=-Jl[i+j*3];
        }
        ind[nv++]=i;
    }
    if (R) {
        for (i=0;i<nv;i++) {
            for (j=0;j<nv;j++) {
                R[i+j*nv]=R0[ind[i]+ind[j]*3];
            }
        }
    }
    trace(3,"v=\n"); tracemat(3,v, 1,nv,12,6);
    trace(3,"H=\n"); tracemat(3,H, 3,nv,12,6);
    trace(3,"R=\n"); tracemat(3,R,nv,nv,12,6);
    return nv;
}
/* find pose measurement data given solution time----------------------------*/
static int findpose(posebuf_t *posebuf,gtime_t time)
{
    int i;
    if (posebuf->dt==0.0&&posebuf->n>10) {
        for (i=0;i<10;i++) {
            posebuf->dt+=timediff(posebuf->data[i+1].time,posebuf->data[i].time);
        }
        posebuf->dt/=10.0;
    }
    for (i=posebuf->start;i<posebuf->end;i++) {
        if (fabs(timediff(posebuf->data[i].time,time))<DTTOL) {
            posebuf->start=i; return i;
        }
    }
    return -1;
}
/* initial covariance of estimate states-------------------------------------*/
static void initP(double *P,int nx)
{
    setzero(P,nx,nx); int i; for (i=0;i<3;i++) P[i+i*nx]=VARPOSE;
}
/* add system noise to covariance--------------------------------------------*/
static void addQ(double *P,int nx)
{
    int i; for (i=0;i<nx;i++) P[i+i*nx]+=SQR(5E-5*D2R);
}
/* build ekf filter matrix---------------------------------------------------*/
static int calibfilt(const insopt_t *opt,const solbuf_t *solbuf,
                     posebuf_t *posebuf,gtime_t ts,gtime_t te,
                     double *C)
{
    double *v,*H,*R,*P,*x,Cvb[9];
    double var,rpy[3];
    int i,j,k,info=0,nv=0;

    trace(3,"calibfilt:\n");

    v=zeros(3,1); H=zeros(3,3);
    R=zeros(3,3); P=zeros(3,3);
    x=zeros(1,3);

    matcpy(Cvb,C,3,3);
    initP(P,3);

    for (i=0;i<solbuf->n;i++) {

#if ADD_NOISE
        addQ(P,3);
#endif
        if (!screent(solbuf->data[i].time,ts,te,0.0)) continue;

        /* only use fixed solution */
        if (solbuf->data[i].stat!=SOLQ_FIX) {
            continue;
        }
        /* check velocity valid */
        if (norm(solbuf->data[i].rr+3,3)<MINVEL) {
            continue;
        }
        /* search pose measurement */
        if ((j=findpose(posebuf,solbuf->data[i].time))<0) {
            continue;
        }
        /* check variance valid */
        for (var=0.0,k=0;k<3;k++) {
            var+=posebuf->data[j].var[k];
        }
        if (var<=0.0||SQRT(var)>15.0*D2R) {
            continue;
        }
        /* measurement matrix */
        if (!(nv=build_ekf(opt,&solbuf->data[i],&posebuf->data[j],
                           Cvb,v,H,R))) {
            continue;
        }
        /* ekf filter */
        setzero(x,1,3);
        if (filter(x,P,H,v,R,3,nv)) {
            continue;
        }
        /* close loop for estimated states */
        corratt(x,Cvb);
    }
    /* check estimated states */
    if (norm(x,3)<1E-4) {

        matcpy(C,Cvb,3,3); /* output estimated states */

        dcm2rpy(Cvb,rpy);
        trace(3,"rpy=\n");
        tracemat(3,rpy,1,3,12,6);

        info=1;
    }
    else {
        seteye(C,3); info=0;
    }
    free(v); free(H); free(R);
    free(x); free(P);
    return info;
}
/* post process for calibration dual ant.------------------------------------*/
static int post_proc(const insopt_t *opt,const solbuf_t *solbuf,
                     posebuf_t *posebuf,gtime_t ts,gtime_t te,
                     const double *Cvb,double *ave,double *std)
{
    double Cbn[9],Cnb[9],Cvn[9],Cvn0[9],Rz[9],Ry[9];
    double Cnv[9],Cnv0[9],phi[3],phi0[3],att[3];
    double *v;
    int i,j,n=0;

    trace(3,"post_proc:\n");

    v=zeros(1,solbuf->n);
    posebuf->start=0;

    for (i=0;i<solbuf->n;i++) {
        if (!screent(solbuf->data[i].time,ts,te,0.0)) continue;

        /* only use fixed solution */
        if (solbuf->data[i].stat!=SOLQ_FIX) {
            continue;
        }
        /* check velocity valid */
        if (norm(solbuf->data[i].rr+3,3)<MINVEL) {
            continue;
        }
        for (j=0;j<3;j++) {
            att[j]=solbuf->data[i].att[j]*D2R;
        }
        rpy2dcm(att,Cnb);
        matt(Cnb,3,3,Cbn);

        matmul("NN",3,3,3,1.0,Cbn,Cvb,0.0,Cvn);

        /* search pose measurement */
        if ((j=findpose(posebuf,solbuf->data[i].time))<0) {
            continue;
        }
        if (norm(posebuf->data[j].var,3)<=0.0) {
            continue;
        }
        Ry(-posebuf->data[j].rpy[1],Ry);
        Rz(-posebuf->data[j].rpy[2],Rz);
        matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cvn0);

        matt(Cvn ,3,3,Cnv );
        matt(Cvn0,3,3,Cnv0);
        dcm2rpy(Cnv ,phi );
        dcm2rpy(Cnv0,phi0);

        v[n]=SQRT(SQR(phi[0]-phi0[0])+SQR(phi[1]-phi0[1])+SQR(phi[2]-phi0[2]));
        if (fabs(v[n])>VARPOSE) {
            continue;
        }
        n++;
    }
    if (std) stds(v,n);
    if (ave) {
        for (*ave=0,i=0;i<n;i++) *ave+=v[i];
        *ave=*ave/n;
    }
    free(v);
    return n/posebuf->n>0.6;
}
/* initial misalignment from pose measurement--------------------------------
 * args:  insopt_t *opt        I  ins options
 *        solbuf_t *solbuf     I  ins solution buffer data
 *        pose_meas_t *posebuf I  dual antennas pose measurement data
 *        gtime_t ts,te        I  start and end time fro process
 *        double *Cvb          O  misalignment from v-frame to b-frame
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
static int initmisali(const insopt_t *opt,const solbuf_t *solbuf,
                      posebuf_t *posebuf,gtime_t ts,gtime_t te,
                      double *Cvb)
{
    double *val,att[3]={0},Ry[9],Rz[9],Cvn[9];
    double Cvb0[9],Cnb[9],var;
    int i,j,k,n=0;

    trace(3,"initmisali:\n");

    val=zeros(3,MINNUMMEAS);

    for (i=0;i<solbuf->n;i++) {
        if (!screent(solbuf->data[i].time,ts,te,0.0)) continue;

        /* only use fixed solution */
        if (solbuf->data[i].stat!=SOLQ_FIX) {
            continue;
        }
        /* check velocity valid */
        if (norm(solbuf->data[i].rr+3,3)<MINVEL) {
            continue;
        }
        /* search pose measurement */
        if ((j=findpose(posebuf,solbuf->data[i].time))<0) {
            continue;
        }
        /* check variance valid */
        for (var=0.0,k=0;k<3;k++) {
            var+=posebuf->data[j].var[k];
        }
        if (var<=0.0||SQRT(var)>90.0*D2R) continue;
        Ry(-posebuf->data[j].rpy[1],Ry);
        Rz(-posebuf->data[j].rpy[2],Rz);
        matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cvn);

        for (j=0;j<3;j++) {
            att[j]=solbuf->data[i].att[j]*D2R;
        }
        rpy2dcm(att,Cnb);
        matmul("NN",3,3,3,1.0,Cnb,Cvn,0.0,Cvb0);

        dcm2rpy(Cvb0,&val[3*n++]);
        if (n>=MINNUMMEAS) break;
    }
    for (i=0;i<n;i++) {
        att[0]+=val[3*i+0];
        att[1]+=val[3*i+1];
        att[2]+=val[3*i+2];
    }
    for (i=0;i<3;i++) att[i]/=n;
    rpy2dcm(att,Cvb); free(val);
    return n==MINNUMMEAS;
}
/* calibration dual antennas misalignment to imu body------------------------
 * args:  insopt_t *opt       I  ins options
 *        solbuf_t *solbuf    I  ins solution buffer data
 *        posebuf_t *posebuf  I  dual antennas pose measurement data
 *        gtime_t ts,te       I  start and end time fro process
 *        double *Cvb         O  misalignment from v-frame to b-frame
 *        double *ave,*std    O  average and std of estimated misalignment
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int calibdualant(const insopt_t *opt,const solbuf_t *solbuf,
                        posebuf_t *posebuf,gtime_t ts,gtime_t te,
                        double *Cvb,double *ave,double *std)
{
    trace(3,"calibdualant:\n");

    if (solbuf->n<=0||posebuf->n<=0) return 0;
    rpy2dcm(opt->mis_euler,Cvb);

    /* initial misalignment */
    if (!initmisali(opt,solbuf,posebuf,ts,te,Cvb)) {

        trace(2,"initial misalignment fail\n");
        return 0;
    }
    /* calibration dual antennas misalignment */
    if (calibfilt(opt,solbuf,posebuf,ts,te,Cvb)) {

        /* post process for calibration */
        if (!post_proc(opt,solbuf,posebuf,ts,te,Cvb,ave,std)) {

            trace(2,"post process fail\n");
            return 0;
        }
        return 1;
    }
    else {
        trace(2,"calibration dual ant. fail\n");
        return 0;
    }
}
/* calibration dual antennas misalignment to imu body------------------------
 * args:  insopt_t *opt       I  ins options
 *        solbuf_t *solbuf    I  ins solution buffer data
 *        posebuf_t *posebuf  I  dual antennas pose measurement data
 *        gtime_t ts,te       I  start and end time fro process
 *        double *Cvb         O  misalignment from v-frame to b-frame 
 *        double *std         O  average and std of estimated misalignment
 * note: this method isn't based on ekf filter
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int calibdualantx(const insopt_t *opt,const solbuf_t *solbuf,
                         posebuf_t *posebuf,gtime_t ts,gtime_t te,
                         double *Cvb,double *std)
{
    double *val,Cvb0[9],phi[3],Ry[9],Rz[9],var;
    double Cvn[9],att[3],Cnb[9];
    int i,j,k,n=0;

    trace(3,"calibdualantx:\n");

    if (solbuf->n<=0||posebuf->n<=0) return 0;
    rpy2dcm(opt->mis_euler,Cvb);

    val=zeros(3,posebuf->n);
    for (i=0;i<solbuf->n;i++) {
        if (!screent(solbuf->data[i].time,ts,te,0.0)) continue;

        /* only use fixed solution */
        if (solbuf->data[i].stat!=SOLQ_FIX) {
            continue;
        }
        /* check velocity valid */
        if (norm(solbuf->data[i].rr+3,3)<MINVEL) {
            continue;
        }
        /* search pose measurement */
        if ((j=findpose(posebuf,solbuf->data[i].time))<0) {
            continue;
        }
        /* check variance valid */
        for (var=0.0,k=0;k<3;k++) {
            var+=posebuf->data[j].var[k];
        }
        if (var<=0.0||SQRT(var)>10.0*D2R) {
            continue;
        }
        Ry(-posebuf->data[j].rpy[1],Ry);
        Rz(-posebuf->data[j].rpy[2],Rz);
        matmul("NN",3,3,3,1.0,Rz,Ry,0.0,Cvn);

        for (k=0;k<3;k++) {
            att[k]=solbuf->data[i].att[k]*D2R;
        }
        rpy2dcm(att,Cnb);
        matmul("NN",3,3,3,1.0,Cnb,Cvn,0.0,Cvb0);

        dcm2rpy(Cvb0,phi);
        so3_log(Cvb0,phi,NULL);

        val[3*n+0]=phi[0];
        val[3*n+1]=phi[1];
        val[3*n+2]=phi[2]; n++;
    }
    setzero(phi,1,3);
    for (i=0;i<n;i++) {
        for (j=0;j<3;j++) phi[j]+=val[3*i+j];
    }
    if (n) {
        phi[0]/=n; phi[1]/=n; phi[2]/=n;
        so3_exp(phi,Cvb);
    }
    /* post process for calibration */
    if (!post_proc(opt,solbuf,posebuf,ts,te,Cvb,NULL,std)) {

        free(val);
        trace(2,"post process fail\n");
        return 0;
    }
    free(val);
    return 1;
}
