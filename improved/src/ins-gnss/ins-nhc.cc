/*-----------------------------------------------------------------------------
* ins-nhc.cc : non-holonomic constraint for ins navigation
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
* history : 2017/11/11 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants ----------------------------------------------------------------*/
#define MAXVEL      0.5          /* max velocity for using non-holonomic constraint */
#define VARVEL      SQR(0.05)    /* initial variance of ins velocity ((m/s)^2) */

/* jacobian of perturb rotation wrt. perturb euler angles--------------------*/
static void jacobian_prot_pang(const double *Cbe,double *S)
{
    double rpy[3]={0};
    c2rpy(Cbe,rpy);
    S[0]= cos(rpy[1])*cos(rpy[2]); S[3]=sin(rpy[2]); S[6]=0.0;
    S[1]=-cos(rpy[1])*sin(rpy[2]); S[4]=cos(rpy[2]); S[7]=0.0;
    S[2]= sin(rpy[1]);             S[5]=0.0;         S[8]=1.0;
}
/* measurement sensitive-matrix for non-holonomic----------------------------*/
static int bldnhc(const insopt_t *opt,const imud_t *imu,const double *Cbe,
                  const double *ve,int nx,double *v,double *H,double *R)
{
    int i,nv,IA,IV;
    double C[9],T[9],vb[3],r[2],S[9];
    double T1[9];

    trace(3,"bldnhc:\n");

    IA=xiA(opt); IV=xiV(opt);

    /* velocity in body-frame */
    matmul("TN",3,1,3,1.0,Cbe,ve,0.0,vb);

    skewsym3(ve,C);
    matmul("TN",3,3,3,1.0,Cbe,C,0.0,T);
    matt(Cbe,3,3,C);

#if UPD_IN_EULER
    jacobian_prot_pang(Cbe,S);
    matcpy(T1,T,3,3);
    matmul("NN",3,3,3,1.0,T1,S,0.0,T);
#endif
    /* build residual vector */
    for (nv=0,i=1;i<3;i++) {

        /* check velocity measurement */
        if (fabs(vb[i])>MAXVEL) {
            trace(2,"too large velocity measurement\n");
            continue;
        }
        /* check gyro measurement */
        if (fabs(norm(imu->gyro,3))>30.0*D2R) {
            trace(2,"too large vehicle turn\n");
            continue;
        }
        H[IA+nv*nx]=T[i]; H[IA+1+nv*nx]=T[i+3]; H[IA+2+nv*nx]=T[i+6];
        H[IV+nv*nx]=C[i]; H[IV+1+nv*nx]=C[i+3]; H[IV+2+nv*nx]=C[i+6];
        
        v[nv  ]=vb[i];
        r[nv++]=VARVEL;
    }
    for (i=0;i<nv;i++) R[i+i*nv]=r[i];
    return nv;
}
/* close loop for non-holonomic constraint-----------------------------------*/
extern void clp(insstate_t *ins,const insopt_t *opt,const double *x)
{
    int i,iba,nba,ibg,nbg,isg,nsg,isa,nsa;
    double *I=eye(3),fibc[3],omgbc[3],ang[3];

    iba=xiBa(opt); nba=xnBa(opt);
    ibg=xiBg(opt); nbg=xnBg(opt);
    isg=xiSg(opt); nsg=xnSg(opt);
    isa=xiSa(opt); nsa=xnSa(opt);

#if UPD_IN_EULER
    c2rpy(ins->Cbe,ang);
    for (i=0;i<3;i++) ang[i]=ang[i]+x[i];
    rpy2c(ang,ins->Cbe);
#else
    /* close-loop attitude correction */
    if (x[0]!=DISFLAG) corratt(x,ins->Cbe);
#endif

    /* close-loop velocity and position correction */
    ins->ve[0]-=x[xiV(opt)+0];
    ins->ve[1]-=x[xiV(opt)+1];
    ins->ve[2]-=x[xiV(opt)+2];

    ins->re[0]-=x[xiP(opt)+0];
    ins->re[1]-=x[xiP(opt)+1];
    ins->re[2]-=x[xiP(opt)+2];

    /* close-loop accl and gyro bias */
    if (nba&&x[iba]!=DISFLAG) {
        ins->ba[0]+=x[iba+0];
        ins->ba[1]+=x[iba+1];
        ins->ba[2]+=x[iba+2];
    }
    if (nbg&&x[ibg]!=DISFLAG) {
        ins->bg[0]+=x[ibg+0];
        ins->bg[1]+=x[ibg+1];
        ins->bg[2]+=x[ibg+2];
    }
    /* close-loop residual scale factors of gyro and accl */
    if (nsg&&x[isg]!=DISFLAG) {
        for (i=isg;i<isg+nsg;i++) ins->Mg[i-isg+(i-isg)*3]+=x[i];
    }
    if (nsa&&x[isa]!=DISFLAG) {
        for (i=isa;i<isa+nsa;i++) ins->Ma[i-isa+(i-isa)*3]+=x[i];
    }
#if 1
    /* correction imu accl and gyro measurements data */
    ins_errmodel2(ins->fb0,ins->omgb0,
                  ins->Ma,ins->Mg,
                  ins->ba,ins->bg,ins->Gg,
                  fibc,omgbc);

    /* correction imu-body accelerometer */
    getaccl(fibc,ins->Cbe,ins->re,ins->ve,ins->ae);
#endif
    free(I);
}
/* using non-holonomic constraint for ins navigation---------------------------
 * args    :  insstate_t* ins  IO  ins state
 *            insopt_t* opt    I   ins options
 *            imud_t* imu      I   imu measurement data
 * return  : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int nhc(insstate_t *ins,const insopt_t *opt,const imud_t *imu)
{
    int nx=ins->nx,info=0,nv;
    double *H,*v,*R,*x;

    trace(3,"nhc:\n");

    H=zeros(2,nx); R=zeros(2,2);
    v=zeros(2,1); x=zeros(1,nx);

    nv=bldnhc(opt,imu,ins->Cbe,ins->ve,nx,v,H,R);
    if (nv>0) {

        /* kalman filter */
        info=filter(x,ins->P,H,v,R,nx,nv);

        /*  check ok? */
        if (info) {
            trace(2,"non-holonomic constraint filter fail\n");
            info=0;
        }
        else {
            /* solution ok */
            ins->stat=INSS_NHC; info=1;
            clp(ins,opt,x);
            trace(3,"use non-holonomic constraint ok\n");
        }
    }
    free(H); free(v);
    free(R); free(x);
    return info;
}