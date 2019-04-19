/*------------------------------------------------------------------------------
* ins-zvu.cc : zero velocity update for ins navigation
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
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define MAXVEL      0.1           /* max velocity for using non-holonomic constraint */
#define MAXGYRO     (10.0*D2R)    /* max rotation speed for using non-holonomic constraint */
#define VARVEL      SQR(0.05)     /* initial variance of receiver vel ((m/s)^2) */
#define MINZC       15            /* min count for zero velocity update once */

/* detect static imu measurement---------------------------------------------*/
extern int detstc(const imud_t *imu,int n,const insopt_t *opt,const double *pos)
{
    switch (opt->detst) {
        case IMUDETST_GLRT: return detstatic_GLRT(imu,n,opt,pos);
        case IMUDETST_MV  : return detstatic_MV  (imu,n,opt);
        case IMUDETST_MAG : return detstatic_MAG (imu,n,opt,pos);
        case IMUDETST_ARE : return detstatic_ARE (imu,n,opt);
        case IMUDETST_ALL :
            return detstatic_GLRT(imu,n,opt,pos)&&detstatic_MV (imu,n,opt)&&
                   detstatic_MAG (imu,n,opt,pos)&&detstatic_ARE(imu,n,opt);
        default:
            return detstatic_GLRT(imu,n,opt,pos);
    }
}
/* zero velocity update for ins navigation -----------------------------------
 * args    :  insstate_t *ins  IO  ins state
 *            insopt_t *opt    I   ins options
 *            imud_t *imu      I   imu measurement data
 *            int flag         I   static flag (1: static, 0: motion)
 * return  : 1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int zvu(insstate_t *ins,const insopt_t *opt,const imud_t *imu,int flag)
{
    int nx=ins->nx,info=0;
    static int nz=0;
    double *x,*H,*R,*v,I[9]={-1,0,0,0,-1,0,0,0,-1};

    trace(3,"zvu:\n");

    flag&=nz++>MINZC?nz=0,true:false;

    if (!flag) return info;

    x=zeros(1,nx); H=zeros(3,nx);
    R=zeros(3,3); v=zeros(3,1);

    /* sensitive matrix */
    asi_blk_mat(H,3,nx,I,3,3,0,3);

    /* variance matrix */
    R[0]=R[4]=R[8]=VARVEL;

    v[0]=ins->ve[0];
    v[1]=ins->ve[1];
    v[2]=ins->ve[2]; /* residual vector */

    if (norm(v,3)<MAXVEL&&norm(imu->gyro,3)<MAXGYRO) {

        /* ekf filter */
        info=filter(x,ins->P,H,v,R,nx,3);

        /* solution fail */
        if (info) {
            trace(2,"zero velocity update filter error\n");
            info=0;
        }
        else {
            /* solution ok */
            ins->stat=INSS_ZVU;
            info=1;
            clp(ins,opt,x);
            trace(3,"zero velocity update ok\n");
        }
    }
    free(x); free(H);
    free(R); free(v);
    return info;
}

