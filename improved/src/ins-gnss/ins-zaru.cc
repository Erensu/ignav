/*------------------------------------------------------------------------------
* ins-zaru.cc : zero angular rate update for ins navigation
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
#define MAXVEL      0.1               /* max velocity for using zero angular rate update */
#define MAXGYRO     (5.0*D2R)         /* max rotation speed for zero angular rate update */
#define VARARE      SQR(1.0*D2R)      /* initial variance of angular rate (rad/s)^2) */
#define MINZAC      100               /* min count for zero velocity update once */

/* zero angular rate update for ins navigation---------------------------------
 * args   :  insstate_t *ins  IO  ins state
 *           insopt_t *opt    I   ins options
 *           imudata_t *imu   I   imu measurement data
 *           int flag         I   static flag
 * return :  1 (ok) or 0 (fail)
 * ---------------------------------------------------------------------------*/
extern int zaru(insstate_t *ins,const insopt_t *opt,const imud_t *imu,int flag)
{
    int info=0,nx=ins->nx;
    static int nz=0;
    double *v,*x,*H,*R,I[9]={-1,0,0,0,-1,0,0,0,-1};

    trace(3,"zaru:\n");

    flag&=nz++>MINZAC?nz=0,true:false;

    if (flag==0||opt->bgopt!=INS_BGEST) return 0;

    x=zeros(1,nx); v=zeros(3,1);
    H=zeros(3,nx); R=zeros(3,3);

    /* sensitive matrix */
    asi_blk_mat(H,3,nx,I,3,3,0,xiBg(opt));

    /* variance matrix */
    R[0]=R[4]=R[8]=VARARE;

    v[0]=-imu[0].gyro[0];
    v[1]=-imu[0].gyro[1];
    v[2]=-imu[0].gyro[2]; /* residual vector */

    if (norm(v,3)<MAXGYRO&&norm(ins->ve,3)<MAXVEL) {

        /* ekf filter */
        info=filter(x,ins->P,H,v,R,nx,3);

        /* solution fail */
        if (info) {
            trace(2,"zero velocity update filter error \n");
            info=0;
        }
        else {
            /* close loop for estimate state */
            ins->stat=INSS_ZARU;
            info=1;
            clp(ins,opt,x);
            trace(3,"zero angular rate update ok\n");
        }
    }
    free(x); free(H);
    free(R); free(v);
    return info;
}