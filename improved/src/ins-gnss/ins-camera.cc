/*------------------------------------------------------------------------------
* ins-camera.cc : camera common functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*    [5] Weiss,Real-Time Metric State Estimation for Modular Vision-Inertial
*        Systems
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/05/06 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* distorts a point on the unit plane (in camera coordinates) according to the
 * radtan distortion model
 * args  :  cam_t *cam  I  camera model
 *          double *in  I  undistorted point coordinates on the unit plane
 *                         (in camera coordinates)
 *          double *out O  distorted point coordinates on the unit plane
 *                         (in camera coordinates)
 *          double *J   O  jacobian matrix of the distortion process(NULL: on output)
 * return: none
 * --------------------------------------------------------------------------*/
extern void distortradtan(const cam_t *cam,const double *in,double *out,double *J)
{
    double mx2,my2,mxy,r2,rad;

    mx2=in[0]*in[0];
    my2=in[1]*in[1];
    mxy=in[0]*in[1];
    r2=mx2+my2;
    rad=cam->k1*r2+cam->k2*r2*r2;

    if (out) {
        out[0]=in[0]+in[0]*rad+2.0*cam->p1*mxy+cam->p2*(r2+2.0*mx2);
        out[1]=in[1]+in[1]*rad+2.0*cam->p2*mxy+cam->p1*(r2+2.0*my2);
    }
    if (J) {
        J[0]=1.0+rad+2.0*cam->k1*mx2+cam->k2*r2*4.0*mx2+
             2.0*cam->p1*in[1]+6.0*cam->p2*in[0];
        J[1]=2.0*cam->k1*mxy+4.0*cam->k2*r2*in[0]*in[1]+
             2.0*cam->p1*in[0]+2.0*cam->p2*in[1];
        J[2]=J[1];
        J[3]=1.0+rad+2.0*cam->k1*my2+cam->k2*r2*4.0*my2+
             6.0*cam->p1*in[1]+2.0*cam->p2*in[0];
    }
}
/* apply undistortion to recover a point in the normalized image plane-------
 * args  :  cam_t *cam  I  camera model
 *          double *in  I  the distorted point. after the function, this point
 *                         is in the normalized image plane
 *          double *out O  undistorted point coordinates on the unit plane
 *          double *J   O  the jacobian of the undistortion function with
 *                         respect to small changes in the input point
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int undistortradtan(const cam_t *cam,const double *in,double *out,
                           double *J)
{
    double ybar[2],ytmp[2],F[4],e[2],du[2];
    double T[4],E[2];
    int i,flag=0;

    matcpy(ybar,in,1,2);
    for (i=0;i<5;i++) { /* iterations for undistorted */

        matcpy(ytmp,ybar,1,3); distortradtan(cam,ytmp,ytmp,F);

        for (i=0;i<2;i++) {
            e[i]=in[i]-ytmp[i];
        }
        matmul("TN",2,2,2,1.0,F,F,0.0,T);
        matmul("TN",2,1,2,1.0,F,e,0.0,E);
        if (matinv(T,2)) {
            continue;
        }
        matmul("NN",2,1,2,1.0,T,E,0.0,du);

        for (i=0;i<2;i++) ybar[i]+=du[i];
        if (norm(e,2)<1E-8) {
            flag=1;
            break;
        }
    }
    if (out) {
        matcpy(out,ybar,1,2);
    }
    if (J) {
        distortradtan(cam,out,NULL,J);
        matinv(J,2);
    }
    return flag;
}