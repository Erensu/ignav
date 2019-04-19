/*------------------------------------------------------------------------------
* se3.cc : Lie Groups functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Ethan Eade, Lie Groups for Computer Vision
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/05/28 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define SMALL_EPS         1E-10
#define quatnorm(q)       quat_normalize_self(q)

/* type define----------------------------------------------------------------*/
typedef struct {          /* so3 struct data type */
    double unitq[4];      /* unit quaternion for handling so3 */
    double t[3];          /* translation */
} se3_t;

/* set quaternion-------------------------------------------------------------*/
static void se3_setq(se3_t *se3,const double *q)
{
    se3->unitq[0]=q[0]; /* w */
    se3->unitq[1]=q[1]; /* x */
    se3->unitq[2]=q[2]; /* y */
    se3->unitq[3]=q[3]; /* z */
}
/* set Lie Groups SE(3)-------------------------------------------------------*/
static void set_se3(const double *R,const double *t,se3_t *se3)
{
    quat_t q;
    dcm2quat(R,&q);

    quatnorm(&q);
    se3->unitq[0]=q.w; se3->unitq[1]=q.x;
    se3->unitq[2]=q.y; se3->unitq[3]=q.z;

    se3->t[0]=t[0];
    se3->t[1]=t[1];
    se3->t[2]=t[2];
}
/* se3 log--------------------------------------------------------------------
 * args:    double *R      I  transformation matrix
 *          double *t      I  translation
 *          double *omega  O  6x1 vector
 * return: none
 * ---------------------------------------------------------------------------*/
extern void se3_log(const double *R,const double *t,double *omega)
{
    static double I[9]={1,0,0,0,1,0,0,0,1};
    double theta,Omega[9],Vi[9];
    int i,j;

    trace(3,"se3_log:\n");

    so3_log(R,omega+3,&theta);

    so3_hat(omega+3,Omega);
    matmul("NN",3,3,3,1.0,Omega,Omega,0.0,Vi);

    if (theta<SMALL_EPS) {

        for (i=0;i<3;i++) for (j=0;j<3;j++) {
            Vi[i+3*j]=I[i+3*j]-0.5*Omega[i+3*j]+1.0/12.0*Vi[i+3*j];
        }
        matmul("NN",3,1,3,1.0,Vi,t,0.0,omega);
    }
    else {

        for (i=0;i<3;i++) for (j=0;j<3;j++) {
            Vi[i+3*j]=I[i+3*j]-0.5*Omega[i+3*j]+
                    (1.0-theta/(2.0*tan(theta/2.0)))/(theta*theta)
                    *Vi[i+3*j];
        }
        matmul("NN",3,1,3,1.0,Vi,t,0.0,omega);
    }
    return;
}
/* se3 exp--------------------------------------------------------------------
 * args:    double *omg  I  6x1 vector
 *          double *R    O  transformation matrix
 *          double *t    O  translation
 * return: none
 * ---------------------------------------------------------------------------*/
extern void se3_exp(const double *omg,double *R,double *t)
{
    static double I[9]={1,0,0,0,1,0,0,0,1};
    double Omega[9],V[9];
    double Omega_sq[9],theta,theta_sq,a,b;
    int i,j;

    trace(3,"se3_exp:\n");

    so3_exp(omg+3,    R);
    so3_hat(omg+3,Omega);
    matmul("NN",3,3,3,1.0,Omega,Omega,0.0,Omega_sq);

    theta=norm(omg+3,3);

    if (theta<SMALL_EPS) {
        matcpy(V,R,3,3);
    }
    else {
        theta_sq=SQR(theta);

        b=(theta-sin(theta))/(theta_sq*theta);
        a=(1.0  -cos(theta))/(theta_sq);

        for (i=0;i<3;i++) for (j=0;j<3;j++) {
            V[i+3*j]=I[i+3*j]+a*Omega[i+3*j]+b*Omega_sq[i+3*j];
        }
    }
    matmul("NN",3,1,3,1.0,V,omg,0.0,t);
    return;
}
/* se3 hat--------------------------------------------------------------------*/
extern void se3_hat(const double *omg,double *Omg)
{
    double Omg_[9];
    int i,j;

    trace(3,"se3_hat:\n");

    setzero(Omg,4,4);

    so3_hat(omg+3,Omg_);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        Omg[4*i+j]=Omg_[3*i+j];
    }
    Omg[4*3+0]=omg[0];
    Omg[4*3+1]=omg[1];
    Omg[4*3+2]=omg[2];
}
/* se3 vee--------------------------------------------------------------------*/
extern void se3_vee(const double *Omg,double *omg)
{
    double T[9];
    int i,j;

    trace(3,"se3_vee:\n");

    omg[0]=Omg[4*3+0];
    omg[1]=Omg[4*3+1];
    omg[2]=Omg[4*3+2];

    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        T[3*i+j]=Omg[4*i+j];
    }
    so3_log(T,omg+3,NULL);
}

