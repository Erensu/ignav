/*------------------------------------------------------------------------------
* so3.cc : Lie Groups functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Ethan Eade, Lie Groups for Computer Vision
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/03/15 1.0 new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants -----------------------------------------------------------------*/
#define SMALL_EPS         1E-10
#define quatnorm(q)       quat_normalize_self(q)

/* type define----------------------------------------------------------------*/
typedef struct {          /* so3 struct data type */
    double unitq[4];      /* unit quaternion for handling so3 */
} so3_t;

/* set quaternion-------------------------------------------------------------*/
static void so3_setq(so3_t *so3,const double *q)
{
    so3->unitq[0]=q[0]; /* w */
    so3->unitq[1]=q[1]; /* x */
    so3->unitq[2]=q[2]; /* y */
    so3->unitq[3]=q[3]; /* z */
}
/* initial so3----------------------------------------------------------------*/
static void so3_init(so3_t *so3)
{
    so3->unitq[0]=1.0; so3->unitq[1]=so3->unitq[2]=so3->unitq[3]=0.0;
}
/* so3 exponential map interface----------------------------------------------*/
static void exp_theta(so3_t *so3,const double *omega,double *theta)
{
    double half_theta,imag_factor,real_factor;
    if (theta) {
        *theta=norm(omega,3); half_theta=0.5**theta;
    }
    else return;

    real_factor=cos(half_theta);
    if ((*theta)<SMALL_EPS) {
        double theta_sq =(*theta)*(*theta);
        double theta_po4=theta_sq*theta_sq;
        imag_factor=0.5-0.0208333*theta_sq+0.000260417*theta_po4;
    }
    else {
        double sin_half_theta=sin(half_theta);
        imag_factor=sin_half_theta/(*theta);
    }
    so3->unitq[0]=real_factor;
    so3->unitq[1]=imag_factor*omega[0];
    so3->unitq[2]=imag_factor*omega[1];
    so3->unitq[3]=imag_factor*omega[2];
}
/* so3 log map interface------------------------------------------------------*/
static void log_theta(so3_t *so3,double *omega,double *theta)
{
    double tmp[3]; matcpy(tmp,so3->unitq+1,1,3);
    double n=norm(tmp,3);
    double w=so3->unitq[0];
    double sw=w*w;
    double ta;
    
    if (n<SMALL_EPS) {
        /* if quaternion is normalized and n=1, then w should be 1;
         * w=0 should never happen here */
        ta=2.0/w-2.0*(n*n)/(w*sw);
    }
    else {
        if (fabs(w)<SMALL_EPS) {
            if (w>0) ta=M_PI/n;
            else {
                ta=-M_PI/n;
            }
        }
        ta=2*atan(n/w)/n;
    }
    if (theta) *theta=ta*n;
    omega[0]=ta*so3->unitq[1];
    omega[1]=ta*so3->unitq[2];
    omega[2]=ta*so3->unitq[3];
}
/* so3 exponential map--------------------------------------------------------*/
extern void so3_exp(const double *omega,double *C)
{
    so3_t so3; quat_t q;
    double theta=0.0;

    so3_init(&so3); exp_theta(&so3,omega,&theta);
    q.w=so3.unitq[0]; q.x=so3.unitq[1];
    q.y=so3.unitq[2]; q.z=so3.unitq[3];
    quat2dcm(&q,C);
}
/* so3 logarithm map----------------------------------------------------------*/
extern void so3_log(const double *C,double *omega,double *theta)
{
    so3_t so3; quat_t q;

    so3_init(&so3); dcm2quat(C,&q);
    quatnorm(&q);
    so3.unitq[0]=q.w; so3.unitq[1]=q.x;
    so3.unitq[2]=q.y; so3.unitq[3]=q.z;
    log_theta(&so3,omega,theta);
}
/* lie bracket----------------------------------------------------------------*/
extern void liebracket(const double *omega1,const double *omega2,double *omega)
{
    cross3(omega1,omega2,omega);
}
/* hat operation--------------------------------------------------------------*/
extern void so3_hat(const double *omega,double *Omg)
{
    skewsym3(omega,Omg);
}
/* vee operation--------------------------------------------------------------*/
extern void so3_vee(const double *Omg,double *omega)
{
    assert(fabs(Omg[5]+Omg[7])<SMALL_EPS);
    assert(fabs(Omg[6]+Omg[2])<SMALL_EPS);
    assert(fabs(Omg[1]+Omg[3])<SMALL_EPS);
    omega[0]=Omg[5]; omega[1]=Omg[6]; omega[2]=Omg[1];
}