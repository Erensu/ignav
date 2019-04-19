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
#include "carvig.h"

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
/* so3 jacobian matrix--------------------------------------------------------*/
extern void so3jac(const double *phi,double *Jri)
{
    double I[9]={1,0,0,0,1,0,0,0,1},W[9],W2[9];
    double n=norm(phi,3);
    int i;

    if (n<=1E-8) {seteye(Jri,3); return;}
    skewsym3(phi,W);
    matmul("NN",3,3,3,1.0,W,W,0.0,W2);

    for (i=0;i<9;i++) {
        Jri[i]=I[i]+0.5*W[i]+(1.0/SQR(n)+(1.0+cos(n))/(2.0*n*sin(n)))*W2[i];
    }
}
/* so3 exp map----------------------------------------------------------------*/
extern void so3exp(const double *phi,double *C)
{
    double wx,wy,wz,theta=norm(phi,3);
    wx=phi[0]/theta;
    wy=phi[1]/theta;
    wz=phi[2]/theta;
    if (fabs(theta)<1E-7) {seteye(C,3); return;}
    double wwTxx=wx*wx,wwTyy=wy*wy,wwTzz=wz*wz;
#ifndef NDEBUG
    double l_n=wwTxx+wwTyy+wwTzz;
    if (fabs(l_n-1.0)>1E-9) {
        fprintf(stderr,"rodriguez: length of n should be 1\n");
        return;
    }
#endif
    double c=cos(theta),s=sin(theta),c_1=1-c;
    double swx=wx*s,swy=wy*s,swz=wz*s;
    double C00=c_1*wwTxx,C01=c_1*wx*wy,C02=c_1*wx*wz;
    double               C11=c_1*wwTyy,C12=c_1*wy*wz;
    double                             C22=c_1*wwTzz;

    C[0]=   c+C00; C[3]=-swz+C01; C[6]= swy+C02;
    C[1]= swz+C01; C[4]=   c+C11; C[7]=-swx+C12;
    C[2]=-swy+C02; C[5]= swx+C12; C[8]=   c+C22;
}
/* so3 log map----------------------------------------------------------------*/
extern void so3log(const double *C,double *phi)
{
    double tr=C[0]+C[4]+C[8],magnitude,theta,tr_3=tr-3.0;

    if (fabs(tr+1.0)<1E-10) {
        if (fabs(C[8]+1.0)>1E-10) {
            phi[0]=(PI/(SQRT(2.0+2.0*C[8])))*C[6];
            phi[1]=(PI/(SQRT(2.0+2.0*C[8])))*C[7];
            phi[2]=(PI/(SQRT(2.0+2.0*C[8])))*(1.0+C[8]);
        }
        else if (fabs(C[4]+1.0)>1E-10) {
            phi[0]=(PI/(SQRT(2.0+2.0*C[4])))*C[3];
            phi[1]=(PI/(SQRT(2.0+2.0*C[4])))*(1.0+C[4]);
            phi[2]=(PI/(SQRT(2.0+2.0*C[4])))*C[5];
        }
        else {
            phi[0]=(PI/(SQRT(2.0+2.0*C[0])))*(1.0+C[0]);
            phi[1]=(PI/(SQRT(2.0+2.0*C[0])))*C[1];
            phi[2]=(PI/(SQRT(2.0+2.0*C[0])))*C[2];
        }
    }
    else {
        if (tr_3<-1E-7) {
            theta=acos((tr-1.0)/2.0);
            magnitude= theta/(2.0*sin(theta));
        }
        else {
            magnitude=0.5-tr_3*tr_3/12.0;
        }
        phi[0]=magnitude*(C[5]-C[7]);
        phi[1]=magnitude*(C[6]-C[2]);
        phi[2]=magnitude*(C[1]-C[3]);
    }
}