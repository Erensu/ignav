/*------------------------------------------------------------------------------
* ins-sim.cc : ins simulation functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/10/13 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* generates a uniform distributed random number between min and max --------*/
static double getuniform(double min,double max)
{
    return 1.0*rand()/RAND_MAX*(max-min)+min;
}
/* creates gaussian distributed random numbers (Box-Müller method)-----------*/
static double getgaussian(double std)
{
    if (std<0.0) std=-std;
    double x1,x2,w,y1;
    do {
        x1=getuniform(-1.0,1.0); x2=getuniform(-1.0,1.0); w=x1*x1+x2*x2;
    } while (w>=1.0);
    w=sqrt((-2.0*log(w))/w); y1=x1*w; return std*y1;
}
/* imu sensor outputs on static base simulation-------------------------------
 * args    :  double *rpy       I  initial attitude (roll,picth,yaw)
 *            double *pos       I  ins position in ned-frame (rad,m)
 *            double ts         I  sampling interval (s)
 *            double T          I  total sampling simulation time
 *            imu_err_t *err    I  imu error model
 *            imu_t *data       O  output imu measurements (rad/s,m/s^2)
 * return : numbers of imu measurement data
 * --------------------------------------------------------------------------*/
extern int sim_imu_static(const double *rpy,const double *pos,const double ts,
                          const double T,const imu_err_t *err,imu_t *data)
{
    int i,j,n;
    double Cnb[9],wnie[3]={0},wbie[3]={0},gn[3],gb[3];
    gtime_t t0={0};

    rpy2dcm(rpy,Cnb);

    /* convert earth rotation to body-frame */
    wnie[0]=OMGE*cos(pos[0]); wnie[2]=-OMGE*sin(pos[0]);
    matmul("TN",3,1,3,1.0,Cnb,wnie,0.0,wbie);

    /* gravity in ned-frame and body-frame */
    gravity_ned(pos,gn);
    matmul("TN",3,1,3,1.0,Cnb,gn,0.0,gb);

    n=ROUND(T/ts);
    data->n=data->nmax=n;
    data->data=(imud_t*)malloc(sizeof(imud_t)*n);

    /* generate imu static simulation measurements */
    for (i=0;i<n;i++) {
        for (j=0;j<3;j++) {
            data->data[i].time=timeadd(t0,i*ts);
            data->data[i].gyro[j]=wbie[j]+err->bg[j]+err->wgn[j]/sqrt(ts)*getgaussian(1.0);
            data->data[i].accl[j]=-gn [j]+err->ba[j]+err->wan[j]/sqrt(ts)*getgaussian(1.0);
        }
    }
    data->coor=IMUCOOR_FRD; data->decfmt=IMUDECFMT_RATE;
    return n;
}