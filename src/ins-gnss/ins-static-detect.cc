/*------------------------------------------------------------------------------
* ins-static-dector.cc : ins detect static imu measurement functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/11/03 1.0 new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* add a span time -----------------------------------------------------------*/
static int addspan(spana_t *span,const span_t *data)
{
    span_t *span_data;

    if (span->nmax<=span->n) {
        if (span->nmax<=0) span->nmax=64; else span->nmax*=2;
        if (!(span_data=(span_t *)realloc(span->tt,sizeof(span_t)*span->nmax))) {
            trace(1,"addspan: memalloc error n=%dx%d\n",sizeof(span_t),span->nmax);
            free(span->tt); span->tt=NULL; span->n=span->nmax=0;
            return -1;
        }
        span->tt=span_data;
    }
    span->tt[span->n++]=*data;
    return 1;
}
/* static detector for per accl axis,this function check accl variance--------*/
static int detstaticaxis(const imud_t *imu,int n,int axis,double thres)
{
    int i;
    double mean,var,*a=mat(n,1);

    for (mean=0.0,i=0;i<n;i++) mean+=imu[i].accl[axis]; mean/=n;
    for (i=0;i<n;i++) a[i]=imu[i].accl[axis]-mean;
    
    matmul("NT",1,1,n,1.0/(n-1),a,a,0.0,&var);

    /* detect static for imu measurement data */
    for (i=0;i<n;i++) {
        if (fabs(a[i]/SQRT(var))<thres) continue; free(a); return 0;
    }
    free(a); return 1;
}
/* static detector for gravity in ned-frame----------------------------------*/
static int detstaticg(const insstate_t *ins,const imud_t *imu,int n,double thres)
{
    int i; double gn[3];

    gravity_ned(ins->rn,gn);
    for (i=0;i<n;i++) {
        if (fabs(norm(imu[i].accl,3)-norm(gn,3))<thres) continue;
        return 0;
    }
    return 1;
}
/* static detector by gyro measurement --------------------------------------*/
static int detstaticgyro(const imud_t *imu,int n,const double *thres)
{
    int i; for (i=0;i<n;i++) {
        if (fabs(imu[i].gyro[0])<thres[0]&&
            fabs(imu[i].gyro[1])<thres[1]&&
            fabs(imu[i].gyro[2])<thres[2]) continue;
        return 0;
    }
    return 1;
}
/* static smooth window detector ----------------------------------------------
 * args   :   imud_t* imu      I  imu measurement data
 *            int n            I  number of imu measurement data
 *            insopt_t *insopt I  thres of static detect
 * return : 1: static,0:motion
 * ---------------------------------------------------------------------------*/
extern int detstaticw(const insstate_t *ins,const imud_t *imu,int n,
                      const insopt_t *opt)
{
    return detstaticaxis(imu,n,0,opt->zvopt.athres[0])&&
           detstaticaxis(imu,n,1,opt->zvopt.athres[1])&&
           detstaticaxis(imu,n,2,opt->zvopt.athres[2])&&
           detstaticgyro(imu,n,opt->zvopt.gyrothres)&&
           detstaticg(ins,imu,n,opt->zvopt.gthres);
}
/* static imu measurement data detector---------------------------------------
 * args   :  imu_t *imu       I  imu measurement data
 *           insstate_t *ins  I  ins state
 *           insopt_t *opt    I  ins options
 *           spana_t *span    O  static imu data time internal
 *           double ts        I  time span of static detect (0.0: detect all)
 * return : number of detect static imu measurement internal
 * --------------------------------------------------------------------------*/
extern int detstatic(const imu_t *imu,const insstate_t *ins,const insopt_t *opt,
                     spana_t *span,double ts)
{
    int i=0,j,k,flag=0,ws=opt->zvopt.ws+1;
    gtime_t te;
    span_t sp={0};
    imud_t *pimu=(imud_t*)malloc(sizeof(imud_t)*ws);

    trace(3,"detstatic:\n");

    te=ts<=0.0?imu->data[imu->n-1].time:timeadd(imu->data[0].time,ts);

    /* first epoch */
    if (sp.ts.time==0.0) {
        sp.ts=imu->data[0].time;
    }
    while (1) {

        if (timediff(imu->data[i].time,te)>0.0||i>imu->n) break;
        if (i++<ws) continue;

        for (k=0,j=i-ws/2;j<i+ws/2;j++) pimu[k++]=imu->data[j];

        /* static detector fot imu measurement */
        if (detstaticw(ins,pimu,k,opt)) flag=1; else flag=0;

        switch (flag) {
            case 0: addspan(span,&sp);
                    sp.ts=imu->data[i].time; sp.n=0; sp.te=sp.ts; break;
            case 1: sp.te=imu->data[i].time; sp.n++;              break;
        }
    }
    for (k=0,i=0;i<span->n;i++) {
        if (timediff(span->tt[i].te,
                     span->tt[i].ts)>opt->zvopt.mt) {
            k++; continue;
        }
        span->tt[i].n=0;
    }
    free(pimu);
    return k;
}
/* runs the generalized likelihood test for detect static imu measurement----
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 *           double *pos     I  ins position (lat,lon,h)
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * -------------------------------------------------------------------------*/
extern int detstatic_GLRT(const imud_t *imu,int n,const insopt_t *opt,
                          const double *pos)
{
    int i,j;
    double gn[3],ym[3]={0},tmp[3];
    double T=0,sg=opt->zvopt.sig_g,sa=opt->zvopt.sig_a;

    trace(3,"detstatic_GLRT:\n");

    gravity_ned(pos,gn);

    for (i=0;i<3;i++) {
        for (j=0;j<n;j++) ym[i]+=imu[j].accl[i]; ym[i]/=n;
    }
    /* detect zero velocity */
    for (i=0;i<n;i++) {
        for (j=0;j<3;j++) tmp[j]=imu[i].accl[j]-norm(gn,3)/norm(ym,3)*ym[j];
        T+=SQR(norm(imu[i].gyro,3))/SQR(sg)+SQR(norm(tmp,3))/SQR(sa);
    }
    T/=n;

    trace(3,"T=%6.4lf,gamma=%6.4lf\n",T,opt->zvopt.gamma[0]);
    return T<opt->zvopt.gamma[0];
}
/* runs the acceleration moving variance detector ---------------------------
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * -------------------------------------------------------------------------*/
extern int detstatic_MV(const imud_t *imu,int n,const insopt_t *opt)
{
    int i,j;
    double ym[3],tmp[3],T=0.0;

    trace(3,"detstatic_MV:\n");

    for (i=0;i<3;i++) {
        for (j=0;j<n;j++) ym[i]+=imu[j].accl[i]; ym[i]/=n;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<3;j++) tmp[j]=imu[i].accl[j]-ym[j];
        T+=SQR(norm(tmp,3));
    }
    T/=(SQR(opt->zvopt.sig_a)*n);

    trace(3,"T=%6.4lf,gamma=%6.4lf\n",T,opt->zvopt.gamma[1]);
    return T<opt->zvopt.gamma[1];
}
/* runs the acceleration magnitude detector ---------------------------------
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 *           double *pos     I  ins position (lat,lon,h)
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * -------------------------------------------------------------------------*/
extern int detstatic_MAG(const imud_t *imu,int n,const insopt_t *opt,
                         const double *pos)
{
    int i,j;
    double gn[3],sa2=SQR(opt->zvopt.sig_a),T;

    trace(3,"detstatic_MAG:\n");

    gravity_ned(pos,gn);

    for (T=0.0,i=0;i<n;i++) {
        T+=SQR(norm(gn,3)-norm(imu[i].accl,3));
    }
    T/=(sa2*n);

    trace(3,"T=%6.4lf,gamma=%6.4lf\n",T,opt->zvopt.gamma[2]);
    return T<opt->zvopt.gamma[2];
}
/* runs the angular rate energy detector -------------------------------------
 * args   :  imud_t *imu     I  imu measurement data
 *           int n           I  number of imu measurement
 *           insopt_t *opt   I  ins options
 * return : 1: zero velocity,0: moving
 * note : n means windows size for static detector
 * --------------------------------------------------------------------------*/
extern int detstatic_ARE(const imud_t *imu,int n,const insopt_t *opt)
{
    int i; double T=0.0,sg2=SQR(opt->zvopt.sig_g);

    trace(3,"detstatic_ARE:\n");

    for (i=0;i<n;i++) T+=SQR(norm(imu[i].gyro,3)); T/=(sg2*n);

    trace(3,"T=%6.4lf,gamma=%6.4lf\n",T,opt->zvopt.gamma[3]);
    return T<opt->zvopt.gamma[3];
}
/* use odometry measurement to detect static zero velocity-------------------
 * args   :  insopt_t *opt  I  ins zero velocity detect options
 *           odod_t *odo    I  odometry measurement data
 * return : 1: zero velocity, 0: motion
 * --------------------------------------------------------------------------*/
extern int detstatic_ODO(const insopt_t *opt,const odod_t *odo)
{
    int info=0;
    static double dt=0.0,dr=0.0;

    trace(3,"detstatic_odo:\n");

    dt+=odo->dt; dr+=odo->dr;

    if (dt>(opt->zvopt.odt==0.0?
            0.1:opt->zvopt.odt)) {
        if (dr/dt==0.0) {
            trace(3,"detect zero velocity by using odometry\n");
            info=1;
        }
        dr=dt=0.0;
    }
    return info;
}








