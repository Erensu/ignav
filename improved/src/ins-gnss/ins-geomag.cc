/*-----------------------------------------------------------------------------
* ins-geomag.cc : ins-geomag measurement coupled common functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Daehee Won,Performance Improvement of Inertial Navigation System by
*        Using Magnetometer with Vehicle Dynamic Constraints,2015
*    [3] Shin,Accuracy Improvement of Low Cost INS/GNSS for Land Application.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/09/06 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"
#include "geomag.h"

/* constants ----------------------------------------------------------------*/
#define VAR_MAG    (5.0*D2R)       /* variance of magnetic heading */

/* global variance define----------------------------------------------------*/
extern BFieldModel mag_model={{0}};
extern BField      mag={0};

/* reads model coefficients from a file--------------------------------------*/
extern int magmodel(const char *file)
{
    trace(3,"magmodel,file=%s\n:",file);
    return read_model(&mag_model,file);
}
/* tilt compensation to obtain the horizontal magnetic measurements----------*/
static void untilt(const mag_t *data,const double *Cbn,double *mh)
{
    double r=-atan (Cbn[2]/SQRT(1.0-SQR(Cbn[2])));
    double p= atan2(Cbn[5],Cbn[8]);
    double Cbh[6];

    Cbh[0]=cos(p); Cbh[2]=sin(r)*sin(p); Cbh[4]=cos(r)*sin(p);
    Cbh[1]=0.0;    Cbh[3]=cos(r);        Cbh[5]=-sin(r);
    matmul("NN",2,1,3,1.0,Cbh,data->val,0.0,mh);
}
/* compensate for the declination angle--------------------------------------*/
static double undecli(gtime_t time,const double *pos)
{
    static double ep[6]={0};

    time2epoch(time,ep);
    get_field_components(&mag,&mag_model,pos[2],Units::kUnitsMeters,
                         CoordinateSystem::kCoordSysGeodetic,
                         pos[0],pos[1],
                         julday(ep[1],ep[2],ep[0]));
    return mag.d;
}
/* compensate for magnetic data distortions----------------------------------*/
static void undistortmag(const magopt_t *opt,double *mh)
{
    mh[0]=opt->sx*mh[0]+opt->ox;
    mh[1]=opt->sy*mh[1]+opt->oy;
}
/* magnetometer measurement data calibration---------------------------------
 * args:    magopt_t *opt  I  magnetometer options
 *          mag_t *data    I  magnetometer measurement data
 *          double *Cbn    I  magnetometer attitude relative to n-frame
 *          double *mh     I  magnetic in n-frame
 * return : magnetic heading angle (deg)
 * --------------------------------------------------------------------------*/
extern double maghead(const magopt_t *opt,const mag_t *data,
                      const double *Cbn,
                      const double *pos)
{
    double mh_[2]={0},decli,yaw;

    trace(3,"magcalib:\n");

    untilt(data,Cbn,mh_);   /* tilt compensation */
    undistortmag(opt,mh_);  /* magnetometer undistort */

    /* declination of the field */
    decli=undecli(data->time,pos);

    /* magnetic heading angle */
    yaw=atan2(mh_[1],mh_[0]);
    return yaw-decli;
}
/* jacobian of heading measurement wrt. attitude error-----------------------*/
static void jacob_head2att(const double *Cbe,const double *pos,double *dhdatt)
{
    double Cne[9],Cbn[9];

    ned2xyz(pos,Cne);
    matmul("TN",3,3,3,1.0,Cne,Cbe,0.0,Cbn);

    dhdatt[0]=(Cbn[0]*(Cne[5]*Cbe[1]-Cne[4]*Cbe[2])-
               Cbn[1]*(Cne[2]*Cbe[1]-Cne[1]*Cbe[2]))/
              (SQR(Cbn[0])+SQR(Cbn[1]));

    dhdatt[1]=(Cbn[0]*(-Cne[5]*Cbe[0]+Cne[3]*Cbe[2])-
               Cbn[1]*(-Cne[2]*Cbe[1]+Cne[0]*Cbe[2]))/
              (SQR(Cbn[0])+SQR(Cbn[1]));

    dhdatt[0]=(Cbn[0]*(Cne[4]*Cbe[0]-Cne[3]*Cbe[1])-
               Cbn[1]*(Cne[1]*Cbe[0]-Cne[0]*Cbe[1]))/
              (SQR(Cbn[0])+SQR(Cbn[1]));
}
/* ekf filter for magnetometer measurement data------------------------------*/
static int magfilt(insstate_t *ins,const insopt_t *opt,const mag_t *data)
{
    double *H,*x,*P,R,v,Cbn[9],pos[3],Cne[9];
    double magh,yaw;

    trace(3,"magfilt: nx=%d\n",ins->nx);

    H=zeros(1,ins->nx); x=ins->x; P=ins->P;

    ecef2pos(ins->re,pos);
    ned2xyz(pos,Cne);
    matmul("TN",3,3,3,1.0,Cne,ins->Cbe,0.0,Cbn);

    /* magnetic heading */
    magh=maghead(&opt->magopt,data,Cbn,pos);
    yaw =atan2(Cbn[1],Cbn[0])*D2R;

    /* jaobian of magnetic heading */
    jacob_head2att(ins->Cbe,pos,H+xiA(opt));

    v=(NORMANG(magh)-NORMANG(yaw))*D2R;
    R=SQR(VAR_MAG);

    if (filter(x,P,H,&v,&R,ins->nx,1)) {  /* update */
        trace(2,"filter error\n");
        free(H);
        return 0;
    }
    /* close loop cor. for ins states */
    clp(ins,opt,x);
    free(H); return 1;
}
/* update ins states using magnetometer measurement data---------------------
 * args:    insstate_t *ins  IO  ins states
 *          insopt_t *opt    I   ins options
 *          mag_t    *data   I   magnetometer measurement data
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int magnetometer(insstate_t *ins,const insopt_t *opt,
                        const mag_t *data)
{
    trace(3,"maghead:\n");

    if (norm(data->val,3)<=0.0) {
        trace(2,"no valid magnetometer data\n");
        return 0;
    }
    if (!magfilt(ins,opt,data)) {
        trace(2,"magnetic head update fail\n");
        return 0;
    }
    ins->stat=INSS_MAGH;
    return 1;
}
