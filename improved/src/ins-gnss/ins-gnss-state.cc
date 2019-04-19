/*-----------------------------------------------------------------------------
* ins-gnss-state.cc : ins-gnss coupled estimated states interface
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
* history : 2017/02/06 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* get number of position states---------------------------------------------*/
extern int xnP(const insopt_t *opt) {return 3;}
/* get number of velocity states---------------------------------------------*/
extern int xnV(const insopt_t *opt) {return 3;}
/* get number of attitude states---------------------------------------------*/
extern int xnA(const insopt_t *opt) {return 3;}
/* get number of accl/gyro bias states---------------------------------------*/
extern int xnBa(const insopt_t *opt) {return ((opt)->baopt==INS_BAEST?3:0);}
extern int xnBg(const insopt_t *opt) {return ((opt)->bgopt==INS_BGEST?3:0);}
/* get number of time synchronization error states---------------------------*/
extern int xnDt(const insopt_t *opt) {return ((opt)->estdt==INS_DTEST?1:0);}

/* get number of residual scale factors of gyroscopes------------------------*/
extern int xnSg(const insopt_t *opt) {return ((opt)->estsg==INS_SGEST?3:0);}
/* get number of residual scale factors of accl------------------------------*/
extern int xnSa(const insopt_t *opt) {return ((opt)->estsa==INS_SAEST?3:0);}

/* get number of non-orthogonal between sensor axes for gyro/accl------------*/
extern int xnRg(const insopt_t *opt) {return ((opt)->estrg==INS_RGEST?6:0);}
extern int xnRa(const insopt_t *opt) {return ((opt)->estra==INS_RAEST?6:0);}

/* get number of lever arm for body to ant.----------------------------------*/
extern int xnLa(const insopt_t *opt) {return ((opt)->estlever==INS_LAEST?3:0);}

/* get number of odometry scale factor---------------------------------------*/
extern int xnOs(const insopt_t *opt) {return ((opt)->estodos?1:0);}
/* get number of odometry lever arm from imu body frame----------------------*/
extern int xnOl(const insopt_t *opt) {return ((opt)->estodol?3:0);}
/* get number of odometry misalignment of body frame and rear frame----------*/
extern int xnOa(const insopt_t *opt) {return ((opt)->estodoa?3:0);}
/* get number of misalignment estimate states of camera to imu body----------*/
extern int xnCm(const insopt_t *opt) {return ((opt)->estcama?3:0);}
/* get number of lever arm states of camera to imu body----------------------*/
extern int xnCla(const insopt_t *opt) {return ((opt)->estcaml?3:0);}

/* get number of camera calibration parameters (fx,fy,cx,cy)-----------------*/
extern int xnCfo(const insopt_t *opt) {return ((opt)->estcam_fo?4:0);}

/* get number of camera calibration parameters (k1,k2,p1,p2)-----------------*/
extern int xnCkp(const insopt_t *opt) {return ((opt)->estcam_kp?4:0);}

/* get number of misalignemnt estimate states of b-frame to v-frame----------*/
extern int xnVm(const insopt_t *opt)
{
    return (opt)->gopt?(((prcopt_t*)(opt)->gopt)->mode==PMODE_INS_LGNSS||
                        ((prcopt_t*)(opt)->gopt)->mode==PMODE_INS_TGNSS)
                       &&opt->estmisv?3:0:0;
}
/* get number of close-loop correction states--------------------------------*/
extern int xnCl(const insopt_t *opt)
{
    return xnP (opt)+xnV (opt)+xnA (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt);
}
/* get number of observation data frequency----------------------------------*/
extern int xnIF(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->ionoopt==IONOOPT_IFLC?1:
                        ((prcopt_t*)opt->gopt)->nf:0);
}
/* get number of receiver clock bias (non-close-loop correction states)------*/
extern int xnRc(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->mode==PMODE_INS_TGNSS
           &&(opt)->tc<=INSTC_RTK?4:0:0);
}
/* get number of receiver clock drift (non-close-loop correction states) ----*/
extern int xnRr(const insopt_t *opt)
{
    return ((opt)->dopp?1:0);
}
/* get number of iono,trop,phase bias and h/w bias states are non-close-loop
 * correction states
 * --------------------------------------------------------------------------*/
extern int xnI(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->ionoopt!=IONOOPT_EST?0:MAXSAT:0);
}
extern int xnT(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->tropopt<TROPOPT_EST?0:
                       (((prcopt_t*)(opt)->gopt)->tropopt<TROPOPT_ESTG?2:6):0);
}
extern int xnL(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->glomodear!=2?0:NFREQGLO:0);
}
extern int xnD(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->nf>=3?1:0:0);
}
extern int xnB(const insopt_t *opt)
{
    return ((opt)->gopt?((prcopt_t*)(opt)->gopt)->mode==PMODE_INS_TGNSS
           &&((opt)->tc==INSTC_PPK||(opt)->tc==INSTC_RTK)?(MAXSAT*xnIF(opt)):0:0);
}
extern int xnRx(const insopt_t *opt)
{
    return (xnP (opt)+xnV (opt)+xnA (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
            xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
            xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt)+xnRr(opt)+
            xnI (opt)+xnT(opt)+xnL(opt)+xnD(opt));
}
/* get number of all states--------------------------------------------------*/
extern int xnX(const insopt_t *opt) {return xnRx(opt)+xnB(opt);}

extern int xiA (const insopt_t *opt) {return 0;}    /* attitude states index */
extern int xiV (const insopt_t *opt) {return 3;}    /* velocity states index */
extern int xiP (const insopt_t *opt) {return 6;}    /* position states index */
extern int xiBa(const insopt_t *opt) {return 9;}    /* accl bias state index */

/* get index of gyro bias states---------------------------------------------*/
extern int xiBg(const insopt_t *opt) {return xnA(opt)+xnV(opt)+xnP(opt)+xnBa(opt);}

/* get index of time synchronization error for IMU and GNSS user equipment---*/
extern int xiDt(const insopt_t *opt)
{
    return xnA(opt)+xnV(opt)+xnP(opt)+xnBa(opt)+xnBg(opt);
}
/* get index of residual scale factors of gyroscopes-------------------------*/
extern int xiSg(const insopt_t *opt)
{
    return xnA(opt)+xnV(opt)+xnP(opt)+xnBa(opt)+xnBg(opt)+xnDt(opt);
}
/* get index of residual scale factors of accl state-------------------------*/
extern int xiSa(const insopt_t *opt)
{
    return xnA(opt)+xnV(opt)+xnP(opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt);
}
/* get index of non-orthogonal between sensor axes for gyro states-----------*/
extern int xiRg(const insopt_t *opt)
{
    return xnA (opt)+xnV(opt)+xnP(opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt);
}
/* get index of non-orthogonal between sensor axes for accl states-----------*/
extern int xiRa(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP(opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt);
}
/* get index lever arm states------------------------------------------------*/
extern int xiLa(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt);
}
/* get index of odometry scale factor states---------------------------------*/
extern int xiOs(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt);
}
/* get index of odometry lever arm states------------------------------------*/
extern int xiOl(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt);
}
/* get index of odometry misalignment of body frame and rear frame states----*/
extern int xiOa(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt);
}
/* get index of misalignment of camera to imu body---------------------------*/
extern int xiCm(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt);
}
/* get index of lever arm of camera to imu ----------------------------------*/
extern int xiCl(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt);
}
/* get index of camera calibration (fx,fy,cx,cy)------------------------------*/
extern int xiCfo(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt);
}
/* get index of camera calibration (k1,k2,p1,p3)-----------------------------*/
extern int xiCkp(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt);
}
/* get index of misalignment of b-frame to v-frame---------------------------*/
extern int xiVm(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt);
}
/* get index of receiver clock bias state------------------------------------*/
extern int xiRc(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt);
}
/* get index of receiver clock drift-----------------------------------------*/
extern int xiRr(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt);
}
/* get index of ionos (s:satellite no)---------------------------------------*/
extern int xiIo(const insopt_t *opt,int s)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt)+xnRr(opt)+(s-1);
}
/* get index of tropos (r:0=rov,1:ref)---------------------------------------*/
extern int xiTr(const insopt_t *opt,int r)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt)+xnRr(opt)+
           xnI (opt)+xnT(opt)/2*r;
}
/* get index of receiver h/w bias--------------------------------------------*/
extern int xiLl(const insopt_t *opt,int f)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt)+xnRr(opt)+
           xnI(opt)+xnT(opt)+f;
}
/* get index of L5-receiver-dcb ---------------------------------------------*/
extern int xiDl(const insopt_t *opt)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt)+xnRr(opt)+
           xnI(opt)+xnT(opt)+xnL(opt);
}
/* get index of phase bias (s:satno,f:freq)----------------------------------*/
extern int xiBs(const insopt_t *opt,int s,int f)
{
    return xnA (opt)+xnV (opt)+xnP (opt)+xnBa(opt)+xnBg(opt)+xnDt(opt)+xnSg(opt)+
           xnSa(opt)+xnRg(opt)+xnRa(opt)+xnLa(opt)+xnOs(opt)+xnOl(opt)+xnOa(opt)+
           xnCm(opt)+xnCla(opt)+xnCfo(opt)+xnCkp(opt)+xnVm(opt)+xnRc(opt)+xnRr(opt)+
           xnI(opt)+xnT(opt)+xnL(opt)+xnD(opt)+MAXSAT*f+s-1;
}
/* disable ins states when filter--------------------------------------------
 * args   :  insopt_t *opt    I   ins update options
 *           int index        I   index of unused ins states
 *           insstate_t *ins  IO  ins states
 *           double *x        IO  ins estimated states
 * return 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int unusex(const insopt_t *opt,int index,insstate_t *ins,double *x)
{
    if (index>=xnCl(opt)) return 0; x[index]=DISFLAG; return 1;
}