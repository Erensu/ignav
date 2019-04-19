/*------------------------------------------------------------------------------
* ins-vo.cc : ins-visual odometry loosely coupled common functions
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
* history : 2018/03/15 1.0 new
*-----------------------------------------------------------------------------*/
#include <vo-mono.h>
#include <carvig.h>

#if ENAOPENCV
#include <cxcore.h>
#include <highgui.h>
#endif

/* constants ----------------------------------------------------------------*/
#define MAXTIMEDIFF   10.0                /* max time difference for estimate camera motion */

static int first=1;                       /* flag of first time to start visual odometry */
static VisualOdometryMono *Mono=NULL;     /* monocular camera visual odometry estimator */

/* reset monocular camera motion estimator----------------------------------*/
extern void resetmonoa()
{
    first=1; if (Mono) delete Mono;
}
/* free monocular camera motion estimator------------------------------------*/
extern void freemonoa()
{
    first=1; if (Mono) delete Mono;
}
/* estimate monocular camera motion for input successive image---------------
 * args  :  voopt_t *opt  I  visual odometry options
 *          img_t *img    I  first image raw data
 *          double *Tr    O  transformation matrix
 *          double *ratio O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
static int estmotionmonoa(const voopt_t *opt,const img_t *img,double *Tr,double *ratio)
{
    Matrix tr;
    static gtime_t pt={0};
    int i,j,dims[3]={img->w,img->h,img->w};
    double tt=0.0;

    trace(3,"estmotionmono: time=%s\n",time_str(img->time,4));

    if (!pt.time) pt=img->time;
    else {
        tt=timediff(img->time,pt); pt=img->time;
    }
    seteye(Tr,4);
    if (img->data==NULL) {
        trace(2,"no valid image raw data\n");
        return 0;
    }
    /* start visual odometry program */
    if (first) {
        if (Mono) delete Mono;

        /* set most important visual odometry parameters */
        VisualOdometryMono::parameters para;

        /* visual odometry parameters */
        para.calib.f =opt->calib.f;
        para.calib.cu=opt->calib.cu;
        para.calib.cv=opt->calib.cv;
        para.inlier_threshold=opt->inlier_thres;
        para.motion_threshold=opt->motion_thres;
        para.ransac_iters    =opt->ransac_iters;
        para.height          =opt->height;

        /* init visual odometry */
        Mono=new VisualOdometryMono(para);
        first=0;
        goto estmono;
    }
    if (tt>MAXTIMEDIFF||(pt.time&&fabs(tt)<1E-5)) {
        trace(2,"too large time difference of two image\n");
        return 0;
    }
estmono:
    /* estimate camera motion */
    if (!Mono->process(img->data,dims,false)) {
        trace(2,"estimate monocular camera motion fail\n");
        return 0;
    }
    tr=Mono->getMotion();
    if (ratio) {
        *ratio=(double)Mono->getNumberOfInliers()/(double)Mono->getNumberOfMatches();
    }
    if (tr.m==0||tr.n==0) {
        trace(2,"estimate monocular camera motion fail\n");
        return 0;
    }
    /* output estimate result */
    if (Tr) {
        for (i=0;i<4;i++) for (j=0;j<4;j++) {
            Tr[i+j*4]=tr.val[i][j];
        }
        trace(3,"Tr=\n");
        tracemat(3,Tr,4,4,12,8);
    }
    return 1;
}
/* visual odometry estimator-------------------------------------------------
 * args  :  voopt_t *opt  I  visual odometry options
 *          img_t *img    I  first image raw data
 *          double *dTr   O  transformation matrix
 *          double *ratio O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int estvo(const voopt_t *opt,const img_t *img,double *dTr,double *ratio)
{
    return estmotionmonoa(opt,img,dTr,ratio);
}




