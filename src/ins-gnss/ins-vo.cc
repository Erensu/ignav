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
#include <navlib.h>

#if ENAOPENCV
#include <cxcore.h>
#include <highgui.h>
#endif

/* rotation parameters and translation parameters convert to transformation---
 * arg   : double *R  I  rotation matrix
 *         double *t  I  translation matrix
 *         double *T  O  transformation matrix
 * return: none
 * --------------------------------------------------------------------------*/
extern void rt2tf(const double *R,const double *t,double *T)
{
    setzero(T,4,4);
    T[0]=R[0]; T[4]=R[3]; T[ 8]=R[6]; T[12]=t[0];
    T[1]=R[1]; T[5]=R[4]; T[ 9]=R[7]; T[13]=t[1];
    T[2]=R[2]; T[6]=R[5]; T[10]=R[8]; T[14]=t[2]; T[15]=1.0;
}
/* transform matrix convert to rotation and translation parameters-----------*/
extern void tf2rt(const double *T,double *R,double *t)
{
    seteye(R,3); setzero(t,1,3);
    R[0]=T[0 ]; R[3]=T[4 ]; R[6]=T[8 ];
    R[1]=T[1 ]; R[4]=T[5 ]; R[7]=T[9 ];
    R[2]=T[2 ]; R[5]=T[6 ]; R[8]=T[10];
    t[0]=T[12]; t[1]=T[13]; t[2]=T[14];
}
/* reset monocular camera motion estimator----------------------------------*/
extern void resetmonoa()
{

}
/* free monocular camera motion estimator------------------------------------*/
extern void freemonoa()
{

}
/* reset stereo camera motion estimator--------------------------------------*/
extern void resetstereo()
{

}
/* free stereo camera motion estimator---------------------------------------*/
extern void freestereo()
{
    
}
/* estimate monocular camera motion for input successive image---------------
 * args  :  voaid_t *opt  I  visual odometry options
 *          img_t *img    I  first image raw data
 *          double *Tr    O  transformation matrix
 *          double *vc    O  camera velocity in vision-frame
 *          double *var   O  motion estimate variance
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int estmotionmonoa(const voaid_t *opt,img_t *img,double *Tr,double *vc,
                          double *var)
{
    return 1;
}
/* estimate monocular camera motion for input two successive image -----------
 * args  :  voaid_t *opt  I  visual odometry options
 *          img_t *img1   I  first image raw data
 *          img_t *img2   I  second image raw data
 *          double *Tr    O  transformation matrix
 *          double *var   O  estimate variance
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int estmotionmono(const voaid_t *opt,img_t *img1,img_t *img2,double *Tr,
                         double *var)
{
    return 1;
}
/* estimate stereo camera motion for input successive image -----------------
 * args  :  voaid_t *opt  I  visual odometry options
 *          img_t *img    I  first image raw data
 *          double *Tr    O  transformation matrix
 *          double *vc    O  camera velocity in vision-frame
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int estmotionstereo(const voaid_t *opt,img_t *img,double *Tr,double *vc,
                           double *var)
{
    return 1;
}
#if ENAOPENCV
/* image buffer convert to Image struct--------------------------------------
 * args  :  img_t *img  I  image gray level measurement data
 *          Image **Img O  image rgb format measurement data
 * return: none
 * --------------------------------------------------------------------------*/
static int img2Img(const img_t *img,IplImage **Img)
{
    int step,i,j; uchar* pimg=NULL;
    CvSize s0={0};

    trace(3,"img2Img:\n");

    if (img==NULL||Img==NULL) {
        return 0;
    }
    s0.width=img->w; s0.height=img->h;
    *Img=cvCreateImage(s0,8,3);
    step=(*Img)->widthStep/sizeof(uchar);

    pimg=(uchar *)(*Img)->imageData;

    for (i=0;i<(*Img)->height;i++) {

        for (j=0;j<(*Img)->width;j++) {
            pimg[i*step+3*j+0]=img->data[i*img->w+j];
            pimg[i*step+3*j+1]=img->data[i*img->w+j];
            pimg[i*step+3*j+2]=img->data[i*img->w+j];
        }
    }
    return 1;
}
/* display image-------------------------------------------------------------*/
extern void dipsplyimg(const img_t *img)
{
    if (img==NULL) return;
    IplImage *Img=NULL;

    if (!img2Img(img,&Img)) {
        trace(2,"image data invalid\n");
        return;
    }
    IplImage *pImage=Img;
    cvShowImage("Image Display",pImage);
    cvWaitKey(1);
}
#endif
/* visual odometry estimator-------------------------------------------------
 * args  :  voaid_t *opt  I  visual odometry options
 *          img_t *img    I  first image raw data
 *          double *Tr    O  transformation matrix
 *          double *vc    O  camera velocity in vision-frame
 *          double *var   O  camera pose estimate variance
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int estvo(const voaid_t *opt,img_t *img,double *Tr,double *vc,double *var)
{
    trace(3,"estvo:\n");
    return 0;
}



