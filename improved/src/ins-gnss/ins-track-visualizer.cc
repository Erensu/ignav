/*------------------------------------------------------------------------------
 * ins-track.cc : feature points tracking functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Bruce D. Lucas and Takeo Kanade. An Iterative Image Registration Technique
 *        with an Application to Stereo Vision. International Joint Conference on
 *        Artificial Intelligence, pages 674-679, 1981.
 *    [3] Carlo Tomasi and Takeo Kanade. Detection and Tracking of Point Features.
 *        Carnegie Mellon University Technical Report CMU-CS-91-132, April 1991.
 *    [4] Jianbo Shi and Carlo Tomasi. Good Features to Track. IEEE Conference on
 *        Computer Vision and Pattern Recognition, pages 593-600, 1994.
 *    [5] Stan Birchfield. Derivation of Kanade-Lucas-Tomasi Tracking Equation.
 *        Unpublished, January 1997.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/10/25 1.0 new
 *-----------------------------------------------------------------------------*/
#include <carvig.h>
#include <cv.h>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

/* image buffer convert to Image struct--------------------------------------
 * args  :  img_t *img  I  image gray level measurement data
 *          Image **Img O  image rgb format measurement data
 * return: none
 * --------------------------------------------------------------------------*/
static int img2Img(const img_t *img,cv::Mat &Img)
{
    int i,j;
    for (i=0;i<Img.rows;i++) for (j=0;j<Img.cols;j++) {
            Img.at<cv::Vec3b>(i,j)[0]=img->data[i*img->w+j];
            Img.at<cv::Vec3b>(i,j)[1]=img->data[i*img->w+j];
            Img.at<cv::Vec3b>(i,j)[2]=img->data[i*img->w+j];
        }
    return 1;
}
/* draw tracks-----------------------------------------------------------------
 * args:    track_t *track  I  tracking data
 *          voopt_t *opt    I  options
 * return: none
 * ---------------------------------------------------------------------------*/
extern void drawtrack(const track_t *track,const voopt_t *opt)
{
    char allName[126][126];
    int i,j;
    img_t *imgp;

    trace(3,"drawtrack: n=%d\n",track->n);

    for (i=0;i<track->n;i++) { /* do for all tracks */
                                                                                                          
        /* for each track. */
        for (j=0;j<track->data[i].n;j++) {

            if (!(imgp=getimgdata(track->data[i].data[j].time))) continue;

            cv::Mat outImage(imgp->h,imgp->w,CV_8UC3);
            img2Img(imgp,outImage);

            cv::Point point;
            point.x=(int)track->data[i].data[j].u;
            point.y=(int)track->data[i].data[j].v;

            cv::circle(outImage,point,2,cv::Scalar(0,0,255),1);
            cv::circle(outImage,point,8,cv::Scalar(0,0,255),1);

            sprintf(allName[j],"image_%d_%d",track->data[i].uid,track->data[i].first_frame+j);
            cv::namedWindow(allName[j],CV_WINDOW_AUTOSIZE);
            cv::imshow(allName[j],outImage);                                                                                                                                  
            cv::waitKey(0);                                                                                                                                                     
        }
        for (j=0;j<track->data[i].n;j++) {
            cv::destroyWindow(allName[j]);
        }
    }
}
/* draw tracks-----------------------------------------------------------------
 * args:    trackd_t *track  I  tracking data
 *          voopt_t *opt     I  options
 * return: none
 * ---------------------------------------------------------------------------*/
extern void drawtrackd(const trackd_t *track,const voopt_t *opt)
{
    static int first=1;
    int j;
    img_t *imgp;

    trace(3,"drawtrackd: n=%d\n",track->n);

    if (first) {
        cv::namedWindow("Match-Image",CV_WINDOW_AUTOSIZE);
        first=0;
    }
    /* for each track. */
    for (j=0;j<track->n;j++) {

        if (!(imgp=getimgdata(track->data[j].time))) continue;

        cv::Mat outImage(imgp->h,imgp->w,CV_8UC3);
        img2Img(imgp,outImage);

        cv::Point point;
        point.x=(int)track->data[j].u;
        point.y=(int)track->data[j].v;

        cv::circle(outImage,point,2,cv::Scalar(0,0,255),1);
        cv::circle(outImage,point,8,cv::Scalar(0,0,255),1);
        cv::imshow("Match-Image",outImage);
        cv::waitKey(0);
    }
}
/* display image-------------------------------------------------------------*/
extern void dipsplyimg(const img_t *img)
{
    static int first=1;
    static char text[126];

    if (img==NULL) return;
    if (first) {
        cv::namedWindow("Display-Image",CV_WINDOW_AUTOSIZE);
        first=0;
    }
    cv::Mat outImage(img->h,img->w,CV_8UC3);
    img2Img(img,outImage);
    sprintf(text,"id=%d",img->id);

    cv::putText(outImage,text,cv::Point(8,15),cv::FONT_HERSHEY_PLAIN,1.0,cv::Scalar(0,255,255));
    cv::imshow("Display-Image",outImage);
    cv::waitKey(3);
}
/*----------------------------------------------------------------------------
 *　Combines two images by scacking one on right of the other
 *
 *　@param img1 left image
 *　@param img2 right image
 *
 *　@return Returns the image resulting from stacking \a img1 on top if \a img2
 *----------------------------------------------------------------------------*/
static IplImage* stack_imgs(IplImage* img1,IplImage* img2)
{
    IplImage* stacked=cvCreateImage(cvSize(MAX(img1->height,img2->height),img1->width+img2->width),
                                    IPL_DEPTH_8U,3);
    cvZero(stacked);
    cvSetImageROI(stacked,cvRect(0,0,img1->width,img1->height));
    cvAdd(img1,stacked,stacked,NULL);
    cvSetImageROI(stacked,cvRect(img2->width,0,img2->width,img2->height));
    cvAdd(img2,stacked,stacked,NULL);
    cvResetImageROI(stacked);
    return stacked;
}
/* draw all track data for debugs--------------------------------------------*/
extern void drawalltrack(const track_t *track)
{
    static int first=1;
    img_t *imgp1,*imgp2;
    gtime_t t1={0},t2={0};
    int i;
    char text[256];

    if (track->nnew) {
        t1=track->data[track->newtrack[0]].ts;
        t2=track->data[track->newtrack[0]].te;
    }
    else if (track->nupd) {
        t1=track->data[track->updtrack[0]].data[0].time;
        t2=track->data[track->updtrack[0]].data[track->data[track->updtrack[0]].n-1].time;
    }
    if (!(imgp1=getimgdata(t1))||!(imgp2=getimgdata(t2))) {
        return;
    }
    cv::Mat I1(imgp1->h,imgp1->w,CV_8UC3);
    img2Img(imgp1,I1);
    sprintf(text,"id=%d",imgp1->id);
    cv::putText(I1,text,cv::Point(8,15),cv::FONT_HERSHEY_PLAIN,1.0,cv::Scalar(0,255,255));

    cv::Mat I2(imgp2->h,imgp2->w,CV_8UC3);
    img2Img(imgp2,I2);
    sprintf(text,"id=%d",imgp2->id);
    cv::putText(I2,text,cv::Point(8,15),cv::FONT_HERSHEY_PLAIN,1.0,cv::Scalar(0,255,255));

    cv::Mat I(MAX(imgp1->h,imgp2->h),imgp1->w+imgp2->w,CV_8UC3);
    cv::Mat ROI1=I(Rect(0,0,imgp1->w,MAX(imgp1->h,imgp2->h)));
    I1.copyTo(ROI1);
    cv::Mat ROI2=I(Rect(imgp1->w,0,imgp2->w,MAX(imgp1->h,imgp2->h)));
    I2.copyTo(ROI2);

    cv::Mat II(MAX(imgp1->h,imgp2->h),imgp1->w+imgp2->w,CV_8UC3);
    I.copyTo(II);

    if (first) {
        cv::namedWindow("All-Match-Image-Display",CV_WINDOW_AUTOSIZE);
        first=0;
    }
    cv::Mat overlay;
    for (i=0;i<track->nnew;i++) {

        cv::Point p1;
        p1.x=(int)track->data[track->newtrack[i]].data[0].u;
        p1.y=(int)track->data[track->newtrack[i]].data[0].v;
        cv::circle(I,p1,2,cv::Scalar(0,0,255),1);

        cv::Point p2;
        p2.x=(int)track->data[track->newtrack[i]].data[1].u+imgp1->w;
        p2.y=(int)track->data[track->newtrack[i]].data[1].v;
        cv::circle(I,p2,2,cv::Scalar(0,255,0),1);
#if 0
        I.copyTo(overlay);

        cv::line(overlay,p1,p2,CV_RGB(255,0,255),1,8,0);
        cv::addWeighted(overlay,0.2,I,0.8,0,I);
#endif
    }
    for (i=0;i<track->nupd;i++) {
        cv::Point p1;
        p1.x=(int)track->data[track->updtrack[i]].data[track->data[track->updtrack[i]].n-2].u;
        p1.y=(int)track->data[track->updtrack[i]].data[track->data[track->updtrack[i]].n-2].v;
        cv::circle(I,p1,2,cv::Scalar(0,0,255),1);

        cv::Point p2;
        p2.x=(int)track->data[track->updtrack[i]].data[track->data[track->updtrack[i]].n-1].u+imgp1->w;
        p2.y=(int)track->data[track->updtrack[i]].data[track->data[track->updtrack[i]].n-1].v;
        cv::circle(I,p2,2,cv::Scalar(0,255,0),1);
#if 0
        I.copyTo(overlay);

        cv::line(overlay,p1,p2,CV_RGB(255,0,0),1,8,0);
        cv::addWeighted(overlay,0.2,I,0.8,0,I);
#endif
    }
    cv::imshow("All-Match-Image-Display",I);
    cv::waitKey(3);
}


