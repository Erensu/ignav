/*------------------------------------------------------------------------------
 * ins-vo-opencv-interface.cc : openCV interface common functions
 *
 * reference :
 *    [1] https://www.opencv.org/
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"
#include "relative-pose.h"

/* recover camera pose from essential matrix-----------------------------------*/
namespace cv {
    void decomposeEssentialMat(InputArray _E,OutputArray _R1,OutputArray _R2,
                               OutputArray _t)
    {

        Mat E=_E.getMat().reshape(1,3);
        CV_Assert(E.cols==3&&E.rows==3);

        Mat D,U,Vt;
        SVD::compute(E,D,U,Vt);

        if (determinant(U) <0) U *=-1.;
        if (determinant(Vt)<0) Vt*=-1.;

        Mat W=(Mat_<double>(3,3)<<0,1,0,-1,0,0,0,0,1);
        W.convertTo(W,E.type());

        Mat R1,R2,t;
        R1=U*W*Vt;
        R2=U*W.t()*Vt;
        t=U.col(2)*1.0;

        R1.copyTo(_R1);
        R2.copyTo(_R2);
        t.copyTo(_t);
    }
    int recoverPose(InputArray E,InputArray _points1,InputArray _points2,InputArray _cameraMatrix,
                    OutputArray _R,
                    OutputArray _t,
                    InputOutputArray _mask)
    {

        Mat points1,points2,cameraMatrix;
        _points1.getMat().convertTo(points1,CV_64F);
        _points2.getMat().convertTo(points2,CV_64F);
        _cameraMatrix.getMat().convertTo(cameraMatrix,CV_64F);

        int npoints=points1.checkVector(2);
        CV_Assert(npoints>=0&&points2.checkVector(2)==npoints&&points1.type()==points2.type());
        CV_Assert(cameraMatrix.rows==3&&cameraMatrix.cols==3&&cameraMatrix.channels()==1);

        if (points1.channels()>1) {
            points1=points1.reshape(1,npoints);
            points2=points2.reshape(1,npoints);
        }

        double fx=cameraMatrix.at<double>(0,0);
        double fy=cameraMatrix.at<double>(1,1);
        double cx=cameraMatrix.at<double>(0,2);
        double cy=cameraMatrix.at<double>(1,2);

        points1.col(0)=(points1.col(0)-cx)/fx;
        points2.col(0)=(points2.col(0)-cx)/fx;
        points1.col(1)=(points1.col(1)-cy)/fy;
        points2.col(1)=(points2.col(1)-cy)/fy;

        points1=points1.t();
        points2=points2.t();

        Mat R1,R2,t;
        decomposeEssentialMat(E,R1,R2,t);
        Mat P0=Mat::eye(3,4,R1.type());
        Mat P1(3,4,R1.type()),P2(3,4,R1.type()),P3(3,4,R1.type()),P4(3,4,R1.type());
        P1(Range::all(),Range(0,3))=R1*1.0; P1.col(3)= t*1.0;
        P2(Range::all(),Range(0,3))=R2*1.0; P2.col(3)= t*1.0;
        P3(Range::all(),Range(0,3))=R1*1.0; P3.col(3)=-t*1.0;
        P4(Range::all(),Range(0,3))=R2*1.0; P4.col(3)=-t*1.0;

        /* Do the cheirality check.
         * Notice here a threshold dist is used to filter
         * out far away points (i.e. infinite points) since
         *there depth may vary between postive and negtive.
         * */
        double dist=50.0;
        Mat Q;
        triangulatePoints(P0,P1,points1,points2,Q);
        Mat mask1=Q.row(2).mul(Q.row(3))>0;
        Q.row(0)/=Q.row(3);
        Q.row(1)/=Q.row(3);
        Q.row(2)/=Q.row(3);
        Q.row(3)/=Q.row(3);
        mask1=(Q.row(2)<dist)&mask1;
        Q=P1*Q;
        mask1=(Q.row(2)>0)&mask1;
        mask1=(Q.row(2)<dist)&mask1;

        triangulatePoints(P0,P2,points1,points2,Q);
        Mat mask2=Q.row(2).mul(Q.row(3))>0;
        Q.row(0)/=Q.row(3);
        Q.row(1)/=Q.row(3);
        Q.row(2)/=Q.row(3);
        Q.row(3)/=Q.row(3);
        mask2=(Q.row(2)<dist)&mask2;
        Q=P2*Q;
        mask2=(Q.row(2)>0)&mask2;
        mask2=(Q.row(2)<dist)&mask2;

        triangulatePoints(P0, P3, points1, points2, Q);
        Mat mask3=Q.row(2).mul(Q.row(3))>0;
        Q.row(0)/=Q.row(3);
        Q.row(1)/=Q.row(3);
        Q.row(2)/=Q.row(3);
        Q.row(3)/=Q.row(3);
        mask3=(Q.row(2)<dist)&mask3;
        Q=P3*Q;
        mask3=(Q.row(2)>0)&mask3;
        mask3=(Q.row(2)<dist)&mask3;

        triangulatePoints(P0,P4,points1,points2,Q);
        Mat mask4=Q.row(2).mul(Q.row(3))>0;
        Q.row(0)/=Q.row(3);
        Q.row(1)/=Q.row(3);
        Q.row(2)/=Q.row(3);
        Q.row(3)/=Q.row(3);
        mask4=(Q.row(2)<dist) & mask4;
        Q=P4*Q;
        mask4=(Q.row(2)>0)&mask4;
        mask4=(Q.row(2)<dist)&mask4;

        mask1=mask1.t();
        mask2=mask2.t();
        mask3=mask3.t();
        mask4=mask4.t();

        /* If _mask is given, then use it to filter outliers. */
        if (!_mask.empty()) {
            Mat mask=_mask.getMat();
            CV_Assert(mask.size()==mask1.size());
            bitwise_and(mask,mask1,mask1);
            bitwise_and(mask,mask2,mask2);
            bitwise_and(mask,mask3,mask3);
            bitwise_and(mask,mask4,mask4);
        }
        if (_mask.empty()&&_mask.needed()) {
            _mask.create(mask1.size(),CV_8U);
        }
        CV_Assert(_R.needed()&& _t.needed());
        _R.create(3,3,R1.type());
        _t.create(3,1,t.type());

        int good1=countNonZero(mask1);
        int good2=countNonZero(mask2);
        int good3=countNonZero(mask3);
        int good4=countNonZero(mask4);

        if (good1>=good2&&good1>=good3&&good1>=good4) {
            R1.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask1.copyTo(_mask);
            return good1;
        }
        else if (good2>=good1&&good2>=good3&&good2>=good4) {
            R2.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask2.copyTo(_mask);
            return good2;
        }
        else if (good3>=good1&&good3>=good2&&good3>=good4) {
            t=-t;
            R1.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask3.copyTo(_mask);
            return good3;
        }
        else {
            t=-t;
            R2.copyTo(_R);
            t.copyTo(_t);
            if (_mask.needed()) mask4.copyTo(_mask);
            return good4;
        }
    }
}
/* solve R|t using opencv------------------------------------------------------
 * args   :  voopt_t *opt      I  visual odometry options
 *           match_set_t *mf   I  feature list
 *           double *dT        O  rotation and translation
 *           double *ratio     O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * ----------------------------------------------------------------------------*/
extern int solveRt(const voopt_t *opt,const insstate_t *ins,const match_set_t *mf,
                   double *dT,double *ratio)
{
    vector<cv::Point2f> ll,rr;
    double R[9],t[3],Rt[9];

    trace(3,"solveRt:\n");

    if (mf->n<8) {
        trace(3,"solve R|t using opencv fail\n");
        return 0;
    }
    for (int i=0;i<mf->n;i++) {
        ll.push_back(cv::Point2f(float((mf->data[i].up-ins->ox)/ins->fx),float((mf->data[i].vp-ins->oy)/ins->fy)));
        rr.push_back(cv::Point2f(float((mf->data[i].uc-ins->ox)/ins->fx),float((mf->data[i].vc-ins->oy)/ins->fy)));
    }
    cv::Mat mask;
    cv::Mat E=cv::findFundamentalMat(ll,rr,cv::FM_RANSAC,0.3/460,0.99,mask);
    cv::Mat cameraMatrix=(cv::Mat_<double>(3,3)<<1,0,0,0,1,0,0,0,1);

    cv::Mat rot,trans;
    int inlier_cnt=cv::recoverPose(E,ll,rr,cameraMatrix,rot,trans,mask);

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) R[j+i*3]=rot.at<double>(j,i);
        t[i]=trans.at<double>(i,0);
    }
    rt2tf(R,t,dT);
    if (ratio) {
        *ratio=(double)inlier_cnt/mf->n;
    }
    trace(3,"dT=\n");
    tracemat(3,dT,4,4,12,5);

    if (inlier_cnt>12) return 1;
    else {
        trace(3,"solve R|t using opencv fail\n");
        return 0;
    }
}
/* solve R|t using opencv------------------------------------------------------
 * args   :  voopt_t *opt      I  visual odometry options
 *           match_set_t *mf   I  feature list
 *           double *dT        O  rotation and translation
 *           double *ratio     O  ratio of inliers
 * return : 1 (ok) or 0 (fail)
 * ----------------------------------------------------------------------------*/
extern int solveRt5p(const voopt_t *opt,const insstate_t *ins,const match_set_t *mf,
                     double *dT,double *ratio)
{
    vector<cv::Point2f> ll,rr;
    Point2d pp={0};
    double R[9],t[3],Rt[9];

    trace(3,"solveRt5p:\n");

    if (mf->n<8) {
        trace(3,"solve R|t using opencv fail\n");
        return 0;
    }
    for (int i=0;i<mf->n;i++) {
        ll.push_back(cv::Point2f(float((mf->data[i].up-ins->ox)/ins->fx),float((mf->data[i].vp-ins->oy)/ins->fy)));
        rr.push_back(cv::Point2f(float((mf->data[i].uc-ins->ox)/ins->fx),float((mf->data[i].vc-ins->oy)/ins->fy)));
    }
    cv::Mat mask;
    cv::Mat E=findEssentialMat5p(ll,rr,1.0,pp,cv::FM_RANSAC,0.99,0.99,mask);

    cv::Mat rot,trans;
    int inlier_cnt=recoverPose5p(E,ll,rr,rot,trans,1.0,pp,mask);

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) R[j+i*3]=rot.at<double>(j,i);
        t[i]=trans.at<double>(i,0);
    }
    rt2tf(R,t,dT); matinv(dT,4);
    if (ratio) {
        *ratio=(double)inlier_cnt/mf->n;
    }
    trace(3,"dT=\n");
    tracemat(3,dT,4,4,12,5);

    if (inlier_cnt>12) return 1;
    else {
        trace(3,"solve R|t using opencv fail\n");
        return 0;
    }
}
