#pragma once

#include <opencv2/opencv.hpp>
#include <vector>
#include <iostream>
#include <ctime>
#include <map>
#include <carvig.h>
using namespace std;
using namespace cv;

#ifdef USE_GPU
#include <opencv2/cudafeatures2d.hpp>
	using cuda::GpuMat;
#endif

#define THRESH_FACTOR 6

#define add_data_func_delc(data_type,add_data_type)            \
    static int add_##data_type##_##add_data_type(              \
                             data_type *ivec,                  \
                             const add_data_type* n)           \
    {                                                          \
        add_data_type *obs_data;                               \
        if (ivec->nmax<=ivec->n) {                             \
            if (ivec->nmax<=0) ivec->nmax=10*2;                \
            else ivec->nmax*=2;                                \
            if (!(obs_data=(add_data_type*)realloc(ivec->data, \
                sizeof(add_data_type)*ivec->nmax))) {          \
                free(ivec->data);                              \
                ivec->data=NULL; ivec->n=ivec->nmax=0;         \
            }                                                  \
            ivec->data=obs_data;                               \
        }                                                      \
        ivec->data[ivec->n++]=*n;                              \
        return 1;                                              \
    }                                                          \

add_data_func_delc(match_set,match_point)

class gms_matcher
{
public:
    // OpenCV Keypoints & Correspond Image Size & Nearest Neighbor Matches
    gms_matcher(const vector<KeyPoint> &vkp1, const Size size1, const vector<KeyPoint> &vkp2, const Size size2, const vector<DMatch> &vDMatches)
    {
        // Input initialize
        NormalizePoints(vkp1, size1, mvP1);
        NormalizePoints(vkp2, size2, mvP2);
        mNumberMatches = vDMatches.size();
        ConvertMatches(vDMatches, mvMatches);

        // Grid initialize
        mGridSizeLeft = Size(20, 20);
        mGridNumberLeft = mGridSizeLeft.width * mGridSizeLeft.height;
    };
    ~gms_matcher() {};

private:
    // 8 possible rotation and each one is 3 X 3
    const int mRotationPatterns[8][9] = {
            {1,2,3,
            4,5,6,
            7,8,9},

            {4,1,2,
            7,5,3,
            8,9,6},

            {7,4,1,
            8,5,2,
            9,6,3},

            {8,7,4,
            9,5,1,
            6,3,2},

            {9,8,7,
            6,5,4,
            3,2,1},

            {6,9,8,
            3,5,7,
            2,1,4},

            {3,6,9,
            2,5,8,
            1,4,7},

            {2,3,6,
            1,5,9,
            4,7,8}
    };
    // 5 level scales
    const double mScaleRatios[5] = { 1.0, 1.0 / 2, 1.0 / sqrt(2.0), sqrt(2.0), 2.0 };

    // Normalized Points
    vector<Point2f> mvP1, mvP2;

    // Matches
    vector<pair<int, int>> mvMatches;

    // Number of Matches
    size_t mNumberMatches;

    // Grid Size
    Size mGridSizeLeft, mGridSizeRight;
    int mGridNumberLeft;

    // Index	  : left grid idx
    // map.first  : right grid idx
    // map.second : how many matches from idx_left to idx_right
    vector<map<int, int>> mMotionStatistics;

    //
    vector<int> mNumberPointsInPerCellLeft;

    // Inldex  : grid_idx_left
    // Value   : grid_idx_right
    vector<int> mCellPairs;

    // Every Matches has a cell-pair
    // first  : grid_idx_left
    // second : grid_idx_right
    vector<pair<int, int> > mvMatchPairs;

    // Inlier Mask for output
    vector<bool> mvbInlierMask;

public:
    // Get Inlier Mask
    // Return number of inliers
    int GetInlierMask(vector<bool> &vbInliers, bool WithScale = false, bool WithRotation = false);

private:

    // Normalize Key Points to Range(0 - 1)
    void NormalizePoints(const vector<KeyPoint> &kp, const Size &size, vector<Point2f> &npts) {
        const size_t numP = kp.size();
        const int width   = size.width;
        const int height  = size.height;
        npts.resize(numP);

        for (size_t i = 0; i < numP; i++)
        {
            npts[i].x = kp[i].pt.x / width;
            npts[i].y = kp[i].pt.y / height;
        }
    }
    // Convert OpenCV DMatch to Match (pair<int, int>)
    void ConvertMatches(const vector<DMatch> &vDMatches, vector<pair<int, int>> &vMatches) {
        vMatches.resize(mNumberMatches);
        for (size_t i = 0; i < mNumberMatches; i++)
        {
            vMatches[i] = pair<int, int>(vDMatches[i].queryIdx, vDMatches[i].trainIdx);
        }
    }
    int GetGridIndexLeft(const Point2f &pt, int type) {
        int x, y;

        if (type == 1) {
            x = floor(pt.x * mGridSizeLeft.width);
            y = floor(pt.y * mGridSizeLeft.height);
        }
        if (type == 2) {
            x = floor(pt.x * mGridSizeLeft.width + 0.5);
            y = floor(pt.y * mGridSizeLeft.height);
        }
        if (type == 3) {
            x = floor(pt.x * mGridSizeLeft.width);
            y = floor(pt.y * mGridSizeLeft.height + 0.5);
        }
        if (type == 4) {
            x = floor(pt.x * mGridSizeLeft.width + 0.5);
            y = floor(pt.y * mGridSizeLeft.height + 0.5);
        }
        if (x >= mGridSizeLeft.width || y >= mGridSizeLeft.height)
        {
            return -1;
        }
        return x + y * mGridSizeLeft.width;
    }

    int GetGridIndexRight(const Point2f &pt) {
        int x = floor(pt.x * mGridSizeRight.width);
        int y = floor(pt.y * mGridSizeRight.height);

        return x + y * mGridSizeRight.width;
    }
    // Assign Matches to Cell Pairs
    void AssignMatchPairs(int GridType);

    // Verify Cell Pairs
    void VerifyCellPairs(int RotationType);

    // Get Neighbor 9
    vector<int> GetNB9(const int idx, const Size& GridSize) {
        vector<int> NB9(9, -1);

        int idx_x = idx % GridSize.width;
        int idx_y = idx / GridSize.width;

        for (int yi = -1; yi <= 1; yi++)
        {
            for (int xi = -1; xi <= 1; xi++)
            {
                int idx_xx = idx_x + xi;
                int idx_yy = idx_y + yi;

                if (idx_xx < 0 || idx_xx >= GridSize.width || idx_yy < 0 || idx_yy >= GridSize.height)
                    continue;

                NB9[xi + 4 + yi * 3] = idx_xx + idx_yy * GridSize.width;
            }
        }
        return NB9;
    }
    // Run
    int run(int RotationType, int Scale);
};
int gms_matcher::GetInlierMask(vector<bool> &vbInliers, bool WithScale, bool WithRotation) {

    int max_inlier = 0;

    if (!WithScale && !WithRotation)
    {
        max_inlier = run(1, 0);
        vbInliers = mvbInlierMask;
        return max_inlier;
    }

    if (WithRotation && WithScale)
    {
        for (int Scale = 0; Scale < 5; Scale++)
        {
            for (int RotationType = 1; RotationType <= 8; RotationType++)
            {
                int num_inlier = run(RotationType, Scale);

                if (num_inlier > max_inlier)
                {
                    vbInliers = mvbInlierMask;
                    max_inlier = num_inlier;
                }
            }
        }
        return max_inlier;
    }
    if (WithRotation && !WithScale)
    {

        for (int RotationType = 1; RotationType <= 8; RotationType++)
        {
            int num_inlier = run(RotationType, 0);

            if (num_inlier > max_inlier)
            {
                vbInliers = mvbInlierMask;
                max_inlier = num_inlier;
            }
        }
        return max_inlier;
    }
    if (!WithRotation && WithScale)
    {
        for (int Scale = 0; Scale < 5; Scale++)
        {
            int num_inlier = run(1, Scale);

            if (num_inlier > max_inlier)
            {
                vbInliers = mvbInlierMask;
                max_inlier = num_inlier;
            }

        }
        return max_inlier;
    }
    return max_inlier;
}
void gms_matcher::AssignMatchPairs(int GridType) {

    for (size_t i = 0; i < mNumberMatches; i++)
    {
        Point2f &lp = mvP1[mvMatches[i].first];
        Point2f &rp = mvP2[mvMatches[i].second];

        const int lgidx = mvMatchPairs[i].first = GetGridIndexLeft(lp, GridType);
        const int rgidx = mvMatchPairs[i].second = GetGridIndexRight(rp);

        mMotionStatistics[lgidx][rgidx] ++;
        mNumberPointsInPerCellLeft[lgidx]++;
    }
}
void gms_matcher::VerifyCellPairs(int RotationType) {

    const int *CurrentRP = mRotationPatterns[RotationType - 1];

    for (int i = 0; i < mGridNumberLeft; i++)
    {
        if (mMotionStatistics[i].empty())
        {
            mCellPairs[i] = -1;
            continue;
        }
        int max_number = 0;
        for (auto &p : mMotionStatistics[i])
        {
            if (p.second > max_number) {
                mCellPairs[i] = p.first;
                max_number = p.second;
            }
        }
        int idx_grid_rt = mCellPairs[i];

        vector<int> NB9_lt = GetNB9(i, mGridSizeLeft);
        vector<int> NB9_rt = GetNB9(idx_grid_rt, mGridSizeRight);

        int score = 0;
        double thresh = 0;
        int numpair = 0;

        for (size_t j = 0; j < 9; j++)
        {
            int ll = NB9_lt[j];
            int rr = NB9_rt[CurrentRP[j] - 1];
            if (ll == -1 || rr == -1)	continue;

            score += mMotionStatistics[ll][rr];
            thresh += mNumberPointsInPerCellLeft[ll];
            numpair++;
        }
        thresh = THRESH_FACTOR * sqrt(thresh / numpair);

        if (score < thresh)
            mCellPairs[i] = -2;
    }
}
int gms_matcher::run(int RotationType, int Scale) {

    mvbInlierMask.assign(mNumberMatches, false);
    for (int GridType = 1; GridType <= 4; GridType++)
    {
        // Set Scale
        mGridSizeRight.width  = mGridSizeLeft.width  * mScaleRatios[Scale];
        mGridSizeRight.height = mGridSizeLeft.height * mScaleRatios[Scale];

        // initialize
        mMotionStatistics.assign(mGridNumberLeft, map<int, int>());
        mCellPairs.assign(mGridNumberLeft, -1);
        mNumberPointsInPerCellLeft.assign(mGridNumberLeft, 0);
        mvMatchPairs.assign(mNumberMatches, pair<int, int>(0, 0));

        AssignMatchPairs(GridType);
        VerifyCellPairs(RotationType);

        // Mark inliers
        for (size_t i = 0; i < mNumberMatches; i++)
        {
            if (mCellPairs[mvMatchPairs[i].first] == mvMatchPairs[i].second)
            {
                mvbInlierMask[i] = true;
            }
        }
    }
    int num_inlier=0;
    for (int i=0;i<mvbInlierMask.size();i++) {
        if (mvbInlierMask[i]) num_inlier++;
    }
    return num_inlier;
}
// utility
inline Mat DrawInlier(Mat &src1, Mat &src2, vector<KeyPoint> &kpt1, vector<KeyPoint> &kpt2, vector<DMatch> &inlier, int type) {
    const int height = max(src1.rows, src2.rows);
    const int width = src1.cols + src2.cols;
    Mat output(height, width, CV_8UC3, Scalar(0, 0, 0));
    src1.copyTo(output(Rect(0, 0, src1.cols, src1.rows)));
    src2.copyTo(output(Rect(src1.cols, 0, src2.cols, src2.rows)));

    if (type == 1)
    {
        for (size_t i = 0; i < inlier.size(); i++)
        {
            Point2f left = kpt1[inlier[i].queryIdx].pt;
            Point2f right = (kpt2[inlier[i].trainIdx].pt + Point2f((float)src1.cols, 0.f));
            line(output, left, right, Scalar(0, 255, 255));
        }
    }
    else if (type == 2)
    {
        for (size_t i = 0; i < inlier.size(); i++)
        {
            Point2f left = kpt1[inlier[i].queryIdx].pt;
            Point2f right = (kpt2[inlier[i].trainIdx].pt + Point2f((float)src1.cols, 0.f));
            line(output, left, right, Scalar(255, 0, 0));
        }
        for (size_t i = 0; i < inlier.size(); i++)
        {
            Point2f left = kpt1[inlier[i].queryIdx].pt;
            Point2f right = (kpt2[inlier[i].trainIdx].pt + Point2f((float)src1.cols, 0.f));
            circle(output, left, 1, Scalar(0, 255, 255), 2);
            circle(output, right, 1, Scalar(0, 255, 0), 2);
        }
    }
    return output;
}
void GmsMatch(Mat &img1,vector<KeyPoint> &kp1,vector<KeyPoint> &kp2, Mat &img2,vector<DMatch>& matches_gms){
    Mat d1, d2;
    vector<DMatch> matches_all;

    Ptr<FeatureDetector> detector = FeatureDetector::create("BRISK");
    Ptr<DescriptorExtractor> descriptor_extractor = DescriptorExtractor::create("BRISK");
    Ptr<DescriptorMatcher> descriptor_matcher = DescriptorMatcher::create("BruteForce");

    detector->detect(img1,kp1);
    detector->detect(img2,kp2);

    Mat descriptors1,descriptors2;
    descriptor_extractor->compute(img1,kp1,descriptors1);
    descriptor_extractor->compute(img2,kp2,descriptors2);

    descriptor_matcher->match(descriptors1, descriptors2, matches_all);

    // GMS filter
    int num_inliers = 0;
    std::vector<bool> vbInliers;
    gms_matcher gms(kp1,img1.size(), kp2,img2.size(), matches_all);
    num_inliers = gms.GetInlierMask(vbInliers, false, true);

    cout << "Get total " << num_inliers << " matches." << endl;

    // draw matches
    for (size_t i = 0; i < vbInliers.size(); ++i)
    {
        if (vbInliers[i] == true)
        {
            matches_gms.push_back(matches_all[i]);
        }
    }
    Mat show = DrawInlier(img1, img2, kp1, kp2, matches_gms, 1);
    imshow("show", show);
    waitKey(2);
}
static int img2Img(const img_t *img,Mat &Img)
{
    int i,j;
    for (i=0;i<Img.rows;i++) for (j=0;j<Img.cols;j++) {
            Img.at<cv::Vec3b>(i,j)[0]=img->data[i*img->w+j];
            Img.at<cv::Vec3b>(i,j)[1]=img->data[i*img->w+j];
            Img.at<cv::Vec3b>(i,j)[2]=img->data[i*img->w+j];
        }
    return 1;
}
extern int getgmsmatches(const voopt_t *opt,match_set *match,const img_t *Ip,const img_t *Ic)
{
    if (Ip==NULL||Ic==NULL) return 0;
    
    vector<DMatch> matches_gms;
    Mat I1(Ip->h,Ip->w,CV_8UC3),I2(Ip->h,Ip->w,CV_8UC3);
    match_point point;
    vector<KeyPoint> kp1,kp2;

    img2Img(Ip,I1);
    img2Img(Ic,I2);
    GmsMatch(I1,kp1,kp2,I2,matches_gms);

    match->n=0;
    for (int i=0;i<matches_gms.size();i++)
    {
        point.up=kp1[matches_gms[i].queryIdx].pt.x;
        point.vp=kp1[matches_gms[i].queryIdx].pt.y;
        point.uc=kp2[matches_gms[i].trainIdx].pt.x;
        point.vc=kp2[matches_gms[i].trainIdx].pt.y;

        add_match_set_match_point(match,&point);
    }
    return match->n;
}