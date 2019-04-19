
#ifndef VISO_MONO_H
#define VISO_MONO_H

#include "vo.h"

class VisualOdometryMono : public VisualOdometry {

public:

    /* monocular-specific parameters (mandatory: height,pitch) */
    struct parameters : public VisualOdometry::parameters {
        double                      height;           /* camera height above ground (meters) */
        double                      pitch;            /* camera pitch (rad, negative=pointing down) */
        int32_t                     ransac_iters;     /* number of RANSAC iterations */
        double                      inlier_threshold; /* fundamental matrix inlier threshold */
        double                      motion_threshold; /* directly return false on small motions */
        parameters () {
            height          =1.0;
            pitch           =0.0;
            ransac_iters    =2000;
            inlier_threshold=0.00001;
            motion_threshold=100.0;
        }
    };

    /* constructor, takes as inpute a parameter structure */
    VisualOdometryMono(parameters param);  
  
    /* deconstructor */
    ~VisualOdometryMono();
  
    /* process a new image, pushs the image back to an internal ring buffer.
     * valid motion estimates are available after calling process for two times.
     * inputs: I ......... pointer to rectified image (uint8, row-aligned)
     *         dims[0] ... width of I
     *         dims[1] ... height of I
     *         dims[2] ... bytes per line (often equal to width)
     *         replace ... replace current image with I, without copying last current
     *                     image to previous image internally. this option can be used
     *                     when small/no motions are observed to obtain Tr_delta wrt
     *                     an older coordinate system / time step than the previous one.
     * output: returns false if motion too small or an error occured
     * */
    bool process (uint8_t *I,int32_t* dims,bool replace=false);

private:

    template<class T> struct idx_cmp {
        idx_cmp(const T arr) : arr(arr) {}
        bool operator()(const size_t a, const size_t b) const {return arr[a]<arr[b];}
        const T arr;
    };

    std::vector<double>  estimateMotion(std::vector<Matcher::p_match> p_matched);
    Matrix               smallerThanMedian (Matrix &X,double &median);
    bool                 normalizeFeaturePoints (std::vector<Matcher::p_match> &p_matched,Matrix &Tp,Matrix &Tc);
    void                 fundamentalMatrix (const std::vector<Matcher::p_match> &p_matched,const std::vector<int32_t> &active,Matrix &F);
    void                 EtoRt(Matrix &E,Matrix &K,std::vector<Matcher::p_match> &p_matched,Matrix &X,Matrix &R,Matrix &t);
    int32_t              triangulateChieral (std::vector<Matcher::p_match> &p_matched,Matrix &K,Matrix &R,Matrix &t,Matrix &X);
    std::vector<int32_t> getInlier (std::vector<Matcher::p_match> &p_matched,Matrix &F);
  
    /* parameters */
    parameters param;
};

#endif

