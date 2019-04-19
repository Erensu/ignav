/*------------------------------------------------------------------------------
 * ins-flow-track.cc : klt feature point track functions
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
 * history : 2017/10/17 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants -----------------------------------------------------------------*/
#define KLT_BOOL int
#define KLT_THRES            5
#define MAX_KERNEL_WIDTH 	 71
#define OUTPUT_PPM           1

#define SWAP_ME(X,Y) {temp=(X);(X)=(Y);(Y)=temp;}

/* type definitions ----------------------------------------------------------*/
typedef struct  {                /* float image data type */
    int ncols;                   /* number of image data rows */
    int nrows;                   /* number of image data cols */
    float *data;                 /* float image data */
} floatimg_t,*pfloatimg_t;

typedef struct {                 /* feature point data type */
    double time;                 /* current feature point track time */
    float x,y;                   /* image coordinate of feature point */
    int valid;                   /* valid flag for feature point */
    int status;                  /* feature point track status */
    unsigned long int id;        /* feature point id */
    unsigned long int count;     /* feature track counts */
    double ts;                   /* start track time */

    /* for affine mapping */
    pfloatimg_t aff_img;
    pfloatimg_t aff_img_gradx;
    pfloatimg_t aff_img_grady;
    float aff_x,aff_y;
    float aff_Axx,aff_Ayx;
    float aff_Axy,aff_Ayy;
} feature_t,*pfeature_t;

typedef struct {                 /* KLT feature point track options type */
    int mindist;		    	 /* min distance b/w features */
    int window_width,window_height;
    KLT_BOOL sequentialMode;	 /* whether to save most recent image to save time */

    /* can set to TRUE manually, but don't set to */
    /* FALSE manually */
    KLT_BOOL smoothBeforeSelecting;	/* whether to smooth image before */
    KLT_BOOL writeInternalImages;	/* selecting features: whether to write internal images */
    KLT_BOOL lighting_insensitive;  /* tracking features:whether to normalize for gain and bias (not in original algorithm) */

    int min_eigenvalue;		        /* smallest eigenvalue allowed for selecting */
    float min_determinant;	        /* th for determining lost */
    float min_displacement;	        /* th for stopping tracking when pixel changes little */
    int max_iterations;		        /* th for stopping tracking when too many iterations */
    float max_residue;		        /* th for stopping tracking when residue is large */
    float grad_sigma;
    float smooth_sigma_fact;
    float pyramid_sigma_fact;
    float step_factor;              /* size of Newton steps; 2.0 comes from equations, 1.0 seems to avoid overshooting */
    int nSkippedPixels;		        /* # of pixels skipped when finding features */
    int borderx;			        /* border in which features will not be found */
    int bordery;
    int nPyramidLevels;	          	/* computed from search_ranges */
    int subsampling;		        /* 		" */

    /* for affine mapping */
    int affine_window_width,affine_window_height;
    int affineConsistencyCheck; /* whether to evaluates the consistency of features with affine mapping
                                   -1 = don't evaluates the consistency
                                    0 = evaluates the consistency of features with translation mapping
                                    1 = evaluates the consistency of features with similarity mapping
                                    2 = evaluates the consistency of features with affine mapping
                                */
    int affine_max_iterations;
    float affine_max_residue;
    float affine_min_displacement;
    float affine_max_displacement_differ; /* th for the difference between the displacement calculated
                                             by the affine tracker and the frame to frame tracker in pel
                                             */
} kltopt_t;

typedef struct {                    /* feature point track option type */
    kltopt_t kltopt;                /* KLT options */
    int track_nc,track_nf;          /* tracking length of camera list/max number of tracking feature points */
    int wth,hgt;                    /* width and height of image data */
} trackopt_t;

typedef struct  {                   /* KLT feature track context type */
    trackopt_t trackopt;            /* feature point track options */

    unsigned long int last_id;      /* last feature point id */

    /* User must not touch these for KLT */
    void *pyramid_last;
    void *pyramid_last_gradx;
    void *pyramid_last_grady;
    unsigned char* img[2];          /* store pointer of image buffer data {0: precious, 1: current} */
} tracking_context_t;

typedef struct  {                    /* feature list data type */
    int n;                           /* number of feature points */
    pfeature_t *feature;             /* feature list data */
} featurelist_t,*pfeaturelist_t;

typedef struct  {
    int subsampling;
    int nLevels;
    pfloatimg_t *img;
    int *ncols,*nrows;
} pyramid_t,*ppyramid_t;

/* type definitions ----------------------------------------------------------*/
typedef struct  {
    int width;
    float data[MAX_KERNEL_WIDTH];
} ConvolutionKernel;
typedef enum {SELECTING_ALL, REPLACING_SOME} selectionMode;

/* global variable -----------------------------------------------------------*/
static const int mindist                  =10;
static const int window_size              =7;
static const int min_eigenvalue           =1;
static const float min_determinant        =0.01;
static const float min_displacement       =0.1;
static const int max_iterations           =10;
static const float max_residue            =10.0;
static const float grad_sigma             =1.0;
static const float smooth_sigma_fact      =0.1;
static const float pyramid_sigma_fact     =0.9;
static const float step_factor            =1.0;
static const KLT_BOOL sequentialMode      =false;
static const KLT_BOOL lighting_insensitive=false;

/* for affine mapping*/
static const int affineConsistencyCheck    =1;
static const int affine_window_size        =15;
static const int affine_max_iterations     =10;
static const float affine_max_residue      =10.0;
static const float affine_min_displacement =0.02;
static const float affine_max_displacement_differ=1.5;

static const KLT_BOOL smoothBeforeSelecting=true;
static const KLT_BOOL writeInternalImages  =false;
static const int search_range              =15;
static const int nSkippedPixels            =0;
static float sigma_last                    =-10.0;

/* Kernels */
static ConvolutionKernel gauss_kernel;
static ConvolutionKernel gaussderiv_kernel;

/* KLT track context */
static tracking_context_t *kltrack=NULL;

/*write float image to PGM ---------------------------------------------------
* args:    pfloatimg_t *img  IO  float image data
*          char *filename    O   write image data file path
* return: none
* ---------------------------------------------------------------------------*/
static void writeFloatImageToPGM(pfloatimg_t img,char *filename)
{
    int npixs =img->ncols*img->nrows;
    float mmax=-999999.9f,mmin=999999.9f;
    float fact;
    float *ptr;
    unsigned char *byteimg,*ptrout;
    int i;

    /* calculate minimum and maximum values of float image */
    ptr=img->data;
    for (i=0;i<npixs;i++)  {
        mmax=MAX(mmax,*ptr);
        mmin=MIN(mmin,*ptr);
        ptr++;
    }
    /* allocate memory to hold converted image */
    byteimg=(unsigned char *)malloc(npixs*sizeof(unsigned char));

    /* Convert image from float to uchar */
    fact=255.0f/(mmax-mmin);
    ptr =img->data;
    ptrout=byteimg;
    for (i=0;i<npixs;i++) {
        *ptrout++=(unsigned char)((*ptr++-mmin)*fact);
    }
    /* write uchar image to PGM */
    pgmWriteFile(filename,byteimg,img->ncols,img->nrows);

    /* free memory */
    free(byteimg);
}
static void writeFeatureListToPPM(pfeaturelist_t featurelist,
                                  unsigned char *greyimg,int ncols,int nrows,
                                  char *filename)
{
    int nbytes=ncols*nrows*sizeof(char);
    unsigned char *redimg,*grnimg,*bluimg;
    int offset;
    int x,y,xx,yy;
    int i;

    /* Allocate memory for component images */
    redimg=(unsigned char *)malloc(nbytes);
    grnimg=(unsigned char *)malloc(nbytes);
    bluimg=(unsigned char *)malloc(nbytes);
    if (redimg==NULL||grnimg==NULL||bluimg==NULL) return;

    memcpy(redimg,greyimg,nbytes);
    memcpy(grnimg,greyimg,nbytes);
    memcpy(bluimg,greyimg,nbytes);

    /* Overlay features in red */
    for (i=0;i<featurelist->n;i++)
        if (featurelist->feature[i]->valid>=0)  {
            x=(int)(featurelist->feature[i]->x+0.5);
            y=(int)(featurelist->feature[i]->y+0.5);
            for (yy=y-1;yy<=y+1;yy++)
                for (xx=x-1;xx<=x+1;xx++)
                    if (xx>=0&&yy>=0&&xx<ncols&&yy<nrows)  {
                        offset=yy*ncols+xx;
                        *(redimg+offset)=255;
                        *(grnimg+offset)=0;
                        *(bluimg+offset)=0;
                    }
        }

    /* Write to PPM file */
    ppmWriteFileRGB(filename,redimg,grnimg,bluimg,ncols,nrows);

    /* Free memory */
    free(redimg);
    free(grnimg);
    free(bluimg);
}
/*----------------------------------------------------------------------------
 * _quicksort
 * Replacement for qsort().  Computing time is decreased by taking
 * advantage of specific knowledge of our array (that there are
 * three ints associated with each point).
 *
 * This routine generously provided by
 *      Manolis Lourakis <lourakis@csi.forth.gr>
 *
 * NOTE: The results of this function may be slightly different from
 * those of qsort().  This is due to the fact that different sort
 * algorithms have different behaviours when sorting numbers with the
 * same value: Some leave them in the same relative positions in the
 * array, while others change their relative positions. For example,
 * if you have the array [c d b1 a b2] with b1=b2, it may be sorted as
 * [a b1 b2 c d] or [a b2 b1 c d].
 *----------------------------------------------------------------------------*/
#define SWAP3(list, i, j)               \
{                                       \
     register int *pi, *pj, tmp;        \
     pi=list+3*(i); pj=list+3*(j);      \
                                        \
     tmp=*pi;                           \
     *pi++=*pj;                         \
     *pj++=tmp;                         \
                                        \
     tmp=*pi;                           \
     *pi++=*pj;                         \
     *pj++=tmp;                         \
                                        \
     tmp=*pi;                           \
     *pi=*pj;                           \
     *pj=tmp;                           \
}
static void quicksort(int *pointlist, int n)
{
    unsigned int i,j,ln,rn;

    while (n>1) {
        SWAP3(pointlist,0,n/2);
        for (i=0,j=n;;) {
            do --j;
            while (pointlist[3*j+2]<pointlist[2]);
            do
                ++i;
            while (i<j&&pointlist[3*i+2]>pointlist[2]);
            if (i>=j) break;
            SWAP3(pointlist,i,j);
        }
        SWAP3(pointlist,j,0);
        ln=j;
        rn=n-++j;
        if (ln<rn) {
            quicksort(pointlist,ln);
            pointlist+=3*j;
            n=rn;
        }
        else {
            quicksort(pointlist+3*j,rn);
            n=ln;
        }
    }
}
#undef SWAP3
/*----------------------------------------------------------------------------
 * Used by qsort (in _KLTSelectGoodFeatures) to determine
 * which feature is better.
 * By switching the '>' with the '<', qsort is fooled into sorting
 * in descending order.
 *----------------------------------------------------------------------------*/
#ifdef KLT_USE_QSORT
static int _comparePoints(const void *a, const void *b)
{
    int v1=*(((int *)a)+2);
    int v2=*(((int *)b)+2);

    if      (v1>v2)  return -1;
    else if (v1<v2)  return 1;
    else return 0;
}
#endif
/*----------------------------------------------------------------------------*/
static void fillFeaturemap(int x, int y,unsigned char *featuremap,
                           int mindist,int ncols,int nrows)
{
    register int ix,iy;

    for (iy=y-mindist;iy<=y+mindist;iy++)
        for (ix=x-mindist;ix<=x+mindist;ix++) {
            if (ix>=0&&ix<ncols&&iy>=0&&iy<nrows) featuremap[iy*ncols+ix]=1;
        }
}
/*----------------------------------------------------------------------------
 * Removes features that are within close proximity to better features.
 *
 * INPUTS
 * featurelist:  A list of features.  The nFeatures property
 *               is used.
 *
 * OUTPUTS
 * featurelist:  Is overwritten.  Nearby "redundant" features are removed.
 *               Writes -1's into the remaining elements.
 *
 * RETURNS
 * The number of remaining features.
 *----------------------------------------------------------------------------*/
static void enforceMinimumDistance(
        int *pointlist,              /* featurepoints */
        int npoints,                 /* number of featurepoints */
        pfeaturelist_t featurelist,  /* features */
        int ncols, int nrows,        /* size of images */
        int mindist,                 /* min. dist b/w features */
        int min_eigenvalue,          /* min. eigenvalue */
        KLT_BOOL overwriteAllFeatures,
        unsigned long int *last_id)  /* last feature point id */
{
    register int indx;         /* Index into features */
    register int x,y,val;      /* Location and trackability of pixel under consideration */
    unsigned char *featuremap; /* Boolean array recording proximity of features */
    int *ptr;

    trace(3,"enforceMinimumDistance:\n");

    /* Cannot add features with an eigenvalue less than one */
    if (min_eigenvalue<1) min_eigenvalue=1;

    /* Allocate memory for feature map and clear it */
    featuremap=(unsigned char *)malloc(ncols*nrows*sizeof(unsigned char));
    memset(featuremap,0,ncols*nrows);

    /* Necessary because code below works with (mindist-1) */
    mindist--;

    /* If we are keeping all old good features, then add them to the featuremap */
    if (!overwriteAllFeatures)
        for (indx=0;indx<featurelist->n;indx++)
            if (featurelist->feature[indx]->status>=0) {
                x=(int)featurelist->feature[indx]->x;
                y=(int)featurelist->feature[indx]->y;
                fillFeaturemap(x,y,featuremap,mindist,ncols,nrows);
            }

    /* For each feature point, in descending order of importance, do ... */
    ptr =pointlist;
    indx=0;
    while (1) {

        /* If we can't add all the points, then fill in the rest
           of the featurelist with -1's */
        if (ptr>=pointlist+3*npoints) {
            while (indx<featurelist->n) {
                if (overwriteAllFeatures||featurelist->feature[indx]->valid<0) {
                    featurelist->feature[indx]->x=-1;
                    featurelist->feature[indx]->y=-1;
                    featurelist->feature[indx]->valid  =KLT_NOT_FOUND;
                    featurelist->feature[indx]->aff_img=NULL;
                    featurelist->feature[indx]->aff_img_gradx=NULL;
                    featurelist->feature[indx]->aff_img_grady=NULL;
                    featurelist->feature[indx]->aff_x=-1.0f;
                    featurelist->feature[indx]->aff_y=-1.0f;
                    featurelist->feature[indx]->aff_Axx=1.0;
                    featurelist->feature[indx]->aff_Ayx=0.0;
                    featurelist->feature[indx]->aff_Axy=0.0;
                    featurelist->feature[indx]->aff_Ayy=1.0;
                }
                indx++;
            }
            break;
        }
        x  =*ptr++;
        y  =*ptr++;
        val=*ptr++;

        /* Ensure that feature is in-bounds */
        assert(x>=0);
        assert(x<ncols);
        assert(y>=0);
        assert(y<nrows);

        /* find index of lose feature point */
        while (!overwriteAllFeatures&&indx<featurelist->n
               &&featurelist->feature[indx]->status>=0) {
            indx++;
        }
        /* we only want n feature points */
        if (indx>=featurelist->n) break;

        /* If no neighbor has been selected, and if the minimum
           eigenvalue is large enough, then add feature to the current list */
        if (!featuremap[y*ncols+x]&&val>=min_eigenvalue)  {
            featurelist->feature[indx]->x    =(float)x;
            featurelist->feature[indx]->y    =(float)y;
            featurelist->feature[indx]->valid=val;
            featurelist->feature[indx]->status=KLT_INIT;
            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;
            featurelist->feature[indx]->aff_x=-1.0f;
            featurelist->feature[indx]->aff_y=-1.0f;
            featurelist->feature[indx]->aff_Axx=1.0;
            featurelist->feature[indx]->aff_Ayx=0.0;
            featurelist->feature[indx]->aff_Axy=0.0;
            featurelist->feature[indx]->aff_Ayy=1.0;

            featurelist->feature[indx]->id=(*last_id)++;
            indx++;

            /* Fill in surrounding region of feature map, but
               make sure that pixels are in-bounds
               */
            fillFeaturemap(x,y,featuremap,mindist,ncols,nrows);
        }
    }
    /* Free feature map  */
    free(featuremap);
}
/* createArray2D--------------------------------------------------------------*/
static void** createArray2D(int ncols, int nrows, int nbytes)
{
    char **tt;
    register int i;

    tt=(char**)malloc(nrows*sizeof(void*)+ncols*nrows*nbytes);
    if (tt==NULL) {
        trace(2,"createArray2D) Out of memory");
    }
    for (i=0;i<nrows;i++) {
        tt[i]=((char*)tt)+(nrows*sizeof(void *)+i*ncols*nbytes);
    }
    return (void**)tt;
}
/* computeKernels-------------------------------------------------------------*/
static void computeKernels(float sigma,ConvolutionKernel *gauss,
                           ConvolutionKernel *gaussderiv)
{
    const float factor=0.01;   /* for truncating tail */
    register int i;

    trace(3,"computeKernels:\n");

    assert(MAX_KERNEL_WIDTH%2==1);
    assert(sigma>=0.0);

    /* Compute kernels, and automatically determine widths */
    {
        const int hw=MAX_KERNEL_WIDTH/2;
        float max_gauss=1.0,max_gaussderiv=(float)(sigma*exp(-0.5));

        /* Compute gauss and deriv */
        for (i=-hw;i<=hw;i++) {
            gauss->data[i+hw]     =(float)exp(-i*i/(2*sigma*sigma));
            gaussderiv->data[i+hw]=-i*gauss->data[i+hw];
        }
        /* Compute widths */
        gauss->width=MAX_KERNEL_WIDTH;
        for (i=-hw;fabs(gauss->data[i+hw]/max_gauss)<factor;
             i++,gauss->width-=2);

        gaussderiv->width=MAX_KERNEL_WIDTH;
        for (i=-hw;fabs(gaussderiv->data[i+hw]/max_gaussderiv)<factor;
             i++,gaussderiv->width-=2);

        if (gauss->width==MAX_KERNEL_WIDTH||gaussderiv->width==MAX_KERNEL_WIDTH) {

            trace(2,"(_computeKernels) MAX_KERNEL_WIDTH %d is too small for "
                    "a sigma of %f",MAX_KERNEL_WIDTH,sigma);
        }
    }
    /* Shift if width less than MAX_KERNEL_WIDTH */
    for (i=0;i<gauss->width;i++) {
        gauss->data[i]=gauss->data[i+(MAX_KERNEL_WIDTH-gauss->width)/2];
    }
    for (i=0;i<gaussderiv->width;i++) {
        gaussderiv->data[i]=gaussderiv->data[i+(MAX_KERNEL_WIDTH-gaussderiv->width)/2];
    }
    /* Normalize gauss and deriv */
    {
        const int hw=gaussderiv->width/2;
        float den;

        den=0.0;
        for (i=0;i<gauss->width;i++) den+=gauss->data[i];
        for (i=0;i<gauss->width;i++) gauss->data[i]/=den;
        den=0.0;
        for (i=-hw;i<=hw;i++) den-=i*gaussderiv->data[i+hw];
        for (i=-hw;i<=hw;i++) gaussderiv->data[i+hw]/=den;
    }
    sigma_last=sigma;
    return;
}
/*----------------------------------------------------------------------------*/
static void KLTGetKernelWidths(float sigma,int *gauss_width,
                               int *gaussderiv_width)
{
    trace(3,"KLTGetKernelWidths:\n");

    computeKernels(sigma,&gauss_kernel,&gaussderiv_kernel);
    *gauss_width=gauss_kernel.width;
    *gaussderiv_width=gaussderiv_kernel.width;
}
/*----------------------------------------------------------------------------*/
static float KLTComputeSmoothSigma(tracking_context_t *tc)
{
    return tc->trackopt.kltopt.smooth_sigma_fact*MAX(tc->trackopt.kltopt.window_width,
                                                     tc->trackopt.kltopt.window_height);
}
/*----------------------------------------------------------------------------*/
static float pyramidSigma(tracking_context_t* tc)
{
    return tc->trackopt.kltopt.pyramid_sigma_fact*tc->trackopt.kltopt.subsampling;
}
/*----------------------------------------------------------------------------*/
static void klt_UpdateTCBorder(tracking_context_t* tc)
{
    register float val;
    register int pyramid_gauss_hw;
    register int smooth_gauss_hw;
    register int gauss_width, gaussderiv_width;
    register int num_levels=tc->trackopt.kltopt.nPyramidLevels;
    register int n_invalid_pixels;
    register int window_hw;
    register int ss=tc->trackopt.kltopt.subsampling;
    register int ss_power;
    register int border;
    register int i;

    trace(3,"klt_UpdateTCBorder:\n");

    window_hw=MAX(tc->trackopt.kltopt.window_width,
                  tc->trackopt.kltopt.window_height)/2;

    /* Find widths of convolution windows */
    KLTGetKernelWidths(KLTComputeSmoothSigma(tc),&gauss_width,&gaussderiv_width);
    smooth_gauss_hw=gauss_width/2;

    KLTGetKernelWidths(pyramidSigma(tc),&gauss_width,&gaussderiv_width);
    pyramid_gauss_hw=gauss_width/2;

    /* Compute the # of invalid pixels at each level of the pyramid.
       n_invalid_pixels is computed with respect to the ith level
       of the pyramid.  So, e.g., if n_invalid_pixels = 5 after
       the first iteration, then there are 5 invalid pixels in
       level 1, which translated means 5*subsampling invalid pixels
       in the original level 0.
       */
    n_invalid_pixels=smooth_gauss_hw;
    for (i=1;i<num_levels;i++)  {
        val=((float) n_invalid_pixels+pyramid_gauss_hw)/ss;
        n_invalid_pixels=(int)(val+0.99);  /* Round up */
    }
    /* ss_power = ss^(num_levels-1) */
    ss_power=1;
    for (i=1;i<num_levels;i++) ss_power*=ss;

    /* Compute border by translating invalid pixels back into */
    /* original image */
    border=(n_invalid_pixels+window_hw)*ss_power;

    tc->trackopt.kltopt.borderx=border;
    tc->trackopt.kltopt.bordery=border;
    return;
}
/* KLT change TCP-Pyram ID----------------------------------------------------
 * args:    tracking_context_t *tc  IO  klt tracking context
 *          int search_range        I   klt search range
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt_changeTCPyramid(tracking_context_t* tc,int search_range)
{
    register float window_halfwidth,subsampling;

    trace(3,"klt_changeTCPyramid:\n");

    window_halfwidth=MIN(tc->trackopt.kltopt.window_width,tc->trackopt.kltopt.window_height)/2.0f;
    subsampling=((float)search_range)/window_halfwidth;

    if (subsampling<1.0) {		/* 1.0 = 0+1 */
        tc->trackopt.kltopt.nPyramidLevels=1;
    }
    else if (subsampling<=3.0) {	/* 3.0 = 2+1 */
        tc->trackopt.kltopt.nPyramidLevels=2;
        tc->trackopt.kltopt.subsampling=2;
    }
    else if (subsampling<= 5.0) {	/* 5.0 = 4+1 */
        tc->trackopt.kltopt.nPyramidLevels=2;
        tc->trackopt.kltopt.subsampling=4;
    }
    else if (subsampling<=9.0) {	/* 9.0 = 8+1 */
        tc->trackopt.kltopt.nPyramidLevels=2;
        tc->trackopt.kltopt.subsampling = 8;
    }
    else {
        /* The following lines are derived from the formula:
           search_range =
           window_halfwidth * \sum_{i=0}^{nPyramidLevels-1} 8^i,
           which is the same as:
           search_range =
           window_halfwidth * (8^nPyramidLevels - 1)/(8 - 1).
           Then, the value is rounded up to the nearest integer.
           */
        float val=(float)(log(7.0*subsampling+1.0)/log(8.0));
        tc->trackopt.kltopt.nPyramidLevels=(int)(val+0.99);
        tc->trackopt.kltopt.subsampling=8;
    }
    return;
}
/* create klt tracking context------------------------------------------------*/
static tracking_context_t* klt_create_tracking_context()
{
    tracking_context_t *tc;

    trace(3,"klt_create_tracking_context:\n");

    /* allocate memory */
    tc=(tracking_context_t*)malloc(sizeof(tracking_context_t));

    /* set values to default values */
    tc->trackopt.kltopt.mindist              =mindist;
    tc->trackopt.kltopt.window_width         =window_size;
    tc->trackopt.kltopt.window_height        =window_size;
    tc->trackopt.kltopt.sequentialMode       =sequentialMode;
    tc->trackopt.kltopt.smoothBeforeSelecting=smoothBeforeSelecting;
    tc->trackopt.kltopt.writeInternalImages  =writeInternalImages;
    tc->trackopt.kltopt.lighting_insensitive =lighting_insensitive;
    tc->trackopt.kltopt.min_eigenvalue       =min_eigenvalue;
    tc->trackopt.kltopt.min_determinant      =min_determinant;
    tc->trackopt.kltopt.max_iterations       =max_iterations;
    tc->trackopt.kltopt.min_displacement     =min_displacement;
    tc->trackopt.kltopt.max_residue          =max_residue;
    tc->trackopt.kltopt.grad_sigma           =grad_sigma;
    tc->trackopt.kltopt.smooth_sigma_fact    =smooth_sigma_fact;
    tc->trackopt.kltopt.pyramid_sigma_fact   =pyramid_sigma_fact;
    tc->trackopt.kltopt.step_factor          =step_factor;
    tc->trackopt.kltopt.nSkippedPixels       =nSkippedPixels;
    tc->pyramid_last         =NULL;
    tc->pyramid_last_gradx   =NULL;
    tc->pyramid_last_grady   =NULL;
    tc->last_id=0;

    /* for affine mapping */
    tc->trackopt.kltopt.affineConsistencyCheck        =affineConsistencyCheck;
    tc->trackopt.kltopt.affine_window_width           =affine_window_size;
    tc->trackopt.kltopt.affine_window_height          =affine_window_size;
    tc->trackopt.kltopt.affine_max_iterations         =affine_max_iterations;
    tc->trackopt.kltopt.affine_max_residue            =affine_max_residue;
    tc->trackopt.kltopt.affine_min_displacement       =affine_min_displacement;
    tc->trackopt.kltopt.affine_max_displacement_differ=affine_max_displacement_differ;

    /* Change nPyramidLevels and subsampling */
    klt_changeTCPyramid(tc,search_range);

    /* Update border, which is dependent upon  */
    /* smooth_sigma_fact, pyramid_sigma_fact, window_size, and subsampling */
    klt_UpdateTCBorder(tc);

    /* Check window size (and correct if necessary) */
    if (tc->trackopt.kltopt.window_width%2!=1) {
        tc->trackopt.kltopt.window_width=tc->trackopt.kltopt.window_width+1;
        trace(3,"Window width must be odd. Changing to %d.\n",tc->trackopt.kltopt.window_width);
    }
    if (tc->trackopt.kltopt.window_height%2!=1) {
        tc->trackopt.kltopt.window_height=tc->trackopt.kltopt.window_height+1;
        trace(3,"Window height must be odd. Changing to %d.\n",
              tc->trackopt.kltopt.window_height);
    }
    if (tc->trackopt.kltopt.window_width<3) {
        tc->trackopt.kltopt.window_width=3;
        trace(3,"Window width must be at least three.\n Changing to %d.\n",
              tc->trackopt.kltopt.window_width);
    }
    if (tc->trackopt.kltopt.window_height<3) {
        tc->trackopt.kltopt.window_height=3;
        trace(3,"Window height must be at least three.\n Changing to %d.\n",
              tc->trackopt.kltopt.window_height);
    }
    return tc;
}
/* create feature points list-------------------------------------------------
 * args:    int nFeatures  I  number of feature points
 * return: pointer of feature points list
 * ---------------------------------------------------------------------------*/
static featurelist_t* klt_create_featurelist(int nFeatures)
{
    featurelist_t* fl;
    pfeature_t first;
    int nbytes=sizeof(featurelist_t)+nFeatures*sizeof(pfeature_t)+
               nFeatures*sizeof(feature_t);
    register int i;

    trace(3,"klt_create_featurelist:\n");

    /* allocate memory for feature list */
    fl=(featurelist_t*)malloc(nbytes);

    /* set parameters */
    fl->n=nFeatures;

    /* set pointers */
    fl->feature=(pfeature_t*)(fl+1);
    first=(pfeature_t)(fl->feature+nFeatures);
    for (i=0;i<nFeatures;i++) {
        fl->feature[i]=first+i;
        fl->feature[i]->aff_img=NULL;  /* initialization fixed by Sinisa Segvic */
        fl->feature[i]->aff_img_gradx=NULL;
        fl->feature[i]->aff_img_grady=NULL;
        fl->feature[i]->count=0;
        fl->feature[i]->id=0; fl->feature[i]->ts=fl->feature[i]->time=0.0;
        fl->feature[i]->valid =0;
        fl->feature[i]->status=KLT_INIT;
    }
    /* return feature list */
    return fl;
}
/* given a pointer to image data (probably unsigned chars), copy data to a
 * float image.
 * args:    unsigned char* img    I  input image data
 *          int nclos,nrows       I  size of image
 *          floatimg_t *floatimg  O  output float image data
 * return: none
 * ---------------------------------------------------------------------------*/
static void tofloatImage(unsigned char *img,int ncols, int nrows,
                         floatimg_t* floatimg)
{
    trace(3,"tofloatImage:\n");

    unsigned char *ptrend=img+ncols*nrows;
    float *ptrout=floatimg->data;

    /* Output image must be large enough to hold result */
    assert(floatimg->ncols>=ncols);
    assert(floatimg->nrows>=nrows);

    floatimg->ncols=ncols;
    floatimg->nrows=nrows;

    while (img<ptrend) *ptrout++=(float)*img++;
}
/* create a float image-------------------------------------------------------
 * args:    int ncol,nrow  I  size of image data
 * return: pointer of image data
 * ---------------------------------------------------------------------------*/
static floatimg_t* create_floatimg(int ncols,int nrows)
{
    floatimg_t* floatimg;
    int nbytes=sizeof(floatimg_t)+ncols*nrows*sizeof(float);

    floatimg=(floatimg_t*)malloc(nbytes);
    if (floatimg==NULL) {
        trace(2,"(create_floatimg) Out of memory");
    }
    floatimg->ncols=ncols;
    floatimg->nrows=nrows;
    floatimg->data =(float*)(floatimg+1);
    return floatimg;
}
/* free float image data------------------------------------------------------
 * args:    floatimg_t *floatimg  IO  input/output float image data
 * return: none
 * ---------------------------------------------------------------------------*/
static void free_floatimg(floatimg_t* floatimg)
{
    trace(3,"free_floatimg:\n");
    free(floatimg); return;
}
/*----------------------------------------------------------------------------*/
static void convolveImageHoriz(floatimg_t* imgin,ConvolutionKernel kernel,
                               floatimg_t* imgout)
{
    float *ptrrow=imgin->data;           /* Points to row's first pixel */
    register float *ptrout=imgout->data, /* Points to next output pixel */
            *ppp;

    register float sum;
    register int radius=kernel.width/2;
    register int ncols =imgin->ncols,nrows=imgin->nrows;
    register int i,j,k;

    /* Kernel width must be odd */
    assert(kernel.width%2==1);

    /* Must read from and write to different images */
    assert(imgin!=imgout);

    /* Output image must be large enough to hold result */
    assert(imgout->ncols>=imgin->ncols);
    assert(imgout->nrows>=imgin->nrows);

    /* For each row, do ... */
    for (j=0;j<nrows;j++)  {

        /* Zero leftmost columns */
        for (i=0;i<radius;i++) *ptrout++=0.0;

        /* Convolve middle columns with kernel */
        for (;i<ncols-radius;i++)  {
            ppp=ptrrow+i-radius;
            sum=0.0;
            for (k=kernel.width-1;k>=0;k--) sum+=*ppp++*kernel.data[k];
            *ptrout++ = sum;
        }
        /* Zero rightmost columns */
        for (;i<ncols;i++) *ptrout++=0.0;
        ptrrow+=ncols;
    }
    return;
}
/*----------------------------------------------------------------------------*/
static void convolveImageVert(floatimg_t* imgin,ConvolutionKernel kernel,
                              floatimg_t* imgout)
{
    float *ptrcol=imgin->data;            /* Points to row's first pixel */
    register float *ptrout=imgout->data,  /* Points to next output pixel */
            *ppp;
    register float sum;
    register int radius=kernel.width/2;
    register int ncols =imgin->ncols,nrows=imgin->nrows;
    register int i,j,k;

    /* Kernel width must be odd */
    assert(kernel.width%2==1);

    /* Must read from and write to different images */
    assert(imgin!=imgout);

    /* Output image must be large enough to hold result */
    assert(imgout->ncols>=imgin->ncols);
    assert(imgout->nrows>=imgin->nrows);

    /* For each column, do ... */
    for (i=0;i<ncols;i++)  {

        /* Zero topmost rows */
        for (j=0;j<radius;j++) {
            *ptrout=0.0; ptrout+=ncols;
        }
        /* Convolve middle rows with kernel */
        for (;j<nrows-radius;j++) {
            ppp=ptrcol+ncols*(j-radius);
            sum=0.0;
            for (k=kernel.width-1;k>=0;k--) {
                sum+=*ppp*kernel.data[k];
                ppp+=ncols;
            }
            *ptrout=sum; ptrout+=ncols;
        }
        /* Zero bottommost rows */
        for (;j<nrows;j++) {
            *ptrout=0.0; ptrout+=ncols;
        }
        ptrcol++;
        ptrout-=nrows*ncols-1;
    }
    return;
}
/*----------------------------------------------------------------------------*/
static void convolveSeparate(floatimg_t* imgin,
                             ConvolutionKernel horiz_kernel,
                             ConvolutionKernel vert_kernel,
                             floatimg_t* imgout)
{
    /* Create temporary image */
    floatimg_t* tmpimg;
    tmpimg=create_floatimg(imgin->ncols, imgin->nrows);

    /* do convolution */
    convolveImageHoriz(imgin,horiz_kernel,tmpimg);
    convolveImageVert (tmpimg,vert_kernel,imgout);

    /* free memory */
    free_floatimg(tmpimg);
    return;
}
/* compute  of image----------------------------------------------------------
 * args:    floatimg_t* img    I  input float image data
 *          float sigma        I  sigma for process image data
 *          floatimg_t* gradx  O  x-gradients of input image data
 *          floatimg_t* grady  O  y-gradients of input image data
 * return: none
 * ---------------------------------------------------------------------------*/
static void compute_gradients(floatimg_t *img, float sigma, floatimg_t *gradx,
                              floatimg_t *grady)
{
    /* Output images must be large enough to hold result */
    assert(gradx->ncols>=img->ncols);
    assert(gradx->nrows>=img->nrows);
    assert(grady->ncols>=img->ncols);
    assert(grady->nrows>=img->nrows);

    /* Compute kernels, if necessary */
    if (fabs(sigma-sigma_last)>0.05) {
        computeKernels(sigma,&gauss_kernel,&gaussderiv_kernel);
    }
    convolveSeparate(img,gaussderiv_kernel,gauss_kernel,gradx);
    convolveSeparate(img,gauss_kernel,gaussderiv_kernel,grady);
    return;
}
/* free pyramid---------------------------------------------------------------
 * args:    pyramid_t*  IO  pyramid
 * return: none
 * ---------------------------------------------------------------------------*/
static void free_pyramid(pyramid_t *pyramid)
{
    register int i;

    /* free images */
    for (i=0;i<pyramid->nLevels;i++) {
        free_floatimg(pyramid->img[i]);
    }
    /* free structure */
    free(pyramid); return;
}
/* free KLT tracking context--------------------------------------------------
 * args:    tracking_context_t* context  IO  input track context
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt_free_track_context(tracking_context_t *tc)
{
    trace(3,"klt_free_track_context:\n");

    if (tc->pyramid_last)  {
        free_pyramid((ppyramid_t)tc->pyramid_last);
    }
    if (tc->pyramid_last_gradx) {
        free_pyramid((ppyramid_t)tc->pyramid_last_gradx);
    }
    if (tc->pyramid_last_grady) {
        free_pyramid((ppyramid_t)tc->pyramid_last_grady);
    }
    free(tc);
}
/* free feature points list---------------------------------------------------
 * args:    featurelist_t *flist  IO  input feature list
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt_free_featurelist(featurelist_t *flist)
{
    /* for affine mapping */
    int indx;

    trace(3,"klt_free_featurelist:\n");

    for (indx=0;indx<flist->n;indx++) {
        /* free image and gradient  */
        free_floatimg(flist->feature[indx]->aff_img);
        free_floatimg(flist->feature[indx]->aff_img_gradx);
        free_floatimg(flist->feature[indx]->aff_img_grady);
        flist->feature[indx]->aff_img = NULL;
        flist->feature[indx]->aff_img_gradx = NULL;
        flist->feature[indx]->aff_img_grady = NULL;
    }
    free(flist);
}
/*----------------------------------------------------------------------------
 * Given the three distinct elements of the symmetric 2x2 matrix
 *                     [gxx gxy]
 *                     [gxy gyy],
 * Returns the minimum eigenvalue of the matrix.
 *----------------------------------------------------------------------------*/
static float minEigenvalue(float gxx,float gxy,float gyy)
{
    return (float)((gxx+gyy-sqrt((gxx-gyy)*(gxx-gyy)+4*gxy*gxy))/2.0);
}
/*----------------------------------------------------------------------------*/
static void sortPointList(int *pointlist,int npoints)
{
#ifdef KLT_USE_QSORT
    qsort(pointlist,npoints,3*sizeof(int),_comparePoints);
#else
    quicksort(pointlist,npoints);
#endif
}
/* get smooth sigma value-----------------------------------------------------
 * args:    tracking_context_t* tc  IO  tracking context
 * return: sigma value
 * ---------------------------------------------------------------------------*/
static float ComputeSmoothSigma(tracking_context_t* tc)
{
    return (tc->trackopt.kltopt.smooth_sigma_fact*MAX(tc->trackopt.kltopt.window_width,
                                                      tc->trackopt.kltopt.window_height));
}
/* image smooth---------------------------------------------------------------
 * args:    pfloatimg_t *img    I  input image data
 *          float sigma         I  sigma for smoothing image
 *          pfloatimg_t smooth  O  smoothed image
 * return: none
 * ---------------------------------------------------------------------------*/
static void compute_smoothedimg(pfloatimg_t img,float sigma,
                                pfloatimg_t smooth)
{
    trace(3,"compute_smoothedimg:\n");

    /* output image must be large enough to hold result */
    assert(smooth->ncols>=img->ncols);
    assert(smooth->nrows>=img->nrows);

    /* compute kernel, if necessary; gauss_deriv is not used */
    if (fabs(sigma-sigma_last)>0.05) {
        computeKernels(sigma,&gauss_kernel,&gaussderiv_kernel);
    }
    convolveSeparate(img,gauss_kernel,
                     gauss_kernel,smooth);
}
/*----------------------------------------------------------------------------*/
static void KLTSelectGoodFeatures(tracking_context_t* tc, unsigned char *img,
                                  int ncols,int nrows,
                                  pfeaturelist_t featurelist,
                                  selectionMode mode)
{
    pfloatimg_t floatimg,gradx,grady;
    int window_hw,window_hh;
    int *pointlist;
    int npoints=0;
    KLT_BOOL overwriteAllFeatures=(mode==SELECTING_ALL)?
                                  true:false;
    KLT_BOOL floatimages_created=false;

    trace(3,"KLTSelectGoodFeatures:\n");

    window_hw=tc->trackopt.kltopt.window_width/2;
    window_hh=tc->trackopt.kltopt.window_height/2;

    /* Create pointlist, which is a simplified version of a featurelist, */
    /* for speed.  Contains only integer locations and values. */
    pointlist=(int*)malloc(ncols*nrows*3*sizeof(int));

    /* Create temporary images, etc. */
    if (mode==REPLACING_SOME&&
        tc->trackopt.kltopt.sequentialMode&&tc->pyramid_last!=NULL) {
        floatimg=((ppyramid_t) tc->pyramid_last)->img[0];
        gradx=((ppyramid_t)tc->pyramid_last_gradx)->img[0];
        grady=((ppyramid_t)tc->pyramid_last_grady)->img[0];
        assert(gradx!=NULL);
        assert(grady!=NULL);
    }
    else {
        floatimages_created=true;
        floatimg=create_floatimg(ncols,nrows);
        gradx   =create_floatimg(ncols,nrows);
        grady   =create_floatimg(ncols,nrows);
        if (tc->trackopt.kltopt.smoothBeforeSelecting) {
            pfloatimg_t tmpimg;
            tmpimg=create_floatimg(ncols,nrows);
            tofloatImage(img,ncols,nrows,tmpimg);

            compute_smoothedimg(tmpimg,ComputeSmoothSigma(tc),floatimg);
            free_floatimg(tmpimg);
        }
        else tofloatImage(img,ncols,nrows,floatimg);

        /* Compute gradient of image in x and y direction */
        compute_gradients(floatimg,tc->trackopt.kltopt.grad_sigma,gradx,grady);
    }
    /* write internal images */
    if (tc->trackopt.kltopt.writeInternalImages)  {
        writeFloatImageToPGM(floatimg,(char*)"kltimg_sgfrlf.pgm");
        writeFloatImageToPGM(gradx,(char*)"kltimg_sgfrlf_gx.pgm");
        writeFloatImageToPGM(grady,(char*)"kltimg_sgfrlf_gy.pgm");
    }
    /* compute trackability of each image pixel as the minimum
       of the two eigenvalues of the Z matrix */
    {
        register float gx,gy;
        register float gxx,gxy, gyy;
        register int xx,yy;
        register int *ptr;
        float val;
        unsigned int limit=1;
        int borderx=tc->trackopt.kltopt.borderx;	/* Must not touch cols */
        int bordery=tc->trackopt.kltopt.bordery;	/* lost by convolution */
        int x,y;
        int i;

        if (borderx<window_hw) borderx=window_hw;
        if (bordery<window_hh) bordery=window_hh;

        /* Find largest value of an int */
        for (i=0;i<sizeof(int);i++) limit*=256;
        limit=limit/2-1;

        /* For most of the pixels in the image, do ... */
        ptr=pointlist;
        for (y=bordery;y<nrows-bordery;y+=tc->trackopt.kltopt.nSkippedPixels+1)
            for (x=borderx;x<ncols-borderx;x+=tc->trackopt.kltopt.nSkippedPixels+1)  {

                /* Sum the gradients in the surrounding window */
                gxx=0; gxy=0; gyy=0;
                for (yy=y-window_hh;yy<=y+window_hh;yy++)
                    for (xx=x-window_hw;xx<=x+window_hw;xx++)  {
                        gx=*(gradx->data+ncols*yy+xx);
                        gy=*(grady->data+ncols*yy+xx);
                        gxx+=gx*gx;
                        gxy+=gx*gy;
                        gyy+=gy*gy;
                    }
                /* store the trackability of the pixel as the minimum
                   of the two eigenvalues */
                *ptr++=x;
                *ptr++=y;
                val=minEigenvalue(gxx,gxy,gyy);
                if (val>limit) {
                    trace(2,"minimum eigenvalue %f is "
                            "greater than the capacity of an int; setting "
                            "to maximum value",val);
                    val=(float)limit;
                }
                *ptr++=(int)val;
                npoints++;
            }
    }
    /* Sort the features  */
    sortPointList(pointlist,npoints);

    /* Check tc->mindist */
    if (tc->trackopt.kltopt.mindist<0)  {
        trace(2,"tracking context field tc->mindist is negative (%d); setting to zero",
              tc->trackopt.kltopt.mindist);
        tc->trackopt.kltopt.mindist=0;
    }
    /* enforce minimum distance between features */
    enforceMinimumDistance(
            pointlist,
            npoints,
            featurelist,
            ncols,nrows,
            tc->trackopt.kltopt.mindist,
            tc->trackopt.kltopt.min_eigenvalue,
            overwriteAllFeatures,&tc->last_id);

    /* Free memory */
    free(pointlist);
    if (floatimages_created)  {
        free_floatimg(floatimg);
        free_floatimg(gradx);
        free_floatimg(grady);
    }
}
/* select good feature points from tracking context---------------------------
 * args:    tracking_context_t* tc  IO  tracking context
 *          unsigned char* img      I   input image data
 *          int ncols,nrows         I   size of image
 *          featurelist_t* fl       O   feature points list
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt_select_goodfeatures(tracking_context_t* tc,unsigned char *img,
                                    int ncols,int nrows,
                                    featurelist_t *fl)
{
    trace(3,"klt_select_goodfeatures:\n");
    KLTSelectGoodFeatures(tc,img,ncols,nrows,fl,SELECTING_ALL);
}
/*----------------------------------------------------------------------------*/
static int count_good_features(pfeaturelist_t fl)
{
    int count=0;
    int i;

    for (i=0;i<fl->n;i++) {
        if (fl->feature[i]->status>=0) count++;
    }
    return count;
}
/* create a pyramid-----------------------------------------------------------
 * args:    int ncols,nrows  I  size of image
 *          int sunsampling  I  step of sub-sample
 *          int nlevels      I  number of pyramid levels
 * return: pointer to a pyramid
 * ---------------------------------------------------------------------------*/
static ppyramid_t klt_create_pyramid(int ncols,int nrows,int subsampling,
                                     int nlevels)
{
    ppyramid_t pyramid;
    int nbytes=sizeof(pyramid_t) +
               nlevels*sizeof(pfloatimg_t *)+
               nlevels*sizeof(int)+
               nlevels*sizeof(int);
    int i;

    if (subsampling!=2&&subsampling!=4&&
        subsampling!=8&&subsampling!=16&&subsampling!=32) {
        trace(2,"Pyramid's subsampling must be either 2, 4, 8, 16, or 32\n");
    }
    /* Allocate memory for structure and set parameters */
    pyramid=(ppyramid_t)malloc(nbytes);
    if (pyramid==NULL) {
        trace(2,"(klt_create_pyramid) Out of memory");
    }
    /* set parameters */
    pyramid->subsampling=subsampling;
    pyramid->nLevels=nlevels;
    pyramid->img  =(pfloatimg_t*)(pyramid+1);
    pyramid->ncols=(int*)(pyramid->img+nlevels);
    pyramid->nrows=(int*)(pyramid->ncols+nlevels);

    /* Allocate memory for each level of pyramid and assign pointers */
    for (i=0;i<nlevels;i++)  {
        pyramid->img[i]  =create_floatimg(ncols,nrows);
        pyramid->ncols[i]=ncols;
        pyramid->nrows[i]=nrows;
        ncols/=subsampling;
        nrows/=subsampling;
    }
    return pyramid;
}
/* compute pyramid of image---------------------------------------------------
 * args:    pfloatimg_t *img     IO  image float data
 *          ppyramid_t *pyramid  IO  pyramid of image data
 *          float sigma_fact     I   sigma fact for processing image data
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt_compute_pyramid(pfloatimg_t img,ppyramid_t pyramid,
                                float sigma_fact)
{
    pfloatimg_t currimg,tmpimg;
    int ncols=img->ncols,nrows=img->nrows;
    int subsampling=pyramid->subsampling;
    int subhalf=subsampling/2;
    float sigma=subsampling*sigma_fact;  /* empirically determined */
    int oldncols;
    int i,x,y;

    if (subsampling!=2&&subsampling!=4&&subsampling!=8&&subsampling!=16&&subsampling!=32){
        trace(2,"Pyramid's subsampling must "
                "be either 2, 4, 8, 16, or 32");
        return;
    }
    assert(pyramid->ncols[0]==img->ncols);
    assert(pyramid->nrows[0]==img->nrows);

    /* Copy original image to level 0 of pyramid */
    memcpy(pyramid->img[0]->data,img->data,ncols*nrows*sizeof(float));

    currimg=img;
    for (i=1;i<pyramid->nLevels;i++) {
        tmpimg=create_floatimg(ncols,nrows);
        compute_smoothedimg(currimg,sigma,tmpimg);

        /* Subsample */
        oldncols=ncols;
        ncols/=subsampling; nrows/=subsampling;
        for (y=0;y<nrows;y++)
            for (x=0;x<ncols;x++) {
                pyramid->img[i]->data[y*ncols+x]=tmpimg->data[(subsampling*y+subhalf)*oldncols+(subsampling*x+subhalf)];
            }
        /* Reassign current image */
        currimg=pyramid->img[i];
        free_floatimg(tmpimg);
    }
    return;
}
/*-----------------------------------------------------------------------------*/
static KLT_BOOL outOfBounds(float x,float y,int ncols,int nrows,int borderx,
                            int bordery)
{
    return (x<borderx||x>ncols-1-borderx||y<bordery||y>nrows-1-bordery);
}
/*----------------------------------------------------------------------------*/
static void am_getSubFloatImage(pfloatimg_t img,float x, float y,
                                pfloatimg_t window)
{
    register int hw=window->ncols/2, hh = window->nrows/2;
    int x0=(int)x;
    int y0=(int)y;
    float *windata=window->data;
    int offset;
    register int i,j;

    assert(x0-hw>=0);
    assert(y0-hh>=0);
    assert(x0+hw<=img->ncols);
    assert(y0+hh<=img->nrows);

    /* copy values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++)  {
            offset=(j+y0)*img->ncols+(i+x0);
            *windata++=*(img->data+offset);
        }
}
/*----------------------------------------------------------------------------*/
static float* allocateFloatWindow(int width,int height)
{
    float* fw;

    fw=(float*)malloc(width*height*sizeof(float));
    if (fw==NULL) {
        trace(2,"(_allocateFloatWindow) Out of memory.\n");
    }
    return fw;
}
/*----------------------------------------------------------------------------*/
static float **am_matrix(long nr, long nc)
{
    float **m;
    int a;
    m=(float**)malloc((size_t)(nr*sizeof(float*)));
    m[0]=(float*)malloc((size_t)((nr*nc)*sizeof(float)));
    for (a=1;a<nr;a++) m[a]=m[a-1]+nc;
    return m;
}
/*----------------------------------------------------------------------------*/
static void am_free_matrix(float **m)
{
    free(m[0]); free(m);
}
/*----------------------------------------------------------------------------
 * Solves the 2x2 matrix equation
 *         [gxx gxy] [dx] = [ex]
 *         [gxy gyy] [dy] = [ey]
 * for dx and dy.
 *
 * Returns KLT_TRACKED on success and KLT_SMALL_DET on failure
 *----------------------------------------------------------------------------*/
static int solveEquation(float gxx, float gxy, float gyy,float ex, float ey,
                         float small,float *dx,float *dy)
{
    float det=gxx*gyy-gxy*gxy;
    if (det<small) return KLT_SMALL_DET;

    *dx=(gyy*ex-gxy*ey)/det;
    *dy=(gxx*ey-gxy*ex)/det;
    return KLT_TRACKED;
}
/*----------------------------------------------------------------------------*/
static float sumAbsFloatWindow(float* fw,int width,int height)
{
    float sum=0.0;
    int w;

    for (;height>0;height--) {
        for (w=0;w<width;w++) sum+=(float)fabs(*fw++);
    }
    return sum;
}
/*----------------------------------------------------------------------------
 * Given a point (x,y) in an image, computes the bilinear interpolated
 * gray-level value of the point in the image.
 *----------------------------------------------------------------------------*/
static float interpolate(float x,float y,pfloatimg_t img)
{
    int xt=(int)x;  /* coordinates of top-left corner */
    int yt=(int)y;
    float ax=x-xt;
    float ay=y-yt;
    float *ptr=img->data+(img->ncols*yt)+xt;

    assert (xt>=0&&yt>=0&&xt<=img->ncols-2&&yt<=img->nrows-2);
    return ((1-ax)*(1-ay)**ptr+ax*(1-ay)**(ptr+1)+
            (1-ax)*ay**(ptr+(img->ncols))+ax*ay**(ptr+(img->ncols)+1));
}
/*----------------------------------------------------------------------------
 * aligns the gradients with the affine transformed window
 *----------------------------------------------------------------------------*/
static void am_getGradientWinAffine(
        pfloatimg_t in_gradx,
        pfloatimg_t in_grady,
        float x, float y,      /* center of window*/
        float Axx, float Ayx , float Axy, float Ayy,  /* affine mapping */
        int width, int height, /* size of window */
        float* out_gradx,      /* output */
        float* out_grady)      /* output */
{
    register int hw=width/2,hh=height/2;
    register int i,j;
    float mi,mj;

    /* Compute values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++)  {
            mi=Axx*i+Axy*j;
            mj=Ayx*i+Ayy*j;
            *out_gradx++=interpolate(x+mi,y+mj,in_gradx);
            *out_grady++=interpolate(x+mi,y+mj,in_grady);
        }
}
/*----------------------------------------------------------------------------*/
static int am_gauss_jordan_elimination(float **a, int n, float **b, int m)
{
    /* re-implemented from Numerical Recipes in C */
    int *indxc,*indxr,*ipiv;
    int i,j,k,l,ll;
    float big,dum,pivinv,temp;
    int col=0;
    int row=0;

    indxc=(int *)malloc((size_t) (n*sizeof(int)));
    indxr=(int *)malloc((size_t) (n*sizeof(int)));
    ipiv =(int *)malloc((size_t) (n*sizeof(int)));
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv[j]!=1)
                for (k=0;k<n;k++) {
                    if (ipiv[k]==0) {
                        if (fabs(a[j][k])>=big) {
                            big=(float) fabs(a[j][k]);
                            row=j;
                            col=k;
                        }
                    } else if (ipiv[k]>1) return KLT_SMALL_DET;
                }
        ++(ipiv[col]);
        if (row != col) {
            for (l=0;l<n;l++) SWAP_ME(a[row][l],a[col][l])
            for (l=0;l<m;l++) SWAP_ME(b[row][l],b[col][l])
        }
        indxr[i]=row;
        indxc[i]=col;
        if (a[col][col]==0.0) return KLT_SMALL_DET;
        pivinv=1.0f/a[col][col];
        a[col][col]=1.0;
        for (l=0;l<n;l++) a[col][l] *= pivinv;
        for (l=0;l<m;l++) b[col][l] *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll!=col) {
                dum=a[ll][col];
                a[ll][col]=0.0;
                for (l=0;l<n;l++) a[ll][l]-=a[col][l]*dum;
                for (l=0;l<m;l++) b[ll][l]-=b[col][l]*dum;
            }
    }
    for (l=n-1;l>=0;l--) {
        if (indxr[l]!=indxc[l])
            for (k=0;k<n;k++) SWAP_ME(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free(ipiv);
    free(indxr);
    free(indxc);

    return KLT_TRACKED;
}
/*----------------------------------------------------------------------------
 * Given two images and the window center in both images,
 * aligns the images with the window and computes the difference
 * between the two overlaid images using the affine mapping.
 *       A =  [ Axx Axy]
 *            [ Ayx Ayy]
 *----------------------------------------------------------------------------*/
static void am_computeIntensityDifferenceAffine(
        pfloatimg_t img1,        /* images */
        pfloatimg_t img2,
        float x1, float y1,      /* center of window in 1st img */
        float x2, float y2,      /* center of window in 2nd img */
        float Axx, float Ayx , float Axy, float Ayy,    /* affine mapping */
        int width, int height,   /* size of window */
        float* imgdiff)          /* output */
{
    register int hw=width/2,hh=height/2;
    float g1,g2;
    register int i,j;
    float mi,mj;

    /* Compute values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++) {
            g1=interpolate(x1+i,y1+j,img1);
            mi=Axx*i+Axy*j;
            mj=Ayx*i+Ayy*j;
            g2=interpolate(x2+mi,y2+mj,img2);
            *imgdiff++=g1-g2;
        }
}
/*----------------------------------------------------------------------------*/
static void am_compute6by6GradientMatrix(float* gradx,float* grady,
                                         int width, int height,float **T)
{
    register int hw=width/2,hh=height/2;
    register int i,j;
    float gx,gy,gxx,gxy,gyy,x,y,xx,xy,yy;

    /* Set values to zero */
    for (j=0;j<6;j++) {
        for (i=j;i<6;i++) {
            T[j][i]=0.0;
        }
    }
    for (j=-hh;j<=hh;j++) {
        for (i=-hw;i<=hw;i++) {
            gx=*gradx++;
            gy=*grady++;
            gxx=gx*gx;
            gxy=gx*gy;
            gyy=gy*gy;
            x=(float) i;
            y=(float) j;
            xx=x*x;
            xy=x*y;
            yy=y*y;

            T[0][0]+=xx*gxx;
            T[0][1]+=xx*gxy;
            T[0][2]+=xy*gxx;
            T[0][3]+=xy*gxy;
            T[0][4]+=x *gxx;
            T[0][5]+=x *gxy;

            T[1][1]+=xx*gyy;
            T[1][2]+=xy*gxy;
            T[1][3]+=xy*gyy;
            T[1][4]+=x *gxy;
            T[1][5]+=x *gyy;

            T[2][2]+=yy*gxx;
            T[2][3]+=yy*gxy;
            T[2][4]+=y *gxx;
            T[2][5]+=y *gxy;

            T[3][3]+=yy*gyy;
            T[3][4]+=y *gxy;
            T[3][5]+=y *gyy;

            T[4][4]+=gxx;
            T[4][5]+=gxy;

            T[5][5]+=gyy;
        }
    }
    for (j=0;j<5;j++) {
        for (i=j+1;i<6;i++) T[i][j]=T[j][i];
    }
}
/*----------------------------------------------------------------------------*/
static void am_compute4by1ErrorVector(float* imgdiff,float* gradx,float* grady,
                                      int width,int height,float **e)
{
    register int hw=width/2,hh=height/2;
    register int i,j;
    register float diff,diffgradx,diffgrady;

    /* Set values to zero */
    for (i=0;i<4;i++) e[i][0]=0.0;

    /* Compute values */
    for (j=-hh;j<=hh;j++) {
        for (i=-hw;i<=hw;i++)  {
            diff = *imgdiff++;
            diffgradx=diff*(*gradx++);
            diffgrady=diff*(*grady++);
            e[0][0]+=diffgradx*i+diffgrady*j;
            e[1][0]+=diffgrady*i-diffgradx*j;
            e[2][0]+=diffgradx;
            e[3][0]+=diffgrady;
        }
    }
    for (i=0;i<4;i++) e[i][0]*=0.5;
}
/*----------------------------------------------------------------------------*/
static void am_compute6by1ErrorVector(float* imgdiff,float* gradx,float* grady,
                                      int width,int height,float **e)
{
    register int hw=width/2,hh=height/2;
    register int i,j;
    register float diff,diffgradx,diffgrady;

    /* Set values to zero */
    for (i =0;i<6;i++) e[i][0] = 0.0;

    /* Compute values */
    for (j=-hh;j<=hh;j++) {
        for (i=-hw;i<=hw;i++)  {
            diff=*imgdiff++;
            diffgradx=diff*(*gradx++);
            diffgrady=diff*(*grady++);
            e[0][0]+=diffgradx*i;
            e[1][0]+=diffgrady*i;
            e[2][0]+=diffgradx*j;
            e[3][0]+=diffgrady*j;
            e[4][0]+=diffgradx;
            e[5][0]+=diffgrady;
        }
    }
    for (i=0;i<6;i++) e[i][0]*=0.5;
}
/*----------------------------------------------------------------------------
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference
 * between the two overlaid images; normalizes for overall gain and bias.
 *----------------------------------------------------------------------------*/
static void computeIntensityDifferenceLightingInsensitive(
        pfloatimg_t img1,       /* images */
        pfloatimg_t img2,
        float x1, float y1,     /* center of window in 1st img */
        float x2, float y2,     /* center of window in 2nd img */
        int width, int height,  /* size of window */
        float* imgdiff)         /* output */
{
    register int hw = width/2, hh = height/2;
    float g1, g2, sum1_squared = 0, sum2_squared = 0;
    register int i, j;

    float sum1 = 0, sum2 = 0;
    float mean1, mean2,alpha,belta;
    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            sum1 += g1;    sum2 += g2;
            sum1_squared += g1*g1;
            sum2_squared += g2*g2;
        }
    mean1=sum1_squared/(width*height);
    mean2=sum2_squared/(width*height);
    alpha = (float) sqrt(mean1/mean2);
    mean1=sum1/(width*height);
    mean2=sum2/(width*height);
    belta = mean1-alpha*mean2;

    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            *imgdiff++ = g1- g2*alpha-belta;
        }
}
/*----------------------------------------------------------------------------
 * Given two gradients and the window center in both images,
 * aligns the gradients wrt the window and computes the sum of the two
 * overlaid gradients; normalizes for overall gain and bias.
 *----------------------------------------------------------------------------*/
static void computeGradientSumLightingInsensitive(
        pfloatimg_t gradx1,      /* gradient images */
        pfloatimg_t grady1,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        pfloatimg_t img1,        /* images */
        pfloatimg_t img2,

        float x1, float y1,      /* center of window in 1st img */
        float x2, float y2,      /* center of window in 2nd img */
        int width, int height,   /* size of window */
        float* gradx,            /* output */
        float* grady)
{
    register int hw = width/2, hh = height/2;
    float g1, g2, sum1_squared = 0, sum2_squared = 0;
    register int i, j;

    float sum1 = 0, sum2 = 0;
    float mean1, mean2, alpha;
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            sum1_squared += g1;    sum2_squared += g2;
        }
    mean1 = sum1_squared/(width*height);
    mean2 = sum2_squared/(width*height);
    alpha = (float) sqrt(mean1/mean2);

    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, gradx1);
            g2 = interpolate(x2+i, y2+j, gradx2);
            *gradx++ = g1 + g2*alpha;
            g1 = interpolate(x1+i, y1+j, grady1);
            g2 = interpolate(x2+i, y2+j, grady2);
            *grady++ = g1+ g2*alpha;
        }
}
/*----------------------------------------------------------------------------
 * Given two images and the window center in both images,
 * aligns the images wrt the window and computes the difference
 * between the two overlaid images.
 *----------------------------------------------------------------------------*/
static void computeIntensityDifference(
        pfloatimg_t img1,       /* images */
        pfloatimg_t img2,
        float x1, float y1,     /* center of window in 1st img */
        float x2, float y2,     /* center of window in 2nd img */
        int width, int height,  /* size of window */
        float* imgdiff)         /* output */
{
    register int hw = width/2, hh = height/2;
    float g1, g2;
    register int i, j;

    /* Compute values */
    for (j = -hh ; j <= hh ; j++)
        for (i = -hw ; i <= hw ; i++)  {
            g1 = interpolate(x1+i, y1+j, img1);
            g2 = interpolate(x2+i, y2+j, img2);
            *imgdiff++ = g1 - g2;
        }
}
/*----------------------------------------------------------------------------
 * Given two gradients and the window center in both images,
 * aligns the gradients wrt the window and computes the sum of the two
 * overlaid gradients.
 *----------------------------------------------------------------------------*/
static void computeGradientSum(
        pfloatimg_t gradx1,      /* gradient images */
        pfloatimg_t grady1,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        float x1, float y1,      /* center of window in 1st img */
        float x2, float y2,      /* center of window in 2nd img */
        int width, int height,   /* size of window */
        float* gradx,            /* output */
        float* grady)            /*   " */
{
    register int hw=width/2,hh=height/2;
    float g1,g2;
    register int i,j;

    /* Compute values */
    for (j=-hh;j<=hh;j++)
        for (i=-hw;i<=hw;i++) {
            g1=interpolate(x1+i,y1+j,gradx1);
            g2=interpolate(x2+i,y2+j,gradx2);
            *gradx++=g1+g2;
            g1=interpolate(x1+i,y1+j,grady1);
            g2=interpolate(x2+i,y2+j,grady2);
            *grady++=g1+g2;
        }
}
/*----------------------------------------------------------------------------*/
static void compute2by1ErrorVector(
        float* imgdiff,
        float* gradx,
        float* grady,
        int width,         /* size of window */
        int height,
        float step_factor, /* 2.0 comes from equations, 1.0 seems to avoid overshooting */
        float *ex,         /* return values */
        float *ey)
{
    register float diff;
    register int i;

    /* Compute values */
    *ex=0; *ey=0;
    for (i=0;i<width*height;i++)  {
        diff=*imgdiff++;
        *ex+=diff*(*gradx++);
        *ey+=diff*(*grady++);
    }
    *ex*=step_factor;
    *ey*=step_factor;
}
/*----------------------------------------------------------------------------*/
static void compute2by2GradientMatrix(
        float* gradx,
        float* grady,
        int width,   /* size of window */
        int height,
        float *gxx,  /* return values */
        float *gxy,
        float *gyy)

{
    register float gx,gy;
    register int i;

    /* Compute values */
    *gxx=0.0;*gxy=0.0;*gyy=0.0;
    for (i=0;i<width*height;i++) {
        gx=*gradx++;
        gy=*grady++;
        *gxx+=gx*gx;
        *gxy+=gx*gy;
        *gyy+=gy*gy;
    }
}
static void am_compute4by4GradientMatrix(
        float* gradx,
        float* grady,
        int width,   /* size of window */
        int height,
        float **T)   /* return values */
{
    register int hw=width/2,hh=height/2;
    register int i,j;
    float gx,gy,x,y;

    /* Set values to zero */
    for (j=0;j<4;j++) {
        for (i=0;i<4;i++) {
            T[j][i]=0.0;
        }
    }
    for (j=-hh;j<=hh;j++) {
        for (i=-hw;i<=hw;i++) {
            gx=*gradx++;
            gy=*grady++;
            x=(float)i;
            y=(float)j;
            T[0][0]+=(x*gx+y*gy)*(x*gx+y*gy);
            T[0][1]+=(x*gx+y*gy)*(x*gy-y*gx);
            T[0][2]+=(x*gx+y*gy)*gx;
            T[0][3]+=(x*gx+y*gy)*gy;

            T[1][1]+=(x*gy-y*gx)*(x*gy-y*gx);
            T[1][2]+=(x*gy-y*gx)*gx;
            T[1][3]+=(x*gy-y*gx)*gy;

            T[2][2]+=gx*gx;
            T[2][3]+=gx*gy;

            T[3][3]+=gy*gy;
        }
    }
    for (j=0;j<3;j++) {
        for (i=j+1;i<4;i++) {
            T[i][j]=T[j][i];
        }
    }
}
/*----------------------------------------------------------------------------
 * Tracks a feature point from the image of first occurrence to the actual image.
 *
 * RETURNS
 * KLT_SMALL_DET or KLT_LARGE_RESIDUE or KLT_OOB if feature is lost,
 * KLT_TRACKED otherwise.
 *----------------------------------------------------------------------------*/
static int am_trackFeatureAffine(
        float x1,  /* location of window in first image */
        float y1,
        float *x2, /* starting location of search in second image */
        float *y2,
        pfloatimg_t img1,
        pfloatimg_t gradx1,
        pfloatimg_t grady1,
        pfloatimg_t img2,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        int width,           /* size of window */
        int height,
        float step_factor,   /* 2.0 comes from equations, 1.0 seems to avoid overshooting */
        int max_iterations,
        float small,         /* determinant threshold for declaring KLT_SMALL_DET */
        float th,            /* displacement threshold for stopping  */
        float th_aff,
        float max_residue,   /* residue threshold for declaring KLT_LARGE_RESIDUE */
        int lighting_insensitive,  /* whether to normalize for gain and bias */
        int affine_map,      /* whether to evaluates the consistency of features with affine mapping */
        float mdd,           /* difference between the displacements */
        float *Axx, float *Ayx,
        float *Axy, float *Ayy) /* used affine mapping */
{
    float *imgdiff,*gradx,*grady;
    float gxx,gxy,gyy,ex,ey,dx,dy;
    int iteration=0;
    int status=0;
    int hw=width/2;
    int hh=height/2;
    int nc1=img1->ncols;
    int nr1=img1->nrows;
    int nc2=img2->ncols;
    int nr2=img2->nrows;
    float **a;
    float **T;
    float one_plus_eps=1.001f;   /* To prevent rounding errors */
    float old_x2=*x2;
    float old_y2=*y2;
    KLT_BOOL convergence=false;

    trace(3,"am_trackFeatureAffine:\n");

    /* Allocate memory for windows */
    imgdiff=allocateFloatWindow(width,height);
    gradx  =allocateFloatWindow(width,height);
    grady  =allocateFloatWindow(width,height);
    T=am_matrix(6,6);
    a=am_matrix(6,1);

    /* Iteratively update the window position */
    do  {
        if(!affine_map) {
            /* pure translation tracker */

            /* If out of bounds, exit loop */
            if (x1-hw <0.0f||nc1-( x1+hw)<one_plus_eps||
                *x2-hw<0.0f||nc2-(*x2+hw)<one_plus_eps||
                y1-hh <0.0f||nr1-( y1+hh)<one_plus_eps||
                *y2-hh<0.0f||nr2-(*y2+hh)<one_plus_eps) {
                status=KLT_OOB;
                break;
            }
            /* Compute gradient and difference windows */
            if (lighting_insensitive) {
                computeIntensityDifferenceLightingInsensitive(img1, img2,x1,y1,*x2,*y2,
                                                              width,height,imgdiff);
                computeGradientSumLightingInsensitive(gradx1,grady1,gradx2,grady2,
                                                      img1, img2,x1,y1,*x2,*y2,
                                                      width,height,gradx,grady);
            } else {
                computeIntensityDifference(img1,img2,x1,y1,*x2,*y2,
                                           width,height,imgdiff);
                computeGradientSum(gradx1,grady1,gradx2,grady2,
                                   x1,y1,*x2,*y2,width,height,gradx,grady);
            }
            /* Use these windows to construct matrices */
            compute2by2GradientMatrix(gradx,grady,width,height,
                                      &gxx,&gxy,&gyy);
            compute2by1ErrorVector(imgdiff,gradx,grady,width,height,step_factor,
                                   &ex,&ey);

            /* Using matrices, solve equation for new displacement */
            status=solveEquation(gxx,gxy,gyy,ex,ey,small,&dx,&dy);

            convergence=(fabs(dx)<th&&fabs(dy)<th);

            *x2+=dx;
            *y2+=dy;

        }
        else {
            /* affine tracker */
            float ul_x=*Axx*(-hw)+*Axy*  hh+*x2;  /* upper left corner */
            float ul_y=*Ayx*(-hw)+*Ayy*  hh+*y2;
            float ll_x=*Axx*(-hw)+*Axy*(-hh)+*x2;  /* lower left corner */
            float ll_y=*Ayx*(-hw)+*Ayy*(-hh)+*y2;
            float ur_x=*Axx*hw+*Axy*  hh +*x2;  /* upper right corner */
            float ur_y=*Ayx*hw+*Ayy*  hh +*y2;
            float lr_x=*Axx*hw+*Axy*(-hh)+*x2;  /* lower right corner */
            float lr_y=*Ayx*hw+*Ayy*(-hh)+*y2;

            /* If out of bounds, exit loop */
            if ( x1-hw<0.0f||nc1-(x1+hw)<one_plus_eps||
                 y1-hh<0.0f||nr1-(y1+hh)<one_plus_eps||
                 ul_x <0.0f||nc2-(ul_x )<one_plus_eps||
                 ll_x <0.0f||nc2-(ll_x )<one_plus_eps||
                 ur_x <0.0f||nc2-(ur_x )<one_plus_eps||
                 lr_x <0.0f||nc2-(lr_x )<one_plus_eps||
                 ul_y <0.0f||nr2-(ul_y )<one_plus_eps||
                 ll_y <0.0f||nr2-(ll_y )<one_plus_eps||
                 ur_y <0.0f||nr2-(ur_y )<one_plus_eps||
                 lr_y <0.0f||nr2-(lr_y )<one_plus_eps) {
                status=KLT_OOB;
                break;
            }
            am_computeIntensityDifferenceAffine(img1,img2,x1,y1,*x2,*y2,*Axx,*Ayx,*Axy,*Ayy,
                                                width,height,imgdiff);

            am_getGradientWinAffine(gradx2,grady2,*x2,*y2,*Axx,*Ayx,*Axy,*Ayy,
                                    width, height,gradx,grady);

            switch (affine_map) {
                case 1:
                    am_compute4by1ErrorVector(imgdiff, gradx,grady,width, height,a);
                    am_compute4by4GradientMatrix(gradx,grady,width,height,T);
                    status=am_gauss_jordan_elimination(T,4,a,1);

                    *Axx+=a[0][0];
                    *Ayx+=a[1][0];
                    *Ayy=*Axx;
                    *Axy=-(*Ayx);

                    dx=a[2][0];
                    dy=a[3][0];

                    break;
                case 2:
                    am_compute6by1ErrorVector(imgdiff,gradx,grady,width,height,a);
                    am_compute6by6GradientMatrix(gradx,grady,width,height,T);

                    status=am_gauss_jordan_elimination(T,6,a,1);

                    *Axx+=a[0][0];
                    *Ayx+=a[1][0];
                    *Axy+=a[2][0];
                    *Ayy+=a[3][0];

                    dx=a[4][0];
                    dy=a[5][0];
                    break;
            }
            *x2+=dx;
            *y2+=dy;

            /* old upper left corner - new upper left corner */
            ul_x-=*Axx*(-hw)+*Axy*hh+*x2;
            ul_y-=*Ayx*(-hw)+*Ayy*hh+*y2;

            /* old lower left corner - new lower left corner */
            ll_x-=*Axx*(-hw)+*Axy*(-hh)+*x2;
            ll_y-=*Ayx*(-hw)+*Ayy*(-hh)+*y2;

            /* old upper right corner - new upper right corner */
            ur_x-=*Axx*hw+*Axy*hh+*x2;
            ur_y-=*Ayx*hw+*Ayy*hh+*y2;

            /* old lower right corner - new lower right corner */
            lr_x-=*Axx*hw+*Axy*(-hh)+*x2;
            lr_y-=*Ayx*hw+*Ayy*(-hh)+*y2;

            convergence=(fabs(dx)<th&&fabs(dy)<th&&
                         fabs(ul_x)<th_aff&&fabs(ul_y)<th_aff&&
                         fabs(ll_x)<th_aff&&fabs(ll_y)<th_aff&&
                         fabs(ur_x)<th_aff&&fabs(ur_y)<th_aff&&
                         fabs(lr_x)<th_aff&&fabs(lr_y)<th_aff);
        }
        if (status==KLT_SMALL_DET) break;
        iteration++;

    } while(!convergence&&iteration<max_iterations);
    /*}  while ((fabs(dx)>=th||fabs(dy)>=th||(affine_map&&iteration<8))&&iteration<max_iterations); */
    am_free_matrix(T);
    am_free_matrix(a);

    /* Check whether window is out of bounds */
    if (*x2-hw<0.0f||nc2-(*x2+hw)<one_plus_eps ||
        *y2-hh<0.0f||nr2-(*y2+hh)<one_plus_eps)
        status=KLT_OOB;

    /* Check whether feature point has moved to much during iteration*/
    if ((*x2-old_x2)>mdd||(*y2-old_y2)>mdd )
        status = KLT_OOB;

    /* Check whether residue is too large */
    if (status==KLT_TRACKED)  {
        if (!affine_map){
            computeIntensityDifference(img1,img2,x1,y1,*x2,*y2,
                                       width,height,imgdiff);
        } else {
            am_computeIntensityDifferenceAffine(img1,img2,x1,y1,*x2,*y2,*Axx,*Ayx,*Axy,*Ayy,
                                                width,height,imgdiff);
        }
        if (sumAbsFloatWindow(imgdiff,width,height)/(width*height)>max_residue) {
            status=KLT_LARGE_RESIDUE;
        }
    }
    /* Free memory */
    free(imgdiff); free(gradx); free(grady);

    /* Return appropriate value */
    return status;
}
/*----------------------------------------------------------------------------
 * Tracks a feature point from one image to the next.
 *
 * RETURNS
 * KLT_SMALL_DET if feature is lost,
 * KLT_MAX_ITERATIONS if tracking stopped because iterations timed out,
 * KLT_TRACKED otherwise.
 *----------------------------------------------------------------------------*/
static int trackFeature(
        float x1,           /* location of window in first image */
        float y1,
        float *x2,          /* starting location of search in second image */
        float *y2,
        pfloatimg_t img1,
        pfloatimg_t gradx1,
        pfloatimg_t grady1,
        pfloatimg_t img2,
        pfloatimg_t gradx2,
        pfloatimg_t grady2,
        int width,           /* size of window */
        int height,
        float step_factor,   /* 2.0 comes from equations, 1.0 seems to avoid overshooting */
        int max_iterations,
        float small,         /* determinant threshold for declaring KLT_SMALL_DET */
        float th,            /* displacement threshold for stopping               */
        float max_residue,   /* residue threshold for declaring KLT_LARGE_RESIDUE */
        int lighting_insensitive)  /* whether to normalize for gain and bias */
{
    float* imgdiff,*gradx,*grady;
    float gxx,gxy,gyy,ex,ey,dx,dy;
    int iteration=0;
    int status;
    int hw=width/2;
    int hh=height/2;
    int nc=img1->ncols;
    int nr=img1->nrows;
    float one_plus_eps=1.001f;   /* To prevent rounding errors */

    /* Allocate memory for windows */
    imgdiff=allocateFloatWindow(width,height);
    gradx  =allocateFloatWindow(width,height);
    grady  =allocateFloatWindow(width,height);

    /* Iteratively update the window position */
    do {
        /* If out of bounds, exit loop */
        if (x1-hw<0.0f||nc-( x1+hw)<one_plus_eps||
            *x2-hw<0.0f||nc-(*x2+hw)<one_plus_eps||
            y1-hh<0.0f||nr-( y1+hh)<one_plus_eps||
            *y2-hh<0.0f||nr-(*y2+hh)<one_plus_eps) {
            status=KLT_OOB;
            break;
        }
        /* Compute gradient and difference windows */
        if (lighting_insensitive) {
            computeIntensityDifferenceLightingInsensitive(img1,img2,x1,y1,*x2,*y2,
                                                          width,height,imgdiff);
            computeGradientSumLightingInsensitive(gradx1,grady1,gradx2,grady2,
                                                  img1,img2,x1,y1,*x2,*y2,width,height,gradx,grady);
        } else {
            computeIntensityDifference(img1, img2,x1,y1,*x2,*y2,
                                       width,height,imgdiff);
            computeGradientSum(gradx1,grady1,gradx2,grady2,
                               x1,y1,*x2,*y2,width,height,gradx,grady);
        }
        /* Use these windows to construct matrices */
        compute2by2GradientMatrix(gradx,grady,width,height,
                                  &gxx,&gxy,&gyy);
        compute2by1ErrorVector(imgdiff,gradx,grady,width,height,step_factor,
                               &ex,&ey);

        /* Using matrices, solve equation for new displacement */
        status=solveEquation(gxx,gxy,gyy,ex,ey,small,&dx,&dy);
        if (status==KLT_SMALL_DET) break;

        *x2+=dx;
        *y2+=dy;
        iteration++;

    } while((fabs(dx)>=th||fabs(dy)>=th)&&iteration<max_iterations);

    /* Check whether window is out of bounds */
    if (*x2-hw<0.0f||nc-(*x2+hw)<one_plus_eps||
        *y2-hh<0.0f||nr-(*y2+hh)<one_plus_eps)
        status=KLT_OOB;

    /* Check whether residue is too large */
    if (status == KLT_TRACKED)  {
        if (lighting_insensitive)
            computeIntensityDifferenceLightingInsensitive(img1,img2,x1,y1,*x2,*y2,
                                                          width,height,imgdiff);
        else
            computeIntensityDifference(img1,img2,x1,y1,*x2,*y2,
                                       width,height,imgdiff);
        if (sumAbsFloatWindow(imgdiff,width,height)/(width*height)>max_residue)
            status=KLT_LARGE_RESIDUE;
    }
    /* Free memory */
    free(imgdiff); free(gradx); free(grady);

    /* Return appropriate value */
    if (status==KLT_SMALL_DET) {
        return KLT_SMALL_DET;
    }
    else if (status==KLT_OOB          ) return KLT_OOB;
    else if (status==KLT_LARGE_RESIDUE) return KLT_LARGE_RESIDUE;
    else if (iteration>=max_iterations) return KLT_MAX_ITERATIONS;
    else return KLT_TRACKED;
}
/* KLT tracking feature points-------------------------------------------------
 * tracks feature points from one image to the next.
 * args:    tracking_context_t *tc  IO  tracking context
 *          unsigned char* img1     I   first image data
 *          unsigned char* img2     I   second image data
 *          int ncols,nrows         I   size if image
 *          pfeaturelist_t fl       O   tracked feature points list
 * return: none
 * ---------------------------------------------------------------------------*/
static void klt(tracking_context_t* tc,unsigned char *img1,unsigned char *img2,
                int ncols,int nrows,pfeaturelist_t featurelist)
{
    pfloatimg_t tmpimg,floatimg1,floatimg2;
    ppyramid_t pyramid1,pyramid1_gradx,pyramid1_grady,pyramid2,pyramid2_gradx,pyramid2_grady;
    float subsampling=(float)tc->trackopt.kltopt.subsampling;
    float xloc,yloc,xlocout,ylocout;
    int val=-1;
    int indx,r;
    KLT_BOOL floatimg1_created=false;
    int i;

    trace(3,"(KLT) Tracking %d features in a %d by %d image...  \n",
          count_good_features(featurelist),ncols,nrows);

    /* Create temporary image */
    tmpimg=create_floatimg(ncols,nrows);

    /* Process first image by converting to float, smoothing, computing */
    /* pyramid, and computing gradient pyramids
     * */
    if (tc->trackopt.kltopt.sequentialMode&&tc->pyramid_last!=NULL) {
        pyramid1=(ppyramid_t)tc->pyramid_last;
        pyramid1_gradx=(ppyramid_t)tc->pyramid_last_gradx;
        pyramid1_grady=(ppyramid_t)tc->pyramid_last_grady;
        if (pyramid1->ncols[0]!=ncols||pyramid1->nrows[0]!=nrows)
            trace(2,"size of incoming image (%d by %d) "
                          "is different from size of previous image (%d by %d)\n",
                  ncols,nrows,pyramid1->ncols[0],pyramid1->nrows[0]);

        assert(pyramid1_gradx!=NULL);
        assert(pyramid1_grady!=NULL);
    }
    else  {
        floatimg1_created=true;
        floatimg1=create_floatimg(ncols,nrows);
        tofloatImage(img1,ncols,nrows,tmpimg);

        compute_smoothedimg(tmpimg,ComputeSmoothSigma(tc),floatimg1);
        pyramid1=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->trackopt.kltopt.nPyramidLevels);

        klt_compute_pyramid(floatimg1,pyramid1,tc->trackopt.kltopt.pyramid_sigma_fact);
        pyramid1_gradx=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->trackopt.kltopt.nPyramidLevels);
        pyramid1_grady=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->trackopt.kltopt.nPyramidLevels);

        for (i=0;i<tc->trackopt.kltopt.nPyramidLevels;i++)
            compute_gradients(pyramid1->img[i],tc->trackopt.kltopt.grad_sigma,
                              pyramid1_gradx->img[i],
                              pyramid1_grady->img[i]);
    }
    /* Do the same thing with second image */
    floatimg2=create_floatimg(ncols,nrows);
    tofloatImage(img2,ncols,nrows,tmpimg);

    compute_smoothedimg(tmpimg,ComputeSmoothSigma(tc),floatimg2);
    pyramid2=klt_create_pyramid(ncols,nrows,(int)subsampling,
                                tc->trackopt.kltopt.nPyramidLevels);

    klt_compute_pyramid(floatimg2,pyramid2,tc->trackopt.kltopt.pyramid_sigma_fact);

    pyramid2_gradx=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->trackopt.kltopt.nPyramidLevels);
    pyramid2_grady=klt_create_pyramid(ncols,nrows,(int)subsampling,tc->trackopt.kltopt.nPyramidLevels);

    for (i=0;i<tc->trackopt.kltopt.nPyramidLevels;i++) {
        compute_gradients(pyramid2->img[i],tc->trackopt.kltopt.grad_sigma,
                          pyramid2_gradx->img[i],
                          pyramid2_grady->img[i]);
    }
    /* For each feature, do ... */
    for (indx=0;indx<featurelist->n;indx++) {

        /* Only track features that are not lost */
        xloc=featurelist->feature[indx]->x;
        yloc=featurelist->feature[indx]->y;

        /* Transform location to coarsest resolution */
        for (r=tc->trackopt.kltopt.nPyramidLevels-1;r>=0;r--) {
            xloc/=subsampling;
            yloc/=subsampling;
        }
        xlocout=xloc; ylocout=yloc;

        /* Beginning with coarsest resolution, do ... */
        for (r=tc->trackopt.kltopt.nPyramidLevels-1;r>=0;r--) {

            /* Track feature at current resolution */
            xloc*=subsampling;
            yloc*=subsampling;
            xlocout*=subsampling; ylocout*=subsampling;

            val=trackFeature(xloc,yloc,&xlocout,&ylocout,
                             pyramid1->img[r],
                             pyramid1_gradx->img[r],pyramid1_grady->img[r],
                             pyramid2->img[r],
                             pyramid2_gradx->img[r],pyramid2_grady->img[r],
                             tc->trackopt.kltopt.window_width,tc->trackopt.kltopt.window_height,
                             tc->trackopt.kltopt.step_factor,
                             tc->trackopt.kltopt.max_iterations,
                             tc->trackopt.kltopt.min_determinant,
                             tc->trackopt.kltopt.min_displacement,
                             tc->trackopt.kltopt.max_residue,
                             tc->trackopt.kltopt.lighting_insensitive);

            if (val==KLT_SMALL_DET||val==KLT_OOB) {
                break;
            }
        }
        /* Record feature */
        if (val==KLT_OOB) {
            featurelist->feature[indx]->x=-1.0f;
            featurelist->feature[indx]->y=-1.0f;
            featurelist->feature[indx]->status=KLT_OOB;
            if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
            if (featurelist->feature[indx]->aff_img_gradx) {
                free_floatimg(featurelist->feature[indx]->aff_img_gradx);
            }
            if (featurelist->feature[indx]->aff_img_grady) {
                free_floatimg(featurelist->feature[indx]->aff_img_grady);
            }
            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;

        }
        else if (outOfBounds(xlocout,ylocout,ncols,nrows,tc->trackopt.kltopt.borderx,tc->trackopt.kltopt.bordery)) {
            featurelist->feature[indx]->x=-1.0f;
            featurelist->feature[indx]->y=-1.0f;
            featurelist->feature[indx]->status=KLT_OOB;
            if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
            if (featurelist->feature[indx]->aff_img_gradx) {
                free_floatimg(featurelist->feature[indx]->aff_img_gradx);
            }
            if (featurelist->feature[indx]->aff_img_grady) {
                free_floatimg(featurelist->feature[indx]->aff_img_grady);
            }
            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;
        }
        else if (val==KLT_SMALL_DET)  {
            featurelist->feature[indx]->x=-1.0f;
            featurelist->feature[indx]->y=-1.0f;
            featurelist->feature[indx]->status=KLT_SMALL_DET;
            if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
            if (featurelist->feature[indx]->aff_img_gradx) {
                free_floatimg(featurelist->feature[indx]->aff_img_gradx);
            }
            if (featurelist->feature[indx]->aff_img_grady) {
                free_floatimg(featurelist->feature[indx]->aff_img_grady);
            }
            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;
        }
        else if (val==KLT_LARGE_RESIDUE) {
            featurelist->feature[indx]->x=-1.0f;
            featurelist->feature[indx]->y=-1.0f;
            featurelist->feature[indx]->status=KLT_LARGE_RESIDUE;

            if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
            if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
            if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);

            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;
        }
        else if (val==KLT_MAX_ITERATIONS) {
            featurelist->feature[indx]->x=-1.0f;
            featurelist->feature[indx]->y=-1.0f;
            featurelist->feature[indx]->status=KLT_MAX_ITERATIONS;
            if (featurelist->feature[indx]->aff_img) free_floatimg(featurelist->feature[indx]->aff_img);
            if (featurelist->feature[indx]->aff_img_gradx) free_floatimg(featurelist->feature[indx]->aff_img_gradx);
            if (featurelist->feature[indx]->aff_img_grady) free_floatimg(featurelist->feature[indx]->aff_img_grady);
            featurelist->feature[indx]->aff_img=NULL;
            featurelist->feature[indx]->aff_img_gradx=NULL;
            featurelist->feature[indx]->aff_img_grady=NULL;
        }
        else {
            featurelist->feature[indx]->x=xlocout;
            featurelist->feature[indx]->y=ylocout;
            featurelist->feature[indx]->status=KLT_TRACKED;
            featurelist->feature[indx]->count++;

            /*for affine mapping*/
            if (tc->trackopt.kltopt.affineConsistencyCheck>=0&&val==KLT_TRACKED) {
                /* add border for interpolation */
                int border=2;

                if (!featurelist->feature[indx]->aff_img){

                    /* save image and gradient for each feature at finest resolution after first successful track */
                    featurelist->feature[indx]->aff_img=create_floatimg((tc->trackopt.kltopt.affine_window_width+border),(tc->trackopt.kltopt.affine_window_height+border));
                    featurelist->feature[indx]->aff_img_gradx=create_floatimg((tc->trackopt.kltopt.affine_window_width+border),(tc->trackopt.kltopt.affine_window_height+border));
                    featurelist->feature[indx]->aff_img_grady=create_floatimg((tc->trackopt.kltopt.affine_window_width+border),(tc->trackopt.kltopt.affine_window_height+border));

                    am_getSubFloatImage(pyramid1->img[0],xloc,yloc,featurelist->feature[indx]->aff_img);
                    am_getSubFloatImage(pyramid1_gradx->img[0],xloc,yloc,featurelist->feature[indx]->aff_img_gradx);
                    am_getSubFloatImage(pyramid1_grady->img[0],xloc,yloc,featurelist->feature[indx]->aff_img_grady);

                    featurelist->feature[indx]->aff_x=xloc-(int)xloc+(tc->trackopt.kltopt.affine_window_width+border)/2;
                    featurelist->feature[indx]->aff_y=yloc-(int)yloc+(tc->trackopt.kltopt.affine_window_height+border)/2;;
                }
                else {
                    /* affine tracking */
                    val=am_trackFeatureAffine(featurelist->feature[indx]->aff_x,featurelist->feature[indx]->aff_y,
                                              &xlocout, &ylocout,
                                              featurelist->feature[indx]->aff_img,
                                              featurelist->feature[indx]->aff_img_gradx,
                                              featurelist->feature[indx]->aff_img_grady,
                                              pyramid2->img[0],
                                              pyramid2_gradx->img[0], pyramid2_grady->img[0],
                                              tc->trackopt.kltopt.affine_window_width,tc->trackopt.kltopt.affine_window_height,
                                              tc->trackopt.kltopt.step_factor,
                                              tc->trackopt.kltopt.affine_max_iterations,
                                              tc->trackopt.kltopt.min_determinant,
                                              tc->trackopt.kltopt.min_displacement,
                                              tc->trackopt.kltopt.affine_min_displacement,
                                              tc->trackopt.kltopt.affine_max_residue,
                                              tc->trackopt.kltopt.lighting_insensitive,
                                              tc->trackopt.kltopt.affineConsistencyCheck,
                                              tc->trackopt.kltopt.affine_max_displacement_differ,
                                              &featurelist->feature[indx]->aff_Axx,
                                              &featurelist->feature[indx]->aff_Ayx,
                                              &featurelist->feature[indx]->aff_Axy,
                                              &featurelist->feature[indx]->aff_Ayy
                    );
                    featurelist->feature[indx]->status=val;
                    if (val!=KLT_TRACKED) {
                        featurelist->feature[indx]->x=-1.0f;
                        featurelist->feature[indx]->y=-1.0f;
                        featurelist->feature[indx]->aff_x=-1.0f;
                        featurelist->feature[indx]->aff_y=-1.0f;

                        /* free image and gradient for lost feature */
                        free_floatimg(featurelist->feature[indx]->aff_img);
                        free_floatimg(featurelist->feature[indx]->aff_img_gradx);
                        free_floatimg(featurelist->feature[indx]->aff_img_grady);
                        featurelist->feature[indx]->aff_img=NULL;
                        featurelist->feature[indx]->aff_img_gradx=NULL;
                        featurelist->feature[indx]->aff_img_grady=NULL;
                    }
                    else {
                        featurelist->feature[indx]->x=xlocout;
                        featurelist->feature[indx]->y=ylocout;
                    }
                }
            }
        }
    }
    if (tc->trackopt.kltopt.sequentialMode) {
        tc->pyramid_last=pyramid2;
        tc->pyramid_last_gradx=pyramid2_gradx;
        tc->pyramid_last_grady=pyramid2_grady;
    }
    else {
        free_pyramid(pyramid2);
        free_pyramid(pyramid2_gradx);
        free_pyramid(pyramid2_grady);
    }
    /* Free memory */
    free_floatimg(tmpimg);
    if (floatimg1_created) free_floatimg(floatimg1);

    free_floatimg(floatimg2);
    free_pyramid(pyramid1);
    free_pyramid(pyramid1_gradx);
    free_pyramid(pyramid1_grady);

    trace(3,"\n\t%d features successfully tracked.\n",
          count_good_features(featurelist));
    return;
}
/* initial klt track---------------------------------------------------------*/
extern void initklt()
{
    trace(3,"initklt:\n");
    kltrack=klt_create_tracking_context();
}
/* free klt track------------------------------------------------------------*/
extern void freeklt()
{
    trace(3,"freeklt:\n");
    klt_free_track_context(kltrack);
}
/* get feature point klt track status-----------------------------------------
 * args:    match_point_t *matchp  IO  matched feature point
 *          img_t* pimg,cimg       I   precious/current image data
 *          voopt_t *opt           I   visual odometry options
 * return: status (1: ok, 0: fail)
 * ---------------------------------------------------------------------------*/
extern int kltstatus(match_point_t *matchp,const img_t *pimg,const img_t *cimg,
                     const voopt_t *opt)
{
    pfeaturelist_t fl;
    unsigned char *img1,*img2;
    int size=sizeof(unsigned char)*cimg->h*cimg->w;

    trace(3,"kltstatus:\n");

    if (kltrack==NULL) initklt();

    /* initial workspace */
    fl=klt_create_featurelist(1);
    fl->feature[0]->x=matchp->up; fl->feature[0]->y=matchp->vp;

    img1=(unsigned char*)malloc(size);
    img2=(unsigned char*)malloc(size);

    memcpy(img1,pimg->data,size);
    memcpy(img2,cimg->data,size);
#if OUTPUT_PPM
    /* just for debug */
    writeFeatureListToPPM(fl,img1,pimg->w,pimg->h,"./klt_feat1.ppm");
#endif
    /* klt track status */
    klt(kltrack,img1,img2,cimg->w,cimg->h,fl);

#if OUTPUT_PPM
    /* just for debug */
    writeFeatureListToPPM(fl,img2,cimg->w,cimg->h,"./klt_feat2.ppm");
#endif
    matchp->kltstat=fl->feature[0]->status;

    /* check feature point track */
    if (SQRT(SQR(fl->feature[0]->x-matchp->uc)+SQR(fl->feature[0]->y-matchp->vc))
        >KLT_THRES) {
        matchp->kltstat=KLT_LARGE_RESIDUE;

        trace(2,"klt track statuc check fail\n");
        goto exit;
    }
exit:
    klt_free_featurelist(fl);
    free(img1); free(img2);
    return 1;
}
