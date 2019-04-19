/*-----------------------------------------------------------------------------
 * ins-img-filter.cc : feature match functions
 *
 * reference :
 *    [1] Kitt B , Geiger A , Lategahn H . Visual odometry based on stereo image
 *    sequences with RANSAC-based outlier rejection scheme[C]// Intelligent
 *    Vehicles Symposium. IEEE, 2010.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/11/18 1.0 new
 *----------------------------------------------------------------------------*/
#include <carvig.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <triangle.h>

/* constants ----------------------------------------------------------------*/
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

#define TRACR_FEAT_POINTS    0       /* output feature points for debugs */
#define MATCH_RATIO          5.0     /* ratio of min-cost to second-min-cost */
#define USE_NEW_BUCKET       1       /* use new bucket feature exctract method */
#define CHK_MATCH_ROI        1       /* check match ROI options */
#define FAST_CORNERS         0       /* FAST corner detect */
#define SCORE_THRESHOLD      20.0    /* Shi-Thoma score threshold */

/* type definitions ----------------------------------------------------------*/
typedef struct maximum {             /* structure for storing interest points */
    int32_t u;                       /* u-coordinate */
    int32_t v;                       /* v-coordinate */
    int32_t val;                     /* value */
    int32_t c;                       /* class */
    int32_t d1,d2,d3,d4,d5,d6,d7,d8; /* descriptor */
} maximum_t;

typedef struct maximum_set {         /* store all interest points */
    int n,nmax;                      /* number and max number of interest points */
    maximum_t *data;                 /* interest points data */
} maximum_set_t;

typedef struct range {               /* u/v ranges for matching stage 0-3 */
    float u_min[4];
    float u_max[4];
    float v_min[4];
    float v_max[4];
} range_t;

typedef struct ranges {             /* store all range data */
    int n,nmax;                     /* number and max number of range data */
    range_t *data;                  /* range data */
} ranges_t;

typedef struct delta {
    float val[8];
} delta_t;

typedef struct deltas {
    int n,nmax;
    delta_t *data;
} deltas_t;

typedef struct ivec {                /* int's vector */
    int n,nmax;                      /* number and max numbers of vector */
    int *data;                       /* vector data */
} ivec_t;

typedef struct area {
    int rect[4][2];                  /* rectangle area for detect feature points */
} area_t;

typedef struct hashtable {           /* hash table type */
    int last_idx;                    /* index of track feature for matching */
    UT_hash_handle hh;               /* makes this structure hashable */
} hashtable_t;

/* global variables-----------------------------------------------------------*/
static int32_t dims_c[3]={0},dims_p[3]={0}; /* dims for current processing image data */
static ranges_t mranges={0};                /* store all range data for matching feature points */
static int32_t  margin;

static int32_t *m1p1,*m1c1;                 /* image data ring buffer */
static int32_t *m1p2,*m1c2;
static int32_t n1p1,n1c1;
static int32_t n1p2,n1c2;
static uint8_t *I1p,*I1c;
static uint8_t *I1p_du,*I1c_du;
static uint8_t *I1p_dv,*I1c_dv;
static uint8_t *I1p_du_full,*I1c_du_full;
static uint8_t *I1p_dv_full,*I1c_dv_full;
static int greplace=0;                      /* replace flag (0: no-replace,1: replace) */
static hashtable_t *hash=NULL;              /* hash table storing all track feature index in `struct:track_t' */

/* declaration add data function----------------------------------------------*/
add_data_func_delc(ivec,int)
add_data_func_delc(ranges,range)
add_data_func_delc(match_set,match_point)
add_data_func_delc(maximum_set,maximum)
add_data_func_delc(deltas,delta)

/* add element to hash table---------------------------------------------------*/
static void hash_add(hashtable_t **ht,int last_idx)
{
    struct hashtable *s=NULL;
    HASH_FIND_INT(*ht,&last_idx,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(struct hashtable *)malloc(sizeof(struct hashtable));
        s->last_idx=last_idx;
        HASH_ADD_INT(*ht,last_idx,s);
    }
    else {
        trace(2,"hash element exist\n");
    }
}
/* delete hash table-----------------------------------------------------------*/
static void hash_delete(hashtable_t **ht)
{
    struct hashtable *current,*tmp;

    HASH_ITER(hh,*ht,current,tmp) {
        HASH_DEL(*ht,current);  /* delete; users advances to next */
        free(current);          /* optional- if you want to free  */
    }
    *ht=NULL; /* delete */
}
/* counts of elements in hash table-------------------------------------------*/
static int hash_counts(const struct hashtable *ht)
{
    return HASH_COUNT(ht);
}
/* find track data from hash table--------------------------------------------*/
static struct hashtable *hash_find(hashtable_t **ht,const int last_idx)
{
    struct hashtable *s=NULL;

    HASH_FIND_INT(*ht,&last_idx,s);  /* s: output pointer */
    return s;
}
/* output feature points------------------------------------------------------*/
static void trace_match_points(const match_set_t *match)
{
    register int i;
    for (i=0;i<match->n;i++) {
        trace(3,"ip=%3d  ic=%3d  (up,vp)=(%8.3lf,%8.3lf),(uc,vc)=(%8.3lf,%8.3lf)\n",
              match->data[i].ip,match->data[i].ic,
              match->data[i].up,
              match->data[i].vp,
              match->data[i].uc,match->data[i].vc);
    }
}
/* initial a int's vector-----------------------------------------------------*/
static int init_ivec(ivec_t **ivec,const int n)
{
    register int i;
    *ivec=(ivec_t*)malloc(sizeof(ivec_t)*n);
    if (*ivec) {
        (*ivec)->n=0; (*ivec)->nmax=n;
        for (i=0;i<n;i++) {
            (*ivec)[i].n=(*ivec)[i].nmax=0; (*ivec)[i].data=NULL;
        }
        return 1;
    }
    return 0;
}
/* free int's vector----------------------------------------------------------*/
static void free_ivec(ivec_t *ivec)
{
    if (ivec->data) {
        free(ivec->data); ivec->data=NULL;
    }
    ivec->n=ivec->nmax=0;
}
/* create position/class bin index vectors------------------------------------*/
static void create_index(const matchopt_t *opt,int32_t* m,int32_t n,ivec_t *k,
                         const int32_t ubn,
                         const int32_t vbn)
{
    /* descriptor step size */
    register int step=sizeof(maximum_t)/sizeof(int32_t),i;

    /* for all points do */
    for (i=0;i<n;i++) {

        /* extract coordinates and class */
        register int u=*(m+step*i+0); /* u-coordinate */
        register int v=*(m+step*i+1); /* v-coordinate */
        register int c=*(m+step*i+3); /* class */

        /* compute row and column of bin to which this observation belongs */
        register int ub=MIN((int)floor((float)u/(float)opt->match_binsize),ubn-1);
        register int vb=MIN((int)floor((float)v/(float)opt->match_binsize),vbn-1);

        /* save index */
        if (!add_ivec_int(&k[(c*vbn+vb)*ubn+ub],&i)) {
            continue;
        }
    }
}
/* find matched feature points------------------------------------------------*/
static int find_match(int32_t* m1,const int i1,int32_t* m2,
                      const int step,ivec_t *k2,
                      const int u_bin_num, const int v_bin_num, const int stat_bin,
                      int *min_ind,int stage,bool use_prior,
                      const matchopt_t *opt,int *flag)
{
    /* init and load image coordinates + feature */
    *min_ind        =0;
    register double  min_cost=1E9,scdmin_cost=1E9;
    int32_t u1      =*(m1+step*i1+0);
    int32_t v1      =*(m1+step*i1+1);
    int32_t c       =*(m1+step*i1+3);
    __m128i xmm1    =_mm_load_si128((__m128i*)(m1+step*i1+4));
    __m128i xmm2    =_mm_load_si128((__m128i*)(m1+step*i1+8));
    __m128i xmm3,xmm4;

    register double cost;
    register float u_min,u_max,v_min,v_max;
    register int u_bin,v_bin,i,u2,v2,k2_ind;

    /* restrict search range with prior */
    if (use_prior) {
        u_min=u1+mranges.data[stat_bin].u_min[stage];
        u_max=u1+mranges.data[stat_bin].u_max[stage];
        v_min=v1+mranges.data[stat_bin].v_min[stage];
        v_max=v1+mranges.data[stat_bin].v_max[stage];
    }
    else {
        u_min=u1-opt->match_radius;
        u_max=u1+opt->match_radius;
        v_min=v1-opt->match_radius;
        v_max=v1+opt->match_radius;
    }
    /* bins of interest */
    register int u_bin_min=MIN(MAX((int)floor(u_min/(float)opt->match_binsize),0),u_bin_num-1);
    register int u_bin_max=MIN(MAX((int)floor(u_max/(float)opt->match_binsize),0),u_bin_num-1);
    register int v_bin_min=MIN(MAX((int)floor(v_min/(float)opt->match_binsize),0),v_bin_num-1);
    register int v_bin_max=MIN(MAX((int)floor(v_max/(float)opt->match_binsize),0),v_bin_num-1);

    /* for all bins of interest do */
    for (u_bin=u_bin_min;u_bin<=u_bin_max;u_bin++) {
        for (v_bin=v_bin_min;v_bin<=v_bin_max;v_bin++) {

            k2_ind=(c*v_bin_num+v_bin)*u_bin_num+u_bin;
            for (i=0;i<k2[k2_ind].n;i++) {

                u2=*(m2+step*(k2[k2_ind].data[i])+0);
                v2=*(m2+step*(k2[k2_ind].data[i])+1);

                if (u2>=u_min&&u2<=u_max&&v2>=v_min&&v2<=v_max) {
                    xmm3=_mm_load_si128((__m128i*)(m2+step*(k2[k2_ind].data[i])+4));
                    xmm4=_mm_load_si128((__m128i*)(m2+step*(k2[k2_ind].data[i])+8));

                    xmm3=_mm_sad_epu8 (xmm1,xmm3);
                    xmm4=_mm_sad_epu8 (xmm2,xmm4);
                    xmm4=_mm_add_epi16(xmm3,xmm4);
                    cost=_mm_extract_epi16(xmm4,0)+_mm_extract_epi16(xmm4,4);

                    if (cost<min_cost) {
                        *min_ind=k2[k2_ind].data[i];
                        scdmin_cost=min_cost; min_cost=cost;
                    }
                }
            }
        }
    }
    /* ratio test */
    if (scdmin_cost/min_cost<MATCH_RATIO) *flag+=1;

    /* return index of feature point */
    return *min_ind;
}
/* computes the address offset for coordinates u,v of an image of given width -*/
static inline int get_address_offset_image(int u,int v,int w)
{
    return v*w+u;
}
/* internal function for matching feature points------------------------------*/
static int match_internal(const matchopt_t *opt,match_set_t *mset,int32_t *mp,
                          int32_t *mc,int np,int nc,bool use_prior,int type)
{
    register int step,ub,vb,b,*M,uc,vc,un,vn,nb,flag;
    ivec_t *kp,*kc;
    register int i,ip,i1c2,up,vp,ind;
    match_point p;

    trace(3,"match_features:\n");

    /* descriptor step size (number of int32_t elements in struct) */
    step=sizeof(maximum_t)/sizeof(int32_t);

    /* compute number of bins */
    ub=(int)ceil((float)dims_c[0]/(float)opt->match_binsize);
    vb=(int)ceil((float)dims_c[1]/(float)opt->match_binsize);
    b =4*vb*ub;

    /* create position/class bin index vectors */
    if (!init_ivec(&kp,b)||!init_ivec(&kc,b)) {
        trace(2,"i-vector malloc fail\n");
        return 0;
    }
    M=(int*)calloc(dims_c[0]*dims_c[1],sizeof(int32_t));

    create_index(opt,mp,np,kp,ub,vb);
    create_index(opt,mc,nc,kc,ub,vb);

    /* for all points do */
    for (i=0;i<nc;i++) {

        flag=0;

        /* coordinates in current image */
        uc=*(mc+step*i+0);
        vc=*(mc+step*i+1);

        un=MIN((int32_t)floor((float)uc/(float)opt->match_binsize),ub-1);
        vn=MIN((int32_t)floor((float)vc/(float)opt->match_binsize),vb-1);
        nb=vn*ub+un;

        /* match forward/backward */
        find_match(mc, i,mp,step,kp,ub,vb,nb,&ip  ,0,use_prior,opt,&flag);
        find_match(mp,ip,mc,step,kc,ub,vb,nb,&i1c2,1,use_prior,opt,&flag);

        /* circle closure success */
        if (i1c2!=i) continue;

        /* ratio test */
        if (flag>=1) continue;

        /* extract coordinates */
        up=*(mp+step*ip+0);
        vp=*(mp+step*ip+1);

        /* add match if this pixel isn't matched yet */
        if (*(M+(ind=get_address_offset_image(uc,vc,dims_c[0])))==0) {

            /* match ok */
            p.up=up; p.vp=vp; p.ip=ip;
            p.uc=uc; p.vc=vc; p.ic=i;
            add_match_set_match_point(mset,&p);

            /* set match flag of this feature point */
            *(M+ind)=1;
        }
    }
    for (i=0;i<b;i++) {
        free_ivec(&kp[i]); free_ivec(&kc[i]);
    }
    free(M);
    return mset->n;
}
static void get_half_resolution_dims(const int32_t *dims,int32_t *dims_half)
{
    dims_half[0]=dims[0]/2;
    dims_half[1]=dims[1]/2;
    dims_half[2]=dims_half[0]+15-(dims_half[0]-1)%16;
}
/* create half resolution image-----------------------------------------------*/
static uint8_t* create_halt_resimg(uint8_t *I,const int32_t* dims)
{
    int32_t dims_half[3],v,u;
    get_half_resolution_dims(dims,dims_half);
    unsigned char* I_half=(uint8_t*)_mm_malloc(dims_half[2]*dims_half[1]*sizeof(unsigned char),16);

    for (v=0;v<dims_half[1];v++)
        for (u=0;u<dims_half[0];u++)
            I_half[v*dims_half[2]+u]=(unsigned char)(((int32_t)I[(v*2+0)*dims[2]+u*2+0]+
                                                      (int32_t)I[(v*2+0)*dims[2]+u*2+1]+
                                                      (int32_t)I[(v*2+1)*dims[2]+u*2+0]+
                                                      (int32_t)I[(v*2+1)*dims[2]+u*2+1])/4);
    return I_half;
}
/* in area?-------------------------------------------------------------------*/
static int inarea(const area_t *area,int u,int v)
{
    return u<area->rect[3][0]&&u>area->rect[0][0]&&
           v>area->rect[0][1]&&v<area->rect[3][1];
}
/* free maximum set-----------------------------------------------------------*/
static void free_maximun_set(maximum_set_t *maxset)
{
    if (maxset->data) {
        free(maxset->data); maxset->data=NULL;
    }
    maxset->n=maxset->nmax=0;
}
/* using FAST corner detect corners-------------------------------------------*/
static void fastfeatdetect(const unsigned char *I,const int img_w,const int img_h,
                           short barrier,const matchopt_t *opt,
                           maximum_set *maxima)
{
    int *uv=imat(img_w*img_h*3,1),num=0,i;
    maximum_t pt;

    trace(3,"fastfeatdetect:\n");

    /* FAST corners detect */
    num=fastfeats(I,img_w,img_h,barrier,opt,uv,&num);
    if (num<=0) {
        trace(2,"FAST corners detect fail\n");
        return;
    }
    /* convert to maximum-set */
    free_maximun_set(maxima);
    for (i=0;i<num;i++) {

        pt.u=uv[3*i+0]; pt.v=uv[3*i+1];
        pt.val=uv[3*i+2];
        pt.c=1;
        add_maximum_set_maximum(maxima,&pt);
    }
    free(uv);
}
/* compute feature point descriptor-------------------------------------------*/
static void compute_descriptor(const uint8_t *I_du, const uint8_t *I_dv,
                               const int32_t &bpl, const int32_t &u, const int32_t &v,
                               uint8_t *desc_addr)
{
    /* get address indices */
    int32_t addr_m1=get_address_offset_image(u,v-1,bpl);
    int32_t addr_m3=addr_m1-2*bpl;
    int32_t addr_m5=addr_m3-2*bpl;
    int32_t addr_p1=addr_m1+2*bpl;
    int32_t addr_p3=addr_p1+2*bpl;
    int32_t addr_p5=addr_p3+2*bpl;

    /* compute descriptor */
    uint32_t k=0;
    desc_addr[k++]=I_du[addr_m1-3]; desc_addr[k++]=I_dv[addr_m1-3];
    desc_addr[k++]=I_du[addr_p1-3]; desc_addr[k++]=I_dv[addr_p1-3];
    desc_addr[k++]=I_du[addr_m1-1]; desc_addr[k++]=I_dv[addr_m1-1];
    desc_addr[k++]=I_du[addr_p1-1]; desc_addr[k++]=I_dv[addr_p1-1];
    desc_addr[k++]=I_du[addr_m1+3]; desc_addr[k++]=I_dv[addr_m1+3];
    desc_addr[k++]=I_du[addr_p1+3]; desc_addr[k++]=I_dv[addr_p1+3];
    desc_addr[k++]=I_du[addr_m1+1]; desc_addr[k++]=I_dv[addr_m1+1];
    desc_addr[k++]=I_du[addr_p1+1]; desc_addr[k++]=I_dv[addr_p1+1];
    desc_addr[k++]=I_du[addr_m5-1]; desc_addr[k++]=I_dv[addr_m5-1];
    desc_addr[k++]=I_du[addr_p5-1]; desc_addr[k++]=I_dv[addr_p5-1];
    desc_addr[k++]=I_du[addr_m5+1]; desc_addr[k++]=I_dv[addr_m5+1];
    desc_addr[k++]=I_du[addr_p5+1]; desc_addr[k++]=I_dv[addr_p5+1];
    desc_addr[k++]=I_du[addr_m3-5]; desc_addr[k++]=I_dv[addr_m3-5];
    desc_addr[k++]=I_du[addr_p3-5]; desc_addr[k++]=I_dv[addr_p3-5];
    desc_addr[k++]=I_du[addr_m3+5]; desc_addr[k++]=I_dv[addr_m3+5];
    desc_addr[k++]=I_du[addr_p3+5];
    desc_addr[k++]=I_dv[addr_p3+5];
}
/* compute all feature point descriptor---------------------------------------*/
static void compute_descriptors(uint8_t *I_du, uint8_t *I_dv, const int32_t bpl,
                                maximum_set_t *maxset)
{
    /* loop variables */
    register int u,v,i;
    unsigned char *desc_addr;

    /* for all maxima do */
    for (i=0;i<maxset->n;i++) {
        u=maxset->data[i].u;
        v=maxset->data[i].v;
        desc_addr=(uint8_t*)(&(maxset->data[i].d1));
        compute_descriptor(I_du,I_dv,bpl,u,v,desc_addr);
    }
}
/* compute FAST corners-------------------------------------------------------*/
static void compute_fastfeats(uint8_t *I,const int32_t* dims,
                              int32_t* &max, int32_t &num,
                              uint8_t* &I_du,uint8_t* &I_dv,
                              uint8_t* &I_du_full,uint8_t* &I_dv_full,
                              const matchopt_t *opt)
{
    register int16_t *I_f1;
    register int16_t *I_f2;
    int i,dims_matching[3];

    trace(3,"compute_fastfeats:\n");

    memcpy(dims_matching,dims,3*sizeof(int32_t));

    /* allocate memory for sobel images and filter images */
    if (!opt->half_res) {
        I_du=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);
        I_dv=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);
        I_f1=(int16_t*)_mm_malloc(dims[2]*dims[1]*sizeof(int16_t),16);
        I_f2=(int16_t*)_mm_malloc(dims[2]*dims[1]*sizeof(int16_t),16);

        sobel5x5(I,I_du,I_dv,dims[2],dims[1]);
        blob5x5(I,I_f1,dims[2],dims[1]);
        checkerboard5x5(I,I_f2,dims[2],dims[1]);
    }
    else {
        unsigned char* I_matching=create_halt_resimg(I,dims);
        get_half_resolution_dims(dims,dims_matching);
        I_du=(unsigned char*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(unsigned char*),16);
        I_dv=(unsigned char*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(unsigned char*),16);

        I_f1=(int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);
        I_f2=(int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);

        I_du_full=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);
        I_dv_full=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);

        sobel5x5(I_matching,I_du,I_dv,dims_matching[2],dims_matching[1]);
        sobel5x5(I,I_du_full,I_dv_full,dims[2],dims[1]);
        blob5x5(I_matching,I_f1,dims_matching[2],dims_matching[1]);
        checkerboard5x5(I_matching,I_f2,dims_matching[2],dims_matching[1]);
        free(I_matching);
    }
    /* extract dense maxima via non-maximum suppression */
    maximum_set mset={0};
    fastfeatdetect(I,dims[0],dims[1],dims[0],opt,&mset);
    compute_descriptors(I_du,I_dv,dims[2],&mset);

    /* release filter images */
    free(I_f1);
    free(I_f2);

    /* get number of interest points and init maxima pointer to NULL */
    num=mset.n;
    max=NULL;

    /* return sparse maxima as 16-bytes aligned memory */
    if (num!=0) {
        max=(int32_t*)_mm_malloc(sizeof(maximum)*num,16);
        int32_t k=0;
        for (i=0;i<mset.n;i++) {
            *(max+k++)=mset.data[i].u; *(max+k++)=mset.data[i].v; *(max+k++)=mset.data[i].val;
            *(max+k++)=mset.data[i].c;
            *(max+k++)=mset.data[i].d1;
            *(max+k++)=mset.data[i].d2;
            *(max+k++)=mset.data[i].d3;
            *(max+k++)=mset.data[i].d4;
            *(max+k++)=mset.data[i].d5;
            *(max+k++)=mset.data[i].d6;
            *(max+k++)=mset.data[i].d7;
            *(max+k++)=mset.data[i].d8;
        }
    }
    free_maximun_set(&mset);
}
/* Alexander Neubeck and Luc Van Gool: Efficient Non-Maximum Suppression,-----
 * ICPR'06, algorithm 4
 * ---------------------------------------------------------------------------*/
static void nonMaximumSuppression(int16_t* I_f1,int16_t* I_f2,const int32_t* dims,
                                  maximum_set *maxima,int nms_n,
                                  const matchopt_t *opt,const area_t *area)
{
    /* extract parameters */
    register int width =dims[0];
    register int height=dims[1];
    register int bpl   =dims[2];
    register int n     =nms_n;
    register int tau   =opt->nms_tau;

    maximum_t max;

    /* loop variables */
    register int32_t f1mini,f1minj,f1maxi,f1maxj,f2mini,f2minj,f2maxi,f2maxj;
    register int32_t f1minval,f1maxval,f2minval,f2maxval,currval;
    register int32_t addr;

    for (int i=n+margin;i<width-n-margin;i+=n+1) {
        for (int j=n+margin;j<height-n-margin;j+=n+1) {

            if (area&&!inarea(area,i,j)) continue;

            f1mini=i; f1minj=j; f1maxi=i; f1maxj=j;
            f2mini=i; f2minj=j; f2maxi=i; f2maxj=j;

            addr    =get_address_offset_image(i,j,bpl);
            f1minval=*(I_f1+addr);
            f1maxval=f1minval;
            f2minval=*(I_f2+addr);
            f2maxval=f2minval;

            for (int i2=i; i2<=i+n; i2++) {
                for (int j2=j; j2<=j+n; j2++) {
                    addr   =get_address_offset_image(i2,j2,bpl);
                    currval=*(I_f1+addr);
                    if (currval<f1minval) {
                        f1mini  =i2;
                        f1minj  =j2;
                        f1minval=currval;
                    }
                    else if (currval>f1maxval) {
                        f1maxi  =i2;
                        f1maxj  =j2;
                        f1maxval=currval;
                    }
                    currval=*(I_f2+addr);
                    if (currval<f2minval) {
                        f2mini  =i2;
                        f2minj  =j2;
                        f2minval=currval;
                    }
                    else if (currval>f2maxval) {
                        f2maxi  =i2;
                        f2maxj  =j2;
                        f2maxval=currval;
                    }
                }
            }
            /* f1 minimum */
            for (int i2=f1mini-n;i2<=MIN(f1mini+n,width-1-margin);i2++) {
                for (int j2=f1minj-n;j2<=MIN(f1minj+n,height-1-margin);j2++) {

                    currval=*(I_f1+get_address_offset_image(i2,j2,bpl));
                    if (currval<f1minval&&(i2<i||i2>i+n||j2<j||j2>j+n)) {
                        goto failed_f1min;
                    }
                }
            }
            if (f1minval<=-tau) {
                max.u=f1mini;
                max.v=f1minj;
                max.val=f1minval;
                max.c=0;
                add_maximum_set_maximum(maxima,&max);
            }
failed_f1min:;
            /* f1 maximum */
            for (int i2=f1maxi-n;i2<=MIN(f1maxi+n,width-1-margin);i2++) {
                for (int j2=f1maxj-n;j2<=MIN(f1maxj+n,height-1-margin);j2++) {

                    currval=*(I_f1+get_address_offset_image(i2,j2,bpl));
                    if (currval>f1maxval&&(i2<i||i2>i+n||j2<j||j2>j+n)) {
                        goto failed_f1max;
                    }
                }
            }
            if (f1maxval>=tau) {
                max.u=f1maxi;
                max.v=f1maxj;
                max.val=f1maxval;
                max.c=1;
                add_maximum_set_maximum(maxima,&max);
            }
failed_f1max:;
            /* f2 minimum */
            for (int i2=f2mini-n;i2<=MIN(f2mini+n,width-1-margin);i2++) {
                for (int j2=f2minj-n;j2<=MIN(f2minj+n,height-1-margin);j2++) {

                    currval=*(I_f2+get_address_offset_image(i2,j2,bpl));
                    if (currval<f2minval&&(i2<i||i2>i+n||j2<j||j2>j+n)) {
                        goto failed_f2min;
                    }
                }
            }
            if (f2minval<=-tau) {
                max.u=f2mini;
                max.v=f2minj;
                max.val=f2minval;
                max.c=2;
                add_maximum_set_maximum(maxima,&max);
            }
failed_f2min:;
            /* f2 maximum */
            for (int i2=f2maxi-n;i2<=MIN(f2maxi+n,width-1-margin);i2++) {
                for (int j2=f2maxj-n;j2<=MIN(f2maxj+n,height-1-margin);j2++) {

                    currval=*(I_f2+get_address_offset_image(i2,j2,bpl));
                    if (currval>f2maxval&&(i2<i||i2>i+n||j2<j||j2>j+n))
                        goto failed_f2max;
                }
            }
            if (f2maxval>=tau) {
                max.u=f2maxi;
                max.v=f2maxj;
                max.val=f2maxval;
                max.c=3;
                add_maximum_set_maximum(maxima,&max);
            }
failed_f2max:;
        }
    }
}
/* compute shi-tomasi score---------------------------------------------------*/
static float shitomasiscore(const unsigned char *img, const int img_w,
                            const int img_h,
                            int u, int v)
{
    register float dXX=0.0;
    register float dYY=0.0;
    register float dXY=0.0;
    int halfbox_size=4;
    int box_size=2*halfbox_size;
    int box_area=box_size*box_size;
    int x_min=u-halfbox_size;
    int x_max=u+halfbox_size;
    int y_min=v-halfbox_size;
    int y_max=v+halfbox_size;

    if (x_min<1||x_max>=img_w-1||y_min<1||y_max>=img_h-1) {
        return 0.0;  /* patch is too close to the boundary */
    }
    int stride=img_w;
    for (int y=y_min;y<y_max;++y) {

        const uint8_t* ptr_left  =img+stride*y+x_min-1;
        const uint8_t* ptr_right =img+stride*y+x_min+1;
        const uint8_t* ptr_top   =img+stride*(y-1)+x_min;
        const uint8_t* ptr_bottom=img+stride*(y+1)+x_min;
        for (int x=0;x<box_size;++x,++ptr_left,++ptr_right,++ptr_top,++ptr_bottom) {
            float dx=*ptr_right -*ptr_left;
            float dy=*ptr_bottom-*ptr_top;
            dXX+=dx*dx;
            dYY+=dy*dy;
            dXY+=dx*dy;
        }
    }
    /* find and return smaller eigenvalue: */
    dXX=dXX/(2.0f*box_area);
    dYY=dYY/(2.0f*box_area);
    dXY=dXY/(2.0f*box_area);
    return 0.5f*(dXX+dYY-sqrtf((dXX+dYY)*(dXX+dYY)-4.0f*(dXX*dYY-dXY*dXY)));
}
/* compute new features for current frame-------------------------------------*/
static int compute_features(uint8_t *I,const int32_t* dims,
                            int32_t* &max1, int32_t &num1,
                            int32_t* &max2, int32_t &num2,
                            uint8_t* &I_du,uint8_t* &I_dv,
                            uint8_t* &I_du_full,uint8_t* &I_dv_full,
                            const matchopt_t *opt)
{
    trace(3,"compute_features:\n");

    register int16_t *I_f1;
    register int16_t *I_f2;

    register int32_t dims_matching[3],i,k;
    register int32_t s=1;
    float score;
    
    memcpy(dims_matching,dims,3*sizeof(int32_t));

    /* allocate memory for sobel images and filter images */
    if (!opt->half_res) {
        I_du=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);
        I_dv=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);
        I_f1=(int16_t*)_mm_malloc(dims[2]*dims[1]*sizeof(int16_t),16);
        I_f2=(int16_t*)_mm_malloc(dims[2]*dims[1]*sizeof(int16_t),16);

        sobel5x5(I,I_du,I_dv,dims[2],dims[1]);
        blob5x5(I,I_f1,dims[2],dims[1]);
        checkerboard5x5(I,I_f2,dims[2],dims[1]);
    }
    else {
        unsigned char* I_matching=create_halt_resimg(I,dims);
        get_half_resolution_dims(dims,dims_matching);
        I_du=(unsigned char*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(unsigned char*),16);
        I_dv=(unsigned char*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(unsigned char*),16);

        I_f1=(int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);
        I_f2=(int16_t*)_mm_malloc(dims_matching[2]*dims_matching[1]*sizeof(int16_t),16);

        I_du_full=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);
        I_dv_full=(unsigned char*)_mm_malloc(dims[2]*dims[1]*sizeof(unsigned char*),16);

        sobel5x5(I_matching,I_du,I_dv,dims_matching[2],dims_matching[1]);
        sobel5x5(I,I_du_full,I_dv_full,dims[2],dims[1]);
        blob5x5(I_matching,I_f1,dims_matching[2],dims_matching[1]);
        checkerboard5x5(I_matching,I_f2,dims_matching[2],dims_matching[1]);
        free(I_matching);
    }
    /* extract sparse maxima (1st pass) via non-maximum suppression */
    maximum_set mset1={0};
    if (opt->multi_stage) {
        int nms_n_sparse=opt->nms_n*3;
        if (nms_n_sparse>10) nms_n_sparse=MAX(opt->nms_n,10);

        nonMaximumSuppression(I_f1,I_f2,dims_matching,&mset1,nms_n_sparse,opt,NULL);
        compute_descriptors(I_du,I_dv,dims_matching[2],&mset1);

        /* shi-tomasi score check */
        for (k=0,i=0;i<mset1.n;i++) {
            score=shitomasiscore(I,dims[0],dims[1],mset1.data[i].u,mset1.data[i].v);
            if (score<SCORE_THRESHOLD) {
                continue;
            }
            memcpy(&mset1.data[k++],&mset1.data[i],sizeof(maximum_t));
        }
        mset1.n=k;
    }
    /* extract dense maxima (2nd pass) via non-maximum suppression */
    maximum_set mset2={0};
    nonMaximumSuppression(I_f1,I_f2,dims_matching,&mset2,opt->nms_n,opt,NULL);
    compute_descriptors(I_du,I_dv,dims_matching[2],&mset2);

    /* shi-tomasi score check */
    for (k=0,i=0;i<mset2.n;i++) {
        score=shitomasiscore(I,dims[0],dims[1],mset2.data[i].u,mset2.data[i].v);
        if (score<SCORE_THRESHOLD) {
            continue;
        }
        memcpy(&mset2.data[k++],&mset2.data[i],sizeof(maximum_t));
    }
    mset2.n=k;

    /* release filter images */
    free(I_f1);
    free(I_f2);

    /* get number of interest points and init maxima pointer to NULL */
    num1=mset1.n;
    num2=mset2.n;
    max1=NULL;
    max2=NULL;

    if (opt->half_res) s=2;

    /* return sparse maxima as 16-bytes aligned memory */
    if (num1!=0) {
        max1=(int32_t*)_mm_malloc(sizeof(maximum)*num1,16);
        int32_t k=0;
        for (i=0;i<mset1.n;i++) {
            *(max1+k++)=mset1.data[i].u*s; *(max1+k++)=mset1.data[i].v*s; *(max1+k++)=0;
            *(max1+k++)=mset1.data[i].c;
            *(max1+k++)=mset1.data[i].d1;
            *(max1+k++)=mset1.data[i].d2;
            *(max1+k++)=mset1.data[i].d3;
            *(max1+k++)=mset1.data[i].d4;
            *(max1+k++)=mset1.data[i].d5;
            *(max1+k++)=mset1.data[i].d6;
            *(max1+k++)=mset1.data[i].d7;
            *(max1+k++)=mset1.data[i].d8;
        }
    }
    /* return dense maxima as 16-bytes aligned memory */
    if (num2!=0) {
        max2=(int32_t*)_mm_malloc(sizeof(maximum)*num2,16);
        int32_t k=0;
        for (i=0;i<mset2.n;i++) {
            *(max2+k++)=mset2.data[i].u*s; *(max2+k++)=mset2.data[i].v*s; *(max2+k++)=0;
            *(max2+k++)=mset2.data[i].c;
            *(max2+k++)=mset2.data[i].d1;
            *(max2+k++)=mset2.data[i].d2;
            *(max2+k++)=mset2.data[i].d3;
            *(max2+k++)=mset2.data[i].d4;
            *(max2+k++)=mset2.data[i].d5;
            *(max2+k++)=mset2.data[i].d6;
            *(max2+k++)=mset2.data[i].d7;
            *(max2+k++)=mset2.data[i].d8;
        }
    }
    free_maximun_set(&mset1);
    free_maximun_set(&mset2);
}
/* ---------------------------------------------------------------------------
 * computes features from images and pushes them back to a ringbuffer,
 * which interally stores the features of the current and previous image pair
 * use this function for matching
 * input: I1    .......... pointers to left and right image (row-aligned),
 *                         range [0..255]
 *        dims[0,1] ...... image width and height (both images must be rectified
 *                         and of same size)
 *        dims[2] ........ bytes per line (often equals width)
 *        replace ........ if this flag is set, the current image is overwritten
 *                         with the input images, otherwise the current image is
 *                         first copied to the previous image (ring buffer
 *                         functionality,descriptors need to be computed only once)
 * --------------------------------------------------------------------------*/
static int puchbackimg(uint8_t *I1,int32_t* dims,const int replace,
                       const matchopt_t *opt)
{
    trace(3,"puchback_img:\n");

    /* image dimensions */
    register int width =dims[0];
    register int height=dims[1];
    register int bpl   =dims[2];

    /* sanity check */
    if (width<=0||height<=0||bpl<width||I1==0) {
        trace(2,"error: image dimension mismatch\n");
        return 0;
    }
    if (replace) {
        if (I1c )        _mm_free(I1c );
        if (m1c1)        _mm_free(m1c1); /* maximum points in current image */
        if (m1c2)        _mm_free(m1c2); /* maximum points in current image*/
        if (I1c_du)      _mm_free(I1c_du);
        if (I1c_dv)      _mm_free(I1c_dv);
        if (I1c_du_full) _mm_free(I1c_du_full);
        if (I1c_dv_full) _mm_free(I1c_dv_full);
    }
    else {
        if (I1p )        _mm_free(I1p );
        if (m1p1)        _mm_free(m1p1); /* maximum points in precious image */
        if (m1p2)        _mm_free(m1p2); /* maximum points in precious image */
        if (I1p_du)      _mm_free(I1p_du);
        if (I1p_dv)      _mm_free(I1p_dv);
        if (I1p_du_full) _mm_free(I1p_du_full);
        if (I1p_dv_full) _mm_free(I1p_dv_full);
        m1p1=m1c1; n1p1=n1c1;
        m1p2=m1c2; n1p2=n1c2;
        I1p        =I1c;
        I1p_du     =I1c_du;
        I1p_dv     =I1c_dv;
        I1p_du_full=I1c_du_full;
        I1p_dv_full=I1c_dv_full;
        dims_p[0]  =dims_c[0];
        dims_p[1]  =dims_c[1];
        dims_p[2]  =dims_c[2];
    }
    /* set new dims (bytes per line must be multiple of 16) */
    dims_c[0]=width;
    dims_c[1]=height;
    dims_c[2]=width+15-(width-1)%16;

    /* copy images to byte aligned memory */
    I1c=(uint8_t*)_mm_malloc(dims_c[2]*dims_c[1]*sizeof(uint8_t),16);

    if (dims_c[2]==bpl) {
        memcpy(I1c,I1,dims_c[2]*dims_c[1]*sizeof(uint8_t));
    }
    else {
        for (int32_t v=0;v<height;v++) {
            memcpy(I1c+v*dims_c[2],I1+v*bpl,dims_c[0]*sizeof(uint8_t));
        }
    }
#if FAST_CORNERS
    /* FAST corner detect */
    compute_fastfeats(I1c,dims_c,m1c2,n1c2, /* dense corners */
                      I1c_du,I1c_dv,
                      I1c_du_full,
                      I1c_dv_full,opt);
#else
    /* compute new features for current frame */
    compute_features(I1c,dims_c,m1c1,n1c1,m1c2,n1c2,
                     I1c_du,I1c_dv,
                     I1c_du_full,
                     I1c_dv_full,opt);
#endif
    return 1;
}
/* clear match set data-------------------------------------------------------*/
extern void free_match_set(match_set_t *mset)
{
    trace(3,"free_match_set:\n");

    if (mset->data) {
        free(mset->data); mset->data=NULL;
    }
    mset->n=mset->nmax=0;
}
/* free ranges data-----------------------------------------------------------*/
static void free_ranges(ranges_t *ranges)
{
    if (ranges->data) {
        free(ranges->data); ranges->data=NULL;
    } 
    ranges->n=ranges->nmax=0;
}
/* set ranges value-----------------------------------------------------------*/
static void set_delta_val(delta_t *delta,const float val)
{
    int i; for (i=0;i<8;i++) delta->val[i]=val;
}
/* remove outliers for feature points-----------------------------------------*/
static void remove_outliers(match_set_t *mset,const matchopt_t *opt,int *index,
                            int *num,int *nindex,int *nnum)
{
    float p1_flow_u,p1_flow_v,p2_flow_u,p2_flow_v,p3_flow_u,p3_flow_v;
    int p1,p2,p3,i,*num_support;
    match_set_t p_matched_copy={0};

    /* do we have enough points for outlier removal? */
    if (mset->n<=3) return;

    /* input/output structure for triangulation */
    struct triangulateio in,out;

    /* inputs */
    in.numberofpoints=mset->n;
    in.pointlist=(float*)malloc(in.numberofpoints*2*sizeof(float));

    /* create copy of p_matched, init vector with number of support points
     * and fill triangle point vector for delaunay triangulation
     * */
    num_support=imat(1,mset->n);

    for (i=0;i<mset->n;i++) {
        add_match_set_match_point(&p_matched_copy,&mset->data[i]);
        num_support[i]=0;

        in.pointlist[2*i+0]=mset->data[i].uc;
        in.pointlist[2*i+1]=mset->data[i].vc;
    }
    /* input parameters */
    in.numberofpointattributes=0;
    in.pointattributelist     =NULL;
    in.pointmarkerlist        =NULL;
    in.numberofsegments       =0;
    in.numberofholes          =0;
    in.numberofregions        =0;
    in.regionlist             =NULL;

    /* outputs */
    out.pointlist             =NULL;
    out.pointattributelist    =NULL;
    out.pointmarkerlist       =NULL;
    out.trianglelist          =NULL;
    out.triangleattributelist =NULL;
    out.neighborlist          =NULL;
    out.segmentlist           =NULL;
    out.segmentmarkerlist     =NULL;
    out.edgelist              =NULL;
    out.edgemarkerlist        =NULL;

    /* do triangulation (z=zero-based, n=neighbors, Q=quiet, B=no boundary markers)
     * attention: not using the B switch or using the n switch
     * creates a memory leak (=> use valgrind!)
     * */
    char parameters[]="zQB";
    triangulate(parameters,&in,&out,NULL);

    /* for all triangles do */
    for (i=0;i<out.numberoftriangles;i++) {

        /* extract triangle corner points */
        p1=out.trianglelist[i*3+0];
        p2=out.trianglelist[i*3+1];
        p3=out.trianglelist[i*3+2];

        /* 1. corner disparity and flow */
        p1_flow_u=p_matched_copy.data[p1].uc-p_matched_copy.data[p1].up;
        p1_flow_v=p_matched_copy.data[p1].vc-p_matched_copy.data[p1].vp;

        /* 2. corner disparity and flow */
        p2_flow_u=p_matched_copy.data[p2].uc-p_matched_copy.data[p2].up;
        p2_flow_v=p_matched_copy.data[p2].vc-p_matched_copy.data[p2].vp;

        /* 3. corner disparity and flow */
        p3_flow_u=p_matched_copy.data[p3].uc-p_matched_copy.data[p3].up;
        p3_flow_v=p_matched_copy.data[p3].vc-p_matched_copy.data[p3].vp;

        /* consistency of 1. edge */
        if (fabs(p1_flow_u-p2_flow_u)+fabs(p1_flow_v-p2_flow_v)<opt->outlier_flow_tol) {
            num_support[p1]++;
            num_support[p2]++;
        }
        /* consistency of 2. edge */
        if (fabs(p2_flow_u-p3_flow_u)+fabs(p2_flow_v-p3_flow_v)<opt->outlier_flow_tol) {
            num_support[p2]++;
            num_support[p3]++;
        }
        /* consistency of 3. edge */
        if (fabs(p1_flow_u-p3_flow_u)+fabs(p1_flow_v-p3_flow_v)<opt->outlier_flow_tol) {
            num_support[p1]++;
            num_support[p3]++;
        }
    }
    /* refill matched feature points */
    free_match_set(mset);
    for (i=0,num==NULL?0:*num=0,nnum==NULL?0:*nnum=0;i<in.numberofpoints;i++) {
        if (num_support[i]>=4) {
            add_match_set_match_point(mset,&p_matched_copy.data[i]);

            /* add to hash table */
            if (hash_find(&hash,p_matched_copy.data[i].ip)) {
                if (index) {
                    index[(*num)++]=MAX(mset->n-1,0);
                }
            }
            else {
                if (nindex) nindex[(*nnum)++]=MAX(mset->n-1,0);
            }
        }
    }
    free_match_set(&p_matched_copy);
    free(in.pointlist);
    free(out.pointlist);
    free(out.trianglelist);
    return;
}
/* free deltas----------------------------------------------------------------*/
static void free_deltas(deltas_t *deltas)
{
    if (deltas->data) {
        free(deltas->data); deltas->data=NULL;
    }
    deltas->n=deltas->nmax=0;
}
/* compute prior statistics---------------------------------------------------*/
static void compute_prior_statistics(const matchopt_t *opt,match_set_t *mset,
                                     ranges_t *rngs)
{
    /* compute number of bins */
    register int32_t u_bin_num=(int32_t)ceil((float)dims_c[0]/(float)opt->match_binsize);
    register int32_t v_bin_num=(int32_t)ceil((float)dims_c[1]/(float)opt->match_binsize);
    register int32_t bin_num  =v_bin_num*u_bin_num;
    register int32_t u_bin_min,u_bin_max,v_bin_min,v_bin_max,v_bin,u_bin;
    register int i,j;

    /* number of matching stages */
    int32_t num_stages=2;

    /* allocate bin accumulator memory */
    deltas_t *delta_accu=(deltas_t*)calloc(sizeof(deltas_t),bin_num);

    /* fill bin accumulator */
    delta delta_curr;
    for (i=0;i<mset->n;i++) {
        delta_curr.val[0]=mset->data[i].up-mset->data[i].uc;
        delta_curr.val[1]=mset->data[i].vp-mset->data[i].vc;
        delta_curr.val[2]=mset->data[i].uc-mset->data[i].up;
        delta_curr.val[3]=mset->data[i].vc-mset->data[i].vp;

        u_bin_min=MIN(MAX((int32_t)floor(mset->data[i].uc/(float)opt->match_binsize)-1,0),u_bin_num-1);
        u_bin_max=MIN(MAX((int32_t)floor(mset->data[i].uc/(float)opt->match_binsize)+1,0),u_bin_num-1);
        v_bin_min=MIN(MAX((int32_t)floor(mset->data[i].vc/(float)opt->match_binsize)-1,0),v_bin_num-1);
        v_bin_max=MIN(MAX((int32_t)floor(mset->data[i].vc/(float)opt->match_binsize)+1,0),v_bin_num-1);

        /* add to accumulator */
        for (v_bin=v_bin_min;v_bin<=v_bin_max;v_bin++) {
            for (u_bin=u_bin_min;u_bin<=u_bin_max;u_bin++) {
                add_deltas_delta(&delta_accu[v_bin*u_bin_num+u_bin],&delta_curr);
            }
        }
    }
    /* clear ranges */
    free_ranges(rngs);

    /* for all bins compute statistics */
    for (v_bin=0;v_bin<v_bin_num;v_bin++) {
        for (u_bin=0;u_bin<u_bin_num;u_bin++) {

            /* use full range in case there are no observations */
            delta delta_min;
            delta delta_max;

            set_delta_val(&delta_max,+opt->match_radius);
            set_delta_val(&delta_min,-opt->match_radius);

            /* otherwise determine delta min and delta max */
            if (delta_accu[v_bin*u_bin_num+u_bin].n>0) {

                /* init displacements 'delta' to 'infinite' */
                set_delta_val(&delta_min,+1000000);
                set_delta_val(&delta_max,-1000000);

                /* find minimum and maximum displacements */
                for (i=0;i<delta_accu[v_bin*u_bin_num+u_bin].n;i++) {

                    const delta_t *pdelta=&delta_accu[v_bin*u_bin_num+u_bin].data[i];
                    for (j=0;j<num_stages*2;j++) {
                        if (pdelta->val[j]<delta_min.val[j]) delta_min.val[j]=pdelta->val[j];
                        if (pdelta->val[j]>delta_max.val[j]) delta_max.val[j]=pdelta->val[j];
                    }
                }
            }
            /* set search range for this bin */
            range r;
            for (i=0;i<num_stages;i++) {

                /* bound minimum search range to 20x20 */
                float delta_u=delta_max.val[i*2+0]-delta_min.val[i*2+0];
                if (delta_u<20) {
                    delta_min.val[i*2+0]-=ceil((20-delta_u)/2);
                    delta_max.val[i*2+0]+=ceil((20-delta_u)/2);
                }
                float delta_v=delta_max.val[i*2+1]-delta_min.val[i*2+1];
                if (delta_v<20) {
                    delta_min.val[i*2+1]-=ceil((20-delta_v)/2);
                    delta_max.val[i*2+1]+=ceil((20-delta_v)/2);
                }
                /* set range for this bin */
                r.u_min[i]=delta_min.val[i*2+0];
                r.u_max[i]=delta_max.val[i*2+0];
                r.v_min[i]=delta_min.val[i*2+1];
                r.v_max[i]=delta_max.val[i*2+1];
            }
            add_ranges_range(rngs,&r);
        }
    }
    /* free bin accumulator memory */
    for (i=0;i<bin_num;i++) {
        free_deltas(&delta_accu[i]);
    }
    return;
}
/* compute small descriptor---------------------------------------------------*/
static void compute_small_des(const uint8_t* I_du,const uint8_t* I_dv,
                              const int32_t &bpl,const int32_t &u,const int32_t &v,
                              uint8_t *desc_addr)
{
    /* get address indices */
    int32_t addr2=get_address_offset_image(u,v,bpl);
    int32_t addr1=addr2-bpl;
    int32_t addr0=addr1-bpl;
    int32_t addr3=addr2+bpl;
    int32_t addr4=addr3+bpl;

    /* compute ELAS-descriptor */
    uint32_t k=0;
    desc_addr[k++]=I_du[addr0  ];
    desc_addr[k++]=I_du[addr1-2];
    desc_addr[k++]=I_du[addr1  ];
    desc_addr[k++]=I_du[addr1+2];
    desc_addr[k++]=I_du[addr2-1];
    desc_addr[k++]=I_du[addr2  ];
    desc_addr[k++]=I_du[addr2  ];
    desc_addr[k++]=I_du[addr2+1];
    desc_addr[k++]=I_du[addr3-2];
    desc_addr[k++]=I_du[addr3  ];
    desc_addr[k++]=I_du[addr3+2];
    desc_addr[k++]=I_du[addr4  ];
    desc_addr[k++]=I_dv[addr1  ];
    desc_addr[k++]=I_dv[addr2-1];
    desc_addr[k++]=I_dv[addr2+1];
    desc_addr[k++]=I_dv[addr3  ];
}
/* parabolic fitting----------------------------------------------------------*/
static int parabolic_fit(const uint8_t* I1_du,const uint8_t* I1_dv,const int32_t* dims1,
                         const uint8_t* I2_du,const uint8_t* I2_dv,const int32_t* dims2,
                         const float &u1,const float &v1,
                         float       &u2,float       &v2,
                         double* At,double* AtA,
                         uint8_t* desc_buffer)
{
    double *c,*b,*x,divisor,ddv,ddu;
    register int i,j,cost_curr;

    /* check if parabolic fitting is feasible (descriptors are within margin) */
    if (u2-3<margin||u2+3>dims2[0]-1-margin||
        v2-3<margin||v2+3>dims2[1]-1-margin) {
        return 0;
    }
    /* compute reference descriptor */
    __m128i xmm1,xmm2;
    compute_small_des(I1_du,I1_dv,dims1[2],u1,v1,desc_buffer);
    xmm1=_mm_load_si128((__m128i*)(desc_buffer));

    /* compute cost matrix */
    int32_t cost[49];
    for (int dv=0;dv<7;dv++) {
        for (int du=0;du<7;du++) {
            compute_small_des(I2_du,I2_dv,dims2[2],(int32_t)u2+du-3,(int32_t)v2+dv-3,desc_buffer);
            xmm2=_mm_load_si128((__m128i*)(desc_buffer));
            xmm2=_mm_sad_epu8(xmm1,xmm2);
            cost[dv*7+du]=_mm_extract_epi16(xmm2,0)+_mm_extract_epi16(xmm2,4);
        }
    }
    /* compute minimum */
    int32_t min_ind =0;
    int32_t min_cost=cost[0];
    for (i=1;i<49;i++) {
        if (cost[i]<min_cost) {
            min_ind=i; min_cost=cost[i];
        }
    }
    /* get indices */
    int32_t du=min_ind%7;
    int32_t dv=min_ind/7;

    /* if minimum is at borders => remove this match */
    if (du==0||du==6||dv==0||dv==6) {
        return 0;
    }
    /* solve least squares system */
    c=mat(9,1); b=mat(6,1); x=mat(6,1);
    for (i=-1;i<=+1;i++) {
        for (j=-1;j<=+1;j++) {
            cost_curr=cost[(dv+i)*7+(du+j)]; c[(i+1)*3+(j+1)]=cost_curr;
        }
    }
    matmul("NN",6,1,9,1.0,At,c,0.0,b);
    if (solve("N",AtA,b,6,1,x)) {
        free(b); free(c); free(x);
        return 0;
    }
    /* extract relative coordinates */
    divisor=(x[2]*x[2]-4.0*x[0]*x[1]);
    if (fabs(divisor)<1E-8||fabs(x[2])<1E-8) {
        return 0;
    }
    ddv= (2.0*x[0]*x[4]-x[2]*x[3])/divisor;
    ddu=-(x[4]+2.0*x[1]*ddv)/x[2];
    if (fabs(ddu)>=1.0||fabs(ddv)>=1.0) return 0;

    /* update target */
    u2+=(float)du-3.0+ddu;
    v2+=(float)dv-3.0+ddv;

    /* return true on success */
    free(b); free(c); free(x);
    return 1;
}
/* relocate minimum-----------------------------------------------------------*/
static void relocate_minimum(const uint8_t* I1_du,const uint8_t* I1_dv,const int32_t* dims1,
                             const uint8_t* I2_du,const uint8_t* I2_dv,const int32_t* dims2,
                             const float &u1,const float &v1,
                             float       &u2,float       &v2,
                             uint8_t* desc_buffer)
{
    /* check if parabolic fitting is feasible (descriptors are within margin) */
    if (u2-2<margin||u2+2>dims2[0]-1-margin||
        v2-2<margin||v2+2>dims2[1]-1-margin) {
        return;
    }
    /* compute reference descriptor */
    __m128i xmm1,xmm2;
    compute_small_des(I1_du,I1_dv,dims1[2],(int32_t)u1,(int32_t)v1,desc_buffer);
    xmm1=_mm_load_si128((__m128i*)(desc_buffer));

    /* compute cost matrix */
    int32_t cost[25];
    for (int dv=0;dv<5;dv++) {
        for (int du=0;du<5;du++) {
            compute_small_des(I2_du,I2_dv,dims2[2],(int32_t)u2+du-2,(int32_t)v2+dv-2,desc_buffer);
            xmm2=_mm_load_si128((__m128i*)(desc_buffer));
            xmm2=_mm_sad_epu8(xmm1,xmm2);
            cost[dv*5+du]=_mm_extract_epi16(xmm2,0)+_mm_extract_epi16(xmm2,4);
        }
    }
    /* compute minimum */
    int32_t min_ind =0;
    int32_t min_cost=cost[0];
    for (int32_t i=1;i<25;i++) {
        if (cost[i]<min_cost) {
            min_ind=i; min_cost=cost[i];
        }
    }
    /* update target */
    u2+=(float)(min_ind%5)-2.0;
    v2+=(float)(min_ind/5)-2.0;
}
/* refinement-----------------------------------------------------------------*/
static void refinement(match_set_t *mset,const matchopt_t *opt)
{
    static double A[9*6]={ 1, 0, 1, 1, 0, 1, 1, 0, 1,
                           1, 1, 1, 0, 0, 0, 1, 1, 1,
                           1, 0,-1, 0, 0, 0,-1, 0, 1,
                          -1, 0, 1,-1, 0, 1,-1, 0, 1,
                          -1,-1,-1, 0, 0, 0, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1};
    double *At,*AtA;
    int i;
    unsigned char* desc_buffer=(uint8_t*)_mm_malloc(32*sizeof(unsigned char),16);

    trace(3,"refinement:\n");

    /* copy vector (for refill) */
    match_set_t mset_cpy={0};
    for (i=0;i<mset->n;i++) add_match_set_match_point(&mset_cpy,&mset->data[i]);
    free_match_set(mset);

    At=mat(9,6); AtA=mat(6,6);
    matt(A,6,9,At);
    matmul("NN",6,6,9,1.0,At,A,0.0,AtA);

    unsigned char* I1p_du_fit=I1p_du;
    unsigned char* I1p_dv_fit=I1p_dv;
    unsigned char* I1c_du_fit=I1c_du;
    unsigned char* I1c_dv_fit=I1c_dv;

    if (opt->half_res) {
        I1p_du_fit=I1p_du_full;
        I1p_dv_fit=I1p_dv_full;
        I1c_du_fit=I1c_du_full;
        I1c_dv_fit=I1c_dv_full;
    }
    /* for all matches do */
    for (i=0;i<mset_cpy.n;i++) {

        /* method: flow */
        if (opt->refine) {
            if (!parabolic_fit(I1c_du_fit,I1c_dv_fit,dims_c,I1p_du_fit,I1p_dv_fit,dims_p,
                               mset_cpy.data[i].uc,
                               mset_cpy.data[i].vc,
                               mset_cpy.data[i].up,
                               mset_cpy.data[i].vp,At,AtA,desc_buffer))
                continue;
        }
        else {
            relocate_minimum(I1c_du_fit,I1c_dv_fit,dims_c,I1p_du_fit,I1p_dv_fit,dims_p,
                             mset_cpy.data[i].uc,
                             mset_cpy.data[i].vc,
                             mset_cpy.data[i].up,
                             mset_cpy.data[i].vp,desc_buffer);
        }
        add_match_set_match_point(mset,&mset_cpy.data[i]);
    }
    /* free memory */
    free_match_set(&mset_cpy);
    free(At); free(AtA);
    free(desc_buffer);
}
/* match feature points from image -------------------------------------------*/
static int matchfeatures(const matchopt_t *opt,match_set_t *mset_sparse,
                         match_set_t *mset_dense,int *index,int *num,int *nindex,
                         int *nnum)
{
    trace(3,"match_features:\n");

#if FAST_CORNERS
    if (m1p2==NULL||n1p2==0||m1c2==NULL||n1c2==0) return 0;
#else
    /* flow match feature points */
    if (m1p2==NULL||n1p2==0||m1c2==NULL||n1c2==0) return 0;
    if (opt->multi_stage) {
        if (m1p1==NULL||n1p1==0||m1c1==NULL||n1c1==0) return 0;
    }
#endif
    /* clear old matches */
    free_match_set(mset_sparse);
    free_match_set(mset_dense);

    /* double pass matching */
    if (opt->multi_stage) {

#if FAST_CORNERS
        /* FAST corners match for dense features */
        match_internal(opt,mset_dense,m1p2,m1c2,n1p2,n1c2,false,1);

        if (opt->refine>0) {
            refinement(mset_dense,opt);
        }
        remove_outliers(mset_dense,opt,index,num,
                        nindex,nnum);
#else
        /* 1st pass (sparse matches) */
        match_internal(opt,mset_sparse,m1p1,m1c1,n1p1,n1c1,false,0);
        remove_outliers(mset_sparse,opt,NULL,NULL,
                        NULL,NULL);

        /* compute search range prior statistics (used for speeding up 2nd pass) */
        compute_prior_statistics(opt,mset_sparse,&mranges);

        /* 2nd pass (dense matches) */
        match_internal(opt,mset_dense,m1p2,m1c2,n1p2,n1c2,true,1);

        if (opt->refine>0) {
            refinement(mset_dense,opt);
        }
        remove_outliers(mset_dense,opt,index,num,
                        nindex,nnum);
#endif
    }
    else {
        
#if FAST_CORNERS
        /* FAST corners match for dense features */
        match_internal(opt,mset_dense,m1p2,m1c2,n1p2,n1c2,false,1);

        if (opt->refine>0) {
            refinement(mset_dense,opt);
        }
        remove_outliers(mset_dense,opt,index,num,
                        nindex,nnum);
#else
        /* single pass matching */
        match_internal(opt,mset_dense,m1p2,m1c2,n1p2,n1c2,false,1);
        if (opt->refine>0) {
            refinement(mset_dense,opt);
        }
        remove_outliers(mset_dense,opt,index,num,
                        nindex,nnum);
#endif
    }
    return mset_sparse->n;
}
/* init match ring buffer ----------------------------------------------------*/
static void init_match_buf(matchopt_t *opt)
{
    I1p_du_full=NULL; I1p_dv_full=NULL;
    I1c_du_full=NULL; I1c_dv_full=NULL;

    m1p1  =NULL; n1p1=0;
    m1p2  =NULL; n1p2=0;
    m1c1  =NULL; n1c1=0;
    m1c2  =NULL; n1c2=0;
    I1p   =NULL;
    I1c   =NULL;
    I1p_du=NULL; I1p_dv=NULL;
    I1c_du=NULL; I1c_dv=NULL;

    /* margin needed to compute descriptor + sobel responses */
    margin=8+1;

    /* adjust match radius on half resolution */
    if (opt->half_res) {
        opt->match_radius/=2;
    }
}
/* free match ring buffer-----------------------------------------------------*/
static void free_match_buf()
{
    if (I1p_du_full) _mm_free(I1p_du_full);
    if (I1p_dv_full) _mm_free(I1p_dv_full);
    if (I1c_du_full) _mm_free(I1c_du_full);
    if (I1c_dv_full) _mm_free(I1c_dv_full);

    if (m1p1) free(m1p1);
    if (m1p2) free(m1p2);
    if (m1c1) free(m1c1);
    if (m1c2) free(m1c2);

    if (I1p) _mm_free(I1p);
    if (I1c) _mm_free(I1c);
    if (I1p_du) _mm_free(I1p_du);
    if (I1c_du) _mm_free(I1c_du);
}
/* get random and unique sample of num numbers from 1:N ---------------------
 * args   :  int n        I  number of feature points
 *           int num      I  number of sample feature sample
 *           int *sample  O  index list of sample feature points
 * return : number of sampled feature points
 * --------------------------------------------------------------------------*/
static int getrand(int n,int num,int *sample)
{
    register int i,j,ns,*list;
    list=imat(n,1);

    /* create vector containing all indices */
    for (i=0;i<n;i++) list[i]=i;

    /* add num indices to current sample */
    for (ns=0,i=0;ns<num&&i<n;i++) {
        j=rand()%n;
        if (list[j]<0) continue;
        sample[ns++]=j; /* add sample index */
        list[j]=-1; /* disable this feature point */
    }
    free(list);
    return ns;
}
/* bucket features extract----------------------------------------------------*/
static void bucketfeat(const matchopt_t *opt,const match_set_t *mp_dense,
                       match_set_t *mp_bucket,int maxnf)
{
    register float u_max=0,bucket_w=(float)opt->bucket.w;
    register float v_max=0,bucket_h=(float)opt->bucket.h;
    register int i,u,v,*rlist,j,n,k;

    trace(3,"bucketfeat:\n");

    if (mp_dense->n<=0) {
        trace(2,"no matched feature points\n");
        return;
    }
    for (i=0;i<mp_dense->n;i++) {
        if (mp_dense->data[i].uc>u_max) u_max=mp_dense->data[i].uc;
        if (mp_dense->data[i].vc>v_max) v_max=mp_dense->data[i].vc;
    }
    /* allocate number of buckets needed */
    int bucket_cols=(int)floor(u_max/bucket_w)+1;
    int bucket_rows=(int)floor(v_max/bucket_h)+1;

    /* assign matches to their buckets */
    match_set *buckets=(match_set*)calloc(sizeof(match_set),bucket_cols*bucket_rows);
    for (i=0;i<mp_dense->n;i++) {
        u=(int)floor(mp_dense->data[i].uc/bucket_w );
        v=(int)floor(mp_dense->data[i].vc/bucket_h);

        add_match_set_match_point(&buckets[v*bucket_cols+u],&(mp_dense->data[i]));
    }
    /* refill p_matched from buckets */
    free_match_set(mp_bucket);

    for (i=0;i<bucket_cols*bucket_rows;i++) {

        /* shuffle bucket indices randomly */
        rlist=imat(1,buckets[i].n);

        n=getrand(buckets[i].n,maxnf>buckets[i].n?buckets[i].n:maxnf,rlist);
        for (j=0;j<n;j++) {
#if CHK_MATCH_ROI
            /* check match feature is in ROI? */
            if (!inroi(buckets[i].data[rlist[j]].uc,buckets[i].data[rlist[j]].vc,opt)) continue;
            if (!inroi(buckets[i].data[rlist[j]].up,buckets[i].data[rlist[j]].vp,opt)) continue;
#endif
            /* add match feature */
            add_match_set_match_point(mp_bucket,&buckets[i].data[rlist[j]]);
        }
        free(rlist);
    }
    for (i=0;i<bucket_cols*bucket_rows;i++) {
        free_match_set(buckets);
    }
}
/* bucket features point extract new version----------------------------------*/
static void buketfeatnew(const matchopt_t *opt,const match_set_t *mp_dense,
                         match_set_t *mp_bucket,int maxnf,int *index,int num,
                         int *nindex,int nnum)
{
    float um=0,bw=(float)opt->bucket.w;
    float vm=0,bh=(float)opt->bucket.h;
    int i,u,v,*rlist,j,n,bc,br;

    trace(3,"buketfeatnew:\n");

    /* first time to match feature */
    if (hash_counts(hash)<=0) {
        bucketfeat(opt,mp_dense,mp_bucket,maxnf);

        /* add match index to hash table */
        for (i=0;i<mp_bucket->n;i++) hash_add(&hash,mp_bucket->data[i].ic);
        return;
    }
    if (num==0) return;

    /* refill p_matched from buckets */
    free_match_set(mp_bucket);

    /* remove all elements in hash table */
    hash_delete(&hash);

    /* add matches from hash table */
    for (i=0;i<num;i++) {

#if CHK_MATCH_ROI
        /* check is in ROI? */
        if (!inroi(mp_dense->data[index[i]].uc,mp_dense->data[index[i]].vc,opt)) continue;
        if (!inroi(mp_dense->data[index[i]].up,mp_dense->data[index[i]].vp,opt)) continue;
#endif
        /* add match feature point */
        add_match_set_match_point(mp_bucket,&mp_dense->data[index[i]]);

        /* update hash table */
        hash_add(&hash,mp_dense->data[index[i]].ic);
    }
#if TRACR_FEAT_POINTS
    trace_match_points(mp_bucket);
#endif
    /* add new matches to hash table */
    if (nnum==0) return;
    for (i=0;i<nnum;i++) {
        if (mp_dense->data[nindex[i]].uc>um) um=mp_dense->data[nindex[i]].uc;
        if (mp_dense->data[nindex[i]].vc>vm) vm=mp_dense->data[nindex[i]].vc;
    }
    /* buckets fill matches */
    bc=(int)floor(um/bw)+1;
    br=(int)floor(vm/bh)+1;

    /* assign matches to their buckets */
    match_set *buckets=(match_set*)calloc(sizeof(match_set),bc*br);
    for (i=0;i<nnum;i++) {
        u=(int)floor(mp_dense->data[nindex[i]].uc/bw);
        v=(int)floor(mp_dense->data[nindex[i]].vc/bh);

        add_match_set_match_point(&buckets[v*bc+u],&(mp_dense->data[nindex[i]]));
    }
    /* refill matches from buckets */
    for (i=0;i<bc*br;i++) {

        /* shuffle bucket indices randomly */
        rlist=imat(1,buckets[i].n);

        n=getrand(buckets[i].n,maxnf>buckets[i].n?buckets[i].n:maxnf/2,rlist);
        for (j=0;j<n;j++) {

#if CHK_MATCH_ROI
            /* check is in ROI? */
            if (!inroi(buckets[i].data[rlist[j]].uc,buckets[i].data[rlist[j]].vc,opt)) continue;
            if (!inroi(buckets[i].data[rlist[j]].up,buckets[i].data[rlist[j]].vp,opt)) continue;
#endif
            add_match_set_match_point(mp_bucket,&buckets[i].data[rlist[j]]);

            /* add to hash table */
            hash_add(&hash,buckets[i].data[rlist[j]].ic);
        }
        free(rlist);
    }
    for (i=0;i<bc*br;i++) {
        free_match_set(buckets);
    }
}
/* add a feature to a frame---------------------------------------------------*/
static void hash_addfeat(feature** feats,int uid,double u,double v,gtime_t time)
{
    feature *s=NULL;
    HASH_FIND_INT(*feats,&uid,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(feature*)malloc(sizeof(feature));
        s->time=time;
        s->uid=uid;
        s->valid=1;
        s->u=u;
        s->v=v;
        HASH_ADD_INT(*feats,uid,s);  /* id: name of key field */
    }
}
/* copy feature point hash table in one image----------------------------------*/
static void hashcp(img_t *out,const img_t *in)
{
    feature *current,*tmp;
    HASH_ITER(hh,out->feat,current,tmp) {

        /* delete feature point */
        HASH_DEL(out->feat,current);
        free(current);
    }
    out->feat=NULL;
    HASH_ITER(hh,in->feat,current,tmp) {
        hash_addfeat(&out->feat,current->uid,current->u,current->v,current->time);
    }
}
/* copy image data to match_t struct-------------------------------------------*/
static void backupimg(const img_t *data,match_t *match)
{
    int size=data->h*data->w*sizeof(unsigned char);

    memcpy(match->Ip.data,match->Ic.data,size);
    memcpy(match->Ic.data,data->data,size);

    match->Ip.w   =match->Ic.w;
    match->Ip.h   =match->Ic.h;
    match->Ip.id  =match->Ic.id;
    match->Ip.time=match->Ic.time;

    match->Ic.w   =data->w;
    match->Ic.h   =data->h;
    match->Ic.id  =data->id;
    match->Ic.time=data->time;

    hashcp(&match->Ip,&match->Ic);
    hashcp(&match->Ic,data);
}
/* initial match struct-------------------------------------------------------
 * args:    match_t *match  IO  match struct
 *          matchopt_t *opt I   match options
 * return: none
 * ---------------------------------------------------------------------------*/
extern void init_match(match_t *match,const matchopt_t *opt)
{
    gtime_t t0={0};

    trace(3,"init_match:\n");

    match->opt=*opt;
    free_match_set(&match->mp_bucket);
    free_match_set(&match->mp_dense );
    free_match_set(&match->mp_sparse);

    if (match->Ip.data==NULL) {
        match->Ip.data=(unsigned char*)malloc(opt->img_w*opt->img_h*sizeof(unsigned char));
    }
    if (match->Ic.data==NULL) {
        match->Ic.data=(unsigned char*)malloc(opt->img_w*opt->img_h*sizeof(unsigned char));
    }
    match->Ic.w=opt->img_w;
    match->Ic.h=opt->img_h; match->Ic.time=t0;
    match->Ic.feat=NULL;

    match->Ip.w=opt->img_w;
    match->Ip.h=opt->img_h; match->Ip.time=t0;
    match->time=t0;
    match->Ip.feat=NULL;

    init_match_buf(&match->opt);
}
/* delete hash table-----------------------------------------------------------*/
static void hash_rmimgfeat(feature **ht)
{
    struct feature *current,*tmp;

    HASH_ITER(hh,*ht,current,tmp) {
        HASH_DEL(*ht,current);  /* delete; users advances to next */
        free(current);          /* optional- if you want to free  */
    }
    *ht=NULL; /* delete */
}
/* free match struct----------------------------------------------------------*/
extern void free_match(match_t *match)
{
    trace(3,"free_match:\n");
    if (match->Ic.data) free(match->Ic.data); match->Ic.data=NULL;
    if (match->Ip.data) free(match->Ip.data); match->Ip.data=NULL;

    free_match_set(&match->mp_bucket);
    free_match_set(&match->mp_sparse);
    free_match_set(&match->mp_dense );

    free_match_buf();
    hash_delete(&hash);

    hash_rmimgfeat(&match->Ip.feat);
    hash_rmimgfeat(&match->Ic.feat);
}
/* match feature point--------------------------------------------------------*/
static int matchimgfeat(match_set_t *pset,const img_t *pimg,const img_t *cimg)
{
    struct feature *current,*tmp,*pf=NULL;
    match_point pt={0};

    if (pimg->feat==NULL||cimg->feat==NULL) return 0;

    pset->n=0;
    HASH_ITER(hh,cimg->feat,current,tmp) {

        HASH_FIND_INT(pimg->feat,&current->uid,pf);
        if (pf==NULL) continue;
        pt.up=(float)pf->u;
        pt.vp=(float)pf->v;
        pt.uc=(float)current->u;
        pt.vc=(float)current->v;

        add_match_set_match_point(pset,&pt);
    }
    return pset->n;
}
/* match feature points from input image data---------------------------------
 * args:    match_t *m         IO  match struct data
 *          img_t *img         I   input image data
 * return: number of matched feature points
 * ---------------------------------------------------------------------------*/
extern int matchfeats(match_t *pmatch,const img_t *img)
{
    int dims[3],*index,num,*nindex,nnum;

    trace(3,"match: time=%s\n",time_str(img->time,4));

    /* ins-gnss-vo simulator data */
    if (img->feat) {
        pmatch->pt  =pmatch->time; 
        pmatch->time=img->time; 

        backupimg(img,pmatch);
        matchimgfeat(&pmatch->mp_bucket,&pmatch->Ip,&pmatch->Ic);

        /* numer of feature points */
        nnum=HASH_COUNT(img->feat);
        return nnum;
    }
    index =imat(1,img->h*img->w/2);
    nindex=imat(1,img->h*img->w/2);

    dims[0]=img->w; dims[1]=img->h;
    dims[2]=img->w;

    pmatch->pt  =pmatch->time;
    pmatch->time=img->time; /* update time */
    
    /* backup image buffer data */
    backupimg(img,pmatch);

    /* input image data */
    if (!puchbackimg(img->data,dims,greplace,&pmatch->opt)) {
        trace(2,"match feature fail\n");
        return 0;
    }
    /* match feature points */
    matchfeatures(&pmatch->opt,&pmatch->mp_sparse,&pmatch->mp_dense,index,&num,nindex,&nnum);

#if USE_NEW_BUCKET
    /* bucket features extract */
    buketfeatnew(&pmatch->opt,&pmatch->mp_dense,&pmatch->mp_bucket,
                 pmatch->opt.bucket.nmax,
                 index,num,nindex,nnum);
#else
    /* bucket features extract */
    bucketfeat(&pmatch->opt,&pmatch->mp_dense,&pmatch->mp_bucket,
               pmatch->opt.bucket.nmax);
#endif

#if TRACR_FEAT_POINTS
    trace_match_points(&pmatch->mp_bucket);
#endif
    /* check match ok? */
    if (pmatch->mp_dense.n<=0||pmatch->mp_bucket.n<=0) {

        free(index); free(nindex);
        trace(2,"match feature fail\n");
        return 0;
    }
    free(index); free(nindex);
    return pmatch->mp_bucket.n;
}
/* initial image struct data-------------------------------------------------*/
extern int initimg(img_t *data,int w,int h,gtime_t time)
{
    trace(3,"initimg:\n");
    if (w<=0||h<=0) return 0;

    if (data->data) free(data->data); data->data=NULL;
    if (!(data->data=(uint8_t*)malloc(sizeof(uint8_t)*w*h))) {
        trace(2,"initial ima fail\n");
        return 0;
    }
    data->time=time; data->w=w; data->h=h;
    data->id=0;
    data->feat=NULL;
    return 1;
}
/* free image struct data----------------------------------------------------*/
extern void freeimg(img_t *data)
{
    trace(3,"freeimg:\n");

    if (data->data) {
        free(data->data); data->data=NULL;
    }
    data->w=data->h=0;
    hash_rmimgfeat(&data->feat);
    data->feat=NULL;
}


