/*--------------------------------------------------------------------------------
 * ins-gnss-vo-msckf.cc : ins-gnss-vo coupled common functions
 *
 * reference :
 *    [01] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *         Navigation System, Artech House, 2008
 *    [02] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *         for IMU calibration without external equipments,2014.
 *    [03] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *         INS 2008.
 *    [04] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [05] Li M, Mourikis A I. High-precision, consistent EKF-based visual–inertial
 *         odometry[J].International Journal of Robotics Research,2013,32(6):690-711.
 *    [06] Monocular Visual Inertial Odometry on a Mobile Device.
 *    [07] Mourikis A / , Roumeliotis S / . A Multi-State Constraint Kalman Filter
 *         for Vision-aided Inertial Navigation[C]// IEEE International Conference
 *         on Robotics & Automation. IEEE, 2007.
 *    [08] Forster C , Carlone L , Dellaert F , et al. On-Manifold Preintegration
 *         for Real-Time Visual-Inertial Odometry[J]. IEEE Transactions on Robotics,
 *         2015, 33(1):1-21.
 *    [09] Observability-constrained vision-aided inertial navigation.
 *    [10] Li M ,Mourikis A I. Improving the accuracy of EKF-based visual-inertial
 *         odometry[C]// IEEE International Conference on Robotics & Automation.
 *         IEEE, 2012.
 *    [11] Li M ,Mourikis A I. High-precision, consistent EKF-based visual-inertial
 *         odometry[M].Sage Publications, Inc. 2013.
 *    [12] Pizzoli M , Forster C , Scaramuzza D . REMODE: Probabilistic, Monocular
 *         Dense Reconstruction in Real Time[C]// IEEE International Conference on
 *         Robotics and Automation (ICRA), Hong Kong, 2014. IEEE, 2014.
 *    [13] Li M , Kim B H , Mourikis A I . Real-time motion tracking on a cellphone
 *         using inertial sensing and a rolling-shutter camera[C]// IEEE International
 *         Conference on Robotics & Automation. IEEE, 2013.
 *    [14] Monocular Visual-Inertial SLAM and Self Calibration for Long Term Autonomy
 *    [15] Lupton T , Sukkarieh S . Visual-Inertial-Aided Navigation for High-Dynamic
 *         Motion in Built Environments Without Initial Conditions[J]. IEEE Transactions
 *         on Robotics, 2012, 28(1):61-76.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants-------------------------------------------------------------------*/
#define MIN_TRACK_LEN               2         /* min length of tracking data */
#define MAX_TRACK_LEN               50        /* max length of tracking data */
#define VAR_FEAT                    SQR(3.0)  /* variance of feature point in image coordinates */
#define SWAP(type,x,y)              {type tmp; tmp=x; x=y; y=tmp;}
#define MAXITER                     10        /* max iteration for Gauss Newton optimization */
#define MAX_FEAT_RES                10.0      /* max residual for a feature point measurement */
#define MAX_GNCOST_NORM             0.03      /* set to inf to allow any triangulation, no matter how bad */
#define OUT_DETECT                  1         /* Mahalanobis gating test for the residual */
#define USE_BUCKET_FEATS            1         /* use bucket matched feature points data */
#define POST_VALIDATION             0         /* post-fit residuals for validation solution */
#define BUCKET_WIDTH                40        /* width of bucket for bucketing matches in pixel */
#define BUCKET_HEIGHT               40        /* height of bucket for bucketing matches in pixel */
#define BUCKET_ROBUST               0         /* bucket matches for robust estimate states */
#define MAXVOGDOP                   30.0      /* reject threshold of gdop */
#define MAX_RANSAC_ITER             100       /* max number of RANSAC for features measurement updates */
#define RATIO_THRESHOLD             0.6       /* ratio threshold for epipolar constraint check */
#define DO_QR_DCMP                  0         /* do QR decomposition */
#define DO_ANOTHER_UPD              1         /* use another method to update feature points */
#define DO_MONO_VO_EST              1         /* monocular visual odometry estimate */
#define DO_LM_METHOD                1         /* apply Levenberg-Marquart method to solve for the 3d position */

/* type definitions ----------------------------------------------------------*/
typedef struct hashtable {                    /* hash table type */
    int id,index;                             /* id of track feature point/index in tracks data */
    const trackd_t *ptrk;                     /* pointer of this track feature point */
    UT_hash_handle hh;                        /* makes this structure hashable */
} hashtable_t;

typedef struct cams {                         /* camera pose states type */
    gtime_t time,imgt;                        /* timestamp of ins states/image data */
    long int find;                            /* frame index corresponding to camera */
    double pnull[3],Cnull[9];                 /* these two variables should have the same physical
                                               * interpretation with `orientation` and `position`.
                                               * There two variables are used to modify the measurement
                                               * Jacobian matrices to make the observability matrix
                                               * have proper null space.
                                               */
    double re[3],Cce[9];                      /* camera position/velocity/attitude states in ecef */
    hashtable *trackfeat;                     /* tracking feature points */
} cams_t;

typedef struct vofilt {                       /* vo filter workspace type */
    int n,nmax;                               /* number and max number of history ins states */
    int nx;                                   /* number error states of current sliding windows */
    double *Px,r1,r2,r3;                      /* current ins states covariance included camera pose in sliding window */
    cams_t *data;                             /* all tracked camera states */
} vofilt_t;

typedef struct bucket {                       /* track features buckets type */
    hashtable *trackfeat;                     /* track features data */
} bucket_t;

typedef struct configLM {                     /* configuration parameters for 3d feature position optimization */
    double tth;                               /* translation threshold */
    double huber;                             /* huber epsilon */
    double ep;                                /* estimation precision */
    double initial_damping;                   /* initial damping */
    int outloop_maxiter,innloop_maxiter;      /* outer/inner loop max iteration */
} configLM_t;

/* global variables------------------------------------------------------------*/
static vofilt_t vofilt={0};                   /* vo aid filter workspace */
static track_t  tracks={0};                   /* all tracking feature points data */
static match_t  matchs={0};                   /* match feature points data */
static const insstate_t *inss=NULL;           /* pointer to ins states */
static const insopt_t   *opts=NULL;           /* pointer to ins options */
static const configLM LMconfig={-1.0,
                                0.01,
                                5E-7,
                                1E-3,
                                5,
                                5};           /* configuration parameters for 3d feature position optimization */

/* add a track feature to hash table------------------------------------------*/
static void hash_add_feature(struct hashtable **ht,const int id,const trackd_t *ptk)
{
    struct hashtable *s=NULL;
    HASH_FIND_INT(*ht,&id,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(struct hashtable *)malloc(sizeof *s);
        s->id=id;
        s->ptrk=ptk;
        HASH_ADD_INT(*ht,id,s);  /* id: name of key field */
    }
}
static void hash_add_featidx(struct hashtable **ht,const int id,const int idx,
                             const trackd_t *ptk)
{
    struct hashtable *s=NULL;
    HASH_FIND_INT(*ht,&id,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(struct hashtable *)malloc(sizeof(struct hashtable));
        s->index=idx;
        s->id=id;
        s->ptrk=ptk;
        HASH_ADD_INT(*ht,id,s);  /* id: name of key field */
    }
}
/* find a track feature from hash table---------------------------------------*/
static hashtable* hash_find_feature(const struct hashtable *ht,const int id)
{
    struct hashtable *s=NULL;

    HASH_FIND_INT(ht,&id,s);  /* s: output pointer */
    return s;
}
/* remove element from hash table---------------------------------------------*/
static void hash_delete_feature(struct hashtable **ht,const int id)
{
    struct hashtable *s=NULL;
    if (!(s=hash_find_feature(*ht,id))) {
        trace(2,"no feature to delete\n");
        return;
    }
    HASH_DEL(*ht,s);  /* user: pointer to delete */
    free(s);
}
/* delete hash table----------------------------------------------------------*/
static void hash_destroy(struct hashtable **ht)
{
    struct hashtable *current,*tmp;

    HASH_ITER(hh,*ht,current,tmp) {
        HASH_DEL(*ht,current);  /* delete; users advances to next */
        free(current);          /* optional- if you want to free  */
    }
}
/* create a hash table--------------------------------------------------------*/
static struct hashtable *hashtable()
{
    return NULL;
}
/* counts of elements in hash table-------------------------------------------*/
static int hash_counts(const struct hashtable *ht)
{
    return HASH_COUNT(ht);
}
/* get hash table element by index--------------------------------------------*/
static struct hashtable *hash_index(struct hashtable **ht,int index)
{
    struct hashtable *current,*tmp;
    int i=0;

    HASH_ITER(hh,*ht,current,tmp) {
        if (i++==index) return current;
    }
    return NULL;
}
/* initial vo-aid-------------------------------------------------------------*/
extern void initvoaid(insopt_t *opt)
{
    trace(3,"initvoaid:\n");

    opt->voopt.match.f =opt->voopt.calib.f;
    opt->voopt.match.fu=opt->voopt.calib.fu;
    opt->voopt.match.fv=opt->voopt.calib.fv;

    opt->voopt.match.cu=opt->voopt.calib.cu;
    opt->voopt.match.cv=opt->voopt.calib.cv;

    init_match(&matchs,&opt->voopt.match);
    inittrackimgbuf(&opt->voopt);

    vofilt.nx=xnX(opt);
    vofilt.Px=zeros(vofilt.nx,vofilt.nx); vofilt.n=vofilt.nmax=0;
    return;
}
/* free vo-aid----------------------------------------------------------------*/
extern void freevoaid()
{
    trace(3,"freevoaid:\n");
    int i;

    freetrackset(&tracks);
    free_match(&matchs);
    if (vofilt.Px) free(vofilt.Px);

    for (i=0;i<vofilt.n;i++) {
        hash_destroy(&vofilt.data[i].trackfeat);
    }
    if (vofilt.data) {
        free(vofilt.data);
    }
    freetrackimgbuf(); return;
}
/* get camera pose using current ins stats------------------------------------*/
static void campose(const insstate_t *ins,const insopt_t *opt,const img_t *img,cams_t *cams,
                    double *Jp,double *Ja,double *Jv)
{
    static double I[9]={1,0,0,0,1,0,0,0,1};
    double T[3],W[9];
    int i,j;

    static int ip=xiP(opt),ia=xiA(opt),ilc=xiCl (opt);
    static int np=xnP(opt),na=xnA(opt),nlc=xnCla(opt);
    static int ila=xiCm(opt),nla=xnCm(opt);

    /* compute camera pose and velocity */
    if (cams) {
        matmul("NT",3,3,3,1.0,ins->Cbe,ins->Cbc,0.0,cams->Cce);
        matmul("NN",3,1,3,1.0,ins->Cbe,ins->lbc,0.0,T);
        for (i=0;i<3;i++) {
            cams->re[i]=ins->re[i]+T[i];
        }
        /* update camera states nullspace */
        matcpy(cams->pnull,cams->re ,3,1);
        matcpy(cams->Cnull,cams->Cce,3,3);

        cams->time=ins->time;
        cams->imgt=img->time;
    }
    /* jacobian for camera position */
    if (Jp) {
        skewsym3(ins->lbc,W);
        matmul("NN",3,3,3,1.0,ins->Cbe,W,0.0,T);

        /* jacobians wrt. camera attitude */
        for (i=0;i<3;i++) {
            for (j=ia;j<na+ia;j++) Jp[i+j*3]=T[i+(j-ia)*3];
        }
        /* jacobians wrt. camera position */
        for (i=0;i<3;i++) {
            for (j=ip;j<np+ip;j++) Jp[i+j*3]=I[i+(j-ip)*3];
        }
        /* jacobians wrt. camera lever arm  */
        if (opt->estcaml) {
            for (i=0;i<3;i++) {
                for (j=ilc;j<ilc+nlc;j++) {
                    Jp[i+j*3]=ins->Cbe[i+(j-ilc)*3];
                }
            }
        }
        trace(3,"Jp=\n");
        tracemat(3,Jp,3,ins->nx,12,6);
    }
    /* jacobian for camera attitude */
    if (Ja) {
        /* wrt. camera attitude */
        for (i=0;i<3;i++) for (j=ia;j<ia+na;j++) {
                Ja[i+j*3]=ins->Cbc[i+j*3];
            }

        /* wrt. camera-ins misalignment */
        if (opt->estcama) {
            for (i=0;i<3;i++) {
                for (j=ila;j<ila+nla;j++) Ja[i+j*3]=I[i+j*3];
            }
        }
        trace(3,"Ja=\n");
        tracemat(3,Ja,3,ins->nx,12,6);
    }
}
/* add a ins state to filter--------------------------------------------------*/
static int addinstates(vofilt_t *filt,const insstate_t *in,const insopt_t *opt,
                       const img_t *img,int cur_fid)
{
    cams_t cams,*p;
    trace(3,"addinstat: time=%s\n",time_str(in->time,4));

    /* compute camera pose */
    campose(in,opt,img,&cams,NULL,NULL,NULL);

    /* create hash table */
    cams.trackfeat=hashtable();
    cams.find=cur_fid;

    /* add ins states */
    if (filt->nmax<=filt->n) {
        if (filt->nmax<=0) filt->nmax=16*2; else filt->nmax*=2;
        if (!(p=(cams_t*)realloc(filt->data,sizeof(cams_t)*filt->nmax))) {

            /* add ins fail */
            filt->n=filt->nmax=0; free(filt->data);
            filt->data=NULL;
            return -1;
        }
        filt->data=p;
        int nxc=filt->nmax*6+in->nx;
        if (!resize(&filt->Px,filt->nx,filt->nx,nxc,nxc)) {
            trace(2,"Py resize matrix fail\n");

            free(filt->Px); filt->Px=NULL;
            return -1;
        }
    }
    filt->data[filt->n++]=cams;
    return 1;
}
/* update the state correspondes to null space ------------------------------*/
static void updnullstate(insstate_t *ins,const imud_t *data)
{
    matcpy(ins->pnull,ins->re ,1,3);
    matcpy(ins->vnull,ins->ve ,1,3);
    matcpy(ins->Cnull,ins->Cbe,3,3);
}
/* modify the transition matrix-----------------------------------------------*/
static void mdfyF(insstate_t *ins,const insopt_t *opt,double dt,double *Fm)
{
    double Rk[9],g[3],u[3],s[3],a;
    double w1[3],w2[3],A1[9],A2[9],dv[3],dp[3];
    double t[3],T[9],W[9];
    int i,j;

    /* initial nullspace */
    if (norm(ins->Cnull,3)==0.0) matcpy(ins->Cnull,ins->pCbe,3,3);
    if (norm(ins->pnull,3)==0.0) matcpy(ins->pnull,ins->pins,1,3);
    if (norm(ins->vnull,3)==0.0) matcpy(ins->vnull,ins->pins+3,1,3);

    /* for attitude */
    matmul("TN",3,3,3,1.0,ins->Cbe,ins->Cnull,0.0,Rk);
    pregrav(ins->pnull,g);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) Fm[xiA(opt)+i+(xiA(opt)+j)*ins->nx]=Rk[i+j*3];
    }
    matmul("NN",3,1,3,1.0,Rk,g,0.0,u);
    matmul("TN",1,1,3,1.0,u,u,0.0,&a);
    s[0]=u[0]/a;
    s[1]=u[1]/a;
    s[2]=u[2]/a;

    /* for velocity */
    dv[0]=ins->vnull[0]-ins->ve[0];
    dv[1]=ins->vnull[1]-ins->ve[1];
    dv[2]=ins->vnull[2]-ins->ve[2];
    skewsym3(dv,W);
    matmul("NN",3,1,3,1.0,W,g,0.0,w1);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            A1[i+j*3]=Fm[xiV(opt)+i+(xiA(opt)+j)*ins->nx];
        }
    }
    matmul("NN",3,1,3,1.0,A1,u,0.0,t);
    t[0]=t[0]-w1[0];
    t[1]=t[1]-w1[1];
    t[2]=t[2]-w1[2];

    matmul("NT",3,3,1,1.0,t,s,0.0,T);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            Fm[xiV(opt)+i+(xiA(opt)+j)*ins->nx]=A1[i+j*3]-T[i+j*3];
        }
    }
    /* for position */
    dp[0]=ins->pnull[0]+ins->vnull[0]*dt-ins->re[0];
    dp[1]=ins->pnull[1]+ins->vnull[1]*dt-ins->re[1];
    dp[2]=ins->pnull[2]+ins->vnull[2]*dt-ins->re[2];
    skewsym3(dp,W);
    matmul("NN",3,1,3,1.0,W,g,0.0,w2);

    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            A2[i+j*3]=Fm[xiP(opt)+i+(xiA(opt)+j)*ins->nx];
        }
    }
    matmul("NN",3,1,3,1.0,A2,u,0.0,t);
    t[0]=t[0]-w2[0];
    t[1]=t[1]-w2[1];
    t[2]=t[2]-w2[2];

    matmul("NT",3,3,1,1.0,t,s,0.0,T);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            Fm[xiP(opt)+i+(xiA(opt)+j)*ins->nx]=A2[i+j*3]-T[i+j*3];
        }
    }
    /* update the state correspondes to null space */
    updnullstate(ins,NULL);

    trace(3,"Fm=\n");
    tracemat(3,Fm,ins->nx,ins->nx,12,5);
}
/* the transition matrix------------------------------------------------------*/
static void statesF(const insopt_t *opt,const insstate_t *ins,double *F,double dt)
{
    static double I[9]={1,0,0,0,1,0,0,0,1};
    double W[9],T[9],fe[3];
    double *F0;
    int i,j,nx=ins->nx;

    int ibg=xiBg(opt),nbg=xnBg(opt);
    int iba=xiBa(opt),nba=xnBa(opt);
    int ixa=xiA(opt),nxa=xnA(opt);
    int ixv=xiV(opt),nxv=xnV(opt);
    int ixp=xiP(opt),nxp=xnP(opt);

    F0=zeros(nx,nx);

    /* for ins attitude */
    skewsym3(ins->omgb,W);
    for (i=0;i<3;i++) {
        for (j=ixa;j<ixa+nxa;j++) F0[i+j*nx]= W[i+(j-ixa)*3];
        for (j=ibg;j<ibg+nbg;j++) F0[i+j*nx]=-I[i+(j-ibg)*3];
    }
    /* for velocity */
    matmul("NN",3,1,3,1.0,ins->Cbe,ins->fb,0.0,fe);
    skewsym3(fe,W);
    matmul("NN",3,3,3,1.0,ins->Cbe,W,0.0,T);
    for (i=ixv;i<ixv+nxv;i++) {
        for (j=iba;j<iba+nba;j++) F0[i+j*nx]=-ins->Cbe[i-ixv+(j-iba)*3];
        for (j=ixa;j<ixa+nxa;j++) F0[i+j*nx]=-T[i-ixv+(j-ixa)*3];
    }
    /* for position */
    for (i=ixp;i<ixp+nxp;i++) {
        for (j=ixv;j<ixv+nxv;j++) {
            F0[i+j*nx]=I[i-ixp+(j-ixv)*3];
        }
    }
    trace(3,"F0=\n");
    tracemat(3,F0,nx,nx,12,5);

    for (i=0;i<nx*nx;i++) F0[i]*=dt;

    expmat(F0,nx,F);

    trace(3,"F(1)=\n");
    tracemat(3,F,nx,nx,12,5);
    free(F0);
}
/* noise transition matrix----------------------------------------------------*/
static void noiseG(const insopt_t *opt,const insstate_t *ins,double *G)
{
    static double I[9]={1,0,0,0,1,0,0,0,1};
    int i,j,nx=ins->nx;

    int ibg=xiBg(opt),nbg=xnBg(opt);
    int iba=xiBa(opt),nba=xnBa(opt);
    int ixa=xiA(opt) ,nxa=xnA(opt);
    int ixv=xiV(opt) ,nxv=xnV(opt);

    /* for attitude */
    for (i=ixa;i<ixa+nxa;i++) {
        for (j=0;j<3;j++) G[i+j*nx]=-I[i-ixa+j*3];
    }
    /* for bg. */
    for (i=ibg;i<ibg+nbg;i++) {
        for (j=3;j<6;j++) G[i+j*nx]=I[i-ibg+(j-3)*3];
    }
    /* for velocity */
    for (i=ixv;i<ixv+nxv;i++) {
        for (j=6;j<9;j++) {
            G[i+j*nx]=-ins->Cbe[i-ixv+(j-6)*3];
        }
    }
    /* for ba. */
    for (i=iba;i<iba+nba;i++) {
        for (j=9;j<12;j++) G[i+j*nx]=I[i-iba+(j-9)*3];
    }
    trace(3,"G=\n");
    tracemat(3,G,nx,4*3,12,6);
}
/* update states covariance---------------------------------------------------*/
static void updallcov(const insopt_t *opt,const double *P0,const double *F,
                      const double *G, double dt, double *P,
                      int nx,int ns)
{
    double *Qs,*Ts,*Q,*P1;
    int i,j;

    trace(3,"P(0)=\n");
    tracemat(3,P,nx,nx,12,6);

    Qs=zeros(ns,ns); Q=zeros(nx,nx);
    Ts=zeros(nx,ns);
    P1=zeros(nx,nx);

    for (i=0;i<3 ;i++) Qs[i+i*ns]=opt->psd.gyro*fabs(dt);
    for (i=6;i<9 ;i++) Qs[i+i*ns]=opt->psd.accl*fabs(dt);
    for (i=9;i<12;i++) {
        Qs[i+i*ns]=opt->psd.ba*fabs(dt);
    }
    for (i=3;i<6 ;i++) {
        Qs[i+i*ns]=opt->psd.bg*fabs(dt);
    }
    matmul("NN",nx,ns,nx,1.0,F,G,0.0,Ts);
    matmul33("NNT",Ts,Qs,Ts,nx,ns,ns,nx,Q);

    trace(3,"Q=\n");
    tracemat(3,Q,nx,nx,12,5);

    matmul33("NNT",F,P0,F,nx,nx,nx,nx,P1);
    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) P[i+j*nx]=P1[i+j*nx]+Q[i+j*nx];
    }
    trace(3,"P(1)=\n");
    tracemat(3,P,nx,nx,12,6);
    free(Qs); free(P1);
    free(Ts); free(Q );
}
static void Pe2Pb(const insstate_t *ins,const insopt_t *opt,double *P)
{
    double *Pa0=zeros(3,3),*Pa1=zeros(3,3);
    int i,j,nx=ins->nx;

    matcpy(P,ins->P,nx,nx);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) Pa0[i+j*3]=ins->P[xiA(opt)+i+(xiA(opt)+j)*nx];
    }
    matmul33("TNN",ins->Cbe,Pa0,ins->Cbe,3,3,3,3,Pa1);
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            P[xiA(opt)+i+(xiA(opt)+j)*nx]=Pa1[i+j*3];
        }
    }
    free(Pa0);
    free(Pa1);
}
/* propagate covariance matrix of ins-camera pose in sliding windows----------*/
static int propagate(vofilt_t *filt,insstate_t *ins,const insopt_t *opt)
{
    double *T,*Pp,*F,*G,*Pip,*Pi;
    int i,j,nx=filt->nx;
    
    trace(3,"propagate: nx=%d\n",filt->nx);

    /* first epoch to initial covariance */
    if (filt->n==0||filt->nx==ins->nx) {
        Pip=zeros(nx,nx);

        Pe2Pb(ins,opt,Pip);
        for (i=0;i<nx;i++) for (j=0;j<nx;j++) filt->Px[i+j*nx]=Pip[i+j*nx];
        free(Pip);
        return 1;
    }
    /* modify the transition matrix */
    F=zeros(ins->nx,ins->nx);
    G=zeros(ins->nx,12);
    Pp=mat(nx,nx); T=mat(nx,nx);
    Pip=zeros(nx,nx);
    Pi =zeros(nx,nx);

    statesF(opt,ins,F,ins->dt);
    mdfyF(ins,opt,ins->dt,F);
    noiseG(opt,ins,G);

    /* update states covariance */
    Pe2Pb(ins,opt,Pip);
    updallcov(opt,Pip,F,G,ins->dt,Pi,ins->nx,12);
    for (i=0;i<nx-ins->nx;i++) {
        for (j=0;j<ins->nx;j++) T[i+j*(nx-ins->nx)]=filt->Px[ins->nx+i+j*nx];
    }
    /* update camera states covariance */
    matmul("NT",nx-ins->nx,ins->nx,ins->nx,1.0,T,F,0.0,Pp);
    for (i=0;i<nx-ins->nx;i++) {
        for (j=0;j<ins->nx;j++) {
            filt->Px[(ins->nx+i)*nx+j]=filt->Px[ins->nx+i+j*nx]=Pp[i+j*(nx-ins->nx)];
        }
    }
    /* update ins states covariance */
    for (i=0;i<ins->nx;i++) {
        for (j=0;j<ins->nx;j++) {
            filt->Px[i+j*nx]=Pi[i+j*ins->nx];
        }
    }
    trace(3,"Px=\n");
    tracemat(3,filt->Px,nx,nx,12,6);

    free(T); free(G );
    free(F); free(Pi);
    free(Pip); free(Pp);
    return 1;
}
/* add a ins state to vo-filter workspace-------------------------------------*/
static void augstates(vofilt_t *filt,const insstate_t *ins,const insopt_t *opt,
                      const img_t *img,int cur_fid)
{
    double *J,*Jp,*Ja,*P,*I;
    int i,j,nx=ins->nx;
    int nxp=filt->nx,nxc=filt->nx+6;

    trace(3,"augstates:\n");

    /* first camera states */
    if (filt->n==0||filt->nx==ins->nx) {
        for (i=0;i<nx;i++) for (j=0;j<nx;j++) filt->Px[i+j*nx]=ins->P[i+j*nx];
    }
    /* add a camera pose to sliding window */
    if (addinstates(filt,ins,opt,img,cur_fid)<0) {
        trace(2,"add camera pose fail\n");
        return;
    }
    P=zeros(nxc,nxc); J=zeros(nxc,nxp);
    Jp=zeros(3,nx); Ja=zeros(3,nx);
    I=eye(nxp);

    /* jacobians matrix */
    campose(ins,opt,img,NULL,Jp,Ja,NULL);

    asi_blk_mat(J,nxc,nxp,Ja,3,nx,nxp+0,0);
    asi_blk_mat(J,nxc,nxp,Jp,3,nx,nxp+3,0);
    asi_blk_mat(J,nxc,nxp,I,nxp,nxp,0,0);
    matmul33("NNT",J,filt->Px,J,nxc,nxp,nxp,nxc,P);

    /* augment covariance matrix */
    for (i=0;i<nxc;i++) {
        for (j=0;j<nxc;j++) filt->Px[i+j*nxc]=P[i+j*nxc];
    }
    /* update nx of `Px' */
    filt->nx=nxc;

    free(J); free(Jp); free(Ja);
    free(P); free(I);
}
/* clone camera---------------------------------------------------------------*/
static void clonecamera(const cams_t *in,cams_t *out)
{
    struct hashtable *current,*tmp;
    if (in==NULL) return;
    out->find=in->find;
    out->imgt=in->imgt;
    out->time=in->time;

    matcpy(out->Cce,in->Cce,3,3);
    matcpy(out->re ,in->re ,3,1);

    hash_destroy(&out->trackfeat);
    HASH_ITER(hh,in->trackfeat,current,tmp) {
        hash_add_featidx(&out->trackfeat,current->id,current->index,current->ptrk);
    }
    return;
}
/* handle of remove camera pose-----------------------------------------------*/
static void handlerm(vofilt_t *filt,const insstate_t *ins,const insopt_t *opt,
                    const int idx)
{
    int i,j,nxc,nxp=filt->nx,m,n;
    double *Pp;

    Pp=mat(nxp,nxp);

    trace(3,"remove camera id=%4ld\n",filt->data[idx].find);

    for (i=idx;i<filt->n-1;i++) {
        clonecamera(&filt->data[i+1],&filt->data[i]);
    }
    if (idx==filt->n-1) {
        hash_destroy(&filt->data[idx].trackfeat);
        filt->data[idx].find=-1;
    }
    matcpy(Pp,filt->Px,nxp,nxp);
    nxc=ins->nx+6*(--filt->n);
    for (n=0,i=0;i<nxp;i++) {
        if (i>=ins->nx+6*idx&&i<=ins->nx+6*idx+5) continue;
        for (m=0,j=0;j<nxp;j++) {
            if (j>=ins->nx+6*idx&&j<=ins->nx+6*idx+5) {
                continue;
            }
            filt->Px[n+m++*nxc]=Pp[i+j*nxp];
        }
        n++;
    }
    filt->nx=nxc;
    free(Pp);
}
/* trace camera tracking features---------------------------------------------*/
static void tracecamera(const vofilt_t *vo)
{
    int i,nf;
    trace(3,"camera tracking features:\n");

    for (i=0;i<vo->n;i++) {
        nf=hash_counts(vo->data[i].trackfeat);
        fprintf(stderr,"cam-id=%3ld  nf=%3d\n",vo->data[i].find,nf);
    }
}
/* remove camera pose from sliding window-------------------------------------*/
static int rmcamera(vofilt_t *filt,const insstate_t *ins,const insopt_t *opt)
{
    trace(3,"rmcamera:\n");
    int i;
    tracecamera(filt);

    for (i=0;i<filt->n;i++) {
        if (i>=0&&hash_counts(filt->data[i].trackfeat)) continue;
        hash_destroy(&filt->data[i].trackfeat);

        /* remove camera */
        handlerm(filt,ins,opt,i);
        i--;
    }
    return 1;
}
/* update all track data------------------------------------------------------*/
static int updatetrack(const insopt_t *opt,const img_t *img)
{
    match_set *pset=&matchs.mp_dense;

#if USE_BUCKET_FEATS
    pset=&matchs.mp_bucket;
#endif
    /* update track data */
    if (!match2track(pset,matchs.pt,matchs.time,img->id,&matchs.Ip,&matchs.Ic,
                     &opt->voopt,&tracks)) {

        trace(2,"update track fail\n");
        return 0;
    }
    return 1;
}
/* get camera pointer through given `uid'-------------------------------------*/
static cams_t* getcampointor(int uid)
{
    cams_t *pcam=NULL;
    int i;
    for (i=0;i<vofilt.n;i++) {
        if (vofilt.data[i].find==uid) {pcam=&vofilt.data[i]; break;}
    }
    return pcam;
}
/* trace tracking length of all feature points--------------------------------*/
static void tracetracklen(const track_t *track)
{
    static int ttable[MAX_TRACK_LEN][MAXBUFF];
    static int ntt[MAX_TRACK_LEN];             /* track feature table which store number and index of
                                                * tracking `min-track-len' to `max-track-len'
                                                * */
    trace(3,"tracetracklen:\n");

    int i,j;
    for (i=0;i<MAX_TRACK_LEN;i++) ntt[i]=0;

    /* get number of tracking lost */
    for (i=0;i<track->nlost;i++) {
        j=(tracks.data[tracks.losttrack[i]].n-MIN_TRACK_LEN)%MAX_TRACK_LEN;
        ttable[j][ntt[j]++]=i;
    }
    trace(3,"track lost: %d\n",tracks.nlost);

    trace(3,"track table:\n");
    for (i=0;i<MAX_TRACK_LEN;i++) {
        trace(3,"track len=%2d:  %4d\n",i+MIN_TRACK_LEN,ntt[i]);
    }
}
/* add current tracking feature points to current camera----------------------*/
static int addcurrenttrack(cams_t *cam,const int pre_id)
{
    cams_t *pcam=NULL;
    int i,j;

    trace(3,"add new tracks: %d\n",tracks.nnew);
    trace(3,"update tracks: %d\n",tracks.nupd);

    /* get number of tracking lost */
    for (i=0,tracks.nlost=0;i<tracks.n;i++) {
        if (tracks.data[i].flag==TRACK_LOST&&tracks.data[i].last_idx!=-1&&tracks.data[i].n>=MIN_TRACK_LEN) {
            tracks.losttrack[tracks.nlost++]=i;
        }
    }
    trace(3,"track lost: %d\n",tracks.nlost);

    /* exceed max track length */
    for (i=0,tracks.nexceed=0;i<tracks.n;i++) {
        if (tracks.data[i].n>=MAX_TRACK_LEN&&tracks.data[i].flag==TRACK_UPDATED) tracks.exceedtrack[tracks.nexceed++]=i;
    }
    trace(3,"track exceed: %d\n",tracks.nexceed);

    /* add new track to camera tracking-feature list */
    for (i=0,pcam=getcampointor(pre_id);i<tracks.nnew&&pcam;i++) {

        j=tracks.newtrack[i];
        hash_add_feature(&cam->trackfeat,tracks.data[j].uid,&tracks.data[j]);

        if (hash_find_feature(pcam->trackfeat,tracks.data[j].uid)) continue;
        hash_add_feature(&pcam->trackfeat,tracks.data[j].uid,
                         &tracks.data[j]);
    }
    /* update old track */
    for (i=0;i<tracks.nupd;i++) {
        j=tracks.updtrack[i];
        hash_add_feature(&cam->trackfeat,tracks.data[j].uid,
                         &tracks.data[j]);
    }
    tracetracklen(&tracks);
    return tracks.nnew||tracks.nupd;
}
/* update current camera-frame track data-------------------------------------*/
static int updatecamtrack(const insopt_t *opt,gtime_t time,int pre_id,int cur_ind)
{
    cams_t *cam;
    trace(3,"updatecamtrack: \n");

    /* add current track features */
    if ((cam=getcampointor(cur_ind))&&addcurrenttrack(cam,pre_id)) return 1;
    return 0;
}
/* distort feature point------------------------------------------------------*/
static void distortfeat(const insstate_t *ins,const double *pfn,double *pfnd,
                        double *J)
{
    cam_t camp={0};
    camp.k1=ins->k1; camp.k2=ins->k2;
    camp.p1=ins->p1;
    camp.p2=ins->p2;
    distortradtan(&camp,pfn,pfnd,J);
}
/* calculates measurement variance for a feature point------------------------*/
static void featureR(const trackd *ftrack,const feature *feat,const cams_t *cam,
                     const insstate_t *ins,double *R,double *r)
{
    trace(3,"featureR:\n");
#if 0
    /* other method to determinate feature variance */
#else
    /* default feature variance */
    if (R) {
        R[0]=VAR_FEAT; R[3]=VAR_FEAT;
    }
    if (r) {
        r[0]=VAR_FEAT; r[1]=VAR_FEAT;
    }
#endif
}
/* calculates measurement residual for a feature point------------------------*/
static void featurev(const trackd *ftrack,const feature *feat,const cams_t *cam,
                     const insstate_t *ins,double *v)
{
    double pf[3],pfn[3],pfnd[3];
    double dp[3],uv[2];
    int i;

    trace(3,"featurev:\n");

    /* feature position in c-frame */
    for (i=0;i<3;i++) dp[i]=ftrack->xyz[i]-cam->re[i];
    matmul3v("T",cam->Cce,dp,pf);

    /* distorts feature point */
    for (i=0;i<3;i++) pfn[i]=pf[i]/pf[2];
    distortfeat(ins,pfn,pfnd,NULL);

    /* project into image */
    uv[0]=ins->fx*pfnd[0]+ins->ox;
    uv[1]=ins->fy*pfnd[1]+ins->oy;

    if (v) {
        v[0]=feat->u-uv[0];
        v[1]=feat->v-uv[1];
        trace(3,"v=\n"); tracemat(3,v,1,2,12,6);
    }
}
/* modify the measurement jacobian to ensure observability constrain----------*/
static void mdfyH(const cams_t *cam,const trackd_t *feat,double *Hmf,double *Hmx)
{
    double g[3],u[6],dp[3],W[9],s[6],a,b[4],B[2*6];
    int i,j;

    if (Hmf==NULL||Hmx==NULL) return;

    pregrav(cam->pnull,g);
    matmul("NN",3,1,3,1.0,cam->Cnull,g,0.0,u);

    for (i=0;i<3;i++) {
        dp[i]=feat->xyz[i]-cam->re[i];
    }
    skewsym3(dp,W);

    matmul("NN",3,1,3,1.0,W,g,0.0,u+3);
    matmul("TN",1,1,3,1.0,u,u,0.0,&a);
    for (i=0;i<6;i++) {
        s[i]=u[i]/a;
    }
    matmul("NN",2,1,6,1.0,Hmx,u,0.0,b);
    matmul("NT",2,6,1,1.0,b,s,0.0,B);

    /* for camera attitude */
    for (i=0;i<2;i++) {
        for (j=0;j<6;j++) Hmx[i+j*2]=Hmx[i+j*2]-B[i+2*j];
    }
    /* for feature position */
    for (i=0;i<2;i++) {
        for (j=0;j<3;j++) Hmf[i+j*2]=-Hmx[i+(j+3)*2];
    }
    trace(3,"Hmx(1)=\n");
    tracemat(3,Hmx,2,6,15,4);

    trace(3,"Hmf(1)=\n");
    tracemat(3,Hmf,2,3,15,4);
}
/* calculates measurement jacobians matrix for a feature point----------------*/
static void featureH(const trackd_t *trackf,const cams_t *cam,
                     const insstate_t *ins,double *Hf,double *Hxc,
                     double *Hfo,
                     double *Hkp)
{
    double pf[3],dp[3],pfn[3];
    double J1[9],J2[6]={0},J3[4]={0},J4[9],J5[9],J6[4],W[9],pfnd[3],dpc[3];
    double Jac[9],Jpc[9];
    double T[9],r2;
    int i,j;

    trace(3,"featureH:\n");

    /* feature position in c-frame */
    for (i=0;i<3;i++) dp[i]=trackf->xyz[i]-cam->re[i];
    matmul3v("T",cam->Cce,dp,pf);

    /* distorts feature point */
    for (i=0;i<3;i++) pfn[i]=pf[i]/pf[2]; distortfeat(ins,pfn,pfnd,J6);

    J2[0]=1.0/pf[2]; J2[3]=1.0/pf[2];
    J2[4]=-pf[0]/SQR(pf[2]);
    J2[5]=-pf[1]/SQR(pf[2]);
    J3[0]=ins->fx;
    J3[3]=ins->fy;
    matmul33("NNN",J3,J6,J2,2,2,2,3,T);

    /* Hf: 2*3 for a feature point */
    if (Hf) {
        matcpy(J1,cam->Cce,3,3);
        matmul("NT",2,3,3,1.0,T,J1,0.0,Hf);

        trace(3,"Hf=\n");
        tracemat(3,Hf,2,3,12,6);
    }
    /* Hxc: 2*9 for a camera pose */
    matmul("TN",3,1,3,1.0,cam->Cce,dp,0.0,dpc);
    skewsym3(dpc,J4);
    matmul("NN",2,3,3,1.0,T,J4,0.0,Jac);

    matt(cam->Cce,3,3,J5); for (i=0;i<9;i++) J5[i]=-J5[i];
    matmul("NT",2,3,3,1.0,T,J5,0.0,Jpc);
    if (Hxc) {
        for (i=0;i<2;i++) {
            for (j=0;j<3;j++) Hxc[i+j*2]=Jac[i+(j-0)*2];
            for (j=3;j<6;j++) Hxc[i+j*2]=Jpc[i+(j-3)*2];
        }
        trace(3,"Hxc=\n");
        tracemat(3,Hxc,2,6,12,6);
    }
    /* Hfo: 2*4 for a camera calibration parameters: fx,fy,ox,oy */
    if (Hfo) {
        Hfo[0]=pfnd[0]; Hfo[4]=1.0;
        Hfo[3]=pfnd[1]; Hfo[7]=1.0;

        trace(3,"Hfo=\n"); tracemat(3,Hfo,2,4,12,6);
    }
    /* Hkp: 2*4 for a camera calibration parameters: k1,k2,p1,p2 */
    if (Hkp) {                                  
        r2=SQR(pfn[0])+SQR(pfn[1]);

        Hkp[0]=pfn[0]*r2; Hkp[2]=pfn[0]*SQR(r2);
        Hkp[4]=2.0*pfn[0]*pfn[1];
        Hkp[6]=r2+2.0*SQR(pfn[0]);

        Hkp[1]=pfn[1]*r2; Hkp[3]=pfn[1]*SQR(r2);
        Hkp[5]=r2+2.0*SQR(pfn[1]);
        Hkp[7]=2.0*pfn[0]*pfn[1];

        matcpy(T,Hkp,2,4);
        matmul("NN",2,4,2,1.0,J3,T,0.0,Hkp);

        trace(3,"Hkp=\n");
        tracemat(3,Hkp,2,4,12,6);
    }
    /* modify the measurement jacobian */
    mdfyH(cam,trackf,Hf,Hxc);
}
/* find all camera poses where feature points are observed--------------------*/
static int findallcamera(const trackd_t *feat,const insopt_t *opt,gtime_t time,
                         int *index)
{
    int i,k;
    for (i=0,k=0;i<vofilt.n;i++) {
        if (!hash_find_feature(vofilt.data[i].trackfeat,feat->uid)) continue;
        if (vofilt.data[i].find==-1) continue;
        index[k++]=i;
    }
    return k;
}
/* get camera calibration matrix----------------------------------------------*/
static void getK(const insstate_t *ins,double *K)
{
    /* camera calibration matrix */
    K[0]=ins->fx; K[4]=ins->fy;
    K[6]=ins->ox; K[7]=ins->oy;
    K[8]=1.0;
}
/* compute relative transformation between two frame -------------------------*/
static void reltf(gtime_t ts,gtime_t te,const cams_t *cam1,const cams_t *cam2,
                  double *C,double *t)
{
    double T1[16],T2[16],dT[16];

    rt2tf(cam1->Cce,cam1->re,T1);
    rt2tf(cam2->Cce,cam2->re,T2);
    if (!matinv(T1,4)) {
        matmul("NN",4,4,4,1.0,T1,T2,0.0,dT);
    }
    else {
        seteye(dT,4);
    }
    matinv(dT,4);
    tf2rt(dT,C,t);
}
/* compute camera project matrix----------------------------------------------*/
static int prjmatrix(const double *C,const double *t,const double *K,double *P)
{
    double T[16],M[12];
    int i,j;

    rt2tf(C,t,T);
    for (i=0;i<3;i++) for (j=0;j<4;j++) M[i+j*3]=T[i+j*4];
    matmul("NN",3,4,3,1.0,K,M,0.0,P);
    return 1;
}
/* initial feature point position using chieral method------------------------*/
static int initfeatpos_chieral(const cams_t *cam1,const cams_t *cam2,const insstate_t *ins,
                               const double *uv1,const double *uv2,
                               const double *K,double *pf)
{
    double C[9],t[3],P1[12],P2[12],T1[16],T2[16],dp[3];
    double J[16],V[16],x1[3],x2[3],X[4];
    int j;

    trace(3,"initfeatpos_chieral:\n");

    /* relative camera transformation */
    rt2tf(cam1->Cce,cam1->re,T1);
    rt2tf(cam2->Cce,cam2->re,T2);

    /* camera project matrix */
    if (matinv(T1,4)) seteye(T1,4); tf2rt(T1,C,t); /* first frame */
    prjmatrix(C,t,K,P1);

    if (matinv(T2,4)) seteye(T2,4); tf2rt(T2,C,t); /* second frame */
    prjmatrix(C,t,K,P2);

    for (j=0;j<4;j++) {
        J[0+4*j]=P1[2+3*j]*uv1[0]-P1[0+3*j];
        J[1+4*j]=P1[2+3*j]*uv1[1]-P1[1+3*j];
        J[2+4*j]=P2[2+3*j]*uv2[0]-P2[0+3*j];
        J[3+4*j]=P2[2+3*j]*uv2[1]-P2[1+3*j];
    }
    if (!svd(J,4,4,NULL,NULL,V)) return 0;

    /* return false if this point is at infinity */
    if (fabs(V[3+3*4])<1E-10) {
        trace(2,"feature point is at infinity\n");
        return 0;
    }
    for (j=0;j<4;j++) X[j]=V[j+3*4];
    matmul("NN",3,1,4,1.0,P1,X,0.0,x1);
    matmul("NN",3,1,4,1.0,P2,X,0.0,x2);

    if (x1[2]*X[3]<0.0||x2[2]*X[3]<0.0) return 0;
    for (j=0;j<3;j++) {
        pf[j]=V[j+3*4]/V[3+3*4];
    }
    for (j=0;j<3;j++) {
        dp[j]=pf[j]-cam1->re[j];
    }
    matmul("TN",3,1,3,1.0,cam1->Cce,dp,0.0,pf);
    return 1;
}
/* estimate feature point position in ecef------------------------------------*/
static int estfeatpos(trackd_t *feat,const insopt_t *opt,const insstate_t *ins,const cams_t *cams,
                      gtime_t time,int *index,double *pf)
{
    static double uv1[2],uv2[2];
    double K[9]={0};
    const cams_t *cam1,*cam2;
    int k=0;

    trace(3,"estfeatpos:\n");

    if (!(k=findallcamera(feat,opt,time,index))||(k!=feat->n)) {
        trace(2,"no found camera observed this feature\n");
        return 0;
    }
    if (k<MIN_TRACK_LEN) return 0;
    getK(ins,K);

    /* camera pointer */
    cam1=&cams[index[0]]; cam2=&cams[index[k-1]];

    /* feature point measurement data */
    uv1[0]=feat->data[        0].u; uv1[1]=feat->data[0].v;
    uv2[0]=feat->data[feat->n-1].u;
    uv2[1]=feat->data[feat->n-1].v;

    /* initial feature point position */
    if (!initfeatpos_chieral(cam1,cam2,ins,uv1,uv2,K,pf)) {
        trace(2,"initial feature position fail\n");
        return 0;
    }
    trace(3,"feature position: %4d  %8.4lf  %8.4lf  %8.4lf\n",feat->uid,pf[0],pf[1],pf[2]);
    feat->flag=TRACK_INIT_POS_OK;
    return k;
}
/* get feature point measurement data-----------------------------------------*/
static int getfeatmeas(const feature *feat,const insstate_t *ins,double *pfn)
{
    double uv[3],K[9]={0},pf[3];
    cam_t cam;

    /* camera calibration parameters */
    cam.k1=ins->k1; cam.k2=ins->k2;
    cam.p1=ins->p1; cam.p2=ins->p2;
    getK(ins,K);

    uv[0]=feat->u; uv[1]=feat->v; pfn[2]=uv[2]=1.0;

    if (matinv(K,3)) return 0;
    matmul3v("N",K,uv,pf);

    /* undistort feature point */
    undistortradtan(&cam,pf,pfn,NULL);
    return 1;
}
/* jacobians of measurement wrt. (alpha,beta,rho)-----------------------------*/
static void jacobian(const double *C,const double *t,const double *h,double *J)
{
    J[0]=-C[0]/h[2]+(h[0]/SQR(h[2]))*C[2]; /* for alpha. */
    J[1]=-C[1]/h[2]+(h[1]/SQR(h[2]))*C[2];

    J[2]=-C[3]/h[2]+(h[0]/SQR(h[2]))*C[5]; /* for beta. */
    J[3]=-C[4]/h[2]+(h[1]/SQR(h[2]))*C[5];

    J[4]=-t[0]/h[2]+(h[0]/SQR(h[2]))*t[2]; /* for rho. */
    J[5]=-t[1]/h[2]+(h[1]/SQR(h[2]))*t[2];
}
/* iteration for Gauss Newton optimization------------------------------------*/
static int itergnop(const int* index,int k,trackd_t *feat,const insopt_t *opt,
                    const insstate_t *ins,const double *x,double *Jprev,
                    double *dx,double *Jderiv)
{
    double *E,*W,*v,C[9],t[3],zhat[3],T[3],h[3],J[6],*w;
    double Jnew,*EWE,*Wi,*Ev,*Et;
    cams_t *cam1,*cam2;
    int i,j,l,n;

    trace(3,"itergnop:\n");

    W =zeros(2*k,2*k); E=zeros(2*k,3); v=zeros(2*k,1);
    Wi=zeros(2*k,2*k); Et=zeros(2*k,3);
    w =zeros(2*k,  1);

    for (n=0,j=0,cam1=&vofilt.data[index[0]];j<k;j++) {

        /* relative camera transformation */
        cam2=&vofilt.data[index[j]];
        reltf(cam1->time,cam2->time,cam1,cam2,C,t);

        /* for the weight matrix */
        w[2*n+0]=VAR_FEAT/SQR(ins->fx);
        w[2*n+1]=VAR_FEAT/SQR(ins->fy);

        /* feature observation */
        if (!getfeatmeas(&feat->data[j],ins,zhat)) continue;
        T[0]=x[0]; T[1]=x[1]; T[2]=1.0;

        matmul3v("N",C,T,h);
        for (i=0;i<3;i++) h[i]=h[i]+x[2]*t[i];

        /* form the error vector */
        v[2*n+0]=zhat[0]-h[0]/h[2];
        v[2*n+1]=zhat[1]-h[1]/h[2];

        /* form the jacobians matrix */
        jacobian(C,t,h,J);

        for (i=0;i<2;i++) {
            for (l=0;l<3;l++) Et[(2*n+i)*3+l]=J[i+l*2];
        }
        n++; /* number of valid observation */
    }
    /* check is valid */
    if (!n) {
        free(v); free(W); free(E); free(Wi); free(Et);
        return 0;
    }
    for (j=0;j<n;j++) {
        W[2*j+2*j*(2*n)]=w[2*j];
        W[2*j+1+(2*j+1)*(2*n)]=w[2*j+1];
    }
    matt(Et,2*n,3,E);

    /* calculate the cost function */
    matcpy(Wi,W,2*n,2*n);
    if (matinv(Wi,2*n)) {
        free(v); free(W); free(E); free(Wi); free(Et); return 0;
    }
    matmul33("TNN",v,Wi,v,1,2*n,2*n,1,&Jnew);

    /* solve */
    EWE=mat(3,3); Ev=mat(3,1);
    matmul33("TNN",E,Wi,E,3,2*n,2*n,3,EWE);
    matmul33("TNN",E,Wi,v,3,2*n,2*n,1,Ev);
    if (matinv(EWE,3)) {
        free(v); free(W); free(E); free(Wi);
        free(EWE); free(Ev); free(Et);
        return 0;
    }
    matmul("NN",3,1,3,1.0,EWE,Ev,0.0,dx);
    for (i=0;i<3;i++) dx[i]=-dx[i];

    /* calculate the cost */
    *Jderiv=fabs((0.5*Jnew-*Jprev)/(0.5*Jnew)); *Jprev=0.5*Jnew;

    free( v); free(  W); free( E);
    free(Wi); free(EWE); free(Ev);
    free(Et);
    return 1;
}
/* calculate the position estimate of the feature using Gauss Newton optimization
 * ---------------------------------------------------------------------------*/
static int calcgnposest(trackd_t *feat,const insopt_t *opt,const insstate_t *ins,const cams_t *cams,
                        gtime_t time,int *index)
{
    int i,k=0;
    double Jprev=1E9,Jderiv,x[3],dx[3],pf[3],Jn;

    trace(3,"calcgnposest:\n");

    /* feature position in first camera frame */
    if (!(k=estfeatpos(feat,opt,ins,cams,time,index,pf))) return 0;

    /* re-parameter for feature position */
    for (i=0;i<2;i++) x[i]=pf[i]/pf[2]; x[2]=1.0/pf[2];
    for (i=0;i<MAXITER;i++) {

        /* iteration for optimization */
        if (!itergnop(index,k,feat,opt,ins,x,&Jprev,dx,&Jderiv)) continue;

        if (Jderiv<1E-5) break;
        else {
            for (i=0;i<3;i++) x[i]+=dx[i];
        }
    }
    if (i==MAXITER) {
        trace(2,"optimization fail\n");
        return 0;
    }
    dx[0]=x[0]/x[2];
    dx[1]=x[1]/x[2];
    dx[2]=1.0/x[2];
    matmul("NN",3,1,3,1.0,vofilt.data[index[0]].Cce,dx,0.0,pf);

    /* update feature */
    for (i=0;i<3;i++) {
        feat->xyz[i]=pf[i]+vofilt.data[index[0]].re[i];
    }
    Jn=Jprev/SQR(k);
    if (Jn>MAX_GNCOST_NORM) feat->flag=TRACK_INIT_POS_FAIL;
    return k;
}
/* compute the initial guess of the feature's 3d position using---------------
 * only two views
 * computed feature position in c1 frame. (vec3)
 * --------------------------------------------------------------------------*/
static void initial_guess(const double *Tc1_c2,const double *z1,const double *z2,
                          double *p)
{
    double z1_[3],m[3],R12[9],t12[3];
    double a[2],b[2],depth;

    z1_[0]=z1[0];
    z1_[1]=z1[1];
    z1_[2]=1.0;

    /* construct a least square problem to solve the depth */
    tf2rt(Tc1_c2,R12,t12);
    matmul("NN",3,1,3,1.0,R12,z1_,0.0,m);
    a[0]=m[0]-m[2]*z2[0];
    a[1]=m[1]-m[2]*z2[1];

    b[0]=z2[0]*t12[2]-t12[0];
    b[1]=z2[1]*t12[2]-t12[1];

    /* solve for the depth */
    depth=dot(a,b,2)/dot(a,a,2);
    p[0]=z1[0]*depth;
    p[1]=z1[1]*depth;
    p[2]=depth;
}
/* get feature point for giving time------------------------------------------*/
static struct feature* getfeature(const trackd_t *track,gtime_t time)
{
    int i;
    for (i=0;i<track->n;i++) if (fabs(timediff(track->data[i].time,time))<1E-5) {
            return &track->data[i];
        }
    return NULL;
}
/* compute the cost of the camera observations--------------------------------*/
static double featcost(const double *T0i,double *x,const double *z)
{
    double xyz[3],R[9],t[3],h[3],zh[3];
    double alpha=x[0];
    double beta =x[1];
    double rho  =x[2];

    xyz[0]=alpha;
    xyz[1]=beta;
    xyz[2]=1.0;

    tf2rt(T0i,R,t); matmul("NN",3,1,3,1.0,R,xyz,0.0,h);
    h[0]=h[0]+rho*t[0];
    h[1]=h[1]+rho*t[1];
    h[2]=h[2]+rho*t[2];

    zh[0]=h[0]/h[2];
    zh[1]=h[1]/h[2];
    return SQR(z[0]-zh[0])+SQR(z[1]-zh[1]);
}
/* compute the Jacobian of the camera observation----------------------------*/
static void LMjacobian(const double *R,const double *t,const double *x,
                       const double *z,double *J,double *r,double *w)
{
    double xyz[3],h[3],*W=zeros(3,3);
    double alpha=x[0];
    double beta =x[1];
    double rho  =x[2];
    int i,j;

    xyz[0]=alpha; xyz[1]=beta; xyz[2]=1.0;

    matmul("NN",3,1,3,1.0,R,xyz,0.0,h);
    h[0]=h[0]+rho*t[0];
    h[1]=h[1]+rho*t[1];
    h[2]=h[2]+rho*t[2];

    /* compute the Jacobian */
    for (i=0;i<3;i++) {
        for (j=0;j<2;j++) W[i+j*3]=R[i+j*3];
    }
    for (i=0;i<3;i++) W[i+2*3]=t[i];
    for (i=0;i<2;i++) {
        for (j=0;j<3;j++) {
            J[i+j*2]=W[i+j*3]/h[2]-W[2+j*3]*h[i]/SQR(h[2]);
        }
    }
    /* compute the residual */
    r[0]=h[0]/h[2]-z[0];
    r[1]=h[1]/h[2]-z[1];

    /* compute the weight based on the residual */
    if (norm(r,2)<LMconfig.huber) *w=1.0;
    else {
        *w=LMconfig.huber/(2*norm(r,2));
    }
    free(W);
}
/* initial feature position in ecef using Levenberg-Marquart method ----------
 * the computed 3d position is used to set the position member variable.
 * note the resulted position is in world frame
 * --------------------------------------------------------------------------*/
static int initfeatposLM(trackd_t *feat,const insopt_t *opt,const insstate_t *ins,const cams_t *cams,
                         gtime_t time,int *index)
{
    /* struct type for camera pose */
    struct campose {
        double C[9]; double t[3];
    };
    campose *camposes=(campose*)malloc(sizeof(campose)*vofilt.n);
    campose *camposet=(campose*)malloc(sizeof(campose)*vofilt.n);
    campose *cam0,*cami;
    feature *f0,*f1;
    double T0[16],Ti[16],dT[16],xyz[3],p0[3],p1[3],I[9]={1,0,0,0,1,0,0,0,1};
    double sol[3],nsol[3],total_cost=0.0;
    double fsol[3];
    int k,i,j,valid=1;

    trace(3,"initfeatposLM:\n");

    if (!(k=findallcamera(feat,opt,time,index))||(k!=feat->n)) {
        trace(2,"no found camera observed this feature\n");

        free(camposet);
        free(camposes);
        return 0;
    }
    if (k<MIN_TRACK_LEN) {
        trace(2,"less than min track length\n");

        free(camposet);
        free(camposes);
        return 0;
    }
    for (i=0;i<k;i++) {
        matcpy(camposes[i].C,cams[index[i]].Cce,3,3);
        matcpy(camposes[i].t,cams[index[i]].re,3,1);
    }
    /* all camera poses should be modified such that it takes a vector from
     * the first camera frame in the buffer to this camera frame
     * */
    for (cam0=&camposes[0],i=0;i<k;i++) {
        rt2tf(cam0->C,cam0->t,T0); cami=&camposes[i];
        rt2tf(cami->C,cami->t,Ti);
        matinv(Ti,4);
        matmul("NN",4,4,4,1.0,Ti,T0,0.0,dT);
        tf2rt(dT,camposet[i].C,camposet[i].t);
    }
    /* generate initial guess */
    f0=getfeature(feat,vofilt.data[index[  0]].imgt);
    f1=getfeature(feat,vofilt.data[index[k-1]].imgt);

    if (!f0||!f1) {
        free(camposet);
        free(camposes);
        return 0;
    }
    getfeatmeas(f0,ins,p0);
    getfeatmeas(f1,ins,p1);
    rt2tf(camposet[k-1].C,camposet[k-1].t,Ti);
    initial_guess(Ti,p0,p1,xyz);

    sol[0]=xyz[0]/xyz[2];
    sol[1]=xyz[1]/xyz[2];
    sol[2]=1.0/xyz[2];

    /* Levenberg-Marquart method to solve for the 3d position */
    double lambd=LMconfig.initial_damping,delta_norm=1E9;
    int inner_loop_count=0;
    int outer_loop_count=0;
    int is_cost_reduced=0;

    /* compute the initial cost */
    for (i=0;i<k;i++) {
        rt2tf(camposet[i].C,camposet[i].t,Ti);
        if (!(f1=getfeature(feat,vofilt.data[index[i]].imgt))) continue;
        getfeatmeas(f1,ins,p1);
        total_cost+=featcost(Ti,sol,p1);
    }
    /* outer loop */
    while (outer_loop_count<LMconfig.outloop_maxiter&&delta_norm>LMconfig.ep) {

        double A[9]={0},b[3]={0},J[6],r[2],z[2],w,delta[3],A_[9],new_cost;
        int info;
        for (i=0;i<k;i++) {

            if (!(f1=getfeature(feat,vofilt.data[index[i]].imgt))) continue;
            getfeatmeas(f1,ins,z);
            LMjacobian(camposet[i].C,camposet[i].t,sol,z,J,r,&w);
            matmul("TN",3,3,2,w*w,J,J,1.0,A);
            matmul("TN",3,1,2,1.0,J,r,1.0,b);
        }
        /* inner loop */
        while (inner_loop_count<LMconfig.innloop_maxiter
               &&!is_cost_reduced) {

            for (i=0;i<3;i++) {
                for (j=0;j<3;j++) A_[i+j*3]=A[i+j*3]+lambd*I[i+j*3];
            }
            info=solve("N",A_,b,3,1,delta);
            if (info) break;

            for (i=0;i<3;i++) {
                nsol[i]=sol[i]-delta[i];
            }
            delta_norm=norm(delta,3);
            new_cost=0.0;

            for (i=0;i<k;i++) {
                rt2tf(camposet[i].C,camposet[i].t,Ti);
                if (!(f1=getfeature(feat,vofilt.data[index[i]].imgt))) continue;
                getfeatmeas(f1,ins,p1);
                
                new_cost+=featcost(Ti,nsol,p1);
            }
            if (new_cost<total_cost) {
                is_cost_reduced=1;
                matcpy(sol,nsol,1,3);
                total_cost=new_cost;
                lambd=MAX(lambd/10.0,1E-10);
            }
            else {
                is_cost_reduced=0;
                lambd=MIN(lambd*10.0,1E12);
            }
            inner_loop_count++;
        }
        inner_loop_count=0;
        outer_loop_count++;
    }
    /* covert the feature position from inverse depth
     * representation to its 3d coordinate
     * */
    fsol[0]=sol[0]/sol[2];
    fsol[1]=sol[1]/sol[2];
    fsol[2]=1.0/sol[2];

    /* check if the solution is valid. Make sure the feature
     * is in front of every camera frame observing it
     * */
    for (i=0;i<k;i++) {
        matcpy(xyz,camposet[i].t,1,3);
        matmul("NN",3,1,3,1.0,camposet[i].C,fsol,1.0,xyz);
        if (xyz[2]<0.0) {
            valid=0;
            break;
        }
    }
    /* convert to ecef */
    matcpy(feat->xyz,camposes[0].t,1,3);
    matmul("NN",3,1,3,1.0,camposes[0].C,fsol,1.0,feat->xyz);
    if (valid) {
        feat->flag=TRACK_INIT_POS_OK;
    }
    else {
        feat->flag=TRACK_INIT_POS_FAIL;
    }
    free(camposet);
    free(camposes);
    return k;
}
/* measurement update for a feature point-------------------------------------*/
static int featmeas(trackd_t *trk,const insstate_t *ins,const insopt_t *opt,
                    const cams_t *cams,
                    gtime_t time,double *H,double *v,double *R,double *Hf,
                    int *index,int *ki)
{
    int i,j,k,nv,nx=vofilt.nx,ixs;
    double *Hc,*Hkp,*Hfo,vc[2],Rc[4]={0},*r,*Hfb;
    const cams_t *pcam=NULL;

    static int ifo=xiCfo(opt),nfo=xnCfo(opt);
    static int ikp=xiCkp(opt),nkp=xnCkp(opt);

#if DO_LM_METHOD
    /*  compute feature position*/
    if (!(*ki=k=initfeatposLM(trk,opt,ins,cams,time,index))) return 0;
    if (trk->flag==TRACK_INIT_POS_FAIL) return 0;
#else
    /*  compute feature position*/
    if (!(*ki=k=calcgnposest(trk,opt,ins,cams,time,index))) return 0;
    if (trk->flag==TRACK_INIT_POS_FAIL) return 0;
#endif
    Hc =zeros(2,6); Hkp=zeros(2,4);
    Hfo=zeros(2,4); Hfb=zeros(2,3);
    r=zeros(2*k,1);

    /* for every camera pose */
    for (i=0,nv=0,*ki=0,pcam=cams;i<k&&pcam;i++) {

        /* residual vector */
        featurev(trk,&trk->data[i],&pcam[index[i]],ins,vc);

        /* jacobians matrix */
        featureH(trk,&pcam[index[i]],ins,
                 Hfb,Hc,Hfo,Hkp);

        /* measurement variance */
        featureR(trk,&trk->data[i],&pcam[index[i]],ins,
                 Rc,NULL);

        /* index of camera pose */
        ixs=ins->nx+index[i]*6;
        if (H) {
            for (j=0;j<nx;j++) H[(2*nv+0)*nx+j]=H[(2*nv+1)*nx+j]=0.0;

            /* for fx,fy,ox,oy */
            for (j=0;j<4&&nfo;j++) {
                H[(2*nv+0)*nx+(ifo+j)]=Hfo[2*j+0];
                H[(2*nv+1)*nx+(ifo+j)]=Hfo[2*j+1];
            }
            /* for k1,k2,p1,p2 */
            for (j=0;j<4&&nkp;j++) {
                H[(2*nv+0)*nx+(ikp+j)]=Hkp[2*j+0];
                H[(2*nv+1)*nx+(ikp+j)]=Hkp[2*j+1];
            }
            /* for camera pose */
            for (j=0;j<6;j++) {
                H[(2*nv+0)*nx+(ixs+j)]=Hc[2*j+0]; /* col-major=>row-major*/
                H[(2*nv+1)*nx+(ixs+j)]=Hc[2*j+1];
            }
        }
        if (v) {
            v[2*nv+0]=vc[0]; v[2*nv+1]=vc[1];
        }
        r[2*nv+0]=Rc[0];
        r[2*nv+1]=Rc[3];
        if (Hf) {
            for (j=0;j<3;j++) {
                Hf[(2*nv+0)*3+j]=Hfb[0+j*2]; /* for feature position */
                Hf[(2*nv+1)*3+j]=Hfb[1+j*2];
            }
        }
        if (fabs(v[2*nv+0])>MAX_FEAT_RES||fabs(v[2*nv+1])>MAX_FEAT_RES) continue;
        index[(*ki)++]=index[i];
        nv++;
    }
    if (nv<=0) {
        trace(2,"no feature point measurement\n");
        free(Hc ); free(Hfo);
        free(Hkp); free(r);
        return 0;
    }
    for (i=0;i<2*nv&&R;i++) {
        for (j=0;j<2*nv;j++) if (i==j) R[i+j*2*nv]=r[i]; else R[i+j*2*nv]=0.0;
    }
    trace(3,"v=\n"); tracemat(3,v,2*nv,1,12,6);

    trace(3,"H=\n");
    tracemat(3,H,nx,2*nv,15,4);

    trace(3,"Hf=\n");
    tracemat(3,Hf,3,2*nv,15,4);

    free(Hc ); free(Hfo);
    free(Hkp); free(r);
    return 2*nv;
}
/* use left nullspace of H to `H/v/R'-----------------------------------------*/
static int nulltrick(const double *H,const double *v,const double *R,
                     int nx,int nv,
                     const double *Hf, int m,int n,
                     double *Ho,double *Ro,double *vo,int *p,int *q)
{
    double *N,*Hot;
    int i,j;

    N=mat(MAX(m,n),MAX(m,n)); Hot=mat(nx,nv); if (!null(Hf,n,m,N,&i,&j)) return 0;

    trace(3,"N=\n"); tracemat(3,N,i,j,12,6);

    matmul("TT",j,nx,i,1.0,N,H,0.0,Hot);
    matmul("TN",j,1,nv,1.0,N,v,0.0,vo);
    matmul33("TNN",N,R,N,j,nv,nv,j,Ro);
    matt(Hot,nx,j,Ho);

    *p=i; *q=j; free(N); free(Hot);
    return 1;
}
/* convert the vector part of a quaternion to a full quaternion.--------------
 * this function is useful to convert delta quaternion which is
 * usually a 3x1 vector to a full quaternion.
 * for more details, check Equation (238) and (239) in "Indirect Kalman
 * filter for 3D Attitude Estimation: A Tutorial for quaternion Algebra"
 * ---------------------------------------------------------------------------*/
static void dtheta2q(const double *dtheta,double *q)
{
    double dq[3]; int i;
    dq[0]=dtheta[0]/2.0;
    dq[1]=dtheta[1]/2.0;
    dq[2]=dtheta[2]/2.0;

    if (norm(dq,3)<1.0) {
        matcpy(q,dq,1,3); q[3]=SQRT(1.0-norm(dq,3));
    }
    else {
        matcpy(q,dq,1,3);
        q[3]=1.0;
        for (i=0;i<4;i++) {
            q[i]/=SQRT(1.0+norm(dq,3));
        }
    }
}
/* quaternion multiplication--------------------------------------------------*/
static void quatmul(const double *q1,const double *q2,double *q)
{
    double nq1,nq2,q1_[4],q2_[4],Q[16],nq;
    int i;

    nq1=norm(q1,4); nq2=norm(q2,4);
    for (i=0;i<4;i++) {
        q1_[i]=q1[i]/nq1;
        q2_[i]=q2[i]/nq2;
    }
    Q[0]= q1_[3]; Q[4]= q1_[2]; Q[ 8]=-q1_[1]; Q[12]=q1_[0];
    Q[1]=-q1_[2]; Q[5]= q1_[3]; Q[ 9]= q1_[0]; Q[13]=q1_[1];
    Q[2]= q1_[1]; Q[6]=-q1_[0]; Q[10]= q1_[3]; Q[14]=q1_[2];
    Q[3]=-q1_[0]; Q[7]=-q1_[1]; Q[11]=-q1_[2]; Q[15]=q1_[3];

    matmul("NN",4,1,4,1.0,Q,q2_,0.0,q);
    nq=norm(q,4);
    for (i=0;i<4;i++) {
        q[i]=q[i]/nq;
    }
}
/* correction attitude--------------------------------------------------------*/
static void correctquat(const double *dx,double *Cec)
{
    double dq[4],q[4],qud[4];
    quat_t qc;
    dcm2quat(Cec,&qc);
    q[0]=qc.x;
    q[1]=qc.y;
    q[2]=qc.z;
    q[3]=qc.w;

    dtheta2q(dx,dq);
    quatmul(dq,q,qud);

    qc.w=qud[3];
    qc.z=qud[2];
    qc.y=qud[1];
    qc.x=qud[0];
    quat2dcm(&qc,Cec);
}
/* update camera pose and ins states------------------------------------------*/
static void updatestat(double *dx,insstate_t *ins,cams_t *cams,int ncam,int nx,
                       const insopt_t *opt)
{
    static int ifo=xiCfo(opt),nfo=xnCfo(opt);
    static int ikp=xiCkp(opt),nkp=xnCkp(opt);
    int i,j,ns=ins->nx;
    double Ceb[9],Cec[9];

    trace(3,"dx(%d,%d)=\n",nx,1);
    tracemat(3,dx,nx,1,12,5);

    /* update ins states: [ba bg p v a]*/
    ins->ba[0]+=dx[xiBa(opt)+0];
    ins->ba[1]+=dx[xiBa(opt)+1];
    ins->ba[2]+=dx[xiBa(opt)+2];
    ins->bg[0]+=dx[xiBg(opt)+0];
    ins->bg[1]+=dx[xiBg(opt)+1];
    ins->bg[2]+=dx[xiBg(opt)+2];

    ins->re[0]+=dx[xiP(opt)+0];
    ins->re[1]+=dx[xiP(opt)+1];
    ins->re[2]+=dx[xiP(opt)+2];
    ins->ve[0]+=dx[xiV(opt)+0];
    ins->ve[1]+=dx[xiV(opt)+1];
    ins->ve[2]+=dx[xiV(opt)+2];
    matt(ins->Cbe,3,3,Ceb); correctquat(dx+xiA(opt),Ceb);
    matt(Ceb,3,3,ins->Cbe);

    /* update camera states */
    for (i=0;i<ncam;i++) {
        matt(cams[i].Cce,3,3,Cec); correctquat(dx+ns+6*i,Cec);
        matt(Cec,3,3,cams[i].Cce);
        
        for (j=0;j<3;j++) {
            cams[i].re[j]+=dx[ns+6*i+3+j];
        }
    }
    /* update camera calibration parameters */
    if (nfo) {
        ins->fx=ins->fx-dx[ifo+0]; ins->fy=ins->fy-dx[ifo+1];
        ins->ox=ins->ox-dx[ifo+2];
        ins->oy=ins->oy-dx[ifo+3];
    }
    if (nkp) {
        ins->k1=ins->k1-dx[ikp+0];
        ins->k2=ins->k2-dx[ikp+1];
        ins->p1=ins->p1-dx[ikp+2];
        ins->p2=ins->p2-dx[ikp+3];
    }
    for (i=0;i<nx;i++) dx[i]=0.0;
}
/* outlier detection for feature point measurement----------------------------*/
static int outdetect(const double *P,const double *Ho,const double *Ro,
                     const double *vo,int nv,int nx)
{
    double *HPHT,*HRP,r;
    int i;

    HPHT=mat(nv,nv); HRP=mat(nv,nv);
    matmul33("TNN",Ho,P,Ho,nv,nx,nx,nv,HPHT);
    for (i=0;i<nv*nv;i++) {
        HRP[i]=Ro[i]+HPHT[i];
    }
    if (matinv(HRP,nv)) {
        free(HPHT); free(HRP); return 0;
    }
    matmul33("TNN",vo,HRP,vo,1,nv,nv,1,&r);
    trace(3,"r=\n"); tracemat(3,&r,1,1,12,6);

    free(HPHT); free(HRP);
    return r>chisqr[nv%100];
}
/* x'Fx check-----------------------------------------------------------------*/
static int xFxchk(const cams_t *cam1,const cams_t *cam2,const insopt_t *opt,
                  const feature *f1,const feature *f2,const insstate_t *ins)
{
    return 1;
}
/* epipolar constraint--------------------------------------------------------*/
static int epolarconstraint(const insstate_t *ins, const cams_t *cams, int ncam,
                            const insopt_t *opt,
                            const track_t *track,
                            const int *fidx, int nfidx, gtime_t time,
                            double *ratio)
{
    return 1;
}
/* post-residuals validation--------------------------------------------------*/
static int postvalsol(const insstate_t *ins,const cams_t *cams,int ncam,const insopt_t *opt,
                      const track_t *track,const double *dx,
                      const int *fidx,int nfidx,double thres,gtime_t time)
{
    int flag=0,iba,nba,ibg,nbg;

    iba=xiBa(opt); nba=xnBa(opt);
    ibg=xiBg(opt); nbg=xnBg(opt);

    /* check estimated states */
    if (     dx[  0]==DISFLAG&&norm(dx+  0,3)>15.0*D2R) flag|=1;
    if (nba&&dx[iba]==DISFLAG&&norm(dx+iba,3)>1E5*Mg2M) flag|=1;
    if (nbg&&dx[ibg]==DISFLAG&&norm(dx+ibg,3)>15.0*D2R) flag|=1;
    if (flag) {
        trace(2,"too large estimated state error\n");
        return 0;
    }
    return 1;
}
/* add updated track features to buckets--------------------------------------*/
static int addupd2buckets(track_t *track,const insopt_t *opt,int bw,int bh,bucket_t *buckets)
{
    int i,j,k,u,v;
    for (k=0,i=0;i<track->nupd;i++) {
        if (track->data[i].flag!=TRACK_UPDATED||track->data[i].flag==TRACK_INIT_POS_FAIL) continue;

        j=track->updtrack[i];
        u=(int)floor(track->data[j].data[0].u/BUCKET_WIDTH);
        v=(int)floor(track->data[j].data[0].v/BUCKET_HEIGHT);
        u=u%bw; v=v%bh;

        hash_add_featidx(&buckets[v*bw+u].trackfeat,track->data[j].uid,j,&track->data[j]);
    }
    return k;
}
/* buckets all match features-------------------------------------------------*/
static int bucketmatch(track_t *track,const insopt_t *opt,bucket_t *buckets,
                       int bw,int bh,int useupd,int tlen)
{
    int i,j,u,v;
    for (i=0;i<track->n;i++) {
        if (track->data[i].flag!=TRACK_LOST) continue;
        if (track->data[i].flag==TRACK_INIT_POS_FAIL) continue;
        if (tracks.data[i].flag==TRACK_NEW) continue;
        if (tracks.data[i].n<MIN_TRACK_LEN) continue;
        if (tlen&&tracks.data[i].n!=tlen) continue;

        u=(int)floor(track->data[i].data[0].u/BUCKET_WIDTH);
        v=(int)floor(track->data[i].data[0].v/BUCKET_HEIGHT);
        u=u%bw; v=v%bh;

        hash_add_featidx(&buckets[v*bw+u].trackfeat,track->data[i].uid,i,
                         &track->data[i]);
    }
    if (useupd) addupd2buckets(track,opt,bw,bh,buckets);

    trace(3,"bucket matches:\n");
    for (i=0;i<bh;i++) {
        for (j=0;j<bw;j++) {
            fprintf(stderr,"%2d  ",hash_counts(buckets[i*bh+j].trackfeat));
        }
        fprintf(stderr,"\n");
    }
    return 1;
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
        j=rand()%n; if (list[j]<0) continue;
        sample[ns++]=j; /* add sample index */
        list[j]=-1; /* disable this feature point */
    }
    free(list);
    return ns;
}
/* re-find new feature--------------------------------------------------------*/
static struct hashtable *refindfeature(const int *rlist,int nr, bucket_t *buc,
                                       int nb,int *idx)
{
    struct hashtable *ht=NULL;
    int i,nidx,nf,fcount=0;
    while (true) {
refind:
        if (fcount++>MAXITER) break;
        for (i=0,*idx=rand()%nb;i<nr;i++) if (rlist[i]==*idx) goto refind;
        if (!(nf=hash_counts(buc[*idx].trackfeat))) continue;
        nidx=rand()%nf;
        if (!(ht=hash_index(&buc[*idx].trackfeat,nidx))) continue;
        return ht;
    }
    return NULL;
}
/* get random features from tracking table------------------------------------*/
static int rndbuckets2len(bucket_t *buckets,int nb,bucket *otrk,int nmax,int tlen)
{
    struct hashtable *ht=NULL;
    int i,idx,n,nf,*rl,m=0;

    if (nmax<=0) return 0;
    rl=imat(nmax,1); n=getrand(nb,nmax,rl);
    for (i=0;i<n;i++) {

        if (!(nf=hash_counts(buckets[rl[i]].trackfeat))) {
            if ((ht=refindfeature(rl,n,buckets,nb,&idx))) {
                hash_add_featidx(&otrk->trackfeat,ht->id,ht->index,ht->ptrk);
                rl[i]=idx; m++;
            }
            continue;
        }
        if (!getrand(nf,1,&idx)) continue;
        if (!(ht=hash_index(&buckets[rl[i]].trackfeat,idx))) {
            continue;
        }
        hash_add_featidx(&otrk->trackfeat,ht->id,
                         ht->index,
                         ht->ptrk);
        m++;
    }
    if (m<=0) {
        trace(2,"get random features fail\n");
    }
    free(rl);
    return m;
}
/* generate random bucket matches---------------------------------------------*/
static int rndbuckets(bucket_t *buckets,int nb,bucket *otrk,int nmax)
{
    struct hashtable *current,*tmp;
    int i,n=nmax,m=0;
    for (i=MAX_TRACK_LEN-1;i>=0;i--) {
        m+=rndbuckets2len(&buckets[i*nb],nb,otrk,n-m,i+2);
        if (m>=n) break;
    }
    trace(3,"generate random bucket matches:\n");
    HASH_ITER(hh,otrk->trackfeat,current,tmp) {
        fprintf(stderr,"feature id=%4d len=%2d index=%5d\n",current->id,current->ptrk->n,current->index);
    }
    return hash_counts(otrk->trackfeat);
}
/* check match feature geometric distribution---------------------------------*/
static int chkgeodistribt(const bucket_t *otrk,const insstate_t *ins)
{
    struct hashtable *current,*tmp;
    double K[9],xyz[3],uv[3],*Q,*H,gdop=1E9;
    int n,nv=0,i;

    if ((n=hash_counts(otrk->trackfeat))<3) return 0;
    getK(ins,K); matinv(K,3);

    H=zeros(3,n); Q=zeros(3,3);
    HASH_ITER(hh,otrk->trackfeat,current,tmp) {

        uv[0]=current->ptrk->data[0].u;
        uv[1]=current->ptrk->data[0].v;
        uv[2]=1.0;

        matmul("NN",3,1,3,1.0,K,uv,0.0,xyz);
        for (i=0;i<3;i++) H[3*nv+i]=xyz[i]/norm(xyz,3);
        nv++;
    }
    matmul("NT",3,3,n,1.0,H,H,0.0,Q);
    if (!matinv(Q,3)) {
        gdop=SQRT(Q[0]+Q[4]+Q[8]); /* gdop */

        trace(3,"gdop=%.4lf\n",gdop);
    }
    free(H); free(Q);
    return gdop<MAXGDOP;
}
/* check match feature geometric distribution---------------------------------*/
static int chkgeodistribtidx(const track_t *track,const insstate_t *ins,
                             const int *idx,int nf,double *dop)
{
    double K[9],xyz[3],uv[3],*Q,*H,gdop=1E9;
    int i,j;

    if (nf<3) return 0;
    getK(ins,K); matinv(K,3);
    H=zeros(3,nf); Q=zeros(3,3);
    for (i=0;i<nf;i++) {
        uv[0]=track->data[idx[i]].data[0].u;
        uv[1]=track->data[idx[i]].data[0].v;
        uv[2]=1.0;

        matmul("NN",3,1,3,1.0,K,uv,0.0,xyz);
        for (j=0;j<3;j++) H[3*i+j]=xyz[j]/norm(xyz,3);
    }
    matmul("NT",3,3,nf,1.0,H,H,0.0,Q);
    if (!matinv(Q,3)) {
        gdop=SQRT(Q[0]+Q[4]+Q[8]); /* gdop */

        trace(3,"gdop=%.4lf\n",gdop);
        if (dop) *dop=gdop;
    }
    free(H); free(Q);
    return gdop<MAXVOGDOP;
}
/* get random feature for processing------------------------------------------*/
static int randfeats(track_t *track,const insopt_t *opt,const insstate_t *ins,
                     int bw,int bh,int **idx,int nfeats,int chkgeo,bucket_t *buckets)
{
    bucket_t otrk={0};
    struct hashtable *ht=NULL;
    int i,n=0;

    if ((n=rndbuckets(buckets,bw*bh,&otrk,nfeats))<3) goto exit;
    if (chkgeo&&!chkgeodistribt(&otrk,ins)) return 0;

    *idx=imat(1,n);

    for (i=0;i<n;i++) {
        ht=hash_index(&otrk.trackfeat,i);
        (*idx)[i]=ht->index;
    }
exit:
    hash_destroy(&otrk.trackfeat);
    return n;
}
/* QR decomposition-----------------------------------------------------------*/
static int doqrdcmp(const double *H,const double *v,const double *R,int nv,int nx,
                    double *Hq,double *vq,double *Rq,int *nq)
{
    static double r0=re_norm(0.96);
    double *Qo,*Ro,*Th,*Q1,g,a;
    int i,j,flag;

    Qo=mat(nv,nv);
    Ro=mat(nv,nx);

    if (!qr(H,nv,nx,Qo,Ro,1)) {
        trace(2,"QR decomposition fail\n");
        free(Qo);
        free(Ro);
        return 0;
    }
    for (i=0;i<nv;i++) {
        for (flag=1,j=0;j<nx;j++) if (fabs(Ro[i+j*nv])>1E-20) {flag=0; break;}
        if (flag) break;
    }
    *nq=i+1;
    Th=mat(*nq,nx);
    Q1=mat(nv,*nq);

    for (i=0;i<nv;i++) {
        for (j=0;j<*nq;j++) Q1[i+j*nv]=Qo[i+j*nv];
    }
    for (i=0;i<*nq;i++) {
        for (j=0;j<nx;j++) Th[i+j**nq]=Ro[i+j*nv];
    }
    if (Hq) {
        for (i=0;i<*nq;i++) {
            for (j=0;j<nx;j++) {
                Hq[i*nx+j]=Th[i+j**nq];
            }
        }
    }
    if (Rq) {
        matmul33("TNN",Q1,R,Q1,*nq,nv,nv,*nq,Rq);
    }
    if (vq) {
        matmul("TN",*nq,1,nv,1.0,Q1,v,0.0,vq);
    }
    g=stds(vq,*nq-1); a=avg(vq,*nq);
    for (i=0;i<*nq;i++) {
        if (fabs(vq[i]-a)/g>r0) Rq[i+i**nq]=SQR(100.0);
    }
    free(Qo);
    free(Ro); free(Th); free(Q1);
    return 1;
}
/* iteration for feature point measurement update-----------------------------*/
static int iterfeatupd(const insopt_t *opt,insstate_t *ins,cams_t *cams,int ncam,
                       double *Pc,gtime_t time,
                       const int *fidx,int fnidx,
                       double *dop,double *r2,double *r3)
{
    double *H,*R,*v,*Ho,*Ro,*vo,*Hf,*x;
    double *Hk,*Rk,*vk,*rk;
    double *Hq,*Rq,*vq;
    double *Hp,*Rp,*vp;
    int i,j,nv,nx=vofilt.nx,n=vofilt.n,index[MAX_TRACK_LEN],ki=0,k,nvp,*idnx;
    int p,q,ok=0,nf,nq;

    trace(3,"iterfeatupd:\n");

    if (fnidx<=0||ncam==0) return 0;

    H=zeros(n*2,nx); R=zeros(n*2,n*2);
    Hf=zeros(n*2,3); v=zeros(n*2,1);

    Ho=mat(n*2,nx); Ro=mat(n*2,n*2);
    vo=mat(n*2,1);
    x=zeros(nx,1);

    Hk=mat(n*2*fnidx,nx); Rk=zeros(n*2*fnidx,n*2*fnidx);
    vk=mat(n*2*fnidx,1);
    rk=mat(n*2*fnidx,1);

    idnx=imat(1,fnidx);

    /* for each feature measurement data */
    for (nv=0,nf=0,k=0;k<fnidx;k++) {
        if (!(nvp=featmeas(&tracks.data[fidx[k]],ins,opt,cams,time,H,v,R,Hf,index,&ki))) continue;

        /* null-space trick */
        if (!nulltrick(H,v,R,nx,nvp,Hf,nvp,3,Ho,Ro,vo,&p,&q)) continue;

        /* stacked matrixs */
        for (i=0;i<q;i++) {
            for (j=0;j<nx;j++) {Hk[(nv+i)*nx+j]=Ho[i*nx+j];}
            vk[nv+i]=vo[i];
            rk[nv+i]=Ro[i*q+i];
        }
        nv=nv+q;
        idnx[nf++]=fidx[k];
    }
    if (nv>3) {
        for (i=0;i<nv;i++) Rk[i+i*nv]=rk[i];

        trace(3,"Hk=\n"); tracemat(3,Hk,nx,nv,12,6);
        trace(3,"vk=\n"); tracemat(3,vk,nv, 1,12,6);
    }
    else {
        trace(2,"no feature measurement data\n");
        free(H ); free(v ); free(R ); free(x);
        free(Hf); free(Ho); free(vo);
        free(Ro); free(Hk); free(Rk);
        free(vk); free(idnx);
        return 0;
    }
    /* match feature geometric distribution */
    if (!chkgeodistribtidx(&tracks,ins,idnx,nf,dop)) {

        trace(2,"geometric distribution check fail\n");
        free(H ); free(v ); free(R ); free(x);
        free(Hf); free(Ho); free(vo);
        free(Ro); free(Hk); free(Rk);
        free(vk); free(idnx);
        return 0;
    }
#if OUT_DETECT
    /* outlier detection */
    if (nv&&outdetect(Pc,Hk,Rk,vk,nv,nx)) {
        trace(2,"outlier detection fail\n");

        free(H ); free(v ); free(R ); free(x);
        free(Hf); free(Ho); free(vo);
        free(Ro); free(Hk); free(Rk);
        free(vk); free(idnx);
        return 0;
    }
#endif

#if DO_QR_DCMP
    Hq=mat(nv,nx); vq=mat(nv,1);
    Rq=mat(nv,nv);

    /* QR decomposition */
    if (nv>nx&&doqrdcmp(Hk,vk,Rk,nv,nx,Hq,vq,Rq,&nq)) {
        Hp=Hq; Rp=Rq; vp=vq;
    }
    else {
        Hp=Hk; Rp=Rk;
        vp=vk; nq=nv;
    }
    /* EKF filter */
    if (filter(x,Pc,Hp,vp,Rp,nx,nq)) {
        trace(2,"EKF filter fail\n");

        free(H ); free(v ); free(R ); free(x);
        free(Hf); free(Ho); free(vo);
        free(Ro); free(Hk); free(Rk);
        free(vk);
        free(Hq); free(Rq); free(vq);
        return 0;
    }
    free(Hq);
    free(Rq); free(vq);
#else
    /* EKF filter */
    if (filter(x,Pc,Hk,vk,Rk,nx,nv)) {
        trace(2,"EKF filter fail\n");

        free(H ); free(v ); free(R ); free(x);
        free(Hf); free(Ho); free(vo);
        free(Ro); free(Hk); free(Rk);
        free(vk); free(idnx);
        return 0;
    }
#endif
    updatestat(x,ins,cams,ncam,nx,opt); /* update all states */

    /* post residual validation */
    if (!postvalsol(NULL,NULL,ncam,opt,NULL,x,NULL,0,4.0,time)) {
        trace(2,"post residual validation fail\n");

        free(H ); free(v ); free(R ); free(x);
        free(Hf); free(Ho); free(vo);
        free(Ro); free(Hk); free(Rk);
        free(vk); free(idnx);
        return 0;
    }
#if EPIPOLAR_VALIDATION
    /* epipolar constraint check */
    epolarconstraint(ins,cams,ncam,opt,&tracks,tracks.losttrack,tracks.nlost,time,r2);

    /* solution validation  */
    if (*r2<RATIO_THRESHOLD) {
        ok=0; /* validation fail */
    }
    else ok=1; /* validation ok */
#else
    ok=1;
#endif
    free(H ); free(v ); free(R ); free(x);
    free(Hf); free(Ho); free(vo);
    free(Ro); free(Hk); free(Rk);
    free(vk); free(idnx);
    return ok;
}
/* remove feature from track datas--------------------------------------------*/
static void rmfeature(trackd_t *track,vofilt_t *vo)
{
    int j;
    for (j=0;j<vo->n;j++) {

        if (!hash_find_feature(vo->data[j].trackfeat,track->uid)) continue;
        hash_delete_feature(&vo->data[j].trackfeat,track->uid);
    }
    track->flag=TRACK_FILTER;
    track->last_idx=-1;
}
/* reset ins and camera all states--------------------------------------------*/
static void resetstates(insstate_t *insp,cams_t *cams,double *Pc,const insstate_t *ins,
                        const vofilt_t *vo)
{
    memcpy(cams,vo->data,sizeof(cams_t)*vo->n);
    memcpy(insp,ins,sizeof(insstate_t));
    matcpy(Pc,vo->Px,vo->nx,vo->nx);
}
/* corrections for ins and camera states--------------------------------------*/
static void updinscamstates(const insstate_t *insp,const cams_t *cams,const double *Pc,
                            const double r1,const double r2,const double r3,
                            insstate_t *ins,vofilt_t *vo)
{
    int i,j;
    memcpy(ins,insp,sizeof(insstate_t));

    for (i=0;i<vo->n;i++) {
        matcpy(vo->data[i].Cce,cams[i].Cce,3,3);
        matcpy(vo->data[i].re,cams[i].re,1,3);

        vo->data[i].find=cams[i].find;
        vo->data[i].imgt=cams[i].imgt;
        vo->data[i].time=cams[i].time;
    }
    matcpy(vo->Px,Pc,vo->nx,vo->nx);
    for (i=0;i<ins->nx;i++) {
        for (j=0;j<ins->nx;j++) ins->P[i+j*ins->nx]=Pc[i+j*vo->nx];
    }
    vo->r1=r1; /* validation ratio */
    vo->r2=r2;
    vo->r3=r3;
}
/* feature point measurement update-------------------------------------------*/
static int updatefeatmeas(const insopt_t *opt,insstate_t *ins,vofilt_t *vo,gtime_t time)
{
    static int bw=opt->voopt.match.img_w/BUCKET_WIDTH;
    static int bh=opt->voopt.match.img_h/BUCKET_HEIGHT;
    insstate_t insp,insb;
    cams_t cams[MAX_TRACK_LEN],camb[MAX_TRACK_LEN];
    bucket_t *buckets,*bucketv;
    int i,*fidx=NULL,fnidx,flag,count=0,less=0,upd=0;
    double *Pc,*Pb,r1,r2,r3,r1b=0,r2b=0,r3b=0;

    trace(3,"updatefeatmeas:\n");

    if (tracks.nexceed==0&&tracks.nlost==0) return 0;
    Pc=mat(vo->nx,vo->nx);
    Pb=mat(vo->nx,vo->nx);
    resetstates(&insp,cams,Pc,ins,vo);
    resetstates(&insb,camb,Pb,ins,vo);

    buckets=(bucket*)calloc(sizeof(bucket),(size_t)bw*bh*MAX_TRACK_LEN);
    bucketv=(bucket*)calloc(sizeof(bucket),(size_t)bw*bh);

    for (i=0;i<MAX_TRACK_LEN;i++) bucketmatch(&tracks,opt,&buckets[bw*bh*i],bw,bh,0,i+2);
    bucketmatch(&tracks,opt,bucketv,bw,bh,1,0);

    /* update exceed max track length feature */
    if (iterfeatupd(opt,&insp,cams,vo->n,Pc,time,tracks.exceedtrack,tracks.nexceed,
                    &r1,&r2,&r3)) {
        updinscamstates(&insp,cams,Pc,r1,r2,r3,ins,vo);
    }
    if (tracks.nlost<50) {
        less=1; /* less feature to update */
    }
#if BUCKET_ROBUST

    /* update track lost feature */
    while (true) {
        if (count++==MAX_RANSAC_ITER) break;

        /* random feature points */
        if (!less&&!(fnidx=randfeats(&tracks,opt,&insp,bw,bh,&fidx,48,1,buckets))) continue;

        /* iteration update */
        if (!(flag=iterfeatupd(opt,&insp,cams,vo->n,Pc,time,less?tracks.losttrack:fidx,less?tracks.nlost:fnidx,&r1,&r2,&r3))) {
            trace(2,"iteration update fail\n");
        }
        if (fidx) free(fidx); fidx=NULL;
        if (flag&&r2>r2b&&r3>r3b) {

            memcpy(camb,cams,sizeof(cams_t)*vo->n);
            memcpy(&insb,&insp,sizeof(insstate_t));
            matcpy(Pb,Pc,vo->nx,vo->nx);
            r1b=r1;
            r2b=r2;
            r3b=r3; upd=1; /* update ok */
        }
        if (less) break;

        /* reset states */
        resetstates(&insp,cams,Pc,
                    ins,vo);
    }
#else
    upd=iterfeatupd(opt,&insb,camb,vo->n,Pb,time,tracks.losttrack,tracks.nlost,&r1,&r2,&r3);
#endif
    /* remove exceed max track length/track lost features */
    for (i=0;i<tracks.nexceed;i++) {
        rmfeature(&tracks.data[tracks.exceedtrack[i]],vo);
    }
    for (i=0;i<tracks.nlost;i++) {
        rmfeature(&tracks.data[tracks.losttrack[i]],vo);
    }
    /* remove feature of track length less than min-len */
    for (i=0;i<tracks.nless;i++) {
        rmfeature(&tracks.data[tracks.lesstrack[i]],vo);
    }
    /* update all states */
    if (upd) {
        updinscamstates(&insb,camb,Pb,r1b,r2b,r3b,ins,vo);
    }
    /* free all buckets */
    for (i=0;i<bw*bh*MAX_TRACK_LEN;i++) {
        hash_destroy(&buckets[i].trackfeat);
    }
    for (i=0;i<bw*bh;i++) {
        hash_destroy(&bucketv[i].trackfeat);
    }
    free(Pc); free(Pb);
    free(buckets);
    free(bucketv); return upd;
}
/* another method of update all feature points -------------------------------*/
static int updatefeatalltwo(const insopt_t *opt,const img_t *img)
{
    cams_t *pcam=NULL;
    int i,j;

    trace(3,"updatefeatalltwo:\n");

    if (img==NULL||img->feat==NULL||matchfeats(&matchs,img)<=0) {
        trace(2,"no feature points measurement data\n");
        return 0;
    }
    if (!(pcam=getcampointor(img->id))) return 0;
    if (!updatetrack(opt,img)) return 0;

    /* update new tracks */
    for (i=0;i<tracks.nnew;i++) {
        j=tracks.newtrack[i];
        hash_add_feature(&pcam->trackfeat,tracks.data[j].uid,&tracks.data[j]);
    }
    /* update old tracks */
    for (i=0;i<tracks.nupd;i++) {
        j=tracks.updtrack[i];
        hash_add_feature(&pcam->trackfeat,tracks.data[j].uid,
                         &tracks.data[j]);
    }
    trace(3,"add new tracks: %d\n",tracks.nnew);
    trace(3,"update tracks: %d\n",tracks.nupd);

    /* get number of tracking lost */
    for (i=0,tracks.nlost=0;i<tracks.n;i++) {
        if (tracks.data[i].flag==TRACK_LOST&&tracks.data[i].last_idx!=-1&&tracks.data[i].n>=MIN_TRACK_LEN) {
            tracks.losttrack[tracks.nlost++]=i;
        }
    }
    trace(3,"track lost: %d\n",tracks.nlost);

    /* exceed max track length */
    for (i=0,tracks.nexceed=0;i<tracks.n;i++) {
        if (tracks.data[i].n>=MAX_TRACK_LEN&&tracks.data[i].flag==TRACK_UPDATED) {
            tracks.exceedtrack[tracks.nexceed++]=i;
        }
    }
    trace(3,"track exceed: %d\n",tracks.nexceed);

    for (tracks.nless=0,i=0;i<tracks.n;i++) {
        if (tracks.data[i].flag==TRACK_LOST&&tracks.data[i].n<MIN_TRACK_LEN) tracks.lesstrack[tracks.nless++]=i;
    }
    trace(3,"track less: %d\n",tracks.nless);

    tracetracklen(&tracks);
    return tracks.nnew||tracks.nupd;
}
/* update all feature points--------------------------------------------------*/
static int updatefeatallone(const insopt_t *opt,const img_t *img)
{
    trace(3,"updatefeatallone:\n");

    /* match feature points */
    if (img==NULL||matchfeats(&matchs,img)<=0) {

        trace(2,"match feature points fail\n");
        return 0;
    }
    /* update all track data */
    if (!updatetrack(opt,img)) return 0;
    if (!updatecamtrack(opt,img->time,matchs.Ip.id,img->id)) return 0;
    return 1;
}
/* update.--------------------------------------------------------------------*/
static int updateall(insstate_t *ins,const insopt_t *opt,const img_t *img)
{
    trace(3,"updateall:\n");

    /* augment states with a new camera pose */
    augstates(&vofilt,ins,opt,img,img==NULL?-1:img->id);

#if !DO_ANOTHER_UPD
    if (!updatefeatallone(opt,img)) return 0;
#else
    if (!updatefeatalltwo(opt,img)) return 0; 
#endif

#if DO_MONO_VO_EST
    /* monocular visual odometry estimate */
    ins->vo.status=estmonort(&opt->voopt,&matchs.mp_bucket,ins->vo.dT,NULL);
#endif
    /* update camera/ins states */
    if (updatefeatmeas(opt,ins,&vofilt,img->time)) {
        
        trace(3,"feature measurement update ok\n");
        ins->stat=SOLQ_VO;
    }
    else {
        trace(3,"feature update fail\n");
        return 0;
    }
    /* remove redundant camera */
    rmcamera(&vofilt,ins,opt);
    return 1;
}
/* using visual odometry to aid ins/gnss pose estimating-----------------------
 * args:    insopt_t *opt   I   ins options
 *          insstate_t *ins IO  ins states
 *          imud_t *imu     I   imu measurement data
 *          img_t *img      I   image measurement data
 *          int flag        I   flag-0: propagate ins-camera pose covariance
 *                              flag-1: feature points measurement update
 * return: status (1: ok, 0: fail)
 * ----------------------------------------------------------------------------*/
extern int voigposmsckf(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                        const img_t *img,int flag)
{
    trace(3,"voigposmsckf: time=%s\n",time_str(imu->time,4));

    switch (flag) {
        case 0: return propagate(&vofilt,ins,opt);
        case 1: return propagate(&vofilt,ins,opt)&&updateall(ins,opt,img);
        default: {
            trace(2,"not support mode\n");
        }
    }
    if (inss==NULL) inss=ins;
    if (opts==NULL) opts=opt;
    return 0;
}
