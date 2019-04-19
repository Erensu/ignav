/*--------------------------------------------------------------------------------
 * ins-gnss-vo-sim.cc : ins-gnss-vo multisensor data simulator
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

/* type definitions ----------------------------------------------------------*/
typedef struct {               /* true motion data type */
    gtime_t time;              /* motion time stamp */
    double re[3],ve[3],Cbe[9]; /* position/velocity/attitude matrix */
} motion_t;

typedef struct {
    int n,nmax;                /* number and max number of motion data */
    motion_t *data;            /* all motion data */
} path_t;

typedef struct {               /* feature point measurement data type */
    double u,v;                /* image coordinate */
    double re[3],rc[3];        /* feature position in ecef/camera-frame */
    double pc[3];              /* camera position in ecef */
    int uid;                   /* a unique identifier of this feature point */
    UT_hash_handle hh;         /* makes this structure hashable */
} feat_t;

typedef struct {               /* frame data type*/
    int uid;                   /* id of this frame data */
    double Cce[9],re[3];       /* ture camera pose in ecef */
    gtime_t time;              /* frame time stamp */
    feat_t *feat;              /* feature point in current frame */
} frame_t;

/* constants-------------------------------------------------------------------*/
#define VAR_POS                     SQR(0.01) /* variance of GPS position measurement data */
#define VAR_FEAT                    SQR(3.0)  /* variance of feature point measurement data */
#define MIN_TRACK_LEN               5         /* min length of tracking data */
#define MAX_TRACK_LEN               20        /* max length of tracking data */
#define MIN_RANGE                   -50       /* min range of generate feature point */
#define MAX_RANGE                   50        /* max range of generate feature point */
#define MIN_DEPTH                   0.5       /* min depth of generate feature point */
#define MAX_DEPTH                   30        /* max depth of generate feature point */
#define MAX_NUM_FEAT                500       /* max number of generate feature point */
#define IMG_W                       500       /* image width in pixels */
#define IMG_H                       500       /* image height in pixels */

/* global variables------------------------------------------------------------*/
static int id_seed=1;                         /* generate a new feature id */
static frame_t pf={0},cf={0};                 /* precious frame and current frame data */

/* imu body position transform to gps antenna---------------------------------*/
static void generategps(const double *Cbe,const double *re,const double *lever,
                        double *rr)
{
    double T[3];
    int i;
    matmul3v("N",Cbe,lever,T);
    for (i=0;i<3;i++) {
        rr[i]=re[i]+T[i];
    }
}
/* generates a uniform distributed random number between min and max --------*/
static double getuniform(double min,double max)
{
    return 1.0*rand()/RAND_MAX*(max-min)+min;
}
/* creates gaussian distributed random numbers (Box-Müller method)-----------*/
static double getgaussian(double std)
{
    if (std<0.0) std=-std;
    double x1,x2,w,y1;
    do {
        x1=getuniform(-1.0,1.0); x2=getuniform(-1.0,1.0); w=x1*x1+x2*x2;
    } while (w>=1.0);
    w=sqrt((-2.0*log(w))/w); y1=x1*w; return std*y1;
}
/* add new feature and image data to track------------------------------------*/
static int addmotion(path_t *path,const motion_t *data)
{
    motion_t *obs_data;
    size_t size;
    if (path->nmax<=path->n) {
        if (path->nmax<=0) path->nmax=1024; else path->nmax+=126;
        size=sizeof(motion_t)*path->nmax;
        if (!(obs_data=(motion_t*)realloc(path->data,size))) {
            free(path->data);
            path->data=NULL; path->n=path->nmax=0;
            return -1;
        }
        else {
            path->data=obs_data;
        }
    }
    path->data[path->n++]=*data;
    return 1;
}
/* find feature in hash table-------------------------------------------------*/
static feat_t* hash_findfeat(const feat_t* feats,const int uid)
{
    feat_t *s=NULL;
    HASH_FIND_INT(feats,&uid,s);  /* s: output pointer */
    return s;
}
/* add a feature to a frame---------------------------------------------------*/
static void hash_addfeat(feat_t** feats,int uid,double u,double v,double *re,
                         double *rc,const double *pc)
{
    feat_t *s=NULL;
    int i;
    HASH_FIND_INT(*feats,&uid,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(feat_t*)malloc(sizeof(feat_t));
        s->uid=uid;
        s->u=u;
        s->v=v;
        for (i=0;i<3;i++) s->re[i]=re[i];
        for (i=0;i<3;i++) s->rc[i]=rc[i];
        for (i=0;i<3;i++) s->pc[i]=pc[i];
        HASH_ADD_INT(*feats,uid,s);  /* id: name of key field */
    }
}
/* destory a frame------------------------------------------------------------*/
static void hash_rmframe(feat_t** feats)
{
    feat_t *current,*tmp;

    HASH_ITER(hh,*feats,current,tmp) {
        HASH_DEL(*feats,current); /* delete; users advances to next */
        free(current);            /* optional- if you want to free  */
    }
}
/* random a feature point-----------------------------------------------------*/
static int randfeatpos(gtime_t time,const cam_t *cam,double *fpos,double *uv)
{
    double xyz[3];
    int i;

    fpos[0]=getuniform(MIN_RANGE,MAX_RANGE);
    fpos[1]=getuniform(MIN_RANGE,MAX_RANGE);
    fpos[2]=getuniform(MIN_DEPTH,MAX_DEPTH);

    for (i=0;i<3;i++) {
        xyz[i]=fpos[i]/fpos[2];
    }
#if 1
    uv[0]=xyz[0]*cam->K[0]+cam->K[6]+getgaussian(4);
    uv[1]=xyz[1]*cam->K[4]+cam->K[7]+getgaussian(4);
#else
    uv[0]=xyz[0]*cam->K[0]+cam->K[6];
    uv[1]=xyz[1]*cam->K[4]+cam->K[7];
#endif
    if (uv[0]<=0.0||uv[1]<=0.0) return 0;
    if (uv[0]>=IMG_W||uv[1]>=IMG_H) {
        return 0;
    }
    return 1;
}
/* generate uid of tracking lost----------------------------------------------*/
static int generatelostfeat(feat_t *pf,int nlost,int *idslost,int *idsupd,
                            int *k)
{
    feat_t *current,*tmp;
    int n=0,i=0,*idx,j,m=0,*flag;

    n=HASH_COUNT(pf); if (n<=0) return 0;
    idx=imat(n,1);
    flag=imat(n,1);

    HASH_ITER(hh,pf,current,tmp) {
        idx[i]=current->uid;
        flag[i++]=1;
    }
    for (i=0;m<(nlost>n?n/2:nlost);i++) {
        j=rand()%n;

        if (m>=nlost) break;
        if (flag[j]==0) continue;
        
        idslost[m++]=idx[j];
        flag[j]=0;
    }
    for (*k=0,i=0;i<n;i++) {
        if (flag[i]==1) idsupd[(*k)++]=idx[i];
    }
    free(idx); return m;
}
/* drop feature out of image range--------------------------------------------*/
static int dropfeat(feat_t **feat)
{
    feat_t *current,*tmp;
    int ndrop=0;
    HASH_ITER(hh,*feat,current,tmp) {

        if (current->u>0.0&&current->u<IMG_W&&
            current->v>0.0&&current->v<IMG_H) {
            continue;
        }
        HASH_DEL(*feat,current); /* delete; users advances to next */
        free(current);           /* optional- if you want to free  */
        ndrop++;
    }
    return ndrop;
}
/* generate new feature point-------------------------------------------------*/
static int generatenewfeat(feat_t **pf,const double *Cce,const double *re,
                           const cam_t *cam,gtime_t time,
                           int nmax,int **idsnew)
{
    double fpos[3],uv[2],fre[3];
    int n=0;

    if (nmax==0) return 0;
    if (idsnew) *idsnew=imat(nmax,1);

    while (true) {
        if (!randfeatpos(time,cam,fpos,uv)) continue;
        if (n>=nmax) break;

        matmul("NN",3,1,3,1.0,Cce,fpos,0.0,fre);

        fre[0]+=re[0];
        fre[1]+=re[1];
        fre[2]+=re[2];

        hash_addfeat(pf,id_seed++,uv[0],uv[1],fre,fpos,re);
        if (idsnew) {
            (*idsnew)[n]=id_seed;
        }
        n++;
    }
    return n;
}
/* copy hash table from another hash table------------------------------------*/
static void hashcopy(frame_t *pf,frame_t *cf)
{
    hash_rmframe(&pf->feat); pf->feat=NULL;
    feat_t *current,*tmp;
    HASH_ITER(hh,cf->feat,current,tmp) {
        hash_addfeat(&pf->feat,current->uid,current->u,current->v,current->re,current->rc,current->pc);
    }
    matcpy(pf->Cce,cf->Cce,3,3);
    matcpy(pf->re ,cf->re ,3,1);

    pf->time=cf->time;
    pf->uid =cf->uid;
    hash_rmframe(&cf->feat); cf->feat=NULL;
}
/* generate image feature point-----------------------------------------------*/
static int generatefeat(const double *Cce,const double *re,const cam_t *cam,
                        gtime_t time,int pre_id,int cur_id,
                        frame_t *pf,frame_t *cf)
{
    feat_t *pfeat;
    int i,*idslost,nf,nlost,nupd,nnew,ndrop;
    int *idsupd,*idsnew=NULL;
    double dp[3],rc[3],u,v;

    matcpy(cf->Cce,Cce,3,3);
    matcpy(cf->re,re,3,1);
    cf->uid=cur_id;
    cf->time=time;

    if (!(nf=HASH_COUNT(pf->feat))) {

        generatenewfeat(&pf->feat,cf->Cce,cf->re,cam,time,MAX_NUM_FEAT,&idsnew);
        return 0;
    }
    idslost=imat(nf,1);
    idsupd=imat(nf,1);

    nlost=generatelostfeat(pf->feat,100,idslost,idsupd,&nupd);
    if (nlost<=0&&nupd) {

        free(idslost); free(idsupd);
        return 0;
    }
    nnew=generatenewfeat(&cf->feat,cf->Cce,cf->re,cam,time,100,&idsnew);
    if (nnew<=0) {

        free(idslost); free(idsupd);
        if (idsnew) free(idsnew);
        return 0;
    }
    for (i=0;i<nupd;i++) {
        if (!(pfeat=hash_findfeat(pf->feat,idsupd[i]))) continue;
        dp[0]=pfeat->re[0]-cf->re[0];
        dp[1]=pfeat->re[1]-cf->re[1];
        dp[2]=pfeat->re[2]-cf->re[2];

        matmul("TN",3,1,3,1.0,cf->Cce,dp,0.0,rc);
        u=rc[0]/rc[2]*cam->K[0]+cam->K[6];
        v=rc[1]/rc[2]*cam->K[4]+cam->K[7];

        hash_addfeat(&cf->feat,pfeat->uid,u,v,pfeat->re,pfeat->rc,pfeat->pc);
    }
    ndrop=dropfeat(&cf->feat);
    generatenewfeat(&cf->feat,cf->Cce,cf->re,cam,time,ndrop,NULL);

    trace(3,"currenr frame: uid=%4d  nf=%4d  \n",cf->uid,HASH_COUNT(cf->feat));
    free(idslost); free(idsupd);
    free(idsnew);
    return 1;
}
/* generate frame data--------------------------------------------------------*/
static int generateframe(const double *Cbe,const double *re,int fid,const cam_t *cam,
                         const insopt_t *insopt,gtime_t time,frame_t **of)
{
    static int pre_fid=-1;
    double Cce[9],rc[3],Cbc[9],lbc[3],T[3];
    int i;

    matcpy(lbc,insopt->voopt.lbc,1,3);
    rpy2dcm(insopt->voopt.ebc,Cbc);

    matmul("NT",3,3,3,1.0,Cbe,Cbc,0.0,Cce);
    matmul("NN",3,1,3,1.0,Cbe,lbc,0.0,T);
    for (i=0;i<3;i++) {
        rc[i]=re[i]+T[i];
    }
    if (pre_fid==-1) {
        matcpy(pf.Cce,Cce,3,3);
        matcpy(pf.re ,rc ,3,1);
        pf.time=time;
        pf.uid =fid;

        generatenewfeat(&pf.feat,Cce,rc,cam,time,500,NULL);
    }
    else {
        generatefeat(Cce,rc,cam,time,pre_fid,fid,&pf,&cf);
        hashcopy(&pf,&cf);
    }
    pre_fid=fid;
    if (of) *of=&pf;
    return 1;
}
/* read motion profiles--------------------------------------------------------
 *  args:  char *file    I  motion profiles file path
 *         path_t *path  O  output all motion data
 *  note:
 *
 *     format of motion profiles:
 *         column 1: time (sec)
 *         column 2: latitude (rad)
 *         column 3: longitude (rad)
 *         column 4: height (m)
 *         column 5: north velocity (m/s)
 *         column 6: east velocity (m/s)
 *         column 7: down velocity (m/s)
 *         column 8: roll angle of body w.r.t NED (rad)
 *         column 9: pitch angle of body w.r.t NED (rad)
 *         column 10: yaw angle of body w.r.t NED (rad)
 *
 *  return : status (1: ok,0: fail)
 *  --------------------------------------------------------------------------*/
static int readmotionp(const char *file,path_t *path)
{
    double time,vn[3],rn[3],rpy[3],Cne[9],Cnb[9];
    char buff[1024];
    FILE *fp=NULL;
    motion_t mt;

    if ((fp=fopen(file,"r"))==NULL) return 0;
    while (fgets(buff,sizeof(buff),fp)) {
        if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&time,rn,rn+1,rn+2,
                   vn,vn+1,vn+2,rpy,rpy+1,rpy+2)<10) {
            continue;
        }
        mt.time.time=(time_t)time;
        mt.time.sec=time-(int)time;

        rn[0]*=D2R;
        rn[1]*=D2R;

        pos2ecef(rn,mt.re);
        ned2xyz(rn,Cne);
        matmul("NN",3,1,3,1.0,Cne,vn,0.0,mt.ve);

        rpy[0]*=D2R;
        rpy[1]*=D2R;
        rpy[2]*=D2R;

        rpy2dcm(rpy,Cnb);
        matmul("NT",3,3,3,1.0,Cne,Cnb,0.0,mt.Cbe);

        addmotion(path,&mt);
    }
    fclose(fp);
    return path->n;
}
/* generate a path contains image feature points, GPS and IMU data------------
 * args:    char *file      I  file path of truth trajectory
 *          cam_t *cam      I  camera model
 *          imu_err_t *err  I  IMU error model
 *          insopt_t *opt   I  ins options
 *          char *imufile   I  output IMU measurement data file path
 *          char *gpsfile   I  output GPS measurement data file path
 *          char *vofile    I  output feature points measurement data file path
 * return: status (1: ok,0: fail)
 * --------------------------------------------------------------------------*/
extern int generatepath(const char *file,const cam_t *cam,const imu_err_t *err,
                        const insopt_t *opt,
                        const char *imufile,
                        const char *gpsfile,const char *vofile)
{
    double *pCbe,*pre,*pve,*cCbe,*cre,*cve;
    double fb[3],omgb[3],dt,gps[3],var[3];
    path_t path={0};
    frame_t *pframe;
    gtime_t pt;
    FILE *fp_imu,*fp_vo,*fp_gps;
    int i,j,gps_hz=1,fhz=10,imu_hz=100;
    int ngps=0,nf=0,nimu=0;
    int gps_c,f_c;

    gps_c=imu_hz/gps_hz;
    f_c  =imu_hz/fhz;

    trace(3,"generatepath:\n");

    if (!(fp_imu=fopen(imufile,"w"))) return 0;
    if (!(fp_gps=fopen(gpsfile,"w"))) return 0;
    if (!(fp_vo =fopen(vofile ,"w"))) return 0;

    if (!readmotionp(file,&path)) return 0;
    for (pCbe=path.data[0].Cbe,pre=path.data[0].re,pve=path.data[0].ve,pt=path.data[0].time,
                 i=0;i<path.n;i++) {

        cCbe=path.data[i].Cbe;
        cre =path.data[i].re;
        cve =path.data[i].ve;
        dt=timediff(path.data[i].time,pt);

        /* generate imu data */
        if (!kinematicsecef(dt,cCbe,pCbe,cve,pve,cre,fb,omgb)) continue;

        /* simulate IMU errors */
        simimumeas(fb,omgb,err,dt);

        /* write imu data to file */
        fprintf(fp_imu,"%12.5lf  %15.8lf  %15.8lf  %15.8lf  %15.8lf  %15.8lf  %15.8lf\n",
                time2secs(path.data[i].time),fb[0],fb[1],fb[2],
                omgb[0],omgb[1],omgb[2]);
        fflush(fp_imu);

        if (nimu%gps_c==0) {
            
            /* generate GPS data */
            generategps(cCbe,cre,opt->lever,gps);

            /* add noise to gps position */
            for (j=0;j<3;j++) {
                var[j]=SQRT(VAR_POS)+getgaussian(0.005);
                gps[j]+=getgaussian(SQRT(VAR_POS));
            }
            /* write GPS data to file */
            fprintf(fp_gps,"%12.5lf  %15.8lf  %15.8lf  %15.8lf  %15.8lf  %15.8lf  %15.8lf\n",
                    time2secs(path.data[i].time),
                    gps[0],gps[1],gps[2],
                    var[0],var[1],var[2]);
            fflush(fp_gps);
            ngps++;
        }
        if (nimu%f_c==0) {
            /* generate image feature */
            if (!generateframe(cCbe,cre,nf,cam,opt,path.data[i].time,&pframe)) continue;

            /* write feature point to file */
            fprintf(fp_vo,"###,uid=%d,time=%8.4lf,nf=%d\n",pframe->uid+1,time2secs(pframe->time),HASH_COUNT(pframe->feat));

            feat_t *current,*tmp;
            HASH_ITER(hh,pframe->feat,current,tmp) {

                fprintf(fp_vo,"%6d  %8.4lf  %8.4lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  "
                              "%12.6lf  %12.6lf  %12.6lf\n",
                        current->uid,
                        current->u,
                        current->v,
                        current->re[0],
                        current->re[1],
                        current->re[2],

                        current->rc[0],
                        current->rc[1],
                        current->rc[2],

                        current->pc[0],
                        current->pc[1],
                        current->pc[2]);
            }
            fflush(fp_vo);
            nf++;
        }
        pt=path.data[i].time; pCbe=cCbe; pre=cre; pve=cve;
        nimu++;
    }
    fclose(fp_gps);
    fclose(fp_imu); 
    fclose(fp_vo );

    hash_rmframe(&pf.feat);
    hash_rmframe(&cf.feat);
    return nimu;
}