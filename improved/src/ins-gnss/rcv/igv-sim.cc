/*--------------------------------------------------------------------------------
 * igv-sim.cc : read ins-gnss-vo multisensor simulator data
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2019/01/08 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants/macros ---------------------------------------------------------*/
static const double gpst0[]={1980,1, 6,0,0,0};   /* gps time reference */
static FILE *fp_vo=NULL;                         /* file pointer of feature point measurement data file */

/* sync a newline------------------------------------------------------------*/
static int syncnewline(unsigned char *buff, int nb)
{
    if (buff[nb-1]=='\n'||(buff[nb-2]=='\r'&&buff[nb-1]=='\n')) return 1;
    return 0;
}
/* input ins-gnss-vo multisensor imu simulator data --------------------------
 * args   : raw_t *raw         IO     receiver raw data control struct
 *          unsigned char data I stream data (1 byte)
 * return >1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_igvsim_imu(raw_t *raw,unsigned char data)
{
    double tmp[7]={0};
    gtime_t t0;

    trace(3,"input_igvsim_imu: data=%02x",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (sscanf((char*)raw->buff,
               "%lf %lf %lf %lf %lf %lf %lf\n",
               &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6])<7) {
        return 0;
    }
    raw->imut.n=0;

    t0=epoch2time(gpst0);
    raw->imu.time=timeadd(t0,tmp[0]);

    raw->imu.accl[0]=tmp[1];
    raw->imu.accl[1]=tmp[2];
    raw->imu.accl[2]=tmp[3]; /* accl-measurement (m/s^2) */

    raw->imu.gyro[0]=tmp[4];
    raw->imu.gyro[1]=tmp[5];
    raw->imu.gyro[2]=tmp[6]; /* angular velocities (rad/s) */

    addimudata(&raw->imut,&raw->imu);
    return 4;
}
/* read gnss position measurement data---------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_igvsim_gnss(raw_t *raw, unsigned char data)
{
    int i; double tmp[25]={0};
    gtime_t t0;

    trace(3,"input_igvsim_gnss: data=%02x",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (sscanf((char*)raw->buff,
               "%lf %lf %lf %lf %lf %lf %lf\n",
               &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6])<7) {
        return 0;
    }
    raw->sol.stat=SOLQ_FIX;
    raw->sol.ns  =9;

    t0=epoch2time(gpst0);
    raw->sol.time=timeadd(t0,tmp[0]);

    for (i=0;i<3;i++) {
        raw->sol.rr[i]=tmp[1+i];
        raw->sol.qr[i]=tmp[4+i];
    }
    return 6;
}
/* add a feature to a hash table---------------------------------------------*/
static void hash_addfeat(feature **feat,double u,double v,gtime_t time,int uid)
{
    feature *s=NULL;
    HASH_FIND_INT(*feat,&uid,s);
    if (s==NULL) {

        s=(feature*)malloc(sizeof(feature));
        s->uid=uid;
        s->u=u;
        s->v=v;
        s->time=time;
        HASH_ADD_INT(*feat,uid,s);
    }
    else {
        s->time=time;
        s->uid=uid;
        s->u=u;
        s->v=v;
    }
}
/* add a frame to hash table--------------------------------------------------*/
static void hash_addframe(img_t *img,img_t **hash)
{
    feature *current,*tmp;

    img_t *imgp=(img_t*)malloc(sizeof(img_t));
    imgp->w=img->w;
    imgp->h=img->h;

    imgp->id  =img->id;
    imgp->time=img->time;
    imgp->feat=NULL;
    imgp->data=NULL;
    HASH_ITER(hh,img->feat,current,tmp) {
        hash_addfeat(&imgp->feat,current->u,current->v,current->time,current->uid);
    }
    HASH_ADD_INT(*hash,id,imgp);
}
/* delete a frame-------------------------------------------------------------*/
static void hash_rmframe(img_t **hash)
{
    feature *fcurrent,*ftmp;
    img_t *current,*tmp;

    HASH_ITER(hh,*hash,current,tmp) {
        HASH_DEL(*hash,current);
        HASH_ITER(hh,current->feat,fcurrent,ftmp) {

            /* remove feature in frame */
            HASH_DEL(current->feat,fcurrent);
            free(fcurrent);
        }
        free(current); /* delete frame */
    }
    *hash=NULL;
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
/* read feature point measurement data from file------------------------------*/
static int readfeatfile(FILE *fp,img_t **hashimg,int nfmax)
{
    feature *fcurrent,*ftmp;
    raw_t raw={0};
    int flag,count=0;

    init_raw(&raw,STRFMT_IGVSIM_FEAT);
    while (true) {
        flag=input_igvsim_featf(&raw,fp);
        if (flag==-2) break;

        if (flag==11) {
            hash_addframe(&raw.img,hashimg);
            if (++count==nfmax) break;
        }
    }
    hash_rmimgfeat(&raw.img.feat);
    free_raw(&raw);
    return HASH_COUNT(*hashimg);
}
/* read feature points measurement data---------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_igvsim_feat(raw_t *raw, unsigned char data)
{
    static double time,u,v,xyz[3],pc[3];
    static int fid,nf,count,id;
    gtime_t t0;

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (strstr((char*)raw->buff,"###")) {
        sscanf((char*)raw->buff,"###,uid=%d,time=%lf,nf=%d\n",&fid,&time,&nf);
        t0=epoch2time(gpst0);
        raw->img.id=fid;
        raw->img.time=timeadd(t0,time);

        count=0;
        return 0;
    }
    /* decode feature point */
    if (sscanf((char*)raw->buff,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",
               &id,&u,&v,xyz,xyz+1,xyz+2,pc,pc+1,pc+2)<9) {
        return 0;
    }
    /* add to hash table */
    hash_addfeat(&raw->img.feat,u,v,raw->img.time,id);
    if (++count==nf) {count=0; return 11;}
    return 0;
}
extern int input_igvsim_featall(raw_t *raw, unsigned char data)
{
    static int input_count=0;
    feature *fcurrent,*ftmp;
    img_t *hash=NULL;

    trace(3,"input_igvsim_featall:\n");

    if (input_count++<10) return 0;
    if (fp_vo==NULL) {
        if (!(fp_vo=fopen(((stream_t*)raw->strp)->path,"r"))) return 0;
    }
    if (!readfeatfile(fp_vo,&hash,1)) return 0;
    HASH_ITER(hh,hash->feat,fcurrent,ftmp) {

        /* add a feature point */
        hash_addfeat(&raw->img.feat,fcurrent->u,fcurrent->v,
                     fcurrent->time,
                     fcurrent->uid);
    }
    input_count=0;

    raw->img.time=hash->time;
    raw->img.id=hash->id;
    hash_rmframe(&hash);
    return 11;
}
/* read feature point measurement data from file------------------------------
 * args   : raw_t  *raw   IO     receiver raw data control struct
 *          FILE   *fp    I      file pointer
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_igvsim_featf(raw_t  *raw,FILE *fp)
{
    int i,data,ret;

    trace(3,"input_igvsim_featf:\n");

    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_igvsim_feat(raw,(unsigned char)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}
/* free igv-sim file pointer--------------------------------------------------*/
extern void freeigvfp()
{
    if (fp_vo) fclose(fp_vo);
}




