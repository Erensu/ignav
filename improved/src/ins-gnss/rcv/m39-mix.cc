/*----------------------------------------------------------------------------
* m39-mix.cc : decode m39/image time tag raw data
*
* version : $Revision:$ $Date:$
* history : 2017/03/08  1.0  new
*----------------------------------------------------------------------------*/
#include <carvig.h>

/* constants/macros ---------------------------------------------------------*/
#define M39SYNC1   0x55        /* m39-mix message start sync code 1 */
#define M39SYNC2   0xAA        /* m39-mix message start sync code 2 */
#define M39SYNC3   0x99        /* m39-mix message start sync code 3 */
#define M39SYNC4   0x12        /* m39-mix message start sync code 4 */
#define M39MSGLEN  27          /* m39-mix message length  */
#define CENTURY    21          /* century of current year */
#define MINCHGC    10          /* min counts of zda-time changing */
#define READIMG    1           /* read image raw data when post-process */

/* check sum-----------------------------------------------------------------*/
static int chksum(const unsigned char *buff, int len, unsigned char data)
{
    int i; unsigned char sum=0;
    for (i=0;i<len;i++) sum+=buff[i]; return sum==data;
}
/* timespec struct convert to seconds----------------------------------------*/
static double tsp2secs(struct timespec tv)
{
    return tv.tv_sec+tv.tv_nsec*1E-9;
}
/* evaluate sow by using pps-------------------------------------------------*/
static int pps2sow(const double pps,const double zda)
{
    static double ppsp=0.0,zdap=0.0;
    static int flag=0;

    if (ppsp==0.0&&zdap==0.0) {ppsp=pps; zdap=zda; return 0;}

    /* pps always early compare with zda */
    if (fabs(ppsp-pps)>0.1) flag=1;
    if (flag==1) {
        if (fabs(zdap-zda)>0.1) flag=2;
    }
    ppsp=pps;
    zdap=zda;
    if (flag==2) {flag=0; return 1;}
    else return 0;
}
/* get image for given frame time--------------------------------------------
 * args:    char *dir          I  directories of image data
 *          uint64_t fts,ftns  I  time of image data
 *          char *imgfile      O  file path of image data
 * return: status (1:ok, 0:fail)
 * --------------------------------------------------------------------------*/
extern int get_m39_img(const char *dir,const uint64_t fts,uint64_t ftns,
                       char *imgfile)
{
    DIR *dp;
    struct dirent *entry;
    struct stat statbuf;
    char file[128];

    if ((dp=opendir(dir))==NULL) {
        fprintf(stderr,"cannot open directory: %s\n",dir);
        return 0;
    }
    chdir(dir);
    sprintf(file,"%ld_%ld",fts,ftns);

    while ((entry=readdir(dp))!=NULL) {
        lstat(entry->d_name,&statbuf);
        if (strstr(entry->d_name,file)) {
            sprintf(imgfile,"%s/%s",dir,entry->d_name);
            closedir(dp); return 1;
        }
    }
    closedir(dp);
    return 0;
}
/* read image data-----------------------------------------------------------*/
static int readimgbuf(raw_t *raw)
{
    prcopt_t *popt=(prcopt_t*)raw->optp;
    char imgpath[1024],*pstr;

    sprintf(imgpath,"%s/%d_%ld_%ld.jpg",popt->monodir,raw->img.id+1,
            raw->m39.fts.tv_sec,
            raw->m39.fts.tv_nsec);

    /* load jpg image file */
    if (!readjpeg(imgpath,raw->m39.time,&raw->img,0)) {
        if (get_m39_img(popt->monodir,raw->m39.fts.tv_sec,raw->m39.fts.tv_nsec,imgpath)) {

            if (!readjpeg(imgpath,raw->m39.time,&raw->img,0)) {
                trace(2,"read jpg image fail: time=%sl\n",time_str(raw->m39.time,4));
                return 0;
            }
            if ((pstr=strrchr(imgpath,'/'))) {
                sscanf(pstr,"/%d_%ld_%ld.jpg",&raw->img.id,
                       &raw->m39.fts.tv_sec,
                       &raw->m39.fts.tv_nsec);
            }
            goto ok;
        }
        else {
            trace(2,"no such file: time=%s\n",
                  time_str(raw->m39.time,4)); return 0;
        }
    }
ok:
#if 0
    /* display image for debug */
    dipsplyimg(&raw->img);
#endif
    return 11;
}
/* decode m39-mix raw data---------------------------------------------------*/
static int decode_m39_mix(raw_t *raw)
{
    static int flag=0,j=0,week;
    static double ppsp=0.0;
    unsigned char *p=raw->buff+4;

    uint64_t ts=0,tn=0;
    gtime_t time;
    double epoch[6],dt;

    /* pps time in tx2-clock */
    ts=p[0]; ts=ts<<8|p[1]; ts=ts<<8|p[2]; ts=ts<<8|p[3];
    tn=p[4]; tn=tn<<8|p[5]; tn=tn<<8|p[6]; tn=tn<<8|p[7];
    raw->m39.pps.tv_sec=ts; raw->m39.pps.tv_nsec=tn;

    /* gps nmea time */
    epoch[0]=p[8 ]; epoch[1]=p[9 ]; epoch[2]=p[10];
    epoch[3]=p[11]; epoch[4]=p[12]; epoch[5]=p[13];
    epoch[0]+=(CENTURY-1)*100.0;
    raw->m39.zda=time=epoch2time(epoch);

    time2gpst(time,&week); /* GPS week */

    /* frame time */
    ts=p[14]; ts=ts<<8|p[15]; ts=ts<<8|p[16]; ts=ts<<8|p[17];
    tn=p[18]; tn=tn<<8|p[19]; tn=tn<<8|p[20]; tn=tn<<8|p[21];
    raw->m39.fts.tv_sec=ts; raw->m39.fts.tv_nsec=tn;

    /* adjust frame time */
    if (!flag) {
        j=pps2sow(tsp2secs(raw->m39.pps),time2gpst(time,NULL));
        raw->m39.sow=time2gpst(time,NULL); 
        if (j&&raw->m39.sowc++>MINCHGC) flag=1;
    }
    if (flag) {
        if (ppsp!=0.0) {
            raw->m39.sow+=ROUND(tsp2secs(raw->m39.pps)-ppsp);
        }
        ppsp=tsp2secs(raw->m39.pps);

        /* get frame time */
        dt=tsp2secs(raw->m39.fts)-tsp2secs(raw->m39.pps);
        raw->m39.time=gpst2time(week,raw->m39.sow+dt);
#if READIMG
        /* read image data */
        return readimgbuf(raw);
#else
        return 11;
#endif
    }
    else return 0; /* fail */
}
/* sync frame----------------------------------------------------------------*/
static int sync_m39_mix(unsigned char *buff, unsigned char data)
{
    buff[0]=buff[1]; buff[1]=buff[2];
    buff[2]=buff[3]; buff[3]=data;
    return buff[0]==M39SYNC1&&buff[1]==M39SYNC2&&
           buff[2]==M39SYNC3&&buff[3]==M39SYNC4;
}
/* read m39 mix raw data (included image/imu)--------------------------------
 * args  :  raw_t *raw         IO  raw struct
 *          unsigned char data I   raw data
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_m39_mix(raw_t *raw, unsigned char data)
{
    trace(5,"input_m39_mix: data=%02x\n",data);

    /* synchronize frame */
    if (raw->nbyte==0) {
        if (sync_m39_mix(raw->buff,data)) raw->nbyte=4;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;
    if (raw->nbyte>MAXRAWLEN) {raw->nbyte=0; return 0;}
    if (raw->nbyte<M39MSGLEN) return 0;

    raw->nbyte=0;

    /* check sum */
    if (!chksum(raw->buff+2,M39MSGLEN-3,data)) {
        return 0;
    }
    /* decode oem4 message */
    return decode_m39_mix(raw);
}
/* read m39-mix raw data from file-------------------------------------------
 * args   : raw_t  *raw   IO     receiver raw data control struct
 *          FILE   *fp    I      file pointer
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_m39_mixf(raw_t *raw,FILE *fp)
{
    int i,data,ret;

    trace(3,"input_m39_mixf:\n");

    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_m39_mix(raw,(unsigned char)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}