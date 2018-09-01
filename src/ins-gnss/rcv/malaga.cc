/*----------------------------------------------------------------------------
* kitti.cc : read the Malaga Urban Dataset
*
* version : $Revision:$ $Date:$
* history : 2018/04/20  1.0  new
*----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants and macros ------------------------------------------------------*/
#define MAXFIELD   64               /* max number of fields in a record */

/* convert gsof message to solution type--------------------------------------*/
static void gsof2sol(const gsof_t *gsof,sol_t *sol)
{
    int i; double Cne[9];

    trace(3,"gsof2sol:\n");

    sol->time=gsof->t;

    pos2ecef(gsof->llh,sol->rr);
    ned2xyz (gsof->llh,Cne);
    matmul("NN",3,1,3,1.0,Cne,gsof->vel,0.0,sol->rr+3);

    for (i=0;i<3;i++) {
        sol->qr[i]=SQR(gsof->sig[i]);
    }
    sol->ns  =gsof->ns;
    sol->stat=gsof->solq;
}
/* sync a newline------------------------------------------------------------*/
static int syncnewline(unsigned char *buff, int nb)
{
    if (buff[nb-1]=='\n'||(buff[nb-2]=='\r'&&buff[nb-1]=='\n')) {
        if (buff[nb-1]=='\n') buff[nb-1]='\0';
        if (buff[nb-2]=='\r'&&buff[nb-1]=='\n') {
            buff[nb-2]='\0';
            buff[nb-1]='\0';
        }
        return 1;
    }
    return 0;
}
/* read gnss position measurement data---------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_malaga_gnss(raw_t *raw, unsigned char data)
{
    static double varpos[4]={0.3,0.5,1.0,3.0};
    int i; double tmp[25]={0};

    trace(3,"input_malaga_gnss: data=%02x\n",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (strstr((char*)raw->buff,"%")) return 0;

    if (sscanf((char*)raw->buff,
               "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
               &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6],&tmp[7],&tmp[8],&tmp[9],
               &tmp[10],&tmp[11],&tmp[12],&tmp[13],&tmp[14],
               &tmp[15],&tmp[16],&tmp[17],&tmp[18],&tmp[19],
               &tmp[20],&tmp[21],&tmp[22],&tmp[23],&tmp[24])<25) {
        return 0;
    }
    raw->gsof.status=(unsigned char)tmp[4];
    raw->gsof.solq  =(unsigned char)tmp[4];
    raw->gsof.ns    =(unsigned char)tmp[5];

    raw->gsof.t.time=(time_t)tmp[0];

    for (i=0;i<3;i++) {
        raw->gsof.llh[i]=tmp[1+i];
    }
    pos2ecef(raw->gsof.llh,raw->gsof.pos);
    for (i=0;i<3;i++) {
        raw->gsof.sig[i]=(float)varpos[raw->gsof.solq-1];
    }
    gsof2sol(&raw->gsof,&raw->sol);
    return 6;
}
/* read imu measurement data--------------------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_malaga_imu(raw_t *raw, unsigned char data)
{
    double tmp[25]={0};

    trace(3,"input_kitti_imu: data=%02x\n",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (strstr((char*)raw->buff,"%")) return 0;

    if (sscanf((char*)raw->buff,
               "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
               &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6],&tmp[7],
               &tmp[8],&tmp[9],&tmp[10],&tmp[11],&tmp[12],&tmp[13],&tmp[14],&tmp[15])<16) {
        return 0;
    }
    raw->imut.n=0;
    raw->imu.time.time=(time_t)tmp[0];
    raw->imu.time.sec =tmp[0]-(time_t)tmp[0];

    raw->imu.accl[0]= tmp[1];
    raw->imu.accl[1]=-tmp[2];
    raw->imu.accl[2]=-tmp[3]; /* accl-measurement (m/s^2) */

    raw->imu.gyro[0]= tmp[6];
    raw->imu.gyro[1]=-tmp[5];
    raw->imu.gyro[2]=-tmp[4]; /* angular velocities (rad/s) */

    addimudata(&raw->imut,&raw->imu);
    return 4;
}
/* enu-attitude convert to ned-attitude--------------------------------------*/
static void convatt(const double *a,double *b)
{
    static double Cbb[9]={1,0,0,0,-1,0,0,0,-1},Cnn[9];
    static double Cnb1[9],Cnb2[9];

    enu2ned(NULL,NULL,Cnn);
    rpy2dcm(a,Cnb1);

    matmul33("NNT",Cbb,Cnb1,Cnn,3,3,3,3,Cnb2);
    dcm2rpy(Cnb2,b);
}
/* read vehicle attitude data------------------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_malga_att(raw_t *raw, unsigned char data)
{
    int i; double tmp[25]={0},att[3];

    trace(3,"input_malga_att: data=%02x\n",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (strstr((char*)raw->buff,"%")) return 0;

    if (sscanf((char*)raw->buff,
               "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
               &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6],&tmp[7],
               &tmp[8],&tmp[9],&tmp[10],&tmp[11],&tmp[12],&tmp[13],&tmp[14],&tmp[15])<16) {
        return 0;
    }
    raw->att.time.time=(time_t)tmp[0];
    raw->att.time.sec =tmp[0]-(time_t)tmp[0];
    
    att[0]=tmp[12];
    att[1]=tmp[11];
    att[2]=tmp[10];
    convatt(att,raw->att.val); /* enu convert to ned */
    for (i=0;i<3;i++) {
        raw->att.std[i]=SQR(1.0*D2R);
    }
    return 33;
}
/* read left and right image data--------------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_malaga_img(raw_t *raw, unsigned char data)
{
    char str[8],path[256];
    unsigned char *tmp=NULL;
    double time=0.0;
    static int flag=0;
    size_t x,y;

    trace(3,"input_malaga_img:\n");

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;

    if (sscanf((char*)raw->buff,"img_CAMERA1_%lf_%s\n",&time,str)<2) {
        return 0;
    }
    raw->img.time.time=(time_t)time;
    raw->img.time.sec =time-(time_t)time;

    if (strstr(str,"left")) { /* left image */

        sprintf(path,"%s/%s",raw->monodir,(char*)raw->buff);

        /* read image raw data */
        if (!readjpeg(path,raw->img.time,&raw->img,0)) {
            trace(2,"no left image raw data\n");
            return 0;
        }
        raw->img.id--;
        flag++;
    }
    else if (strstr(str,"right")) { /* right image */

#if 0
        sprintf(path,"%s/%s",raw->monodir,(char*)raw->buff);

        /* read image raw data */
        if (!readjpeg(path,raw->img.time,&raw->img,1)) {
            trace(2,"no left image raw data\n");
            return 0;
        }
#endif
        raw->img.id--;
        flag++;
    }
    if (flag==2&&!(flag=0)) {
        raw->img.id++; return 11;
    }
    else return 0;
}
