/*----------------------------------------------------------------------------
* euroc.cc : read the EuRoC MAV dataset
*
* version : $Revision:$ $Date:$
* history : 2017/03/20  1.0  new
*----------------------------------------------------------------------------*/
#include <navlib.h>

/* sync a newline------------------------------------------------------------*/
static int syncnewline(unsigned char *buff, int nb)
{
    if (buff[nb-1]=='\n'||(buff[nb-2]=='\r'&&buff[nb-1]=='\n')) return 1;
    return 0;
}
/* read euroc imu measurement data--------------------------------------------
 * args  :  raw_t *raw         IO  raw struct
 *          unsigned char data I   raw data
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_imu_euroc(raw_t *raw, unsigned char data)
{
    int i; double val[7]={0};

    trace(3,"input_imu_euroc: data=%02x\n",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (strstr((char*)raw->buff,"#")) return 0;

    /* decode imu measurement data */
    sscanf((char*)raw->buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
           val,val+1,val+2,val+3,val+4,val+5,val+6);

    val[0]*=1E-9; /* time */
    raw->imu.time.time=(time_t)val[0];
    raw->imu.time.sec =val[0]-(time_t)val[0];

    /* gyro/accl measurement data */
    for (i=0;i<3;i++) {
        raw->imu.gyro[i]=val[1+i]; /* rad/s */
        raw->imu.accl[i]=val[4+i]; /* m/s^2 */
    }
    /* add imu measurement data */
    raw->imut.n=0;
    addimudata(&raw->imut,&raw->imu);
    return 4;
}
/* read euroc image data-----------------------------------------------------
 * args  : raw_t *raw         IO  raw struct
 *         unsigned char data I   raw data
 * return >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_img_euroc(raw_t *raw, unsigned char data)
{
    double time=0.0; char path[126],filename[126];
    unsigned char *tmp=NULL;
    size_t x,y;
    trace(3,"input_img_euroc: data=%02x\n",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    raw->nbyte=0;
    if (strstr((char*)raw->buff,"#")) return 0;
    if (sscanf((char*)raw->buff,"%lf,%s\n",&time,filename)<2) {
        return 0;
    }
    time=time*1E-9;
    /* read image raw data */
    sprintf(path,"%s/%s",raw->monodir,filename);
#if 0
    tmp=io_png_read_u8_gray(path,&x,&y);
#endif
    /* save image raw data */
    raw->img.w=x; raw->img.h=y;
    raw->img.dims[0]=(int)x;
    raw->img.dims[1]=(int)y;
    raw->img.dims[2]=(int)x;
    raw->img.id++; raw->img.flag=0;
    raw->img.data=tmp;
    raw->img.time.time=(time_t)time;
    raw->img.time.sec =time-(time_t)time;
    return 11;
}