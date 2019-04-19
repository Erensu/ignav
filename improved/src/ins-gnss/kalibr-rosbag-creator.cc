/*------------------------------------------------------------------------------
* kalibr-rosbag-creator.cc : create imu-cam rosbag for Kalibr calibrator
*
* reference :
*    [1] https://github.com/ethz-asl/kalibr/wiki/bag-format
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/12/13 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

#define MAX_INDEX          1000

/* copy file------------------------------------------------------------------*/
static int copyfile(const char* source_path,const char *destination_path)
{
    char buffer[1024];
    FILE *in,*out;

    if ((in=fopen(source_path,"r"))==NULL) return 0;
    if ((out=fopen(destination_path,"w"))==NULL) return 0;

    int len;
    while ((len=fread(buffer,1,1024,in))>0) {
        fwrite(buffer,1,len,out);
    }
    fclose(out);
    fclose(in); return 1;
}
/* create imu-cam rosbag for Kalibr calibrator--------------------------------
 * args:    char *datfile  I  dat-format file for synchro image data
 *          char *imufilr  I  imu raw measurement data
 *          char *imgdir   I  image raw data dir
 *          char *output   I  output dir
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int kalibrrosbag(const char *datfile,const char *imufile,const char *imgdir,
                        const char *output)
{
    char imgfile_src[1024],imgfile_des[1024],imupath[1024];
    char camdir[1024];
    int flag,i,n=0; long int timestamp;

    raw_t raw={{0}};
    FILE *fp_dat;

    trace(3,"kalibrrosbag:\n");

    if (!(fp_dat=fopen(datfile,"rb"))) return 0;
    init_raw(&raw,1);

    sprintf(camdir,"%s/cam0",output);
    if (access(camdir,F_OK)==-1) {
        if (mkdir(camdir,0777)!=0) {free_raw(&raw); return 0;}
    }
    /* output image raw data */
    while (true) {
        flag=input_m39_mixf(&raw,fp_dat);
        if (flag==-2) break;
        if (flag==11) {
            if (n++>MAX_INDEX) break;

            if (get_m39_img(imgdir,raw.m39.fts.tv_sec,raw.m39.fts.tv_nsec,imgfile_src)) {
                timestamp=(long int)(time2gpst(raw.m39.time,NULL)*1E9);
                sprintf(imgfile_des,"%s/%ld.jpg",camdir,timestamp);

                copyfile(imgfile_src,imgfile_des);
            }
            else {
                fprintf(stderr,"lost image: %.4lf\n",raw.m39.sow);
            }
        }
    }
    fclose(fp_dat);

    /* output imu raw data */
    imu_t imu={0};
    readimub(imufile,&imu,IMUDECFMT_RATE,IMUFMT_GI310,
             IMUCOOR_RFU,IMUVALFMT_DEG);

    sortimudata(&imu);

    sprintf(imupath,"%s/imu.csv",output);
    if (!(fp_dat=fopen(imupath,"w"))) {
        free_raw(&raw); return 0;
    }
    fprintf(fp_dat,"timestamp,omega_x,omega_y,omega_z,alpha_x,alpha_y,alpha_z\n");
    for (i=0;i<imu.n;i++) {

        /* write imu measurement data */
        timestamp=(long int)(time2gpst(imu.data[i].time,NULL)*1E9);

        fprintf(fp_dat,"%ld,%lf,%lf,%lf,%lf,%lf,%lf\n",timestamp,
                imu.data[i].gyro[0]*D2R,
                imu.data[i].gyro[1]*D2R,
                imu.data[i].gyro[2]*D2R,
                imu.data[i].accl[0],
                imu.data[i].accl[1],
                imu.data[i].accl[2]);
    }
    fclose(fp_dat);
    free_raw(&raw); freeimudata(&imu);
    return 1;
}