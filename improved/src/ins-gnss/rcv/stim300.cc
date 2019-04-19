/*------------------------------------------------------------------------------
* stim300.cc : decode STIM300 imu raw measurement data
*
* version : $Revision:$ $Date:$
* history : 2019/04/14  1.0  new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* read STIM300 IMU raw data--------------------------------------------------
 * args:  char* file  I  STIM300 IMU data file
 *        imu_t *imu  O  imu measurement data
 * return: epochs of imu measurement data
 * --------------------------------------------------------------------------*/
extern int readstim300(const char *file,imu_t *imu)
{
    imud_t imud={0};
    double val[7];
    FILE *handle;

    trace(3,"readstim300:\n");

    if (!(handle=fopen(file,"rb"))) {
        trace(2,"open imu file fail\n");
        return 0;
    }
    while (fread(val,sizeof(double),7,handle)>0) {
        imud.time.time=(time_t)val[0];
        imud.time.sec =val[0]-(time_t)val[0];

        imud.gyro[0]=val[1];
        imud.gyro[1]=val[2];
        imud.gyro[2]=val[3];
        imud.accl[0]=val[4];
        imud.accl[1]=val[5];
        imud.accl[2]=val[6];
        addimudata(imu,&imud);
    }
    fclose(handle);
    return imu->n;
}