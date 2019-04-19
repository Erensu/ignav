/*-----------------------------------------------------------------------------
* lc-fbsm.cc : ins-gnss loosely coupled fwd/bwd fixed-point smoother app.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/10/01 1.0 new
*----------------------------------------------------------------------------*/
#include <carvig.h>

/* receiver options table ----------------------------------------------------*/
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,7:ntripcli,8:ftp,9:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripsvr,11:ntripc_c"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,3:oem3,4:ubx,5:ss2,6:hemis,7:skytraq,8:gw10,9:javad,10:nvs,11:binex,12:rt17,13:sbf,14:cmr,15:tersus,18:sp3,19:rnxclk,20:sbas,21:nmea,22:gsof,23:ublox-evk-m8u,24:ublox-sol,25:m39,26:rinex,27:m39-mix,28:euroc-imu,29:euroc-img,30:karl-img,31:malaga-gnss,32:malaga-imu,33:malaga-img,34:oem6-sol,35:oem6-pose,36:oem6-raw"
#define NMEOPT  "0:off,1:latlon,2:single"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:stat,5:gsif,6:ins"
#define MSGOPT  "0:all,1:rover,2:base,3:corr"

#define OPTSDIR     "."                 /* default config directory */
#define OPTSFILE    "rtkrcv.conf"       /* default config file */
#define MAXSTR      1024                /* max length of a stream */

/* help text -----------------------------------------------------------------*/
static const char *usage[]={
        "usage: lc-rts [-p port][-o file]",
        "options",
        "  -p port    port number for telnet console",
        "  -o file    processing options file",
};
static int strtype[]={                  /* stream types */
        STR_SERIAL,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,
        STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE
};
static char strpath[14][MAXSTR]={"","","","","","","","","","","","","",""}; /* stream paths */
static int strfmt[]={                   /* stream formats */
        STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,STRFMT_UBXM8,STRFMT_UBXSOL,SOLF_LLH,SOLF_NMEA,
        SOLF_LLH,SOLF_LLH
};
static opt_t rcvopts[]={
        {"inpstr1-type",    3,  (void *)&strtype[0],         ISTOPT },
        {"inpstr2-type",    3,  (void *)&strtype[1],         ISTOPT },
        {"inpstr3-type",    3,  (void *)&strtype[2],         ISTOPT },
        {"inpstr4-type",    3,  (void *)&strtype[3],         ISTOPT },
        {"inpstr5-type",    3,  (void *)&strtype[4],         ISTOPT },
        {"inpstr6-type",    3,  (void *)&strtype[5],         ISTOPT },
        {"inpstr7-type",    3,  (void *)&strtype[6],         ISTOPT },
        {"inpstr1-path",    2,  (void *)strpath [0],         ""     },
        {"inpstr2-path",    2,  (void *)strpath [1],         ""     },
        {"inpstr3-path",    2,  (void *)strpath [2],         ""     },
        {"inpstr4-path",    2,  (void *)strpath [3],         ""     },
        {"inpstr5-path",    2,  (void *)strpath [4],         ""     },
        {"inpstr6-path",    2,  (void *)strpath [5],         ""     },
        {"inpstr7-path",    2,  (void *)strpath [6],         ""     },
        {"inpstr1-format",  3,  (void *)&strfmt [0],         FMTOPT },
        {"inpstr2-format",  3,  (void *)&strfmt [1],         FMTOPT },
        {"inpstr3-format",  3,  (void *)&strfmt [2],         FMTOPT },
        {"inpstr4-format",  3,  (void *)&strfmt [3],         FMTOPT },
        {"inpstr5-format",  3,  (void *)&strfmt [4],         FMTOPT },
        {"inpstr6-format",  3,  (void *)&strfmt [5],         FMTOPT },
        {"inpstr7-format",  3,  (void *)&strfmt [6],         FMTOPT },
        {"outstr1-type",    3,  (void *)&strtype[7],         OSTOPT },
        {"outstr2-type",    3,  (void *)&strtype[8],         OSTOPT },
        {"outstr1-path",    2,  (void *)strpath [7],         ""     },
        {"outstr2-path",    2,  (void *)strpath [8],         ""     },
        {"outstr1-format",  3,  (void *)&strfmt [7],         SOLOPT },
        {"outstr2-format",  3,  (void *)&strfmt [8],         SOLOPT },
        {"logstr1-type",    3,  (void *)&strtype[9],         OSTOPT },
        {"logstr2-type",    3,  (void *)&strtype[10],        OSTOPT },
        {"logstr3-type",    3,  (void *)&strtype[11],        OSTOPT },
        {"logstr4-type",    3,  (void *)&strtype[12],        OSTOPT },
        {"logstr5-type",    3,  (void *)&strtype[13],        OSTOPT },
        {"logstr1-path",    2,  (void *)strpath [9],         ""     },
        {"logstr2-path",    2,  (void *)strpath [10],        ""     },
        {"logstr3-path",    2,  (void *)strpath [11],        ""     },
        {"logstr4-path",    2,  (void *)strpath [12],        ""     },
        {"logstr5-path",    2,  (void *)strpath [13],        ""     },

        {"",0,NULL,""}
};
static prcopt_t prcopt;                 /* processing options */
static solopt_t solopt={0};             /* solution options */
static filopt_t filopt={""};            /* file options */
static imu_t       imu={0};             /* imu measurement data */
static gsof_data_t pos={0};             /* position measurement data */

/* print usage ---------------------------------------------------------------*/
static void printusage(void)
{
    int i;
    for (i=0;i<(int)(sizeof(usage)/sizeof(*usage));i++) {
        fprintf(stderr,"%s\n",usage[i]);
    }
    exit(0);
}
/* adjust imu measurement data-----------------------------------------------*/
static void adj_imudata(const prcopt_t *opt,imu_t *imu)
{
    int i;
    for (i=0;i<imu->n;i++) adjustimu(opt,&imu->data[i]);
}
/* read imu measurement data-------------------------------------------------*/
static int readimu(const char *file,int type)
{
    int nimu=0;
    switch (type) {
        case STRFMT_M39: readimub(file,&imu,prcopt.insopt.imudecfmt,
                                  prcopt.insopt.imuformat,
                                  prcopt.insopt.imucoors,
                                  prcopt.insopt.imuvalfmt);
    }
    /* sort imu measurement data */
    nimu=sortimudata(&imu);

    /* adjust imu measurement data */
    adj_imudata(&prcopt,&imu);
    return nimu;
}
/* read position measurement data--------------------------------------------*/
static int readpos(const char *file,int type)
{
    int npos=0;
    switch (type) {
        case STRFMT_GSOF: readgsoff(file,&pos);
    }
    /* sort gsof measurement data */
    npos=sortgsof(&pos);
    return npos;
}
/* FB-SM main------------------------------------------------------------------
* sysnopsis
*     lc-fbsm [-p port] [-o file]
*
* description
*     A command line version of RTS smoother for ins and gnss loosely coupled.
*
* option
*     -p  port number for monitor stream
*     -o  processing options file
*
* --------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int i,port=0;
    char file[1024];

    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-p")&&i+1<argc) port=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-o")&&i+1<argc) strcpy(file,argv[++i]);
        else printusage();
    }
    /* load options file */
    if (!*file) sprintf(file,"%s/%s",OPTSDIR,OPTSFILE);
    resetsysopts();
    if (!loadopts(file,rcvopts)||!loadopts(file,sysopts)||
        !loadopts(file,insopts)) {
        fprintf(stderr,"no options file: %s. defaults used\n",file);
    }
    getsysopts(&prcopt,&solopt,&filopt);

    /* read imu/pos data */
    readimu(strpath[4],strfmt[4]);
    readpos(strpath[3],strfmt[3]);

    /* set forward solution-binary file path */
    set_fwdtmp_file(NULL);

    /* forward/backward combined */
    lcfbsm(&imu,&pos,&prcopt,&solopt,port,strpath[7]);

    freeimudata(&imu);
    freegsofdata(&pos);
    return 0;
}