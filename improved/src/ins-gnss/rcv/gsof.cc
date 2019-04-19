/*------------------------------------------------------------------------------
* gsof.cc: Trimble General Serial Output Format
*
* reference :
*     [1] http://www.trimble.com/OEM_ReceiverHelp/V5.11/en/#GSOFmessages_GSOF.html
*
* version : $Revision:$ $Date:$
* history : 2017/10/27  1.0  new
*-----------------------------------------------------------------------------*/
#include <carvig.h>

#define GSOFHEADLEN   4                       /* number bytes of gsof message header */
#define GSOFMIDLEN    3                       /* number bytes of gsof message information */
#define GSOFMMLEN     2                       /* number bytes of gsof message record header */
#define GSOFTYPE01    1                       /* type of gsof message */
#define GSOFTYPE02    2
#define GSOFTYPE03    3
#define GSOFTYPE06    6
#define GSOFTYPE08    8
#define GSOFTYPE09    9
#define GSOFTYPE10    10
#define GSOFTYPE11    11
#define GSOFTYPE12    12
#define GSOFTYPE38    38
#define GSOFTYPE41    41

/* extract field (big-endian) ------------------------------------------------*/
#define U1(p)       (*((unsigned char *)(p)))
#define I1(p)       (*((signed char *)(p)))

/* used by the decoding routines to grab 4 bytes and pack them into a u32----*/
static unsigned long getu32(unsigned char **data)
{
    unsigned int ret;
    unsigned char *pb;

    pb=(unsigned char*)(&ret)+3;
    *pb--=*(*data)++; *pb--=*(*data)++;
    *pb--=*(*data)++; *pb  =*(*data)++;
    return ret;
}
/* used by the decoding routines to grab 4 bytes and pack them into a float --*/
static float getfloat(unsigned char **data)
{
    float ret;
    unsigned char *pb;

    pb=(unsigned char*)(&ret)+3;
    *pb--=*(*data)++; *pb--=*(*data)++;
    *pb--=*(*data)++; *pb  =*(*data)++;
    return ret;
}
/* used by the decoding routines to grab 8 bytes and pack them into a double--*/
static double getdouble(unsigned char**data)
{
    double ret;
    unsigned char *pb;

    pb=(unsigned char *)(&ret)+7;
    *pb--=*(*data)++; *pb--=*(*data)++;
    *pb--=*(*data)++; *pb--=*(*data)++;
    *pb--=*(*data)++; *pb--=*(*data)++;
    *pb--=*(*data)++; *pb  =*(*data)++;
    return ret;
}
/* used by the decoding routines to grab 2 bytes and pack them into a 16 -----*/
static unsigned short getu16(unsigned char**data)
{
    unsigned short ret;
    unsigned char *pb;

    pb=(unsigned char*)(&ret)+1;
    *pb--=*(*data)++; *pb=*(*data)++;
    return ret;
}
/* decode gsof message type 10 record ----------------------------------------*/
static void decode_type10(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type10 :\n");

    raw->gsof.clk.flags=data[0]; data++;
    raw->gsof.clk.off  =getdouble(&data)/1000.0;
    raw->gsof.clk.foff =getdouble(&data)/1E6;
}
/* decode gsof message type 02 record ----------------------------------------*/
static void decode_type02(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type02 :\n");

    raw->gsof.llh[0]=getdouble(&data);
    raw->gsof.llh[1]=getdouble(&data);
    raw->gsof.llh[2]=getdouble(&data);

    pos2ecef(raw->gsof.llh,raw->gsof.pos);
}
/* decode gsof message type 03 record-----------------------------------------*/
static void decode_type03(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type03 :\n");

    raw->gsof.pos[0]=getdouble(&data);
    raw->gsof.pos[1]=getdouble(&data);
    raw->gsof.pos[2]=getdouble(&data);
}
/* decode gsof message type 09 record-----------------------------------------*/
static void decode_type09(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type09 :\n");

    raw->gsof.dop[0]=getdouble(&data);
    raw->gsof.dop[1]=getdouble(&data);
    raw->gsof.dop[2]=getdouble(&data);
    raw->gsof.dop[3]=getdouble(&data);
}
/* decode gsof message type 01 record-----------------------------------------*/
static void decode_type01(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type01 :\n");
    int week; double sow;

    raw->gsof.ns=data[6];
    if ((data[8]&0x00)==0x00) raw->gsof.solq=SOLQ_SINGLE;
    if ((data[8]&0x01)==0x01) raw->gsof.solq=SOLQ_DGPS  ;
    if ((data[8]&0x03)==0x03) raw->gsof.solq=SOLQ_FLOAT ;
    if ((data[8]&0x07)==0x07) raw->gsof.solq=SOLQ_FIX   ;

    sow=getu32(&data)/1000.0; week=getu16(&data);
    raw->gsof.t=gpst2time(week,sow);
}
/* decode gsof message type 12 record-----------------------------------------*/
static void decode_type11(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type11 :\n");

    raw->gsof.cov[0]=getfloat(&data);
    raw->gsof.cov[1]=getfloat(&data);
    raw->gsof.cov[2]=getfloat(&data);
    raw->gsof.cov[3]=getfloat(&data);
    raw->gsof.cov[4]=getfloat(&data);
    raw->gsof.cov[5]=getfloat(&data);
    raw->gsof.cov[6]=getfloat(&data);
    raw->gsof.cov[7]=getfloat(&data);
}
/* decode gsof message type 12 record-----------------------------------------*/
static void decode_type12(raw_t *raw, unsigned char* data)
{
    trace(3,"decode_type12 :\n");

    raw->gsof.sig[0]=getfloat(&data);
    raw->gsof.sig[1]=getfloat(&data); /* east */
    raw->gsof.sig[2]=getfloat(&data); /* north */
    raw->gsof.sig[3]=getfloat(&data); /* east-north */
    raw->gsof.sig[4]=getfloat(&data); /* up */
    getfloat(&data); getfloat(&data); getfloat(&data);
    raw->gsof.sig[5]=getfloat(&data); /* unit variance */
}
/* decode gsof message type 38 record-----------------------------------------*/
static void decode_type38(raw_t *raw, unsigned char *data)
{
    trace(3,"decode_type38 :\n");
}
/* decode gsof message type 41 record------------------------------------------*/
static void decode_type41(raw_t *raw, unsigned char *data)
{
    trace(3,"decode_type41 :\n");

    int week; double sow; sow=getu32(&data)/1000.0; week=getu16(&data);
    raw->gsof.base.t      =gpst2time(week,sow);
    raw->gsof.base.pos[0] =getdouble(&data);
    raw->gsof.base.pos[1] =getdouble(&data);
    raw->gsof.base.pos[2] =getdouble(&data);
    raw->gsof.base.quality=data[30];
}
/* decode gsof message type 06 record------------------------------------------*/
static void decode_type06(raw_t *raw, unsigned char *data)
{
    trace(3,"decode_type06 :\n");

    raw->gsof.delta[0]=getdouble(&data);
    raw->gsof.delta[1]=getdouble(&data);
    raw->gsof.delta[2]=getdouble(&data);
}
/* decode gsof message type 08 record------------------------------------------*/
static void decode_type08(raw_t *raw, unsigned char *data)
{
    trace(3,"decode_type08 :\n");

    double s,head;
    raw->gsof.velf=data[0]; data++;
    s   =getfloat(&data);
    head=getfloat(&data);
    raw->gsof.vel[0]=s*sin(head); /* east */
    raw->gsof.vel[1]=s*cos(head); /* north */
    raw->gsof.vel[2]=getfloat(&data);  /* up */
}
/* convert gsof message to solution type--------------------------------------*/
static void gsof2sol(const gsof_t *gsof,sol_t *sol)
{
    int i;
    double Cne[9],std[3],e[3];

    trace(3,"gsof2sol:\n");

    sol->time=gsof->t;

    pos2ecef(gsof->llh,sol->rr);
    ned2xyz (gsof->llh,Cne);
    matmul("NN",3,1,3,1.0,Cne,gsof->vel,0.0,sol->rr+3);

    std[0]=gsof->sig[1];
    std[1]=gsof->sig[2];
    std[2]=gsof->sig[4];

    enu2ecef(gsof->llh,std,e);
    for (i=0;i<3;i++) sol->qr[i]=SQR((float)e[i]);

    sol->ns  =(unsigned char)gsof->ns;
    sol->stat=(unsigned char)gsof->solq;
}
/* check gsof header----------------------------------------------------------*/
static int chkhead(raw_t *raw)
{
    return raw->buff[0]==0x02&&raw->buff[2]==0x40;
}
/* check sum of message ------------------------------------------------------*/
static int chksum(raw_t *raw)
{
    unsigned char cs=0;
    int i;
    for (i=1;i<raw->len+4;i++) cs+=raw->buff[i];
    return cs%256==raw->buff[raw->len+4]&&raw->buff[raw->len+5]==0x03;
}
/* decode gsof message--------------------------------------------------------*/
static void decode_gsof(raw_t *raw)
{
    int type,len,cb=0,nt=GSOFHEADLEN+GSOFMIDLEN;

    trace(3,"decode_gsof: type=0x%02X len=%d\n",type,raw->len);
    
    raw->gsof.status=U1(raw->buff+1);
    raw->gsof.no    =U1(raw->buff+4);

    while (cb<raw->len-GSOFMIDLEN&&raw->len>0) {
        type=U1(raw->buff+nt+cb  ); /* message type */
        len =U1(raw->buff+nt+cb+1); /* message length */
        switch (type) {
            case GSOFTYPE01 : decode_type01(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE02 : decode_type02(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE03 : decode_type03(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE06 : decode_type06(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE08 : decode_type08(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE09 : decode_type09(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE10 : decode_type10(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE11 : decode_type11(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE12 : decode_type12(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE38 : decode_type38(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            case GSOFTYPE41 : decode_type41(raw,raw->buff+GSOFMMLEN+nt+cb); break;
            default:
                trace(2,"no support gsof message \n");
                break;
        }
        cb+=(len+2);
    }
}
/* input gsof raw message in backward------------------------------------------*/
static int decode_gsofb(raw_t *raw, unsigned char data)
{
    if (raw->gsofb.n==0) {
        stream_t *strp=(stream_t*)raw->strp;
        readgsoff(strp->path,&raw->gsofb);

        raw->curb=raw->gsofb.n-1; /* index of current raw data */
    }
    if (raw->curb>=0) {
        gsof2sol(&raw->gsofb.data[raw->curb--],&raw->sol);
        return 6; /* ok */
    }
    return 0; /* fail */
}
/* input gsof raw message int forward------------------------------------------*/
static int decode_gsoff(raw_t *raw, unsigned char data)
{
    raw->buff[raw->nbyte++]=data;

    /* synchronize frame */
    if (raw->nbyte<GSOFHEADLEN) return 0;

    /* check gsof message head */
    if (chkhead(raw)) {
        raw->len=U1(raw->buff+3); /* length of message */
    }
    else {
        raw->nbyte=0; return 0; /* no gsof message */
    }
    if (raw->nbyte<raw->len+6) return 0;

    if (!chksum(raw)) {
        tracet(2,"gsof message checksum error msg=%d\n",raw->buff[2]);
        raw->buff[0]='\0';
        raw->nbyte=0;
        return -1;
    }
    /* decode gsof message */
    decode_gsof(raw);

    /* convert gsof message to solution type */
    gsof2sol(&raw->gsof,&raw->sol);

    raw->buff[0]='\0'; raw->nbyte=0;
    return 6;
}
/* input gsof raw message ------------------------------------------------------
 * args   : raw_t *raw         IO receiver raw data control struct
 *          unsigned char data I  stream data (1 byte)
 * return 1:ok ,0:fail
 * ----------------------------------------------------------------------------*/
extern int input_gsof(raw_t *raw, unsigned char data)
{
    trace(5,"input_gsof: data=%02x\n",data);
    if (raw->dire) return decode_gsofb(raw,data);
    else           return decode_gsoff(raw,data);
}
/* input gsof raw message from file --------------------------------------------
* input next gsof raw message from file
* args   : raw_t  *raw   IO     receiver raw data control struct
*          FILE   *fp    I      file pointer
* return : status(-2: end of file, 1: ok 0: fail)
*-----------------------------------------------------------------------------*/
extern int input_gsoff(raw_t *raw,FILE *fp)
{
    int i,data,ret;

    trace(4,"input_gsoff:\n");

    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_gsof(raw,(unsigned char)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}
/* add gps position/velocity measurement data -------------------------------*/
extern int addgmea(gmeas_t *gmeas, const gmea_t *data)
{
    gmea_t *obs_data;

    if (gmeas->nmax<=gmeas->n) {
        if (gmeas->nmax<=0) gmeas->nmax=NPOS*2; else gmeas->nmax*=2;
        if (!(obs_data=(gmea_t *)realloc(gmeas->data,sizeof(gmea_t)*gmeas->nmax))) {
            trace(1,"addgmea: memalloc error n=%dx%d\n",sizeof(gmea_t),gmeas->nmax);
            free(gmeas->data);
            gmeas->data=NULL;
            gmeas->n=gmeas->nmax=0; return -1;
        }
        gmeas->data=obs_data;
    }
    gmeas->data[gmeas->n++]=*data;
    return 1;
}
/* compare gsof message------------------------------------------------------*/
static int cmpgsof(const void *p1, const void *p2)
{
    gsof_t *q1=(gsof_t *)p1,*q2=(gsof_t *)p2;
    double tt=timediff(q1->t,q2->t);
    if (fabs(tt)>DTTOL) return tt<0?-1:1;
}
/* sort and unique gsof observation data --------------------------------------
* sort and unique observation data by time
* args   : gsof_data_t *gsof    IO     observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortgsof(gsof_data_t *gsof)
{
    int i,j,n;

    trace(3,"sortgsof: nobs=%d\n",gsof->n);

    if (gsof->n<=0) return 0;

    qsort(gsof->data,gsof->n,sizeof(gsof_t),cmpgsof);

    /* delete duplicated data */
    for (i=j=0;i<gsof->n;i++) {
        if (timediff(gsof->data[i].t,gsof->data[j].t)!=0.0) {
            gsof->data[++j]=gsof->data[i];
        }
    }
    gsof->n=j+1;

    for (i=n=0;i<gsof->n;i=j,n++) {
        for (j=i+1;j<gsof->n;j++) {
            if (timediff(gsof->data[j].t,gsof->data[i].t)>DTTOL) break;
        }
    }
    return n;
}
/* input gsof raw message from file --------------------------------------------
* input next gw10 raw message from file
* args   : char   *file       I  input gsof file
*          gsof_data_t *gsof  O  output gsof observation data (NULL: no output)
* return : status(1: ok 0: fail,-2 :end of file)
*-----------------------------------------------------------------------------*/
extern int readgsoff(const char *file,gsof_data_t *gsof)
{
    trace(3,"readgsoff :\n");

    int data,siz;
    raw_t raw={{0}};
    FILE *fp;

    if (!(fp=fopen(file,"r"))) {
        trace(2,"gsof file open error \n");
        return 0;
    }
    /* read gsof message from file */
    while (1) {
        if ((data=fgetc(fp))==EOF) break;
        if ((input_gsof(&raw,(unsigned char)data))) {

            if (gsof->n>=gsof->nmax) {
                trace(5,"readgsoff: gsof->n=%d nmax=%d\n",gsof->n,gsof->nmax);

                siz=sizeof(gsof_t)*(gsof->nmax+=4096);
                if (!(gsof->data=(gsof_t *)realloc(gsof->data,siz))) {
                    
                    fprintf(stderr,"readgsoff :memory allocation error\n");
                    free(gsof->data);
                    gsof->n=gsof->nmax=0;
                    break;
                }
            }
            gsof->data[gsof->n++]=raw.gsof;
        }
    }
    fclose(fp);
    return gsof->n>0;
}
/* free gsof measurement data------------------------------------------------*/
extern void freegsofdata(gsof_data_t *data)
{
    trace(3,"freegsofdata:\n");
    if (data->data) {
        free(data->data); data->data=NULL; data->n=data->nmax=0;
    }
}




