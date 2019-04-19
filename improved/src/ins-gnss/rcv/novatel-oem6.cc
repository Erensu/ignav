/*------------------------------------------------------------------------------
* notvatel-oem6.cc : NovAtel OEM6 receiver functions
*
* reference :
*     [1] NovAtel, OEM6 Family Firmware Reference Manual
*
* version : $Revision: 1.2 $ $Date: 2008/07/14 00:05:05 $
* history : 2018/05/24 1.0 new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

#define OEM6SYNC1   0xAA        /* oem4 message start sync code 1 */
#define OEM6SYNC2   0x44        /* oem4 message start sync code 2 */
#define OEM6SYNC3   0x12        /* oem4 message start sync code 3 */
#define OEM6HLEN    28          /* oem4 message header length (bytes) */

#define ID_QZSSIONUTC 1347      /* message id: oem6 qzss ion/utc parameters */
#define ID_QZSSRAWEPHEM 1330    /* message id: oem6 qzss raw ephemeris */
#define ID_QZSSRAWSUBFRAME 1331 /* message id: oem6 qzss raw subframe */
#define ID_RAWSBASFRAME 973     /* message id: oem6 raw sbas frame */
#define ID_GALEPHEMERIS 1122    /* message id: oem6 decoded galileo ephemeris */
#define ID_GALALMANAC 1120      /* message id: oem6 decoded galileo almanac */
#define ID_GALCLOCK 1121        /* message id: oem6 galileo clockinformation */
#define ID_GALIONO  1127        /* message id: oem6 decoded galileo iono corrections */
#define ID_GALFNAVRAWPAGE 1413  /* message id: oem6 raw galileo f/nav paga data */
#define ID_GALINAVRAWWORD 1414  /* message id: oem6 raw galileo i/nav word data */
#define ID_RAWCNAVFRAME 1066    /* message id: oem6 raw cnav frame data */
#define ID_BDSEPHEMERIS 1696    /* message id: oem6 decoded bds ephemeris */
#define ID_BESTPOS 42           /* message id: oem6 decoded best position */
#define ID_BESTVEL 99           /* message id: oem6 decoded best velocity */
#define ID_BESTXYZ 241          /* message id: oem6 decoded best available cartesian position and velocity */
#define ID_HEAD    971          /* message id: oem6 decoded heading information */
#define ID_HEAD2   1335         /* message id: oem6 decoded heading information with multiple rovers */
#define ID_HEADRATE 1698        /* message id: oem6 decoded heading rate information */

/* get fields (little-endian) ------------------------------------------------*/
#define U1(p) (*((unsigned char *)(p)))
#define I1(p) (*((signed char *)(p)))
static unsigned short U2(unsigned char *p) {unsigned short u; memcpy(&u,p,2); return u;}
static unsigned int   U4(unsigned char *p) {unsigned int   u; memcpy(&u,p,4); return u;}
static int            I4(unsigned char *p) {int            i; memcpy(&i,p,4); return i;}
static float          R4(unsigned char *p) {float          r; memcpy(&r,p,4); return r;}
static double         R8(unsigned char *p) {double         r; memcpy(&r,p,8); return r;}

/* adjust weekly rollover of gps time ----------------------------------------*/
static gtime_t adjweek(gtime_t time, double tow)
{
    double tow_p;
    int week;
    tow_p=time2gpst(time,&week);
    if      (tow<tow_p-302400.0) tow+=604800.0;
    else if (tow>tow_p+302400.0) tow-=604800.0;
    return gpst2time(week,tow);
}
/* ura value (m) to ura index ------------------------------------------------*/
static int uraindex(double value)
{
    static const double ura_eph[]={
            2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
            3072.0,6144.0,0.0
    };
    int i;
    for (i=0;i<15;i++) if (ura_eph[i]>=value) break;
    return i;
}
/* decode qzss rawephemb -----------------------------------------------------*/
static int decode_qzssrawephemb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN,*q;
    eph_t eph={0};
    char *msg;
    int i,prn,id,sat;

    trace(3,"decode_qzssrawephemb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+44) {
        trace(2,"oem4 qzssrawephemb length error: len=%d\n",raw->len);
        return -1;
    }
    prn=U4(p);
    id =U4(p+4);

    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," prn=%3d id=%d",prn,id);
    }
    if (!(sat=satno(SYS_QZS,prn))) {
        trace(2,"oem4 qzssrawephemb satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (id<1||3<id) return 0;

    q=raw->subfrm[sat-1]+(id-1)*30;
    for (i=0;i<30;i++) *q++=p[8+i];

    if (id<3) return 0;
    if (decode_frame(raw->subfrm[sat-1]   ,&eph,NULL,NULL,NULL,NULL)!=1||
        decode_frame(raw->subfrm[sat-1]+30,&eph,NULL,NULL,NULL,NULL)!=2||
        decode_frame(raw->subfrm[sat-1]+60,&eph,NULL,NULL,NULL,NULL)!=3) {
        return 0;
    }
    if (!strstr(raw->opt,"-EPHALL")) {
        if (eph.iodc==raw->nav.eph[sat-1].iodc&&
            eph.iode==raw->nav.eph[sat-1].iode) return 0; /* unchanged */
    }
    eph.sat=sat;
    raw->nav.eph[sat-1]=eph;
    raw->ephsat=sat;
    trace(4,"decode_qzssrawephemb: sat=%2d\n",sat);
    return 2;
}
/* decode rawwaasframeb ------------------------------------------------------*/
static int decode_rawwaasframeb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    int i,prn;

    trace(3,"decode_rawwaasframeb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+48) {
        trace(2,"oem4 rawwaasframeb length error: len=%d\n",raw->len);
        return -1;
    }
    prn=U4(p+4);

    if (MINPRNQZS_S<=prn&&prn<=MAXPRNQZS_S) {
        prn+=10; /* QZSS SAIF PRN -> QZSS PRN */
    }
    else if (prn<MINPRNSBS||MAXPRNSBS<prn) return 0;

    raw->sbsmsg.tow=(int)time2gpst(raw->time,&raw->sbsmsg.week);
    raw->sbsmsg.prn=prn;
    for (i=0,p+=12;i<29;i++,p++) raw->sbsmsg.msg[i]=*p;
    return 3;
}
/* decode rawsbasframeb ------------------------------------------------------*/
static int decode_rawsbasframeb(raw_t *raw)
{
    trace(3,"decode_rawsbasframeb: len=%d\n",raw->len);

    /* format same as rawwaasframeb */
    return decode_rawwaasframeb(raw);
}
/* decode qzss rawsubframeb --------------------------------------------------*/
static int decode_qzssrawsubframeb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    eph_t eph={0};
    char *msg;
    int prn,sat;

    trace(3,"decode_qzssrawsubframeb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+44) {
        trace(2,"oem4 qzssrawsubframeb length error: len=%d\n",raw->len);
        return -1;
    }
    prn=U4(p);

    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," prn=%3d",prn);
    }
    if (!(sat=satno(SYS_QZS,prn))) {
        trace(2,"oem4 qzssrawephemb satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (decode_frame(p+12,&eph,NULL,NULL,NULL,NULL)!=1||
        decode_frame(p+42,&eph,NULL,NULL,NULL,NULL)!=2||
        decode_frame(p+72,&eph,NULL,NULL,NULL,NULL)!=3) {
        return 0;
    }
    if (!strstr(raw->opt,"-EPHALL")) {
        if (eph.iodc==raw->nav.eph[sat-1].iodc&&
            eph.iode==raw->nav.eph[sat-1].iode) return 0; /* unchanged */
    }
    eph.sat=sat;
    raw->nav.eph[sat-1]=eph;
    raw->ephsat=sat;
    trace(4,"decode_qzssrawsubframeb: sat=%2d\n",sat);
    return 2;
}
/* decode qzssionutcb --------------------------------------------------------*/
static int decode_qzssionutcb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    int i;

    trace(3,"decode_qzssionutcb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+108) {
        trace(2,"oem4 qzssionutcb length error: len=%d\n",raw->len);
        return -1;
    }
    for (i=0;i<8;i++) raw->nav.ion_qzs[i]=R8(p+i*8);
    raw->nav.utc_qzs[0]=R8(p+72);
    raw->nav.utc_qzs[1]=R8(p+80);
    raw->nav.utc_qzs[2]=U4(p+68);
    raw->nav.utc_qzs[3]=U4(p+64);
    raw->nav.leaps =I4(p+96);
    return 9;
}
/* decode galephemerisb ------------------------------------------------------*/
static int decode_galephemerisb(raw_t *raw)
{
    eph_t eph={0};
    unsigned char *p=raw->buff+OEM6HLEN;
    double tow,sqrtA,af0_fnav,af1_fnav,af2_fnav,af0_inav,af1_inav,af2_inav,tt;
    char *msg;
    int prn,rcv_fnav,rcv_inav,svh_e1b,svh_e5a,svh_e5b,dvs_e1b,dvs_e5a,dvs_e5b;
    int toc_fnav,toc_inav,week,sel_nav=0;

    trace(3,"decode_galephemerisb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+220) {
        trace(2,"oem4 galephemrisb length error: len=%d\n",raw->len);
        return -1;
    }
    prn       =U4(p);   p+=4;
    rcv_fnav  =U4(p)&1; p+=4;
    rcv_inav  =U4(p)&1; p+=4;
    svh_e1b   =U1(p)&3; p+=1;
    svh_e5a   =U1(p)&3; p+=1;
    svh_e5b   =U1(p)&3; p+=1;
    dvs_e1b   =U1(p)&1; p+=1;
    dvs_e5a   =U1(p)&1; p+=1;
    dvs_e5b   =U1(p)&1; p+=1;
    eph.sva   =U1(p);   p+=1+1; /* SISA */
    eph.iode  =U4(p);   p+=4;   /* IODNav */
    eph.toes  =U4(p);   p+=4;
    sqrtA     =R8(p);   p+=8;
    eph.deln  =R8(p);   p+=8;
    eph.M0    =R8(p);   p+=8;
    eph.e     =R8(p);   p+=8;
    eph.omg   =R8(p);   p+=8;
    eph.cuc   =R8(p);   p+=8;
    eph.cus   =R8(p);   p+=8;
    eph.crc   =R8(p);   p+=8;
    eph.crs   =R8(p);   p+=8;
    eph.cic   =R8(p);   p+=8;
    eph.cis   =R8(p);   p+=8;
    eph.i0    =R8(p);   p+=8;
    eph.idot  =R8(p);   p+=8;
    eph.OMG0  =R8(p);   p+=8;
    eph.OMGd  =R8(p);   p+=8;
    toc_fnav  =U4(p);   p+=4;
    af0_fnav  =R8(p);   p+=8;
    af1_fnav  =R8(p);   p+=8;
    af2_fnav  =R8(p);   p+=8;
    toc_inav  =U4(p);   p+=4;
    af0_inav  =R8(p);   p+=8;
    af1_inav  =R8(p);   p+=8;
    af2_inav  =R8(p);   p+=8;
    eph.tgd[0]=R8(p);   p+=8; /* BGD: E5A-E1 (s) */
    eph.tgd[1]=R8(p);         /* BGD: E5B-E1 (s) */
    eph.iodc  =eph.iode;
    eph.svh   =(svh_e5b<<7)|(dvs_e5b<<6)|(svh_e5a<<4)|(dvs_e5a<<3)|
               (svh_e1b<<1)|dvs_e1b;

    /* ephemeris selection (0:INAV,1:FNAV) */
    if      (strstr(raw->opt,"-GALINAV")) sel_nav=0;
    else if (strstr(raw->opt,"-GALFNAV")) sel_nav=1;
    else if (!rcv_inav&&rcv_fnav) sel_nav=1;

    eph.A     =sqrtA*sqrtA;
    eph.f0    =sel_nav?af0_fnav:af0_inav;
    eph.f1    =sel_nav?af1_fnav:af1_inav;
    eph.f2    =sel_nav?af2_fnav:af2_inav;
    eph.code  =sel_nav?2:1; /* data source 1:I/NAV E1B,2:F/NAV E5a-I */

    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," prn=%3d iod=%3d toes=%6.0f",prn,eph.iode,eph.toes);
    }
    if (!(eph.sat=satno(SYS_GAL,prn))) {
        trace(2,"oemv galephemeris satellite error: prn=%d\n",prn);
        return -1;
    }
    tow=time2gpst(raw->time,&week);
    eph.week=week; /* gps-week = gal-week */
    eph.toe=gpst2time(eph.week,eph.toes);

    /* for week-handover problem */
    tt=timediff(eph.toe,raw->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt> 302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=adjweek(eph.toe,sel_nav?toc_fnav:toc_inav);
    eph.ttr=adjweek(eph.toe,tow);

    if (!strstr(raw->opt,"-EPHALL")) {
        if (raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].code==eph.code) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* decode galalmanacb --------------------------------------------------------*/
static int decode_galalmanacb(raw_t *raw)
{
    alm_t alm={0};
    unsigned char *p=raw->buff+OEM6HLEN;
    double dsqrtA,sqrtA=sqrt(29601297.0);
    int prn,rcv_fnav,rcv_inav,svh_e1b,svh_e5a,svh_e5b,ioda;

    trace(3,"decode_galalmanacb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+100) {
        trace(2,"oem4 galephemrisb length error: len=%d\n",raw->len);
        return -1;
    }
    prn     =U4(p);   p+=4;
    rcv_fnav=U4(p)&1; p+=4;
    rcv_inav=U4(p)&1; p+=4;
    svh_e1b =U1(p)&3; p+=1;
    svh_e5a =U1(p)&3; p+=1;
    svh_e5b =U1(p)&3; p+=1+1;
    ioda    =U4(p);   p+=4;
    alm.week=U4(p);   p+=4; /* gst week */
    alm.toas=U4(p);   p+=4;
    alm.e   =R8(p);   p+=8;
    alm.OMGd=R8(p);   p+=8;
    alm.OMG0=R8(p);   p+=8;
    alm.omg =R8(p);   p+=8;
    alm.M0  =R8(p);   p+=8;
    alm.f0  =R8(p);   p+=8;
    alm.f1  =R8(p);   p+=8;
    dsqrtA  =R8(p);   p+=8;
    alm.i0  =(R8(p)+56.0)*D2R;
    alm.svh =(svh_e5b<<7)|(svh_e5a<<4)|(svh_e1b<<1);
    alm.A   =(sqrtA+dsqrtA)*(sqrtA+dsqrtA);

    if (!(alm.sat=satno(SYS_GAL,prn))) {
        trace(2,"oemv galalmanac satellite error: prn=%d\n",prn);
        return -1;
    }
    alm.toa=gst2time(alm.week,alm.toas);
    raw->nav.alm[alm.sat-1]=alm;
    return 0;
}
/* decode galclockb ----------------------------------------------------------*/
static int decode_galclockb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    double a0,a1,a0g,a1g;
    int leaps,tot,wnt,wnlsf,dn,dtlsf,t0g,wn0g;

    trace(3,"decode_galclockb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+64) {
        trace(2,"oem4 galclockb length error: len=%d\n",raw->len);
        return -1;
    }
    a0   =R8(p); p+=8;
    a1   =R8(p); p+=8;
    leaps=I4(p); p+=4;
    tot  =U4(p); p+=4;
    wnt  =U4(p); p+=4;
    wnlsf=U4(p); p+=4;
    dn   =U4(p); p+=4;
    dtlsf=U4(p); p+=4;
    a0g  =R8(p); p+=8;
    a1g  =R8(p); p+=8;
    t0g  =U4(p); p+=4;
    wn0g =U4(p);

    raw->nav.utc_gal[0]=a0;
    raw->nav.utc_gal[1]=a1;
    raw->nav.utc_gal[2]=tot; /* utc reference tow (s) */
    raw->nav.utc_gal[3]=wnt; /* utc reference week */
    return 9;
}
/* decode galionob -----------------------------------------------------------*/
static int decode_galionob(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    double ai[3];
    int i,sf[5];

    trace(3,"decode_galionob: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+29) {
        trace(2,"oem4 galionob length error: len=%d\n",raw->len);
        return -1;
    }
    ai[0]=R8(p); p+=8;
    ai[1]=R8(p); p+=8;
    ai[2]=R8(p); p+=8;
    sf[0]=U1(p); p+=1;
    sf[1]=U1(p); p+=1;
    sf[2]=U1(p); p+=1;
    sf[3]=U1(p); p+=1;
    sf[4]=U1(p);

    for (i=0;i<3;i++) raw->nav.ion_gal[i]=ai[i];
    return 9;
}
/* decode galfnavrawpageb ----------------------------------------------------*/
static int decode_galfnavrawpageb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    unsigned char buff[27];
    int i,sigch,satid,page;

    trace(3,"decode_galfnavrawpageb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+35) {
        trace(2,"oem4 galfnavrawpageb length error: len=%d\n",raw->len);
        return -1;
    }
    sigch=U4(p); p+=4;
    satid=U4(p); p+=4;
    for (i=0;i<27;i++) {
        buff[i]=U1(p); p+=1;
    }
    page=getbitu(buff,0,6);

    trace(3,"%s E%2d FNAV     (%2d) ",time_str(raw->time,0),satid,page);
    traceb(3,buff,27);

    return 0;
}
/* decode galinavrawwordb ----------------------------------------------------*/
static int decode_galinavrawwordb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    unsigned char buff[16];
    gtime_t time=raw->time;
    char *sig;
    int i,sigch,satid,sigtype,type,week=0,tow=0;

    trace(3,"decode_galinavrawwordb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+28) {
        trace(2,"oem4 galinavrawwordb length error: len=%d\n",raw->len);
        return -1;
    }
    sigch  =U4(p); p+=4;
    satid  =U4(p); p+=4;
    sigtype=U4(p); p+=4;

    switch (sigtype) {
        case 10433: sig="E1 "; break;
        case 10466: sig="E5A"; break;
        case 10499: sig="E5B"; break;
        default: sig="???"   ; break;
    }
    for (i=0;i<16;i++) {
        buff[i]=U1(p); p+=1;
    }
    type=getbitu(buff,0,6);
    if (type==0&&getbitu(buff,6,2)==2) {
        week=getbitu(buff, 96,12); /* gst week */
        tow =getbitu(buff,108,20);
        time=gst2time(week,tow);
    }
    trace(3,"%s E%2d INAV-%s (%2d) ",time_str(time,0),satid,sig,type);
    traceb(3,buff,16);

    return 0;
}
/* decode rawcnavframeb ------------------------------------------------------*/
static int decode_rawcnavframeb(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    unsigned char buff[38];
    int i,sigch,prn,frmid;

    trace(3,"decode_rawcnavframeb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+50) {
        trace(2,"oem4 rawcnavframeb length error: len=%d\n",raw->len);
        return -1;
    }
    sigch=U4(p); p+=4;
    prn  =U4(p); p+=4;
    frmid=U4(p); p+=4;

    for (i=0;i<38;i++) {
        buff[i]=U1(p); p+=1;
    }
    trace(3,"%s PRN=%3d FRMID=%2d ",time_str(raw->time,0),prn,frmid);
    traceb(3,buff,38);

    return 0;
}
/* decode bdsephemerisb ------------------------------------------------------*/
static int decode_bdsephemerisb(raw_t *raw)
{
    eph_t eph={0};
    unsigned char *p=raw->buff+OEM6HLEN;
    double ura,sqrtA;
    char *msg;
    int prn,toc;

    trace(3,"decode_bdsephemerisb: len=%d\n",raw->len);

    if (raw->len<OEM6HLEN+196) {
        trace(2,"oem4 bdsephemrisb length error: len=%d\n",raw->len);
        return -1;
    }
    prn       =U4(p);   p+=4;
    eph.week  =U4(p);   p+=4;
    ura       =R8(p);   p+=8;
    eph.svh   =U4(p)&1; p+=4;
    eph.tgd[0]=R8(p);   p+=8; /* TGD1 for B1 (s) */
    eph.tgd[1]=R8(p);   p+=8; /* TGD2 for B2 (s) */
    eph.iodc  =U4(p);   p+=4; /* AODC */
    toc       =U4(p);   p+=4;
    eph.f0    =R8(p);   p+=8;
    eph.f1    =R8(p);   p+=8;
    eph.f2    =R8(p);   p+=8;
    eph.iode  =U4(p);   p+=4; /* AODE */
    eph.toes  =U4(p);   p+=4;
    sqrtA     =R8(p);   p+=8;
    eph.e     =R8(p);   p+=8;
    eph.omg   =R8(p);   p+=8;
    eph.deln  =R8(p);   p+=8;
    eph.M0    =R8(p);   p+=8;
    eph.OMG0  =R8(p);   p+=8;
    eph.OMGd  =R8(p);   p+=8;
    eph.i0    =R8(p);   p+=8;
    eph.idot  =R8(p);   p+=8;
    eph.cuc   =R8(p);   p+=8;
    eph.cus   =R8(p);   p+=8;
    eph.crc   =R8(p);   p+=8;
    eph.crs   =R8(p);   p+=8;
    eph.cic   =R8(p);   p+=8;
    eph.cis   =R8(p);
    eph.A     =sqrtA*sqrtA;
    eph.sva   =uraindex(ura);

    if (raw->outtype) {
        msg=raw->msgtype+strlen(raw->msgtype);
        sprintf(msg," prn=%3d iod=%3d toes=%6.0f",prn,eph.iode,eph.toes);
    }
    if (!(eph.sat=satno(SYS_CMP,prn))) {
        trace(2,"oemv bdsephemeris satellite error: prn=%d\n",prn);
        return -1;
    }
    eph.toe=bdt2gpst(bdt2time(eph.week,eph.toes)); /* bdt -> gpst */
    eph.toc=bdt2gpst(bdt2time(eph.week,toc));      /* bdt -> gpst */
    eph.ttr=raw->time;

    if (!strstr(raw->opt,"-EPHALL")) {
        if (timediff(raw->nav.eph[eph.sat-1].toe,eph.toe)==0.0&&
            raw->nav.eph[eph.sat-1].iode==eph.iode&&
            raw->nav.eph[eph.sat-1].iodc==eph.iodc) return 0; /* unchanged */
    }
    raw->nav.eph[eph.sat-1]=eph;
    raw->ephsat=eph.sat;
    return 2;
}
/* decode best position-------------------------------------------------------*/
static int decode_bestpos(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    char sta_id[8]={0};
    sol_t sol={0};
    int i,sol_stat,pos_stat,datum_id,ns,ns_sol,ns_sol_l1,ns_sol_mul;
    double pos[3],cov[9]={0},Cne[9],cov_xyz[9];
    float undulation,sig[3],diff_age,sol_age;

    trace(3,"decode_bestpos: len=%d\n",raw->len);

    sol_stat=U4(p); p+=4;
    pos_stat=U4(p); p+=4;

    pos[0]  =R8(p); p+=8;
    pos[1]  =R8(p); p+=8;
    pos[2]  =R8(p); p+=8;

    undulation=R4(p); p+=4;
    datum_id  =U4(p); p+=4;
    sig[0]    =R4(p); p+=4;
    sig[1]    =R4(p); p+=4;
    sig[2]    =R4(p); p+=4;

    memcpy(sta_id,p,sizeof(unsigned char)*4); p+=4;

    diff_age  =R4(p); p+=4;
    sol_age   =R4(p); p+=4;
    
    ns        =U1(p); p+=1;
    ns_sol    =U1(p); p+=1;
    ns_sol_l1 =U1(p); p+=1;
    ns_sol_mul=U1(p);

    if (sol_stat==0) {
        switch (pos_stat) {
            case  0: sol.stat=SOLQ_NONE  ; break;
            case  1:
            case  2:
            case 48:
            case 50: sol.stat=SOLQ_FIX   ; break;
            case 16: sol.stat=SOLQ_SINGLE; break;
            case 17: sol.stat=SOLQ_DGPS  ; break;
            case 18: sol.stat=SOLQ_WAAS  ; break;
            case 19: sol.stat=SOLQ_PROP  ; break;
            case 32:
            case 33:
            case 34: sol.stat=SOLQ_FLOAT ; break;
            case 68:
            case 69:
            case 77:
            case 78: sol.stat=SOLQ_PPP   ; break;

            case 70:
            case 71:
            case 72: sol.stat=SOLQ_SINGLE; break;

        }
    }
    else {
        sol.stat=SOLQ_NONE;
    }
    sol.age =sol_age; sol.type=0;
    sol.ns  =ns_sol;
    sol.time=raw->time;

    pos[0]*=D2R; pos[1]*=D2R;
    pos2ecef(pos,sol.rr);

    cov[0]=SQR(sig[0]);
    cov[4]=SQR(sig[1]);
    cov[8]=SQR(sig[2]);
    ned2xyz(pos,Cne);

    matmul33("NNT",Cne,cov,Cne,3,3,3,3,cov_xyz);
    for (i=0;i<3;i++) {
        sol.qr[i]=(float)cov_xyz[i+i*3];
    }
    sol.qr[3]=(float)cov_xyz[0+1*3];
    sol.qr[4]=(float)cov_xyz[1+2*3];
    sol.qr[5]=(float)cov_xyz[2+1*3];
    raw->sol=sol;
    return 6;
}
/* decode best velocity-------------------------------------------------------*/
static int decode_bestvel(raw_t *raw)
{
    trace(3,"decode_bestvel: len=%d\n",raw->len);
    return 0;
}
/* decode best available cartesian position and velocity----------------------*/
static int decode_bestxyz(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    sol_t sol={0};
    int sol_stat,pos_stat,ns,ns_sol,ns_sol_l1,ns_sol_mul;
    int vel_stat,vel_type;
    float v_latency;

    trace(3,"decode_bestxyz: len=%d\n",raw->len);

    sol_stat =U4(p); p+=4; /* solution status */
    pos_stat =U4(p); p+=4; /* position status */

    sol.rr[0]=R4(p); p+=8;
    sol.rr[1]=R8(p); p+=8;
    sol.rr[2]=R8(p); p+=8;

    sol.qr[0]=R4(p); p+=4;
    sol.qr[1]=R4(p); p+=4;
    sol.qr[2]=R4(p); p+=4;

    vel_stat =U4(p); p+=4;
    vel_type =U4(p); p+=4;

    sol.sol_vel.v[0] =R8(p); p+=8;
    sol.sol_vel.v[1] =R8(p); p+=8;
    sol.sol_vel.v[2] =R8(p); p+=8;

    sol.sol_vel.qv[0]=R4(p); p+=4;
    sol.sol_vel.qv[1]=R4(p); p+=4;
    sol.sol_vel.qv[2]=R4(p);

    p+=4; p+=4; /* base station identification */
    v_latency=R4(p); p+=4;

    sol.age   =R4(p); p+=4;
    ns        =U1(p); p+=1;
    ns_sol    =U1(p); p+=1;
    ns_sol_l1 =U1(p); p+=1;
    ns_sol_mul=U1(p);

    if (sol_stat==0) {
        switch (pos_stat) {
            case  0: sol.stat=SOLQ_NONE  ; break;
            case  1:
            case  2:
            case 48:
            case 50: sol.stat=SOLQ_FIX   ; break;
            case 16: sol.stat=SOLQ_SINGLE; break;
            case 17: sol.stat=SOLQ_DGPS  ; break;
            case 18: sol.stat=SOLQ_WAAS  ; break;
            case 19: sol.stat=SOLQ_PROP  ; break;
            case 32:
            case 33:
            case 34: sol.stat=SOLQ_FLOAT ; break;
            case 68:
            case 69:
            case 77:
            case 78: sol.stat=SOLQ_PPP   ; break;

            case 70:
            case 71:
            case 72: sol.stat=SOLQ_SINGLE; break;

        }
    }
    else {
        sol.stat=SOLQ_NONE;
    }
    sol.ns  =ns_sol;
    sol.time=raw->time;
    sol.type=0;

    if (vel_stat==0) {
        switch (vel_type) {
            case 0 : sol.sol_vel.type=0; break;
            case 8 :
            case 17:
            case 18: sol.sol_vel.type=1; break;
            default: sol.sol_vel.type=2; break;
        }
    }
    else {
        sol.sol_vel.type=0;
    }
    sol.sol_vel.time=timeadd(raw->time,-v_latency);
    raw->sol=sol;
    return 6;
}
/* decode heading information-------------------------------------------------*/
static int decode_head(raw_t *raw)
{
    unsigned char *p=raw->buff+OEM6HLEN;
    int stat,type;
    float len;
    pose_meas_t pose={0};

    trace(3,"decode_head: len=%d\n",raw->len);

    stat=U4(p); p+=4;
    type=U4(p); p+=4;
    len =R4(p); p+=4;

    pose.rpy[2]=R4(p)*D2R; p+=4;
    pose.rpy[1]=R4(p)*D2R;

    p+=4; p+=4; /* reserved */
    
    pose.var[2]=SQR(R4(p)*D2R); p+=4;
    pose.var[1]=SQR(R4(p)*D2R);

    if (stat==0) {
        switch (type) {
            case  0: pose.stat=SOLQ_NONE  ; break;
            case  1:
            case  2:
            case 48:
            case 50: pose.stat=SOLQ_FIX   ; break;
            case 16: pose.stat=SOLQ_SINGLE; break;
            case 17: pose.stat=SOLQ_DGPS  ; break;
            case 18: pose.stat=SOLQ_WAAS  ; break;
            case 19: pose.stat=SOLQ_PROP  ; break;
            case 32:
            case 33:
            case 34: pose.stat=SOLQ_FLOAT ; break;
            case 68:
            case 69:
            case 77:
            case 78: pose.stat=SOLQ_PPP   ; break;

            case 70:
            case 71:
            case 72: pose.stat=SOLQ_SINGLE; break;

        }
    }
    else {
        pose.stat=SOLQ_NONE;
    }
    pose.time=raw->time;
    pose.len =len;
    pose.type=POSE_DUAL_ANT;
#if 1
    rpy2dcm(pose.rpy,pose.C); /* convert to dcm */
#endif
    raw->pose=pose;
    return 34;
}
/* decode heading information with multiple rovers----------------------------*/
static int decode_head2(raw_t *raw)
{
    trace(3,"decode_head: len=%d\n",raw->len);

    return decode_head(raw);
}
/* decode heading rate information--------------------------------------------*/
static int decode_headrate(raw_t *raw)
{
    trace(3,"decode_headrate:\n");
    return 0;
}
/* decode oem6 message -------------------------------------------------------*/
static int decode_oem6_sol(raw_t *raw)
{
    double tow;
    int msg,week,type=U2(raw->buff+4);

    trace(3,"decode_oem6_sol: type=%3d len=%d\n",type,raw->len);

    /* check crc32 */
    if (rtk_crc32(raw->buff,raw->len)!=U4(raw->buff+raw->len)) {
        trace(2,"oem6 crc error: type=%3d len=%d\n",type,raw->len);
        return -1;
    }
    msg =(U1(raw->buff+6)>>4)&0x3;
    if (!(week=U2(raw->buff+14))) {
        return -1;
    }
    week=adjgpsweek(week);
    tow =U4(raw->buff+16)*0.001;
    raw->time=gpst2time(week,tow);

    if (raw->outtype) {
        sprintf(raw->msgtype,"OEM6 %4d (%4d): msg=%d %s",type,raw->len,msg,
                time_str(gpst2time(week,tow),2));
    }
    if (msg!=0) return 0; /* message type: 0=binary,1=ascii */

    switch (type) {
        case ID_BESTPOS: return decode_bestpos(raw);
        case ID_BESTVEL: return decode_bestvel(raw);
        case ID_BESTXYZ: return decode_bestxyz(raw);
    }
    return 0;
}
/* decode oem6 pose message -------------------------------------------------*/
static int decode_oem6_pose(raw_t *raw)
{
    double tow;
    int msg,week,type=U2(raw->buff+4);

    trace(3,"decode_oem6_pose: type=%3d len=%d\n",type,raw->len);

    /* check crc32 */
    if (rtk_crc32(raw->buff,raw->len)!=U4(raw->buff+raw->len)) {
        trace(2,"oem6 crc error: type=%3d len=%d\n",type,raw->len);
        return -1;
    }
    msg =(U1(raw->buff+6)>>4)&0x3;
    if (!(week=U2(raw->buff+14))) {
        return -1;
    }
    week=adjgpsweek(week);
    tow =U4(raw->buff+16)*0.001;
    raw->time=gpst2time(week,tow);

    if (raw->outtype) {
        sprintf(raw->msgtype,"OEM6 %4d (%4d): msg=%d %s",type,raw->len,msg,
                time_str(gpst2time(week,tow),2));
    }
    if (msg!=0) return 0; /* message type: 0=binary,1=ascii */

    switch (type) {
        case ID_HEAD    : return decode_head    (raw);
        case ID_HEAD2   : return decode_head2   (raw);
        case ID_HEADRATE: return decode_headrate(raw);
    }
    return 0;
}
/* decode oem6 raw observation data-------------------------------------------*/
extern int decode_oem6_raw(raw_t *raw)
{
    double tow;
    int msg,week,type=U2(raw->buff+4);

    trace(3,"decode_oem6_raw: type=%3d len=%d\n",type,raw->len);

    /* check crc32 */
    if (rtk_crc32(raw->buff,raw->len)!=U4(raw->buff+raw->len)) {
        trace(2,"oem6 crc error: type=%3d len=%d\n",type,raw->len);
        return -1;
    }
    msg =(U1(raw->buff+6)>>4)&0x3;
    if (!(week=U2(raw->buff+14))) {
        return -1;
    }
    week=adjgpsweek(week);
    tow =U4(raw->buff+16)*0.001;
    raw->time=gpst2time(week,tow);

    if (raw->outtype) {
        sprintf(raw->msgtype,"OEM6 %4d (%4d): msg=%d %s",type,raw->len,msg,
                time_str(gpst2time(week,tow),2));
    }
    if (msg!=0) return 0; /* message type: 0=binary,1=ascii */

    switch (type) {
        case ID_RAWSBASFRAME  : return decode_rawsbasframeb  (raw);
        case ID_QZSSRAWEPHEM  : return decode_qzssrawephemb  (raw);
        case ID_QZSSRAWSUBFRAME:return decode_qzssrawsubframeb(raw);
        case ID_QZSSIONUTC    : return decode_qzssionutcb    (raw);
        case ID_GALEPHEMERIS  : return decode_galephemerisb  (raw);
        case ID_GALALMANAC    : return decode_galalmanacb    (raw);
        case ID_GALCLOCK      : return decode_galclockb      (raw);
        case ID_GALIONO       : return decode_galionob       (raw);
        case ID_GALFNAVRAWPAGE: return decode_galfnavrawpageb(raw);
        case ID_GALINAVRAWWORD: return decode_galinavrawwordb(raw);
        case ID_RAWCNAVFRAME  : return decode_rawcnavframeb  (raw);
        case ID_BDSEPHEMERIS  : return decode_bdsephemerisb  (raw);
    }
    return 0;
}
/* sync header ---------------------------------------------------------------*/
static int sync_oem6(unsigned char *buff, unsigned char data)
{
    buff[0]=buff[1]; buff[1]=buff[2]; buff[2]=data;
    return buff[0]==OEM6SYNC1&&buff[1]==OEM6SYNC2&&buff[2]==OEM6SYNC3;
}
/* input oem6 raw data from stream --------------------------------------------
* fetch next novatel oem6 raw data and input a mesasge from stream
* args   : raw_t *raw   IO     receiver raw data control struct
*          unsigned char data I stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 3: input sbas message,
*                  9: input ion/utc parameter
*                  6: pvt solutions for ins/gnss coupled
*                  34: pose measurement data (from camera or dual ant.))
*
* notes  : to specify input options for oem4, set raw->opt to the following
*          option strings separated by spaces.
*
*          -EPHALL : input all ephemerides
*          -GL1P   : select 1P for GPS L1 (default 1C)
*          -GL2X   : select 2X for GPS L2 (default 2W)
*          -RL2C   : select 2C for GLO L2 (default 2P)
*          -EL2C   : select 2C for GAL L2 (default 2C)
*          -GALINAV: use I/NAV for GAL ephemeris
*          -GALFNAV: use F/NAV for GAL ephemeris
*
*-----------------------------------------------------------------------------*/
extern int input_oem6_sol(raw_t *raw, unsigned char data)
{
    trace(5,"input_oem6_sol: data=%02x\n",data);

    /* synchronize frame */
    if (raw->nbyte==0) {
        if (sync_oem6(raw->buff,data)) raw->nbyte=3;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte==10&&(raw->len=U2(raw->buff+8)+OEM6HLEN)>MAXRAWLEN-4) {
        trace(2,"oem6 length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (raw->nbyte<10||raw->nbyte<raw->len+4) return 0;
    raw->nbyte=0;

    /* decode oem6 message */
    return decode_oem6_sol(raw);
}
extern int input_oem6_pose(raw_t *raw, unsigned char data)
{
    trace(5,"input_oem6_pose: data=%02x\n",data);

    /* synchronize frame */
    if (raw->nbyte==0) {
        if (sync_oem6(raw->buff,data)) raw->nbyte=3;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte==10&&(raw->len=U2(raw->buff+8)+OEM6HLEN)>MAXRAWLEN-4) {
        trace(2,"oem6 length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (raw->nbyte<10||raw->nbyte<raw->len+4) return 0;
    raw->nbyte=0;

    /* decode oem6 message */
    return decode_oem6_pose(raw);
}
extern int input_oem6_raw(raw_t *raw, unsigned char data)
{
    trace(5,"input_oem6_raw: data=%02x\n",data);

    /* synchronize frame */
    if (raw->nbyte==0) {
        if (sync_oem6(raw->buff,data)) raw->nbyte=3;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte==10&&(raw->len=U2(raw->buff+8)+OEM6HLEN)>MAXRAWLEN-4) {
        trace(2,"oem6 length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (raw->nbyte<10||raw->nbyte<raw->len+4) return 0;
    raw->nbyte=0;

    /* decode oem6 message */
    return decode_oem6_raw(raw);
}
/* input oem6 raw data from file ----------------------------------------------
* fetch next novatel oem6 raw data and input a message from file
* args   : raw_t  *raw   IO     receiver raw data control struct
*          int    format I      receiver raw data format (STRFMT_???)
*          FILE   *fp    I      file pointer
* return : status(-2: end of file, -1...9: same as above)
*-----------------------------------------------------------------------------*/
extern int input_oem6f_sol(raw_t *raw, FILE *fp)
{
    int i,data;

    trace(4,"input_oem6f_sol:\n");

    /* synchronize frame */
    if (raw->nbyte==0) {
        for (i=0;;i++) {
            if ((data=fgetc(fp))==EOF) return -2;
            if (sync_oem6(raw->buff,(unsigned char)data)) break;
            if (i>=4096) return 0;
        }
    }
    if (fread(raw->buff+3,7,1,fp)<1) return -2;
    raw->nbyte=10;

    if ((raw->len=U2(raw->buff+8)+OEM6HLEN)>MAXRAWLEN-4) {
        trace(2,"oem6 length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (fread(raw->buff+10,raw->len-6,1,fp)<1) return -2;
    raw->nbyte=0;

    /* decode oem6 message */
    return decode_oem6_sol(raw);
}
extern int input_oem6f_pose(raw_t *raw, FILE *fp)
{
    int i,data;

    trace(4,"input_oem6f_pose:\n");

    /* synchronize frame */
    if (raw->nbyte==0) {
        for (i=0;;i++) {
            if ((data=fgetc(fp))==EOF) return -2;
            if (sync_oem6(raw->buff,(unsigned char)data)) break;
            if (i>=4096) return 0;
        }
    }
    if (fread(raw->buff+3,7,1,fp)<1) return -2;
    raw->nbyte=10;

    if ((raw->len=U2(raw->buff+8)+OEM6HLEN)>MAXRAWLEN-4) {
        trace(2,"oem6 length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (fread(raw->buff+10,raw->len-6,1,fp)<1) return -2;
    raw->nbyte=0;

    /* decode oem6 message */
    return decode_oem6_pose(raw);
}
extern int input_oem6f_raw(raw_t *raw, FILE *fp)
{
    int i,data;

    trace(4,"input_oem6f_raw:\n");

    /* synchronize frame */
    if (raw->nbyte==0) {
        for (i=0;;i++) {
            if ((data=fgetc(fp))==EOF) return -2;
            if (sync_oem6(raw->buff,(unsigned char)data)) break;
            if (i>=4096) return 0;
        }
    }
    if (fread(raw->buff+3,7,1,fp)<1) return -2;
    raw->nbyte=10;

    if ((raw->len=U2(raw->buff+8)+OEM6HLEN)>MAXRAWLEN-4) {
        trace(2,"oem6 length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (fread(raw->buff+10,raw->len-6,1,fp)<1) return -2;
    raw->nbyte=0;

    /* decode oem6 message */
    return decode_oem6_raw(raw);
}