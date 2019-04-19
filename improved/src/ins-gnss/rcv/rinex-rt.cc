/*------------------------------------------------------------------------------
* rinex-rt.cc : rinex functions for real-time input stream
*
* reference :
*     [1] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 2.11, December 10, 2007
*     [2] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.00, November 28, 2007
*     [3] IS-GPS-200D, Navstar GPS Space Segment/Navigation User Interfaces,
*         7 March, 2006
*     [4] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 2.12, June 23, 2009
*     [5] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.01, June 22, 2009
*     [6] J.Ray and W.Gurtner, RINEX extentions to handle clock information
*         version 3.02, September 2, 2010
*     [7] RINEX The Receiver Independent Exchange Format Version 3.02,
*         International GNSS Service (IGS), RINEX Working Group and Radio
*         Technical Commission for Maritime Services Special Committee 104
*         (RTCM-SC104), December 10, 2012
*
* version : $Revision:$
* history : 2018/02/13 1.0  new
*-----------------------------------------------------------------------------*/
#include <carvig.h>

/* constants/macros ----------------------------------------------------------*/
#define NUMSYS      6                   /* number of systems */
#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD  1024                /* max head line position */
#define MINFREQ_GLO -7                  /* min frequency number glonass */
#define MAXFREQ_GLO 13                  /* max frequency number glonass */
#define NINCOBS     262144              /* inclimental number of obs data */

static const int navsys[]={             /* satellite systems */
        SYS_GPS,SYS_GLO,SYS_GAL,SYS_QZS,SYS_SBS,SYS_CMP,SYS_IRN,0
};
static const char syscodes[]="GREJSCI"; /* satellite system codes */
static const char obscodes[]="CLDS";    /* obs type codes */
static const char frqcodes[]="1256789"; /* frequency codes */

/* set string without tail space ---------------------------------------------*/
static void setstr(char *dst, const char *src, int n)
{
    char *p=dst;
    const char *q=src;
    while (*q&&q<src+n) *p++=*q++;
    *p--='\0';
    while (p>=dst&&*p==' ') *p--='\0';
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjweek(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-302400.0) return timeadd(t, 604800.0);
    if (tt> 302400.0) return timeadd(t,-604800.0);
    return t;
}
/* adjust time considering week handover -------------------------------------*/
static gtime_t adjday(gtime_t t, gtime_t t0)
{
    double tt=timediff(t,t0);
    if (tt<-43200.0) return timeadd(t, 86400.0);
    if (tt> 43200.0) return timeadd(t,-86400.0);
    return t;
}
/* time string for ver.3 (yyyymmdd hhmmss UTC) -------------------------------*/
static void timestr_rnx(char *str)
{
    gtime_t time;
    double ep[6];
    time=timeget();
    time.sec=0.0;
    time2epoch(time,ep);
    sprintf(str,"%04.0f%02.0f%02.0f %02.0f%02.0f%02.0f UTC",ep[0],ep[1],ep[2],
            ep[3],ep[4],ep[5]);
}
/* satellite to satellite code -----------------------------------------------*/
static int sat2code(int sat, char *code)
{
    int prn;
    switch (satsys(sat,&prn)) {
        case SYS_GPS: sprintf(code,"G%2d",prn-MINPRNGPS+1); break;
        case SYS_GLO: sprintf(code,"R%2d",prn-MINPRNGLO+1); break;
        case SYS_GAL: sprintf(code,"E%2d",prn-MINPRNGAL+1); break;
        case SYS_SBS: sprintf(code,"S%2d",prn-100); break;
        case SYS_QZS: sprintf(code,"J%2d",prn-MINPRNQZS+1); break;
        case SYS_CMP: sprintf(code,"C%2d",prn-MINPRNCMP+1); break;
        case SYS_IRN: sprintf(code,"I%2d",prn-MINPRNIRN+1); break;
        default: return 0;
    }
    return 1;
}
/* decode obs header ---------------------------------------------------------*/
static void decode_obsh(char *buff,int *tsys,char tobs[][MAXOBSTYPE][4],
                        nav_t *nav, sta_t *sta)
{
    /* default codes for unknown code */
    const char *defcodes[]={
            "CWX    ",  /* GPS: L125____ */
            "CC     ",  /* GLO: L12_____ */
            "X XXXX ",  /* GAL: L1_5678_ */
            "CXXX   ",  /* QZS: L1256___ */
            "C X    ",  /* SBS: L1_5____ */
            "X  XX  ",  /* BDS: L1__67__ */
            "  A   A"   /* IRN: L__5___9 */
    };
    double del[3];
    int i,j,k,n,nt,prn,fcn;
    const char *p;
    char *label=buff+60;

    if      (strstr(label,"MARKER NAME"         )) {
        if (sta) setstr(sta->name,buff,60);
    }
    else if (strstr(label,"MARKER NUMBER"       )) { /* opt */
        if (sta) setstr(sta->marker,buff,20);
    }
    else if (strstr(label,"MARKER TYPE"         )) ; /* ver.3 */
    else if (strstr(label,"OBSERVER / AGENCY"   )) ;
    else if (strstr(label,"REC # / TYPE / VERS" )) {
        if (sta) {
            setstr(sta->recsno, buff,   20);
            setstr(sta->rectype,buff+20,20);
            setstr(sta->recver, buff+40,20);
        }
    }
    else if (strstr(label,"ANT # / TYPE"        )) {
        if (sta) {
            setstr(sta->antsno,buff   ,20);
            setstr(sta->antdes,buff+20,20);
        }
    }
    else if (strstr(label,"APPROX POSITION XYZ" )) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) sta->pos[i]=str2num(buff,j,14);
        }
    }
    else if (strstr(label,"ANTENNA: DELTA H/E/N")) {
        if (sta) {
            for (i=0,j=0;i<3;i++,j+=14) del[i]=str2num(buff,j,14);
            sta->del[2]=del[0]; /* h */
            sta->del[0]=del[1]; /* e */
            sta->del[1]=del[2]; /* n */
        }
    }
    else if (strstr(label,"ANTENNA: DELTA X/Y/Z")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: PHASECENTER")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: B.SIGHT XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR AZI")) ; /* opt ver.3 */
    else if (strstr(label,"ANTENNA: ZERODIR XYZ")) ; /* opt ver.3 */
    else if (strstr(label,"CENTER OF MASS: XYZ" )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / # / OBS TYPES" )) { /* ver.3 */
        if (!(p=strchr(syscodes,buff[0]))) {
            trace(2,"invalid system code: sys=%c\n",buff[0]);
            return;
        }
        i=(int)(p-syscodes);
        n=(int)str2num(buff,3,3);
        for (j=nt=0,k=7;j<n;j++,k+=4) {
            if (nt<MAXOBSTYPE-1) setstr(tobs[i][nt++],buff+k,3);
        }
        *tobs[i][nt]='\0';

        /* change beidou B1 code: 3.02 draft -> 3.02/3.03 */
        if (i==5) {
            for (j=0;j<nt;j++) if (tobs[i][j][1]=='2') tobs[i][j][1]='1';
        }
        /* if unknown code in ver.3, set default code */
        for (j=0;j<nt;j++) {
            if (tobs[i][j][2]) continue;
            if (!(p=strchr(frqcodes,tobs[i][j][1]))) continue;
            tobs[i][j][2]=defcodes[i][(int)(p-frqcodes)];
            trace(2,"set default for unknown code: sys=%c code=%s\n",buff[0],
                  tobs[i][j]);
        }
    }
    else if (strstr(label,"WAVELENGTH FACT L1/2")) ; /* opt ver.2 */
    else if (strstr(label,"SIGNAL STRENGTH UNIT")) ; /* opt ver.3 */
    else if (strstr(label,"INTERVAL"            )) ; /* opt */
    else if (strstr(label,"TIME OF FIRST OBS"   )) {
        if      (!strncmp(buff+48,"GPS",3)) *tsys=TSYS_GPS;
        else if (!strncmp(buff+48,"GLO",3)) *tsys=TSYS_UTC;
        else if (!strncmp(buff+48,"GAL",3)) *tsys=TSYS_GAL;
        else if (!strncmp(buff+48,"QZS",3)) *tsys=TSYS_QZS; /* ver.3.02 */
        else if (!strncmp(buff+48,"BDT",3)) *tsys=TSYS_CMP; /* ver.3.02 */
        else if (!strncmp(buff+48,"IRN",3)) *tsys=TSYS_IRN; /* ver.3.03 */
    }
    else if (strstr(label,"TIME OF LAST OBS"    )) ; /* opt */
    else if (strstr(label,"RCV CLOCK OFFS APPL" )) ; /* opt */
    else if (strstr(label,"SYS / DCBS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PCVS APPLIED"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / SCALE FACTOR"  )) ; /* opt ver.3 */
    else if (strstr(label,"SYS / PHASE SHIFTS"  )) ; /* ver.3.01 */
    else if (strstr(label,"GLONASS SLOT / FRQ #")) { /* ver.3.02 */
        if (nav) {
            for (i=0,p=buff+4;i<8;i++,p+=8) {
                if (sscanf(p,"R%2d %2d",&prn,&fcn)<2) continue;
                if (1<=prn&&prn<=MAXPRNGLO) nav->glo_fcn[prn-1]=fcn+8;
            }
        }
    }
    else if (strstr(label,"GLONASS COD/PHS/BIS" )) { /* ver.3.02 */
        if (nav) {
            for (i=0,p=buff;i<4;i++,p+=13) {
                if      (strncmp(p+1,"C1C",3)) nav->glo_cpbias[0]=str2num(p,5,8);
                else if (strncmp(p+1,"C1P",3)) nav->glo_cpbias[1]=str2num(p,5,8);
                else if (strncmp(p+1,"C2C",3)) nav->glo_cpbias[2]=str2num(p,5,8);
                else if (strncmp(p+1,"C2P",3)) nav->glo_cpbias[3]=str2num(p,5,8);
            }
        }
    }
    else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
        if (nav) nav->leaps=(int)str2num(buff,0,6);
    }
    else if (strstr(label,"# OF SALTELLITES"    )) { /* opt */
        /* skip */ ;
    }
    else if (strstr(label,"PRN / # OF OBS"      )) { /* opt */
        /* skip */ ;
    }
}
/* sync a newline------------------------------------------------------------*/
static int syncnewline(unsigned char *buff, int nb)
{
    if (buff[nb-1]=='\n'||(buff[nb-2]=='\r'&&buff[nb-1]=='\n')) return 1;
    return 0;
}
/* decode obs epoch ----------------------------------------------------------*/
static int decode_obsepoch(char *buff,gtime_t *time,int *flag)
{
    int n;

    if ((n=(int)str2num(buff,32,3))<=0) return 0;
    *flag=(int)str2num(buff,31,1);

    if (3<=*flag&&*flag<=5) return n;
    if (buff[0]!='>'||str2time(buff,1,28,time)) {
        return 0;
    }
    return n;
}
/* decode obs data -----------------------------------------------------------*/
static int decode_obsdata(char *buff, int mask,sigind_t *index, obsd_t *obs)
{
    sigind_t *ind;
    double val[MAXOBSTYPE]={0};
    unsigned char lli[MAXOBSTYPE]={0};
    char satid[8]="";
    int i,j,n,m,stat=1,p[MAXOBSTYPE],k[16],l[16];

    strncpy(satid,buff,3);
    obs->sat=(unsigned char)satid2no(satid);

    if (!obs->sat) {
        trace(4,"decode_obsdata: unsupported sat sat=%s\n",satid);
        stat=0;
    }
    else if (!(satsys(obs->sat,NULL)&mask)) {
        stat=0;
    }
    /* read obs data fields */
    switch (satsys(obs->sat,NULL)) {
        case SYS_GLO: ind=index+1; break;
        case SYS_GAL: ind=index+2; break;
        case SYS_QZS: ind=index+3; break;
        case SYS_SBS: ind=index+4; break;
        case SYS_CMP: ind=index+5; break;
        default:      ind=index  ; break;
    }
    for (i=0,j=3;i<ind->n;i++,j+=16) {
        if (stat) {
            val[i]=str2num(buff,j,14)+ind->shift[i];
            lli[i]=(unsigned char)str2num(buff,j+14,1)&3;
        }
    }
    if (!stat) return 0;

    for (i=0;i<NFREQ+NEXOBS;i++) {
        obs->P[i]=obs->L[i]=0.0; obs->D[i]=0.0f;
        obs->SNR[i]=obs->LLI[i]=obs->code[i]=0;
    }
    /* assign position in obs data */
    for (i=n=m=0;i<ind->n;i++) {

        p[i]=ind->pos[i];

        if (ind->type[i]==0&&p[i]==0) k[n++]=i; /* C1? index */
        if (ind->type[i]==0&&p[i]==1) l[m++]=i; /* C2? index */
    }
    /* save obs data */
    for (i=0;i<ind->n;i++) {
        if (p[i]<0||val[i]==0.0) continue;
        switch (ind->type[i]) {
            case 0: obs->P[p[i]]=val[i]; obs->code[p[i]]=ind->code[i]; break;
            case 1: obs->L[p[i]]=val[i]; obs->LLI [p[i]]=lli[i];       break;
            case 2: obs->D[p[i]]=(float)val[i];                        break;
            case 3: obs->SNR[p[i]]=(unsigned char)(val[i]*4.0+0.5);    break;
        }
    }
    return 1;
}
/* set signal index ----------------------------------------------------------*/
static void sigindex(int sys, const char *opt,char tobs[MAXOBSTYPE][4],
                     sigind_t *ind)
{
    const char *p;
    char str[8],*optstr="";
    double shift;
    int i,j,k,n;

    for (i=n=0;*tobs[i];i++,n++) {
        ind->code[i]=obs2code(tobs[i]+1,ind->frq+i);
        ind->type[i]=(p=strchr(obscodes,tobs[i][0]))?(int)(p-obscodes):0;
        ind->pri[i]=getcodepri(sys,ind->code[i],opt);
        ind->pos[i]=-1;

        /* frequency index for beidou */
        if (sys==SYS_CMP) {
            if      (ind->frq[i]==5) ind->frq[i]=2; /* B2 */
            else if (ind->frq[i]==4) ind->frq[i]=3; /* B3 */
        }
    }
    /* parse phase shift options */
    switch (sys) {
        case SYS_GPS: optstr="-GL%2s=%lf"; break;
        case SYS_GLO: optstr="-RL%2s=%lf"; break;
        case SYS_GAL: optstr="-EL%2s=%lf"; break;
        case SYS_QZS: optstr="-JL%2s=%lf"; break;
        case SYS_SBS: optstr="-SL%2s=%lf"; break;
        case SYS_CMP: optstr="-CL%2s=%lf"; break;
        case SYS_IRN: optstr="-IL%2s=%lf"; break;
    }
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str,&shift)<2) continue;
        for (i=0;i<n;i++) {
            if (strcmp(code2obs(ind->code[i],NULL),str)) continue;
            ind->shift[i]=shift;
            trace(2,"phase shift: sys=%2d tobs=%s shift=%.3f\n",sys,
                  tobs[i],shift);
        }
    }
    /* assign index for highest priority code */
    for (i=0;i<NFREQ;i++) {
        for (j=0,k=-1;j<n;j++) {
            if (ind->frq[j]==i+1&&ind->pri[j]&&(k<0||ind->pri[j]>ind->pri[k])) {
                k=j;
            }
        }
        if (k<0) continue;

        for (j=0;j<n;j++) {
            if (ind->code[j]==ind->code[k]) ind->pos[j]=i;
        }
    }
    /* assign index of extended obs data */
    for (i=0;i<NEXOBS;i++) {
        for (j=0;j<n;j++) {
            if (ind->code[j]&&ind->pri[j]&&ind->pos[j]<0) break;
        }
        if (j>=n) break;

        for (k=0;k<n;k++) {
            if (ind->code[k]==ind->code[j]) ind->pos[k]=NFREQ+i;
        }
    }
    for (i=0;i<n;i++) {
        if (!ind->code[i]||!ind->pri[i]||ind->pos[i]>=0) continue;
        trace(4,"reject obs type: sys=%2d, obs=%s\n",sys,tobs[i]);
    }
    ind->n=n;
}
/* clear message buffer------------------------------------------------------*/
static void clearbuff(raw_t *raw)
{
    int i; for (i=0;i<raw->nbyte;i++) raw->buff[i]=0; raw->nbyte=0;
}
/* reset rinex deocde temp variable------------------------------------------*/
static void resetrnx(raw_t *raw)
{
    raw->rinex.start=0;
    raw->rinex.i=raw->rinex.n=raw->rinex.nsat=0;
}
/* set signal index----------------------------------------------------------*/
static void setsigindex(raw_t *raw)
{
    sigindex(SYS_GPS,"",raw->rinex.tobs[0],raw->rinex.index  );
    sigindex(SYS_GLO,"",raw->rinex.tobs[1],raw->rinex.index+1);
    sigindex(SYS_GAL,"",raw->rinex.tobs[2],raw->rinex.index+2);
    sigindex(SYS_QZS,"",raw->rinex.tobs[3],raw->rinex.index+3);
    sigindex(SYS_SBS,"",raw->rinex.tobs[4],raw->rinex.index+4);
    sigindex(SYS_CMP,"",raw->rinex.tobs[5],raw->rinex.index+5);
    sigindex(SYS_IRN,"",raw->rinex.tobs[6],raw->rinex.index+6);
}
/* input rinex data from stream ---------------------------------------------
 * args   : raw_t *raw         IO  receiver raw data control struct
 *          unsigned char data I   stream data (1 byte)
 * return : status (-1: error message, 0: no message, 1: input observation data,
 *                   2: input ephemeris, 3: input sbas message,
 *                   4: imu measurement data
 *                   6: pvt solution data
 *                   9: input ion/utc parameter)
 *                  32: update signal index
 * notes  : now it only support rinex version 3.02
 *----------------------------------------------------------------------------*/
extern int input_rinex(raw_t *raw, unsigned char data)
{
    int j,tsys=TSYS_GPS,flag=0,mask=SYS_ALL;
    char *p=(char*)raw->buff;

    trace(5,"input_rinex: data=%02x\n",data);

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        resetrnx(raw);
        raw->nbyte=0; return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    if (strstr(p,"END OF HEADER")) { /* end header */

        clearbuff(raw);
        setsigindex(raw); raw->rinex.endhead=1;
        return 32;
    }
    if (!raw->rinex.endhead) {
        /* decode rinex obs header */
        decode_obsh(p,&tsys,raw->rinex.tobs,&raw->nav,NULL);
    }
    else {
        /* decode rinex obs data */
        if (!raw->rinex.start
            &&(raw->rinex.nsat=decode_obsepoch(p,&raw->rinex.time,&flag))) {
            raw->rinex.start=1;
            clearbuff(raw);
            return 0;
        }
        else if (flag<=2||flag==6) {

            raw->rinex.obs[raw->rinex.n].time=raw->rinex.time;

            /* decode obs data */
            if (decode_obsdata(p,mask,raw->rinex.index,raw->rinex.obs+raw->rinex.n)
                &&raw->rinex.n<MAXOBS) raw->rinex.n++;
        }
        /* new site or header info follows */
        else if (flag==3||flag==4) {

            /* decode obs header */
            decode_obsh(p,&tsys,raw->rinex.tobs,NULL,NULL);
        }
        if (++raw->rinex.i>=raw->rinex.nsat) {
            
            for (j=0;j<raw->rinex.n;j++) {

                /* utc -> gpst */
                if (tsys==TSYS_UTC) {
                    raw->rinex.obs[j].time=utc2gpst(raw->rinex.obs[j].time);
                }
                raw->obs.data[j]=raw->rinex.obs[j];
            }
            raw->time=raw->rinex.time;
            raw->obs.n=raw->rinex.n;
            resetrnx(raw);
            clearbuff(raw); return 1;
        }
    }
    clearbuff(raw);
    return 0;
}
