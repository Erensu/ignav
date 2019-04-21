/*------------------------------------------------------------------------------
* rtksvr.c : rtk server functions
*
* options : -DWIN32    use WIN32 API
*
* version : $Revision:$ $Date:$
* history : 2009/01/07  1.0  new
*           2009/06/02  1.1  support glonass
*           2010/07/25  1.2  support correction input/log stream
*                            supoort online change of output/log streams
*                            supoort monitor stream
*                            added api:
*                                rtksvropenstr(),rtksvrclosestr()
*                            changed api:
*                                rtksvrstart()
*           2010/08/25  1.3  fix problem of ephemeris time inversion (2.4.0_p6)
*           2010/09/08  1.4  fix problem of ephemeris and ssr squence upset
*                            (2.4.0_p8)
*           2011/01/10  1.5  change api: rtksvrstart(),rtksvrostat()
*           2011/06/21  1.6  fix ephemeris handover problem
*           2012/05/14  1.7  fix bugs
*           2013/03/28  1.8  fix problem on lack of glonass freq number in raw
*                            fix problem on ephemeris with inverted toe
*                            add api rtksvrfree()
*           2014/06/28  1.9  fix probram on ephemeris update of beidou
*           2015/04/29  1.10 fix probram on ssr orbit/clock inconsistency
*           2015/07/31  1.11 add phase bias (fcb) correction
*           2015/12/05  1.12 support opt->pppopt=-DIS_FCB
*           2016/07/01  1.13 support averaging single pos as base position
*           2016/07/31  1.14 fix bug on ion/utc parameters input
*           2016/08/20  1.15 support api change of sendnmea()
*           2016/09/18  1.16 fix server-crash with server-cycle > 1000
*           2016/09/20  1.17 change api rtksvrstart()
*           2016/10/01  1.18 change api rtksvrstart()
*           2016/10/04  1.19 fix problem to send nmea of single solution
*           2016/10/09  1.20 add reset-and-single-sol mode for nmea-request
*           2017/04/11  1.21 add rtkfree() in rtksvrfree()
*----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants ----------------------------------------------------------------*/
#define MIN_INT_RESET   30000      /* mininum interval of reset command (ms) */
#define DTTOLM          0.1        /* threshold of rover and base observation data */
#define REALTIME        0          /* real time process rover observation data */
#define MAXTIMEDIFF     0.5        /* max time difference for suspend input stream */
#define OUTSOLFRQ       50         /* frequency of output ins solutions */

#define NS(i,j,max)     ((((j)-1)%(max)-(i))<0?(((j)-1)%(max)-(i)+(max)):(((j)-1)%(max)-(i)))
#define NE(i,j,max)     MAX(0,(((i)-(j))<0?((i)-(j)+(max)):((i)-(j))))

/* constants/global variables ------------------------------------------------*/
struct obs {int n; obsd_t data[MAXOBS];};
struct imu {int n; imud_t data[MAXIMU];};
struct pvt {int n; sol_t  data[MAXSOL];};
struct img {int n; img_t  data[MAXIMG];};

/* write solution header to output stream ------------------------------------*/
static void writesolhead(stream_t *stream, const solopt_t *solopt)
{
    unsigned char buff[1024];
    int n;
    
    n=outsolheads(buff,solopt);
    strwrite(stream,buff,n);
}
/* save output buffer --------------------------------------------------------*/
static void saveoutbuf(rtksvr_t *svr, unsigned char *buff, int n, int index)
{
    rtksvrlock(svr);
    
    n=n<svr->buffsize-svr->nsb[index]?n:svr->buffsize-svr->nsb[index];
    memcpy(svr->sbuf[index]+svr->nsb[index],buff,n);
    svr->nsb[index]+=n;
    
    rtksvrunlock(svr);
}
/* write solution to output stream -------------------------------------------*/
static void writesol(rtksvr_t *svr, int index)
{
    prcopt_t *opt=&svr->rtk.opt;
    unsigned char buff[MAXSOLMSG+1];
    static int c=0; int i,n;

    tracet(4,"writesol: index=%d\n",index);

    /* update ins update solution status */
    if (opt->mode>=PMODE_INS_UPDATE&&opt->mode<=PMODE_INS_TGNSS) {
        ins2sol(&svr->rtk.ins,&opt->insopt,&svr->rtk.sol);
    }
    /* camera pose relative to initial position */
    if (opt->mode==PMODE_VO&&svr->rtk.sol.stat==SOLQ_VO) {
        for (i=0;i<3;i++) {
            svr->rtk.sol.rr[i+0]=svr->rtk.ins.vo.rc[i];
            svr->rtk.sol.rr[i+3]=svr->rtk.ins.vo.vc[i];
        }
        svr->rtk.sol.time=svr->rtk.ins.vo.time;
        svr->rtk.sol.stat=SOLQ_VO;
    }
    /* write solution status output stream */
    for (i=0;i<2;i++) {
        if (svr->solopt[i].posf==SOLF_STAT) { /* solution status */

            /* output solution status */
            rtksvrlock(svr);
            n=rtkoutstat(&svr->rtk,(char*)buff);
            rtksvrunlock(svr);
        }
        else {
            /* output solution */
            n=outsols(buff,&svr->rtk.sol,svr->rtk.rb,svr->solopt+i,
                      opt->mode>=PMODE_INS_UPDATE?&svr->rtk.ins:NULL,
                      &svr->rtk.opt.insopt,0);
        }
        trace(3,"solution: %s",buff);

        strwrite(svr->stream+i+7,buff,n);

        /* save output buffer */
        saveoutbuf(svr,buff,n,i);

        /* output extended solution */
        n=outsolexs(buff,&svr->rtk.sol,svr->rtk.ssat,svr->solopt+i,0);
        strwrite(svr->stream+i+7,buff,n);

        /* save output buffer */
        saveoutbuf(svr,buff,n,i);
    }
    /* output solution to monitor port */
    const solopt_t *psolopt=NULL;
    if (opt->mode>=PMODE_INS_UPDATE&&opt->mode<=PMODE_INS_TGNSS) {
        psolopt=&solopt_ins_default;
    }
    else if (opt->mode==PMODE_VO) {
        psolopt=&solopt_vo_default;
    }
    else {
        psolopt=&solopt_default;
    }
    if (svr->moni) {
        n=outsols(buff,&svr->rtk.sol,svr->rtk.rb,psolopt,
                  opt->mode>=PMODE_INS_UPDATE?&svr->rtk.ins:NULL,
                  &svr->rtk.opt.insopt,1);
        if (opt->mode>=PMODE_INS_UPDATE&&opt->mode<=PMODE_INS_TGNSS) {
            if (c++>OUTSOLFRQ) {
                strwrite(svr->moni,buff,n); c=0;
            }
        }
        else {                                                      
            strwrite(svr->moni,buff,n);
        }
    }
    /* output ground truth solution to monitor port */
    if (svr->groundtruth) {

        /* search ground truth solution */
        i=findgtsols(&svr->gtsols,svr->rtk.sol.time);
        if (i>=0) {
            n=outgroundtruth(buff,&svr->gtsols.data[i],&solopt_default);
        }
        strwrite(svr->groundtruth,buff,n);
        trace(3,"solution: %s",buff);
    }
    /* save solution buffer */
    if (svr->nsol<MAXSOLBUF) {
        rtksvrlock(svr);
        svr->solbuf[svr->nsol++]=svr->rtk.sol;
        rtksvrunlock(svr);
    }
}
/* update navigation data ----------------------------------------------------*/
static void updatenav(nav_t *nav)
{
    int i,j;
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        nav->lam[i][j]=satwavelen(i+1,j,nav);
    }
}
/* update glonass frequency channel number in raw data struct ----------------*/
static void updatefcn(rtksvr_t *svr)
{
    int i,j,sat,frq;

    trace(3,"updatefcn:\n");
    
    for (i=0;i<MAXPRNGLO;i++) {
        sat=satno(SYS_GLO,i+1);
        
        for (j=0,frq=-999;j<3;j++) {
            if (svr->raw[j].nav.geph[i].sat!=sat) continue;
            frq=svr->raw[j].nav.geph[i].frq;
        }
        if (frq<-7||frq>6) continue;
        
        for (j=0;j<3;j++) {
            if (svr->raw[j].nav.geph[i].sat==sat) continue;
            svr->raw[j].nav.geph[i].sat=sat;
            svr->raw[j].nav.geph[i].frq=frq;
        }
    }
}
/* update observation data----------------------------------------------------*/
static void updateobs(rtksvr_t *svr,obs_t *obs,int index,int iobs)
{
    trace(3,"updateobs (index=%d,iobs=%d):\n");

    int i,n=0; if (iobs<MAXOBSBUF) {
        for (i=0;i<obs->n;i++) {
            if (svr->rtk.opt.exsats[obs->data[i].sat-1]==1||
                !(satsys(obs->data[i].sat,NULL)&svr->rtk.opt.navsys)) continue;
            svr->obs[index][iobs].data[n]=obs->data[i];
            svr->obs[index][iobs].data[n++].rcv=index+1;

            time2gpst(obs->data[i].time,&svr->week);
        }
        svr->obs[index][iobs].n=n;
        sortobs(&svr->obs[index][iobs]);
    }
    svr->nmsg[index][0]++;
}
/* update ephemeris-----------------------------------------------------------*/
static void updateeph(rtksvr_t *svr,nav_t *nav,int sat,int index)
{
    int prn;
    eph_t *eph1,*eph2,*eph3;
    geph_t *geph1,*geph2,*geph3;

    trace(3,"updateeph:\n");

    if (satsys(sat,&prn)!=SYS_GLO) {
        if (!svr->navsel||svr->navsel==index+1) {
            eph1=nav->eph+sat-1;
            eph2=svr->nav.eph+sat-1;
            eph3=svr->nav.eph+sat-1+MAXSAT;
            if (eph2->ttr.time==0||
                (eph1->iode!=eph3->iode&&eph1->iode!=eph2->iode)||
                (timediff(eph1->toe,eph3->toe)!=0.0&&
                 timediff(eph1->toe,eph2->toe)!=0.0)) {
                *eph3=*eph2;
                *eph2=*eph1;
                updatenav(&svr->nav);
            }
        }
        svr->nmsg[index][1]++;
    }
    else {
        if (!svr->navsel||svr->navsel==index+1) {
            geph1=nav->geph+prn-1;
            geph2=svr->nav.geph+prn-1;
            geph3=svr->nav.geph+prn-1+MAXPRNGLO;
            if (geph2->tof.time==0||
                (geph1->iode!=geph3->iode&&geph1->iode!=geph2->iode)) {
                *geph3=*geph2;
                *geph2=*geph1;
                updatenav(&svr->nav);
                updatefcn(svr);
            }
        }
        svr->nmsg[index][6]++;
    }
}
/* update sbas message data---------------------------------------------------*/
static void updatesbas(rtksvr_t *svr, sbsmsg_t *sbsmsg, int index)
{
    int sbssat=svr->rtk.opt.sbassatsel,i;

    trace(3,"updatesbas:\n");

    if (sbsmsg&&(sbssat==sbsmsg->prn||sbssat==0)) {
        if (svr->nsbs<MAXSBSMSG) {
            svr->sbsmsg[svr->nsbs++]=*sbsmsg;
        }
        else {
            for (i=0;i<MAXSBSMSG-1;i++) svr->sbsmsg[i]=svr->sbsmsg[i+1];
            svr->sbsmsg[i]=*sbsmsg;
        }
        sbsupdatecorr(sbsmsg,&svr->nav);
    }
    svr->nmsg[index][3]++;
}
/* update ion/utc parameters--------------------------------------------------*/
static void updateionutc(rtksvr_t *svr, nav_t *nav,int index)
{
    trace(3,"updateionutc:\n");

    int i; if (svr->navsel==0||svr->navsel==index+1) {
        for (i=0;i<8;i++) svr->nav.ion_gps[i]=nav->ion_gps[i];
        for (i=0;i<4;i++) svr->nav.utc_gps[i]=nav->utc_gps[i];
        for (i=0;i<4;i++) svr->nav.ion_gal[i]=nav->ion_gal[i];
        for (i=0;i<4;i++) svr->nav.utc_gal[i]=nav->utc_gal[i];
        for (i=0;i<8;i++) svr->nav.ion_qzs[i]=nav->ion_qzs[i];
        for (i=0;i<4;i++) svr->nav.utc_qzs[i]=nav->utc_qzs[i];
        svr->nav.leaps=nav->leaps;
    }
    svr->nmsg[index][2]++;
}
/* update antenna postion parameters------------------------------------------*/
static void updateant(rtksvr_t *svr,int index)
{
    int i;
    double pos[3],del[3]={0},dr[3];

    trace(3,"updateant:\n");

    if (svr->rtk.opt.refpos==POSOPT_RTCM&&index==1) {
        for (i=0;i<3;i++) {
            svr->rtk.rb[i]=svr->rtcm[1].sta.pos[i];
        }
        /* antenna delta */
        ecef2pos(svr->rtk.rb,pos);
        if (svr->rtcm[1].sta.deltype) { /* xyz */
            del[2]=svr->rtcm[1].sta.hgt;
            enu2ecef(pos,del,dr);
            for (i=0;i<3;i++) {
                svr->rtk.rb[i]+=svr->rtcm[1].sta.del[i]+dr[i];
            }
        }
        else { /* enu */
            enu2ecef(pos,svr->rtcm[1].sta.del,dr);
            for (i=0;i<3;i++) {
                svr->rtk.rb[i]+=dr[i];
            }
        }
    }
    else if (svr->rtk.opt.refpos==POSOPT_RAW&&index==1) {
        for (i=0;i<3;i++) {
            svr->rtk.rb[i]=svr->raw[1].sta.pos[i];
        }
        /* antenna delta */
        ecef2pos(svr->rtk.rb,pos);
        if (svr->raw[1].sta.deltype) { /* xyz */
            del[2]=svr->raw[1].sta.hgt;
            enu2ecef(pos,del,dr);
            for (i=0;i<3;i++) {
                svr->rtk.rb[i]+=svr->raw[1].sta.del[i]+dr[i];
            }
        }
        else { /* enu */
            enu2ecef(pos,svr->raw[1].sta.del,dr);
            for (i=0;i<3;i++) {
                svr->rtk.rb[i]+=dr[i];
            }
        }
    }
    svr->nmsg[index][4]++;
}
/* update dgps correction data------------------------------------------------*/
static void updatedgps(rtksvr_t *svr,int index)
{
    trace(3,"updatedgps:\n"); svr->nmsg[index][5]++;
}
/* update ssr message---------------------------------------------------------*/
static void updatessr(rtksvr_t *svr,int index)
{
    trace(3,"updatessr:\n");

    int i,prn,sys,iode;
    for (i=0;i<MAXSAT;i++) {
        if (!svr->rtcm[index].ssr[i].update) continue;

        /* check consistency between iods of orbit and clock */
        if (svr->rtcm[index].ssr[i].iod[0]!=
            svr->rtcm[index].ssr[i].iod[1]) continue;

        svr->rtcm[index].ssr[i].update=0;

        iode=svr->rtcm[index].ssr[i].iode;
        sys=satsys(i+1,&prn);

        /* check corresponding ephemeris exists */
        if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS) {
            if (svr->nav.eph[i       ].iode!=iode&&
                svr->nav.eph[i+MAXSAT].iode!=iode) {
                continue;
            }
        }
        else if (sys==SYS_GLO) {
            if (svr->nav.geph[prn-1          ].iode!=iode&&
                svr->nav.geph[prn-1+MAXPRNGLO].iode!=iode) {
                continue;
            }
        }
        svr->nav.ssr[i]=svr->rtcm[index].ssr[i];
    }
    svr->nmsg[index][7]++;
}
/* update lex message---------------------------------------------------------*/
static void updatelex(rtksvr_t *svr,int index)
{
    trace(3,"updatelex:\n");
    gtime_t tof;
    lexupdatecorr(&svr->raw[index].lexmsg,&svr->nav,&tof);
    svr->nmsg[index][8]++;
}
/* update error message-------------------------------------------------------*/
static void updateerror(rtksvr_t *svr,int index)
{
    trace(3,"updateerror:\n"); svr->nmsg[index][9]++;
}
/* update imu measurement data------------------------------------------------*/
static void updateimu(rtksvr_t *svr, const imud_t *imu, int iimu)
{
    int week;

    trace(3,"updateimu:\n");

    svr->imu[iimu]=*imu;
    if (time2gpst(imu->time,&week)>0.0&&week==0) {
        svr->imu[iimu].time=timeadd(imu->time,svr->week*604800.0);
    }
}
/* update PVT solutions data--------------------------------------------------*/
static void updatepvt(rtksvr_t *svr, const sol_t *sol ,int isol)
{
    trace(3,"updatepvt:\n");

    time2gpst(sol->time,&svr->week);
    
    if (isol>=MAXSOLBUF) return;
    if (fabs(timediff(svr->pvt[isol-1].time,sol->time))<DTTOL
        &&isol>=1) return;

    svr->pvt[isol]=*sol;
}
/* update image raw data------------------------------------------------------*/
static void updateimg(rtksvr_t *svr, const img_t *img ,int iimg)
{
    trace(3,"updateimg:\n");

    if (iimg>=MAXIMGBUF) return;
    if (fabs(timediff(svr->img[iimg-1].time,img->time))<DTTOL
        &&iimg>=1) return;

    /* update image raw data buffer */
    if (svr->img[iimg].data) {
        if (svr->img[iimg].flag) {
            free(svr->img[iimg].data); svr->img[iimg].data=NULL;
        }
        else {
            free(svr->img[iimg].data);
            svr->img[iimg].data=NULL;
        }
    }
    if (svr->img[iimg].datar) {
        if (svr->img[iimg].flag) {
            free(svr->img[iimg].datar);
            svr->img[iimg].datar=NULL;
        }
        else {
            free(svr->img[iimg].datar);
            svr->img[iimg].datar=NULL;
        }
    }
    svr->img[iimg]=*img;
}
/* update pose measurement data-----------------------------------------------*/
static void updatepose(rtksvr_t *svr,const pose_meas_t *pose,int ipose)
{
    trace(3,"updatepose:\n");
    if (ipose>=MAXPOSEBUF) return;
    if (ipose>=1&&fabs(timediff(svr->pose[ipose-1].time,svr->pose[ipose].time))
                  <DTTOL) {
        return;
    }
    svr->pose[ipose]=*pose;
}
/* solution convert to gnss measurement data----------------------------------*/
static void sol2gnss(const sol_t *sol,gmea_t *gnss)
{
    /* number of valid satellite and solution status */
    gnss->ns=sol->ns; gnss->stat=sol->stat;
    gnss->t=sol->time;

    /* position and velocity */
    matcpy(gnss->pe,sol->rr+0,1,3);
    matcpy(gnss->ve,sol->rr+3,1,3);

    /* covariance of position and velocity */
    int i; for (i=0;i<3;i++) {
        gnss->std[i+0]=SQRT(sol->qr[i]);
        gnss->std[i+3]=SQRT(sol->qv[i]);
    }
}
/* input rover and base observation data--------------------------------------*/
static int inputobs(rtksvr_t *svr,obsd_t *obs)
{
    int i,j=0,k,nr,nb=0,n=0;
    const obs_t *robs=NULL,*bobs=NULL;

    tracet(3,"inputobs:\n");

    /* collect rover observation data */
    for (nr=0,robs=svr->obs[0],i=0;i<robs[svr->syn.rover].n;i++) {
        obs[nr++]=robs[svr->syn.rover].data[i];
    }
    /* check position mode if need base observation */
    switch (svr->rtk.opt.mode) {
        case PMODE_SINGLE    : return nr;
        case PMODE_PPP_KINEMA: return nr;
        case PMODE_PPP_STATIC: return nr;
        case PMODE_PPP_FIXED : return nr;
    }
    /* check observation time sync */
    if (!svr->syn.tali[0]&&!svr->syn.tali[2]) {
        trace(2,"check time alignment fail\n");
        return 0;
    }
    /* search base observation data by rover obs time */
    for (i=0,bobs=svr->obs[1],n=NS(svr->syn.base,svr->syn.nb,MAXOBSBUF);
         i<n+1&&nr;i++) {

        j=svr->syn.base+i;
        if (j>=MAXOBSBUF) j%=MAXOBSBUF;

        if (fabs(timediff(bobs[j].data[0].time,obs[0].time))<DTTOLM) {
            for (nb=0,k=0;k<bobs[i].n;k++) {
                obs[nr+nb++]=bobs[j].data[k];
            }
            /* update base sync index */
            if (nb) svr->syn.base=j; break;
        }
    }
#if REALTIME
    if ((svr->syn.rover+=(nr>0))>=MAXOBSBUF) {
        svr->syn.rover%=MAXOBSBUF;
        trace(3,"rover sync index overflow\n");
    }    
#else
    /* no update rover sync index if no base observation data */
    if ((svr->syn.rover+=(nr>0&&nb>0))>=MAXOBSBUF) {
        svr->syn.rover%=MAXOBSBUF;
        trace(3,"rover sync index overflow\n");
    }
    /* check number of base observation data */
    if (nb<=0) {
        trace(2,"no match base observation data\n");
        return 0;
    }
#endif
    return nr+nb;
}
/* input rover and base observation data for ins tightly coupled-------------*/
static int inputobstc(rtksvr_t *svr,gtime_t time,obsd_t *obs)
{
    tracet(3,"inputobstc: time=%s\n",time_str(time,4));

    int i,j,k,n=0,m=0,flag=0;
    double dt;
    gtime_t tr={0};

    const obs_t *robs=svr->obs[0],*bobs=svr->obs[1];

    /* match rover observation data */
    for (i=0,dt=0.0,n=NS(svr->syn.rover,svr->syn.nr,MAXOBSBUF);
         i<n+1&&svr->syn.nr;i++) {
        j=svr->syn.rover+i;

        if (j>=MAXOBSBUF) j=j%MAXOBSBUF;
        if (dt&&fabs(dt)<fabs(timediff(time,robs[j].data[0].time))) {
            break;
        }
        if (fabs((dt=timediff(time,robs[j].data[0].time)))<DTTOL
            &&robs[j].n!=0) {
            for (k=0,m=0;k<robs[j].n;k++) obs[m++]=robs[j].data[k];
            svr->syn.rover=j; tr=obs[0].time;
            flag=1; break;
        }
    }
    /* match base observation data by rover time */
    if (tr.time) {
        for (i=0,dt=0.0,n=NS(svr->syn.base,svr->syn.nb,MAXOBSBUF);
             i<n+1&&svr->syn.nb;i++) {
            j=svr->syn.base+i;

            if (j>=MAXOBSBUF) j=j%MAXOBSBUF;
            if (dt&&fabs(dt)<fabs(timediff(tr,bobs[j].data[0].time))) {
                break;
            }
            if (fabs((dt=timediff(tr,bobs[j].data[0].time)))<DTTOLM
                &&bobs[j].n!=0) {
                for (k=0;k<bobs[j].n;k++) obs[m++]=bobs[j].data[k];
                svr->syn.base=j; flag=2; break;
            }
        }
    }
#if REALTIME
    /* always input rover observation data despite base */
    return m;
#else
    /* ok when rover and base have input */
    if (flag!=2) return 0; else return m;
#endif
}
/* input imu measurement data------------------------------------------------*/
static int inputimu(rtksvr_t *svr,imud_t *data)
{
    int i,j,k,n=0;

    tracet(3,"inputimu:\n");

    /* check time alignment of input streams */
    if (!svr->syn.tali[1]&&!svr->syn.tali[2]&&!svr->syn.tali[3]) {
        trace(2,"check time alignment fail\n");
        return 0;
    }
    for (i=0,k=0,n=NS(svr->syn.imu,svr->syn.ni,MAXIMUBUF),
        j=svr->syn.imu;i<n+1;i++) {
        j=svr->syn.imu+i;
        if (j>=MAXIMUBUF) j=j%MAXIMUBUF;
        if (k<MAXIMU) data[k++]=svr->imu[j]; else break;
    }
    svr->syn.imu=(j+1)%MAXIMUBUF;
    return k;
}
/* input pvt solution data----------------------------------------------------*/
static int inputpvt(rtksvr_t *svr,gtime_t t0,sol_t *sol)
{
    tracet(3,"inputpvt:\n");

    sol_t sol0={0}; *sol=sol0;
    int i,j,n=0; double dt;

    /* search pvt solution by t0 */
    for (i=0,dt=0.0,n=NS(svr->syn.pvt,svr->syn.ns,MAXSOLBUF);
        i<n+1&&svr->syn.ns;i++) {
        j=svr->syn.pvt+i;

        if (j>=MAXSOLBUF) j=j%MAXSOLBUF;
        if (dt&&fabs(dt)<fabs(timediff(t0,svr->pvt[j].time))) {
            break;
        }
        if (fabs((dt=timediff(t0,svr->pvt[j].time)))<DTTOL
            &&svr->pvt[j].time.time!=0
            &&svr->pvt[j].stat!=SOLQ_NONE) {
            memcpy(sol,&svr->pvt[j],sizeof(sol_t));
            svr->syn.pvt=j;
            return 1; /* ok */
        }
    }
    return 0; /* fail */
}
/* input image raw data-------------------------------------------------------*/
static int inputimg(rtksvr_t *svr,gtime_t t0,img_t **img)
{
    trace(3,"inputimg:\n");
    int i,j,n; double dt;

    /* search image raw data by t0 */
    for (i=0,dt=0.0,n=NS(svr->syn.img,svr->syn.nm,MAXIMGBUF);
         i<n+1&&svr->syn.nm;i++) {
        j=svr->syn.img+i;

        if (j>=MAXIMGBUF) j=j%MAXIMGBUF;
        if (dt&&fabs(dt)<fabs(timediff(t0,svr->img[j].time))) {
            break;
        }
        if (fabs((dt=timediff(t0,svr->img[j].time)))<DTTOL
            &&svr->img[j].flag==0&&svr->img[j].h&&svr->img[j].w
            &&svr->img[j].time.time!=0) {
            *img=&svr->img[j];
            svr->syn.img=j; return 1; /* ok */
        }
    }
    *img=NULL;
    return 0; /* fail */
}
/* input image measurement data-----------------------------------------------*/
static int inputimage(rtksvr_t *svr,img_t *data)
{
    int i,j,k,n=0;

    for (i=0,k=0,n=NS(svr->syn.img,svr->syn.nm,MAXIMGBUF),
                 j=svr->syn.img;i<n+1;i++) {
        j=svr->syn.img+i;
        if (j>=MAXIMGBUF) j=j%MAXIMGBUF;
        if (k<MAXIMG&&svr->img[j].data) {
            data[k++]=svr->img[j]; svr->img[j].flag=1;
        }
        else break;
    }
    svr->syn.img=j;
    return k;
}
/* input pose measurement data------------------------------------------------*/
static int inputpose(rtksvr_t *svr,gtime_t time,pose_meas_t *pose)
{
    double dt=0.0;
    int i,j,np;

    np=NS(svr->syn.pose,svr->syn.np,MAXPOSEBUF);

    for (i=0;i<np+1;i++) {
        j=svr->syn.pose+i;
        if (j>=MAXPOSEBUF) j=j%MAXPOSEBUF;
        if (dt&&fabs(dt)<fabs(timediff(time,svr->pose[j].time))) {
            break;
        }
        if (fabs((dt=timediff(time,svr->pose[j].time)))<DTTOL
            &&svr->pose[j].time.time
            &&svr->pose[j].stat
            &&norm(svr->pose[j].rpy,3)>0.0) {
            pose[0]=svr->pose[j];
            svr->syn.pose=j; return 1;
        }
    }
    return 0;
}
/* decode receiver raw/rtcm data ---------------------------------------------*/
static int decoderaw(rtksvr_t *svr, int index)
{
    obs_t *obs=NULL;
    nav_t *nav=NULL;
    imu_t *imu=NULL;
    sol_t *sol=NULL;
    img_t *img=NULL;
    pose_meas_t *pose=NULL;
    sbsmsg_t *sbsmsg=NULL;
    int i,j,ret=0,sat,k=0,fobs=0;
    
    tracet(4,"decoderaw: index=%d\n",index);
    
    rtksvrlock(svr);
    switch (index) {
        case 0: fobs=svr->syn.nr; svr->syn.nrp=svr->syn.nr; break;
        case 1: fobs=svr->syn.nb; svr->syn.nbp=svr->syn.nb; break;
        case 3: fobs=svr->syn.ns; svr->syn.nsp=svr->syn.ns; break;
        case 4: fobs=svr->syn.ni; svr->syn.nip=svr->syn.ni; break;
        case 5: fobs=svr->syn.nm; svr->syn.nmp=svr->syn.nm; break;
        case 6: fobs=svr->syn.np; svr->syn.npp=svr->syn.np; break;
    }
    for (i=0;i<svr->nb[index];i++) {
        
        /* input rtcm/receiver raw data from stream */
        if (svr->format[index]==STRFMT_RTCM2) {
            ret=input_rtcm2(svr->rtcm+index,svr->buff[index][i]);
            obs=&svr->rtcm[index].obs;
            nav=&svr->rtcm[index].nav;
            sat=svr->rtcm[index].ephsat;
        }
        else if (svr->format[index]==STRFMT_RTCM3) {
            ret=input_rtcm3(svr->rtcm+index,svr->buff[index][i]);
            obs=&svr->rtcm[index].obs;
            nav=&svr->rtcm[index].nav;
            sat=svr->rtcm[index].ephsat;
        }
        else {
            ret=input_raw(svr->raw+index,svr->format[index],svr->buff[index][i]);
            obs=&svr->raw[index].obs;
            imu=&svr->raw[index].imut;
            nav=&svr->raw[index].nav;
            sol=&svr->raw[index].sol;
            img=&svr->raw[index].img;
            sat=svr->raw[index].ephsat;
            sbsmsg=&svr->raw[index].sbsmsg;
            pose=&svr->raw[index].pose;
        }
        /* update cmr rover observations cache */
        if (svr->format[1]==STRFMT_CMR&&index==0&&ret==1) {
            update_cmr(&svr->raw[1],svr,obs);
        }
        /* update observation data */
        if (ret==1&&obs->n) {
            updateobs(svr,obs,index,fobs++%MAXOBSBUF);
            if (k<MAXOBSBUF) k++; else svr->prcout++;
            if (fobs>=MAXOBSBUF&&index==0) svr->syn.of[0]=fobs/MAXOBSBUF;
            if (fobs>=MAXOBSBUF&&index==1) svr->syn.of[1]=fobs/MAXOBSBUF;
        }
        /* update ephemeris */
        if (ret==2) {
            updateeph(svr,nav,sat,index);
        }
        /* update sbas messages data */
        if (ret==3) {
            updatesbas(svr,sbsmsg,index);
        }
        /* update imu measurement data */
        if (ret==4) {
            for (j=0;j<imu->n;j++) {
                adjustimu(&svr->rtk.opt,&imu->data[j]);
                updateimu(svr,&imu->data[j],fobs++%MAXIMUBUF);
                if (k<MAXIMUBUF) k++; else svr->iprcout++;
                if (fobs>=MAXIMUBUF) svr->syn.of[2]=fobs/MAXIMUBUF;
            }
        }
        /* update antenna postion parameters */
        if (ret==5) {
            updateant(svr,index);
        }
        /* update GPS solutions for ins/gnss coupled */
        if (ret==6&&svr->format[index]!=STRFMT_RTCM2&&svr->format[index]!=STRFMT_RTCM3) {
            updatepvt(svr,sol,fobs++%MAXSOLBUF);
            if (k<MAXSOLBUF) k++; else svr->gprcout++;
            if (fobs>=MAXSOLBUF) svr->syn.of[3]=fobs/MAXSOLBUF;
        }
        /* update dgps correction data */
        if (ret==7) {
            updatedgps(svr,index);
        }
        /* update ion/utc parameter */
        if (ret==9) {
            updateionutc(svr,nav,index);
        }
        /* update ssr message data */
        if (ret==10) updatessr(svr,index);

        /* update image raw data */
        if (ret==11) {
            updateimg(svr,img,fobs++%MAXIMGBUF);
            if (fobs>=MAXIMGBUF) svr->syn.of[4]=fobs/MAXIMGBUF;
        }
        /* update lex message data */
        if (ret==31) {
            updatelex(svr,index);
        }
        /* update signal index */
        if (ret==32&&(index==0||index==1)) {
            /* observation signal index for rover and base */
            for (j=0;j<7;j++) svr->rtk.opt.sind[0][j]=svr->raw[0].rinex.index[j];
            for (j=0;j<7;j++) svr->rtk.opt.sind[1][j]=svr->raw[1].rinex.index[j];

            for (j=0;j<7;j++) svr->nav.sind[0][j]=svr->raw[0].rinex.index[j];
            for (j=0;j<7;j++) svr->nav.sind[1][j]=svr->raw[1].rinex.index[j];
        }
        /* update pose measurement */
        if (ret==34) {
            updatepose(svr,pose,fobs++%MAXPOSEBUF);
            if (fobs>=MAXPOSEBUF) {
                svr->syn.of[5]=fobs/MAXPOSEBUF;
            }
        }
        if (ret==-1) updateerror(svr,index);
    }
    /* update time synchronization index */
    switch (index) {
        case 0: svr->syn.nr=fobs; break;
        case 1: svr->syn.nb=fobs; break;
        case 3: svr->syn.ns=fobs; break;
        case 4: svr->syn.ni=fobs; break;
        case 5: svr->syn.nm=fobs; break;
        case 6: svr->syn.np=fobs; break;
    }
    /* update current time of measurement */
    if (fobs) {
        if (index==0) svr->syn.time[0]=svr->obs[0][(fobs-1)%MAXOBSBUF].data[0].time;
        if (index==1) svr->syn.time[1]=svr->obs[1][(fobs-1)%MAXOBSBUF].data[0].time;
        if (index==3) svr->syn.time[3]=svr->pvt[(fobs-1)%MAXSOLBUF].time;
        if (index==4) svr->syn.time[2]=svr->imu[(fobs-1)%MAXIMUBUF].time;
        if (index==5) svr->syn.time[4]=svr->img[(fobs-1)%MAXIMGBUF].time;
        if (index==6) svr->syn.time[5]=svr->pose[(fobs-1)%MAXPOSEBUF].time;
    }
    switch (index) {
        case 0: k=MIN(NE(fobs%MAXOBSBUF,svr->syn.rover,MAXOBSBUF),MAXOBS); break;
        case 1: k=MIN(NE(fobs%MAXOBSBUF,svr->syn.base ,MAXOBSBUF),MAXOBS); break;
        case 3: k=MIN(NE(fobs%MAXSOLBUF,svr->syn.pvt  ,MAXSOLBUF),MAXSOL); break;
        case 4: k=MIN(NE(fobs%MAXIMUBUF,svr->syn.imu  ,MAXIMUBUF),MAXIMU); break;
        case 5: k=MIN(NE(fobs%MAXIMGBUF,svr->syn.img  ,MAXIMGBUF),MAXIMG); break;
        case 6: k=MIN(NE(fobs%MAXPOSEBUF,svr->syn.pose,MAXPOSEBUF),MAXPOSE); break;
    }
    svr->nb[index]=0;
    rtksvrunlock(svr);
    return k;
}
/* decode download file ------------------------------------------------------*/
static void decodefile(rtksvr_t *svr, int index)
{
    nav_t nav={0};
    char file[1024];
    int nb;
    
    tracet(4,"decodefile: index=%d\n",index);
    
    rtksvrlock(svr);
    
    /* check file path completed */
    if ((nb=svr->nb[index])<=2||
        svr->buff[index][nb-2]!='\r'||svr->buff[index][nb-1]!='\n') {
        rtksvrunlock(svr);
        return;
    }
    strncpy(file,(char*)svr->buff[index],nb-2); file[nb-2]='\0';
    svr->nb[index]=0;
    
    rtksvrunlock(svr);
    
    if (svr->format[index]==STRFMT_SP3) { /* precise ephemeris */
        
        /* read sp3 precise ephemeris */
        readsp3(file,&nav,0);
        if (nav.ne<=0) {
            tracet(1,"sp3 file read error: %s\n",file);
            return;
        }
        /* update precise ephemeris */
        rtksvrlock(svr);
        
        if (svr->nav.peph) free(svr->nav.peph);
        svr->nav.ne=svr->nav.nemax=nav.ne;
        svr->nav.peph=nav.peph;
        svr->ftime[index]=utc2gpst(timeget());
        strcpy(svr->files[index],file);
        
        rtksvrunlock(svr);
    }
    else if (svr->format[index]==STRFMT_RNXCLK) { /* precise clock */
        
        /* read rinex clock */
        if (readrnxc(file,&nav)<=0) {
            tracet(1,"rinex clock file read error: %s\n",file);
            return;
        }
        /* update precise clock */
        rtksvrlock(svr);

        if (svr->nav.pclk) free(svr->nav.pclk);
        svr->nav.nc=svr->nav.ncmax=nav.nc;
        svr->nav.pclk=nav.pclk;
        svr->ftime[index]=utc2gpst(timeget());
        strcpy(svr->files[index],file);

        rtksvrunlock(svr);
    }
}
/* periodic command ----------------------------------------------------------*/
static void periodic_cmd(int cycle, const char *cmd, stream_t *stream)
{
    const char *p=cmd,*q;
    char msg[1024],*r;
    int n,period;
    
    for (p=cmd;;p=q+1) {
        for (q=p;;q++) if (*q=='\r'||*q=='\n'||*q=='\0') break;
        n=(int)(q-p); strncpy(msg,p,n); msg[n]='\0';
        
        period=0;
        if ((r=strrchr(msg,'#'))) {
            sscanf(r,"# %d",&period);
            *r='\0';
            while (*--r==' ') *r='\0'; /* delete tail spaces */
        }
        if (period<=0) period=1000;
        if (*msg&&cycle%period==0) {
            strsendcmd(stream,msg);
        }
        if (!*q) break;
	}
}
/* baseline length -----------------------------------------------------------*/
static double baseline_len(const rtk_t *rtk)
{
	double dr[3];
	int i;

	if (norm(rtk->sol.rr,3)<=0.0||norm(rtk->rb,3)<=0.0) return 0.0;

	for (i=0;i<3;i++) {
		dr[i]=rtk->sol.rr[i]-rtk->rb[i];
	}
	return norm(dr,3)*0.001; /* (km) */
}
/* send nmea request to base/nrtk input stream -------------------------------*/
static void send_nmea(rtksvr_t *svr, unsigned int *tickreset)
{
	sol_t sol_nmea={{0}};
	double vel,bl;
	unsigned int tick=tickget();
	int i;

	if (svr->stream[1].state!=1) return;

	if (svr->nmeareq==1) { /* lat-lon-hgt mode */
		sol_nmea.stat=SOLQ_SINGLE;
		sol_nmea.time=utc2gpst(timeget());
		matcpy(sol_nmea.rr,svr->nmeapos,3,1);
		strsendnmea(svr->stream+1,&sol_nmea);
	}
	else if (svr->nmeareq==2) { /* single-solution mode */
		if (norm(svr->rtk.sol.rr,3)<=0.0) return;
		sol_nmea.stat=SOLQ_SINGLE;
		sol_nmea.time=utc2gpst(timeget());
		matcpy(sol_nmea.rr,svr->rtk.sol.rr,3,1);
		strsendnmea(svr->stream+1,&sol_nmea);
	}
	else if (svr->nmeareq==3) { /* reset-and-single-sol mode */

		/* send reset command if baseline over threshold */
		bl=baseline_len(&svr->rtk);
		if (bl>=svr->bl_reset&&(int)(tick-*tickreset)>MIN_INT_RESET) {
			strsendcmd(svr->stream+1,svr->cmd_reset);
			
			tracet(2,"send reset: bl=%.3f rr=%.3f %.3f %.3f rb=%.3f %.3f %.3f\n",
				   bl,svr->rtk.sol.rr[0],svr->rtk.sol.rr[1],svr->rtk.sol.rr[2],
				   svr->rtk.rb[0],svr->rtk.rb[1],svr->rtk.rb[2]);
			*tickreset=tick;
		}
		if (norm(svr->rtk.sol.rr,3)<=0.0) return;
		sol_nmea.stat=SOLQ_SINGLE;
		sol_nmea.time=utc2gpst(timeget());
		matcpy(sol_nmea.rr,svr->rtk.sol.rr,3,1);

		/* set predicted position if velocity > 36km/h */
		if ((vel=norm(svr->rtk.sol.rr+3,3))>10.0) {
			for (i=0;i<3;i++) {
				sol_nmea.rr[i]+=svr->rtk.sol.rr[i+3]/vel*svr->bl_reset*0.8;
			}
		}
		strsendnmea(svr->stream+1,&sol_nmea);

		tracet(3,"send nmea: rr=%.3f %.3f %.3f\n",sol_nmea.rr[0],sol_nmea.rr[1],
			   sol_nmea.rr[2]);
	}
}
/* output solution------------------------------------------------------------*/
static void outrslt(rtksvr_t *svr,gmea_t *gm,int tick,int index)
{
    double tt;
    insstate_t *ins=&svr->rtk.ins;

    trace(3,"outrslt: tick=%d\n",tick);

    /* adjust current time */
    tt=(int)(tickget()-tick)/1000.0+DTTOL;
    timeset(gpst2utc(timeadd(svr->rtk.sol.time,tt)));

    /* update ins solution status */
    if (ins->stat==INSS_LCUD&&gm) {
        ins->ns=gm->ns; ins->gstat=gm->stat;
    }
    else if (ins->stat==INSS_MECH) {
        ins->ns=0;
        ins->gstat=SOLQ_NONE;
    }
    writesol(svr,index); /* write solution */
}
/* update time difference between input stream--------------------------------*/
static void updatetimediff(rtksvr_t *svr)
{
    syn_t *syn=&svr->syn;
    
    trace(3,"updatetimediff:\n");

    if (svr->rtk.opt.mode<=PMODE_FIXED) { /* for gnss position mode */
        if (syn->nr&&syn->nb) {
            syn->dt[0]=time2gpst(syn->time[0],NULL)-time2gpst(syn->time[1],NULL);
        }
    }
    /* ins loosely coupled position mode */
    if (svr->rtk.opt.mode==PMODE_INS_LGNSS) {

        if (svr->rtk.opt.insopt.lcopt==IGCOM_USESOL&&syn->ni&&syn->ns) {
            syn->dt[0]=time2gpst(syn->time[2],NULL)-time2gpst(syn->time[3],NULL);
        }
        if (svr->rtk.opt.insopt.lcopt==IGCOM_USEOBS&&syn->ni) {
            if (syn->nr) {
                syn->dt[0]=time2gpst(syn->time[2],NULL)-time2gpst(syn->time[0],NULL);
            }
            if (syn->nb) {
                syn->dt[1]=time2gpst(syn->time[2],NULL)-time2gpst(syn->time[1],NULL);
            }
        }
        if (svr->rtk.opt.insopt.pose_aid||svr->rtk.opt.insopt.align_dualants) {
            if (syn->np&&syn->ni) {
                syn->dt[2]=time2gpst(syn->time[5],NULL)-
                           time2gpst(syn->time[2],NULL);
            }
        }
    }
    /* ins tightly coupled position mode */
    if (svr->rtk.opt.mode==PMODE_INS_TGNSS) {
        if (syn->ni&&syn->nr) {
            syn->dt[0]=time2gpst(syn->time[2],NULL)-
                       time2gpst(syn->time[0],NULL);
        }
        if (syn->ni&&syn->nb) {
            syn->dt[1]=time2gpst(syn->time[2],NULL)-
                       time2gpst(syn->time[1],NULL);
        }
        if (svr->rtk.opt.insopt.pose_aid||svr->rtk.opt.insopt.align_dualants) {
            if (syn->np&&syn->ni) {
                syn->dt[2]=time2gpst(syn->time[5],NULL)-
                           time2gpst(syn->time[2],NULL);
            }
        }
    }
    /* camera measurement data aid */
    if (syn->nm&&syn->ni) {
        syn->dt[2]=time2gpst(syn->time[2],NULL)-
                   time2gpst(syn->time[4],NULL);
    }
    /* visual odometry and ins loosely coupled */
    if (svr->rtk.opt.mode==PMODE_INS_LVO) {
        syn->dt[0]=time2gpst(syn->time[2],NULL)-
                   time2gpst(syn->time[4],NULL);
    }
}
/* update input stream supend flag--------------------------------------------*/
static int suspend(rtksvr_t *svr,int index)
{
    syn_t *syn=&svr->syn;
    int lc=svr->rtk.opt.insopt.lcopt;
    double d,dT=MAXTIMEDIFF;

    d=svr->rtk.opt.soltype?-1.0:1.0;

    if (svr->rtk.opt.mode<=PMODE_FIXED) { /* gnss position mode */
        switch (index) {
            case 0: if (d*syn->dt[0]> dT) return 1; break;
            case 1: if (d*syn->dt[0]<-dT) return 1; break;
        }
    }
    /* ins loosely coupled position mode */
    if (svr->rtk.opt.mode==PMODE_INS_LGNSS) {
        if (lc==IGCOM_USESOL) { /* for solution data */
            switch (index) {
                case 3: if (d*syn->dt[0]<-dT) return 1; break;
                case 4: if (d*syn->dt[0]> dT) return 1; break;
            }
        }
        if (lc==IGCOM_USEOBS) { /* for observation data */
            switch (index) {
                case 4: if (d*syn->dt[0]> dT) return 1;
                        if (d*syn->dt[1]> dT) return 1; break;
                case 0: if (d*syn->dt[0]<-dT) return 1; break;
                case 1: if (d*syn->dt[1]<-dT) return 1; break;
            }
        }
        if (svr->rtk.opt.insopt.pose_aid||svr->rtk.opt.insopt.align_dualants) {
            switch (index) {
                case 4: if (d*syn->dt[2]<-dT) return 1; break;
                case 6: if (d*syn->dt[2]> dT) return 1; break;
            }
        }
    }
    /* ins tightly coupled position mode */
    if (svr->rtk.opt.mode==PMODE_INS_TGNSS) {
        switch (index) {
            case 4: if (d*syn->dt[0]> dT) return 1;
                    if (d*syn->dt[1]> dT) return 1; break;
            case 0: if (d*syn->dt[0]<-dT) return 1; break;
            case 1: if (d*syn->dt[1]<-dT) return 1; break;
        }
        if (svr->rtk.opt.insopt.pose_aid||svr->rtk.opt.insopt.align_dualants) {
            switch (index) {
                case 4: if (d*syn->dt[2]<-dT) return 1; break;
                case 6: if (d*syn->dt[2]> dT) return 1; break;
            }
        }
    }
    /* visual odometry and ins loosely coupled */
    if (svr->rtk.opt.mode==PMODE_INS_LVO) {
        switch (index) {
            case 4: if (d*syn->dt[0]> dT) return 1; break;
            case 5: if (d*syn->dt[0]<-dT) return 1; break;
        }
    }
    return 0;
}
/* time alignment for rover and base observation data------------------------*/
static int rbobsalign(rtksvr_t *svr)
{
    int i,j;
    obs_t *p1=NULL,*p2=NULL;
    syn_t *psyn=&svr->syn;
    gtime_t t0,t1;

    for (i=0;i<(psyn->of[0]?MAXOBSBUF:psyn->nr)&&!svr->syn.tali[0];i++) {

        /* for base observation data */
        for (j=0;j<(psyn->of[1]?MAXOBSBUF:psyn->nb);j++) {
            p1=&svr->obs[0][i]; p2=&svr->obs[1][j];

            if (!p1->n&&!p2->n) continue;
            t0=p1->data[0].time; t1=p2->data[0].time;

            if (fabs(timediff(t0,t1))>DTTOLM) continue;
            trace(3,"rover and base observation time align ok\n");
            psyn->rover=i;
            psyn->base =j;
            psyn->tali[0]=1; return 1;
        }
        /* align ok */
        if (psyn->tali[0]) break;
    }
    return 0;
}
/* time alignment for solution and imu data----------------------------------*/
static int solimualign(rtksvr_t *svr)
{
    int i,j; double sow1,sow2;
    syn_t *syn=&svr->syn;

    for (i=0;i<(syn->of[2]?MAXIMUBUF:syn->ni)&&!syn->tali[1];i++) {
        sow1=time2gpst(svr->imu[i].time,NULL);

        /* for pvt solution data */
        for (j=0;j<(syn->of[3]?MAXSOLBUF:syn->ns);j++) {

            sow2=time2gpst(svr->pvt[j].time,&svr->week);
            if (sow1==0.0||sow2==0.0||fabs(sow1-sow2)>DTTOL) continue;

            trace(3,"imu and pvt solution time align ok\n");
            syn->imu=i;
            syn->pvt=j;
            syn->tali[1]=1; return 1;
        }
        /* alignment ok */
        if (syn->tali[1]) break;
    }
    return 0;
}
/* time alignment for observation and imu data-------------------------------*/
static int imuobsalign(rtksvr_t *svr)
{
    int i,j,k,n; double sow1,sow2,sow3;
    obs_t *p1=NULL,*p2=NULL;
    syn_t *psyn=&svr->syn;

    /* observation and imu data time alignment */
    n=psyn->of[2]?MAXIMUBUF:psyn->ni;
    p1=svr->obs[0];
    p2=svr->obs[1];

    for (i=0;i<n&&svr->syn.tali[2]!=2;i++) { /* start time alignment */

        sow1=time2gpst(svr->imu[i].time,NULL);

        /* match rover observation */
        for (j=0;j<(psyn->of[0]?MAXOBSBUF:psyn->nr);j++) {
            sow2=time2gpst(p1[j].data[0].time,NULL);
            if (p1[j].n) {
                if (fabs(sow1-sow2)>DTTOL) continue;
            }
            psyn->imu    =i;
            psyn->rover  =j;
            psyn->tali[2]=1;
            break;
        }
        /* match base observation */
        if (psyn->tali[2]==1) {
            for (k=0;k<(psyn->of[1]?MAXOBSBUF:psyn->nb);k++) {
                sow3=time2gpst(p2[k].data[0].time,NULL);
                if (p2[k].n) {
                    if (fabs(sow2-sow3)>DTTOLM) continue;
                }
                else continue;
                psyn->base   =k;
                psyn->tali[2]=2;
                break;
            }
        }
        if (psyn->tali[2]==2) {
            tracet(3,"imu and rover/base align ok\n");
            return 1;
        }
        else psyn->tali[2]=0; /* fail */
    }
    return 0;
}
/* time alignment for imu and image measurement data-------------------------*/
static int imuimgalign(rtksvr_t *svr)
{
    int i,j; double sow1,sow2;
    syn_t *psyn=&svr->syn;

    /* alignment of imu and image measurement data */
    for (i=0;i<(psyn->of[2]?MAXIMUBUF:psyn->ni)&&svr->syn.tali[3]!=1;i++) {
        sow1=time2gpst(svr->imu[i].time,NULL);

        /* match image measurement data */
        for (j=0;j<(psyn->of[4]?MAXIMGBUF:psyn->nm);j++) {
            sow2=time2gpst(svr->img[j].time,NULL);
            if (svr->img[j].flag||svr->img[j].h==0||svr->img[j].w==0) continue;
            if (fabs(sow1-sow2)>DTTOL) continue;

            psyn->imu=i;
            psyn->img=j;
            psyn->tali[3]=1; break;
        }
        if (psyn->tali[3]==1) {
            tracet(3,"imu and image align ok\n");
            return 1;
        }
        else psyn->tali[3]=0; /* fail */
    }
    return 0; /* fail */
}
/* time alignment for measurement data---------------------------------------*/
static int timealign(rtksvr_t *svr)
{
    if (svr->rtk.opt.mode<=PMODE_FIXED) return rbobsalign(svr);
    
    if (svr->rtk.opt.mode==PMODE_INS_LGNSS) {
        if (svr->rtk.opt.insopt.lcopt==IGCOM_USEOBS) return imuobsalign(svr);
        if (svr->rtk.opt.insopt.lcopt==IGCOM_USESOL) return solimualign(svr);
    }
    if (svr->rtk.opt.mode==PMODE_INS_TGNSS) return imuobsalign(svr);
    if (svr->rtk.opt.mode==PMODE_INS_LVO  ) return imuimgalign(svr);
    return 0;
}
/* motion constraint for ins states update-----------------------------------*/
static void motion(const insopt_t *opt,imud_t *imuz,insstate_t *ins,
                   imud_t *imu,int ws)
{
    static int i,nc=0,zf=0; double pos[3];

    trace(3,"motion:\n");

    /* prepare imu data for static detect */
    for (i=0;i<ws-1;i++) imuz[i]=imuz[i+1]; imuz[i]=*imu;

    /* non-holonomic constraint */
    if (opt->nhc&&(nc++>opt->nhz?nc=0,true:false)) {
        nhc(ins,opt,imu);
    }
    /* zero velocity/zero angular rate update */
    ecef2pos(ins->re,pos);
    if (opt->zvu||opt->zaru) {

        /* static imu data detector */
        zf=detstc(imuz,ws,opt,pos);
        if (opt->odo) zf|=detstatic_ODO(opt,&imu->odo);

        /* zero velocity update */
        if (zf&&opt->zvu) {
            zvu(ins,opt,imuz,1);
        }
        /* zero angular rate update */
        if (zf&&opt->zaru) zaru(ins,opt,imuz,1);
    }
}
/* rtk server thread --------------------------------------------------------*/
#ifdef WIN32
static DWORD WINAPI rtksvrthread(void *arg)
#else
static void *rtksvrthread(void *arg)
#endif
{
    rtksvr_t *svr=(rtksvr_t *)arg;
    sol_t sol={0},psol={0};
    prcopt_t *opt=&svr->rtk.opt;
    insopt_t *iopt=&opt->insopt;
    insstate_t *ins=&svr->rtk.ins;
    gmea_t gnss={0};
    imud_t *imuz=NULL;
    img_t **imgt=NULL;

    static obs obss[MAXOBS]={{0}},obsd={0};
    static imu imus={0};
    static img imgs={0};
    static pose_meas_t pose={0};
    static mag_t mag={0};

    unsigned int tick,ticknmea,tick1hz,tickreset;
    unsigned char *p,*q;
    char msg[128];
    int i,j=0,n=0,ws,fobs[7]={0},cycle,cputime,init=0,flag=0;

    tracet(3,"rtksvrthread:\n");

    svr->state=1;
    svr->tick=tickget();
    ticknmea=tick1hz=svr->tick-1000;
    tickreset=svr->tick-MIN_INT_RESET;

    /* static detect window size */
    ws=opt->insopt.zvopt.ws<=0?5:opt->insopt.zvopt.ws;
    if (!(imuz=(imud_t*)malloc(sizeof(imud_t)*ws))) {
        fprintf(stderr,"malloc error\n");
        return NULL;
    }
    if (!(imgt=(img_t**)malloc(sizeof(img_t*)))) {
        fprintf(stderr,"malloc error\n");
        return NULL;
    }
    for (cycle=0;svr->state;cycle++) {
        tick=tickget();
        
        if (svr->pause ) continue;
        if (svr->reinit) init=0;

        for (i=0;i<7;i++) {
            p=svr->buff[i]+svr->nb[i]; q=svr->buff[i]+svr->buffsize;

            /* check suspend input stream */
            if (suspend(svr,i)) continue;

            /* check solution direction */
            if (svr->raw[i].dire==1) {
                svr->nb[i]=1; continue;
            }
            /* read receiver raw/rtcm data from input stream */
            if ((n=strread(svr->stream+i,p,q-p))<=0) {
                continue;    
            }
            /* write receiver raw/rtcm data to log stream */
            strwrite(svr->stream+i+9,p,n);
            svr->nb[i]+=n;
            
            /* save peek buffer */
            rtksvrlock(svr);
            n=n<svr->buffsize-svr->npb[i]?n:svr->buffsize-svr->npb[i];
            memcpy(svr->pbuf[i]+svr->npb[i],p,n);
            svr->npb[i]+=n;
            rtksvrunlock(svr);
        }
        for (i=0;i<7;i++) {
            if (svr->format[i]==STRFMT_SP3||svr->format[i]==STRFMT_RNXCLK) {
                /* decode download file */
                decodefile(svr,i);
            }
            else {
                /* decode receiver raw/rtcm data */
                fobs[i]=decoderaw(svr,i);
            }
        }
        /* update time difference between input stream */
        updatetimediff(svr);
        
        /* time alignment for measurement data */
        if (timealign(svr)) continue;

        /* averaging single base position */
        if (fobs[1]>0&&svr->rtk.opt.rb[0]==0.0&&
            svr->rtk.opt.refpos==POSOPT_SINGLE) {
            if ((svr->rtk.opt.maxaveep<=0||svr->nave<svr->rtk.opt.maxaveep)&&
                pntpos(svr->obs[1][0].data,svr->obs[1][0].n,&svr->nav,&svr->rtk.opt,&sol,
                       NULL,NULL,NULL,msg)) {
                for (i=0;i<6;i++) svr->rtk.opt.rb[i]=sol.rr[i];
            }
        }
        /* input imu measurement data */
        if (fobs[4]) {
            imus.n=inputimu(svr,imus.data);
        }
        /* for rover and base observation data from buffer */
        if (opt->mode<=PMODE_FIXED) for (i=0;i<fobs[0]&&i<MAXOBS;i++) {

            /* input rover and base observation data */
            if (!(obss[i].n=inputobs(svr,obss[i].data))) {
                continue;
            }
            /* carrier-phase bias correction */
            if (svr->nav.nf>0) {
                corr_phase_bias_fcb(obsd.data,n,&svr->nav);
            }
            else if (!strstr(opt->pppopt,"-DIS_FCB")) {
                corr_phase_bias_ssr(obsd.data,n,&svr->nav);
            }
            rtksvrlock(svr);
            rtkpos(&svr->rtk,obss[i].data,obss[i].n,&svr->nav);
            rtksvrunlock(svr);

            /* output results */
            if (opt->mode<=PMODE_PPP_FIXED) {
                if (svr->rtk.sol.stat!=SOLQ_NONE) {
                    outrslt(svr,NULL,tick,i);
                }
            }
            /* if cpu overload,inclement obs outage counter and break */
            if ((int)(tickget()-tick)>=svr->cycle) {
                svr->prcout+=fobs[0]-i-1;
            }
        }
        /* start loosely coupled position */
        if (opt->mode==PMODE_INS_LGNSS) for (i=0;i<imus.n;i++) {

            /* rtk position by observation data */
            if (iopt->lcopt==IGCOM_USEOBS) {

                if ((obsd.n=inputobstc(svr,imus.data[i].time,obsd.data))) {

                    /* rtk positioning */
                    rtksvrlock(svr);
                    rtkpos(&svr->rtk,obsd.data,obsd.n,&svr->nav);

                    if (svr->rtk.sol.stat!=SOLQ_NONE) {
                        sol2gnss(&svr->rtk.sol,&gnss);
                        j=INSUPD_MEAS;
                    }
                    rtksvrunlock(svr);
                }
                else {
                    j=INSUPD_TIME; /* ins mechanization */
                }
            }
            /* input pvt solution data from pvt buffer */
            else if (iopt->lcopt==IGCOM_USESOL) {
                if (inputpvt(svr,imus.data[i].time,&psol)) {
                    sol2gnss(&psol,&gnss);
                    j=INSUPD_MEAS;
                }
                else {
                    /* ins mechanization */
                    j=INSUPD_TIME;
                }
            }
            if (!init) {
                flag=0; if (iopt->align_dualants) {

                    /* initial ins states from dual antennas */
                    if (inputpose(svr,imus.data[i].time,&pose)
                        &&j==INSUPD_MEAS) {
                        flag=insinitdualant(svr,&pose,&psol,&imus.data[i]);
                    }
                }
                else {
                    /* initial ins states from solutions */
                    if (j==INSUPD_MEAS) {
                        flag=insinitrt(svr,&psol,&imus.data[i]);
                    }
                }
                /* initial ins states ok */
                if (flag) {
                    if (svr->reinit) {
                        svr->rtk.ins.stat=INSS_REINIT;
                        svr->reinit=0;
                    }
                    init=1;
                    tracet(3,"ins initial ok\n");
                }
                if (!init) {
                    tracet(2,"ins still initialing\n");
                    continue;
                }
            }
            /* loosely coupled position */
            rtksvrlock(svr);
            lcigpos(iopt,imus.data+i,ins,&gnss,j);

            /* motion constraint update */
            motion(iopt,imuz,ins,&imus.data[i],ws);

            /* odometry velocity aid */
            if (iopt->odo) {
                odo(iopt,&imus.data[i],&imus.data[i].odo,ins);
            }
            /* camera visual odometry aid */
            if (iopt->usecam) {
                /* todo: add functions of visual odometry aid */
            }
            /* dual ant. or camera pose aid */
            if (iopt->pose_aid) {

                flag=inputpose(svr,imus.data[i].time,&pose);

                if (flag) {
                    posefusion(iopt,&pose,ins,INSUPD_MEAS);
                }
            }
            /* magnetometer auxiliary */
            if (iopt->magh) {
                magnetometer(ins,iopt,&mag);
            }
            rtksvrunlock(svr);

            /* output results */
            outrslt(svr,&gnss,tick,i);

            /* if cpu overload,inclement obs outage counter and break */
            if ((int)(tickget()-tick)>=svr->cycle) {
                svr->iprcout+=imus.n-i-1;
            }
        }
        /* tightly coupled position start */
        if (opt->mode==PMODE_INS_TGNSS) for (j=INSUPD_TIME,i=0;i<imus.n;
            i++,j=INSUPD_TIME) {

            /* match observation data */
            if ((obsd.n=inputobstc(svr,imus.data[i].time,obsd.data))) {
                j=INSUPD_MEAS;
            }
            /* initial ins states */
            if (!init||svr->reinit) {

                if (j==INSUPD_MEAS) {
                    if (insinirtobs(svr,obsd.data,obsd.n,&imus.data[i])) {
                        init=1;
                    }
                }
                /* set initialed flag */
                if (init) {
                    tracet(3,"ins initial states ok\n");
                    if (svr->reinit) {
                        svr->rtk.ins.stat=INSS_REINIT; svr->reinit=0;
                    }
                    continue;
                }
                else {
                    tracet(2,"ins still initialing\n");
                    continue;
                }
            }
            /* start tightly coupled position */
            rtksvrlock(svr);
            tcigpos(opt,obsd.data,obsd.n,&svr->nav,&imus.data[i],&svr->rtk,ins,j);
                                       
            /* motion constraint update */
            motion(iopt,imuz,ins,&imus.data[i],ws);

            /* odometry velocity aid */
            if (iopt->odo) {
                odo(iopt,&imus.data[i],&imus.data[i].odo,ins);
            }
            /* doppler measurement aid */
            if (iopt->dopp) {
                doppler(obsd.data,obsd.n,&svr->nav,opt,ins);
            }
            /* magnetometer auxiliary */
            if (iopt->magh) {
                magnetometer(ins,iopt,&mag);
            }
            rtksvrunlock(svr);

            /* output results */
            outrslt(svr,NULL,tick,i);

            /* if cpu overload, inclement obs outage counter and break */
            if ((int)(tickget()-tick)>=svr->cycle) {
                svr->iprcout+=imus.n-i-1;
            }
        }
        /* visual odometry position mode */
        if (opt->mode==PMODE_VO) {
            if (fobs[5]) {
                imgs.n=inputimage(svr,imgs.data); /* input image measurement data */
            }
            /* reset visual odometry if needed */
            if (imgs.n&&ins->vo.time.time&&timediff(imgs.data[0].time,ins->vo.time)>MAXTIMEDIFF) {
                resetmonoa();
                init=0; /* re-initial */
            }
            /* process image one by one */
            for (i=0;i<imgs.n;i++) {

                /* initial camera position and attitude in e-frame */
                if (!init) {
#if 0
                    if (initcamgps(imgs.data[i].time,imgs.data[i].rr,iopt,ins)) {
                        tracet(3,"initial camera pose ok\n");

                        init=1; goto updatevo;
                    }
                    tracet(2,"camera pose still initial\n");
                    continue;
#elif 1
                    /* use gps solutions to initial */
                    for (j=0;j<svr->syn.ns;j++) {
                        if (fabs(timediff(imgs.data[i].time,svr->pvt[j].time))>0.05/2) {
                            continue;
                        }
                        if (initcamgps(svr->pvt[j].time,svr->pvt[j].rr,iopt,ins)) {
                            init=1; goto updatevo;
                        }
                    }
                    continue;
#else
                    /* no external measurement to initial */
                    setzero(ins->vo.rc,1,3);
                    setzero(ins->vo.vc,1,3);
                    seteye(ins->vo.Cce,3);
                    init=1;
#endif
                }
updatevo:
                /* start visual odometry */
                rtksvrlock(svr);
                flag=estvo(&opt->insopt.voopt,&imgs.data[i],svr->rtk.sol.dT,svr->rtk.sol.vc,
                           svr->rtk.sol.var);

                /* update camera pose */
                if (flag) {
                    flag=updatecam(ins,iopt,svr->rtk.sol.dT,imgs.data[i].time);
                }
                svr->rtk.sol.time=imgs.data[i].time;
                svr->rtk.sol.stat=(flag?SOLQ_VO:SOLQ_NONE);
                rtksvrunlock(svr);

                /* output vo results */
                if (flag) outrslt(svr,NULL,tick,i);
            }
        }
        /* reset measurement buffer */
        imus.n=0; imgs.n=0;

        /* send null solution if no solution (1hz) */
        if (svr->rtk.sol.stat==SOLQ_NONE&&(int)(tick-tick1hz)>=1000) {
            writesol(svr,0);
            tick1hz=tick;
        }
        /* write periodic command to input stream */
        for (i=0;i<5;i++) {
            periodic_cmd(cycle*svr->cycle,svr->cmds_periodic[i],svr->stream+i);
        }
        /* send nmea request to base/nrtk input stream */
        if (svr->nmeacycle>0&&(int)(tick-ticknmea)>=svr->nmeacycle) {
            send_nmea(svr,&tickreset);
            ticknmea=tick;
        }
        if ((cputime=(int)(tickget()-tick))>0) svr->cputime=cputime;
        
        /* sleep until next cycle */
        sleepms(svr->cycle-cputime);
    }
    for (i=0;i<MAXSTRRTK;i++) rtksvrclosestr(svr,i);
    for (i=0;i<5;i++) {
        svr->nb[i]=svr->npb[i]=0;
        free(svr->buff[i]); svr->buff[i]=NULL;
        free(svr->pbuf[i]); svr->pbuf[i]=NULL;
        free_raw (svr->raw +i);
        free_rtcm(svr->rtcm+i);
    }
    for (i=0;i<2;i++) {
        svr->nsb[i]=0;
        free(svr->sbuf[i]); svr->sbuf[i]=NULL;
    }
    for (i=0;i<MAXIMGBUF;i++) {
        if (svr->img[i].data) free(svr->img[i].data);
        svr->img[i].data=NULL;
    }
    if (opt->insopt.usecam) freemonoa();
    if (opt->mode==PMODE_VO) freemonoa();
    free(imuz); free(imgt);
    return NULL;
}
/* initialize rtk server -------------------------------------------------------
* initialize rtk server
* args   : rtksvr_t *svr    IO rtk server
* return : status (0:error,1:ok)
*-----------------------------------------------------------------------------*/
extern int rtksvrinit(rtksvr_t *svr)
{
    gtime_t time0={0};
    sol_t  sol0 ={{0}};
    eph_t  eph0 ={0,-1,-1};
    geph_t geph0={0,-1};
    seph_t seph0={0};
    sigind_t sig0={0};
    int i,j;
    
    tracet(3,"rtksvrinit:\n");
    
    svr->state=svr->cycle=svr->nmeacycle=svr->nmeareq=0;
    svr->syn.ws=3;
    for (i=0;i<3;i++) svr->nmeapos[i]=0.0;
    svr->buffsize=0;
    for (i=0;i<5;i++) svr->format[i]=0;
    for (i=0;i<2;i++) svr->solopt[i]=solopt_default;
    svr->navsel=svr->nsbs=svr->nsol=0;
    rtkinit(&svr->rtk,&prcopt_default);
    for (i=0;i<5;i++) svr->nb[i]=0;
    for (i=0;i<2;i++) svr->nsb[i]=0;
    for (i=0;i<5;i++) svr->npb[i]=0;
    for (i=0;i<5;i++) svr->buff[i]=NULL;
    for (i=0;i<2;i++) svr->sbuf[i]=NULL;
    for (i=0;i<5;i++) svr->pbuf[i]=NULL;
    for (i=0;i<MAXSOLBUF;i++) svr->solbuf[i]=sol0;
    for (i=0;i<5;i++) for (j=0;j<15;j++) svr->nmsg[i][j]=0;
    for (i=0;i<3;i++) svr->ftime[i]=time0;
    for (i=0;i<3;i++) svr->files[i][0]='\0';
    svr->moni=NULL;
    svr->tick=0;
    svr->thread=0;
    svr->cputime=svr->prcout=svr->nave=0;
    for (i=0;i<3;i++) svr->rb_ave[i]=0.0;
    
    if (!(svr->nav.eph =(eph_t  *)malloc(sizeof(eph_t )*MAXSAT *2))||
        !(svr->nav.geph=(geph_t *)malloc(sizeof(geph_t)*NSATGLO*2))||
        !(svr->nav.seph=(seph_t *)malloc(sizeof(seph_t)*NSATSBS*2))) {
        tracet(1,"rtksvrinit: malloc error\n");
        return 0;
    }
    for (i=0;i<MAXSAT *2;i++) svr->nav.eph [i]=eph0;
    for (i=0;i<NSATGLO*2;i++) svr->nav.geph[i]=geph0;
    for (i=0;i<NSATSBS*2;i++) svr->nav.seph[i]=seph0;
    svr->nav.n =MAXSAT *2;
    svr->nav.ng=NSATGLO*2;
    svr->nav.ns=NSATSBS*2;
    for (i=0;i<2;i++) for (j=0;j<7;j++) svr->nav.sind[i][j]=sig0;
    for (i=0;i<MAXIMGBUF;i++) {
        for (j=0;j<3;j++) svr->img[i].dims[j]=0;
        svr->img[i].h=svr->img[i].w=0;
        svr->img[i].data=NULL;
    }
    for (i=0;i<3;i++) for (j=0;j<MAXOBSBUF;j++) {
        if (!(svr->obs[i][j].data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))) {
            tracet(1,"rtksvrinit: malloc error\n");
            return 0;
        }
    }
    for (i=0;i<7;i++) {
        memset(svr->raw +i,0,sizeof(raw_t ));
        memset(svr->rtcm+i,0,sizeof(rtcm_t));
    }
    for (i=0;i<MAXSTRRTK;i++) strinit(svr->stream+i);
    
    for (i=0;i<7;i++) *svr->cmds_periodic[i]='\0';
    *svr->cmd_reset='\0';
    svr->bl_reset=10.0;
    initlock(&svr->lock);
    
    return 1;
}
/* free rtk server -------------------------------------------------------------
* free rtk server
* args   : rtksvr_t *svr    IO rtk server
* return : none
*-----------------------------------------------------------------------------*/
extern void rtksvrfree(rtksvr_t *svr)
{
    int i,j;
    
    free(svr->nav.eph ); svr->nav.eph =NULL;
    free(svr->nav.geph); svr->nav.geph=NULL;
    free(svr->nav.seph); svr->nav.seph=NULL;
    for (i=0;i<3;i++) for (j=0;j<MAXOBSBUF;j++) {
        free(svr->obs[i][j].data); svr->obs[i][j].data=NULL;
    }
    if (svr->rtk.ins.rtkp) rtkfree((rtk_t*)svr->rtk.ins.rtkp);
    rtkfree(&svr->rtk);
}
/* lock/unlock rtk server ------------------------------------------------------
* lock/unlock rtk server
* args   : rtksvr_t *svr    IO rtk server
* return : status (1:ok 0:error)
*-----------------------------------------------------------------------------*/
extern void rtksvrlock  (rtksvr_t *svr) {lock  (&svr->lock);}
extern void rtksvrunlock(rtksvr_t *svr) {unlock(&svr->lock);}
/* start rtk server ------------------------------------------------------------
* start rtk server thread
* args   : rtksvr_t *svr    IO rtk server
*          int     cycle    I  server cycle (ms)
*          int     buffsize I  input buffer size (bytes)
*          int     *strs    I  stream types (STR_???)
*                              types[0]=input stream rover
*                              types[1]=input stream base station
*                              types[2]=input stream correction
*                              types[3]=input stream imu
*                              types[4]=input stream solution
*                              types[4]=output stream solution 1
*                              types[5]=output stream solution 2
*                              types[6]=log stream rover
*                              types[7]=log stream base station
*                              types[8]=log stream correction
*          char    **paths  I  input stream paths
*          int     *format  I  input stream formats (STRFMT_???)
*                              format[0]=input stream rover
*                              format[1]=input stream base station
*                              format[2]=input stream correction
*                              format[3]=input stream imu
*                              format[4]=input stream solution
*          int     navsel   I  navigation message select
*                              (0:rover,1:base,2:ephem,3:all)
*          char    **cmds   I  input stream start commands
*                              cmds[0]=input stream rover (NULL: no command)
*                              cmds[1]=input stream base (NULL: no command)
*                              cmds[2]=input stream corr (NULL: no command)
*                              cmds[3]=input stream imu (NULL: no comman)
*                              cmds[4]=input stream solution (NULL: no command)
*          char    **cmds_periodic I input stream periodic commands
*                              cmds[0]=input stream rover (NULL: no command)
*                              cmds[1]=input stream base (NULL: no command)
*                              cmds[2]=input stream corr (NULL: no command)
*                              cmds[3]=input stram imu (NULL: no command)
*          char    **rcvopts I receiver options
*                              rcvopt[0]=receiver option rover
*                              rcvopt[1]=receiver option base
*                              rcvopt[2]=receiver option corr
*                              rcvopt[3]=receiver option imu
*                              rcvopt[4]=receiver option sol
*          int     nmeacycle I nmea request cycle (ms) (0:no request)
*          int     nmeareq  I  nmea request type
*                              (0:no,1:base pos,2:single sol,3:reset and single)
*          double *nmeapos  I  transmitted nmea position (ecef) (m)
*          prcopt_t *prcopt I  rtk processing options
*          solopt_t *solopt I  solution options
*                              solopt[0]=solution 1 options
*                              solopt[1]=solution 2 options
*          stream_t *moni   I  monitor stream (NULL: not used)
*          char   *errmsg   O  error message
* return : status (1:ok 0:error)
*-----------------------------------------------------------------------------*/
extern int rtksvrstart(rtksvr_t *svr, int cycle, int buffsize, int *strs,
                       char **paths, int *formats, int navsel, char **cmds,
                       char **cmds_periodic, char **rcvopts, int nmeacycle,
                       int nmeareq, const double *nmeapos, prcopt_t *prcopt,
                       solopt_t *solopt, stream_t *moni, char *errmsg)
{
    gtime_t time;
    int i,j,rw;
    
    tracet(3,"rtksvrstart: cycle=%d buffsize=%d navsel=%d nmeacycle=%d nmeareq=%d\n",
           cycle,buffsize,navsel,nmeacycle,nmeareq);
    
    if (svr->state) {
        sprintf(errmsg,"server already started");
        return 0;
    }
    strinitcom();
    svr->pause=svr->reinit=0;
    svr->cycle=cycle>1?cycle:1;
    svr->nmeacycle=nmeacycle>1000?nmeacycle:1000;
    svr->nmeareq=nmeareq;
    for (i=0;i<3;i++) svr->nmeapos[i]=nmeapos[i];
    svr->buffsize=buffsize<=0?4096:buffsize;
    for (i=0;i<7;i++) svr->format[i]=formats[i];
    svr->navsel=navsel;
    svr->nsbs=0;
    svr->nsol=0;
    svr->prcout=0;
    rtkfree(&svr->rtk);
    rtkinit(&svr->rtk,prcopt);
    
    if (prcopt->initrst) { /* init averaging pos by restart */
        svr->nave=0;
        for (i=0;i<3;i++) svr->rb_ave[i]=0.0;
    }
    for (i=0;i<7;i++) { /* input/log streams */
        svr->nb[i]=svr->npb[i]=0;
        if (!(svr->buff[i]=(unsigned char*)malloc(buffsize))||
            !(svr->pbuf[i]=(unsigned char*)malloc(buffsize))) {
            tracet(1,"rtksvrstart: malloc error\n");
            sprintf(errmsg,"rtk server malloc error");
            return 0;
        }
        for (j=0;j<15;j++) svr->nmsg[i][j]=0;
        if (i<3) {
            for (j=0;j<MAXOBSBUF;j++) svr->obs[i][j].n=0;
        }
        strcpy(svr->cmds_periodic[i],!cmds_periodic[i]?"":cmds_periodic[i]);
        
        /* initialize receiver raw and rtcm control */
        init_raw(svr->raw+i,formats[i]);
        init_rtcm(svr->rtcm+i);
        
        /* set receiver and rtcm option */
        strcpy(svr->raw [i].opt,rcvopts[i]);
        strcpy(svr->rtcm[i].opt,rcvopts[i]);

        /* connect dgps corrections */
        svr->rtcm[i].dgps=svr->nav.dgps;
    }
    for (i=0;i<2;i++) { /* output peek buffer */
        if (!(svr->sbuf[i]=(unsigned char *)malloc(buffsize))) {
            tracet(1,"rtksvrstart: malloc error\n");
            sprintf(errmsg,"rtk server malloc error");
            return 0;
        }
    }
    /* set solution options */
    for (i=0;i<2;i++) {
        svr->solopt[i]=solopt[i];
    }
    /* set base station position */
    if (prcopt->refpos!=POSOPT_SINGLE) {
        for (i=0;i<6;i++) {
            svr->rtk.rb[i]=i<3?prcopt->rb[i]:0.0;
        }
    }
    updatenav(&svr->nav); /* update navigation data */
    
    /* set monitor stream */
    svr->moni=moni;

    /* open input streams */
    for (i=0;i<14;i++) {
        rw=i<7?STR_MODE_R:STR_MODE_W; if (strs[i]!=STR_FILE) rw|=STR_MODE_W;
        if (!stropen(svr->stream+i,strs[i],rw,paths[i])) {
            sprintf(errmsg,"str%d open error path=%s",i+1,paths[i]);
            for (i--;i>=0;i--) strclose(svr->stream+i);
            return 0;
        }
        /* set initial time for rtcm and raw */
        if (i<7) {
            time=utc2gpst(timeget());
            svr->raw [i].time=strs[i]==STR_FILE?strgettime(svr->stream+i):time;
            svr->rtcm[i].time=strs[i]==STR_FILE?strgettime(svr->stream+i):time;
        }
    }
    /* put option to raw struct */
    for (i=0;i<7;i++) svr->raw[i].optp=&svr->rtk.opt;

    /* put stream pointer to raw struct */
    for (i=0;i<7;i++) svr->raw[i].strp=&svr->stream[i];

    /* set direction of raw data input */
    for (i=0;i<7;i++) {
        svr->raw[i].dire=(unsigned char)svr->rtk.opt.soltype;
    }
    /* mono camera image directory */
    for (i=0;i<7;i++) {
        strcpy(svr->raw[i].monodir,prcopt->monodir);
    }
    /* sync input streams */
    strsync(svr->stream,svr->stream+1);
    strsync(svr->stream,svr->stream+2);
    strsync(svr->stream,svr->stream+3);
    strsync(svr->stream,svr->stream+4);
    
    /* write start commands to input streams */
    for (i=0;i<7;i++) {
        if (!cmds[i]) continue;
        strwrite(svr->stream+i,(unsigned char *)"",0); /* for connect */
        sleepms(100);
        strsendcmd(svr->stream+i,cmds[i]);
    }
    /* write solution header to solution streams */
    for (i=7;i<9;i++) {
        writesolhead(svr->stream+i,svr->solopt+i-7);
    }
    /* create rtk server thread */
#ifdef WIN32
    if (!(svr->thread=CreateThread(NULL,0,rtksvrthread,svr,0,NULL))) {
#else
    if (pthread_create(&svr->thread,NULL,rtksvrthread,svr)) {
#endif
        for (i=0;i<MAXSTRRTK;i++) strclose(svr->stream+i);
        sprintf(errmsg,"thread create error\n");
        return 0;
    }
    return 1;
}
/* stop rtk server -------------------------------------------------------------
* start rtk server thread
* args   : rtksvr_t *svr    IO rtk server
*          char    **cmds   I  input stream stop commands
*                              cmds[0]=input stream rover (NULL: no command)
*                              cmds[1]=input stream base  (NULL: no command)
*                              cmds[2]=input stream ephem (NULL: no command)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtksvrstop(rtksvr_t *svr, char **cmds)
{
    int i;
    
    tracet(3,"rtksvrstop:\n");
    
    /* write stop commands to input streams */
    rtksvrlock(svr);
    for (i=0;i<3;i++) {
        if (cmds[i]) strsendcmd(svr->stream+i,cmds[i]);
    }
    rtksvrunlock(svr);
    
    /* stop rtk server */
    svr->state=0;
    
    /* free rtk server thread */
#ifdef WIN32
    WaitForSingleObject(svr->thread,10000);
    CloseHandle(svr->thread);
#else
    pthread_join(svr->thread,NULL);
#endif
}
/* open output/log stream ------------------------------------------------------
* open output/log stream
* args   : rtksvr_t *svr    IO rtk server
*          int     index    I  output/log stream index
*                              (3:solution 1,4:solution 2,5:log rover,
*                               6:log base station,7:log correction)
*          int     str      I  output/log stream types (STR_???)
*          char    *path    I  output/log stream path
*          solopt_t *solopt I  solution options
* return : status (1:ok 0:error)
*-----------------------------------------------------------------------------*/
extern int rtksvropenstr(rtksvr_t *svr, int index, int str, const char *path,
                         const solopt_t *solopt)
{
    tracet(3,"rtksvropenstr: index=%d str=%d path=%s\n",index,str,path);
    
    if (index<3||index>7||!svr->state) return 0;
    
    rtksvrlock(svr);
    
    if (svr->stream[index].state>0) {
        rtksvrunlock(svr);
        return 0;
    }
    if (!stropen(svr->stream+index,str,STR_MODE_W,path)) {
        tracet(2,"stream open error: index=%d\n",index);
        rtksvrunlock(svr);
        return 0;
    }
    if (index<=4) {
        svr->solopt[index-3]=*solopt;
        
        /* write solution header to solution stream */
        writesolhead(svr->stream+index,svr->solopt+index-3);
    }
    rtksvrunlock(svr);
    return 1;
}
/* close output/log stream -----------------------------------------------------
* close output/log stream
* args   : rtksvr_t *svr    IO rtk server
*          int     index    I  output/log stream index
*                              (3:solution 1,4:solution 2,5:log rover,
*                               6:log base station,7:log correction)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtksvrclosestr(rtksvr_t *svr, int index)
{
    tracet(3,"rtksvrclosestr: index=%d\n",index);
    
    if (!svr->state) return;
    
    rtksvrlock(svr);
    strclose(svr->stream+index);
    rtksvrunlock(svr);
}
/* get observation data status -------------------------------------------------
* get current observation data status
* args   : rtksvr_t *svr    I  rtk server
*          int     rcv      I  receiver (0:rover,1:base,2:ephem)
*          gtime_t *time    O  time of observation data
*          int     *sat     O  satellite prn numbers
*          double  *az      O  satellite azimuth angles (rad)
*          double  *el      O  satellite elevation angles (rad)
*          int     **snr    O  satellite snr for each freq (dBHz)
*                              snr[i][j] = sat i freq j snr
*          int     *vsat    O  valid satellite flag
* return : number of satellites
*-----------------------------------------------------------------------------*/
extern int rtksvrostat(rtksvr_t *svr, int rcv, gtime_t *time, int *sat,
                       double *az, double *el, int **snr, int *vsat)
{
    int i,j,ns;
    
    tracet(4,"rtksvrostat: rcv=%d\n",rcv);
    
    if (!svr->state) return 0;
    rtksvrlock(svr);
    ns=svr->obs[rcv][0].n;
    if (ns>0) {
        *time=svr->obs[rcv][0].data[0].time;
    }
    for (i=0;i<ns;i++) {
        sat [i]=svr->obs[rcv][0].data[i].sat;
        az  [i]=svr->rtk.ssat[sat[i]-1].azel[0];
        el  [i]=svr->rtk.ssat[sat[i]-1].azel[1];
        for (j=0;j<NFREQ;j++) {
            snr[i][j]=(int)(svr->obs[rcv][0].data[i].SNR[j]*0.25);
        }
        if (svr->rtk.sol.stat==SOLQ_NONE||svr->rtk.sol.stat==SOLQ_SINGLE) {
            vsat[i]=svr->rtk.ssat[sat[i]-1].vs;
        }
        else {
            vsat[i]=svr->rtk.ssat[sat[i]-1].vsat[0];
        }
    }
    rtksvrunlock(svr);
    return ns;
}
/* get stream status -----------------------------------------------------------
* get current stream status
* args   : rtksvr_t *svr    I  rtk server
*          int     *sstat   O  status of streams
*          char    *msg     O  status messages
* return : none
*-----------------------------------------------------------------------------*/
extern void rtksvrsstat(rtksvr_t *svr, int *sstat, char *msg)
{
    int i;
    char s[MAXSTRMSG],*p=msg;
    
    tracet(4,"rtksvrsstat:\n");
    
    rtksvrlock(svr);
    for (i=0;i<MAXSTRRTK;i++) {
        sstat[i]=strstat(svr->stream+i,s);
        if (*s) p+=sprintf(p,"(%d) %s ",i+1,s);
    }
    rtksvrunlock(svr);
}
/* mark current position -------------------------------------------------------
* open output/log stream
* args   : rtksvr_t *svr    IO rtk server
*          char    *name    I  marker name
*          char    *comment I  comment string
* return : status (1:ok 0:error)
*-----------------------------------------------------------------------------*/
extern int rtksvrmark(rtksvr_t *svr, const char *name, const char *comment)
{
    char buff[MAXSOLMSG+1],tstr[32],*p,*q;
    double tow,pos[3];
    int i,sum,week;
    
    tracet(4,"rtksvrmark:name=%s comment=%s\n",name,comment);
    
    if (!svr->state) return 0;
    
    rtksvrlock(svr);
    
    time2str(svr->rtk.sol.time,tstr,3);
    tow=time2gpst(svr->rtk.sol.time,&week);
    ecef2pos(svr->rtk.sol.rr,pos);
    
    for (i=0;i<2;i++) {
        p=buff;
        if (svr->solopt[i].posf==SOLF_STAT) {
            p+=sprintf(p,"$MARK,%d,%.3f,%d,%.4f,%.4f,%.4f,%s,%s\n",week,tow,
                       svr->rtk.sol.stat,svr->rtk.sol.rr[0],svr->rtk.sol.rr[1],
                       svr->rtk.sol.rr[2],name,comment);
        }
        else if (svr->solopt[i].posf==SOLF_NMEA) {
            p+=sprintf(p,"$GPTXT,01,01,02,MARK:%s,%s,%.9f,%.9f,%.4f,%d,%s",
                       name,tstr,pos[0]*R2D,pos[1]*R2D,pos[2],svr->rtk.sol.stat,
                       comment);
            for (q=(char *)buff+1,sum=0;*q;q++) sum^=*q; /* check-sum */
            p+=sprintf(p,"*%02X%c%c",sum,0x0D,0x0A);
        }
        else {
            p+=sprintf(p,"%s MARK: %s,%s,%.9f,%.9f,%.4f,%d,%s\n",COMMENTH,
                       name,tstr,pos[0]*R2D,pos[1]*R2D,pos[2],svr->rtk.sol.stat,
                       comment);
        }
        strwrite(svr->stream+i+3,(unsigned char *)buff,p-buff);
        saveoutbuf(svr,(unsigned char *)buff,p-buff,i);
    }
    if (svr->moni) {
        p=buff;
        p+=sprintf(p,"%s MARK: %s,%s,%.9f,%.9f,%.4f,%d,%s\n",COMMENTH,
                   name,tstr,pos[0]*R2D,pos[1]*R2D,pos[2],svr->rtk.sol.stat,
                   comment);
        strwrite(svr->moni,(unsigned char *)buff,p-buff);
    }
    rtksvrunlock(svr);
    return 1;
}
