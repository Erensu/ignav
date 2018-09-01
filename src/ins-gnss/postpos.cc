/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/05/08  1.0  new
*           2008/06/16  1.1  support binary inputs
*           2009/01/02  1.2  support new rtk positioing api
*           2009/09/03  1.3  fix bug on combined mode of moving-baseline
*           2009/12/04  1.4  fix bug on obs data buffer overflow
*           2010/07/26  1.5  support ppp-kinematic and ppp-static
*                            support multiple sessions
*                            support sbas positioning
*                            changed api:
*                                postpos()
*                            deleted api:
*                                postposopt()
*           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
*           2010/12/09  1.7  support qzss lex and ssr corrections
*           2011/02/07  1.8  fix bug on sbas navigation data conflict
*           2011/03/22  1.9  add function reading g_tec file
*           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
*           2011/09/15  1.11 add function reading stec file
*           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
*           2013/03/11  1.13 add function reading otl and erp data
*           2014/06/29  1.14 fix problem on overflow of # of satellites
*           2015/03/23  1.15 fix bug on ant type replacement by rinex header
*                            fix bug on combined filter for moving-base mode
*           2015/04/29  1.16 fix bug on reading rtcm ssr corrections
*                            add function to read satellite fcb
*                            add function to read stec and troposphere file
*                            add keyword replacement in dcb, erp and ionos file
*           2015/11/13  1.17 add support of L5 antenna phase center paramters
*                            add *.stec and *.trp file for ppp correction
*           2015/11/26  1.18 support opt->freqopt(disable L2)
*           2016/01/12  1.19 add carrier-phase bias correction by ssr
*           2016/07/31  1.20 fix error message problem in rnx2rtkp
*           2016/08/29  1.21 suppress warnings
*           2016/10/10  1.22 fix bug on identification of file fopt->blq
*           2017/06/13  1.23 add smoother of velocity solution
*-----------------------------------------------------------------------------*/
#include <navlib.h>

#define MAXPRCDAYS  100          /* max days of continuous processing */
#define MAXINFILE   1000         /* max number of input files */
#define MINEXPIRE   200

/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss={0};        /* receiver antenna parameters */
static pcvs_t pcvsr={0};        /* satellite antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sbs_t sbss={0};          /* sbas messages */
static lex_t lexs={0};          /* lex messages */
static sta_t stas[MAXRCV];      /* station infomation */
static imu_t imu={0};           /* imu measurement data */
static gsof_data_t gsof={0};    /* gsof measurement data for ins-gnss coupled */
static int nepoch=0;            /* number of observation epochs */
static int nimu=0;              /* number of imu measurements epochs */
static int ngsof=0;             /* number of gsof measurement data */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int isbs  =0;            /* current sbas message index */
static int ilex  =0;            /* current lex message index */
static int igsof =0;            /* current gsof message index */
static int iimu  =0;            /* current imu measurement data */
static int revs  =0;            /* analysis direction (0:forward,1:backward) */
static int aborts=0;            /* abort status */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */
static char rtcm_file[1024]=""; /* rtcm data file */
static char rtcm_path[1024]=""; /* rtcm data path */
static rtcm_t rtcm;             /* rtcm control struct */
static FILE *fp_rtcm=NULL;      /* rtcm data file pointer */

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3];
    const char *sep=opt->sep;
    
    trace(3,"outrpos :\n");
    
    if (opt->posf==SOLF_LLH||opt->posf==SOLF_ENU) {
        ecef2pos(r,pos);
        if (opt->degf) {
            deg2dms(pos[0]*R2D,dms1,5);
            deg2dms(pos[1]*R2D,dms2,5);
            fprintf(fp,"%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
                    dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],
                    sep,dms2[2],sep,pos[2]);
        }
        else {
            fprintf(fp,"%13.9f%s%14.9f%s%10.4f",pos[0]*R2D,sep,pos[1]*R2D,
                    sep,pos[2]);
        }
    }
    else if (opt->posf==SOLF_XYZ) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];
    
    trace(3,"outheader: n=%d\n",n);
    
    if (sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        return;
    }
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : RTKLIB ver.%s\n",COMMENTH,VER_RTKLIB);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
        for (i=0;i<n;i++) {
            fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
        }
        for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
        for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;

        if (popt->mode<PMODE_INS_UPDATE) {
            if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        }
        ts=popt->mode<PMODE_INS_UPDATE?obss.data[i].time:imu.data[0      ].time;
        te=popt->mode<PMODE_INS_UPDATE?obss.data[j].time:imu.data[imu.n-1].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) ts=gpst2utc(ts);
        if (sopt->times>=1) te=gpst2utc(te);
        if (sopt->times==2) ts=timeadd(ts,9*3600.0);
        if (sopt->times==2) te=timeadd(te,9*3600.0);
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }
    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED&&popt->mode!=PMODE_MOVEB) {
        fprintf(fp,"%s ref pos   :",COMMENTH);
        outrpos(fp,popt->rb,sopt);
        fprintf(fp,"\n");
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);
    
    outsolhead(fp,sopt,&popt->insopt);
}
/* search next observation data index ----------------------------------------*/
extern int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
extern int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;
    
    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* convert time to index of imu measurement data------------------------------*/
static int time2index(const gtime_t t,const imu_t *imu)
{
    int i; for (i=0;i<imu->n;i++) {
        if (fabs(timediff(imu->data[i].time,t))<DTTOL) return i;
    }
    return -1;
}
/* update rtcm ssr correction ------------------------------------------------*/
static void update_rtcm_ssr(gtime_t time)
{
    char path[1024];
    int i;
    
    /* open or swap rtcm file */
    reppath(rtcm_file,path,time,"","");
    
    if (strcmp(path,rtcm_path)) {
        strcpy(rtcm_path,path);
        
        if (fp_rtcm) fclose(fp_rtcm);
        fp_rtcm=fopen(path,"rb");
        if (fp_rtcm) {
            rtcm.time=time;
            input_rtcm3f(&rtcm,fp_rtcm);
            trace(2,"rtcm file open: %s\n",path);
        }
    }
    if (!fp_rtcm) return;
    
    /* read rtcm file until current time */
    while (timediff(rtcm.time,time)<1E-3) {
        if (input_rtcm3f(&rtcm,fp_rtcm)<-1) break;
        
        /* update ssr corrections */
        for (i=0;i<MAXSAT;i++) {
            if (!rtcm.ssr[i].update||
                rtcm.ssr[i].iod[0]!=rtcm.ssr[i].iod[1]||
                timediff(time,rtcm.ssr[i].t0[0])<-1E-3) continue;
            navs.ssr[i]=rtcm.ssr[i];
            rtcm.ssr[i].update=0;
        }
    }
}
/* synchronization of imu and gsof data---------------------------------------*/
static int sysncimugsof(const imu_t *imu,const gsof_data_t *gsof,int direction)
{
    int i=0,j=0,flag=0;

    for (direction==0?i=0:imu->n-1;direction==0?i<imu->n:i>=0;
         direction==0?i++:i--) {
        for (direction==0?j=0:j=gsof->n;direction==0?j<gsof->n:j>=0;
             direction==0?j++:j--) {
            if (fabs(timediff(imu->data[i].time,
                              gsof->data[j].t))<DTTOL) {
                flag=1; break;
            }
        }
        if (flag) break;
    }
    return direction==0?
           (iimu=i)<imu->n&&(igsof=j)<gsof->n:
           (iimu=i)>=0&&(igsof=j)>=0;
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs,int solq, const prcopt_t *popt)
{
    gtime_t time={0};
    int i,nu,nr,n=0;
    
    trace(3,"infunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n",revs,iobsu,iobsr,isbs);
    
    if (0<=iobsu&&iobsu<obss.n) {
        settime((time=obss.data[iobsu].time));
        if (checkbrk("processing : %s Q=%d",time_str(time,0),solq)) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }
    if (!revs) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
        }
        nr=nextobsf(&obss,&iobsr,2);
        if (nr<=0) {
            nr=nextobsf(&obss,&iobsr,2);
        }
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];
        iobsu+=nu;
        
        /* update sbas corrections */
        while (isbs<sbss.n) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            isbs++;
        }
        /* update lex corrections */
        while (ilex<lexs.n) {
            if (lexupdatecorr(lexs.msgs+ilex,&navs,&time)) {
                if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            }
            ilex++;
        }
        /* update rtcm ssr corrections */
        if (*rtcm_file) {
            update_rtcm_ssr(obs[0].time);
        }
    }
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        }
        else {
            for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)<-DTTOL) break;
        }
        nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr-nr+1+i];
        iobsu-=nu;
        
        /* update sbas corrections */
        while (isbs>=0) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);
            
            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            isbs--;
        }
        /* update lex corrections */
        while (ilex>=0) {
            if (lexupdatecorr(lexs.msgs+ilex,&navs,&time)) {
                if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            }
            ilex--;
        }
    }
    return n;
}
/* input imu measurement data-------------------------------------------------*/
static int inputimu(imud_t *imudata,const prcopt_t* opt,imud_t *imuz,int ws)
{
    int i,k;
    gtime_t time;

    if (0<=iimu&&iimu<imu.n) {
        settime((time=imu.data[iimu].time));
        if (checkbrk("imu measurement data time=%s",time_str(time,4))) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }
    else return 0; /* no imu measurement data */

    /* prepare imu data for static detect */
    for (k=0,i=iimu;k<ws&&i>=0&&i<imu.n;
         k++,revs?i--:i++) {
        imuz[k]=imu.data[i];
    }
    if (!revs) *imudata=imu.data[iimu++]; /* forward input */
    else       *imudata=imu.data[iimu--]; /* backward input */
    return 1;
}
/* input gsof message data----------------------------------------------------*/
static int inputgsof(gsof_t *gsofdata,gtime_t timu,const prcopt_t* opt)
{
    int flag=0,i;
    static double ts=0.5/opt->insopt.hz;
    static gsof_t gsof0={0};
    gtime_t time;

    trace(3,"inputimudata:\n");

    *gsofdata=gsof0;

    if (0<=igsof&&igsof<gsof.n) {
        settime((time=gsof.data[igsof].t));
        if (checkbrk("gsof measurement data time=%s",time_str(time,4))) {
            aborts=1;
            showmsg("aborted"); return -1;
        }
    }
    else return 0;

    /* find the closest gsof message for imu data */
    if (!revs) { /* forward */
        for (i=igsof>5?igsof-5:0;i<gsof.n;i++) {
            if (fabs(timediff(timu,gsof.data[i].t))<ts) {
                flag=1;
                igsof=i;
                *gsofdata=gsof.data[i]; break;
            }
        }
    }
    else { /* backward */

    }
    if (!flag) return 0; return 1;
}
/* exclude imu measurement data-----------------------------------------------*/
extern int excluimudata(const prcopt_t *popt,const imud_t *data)
{
    int i;

    for (i=0;i<16;i++) {
        if (popt->insopt.ext[i][0].time==0.0&&
            popt->insopt.ext[i][1].time==0.0)
            continue;
        if (timediff(data->time,popt->insopt.ext[i][0])>=0&&
            timediff(data->time,popt->insopt.ext[i][1])<=0)
            return 1;
    }
    return 0;
}
/* exclude gsof message data--------------------------------------------------*/
static int outagegsof(const prcopt_t *popt, const gsof_t *data)
{
    int i;

    for (i=0;i<16;i++) {
        if (popt->ext[i][0].time==0.0&&
            popt->ext[i][1].time==0.0)
            continue;
        if (timediff(data->t,popt->insopt.ext[i][0])>=0&&
            timediff(data->t,popt->insopt.ext[i][1])<=0)
            return 1;
    }
    return 0;
}
/* gsof message convert to gnss measurement data------------------------------*/
static int gsof2gnss(const gsof_t *gsofs, gmea_t *gnss_meas)
{
    double std[3],Cne[9];

    gnss_meas->t=gsofs->t;
    
    pos2ecef(gsofs->llh,gnss_meas->pe);
    ned2xyz (gsofs->llh,Cne);
    matmul("NN",3,1,3,1.0,Cne,gsofs->vel,0.0,gnss_meas->ve);

    std[0]=gsofs->sig[1];
    std[1]=gsofs->sig[2];
    std[2]=gsofs->sig[4];
    enu2ecef(gsofs->llh,std,gnss_meas->std);

    return norm(gnss_meas->std,3)!=0.0&&
           norm(gnss_meas->pe ,3)!=0.0;
}
/* using zero velecity imu measurement for ins corse alignment----------------
 * args  :  prcopt_t *popt  I  ins options
 *          insstate_t *ins I  ins state
 *          int *inds       O  zero velocity time start index in imu data
 *          int *inde       O  zero velocity time end index in imu data
 * return : 1: detected,0: no detected
 * --------------------------------------------------------------------------*/
static int detzerovel(const prcopt_t *popt,insstate_t *ins,int *inds,int *inde)
{
    int i,j;
    spana_t span={0};

    trace(3,"detzerovel:\n");

    if (!detstatic(&imu,ins,&popt->insopt,
                  &span,popt->insopt.zvopt.sp)) {
        trace(2,"no zero velocity measurement \n");
        return 0;
    }
    /* choose max time-span of zero velocity */
    for (j=0,i=0;i<span.n;i++) {
        if (span.tt[i].n==0) continue;
        if (timediff(span.tt[i].te,span.tt[i].ts)>=
            timediff(span.tt[j].te,span.tt[j].ts)) {
            j=i;
        }
    }
    /* check synchronization */
    *inds=time2index(span.tt[j].ts,&imu);
    *inde=time2index(span.tt[j].te,&imu);

    if (*inde<=iimu+MINEXPIRE) {
        trace(2,"check synchronization failed \n");
        return 0;
    }
    return 1;
}
/* adjust synchronization for ins and gnss measurement data-------------------*/
static void adjsync(int inde)
{
    int i; for (i=0;i<gsof.n;i++) {
        if (fabs(timediff(gsof.data[i].t,imu.data[inde].time))<DTTOL) {
            igsof=i; break;
        }
    }
    iimu=inde;
}
/* extract position from gsof message by solution status----------------------
 * args  :  double *pos  O  output position
 *          int solq     I  position status from gsof message
 *          int is       I  start index
 *          int ie       I  end index
 *          int * stat   O  position status form gsof message
 * return : number of position
 * --------------------------------------------------------------------------*/
static int gsof2pos(double *pos,int solq,int is,int ie,int *stat)
{
    int i,k;
    *stat=solq;
    for (k=0,i=is<0?0:is;i<=ie<=0?gsof.n:ie;i++) {
        if (gsof.data[i].solq==solq) matcpy(pos+3*k++,gsof.data[i].llh,3,1);
    }
    return k;
}
/* initial ins position,velecity and attitude by gsof message----------------
 * args  :  prcopt_t* opt    I  ins options
 *          insstate_t *ins  I  ins states
 *          int *gstae       O  solution status
 * return : 1 (ok) or 0 (fail)
 * note : this function only use in function: proclcgsof()
 * --------------------------------------------------------------------------*/
static int initinspva(const prcopt_t *popt,insstate_t *ins,int *gstat)
{
    int i,j,k,inds,inde,igs=-1,ige=-1;
    int sq[4]={SOLQ_FIX,SOLQ_FLOAT,SOLQ_DGPS,SOLQ_SINGLE};

    double *gr,Cne[9],grn[3]={0},gre[3]={0},gve[3]={0};
    imud_t imus={0};
    const insopt_t *opt=&popt->insopt;

    trace(3,"initinspva:\n");

    /* static alignment for initial ins states */
    if (opt->alimethod<INSALIGN_VELMATCH&&
            detzerovel(popt,ins,&inds,&inde)) {

        for (i=igsof;i<gsof.n;i++) {
            if (fabs(timediff(imu.data[inds].time,gsof.data[i].t))<DTTOL) igs=i;
            if (fabs(timediff(imu.data[inde].time,gsof.data[i].t))<DTTOL) ige=i;
        }
        gr=mat(ige-igs<=0?1:ige-igs,3);
        
        if (!((k=gsof2pos(gr,SOLQ_FIX   ,igs,ige,gstat))||
              (k=gsof2pos(gr,SOLQ_FLOAT ,igs,ige,gstat))||
              (k=gsof2pos(gr,SOLQ_DGPS  ,igs,ige,gstat))||
              (k=gsof2pos(gr,SOLQ_SINGLE,igs,ige,gstat)))) {
            return 0;
        }
        for (j=0;j<3;j++) {
            for (i=0;i<k;i++) grn[j]+=gr[3*i+j]; grn[j]/=k;
        }
        free(gr);
        
        pos2ecef(grn,gre);
        ned2xyz(grn,Cne);

        estatt(imu.data+inds,inde-inds,ins->Cbn);
        matmul3("NN",Cne,ins->Cbn,ins->Cbe);

        gapv2ipv(gre,gve,ins->Cbe,opt->lever,&imus,ins->re,ins->ve);

        /* alignment ins states using static imu data */
        switch (opt->alimethod) {
            case INSALIGN_CORSE : return coarse_align  (ins,imu.data+inds,inde-inds,opt);
            case INSALIGN_FINE  : return fine_align    (ins,imu.data+inds,inde-inds,opt);
            case INSALIGN_FINEEX: return fine_alignex  (ins,imu.data+inds,inde-inds,opt);
            case INSALIGN_LARGE : return fine_align_lym(ins,imu.data+inds,inde-inds,opt);
            default             : return coarse_align  (ins,imu.data+inds,inde-inds,opt);
        }
    }
    /* default method for initial ins states */
    else if (opt->alimethod==INSALIGN_DEFAULT) {

        /* extend method of default method*/
        if (opt->exvm) {
            for (i=0;i<4;i++) {
                if (cvmalign(&gsof,igsof,&imu,iimu,opt,sq[i],ins)) break;
            }
        }
        else { /* easy method of default method*/
            for (i=0;i<4;i++) {
                if (easyvmali(&gsof,igsof,&imu,iimu,opt,sq[i],ins)) break;
            }
        }
        *gstat=sq[i]; /* solution status */
        return i<4; /* check initial ok? */
    }
    return 0;
}
/* initial position/velocity/attitude for ins/gnss tightly coupled-----------*/
static int initcapv(const obs_t *obs,const nav_t *nav,const imu_t *imu,
                    const prcopt_t *opt,insstate_t *ins,int *iobsu,int *iobsr,
                    int *iimu)
{
    trace(3,"initc: mode=%d\n",opt->insopt.alimethod);

    switch (opt->insopt.alimethod) {
        case INSALIGN_SINGLE: return alisgpos (obs,nav,imu,opt,ins,iobsu,iimu);
        case INSALIGN_PPK   : return alirtkpos(obs,nav,imu,opt,ins,iobsu,iobsr,iimu);
        case INSALIGN_DGPS  : return alirtkpos(obs,nav,imu,opt,ins,iobsu,iobsr,iimu);
        case INSALIGN_RTK   : return alirtkpos(obs,nav,imu,opt,ins,iobsu,iobsr,iimu);
        default             : return alisgpos (obs,nav,imu,opt,ins,iobsu,iimu);
    }
    trace(2,"unknown problem arises\n");
    return 0;
}
/* find observation index for tightly coupling--------------------------------*/
static int fnobs(gtime_t imut,int *iobs,const imu_t *imu,const obs_t *obs)
{
    int i,info=0;

    if (!revs) { /* forward */
        for (i=*iobs-100<0?0:*iobs-50;i<obs->n;i++) {
            if (fabs(timediff(obs->data[i].time,imut))<DTTOL) {
                *iobs=i; info=1; break;
            }
        }
    }
    else { /* backward */
        /* todo: backward implemented */
    }
    return info;
}
/* process positioning with ins and gsof measurements data--------------------*/
static void proclcgsof(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
                        int mode)
{
    gsof_t gsofs={0};
    imud_t imus={0},*imuz;
    rtk_t rtk={{0}};
    gmea_t gnss_meas={0};
    
    int ws,stat=0,gs=SOLQ_NONE,nc=0,flag,zf=0;
    double pos[3];

    trace(3,"procinsgsof : mode=%d\n",mode);

    if (!sysncimugsof(&imu,&gsof,revs)) {
        trace(2,"synchronization of ins and gnss fail\n");
        return;
    }
    rtkinit(&rtk,popt);

    /* initial ins states before process */
    if (revs==0) { /* forward */

        if (!initinspva(popt,&rtk.ins,&gs)) {
            trace(2,"initial ins state fail\n");
            rtkfree(&rtk);
            return;
        }
    }
    else { /* backward */

    }
    /* check rtk position status */
    if (gs<popt->insopt.iisu) {
        trace(2,"initial ins state use rtk %s<%s\n",
              solqstrs[gs],solqstrs[popt->insopt.iisu]);
        rtkfree(&rtk);
        return;
    }
    /* adjust synchronization */
    adjsync(time2index(rtk.ins.time,&imu)<0?0:time2index(rtk.ins.time,&imu));

    flag=popt->insopt.zvu||popt->insopt.zaru;

    /* static detect window size */
    ws=popt->insopt.zvopt.ws<=0?5:popt->insopt.zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* process loosely coupled*/
    while (inputimu(&imus,popt,imuz,ws)) {

        if (excluimudata(popt,&imus)) continue;
        if (inputgsof(&gsofs,imus.time,popt)) {

            /* check gsof measurement data */
            if (outagegsof(popt,&gsofs)) continue;
            if (!gsof2gnss(&gsofs,&gnss_meas)) continue;

            /* ins/gnss loosely coupled */
            stat=lcigpos(&popt->insopt,&imus,&rtk.ins,&gnss_meas,INSUPD_MEAS);

            /* solution status */
            if (stat) {
                rtk.ins.ns=gsofs.ns;
                rtk.ins.gstat=gsofs.solq;
            }
            else {
                /* solution fail */
                rtk.ins.ns=0;
                rtk.ins.gstat=SOLQ_NONE;
            }
        }
        else { /* ins mechanization */
            if (!lcigpos(&popt->insopt,&imus,&rtk.ins,NULL,INSUPD_TIME)) {
                continue;
            }
        }
        /* non-holonomic constraint */
        if (popt->insopt.nhc
            &&(nc++>popt->insopt.nhz?nc=0,true:false)) {
            nhc(&rtk.ins,&popt->insopt,&imus);
        }
        /* zero velocity/zero angular rate update */
        ecef2pos(rtk.ins.re,pos);

        /* static detector */
        if (flag) {

            zf=detstc(imuz,ws,&popt->insopt,pos);
            if (popt->insopt.odo) {
                zf|=detstatic_ODO(&popt->insopt,&imus.odo);
            }
            if (zf) {
                /* zero velocity update */
                if (popt->insopt.zvu) {
                    zvu(&rtk.ins,&popt->insopt,imuz,1);
                }
                /* zero angular rate update */
                if (popt->insopt.zaru) {
                    zaru(&rtk.ins,&popt->insopt,imuz,1);
                }
            }
        }
        /* odometry velocity aid */
        if (popt->insopt.odo) {
            odo(&popt->insopt,&imus,&imus.odo,&rtk.ins);
        }
        if (mode==0) { /* forward/backward */
            outsol(fp,&rtk.sol,rtk.rb,sopt,&rtk.ins,&popt->insopt);
        }
        else if (!revs) { /* combined-forward */
            /* todo: combined-forward solutions */
        }
        else { /* combined-backward */
            /* todo: combined-backward solutions */
        }
    }
    rtkfree(&rtk);
    free(imuz);
}
/* process loosely-coupled with observation data------------------------------*/
static void proclcobs(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
                      int mode)
{
    imud_t imus={0},*imuz;
    rtk_t rtk={{0}};
    gmea_t gmeas={0};
    obsd_t obs[MAXOBS*2]; /* for rover and base */
    int ws,flag=0,nobs,i,n,stat=0,nc=0,zf=0;
    double pos[3];

    trace(3,"proclcobs:\n");

    /* initial ins states */
    rtkinit(&rtk,popt);

    /* initial ins states */
    if (!initcapv(&obss,&navs,&imu,popt,&rtk.ins,
                  &iobsu,&iobsr,&iimu)) {
        trace(2,"initial ins states fail\n");
        rtkfree(&rtk);
        return;
    }
    /* static detect window size */
    ws=popt->insopt.zvopt.ws<=0?5:popt->insopt.zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* loosely coupled process */
    while (inputimu(&imus,popt,imuz,ws)) {

        /* exclude imu data */
        if (excluimudata(popt,&imus)) continue;

        /* match observation for imu measurement data */
        flag=fnobs(imus.time,&iobsu,&imu,&obss);

        if (flag) {
            nobs=inputobs(obs,rtk.sol.stat,popt);

            if (nobs) {
                /* exclude satellites */
                for (i=n=0;i<nobs;i++) {
                    if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                        popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
                }
                if (n<=0) continue;

                /* carrier-phase bias correction */
                if (navs.nf>0) {
                    corr_phase_bias_fcb(obs,n,&navs);
                }
                else if (!strstr(popt->pppopt,"-DIS_FCB")) {
                    corr_phase_bias_ssr(obs,n,&navs);
                }
                /* rtk position */
                if (!rtkpos(&rtk,obs,n,&navs)) continue;

                matcpy(gmeas.pe,rtk.sol.rr+0,1,3);
#if 0
                matcpy(gmeas.ve,rtk.sol.rr+3,1,3);
#endif
                for (i=0;i<3;i++) gmeas.std[i]=SQRT(rtk.sol.qr[i  ]);
                for (i=3;i<6;i++) gmeas.std[i]=SQRT(rtk.sol.qv[i-3]);

                gmeas.t=rtk.sol.time;

                /* loosely-coupled start */
                stat=lcigpos(&popt->insopt,&imus,&rtk.ins,&gmeas,INSUPD_MEAS);

                /* couple ok */
                if (stat) {

                    rtk.ins.ns   =rtk.sol.ns;
                    rtk.ins.age  =rtk.sol.age;
                    rtk.ins.gstat=rtk.sol.stat;
                    rtk.ins.ratio=rtk.sol.ratio;

                    trace(3,"ins/gnss loosely coupled: time=%s\n",time_str(imus.time,4));
                }
                else {
                    /* couple fail */
                    rtk.ins.ns   =0;
                    rtk.ins.age  =0.0;
                    rtk.ins.ratio=0.0;
                    rtk.ins.gstat=SOLQ_NONE;
                }
            }
        }
        else { /* ins mechanization */
            if (!lcigpos(&popt->insopt,&imus,&rtk.ins,NULL,INSUPD_TIME)) {
                continue;
            }
            trace(3,"ins update states: time=%s\n",time_str(imus.time,3));
        }
        /* non-holonomic constraint */
        if (popt->insopt.nhc
            &&(nc++>popt->insopt.nhz?nc=0,true:false)) {
            nhc(&rtk.ins,&popt->insopt,&imus);
        }
        /* zero velocity/zero angular rate update */
        ecef2pos(rtk.ins.re,pos);

        /* static detector */
        if (flag) {

            zf=detstc(imuz,ws,&popt->insopt,pos);
            if (popt->insopt.odo) {
                zf|=detstatic_ODO(&popt->insopt,&imus.odo);
            }
            if (zf) {
                /* zero velocity update */
                if (popt->insopt.zvu) {
                    zvu(&rtk.ins,&popt->insopt,imuz,1);
                }
                /* zero angular rate update */
                if (popt->insopt.zaru) {
                    zaru(&rtk.ins,&popt->insopt,imuz,1);
                }
            }
        }
        /* odometry velocity aid */
        if (popt->insopt.odo) {
            odo(&popt->insopt,&imus,&imus.odo,&rtk.ins);
        }
        if (mode==0) { /* forward/backward */
            outsol(fp,&rtk.sol,rtk.rb,sopt,&rtk.ins,&popt->insopt);
        }
        else if (!revs) { /* combined-forward */
            /* todo: combined-forward solutions */
        }
        else { /* combined-backward */
            /* todo: combined-backward solutions */
        }
    }
    rtkfree(&rtk);
    free(imuz);
}
/* ins/gnss tighly coupled use observation------------------------------------*/
static void proctcpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
                      int mode)
{
    int i,ws,flag=0,nobs,n=0,nc=0,zf=0;
    double pos[3];
    rtk_t rtk;
    imud_t imus,*imuz;
    obsd_t obs[MAXOBS*2]; /* for rover and base */

    trace(3,"proctcpos:\n");

    rtkinit(&rtk,popt);

    /* initial ins states */
    if (!initcapv(&obss,&navs,&imu,&rtk.opt,&rtk.ins,
                  &iobsu,&iobsr,&iimu)) {
        trace(2,"initial ins states fail\n");
        rtkfree(&rtk);
        return;
    }
    /* static detect window size */
    ws=rtk.opt.insopt.zvopt.ws<=0?5:rtk.opt.insopt.zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* tightly coupled process */
    while (inputimu(&imus,popt,imuz,ws)) {

        /* exclude imu data */
        if (excluimudata(popt,&imus)) continue;

        /* match observation for imu measurement data */
        flag=fnobs(imus.time,&iobsu,&imu,&obss);

        if (flag) {
            /* observation data */
            nobs=inputobs(obs,rtk.sol.stat,popt);

            if (nobs) {
                /* exclude satellites */
                for (i=n=0;i<nobs;i++) {
                    if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                        popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
                }
                if (n<=0) continue;

                /* carrier-phase bias correction */
                if (navs.nf>0) {
                    corr_phase_bias_fcb(obs,n,&navs);
                }
                else if (!strstr(popt->pppopt,"-DIS_FCB")) {
                    corr_phase_bias_ssr(obs,n,&navs);
                }
#if 1
                /* tightly coupled */
                tcigpos(&rtk.opt,obs,n,&navs,&imus,&rtk,&rtk.ins,INSUPD_MEAS);
#endif
                /* doppler measurement aid */
                if (popt->insopt.dopp) {
                    doppler(obs,nobs,&navs,&rtk.opt,&rtk.ins);
                }
            }
        }
        else {
            /* ins mechanization */
            tcigpos(&rtk.opt,NULL,n,&navs,&imus,&rtk,&rtk.ins,INSUPD_TIME);
        }
        /* non-holonomic constraint */
        if (popt->insopt.nhc
            &&(nc++>rtk.opt.insopt.nhz?nc=0,true:false)) {
            nhc(&rtk.ins,&rtk.opt.insopt,&imus);
        }
        /* zero velocity/zero angular rate update */
        ecef2pos(rtk.ins.re,pos);

        /* static detector */
        if (flag) {

            zf=detstc(imuz,ws,&rtk.opt.insopt,pos);
            if (popt->insopt.odo) {
                zf|=detstatic_ODO(&rtk.opt.insopt,&imus.odo);
            }
            if (zf) {
                /* zero velocity update */
                if (popt->insopt.zvu) {

                    zvu(&rtk.ins,&rtk.opt.insopt,imuz,1);
                }
                /* zero angular rate update */
                if (popt->insopt.zaru) {

                    zaru(&rtk.ins,&rtk.opt.insopt,imuz,1);
                }
            }
        }
        /* odometry velocity aid */
        if (popt->insopt.odo) {
            odo(&rtk.opt.insopt,&imus,&imus.odo,&rtk.ins);
        }
        if (mode==0) { /* forward/backward */
            outsol(fp,&rtk.sol,rtk.rb,sopt,&rtk.ins,&rtk.opt.insopt);
        }
        else if (!revs) { /* combined-forward */
            /* todo: combined-forward solutions */
        }
        else { /* combined-backward */
            /* todo: combined-backward solutions */
        }
    }
    rtkfree(&rtk);
    free(imuz);
}
/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, const prcopt_t *popt, const solopt_t *sopt,
                    int mode)
{
    gtime_t time={0};
    sol_t sol={{0}};
    rtk_t rtk={{0}};
    obsd_t obs[MAXOBS*2]; /* for rover and base */
    double rb[3]={0};
    int i,nobs,n,solstatic,pri[]={0,1,2,3,4,5,1,6};
    
    trace(3,"procpos : mode=%d\n",mode);
    
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
    rtkinit(&rtk,popt);
    rtcm_path[0]='\0';
    
    while ((nobs=inputobs(obs,rtk.sol.stat,popt))>=0) {
        
        /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
        }
        if (n<=0) continue;
        
        /* carrier-phase bias correction */
        if (navs.nf>0) {
            corr_phase_bias_fcb(obs,n,&navs);
        }
        else if (!strstr(popt->pppopt,"-DIS_FCB")) {
            corr_phase_bias_ssr(obs,n,&navs);
        }
        /* disable L2 */
#if 0
        if (popt->freqopt==1) {
            for (i=0;i<n;i++) obs[i].L[1]=obs[i].P[1]=0.0;
        }
#endif
        if (!rtkpos(&rtk,obs,n,&navs)) continue;
        
        if (mode==0) { /* forward/backward */
            if (!solstatic) {
                outsol(fp,&rtk.sol,rtk.rb,sopt,NULL,&popt->insopt);
            }
            else if (time.time==0||pri[rtk.sol.stat]<=pri[sol.stat]) {
                sol=rtk.sol;
                for (i=0;i<3;i++) rb[i]=rtk.rb[i];
                if (time.time==0||timediff(rtk.sol.time,time)<0.0) {
                    time=rtk.sol.time;
                }
            }
        }
        else if (!revs) { /* combined-forward */
            if (isolf>=nepoch) return;
            solf[isolf]=rtk.sol;
            for (i=0;i<3;i++) rbf[i+isolf*3]=rtk.rb[i];
            isolf++;
        }
        else { /* combined-backward */
            if (isolb>=nepoch) return;
            solb[isolb]=rtk.sol;
            for (i=0;i<3;i++) rbb[i+isolb*3]=rtk.rb[i];
            isolb++;
        }
    }
    if (mode==0&&solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt,NULL,&popt->insopt);
    }
    rtkfree(&rtk);
}
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb)
{
    double dr[3],var[3];
    int i;
    char tstr[32];
    
    trace(3,"valcomb :\n");
    
    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr [i]=solf->rr[i]-solb->rr[i];
        var[i]=solf->qr[i]+solb->qr[i];
    }
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */
        
        time2str(solf->time,tstr,2);
        trace(2,"degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
              tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}
/* combine forward/backward solutions and output results ---------------------*/
static void combres(FILE *fp, const prcopt_t *popt, const solopt_t *sopt)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}};
    double tt,Qf[9],Qb[9],Qs[9],rbs[3]={0},rb[3]={0},rr_f[3],rr_b[3],rr_s[3];
    int i,j,k,solstatic,pri[]={0,1,2,3,4,5,1,6};
    
    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);
    
    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_PPP_STATIC);
    
    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {
        
        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;
        }
        else if (solf[i].stat<solb[j].stat) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (solf[i].stat>solb[j].stat) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        else {
            sols=solf[i];
            sols.time=timeadd(sols.time,-tt/2.0);
            
            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&
                sols.stat==SOLQ_FIX) {
                
                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j)) sols.stat=SOLQ_FLOAT;
            }
            for (k=0;k<3;k++) {
                Qf[k+k*3]=solf[i].qr[k];
                Qb[k+k*3]=solb[j].qr[k];
            }
            Qf[1]=Qf[3]=solf[i].qr[3];
            Qf[5]=Qf[7]=solf[i].qr[4];
            Qf[2]=Qf[6]=solf[i].qr[5];
            Qb[1]=Qb[3]=solb[j].qr[3];
            Qb[5]=Qb[7]=solb[j].qr[4];
            Qb[2]=Qb[6]=solb[j].qr[5];
            
            if (popt->mode==PMODE_MOVEB) {
                for (k=0;k<3;k++) rr_f[k]=solf[i].rr[k]-rbf[k+i*3];
                for (k=0;k<3;k++) rr_b[k]=solb[j].rr[k]-rbb[k+j*3];
                if (smoother(rr_f,Qf,rr_b,Qb,3,rr_s,Qs)) continue;
                for (k=0;k<3;k++) sols.rr[k]=rbs[k]+rr_s[k];
            }
            else {
                if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            }
            sols.qr[0]=(float)Qs[0];
            sols.qr[1]=(float)Qs[4];
            sols.qr[2]=(float)Qs[8];
            sols.qr[3]=(float)Qs[1];
            sols.qr[4]=(float)Qs[5];
            sols.qr[5]=(float)Qs[2];
            
            /* smoother for velocity solution */
            if (popt->dynamics) {
                for (k=0;k<3;k++) {
                    Qf[k+k*3]=solf[i].qv[k];
                    Qb[k+k*3]=solb[j].qv[k];
                }
                Qf[1]=Qf[3]=solf[i].qv[3];
                Qf[5]=Qf[7]=solf[i].qv[4];
                Qf[2]=Qf[6]=solf[i].qv[5];
                Qb[1]=Qb[3]=solb[j].qv[3];
                Qb[5]=Qb[7]=solb[j].qv[4];
                Qb[2]=Qb[6]=solb[j].qv[5];
                if (smoother(solf[i].rr+3,Qf,solb[j].rr+3,Qb,3,sols.rr+3,Qs)) continue;
                sols.qv[0]=(float)Qs[0];
                sols.qv[1]=(float)Qs[4];
                sols.qv[2]=(float)Qs[8];
                sols.qv[3]=(float)Qs[1];
                sols.qv[4]=(float)Qs[5];
                sols.qv[5]=(float)Qs[2];
            }
        }
        if (!solstatic) {
            outsol(fp,&sols,rbs,sopt,NULL,&popt->insopt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            for (k=0;k<3;k++) rb[k]=rbs[k];
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
    }
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt,NULL,&popt->insopt);
    }
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
static void readpreceph(char **infile, int n, const prcopt_t *prcopt,
                        nav_t *nav, sbs_t *sbs, lex_t *lex)
{
    seph_t seph0={0};
    int i;
    char *ext;
    
    trace(2,"readpreceph: n=%d\n",n);
    
    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    nav->nf=nav->nfmax=0;
    sbs->n =sbs->nmax =0;
    lex->n =lex->nmax =0;
    
    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readrnxc(infile[i],nav);
    }
    /* read satellite fcb files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".fcb")||!strcmp(ext,".FCB"))) {
            readfcb(infile[i],nav);
        }
    }
    /* read solution status files for ppp correction */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".stat")||!strcmp(ext,".STAT")||
             !strcmp(ext,".stec")||!strcmp(ext,".STEC")||
             !strcmp(ext,".trp" )||!strcmp(ext,".TRP" ))) {
            pppcorr_read(&nav->pppcorr,infile[i]);
        }
    }
    /* read sbas message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        sbsreadmsg(infile[i],prcopt->sbassatsel,sbs);
    }
    /* read lex message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        lexreadmsg(infile[i],0,lex);
    }
    /* allocate sbas ephemeris */
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
         showmsg("error : sbas ephem memory allocation");
         trace(1,"error : sbas ephem memory allocation");
         return;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;
    
    /* set rtcm file and initialize rtcm struct */
    rtcm_file[0]=rtcm_path[0]='\0'; fp_rtcm=NULL;
    
    for (i=0;i<n;i++) {
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
            strcpy(rtcm_file,infile[i]);
            init_rtcm(&rtcm);
            break;
        }
    }
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs, lex_t *lex)
{
    int i;
    
    trace(3,"freepreceph:\n");
    
    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->fcb ); nav->fcb =NULL; nav->nf=nav->nfmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(sbs->msgs); sbs->msgs=NULL; sbs->n =sbs->nmax =0;
    free(lex->msgs); lex->msgs=NULL; lex->n =lex->nmax =0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;
    
    if (fp_rtcm) fclose(fp_rtcm);
    free_rtcm(&rtcm);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, prcopt_t *prcopt,obs_t *obs,
                      nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;
    
    trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);
    
    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;
    
    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;
        
        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n; 
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("error : no nav data");
        trace(1,"\n");
        return 0;
    }
    /* sort observation data */
    nepoch=sortobs(obs);

    /* observation signal index for rover and base */
    for (i=0;i<2;i++) {
        for (j=0;j<7;j++) prcopt->sind[i][j]=obs->sind[i][j];
        for (j=0;j<7;j++) navs.sind[i][j]=obs->sind[i][j];
    }
    /* delete duplicated ephemeris */
    uniqnav(nav);
    
    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            settspan(ts,te);
        }
    }
    return 1;
}
/* adjust imu measurement data to frd-ned-frame and get imu time---------------
 * args    :  prcopt_t *opt  I   ins options
 *            imu_t *imu     IO  imu measurement data
 * return  : none
 * ---------------------------------------------------------------------------*/
extern void adjimudata(const prcopt_t *opt,imu_t *imu)
{
    int i,j,week,flag=0;
    double gyro[3],accl[3],dt=1.0/opt->insopt.hz,sg,si,so;

    trace(3,"adjimudata:\n");

    /* obtain gps week from gsof message for adjust imu time */
    for (i=0;i<gsof.n&&!flag;i++) {
        sg=time2gpst(gsof.data[i].t,&week);
        for (j=0;j<imu->n;j++) {
            si=time2gpst(imu->data[j].time,NULL);
            if (fabs(si-sg)<=DTTOL) {flag=1; break;}
        }
    }
    /* obtain gps week from observation data */
    for (i=0;i<imu->n&&!flag;i++) {
        si=time2gpst(imu->data[i].time,NULL);
        for (j=0;j<obss.n;j++) {
            so=time2gpst(obss.data[j].time,&week);
            if (fabs(si-so)<=DTTOL) {flag=1; break;}
        }
    }
    if (!flag) {
        trace(2,"imu and gsof measurement data synchro fail\n");
        return;
    }
    /* adjust imu data to frd-ned frame and convert to angular rate/acceleration */
    for (i=0;i<imu->n;i++) {

        /* add gps week to imu time */
        imu->data[i].time=timeadd(imu->data[i].time,week*604800.0);
        
        if (opt->insopt.imucoors==IMUCOOR_RFU) { /* convert to frd-ned-frame */
            matcpy(gyro,imu->data[i].gyro,1,3);
            matcpy(accl,imu->data[i].accl,1,3);
            matmul("NN",3,1,3,1.0,Crf,gyro,0.0,imu->data[i].gyro);
            matmul("NN",3,1,3,1.0,Crf,accl,0.0,imu->data[i].accl);
        }
        if (opt->insopt.imudecfmt==IMUDECFMT_INCR) {
            for (j=0;j<3;j++) {
                imu->data[i].gyro[j]/=dt;
                imu->data[i].accl[j]/=dt; /* convert to rate/acceleration */
            }
        }
        if (opt->insopt.imuvalfmt==IMUVALFMT_DEG) {
            for (j=0;j<3;j++) imu->data[i].gyro[j]*=D2R; /* convert to rad */
        }
    }
}
/* read imu measurements data-------------------------------------------------*/
static int readimudata(char **infile,const int *index, int n,
                       const prcopt_t *prcopt,imu_t *imu)
{
    int i;

    trace(3,"readimudata:\n");

    imu->data=NULL; imu->n=imu->nmax=0;
    
    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;

        if (!strstr(infile[i],"imu")) continue;

        /* read imu measurements data */
        if (!readimu(infile[i],imu,
                     prcopt->insopt.imudecfmt,
                     prcopt->insopt.imuformat,
                     prcopt->insopt.imucoors,
                     prcopt->insopt.imuvalfmt)&&
            !readimub(infile[i],imu,
                      prcopt->insopt.imudecfmt,
                      prcopt->insopt.imuformat,
                      prcopt->insopt.imucoors,
                      prcopt->insopt.imuvalfmt)) continue;
    }
    if (imu->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"\n");
        return 0;
    }
    /* sort imu measurement data */
    nimu=sortimudata(imu);

    /* adjust imu measurement data */
    adjimudata(prcopt,imu);
    
    return 1;
}
/* read gsof message date from file-------------------------------------------*/
static int readgsofs(char **infile,const int *index, int n,
                     const prcopt_t *prcopt,gsof_data_t *gsof)
{
    int i;

    trace(3,"readgsofs :\n");

    if (gsof->data) free(gsof->data);
    gsof->data=NULL; gsof->n=gsof->nmax=0;

    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;

        /* read gsof messages from input files */
        if (!strstr(infile[i],"gsof")) continue;
        
        if (!readgsoff(infile[i],gsof)) continue;
    }
    if (gsof->n<=0) {
        checkbrk("error : no  data");
        trace(1,"\n");
        return 0;
    }
    /* sort gsof measurement data */
    ngsof=sortgsof(gsof);
    return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");
    
    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* free gsof measurement data-------------------------------------------------*/
static void freegsofdata(gsof_data_t *gsof)
{
    trace(3,"freegsofdata:\n");
    if (gsof->data) free(gsof->data); gsof->data=NULL; gsof->n=gsof->nmax=0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];
    
    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);
    
    for (i=0;i<3;i++) ra[i]=0.0;
    
    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {
        
        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */
        
        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,NULL,msg)) continue;
        
        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, char *name, double *r)
{
    FILE *fp;
    char buff[256],sname[256],*p,*q;
    double pos[3];
    
    trace(3,"getstapos: file=%s name=%s\n",file,name);
    
    if (!(fp=fopen(file,"r"))) {
        trace(1,"station position file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if ((p=strchr(buff,'%'))) *p='\0';
        
        if (sscanf(buff,"%lf %lf %lf %s",pos,pos+1,pos+2,sname)<4) continue;
        
        for (p=sname,q=name;*p&&*q;p++,q++) {
            if (toupper((int)*p)!=toupper((int)*q)) break;
        }
        if (!*p) {
            pos[0]*=D2R;
            pos[1]*=D2R;
            pos2ecef(pos,r);
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(1,"no station position: %s %s\n",name,file);
    return 0;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;
    
    trace(3,"antpos  : rcvno=%d\n",rcvno);
    
    if (postype==POSOPT_SINGLE) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            showmsg("error : station pos computation");
            return 0;
        }
    }
    else if (postype==POSOPT_FILE) { /* read from position file */
        name=stas[rcvno==1?0:1].name;
        if (!getstapos(posfile,name,rr)) {
            showmsg("error : no position of %s in %s",name,posfile);
            return 0;
        }
    }
    else if (postype==POSOPT_RINEX) { /* get from rinex header */
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) {
            showmsg("error : no position in rinex header");
            trace(1,"no position position in rinex header\n");
            return 0;
        }
        /* antenna delta */
        if (stas[rcvno==1?0:1].deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
            del[2]+=stas[rcvno==1?0:1].hgt;
            ecef2pos(stas[rcvno==1?0:1].pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}
/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    int i;
    
    trace(3,"openses :\n");
    
    /* read satellite antenna parameters */
    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
        showmsg("error : no sat ant pcv in %s",fopt->satantp);
        trace(1,"sat antenna pcv read error: %s\n",fopt->satantp);
        return 0;
    }
    /* read receiver antenna parameters */
    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
        showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
        trace(1,"rec antenna pcv read error: %s\n",fopt->rcvantp);
        return 0;
    }
    /* use satellite L2 offset if L5 offset does not exists */
    for (i=0;i<pcvs->n;i++) {
        if (norm(pcvs->pcv[i].off[2],3)>0.0) continue;
        matcpy(pcvs->pcv[i].off[2],pcvs->pcv[i].off[1], 3,1);
        matcpy(pcvs->pcv[i].var[2],pcvs->pcv[i].var[1],19,1);
    }
    for (i=0;i<pcvr->n;i++) {
        if (norm(pcvr->pcv[i].off[2],3)>0.0) continue;
        matcpy(pcvr->pcv[i].off[2],pcvr->pcv[i].off[1], 3,1);
        matcpy(pcvr->pcv[i].var[2],pcvr->pcv[i].var[1],19,1);
    }
    return 1;
}
/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"closeses:\n");
    
    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;
    
    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;
    
    /* close solution statistics and debug trace */
    rtkclosestat();
    traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv;
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];
    
    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(3,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
        popt->pcvr[i]=*pcv;
    }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    
    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt)
{
    FILE *fp=stdout;
    
    trace(3,"outhead: outfile=%s n=%d\n",outfile,n);
    
    if (*outfile) {
        createdir(outfile);
        
        if (!(fp=fopen(outfile,"w"))) {
            showmsg("error : open output file %s",outfile);
            return 0;
        }
    }
    /* output header */
    outheader(fp,infile,n,popt,sopt);
    
    if (*outfile) fclose(fp);
    
    return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);
    
    return !*outfile?stdout:fopen(outfile,"a");
}
/* execute processing session ------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                   const solopt_t *sopt, const filopt_t *fopt, int flag,
                   char **infile, const int *index, int n, char *outfile)
{
    FILE *fp;
    prcopt_t popt_=*popt;
    char tracefile[1024],statfile[1024],path[1024];
    const char *ext;
    
    trace(3,"execses : n=%d outfile=%s\n",n,outfile);
    
    /* open debug trace */
    if (flag&&sopt->trace>0) {
        if (*outfile) {
            strcpy(tracefile,outfile);
            strcat(tracefile,".trace");
        }
        else {
            strcpy(tracefile,fopt->trace);
        }
        traceclose();
        traceopen(tracefile);
        tracelevel(sopt->trace);
    }
    /* read ionosphere data file */
    if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
        if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I')) {
            reppath(fopt->iono,path,ts,"","");
            readtec(path,&navs,1);
        }
    }
    /* read erp data */
    if (*fopt->eop) {
        free(navs.erp.data); navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;
        reppath(fopt->eop,path,ts,"","");
        if (!readerp(path,&navs.erp)) {
            showmsg("error : no erp data %s",path);
            trace(2,"no erp data %s\n",path);
        }
    }
    if (popt_.mode==PMODE_INS_LGNSS) { /* loosely coupled */

        if (popt_.insopt.lcopt==IGCOM_USEGSOF) {
            /* read gsof message data from file */
            if (!readgsofs(infile,index,n,&popt_,&gsof)) {
                trace(2,"no gsof message measurement data,ins-gnss coupled solution may fail\n");
                showmsg("error : no gsof data");
                return 0;
            }
        }
        else if (popt_.insopt.lcopt==IGCOM_USEOBS) {
            if (!readobsnav(ts,te,ti,infile,index,n,
                            &popt_,&obss,&navs,stas)) {
                showmsg("no observation or navigation data");
                return 0;
            }
        }
    }
    else if (popt_.mode==PMODE_INS_TGNSS) {  /* tightly coupled */

        if (!readobsnav(ts,te,ti,infile,index,n,
                        &popt_,&obss,&navs,stas)) {
            showmsg("no observation or navigation data");
            return 0;
        }
    }
    else if (popt_.mode<PMODE_INS_UPDATE){ /* only gnss positioning */
        /* read obs and nav data */
        if (!readobsnav(ts,te,ti,infile,index,n,
                        &popt_,&obss,&navs,stas)) {
            return 0;
        }
    }
    /* read imu measurements data */
    if (!readimudata(infile,index,n,&popt_,&imu)) {
        trace(2,"no imu measurement data,ins-gnss coupled solution is disabled\n");
    }
    /* read dcb parameters */
    if (*fopt->dcb) {
        reppath(fopt->dcb,path,ts,"","");
        readdcb(path,&navs,stas);
    }
    /* set antenna paramters */
    if (popt_.mode!=PMODE_SINGLE) {
        setpcv(obss.n>0?obss.data[0].time:timeget(),&popt_,&navs,&pcvss,&pcvsr,
               stas);
    }
    /* read ocean tide loading parameters */
    if (popt_.mode>PMODE_SINGLE&&*fopt->blq) {
        readotl(&popt_,fopt->blq,stas);
    }
    /* rover/reference fixed position */
    if (popt_.mode==PMODE_FIXED) {
        if (!antpos(&popt_,1,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            return 0;
        }
    }
    else if ((PMODE_DGPS<=popt_.mode&&popt_.mode<=PMODE_STATIC)
             ||(popt_.mode==PMODE_INS_TGNSS&&popt_.insopt.tc>=INSTC_DGPS)) {
        if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            freeimudata(&imu);
            freegsofdata(&gsof);
            return 0;
        }
    }
    /* open solution statistics */
    if (flag&&sopt->sstat>0) {
        strcpy(statfile,outfile);
        strcat(statfile,".stat");
        rtkclosestat();
        rtkopenstat(statfile,sopt->sstat);
    }
    /* write header to output file */
    if (flag&&!outhead(outfile,infile,n,&popt_,sopt)) {
        freeobsnav(&obss,&navs);
        freeimudata(&imu);
        freegsofdata(&gsof);
        return 0;
    }
    iobsu=iobsr=isbs=ilex=revs=aborts=0;
    
    if (popt_.mode==PMODE_SINGLE||popt_.soltype==0) { /* forward */
        if ((fp=openfile(outfile))) {
            if      (popt_.mode<PMODE_INS_UPDATE) procpos(fp,&popt_,sopt,0);
            else if (popt_.mode==PMODE_INS_LGNSS) {
                if      (popt_.insopt.lcopt==IGCOM_USEGSOF) proclcgsof(fp,&popt_,sopt,0);
                else if (popt_.insopt.lcopt==IGCOM_USEOBS ) proclcobs (fp,&popt_,sopt,0);
            }
            else if (popt_.mode==PMODE_INS_TGNSS) proctcpos(fp,&popt_,sopt,0);
            fclose(fp);
        }
    }
    else if (popt_.soltype==1) { /* backward */
        if ((fp=openfile(outfile))) {
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1; ilex=lexs.n-1;
            procpos(fp,&popt_,sopt,0);
            fclose(fp);
        }
    }
    else { /* combined */
        solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        rbf=(double *)malloc(sizeof(double)*nepoch*3);
        rbb=(double *)malloc(sizeof(double)*nepoch*3);
        
        if (solf&&solb) {
            isolf=isolb=0;
            procpos(NULL,&popt_,sopt,1); /* forward */
            revs=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1; ilex=lexs.n-1;
            procpos(NULL,&popt_,sopt,1); /* backward */
            
            /* combine forward/backward solutions */
            if (!aborts&&(fp=openfile(outfile))) {
                combres(fp,&popt_,sopt);
                fclose(fp);
            }
        }
        else showmsg("error : memory allocation");
        free(solf);
        free(solb);
        free(rbf);
        free(rbb);
    }
    /* free obs and nav data */
    freeobsnav(&obss,&navs);

    /* free imu measurement data and gsof data*/
    freeimudata(&imu);
    freegsofdata(&gsof);

    return aborts?1:0;
}
/* execute processing session for each rover ---------------------------------*/
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*rov_,*p,*q,s[64]="";
    
    trace(3,"execses_r: n=%d outfile=%s\n",n,outfile);
    
    for (i=0;i<n;i++) if (strstr(infile[i],"%r")) break;
    
    if (i<n) { /* include rover keywords */
        if (!(rov_=(char *)malloc(strlen(rov)+1))) return 0;
        strcpy(rov_,rov);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(rov_); for (;i>=0;i--) free(ifile[i]);
                return 0;
            }
        }
        for (p=rov_;;p=q+1) { /* for each rover */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_rov,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,p,"");
                reppath(outfile,ofile,t0,p,"");
                
                /* execute processing session */
                stat=execses(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile);
            }
            if (stat==1||!q) break;
        }
        free(rov_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        /* execute processing session */
        stat=execses(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile);
    }
    return stat;
}
/* execute processing session for each base station --------------------------*/
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov, const char *base)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*base_,*p,*q,s[64];
    
    trace(3,"execses_b: n=%d outfile=%s\n",n,outfile);
    
    /* read prec ephemeris and sbas data */
    readpreceph(infile,n,popt,&navs,&sbss,&lexs);
    
    for (i=0;i<n;i++) if (strstr(infile[i],"%b")) break;
    
    if (i<n) { /* include base station keywords */
        if (!(base_=(char *)malloc(strlen(base)+1))) {
            freepreceph(&navs,&sbss,&lexs);
            return 0;
        }
        strcpy(base_,base);
        
        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(base_); for (;i>=0;i--) free(ifile[i]);
                freepreceph(&navs,&sbss,&lexs);
                return 0;
            }
        }
        for (p=base_;;p=q+1) { /* for each base station */
            if ((q=strchr(p,' '))) *q='\0';
            
            if (*p) {
                strcpy(proc_base,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,"",p);
                reppath(outfile,ofile,t0,"",p);
                
                stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,ifile,index,n,ofile,rov);
            }
            if (stat==1||!q) break;
        }
        free(base_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile,rov);
    }
    /* free prec ephemeris and sbas data */
    freepreceph(&navs,&sbss,&lexs);
    
    return stat;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise 
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first 
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extention as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .lex,.LEX            : qzss lex message log files
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              .fcb,.FCB            : satellite fcb
*              .imu,.bin            : imu measurement data
*              .gsof                : gsof measurements data from trimble
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, char **infile, int n, char *outfile,
                   const char *rov, const char *base)
{
    gtime_t tts,tte,ttte;
    double tunit,tss;
    int i,j,k,nf,stat=0,week,flag=1,index[MAXINFILE]={0};
    char *ifile[MAXINFILE],ofile[1024],*ext;
    
    trace(3,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);
    
    /* open processing session */
    if (!openses(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;
    
    if (ts.time!=0&&te.time!=0&&tu>=0.0) {
        if (timediff(te,ts)<0.0) {
            showmsg("error : no period");
            closeses(&navs,&pcvss,&pcvsr);
            return 0;
        }
        for (i=0;i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
        }
        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;
        settspan(ts,te);
        tunit=tu<86400.0?tu:86400.0;
        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);
        
        for (i=0;;i++) { /* for each periods */
            tts=gpst2time(week,tss+i*tu);
            tte=timeadd(tts,tu-DTTOL);
            if (timediff(tts,te)>0.0) break;
            if (timediff(tts,ts)<0.0) tts=ts;
            if (timediff(tte,te)>0.0) tte=te;
            
            strcpy(proc_rov ,"");
            strcpy(proc_base,"");
            if (checkbrk("reading    : %s",time_str(tts,0))) {
                stat=1;
                break;
            }
            for (j=k=nf=0;j<n;j++) {
                
                ext=strrchr(infile[j],'.');
                
                if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
                    strcpy(ifile[nf++],infile[j]);
                }
                else {
                    /* include next day precise ephemeris or rinex brdc nav */
                    ttte=tte;
                    if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
                              !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
                        ttte=timeadd(ttte,3600.0);
                    }
                    else if (strstr(infile[j],"brdc")) {
                        ttte=timeadd(ttte,7200.0);
                    }
                    nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,tts,ttte,"","");
                }
                while (k<nf) index[k++]=j;
                
                if (nf>=MAXINFILE) {
                    trace(2,"too many input files. trancated\n");
                    break;
                }
            }
            if (!reppath(outfile,ofile,tts,"","")&&i>0) flag=0;
            
            /* execute processing session */
            stat=execses_b(tts,tte,ti,popt,sopt,fopt,flag,ifile,index,nf,ofile,
                           rov,base);
            
            if (stat==1) break;
        }
        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
    }
    else if (ts.time!=0) {
        for (i=0;i<n&&i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                return -1;
            }
            reppath(infile[i],ifile[i],ts,"","");
            index[i]=i;
        }
        reppath(outfile,ofile,ts,"","");
        
        /* execute processing session */
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,ifile,index,n,ofile,rov,
                       base);
        
        for (i=0;i<n&&i<MAXINFILE;i++) free(ifile[i]);
    }
    else {
        for (i=0;i<n;i++) index[i]=i;
        
        /* execute processing session */
        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,rov,
                       base);
    }
    /* close processing session */
    closeses(&navs,&pcvss,&pcvsr);
    return stat;
}
