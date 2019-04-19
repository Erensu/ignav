/*-----------------------------------------------------------------------------
* ins-gnss-rts.cc : ins-gnss RTS smoother
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
*        INS 2008.
*    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/09/17 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants-----------------------------------------------------------------*/
#define MAXTIMEDIFF  3.0                 /* max time difference for RTS smoother */
#define NP           3                   /* number of position solutions */
#define MINVEL       5.0                 /* min velocity for initial ins states */
#define MAXGYRO      (30.0*D2R)          /* max rotation speed value for initial */
#define MAXDIFF      10.0                /* max time difference between solution */
#define INTKEEPALIVE 1000                /* keep alive interval (ms) */
#define OUTSOLFRQ    100                 /* frequency of output ins solutions */
#define OUT_MONITOR  1                   /* output solution to monitor */
#define SOL_OUTPUT_FILE 0                /* output solution to file */
#define FORWARD_IN_MEMO 0                /* forward solution save in memory, otherwise in file */

typedef struct {                         /* RTS ins solution data type */
    gtime_t time;                        /* ins solution time */
    int nx;                              /* number of error states */
    int ns,stat,gstat;                   /* number of valid satellite/ins solution status/gps status */
    double cCbe[9],cre[3],cve[3],cae[3];
    double cba [3],cbg[3],cMa[9],cMg[9];
    double clever[3];                    /* current ins solution data */

    double pCbe[9],pre[3],pve[3],pae[3];
    double pba [3],pbg[3],pMa[9],pMg[9]; /* predicted ins solution data */
    double plever[3];

    double sCbe[9],sre[3],sve[3],sae[3];
    double sba [3],sbg[3],sMa[9],sMg[9]; /* smoother ins solution data */
    double slever[3];
    double *Pc,*Pp,*Ps;                  /* updated/predicted/smoothed error state covariance */
    double *F;                           /* transmit matrix */
} ins_sol_t;

typedef struct {                         /* ins solution buffer type */
    int n,nmax;                          /* number and max number of solutions */
    ins_sol_t  *data;                    /* ins solution data */
} ins_solbuf_t;

/* constants/global variables -----------------------------------------------*/
static stream_t    moni={0};             /* monitor stream */
static stream_t    frst={0};             /* solution result file */
static prcopt_t prcopt ={0};             /* processing options */
static solopt_t solopt ={0};             /* solution options */
static filopt_t filopt ={""};            /* file options */
static ins_solbuf_t insbuf={0};          /* ins solution buffer */
static ins_sol_t    insol={0};           /* ins solution data for temporary savings */

static int ipos=0;                       /* current gsof message index */
static int iimu=0;                       /* current imu measurement data */
static int keepalive=0;                  /* keep alive flag */
static int timeout  =10000;              /* timeout time (ms) */
static int reconnect=10000;              /* reconnect interval (ms) */
static int week=0;                       /* GPS week */
static char solfile[1024];               /* solution output file path */
static FILE *fp_fwd_sol=NULL;            /* foeward solution file pointer */

/* adjust imu measurement time ----------------------------------------------*/
static void adj_imutime(imud_t *imu,const prcopt_t *opt)
{
    int w=0;
    time2gpst(imu->time,&w);
    if (w==0) {
        imu->time=timeadd(imu->time,week*604800.0);
    }
}
/* input imu measurement data------------------------------------------------*/
static int inputimu(const imu_t *imu,imud_t *imudata,const prcopt_t* opt,
                    imud_t *imuz,int ws)
{
    int i,k;
    if (0<=iimu&&iimu<imu->n) {
        trace(3,"input imu measurement data\n");
    }
    else return 0;

    /* imu data for static detect */
    for (k=0,i=iimu-ws;k<ws&&i>=0&&i<imu->n;k++,i++) {
        /* adjust time */
        adj_imutime(&imu->data[i],opt);

        imuz[k]=imu->data[i];
    }
    adj_imutime(&imu->data[iimu],opt);

    /* forward input */
    *imudata=imu->data[iimu++];
    return 1;
}
/* input position measurement data-------------------------------------------*/
static int inputpos(const gsof_data_t *pos,gsof_t *gsofdata,gtime_t timu,
                    const prcopt_t* opt)
{
    int flag=0,i;
    static double ts=0.5/opt->insopt.hz,sowi,sowp;

    trace(3,"inputpos:\n");

    if (0<=ipos&&ipos<pos->n) {
        trace(3,"input position measurement data\n");
    }
    else return 0;

    sowi=time2gpst(timu,NULL);

    /* find the closest gsof message for imu data */
    for (i=ipos>5?ipos-5:0;i<pos->n;i++) {
        sowp=time2gpst(pos->data[i].t,NULL);
        if (fabs(sowi-sowp)<ts) {
            flag=1;
            ipos=i;
            *gsofdata=pos->data[i];

            /* GPS week */
            time2gpst(gsofdata->t,&week);
            break;
        }
    }
    if (!flag) return 0; return 1;
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
    double cov[9]={0},Cne[9],cove[9];

    if (gsofs->solq==SOLQ_NONE) return 0;
    if (gsofs->sig[2]<0.0||gsofs->sig[1]<0.0||gsofs->sig[4]<0.0) {
        return 0;
    }
    gnss_meas->t=gsofs->t;

    pos2ecef(gsofs->llh,gnss_meas->pe);
    ned2xyz (gsofs->llh,Cne);
    matmul("NN",3,1,3,1.0,Cne,gsofs->vel,0.0,gnss_meas->ve);

    cov[0]=SQR(gsofs->sig[2]==0.0?30.0:gsofs->sig[2]); /* ned */
    cov[4]=SQR(gsofs->sig[1]==0.0?30.0:gsofs->sig[1]);
    cov[8]=SQR(gsofs->sig[4]==0.0?30.0:gsofs->sig[4]);
    cov[1]=cov[3]=(gsofs->sig[3]<0.0?-SQR(gsofs->sig[3]):SQR(gsofs->sig[3]));
    
    matmul33("NNT",Cne,cov,Cne,3,3,3,3,cove);
    matcpy(gnss_meas->covp,cove,3,3);

    gnss_meas->std[0]=SQRT(cove[0]);
    gnss_meas->std[1]=SQRT(cove[4]);
    gnss_meas->std[2]=SQRT(cove[8]);

    gnss_meas->stat=(unsigned char)gsofs->solq;
    gnss_meas->ns  =(unsigned char)gsofs->ns;
    return norm(gnss_meas->std,3)!=0.0&&
           norm(gnss_meas->pe ,3)!=0.0;
}
/* check solution valid------------------------------------------------------*/
static int chksol(const gsof_t *sols)
{
    if (sols->solq==SOLQ_NONE) return 0;
    if (sols->t.time==0) return 0;
    if (norm(sols->pos,3)==0.0) return 0;
    return 1;
}
/* solution convert to velocity----------------------------------------------*/
static void sol2vel(const gsof_t *sol1,const gsof_t *sol2,double *v)
{
    v[0]=(sol1->pos[0]-sol2->pos[0])/timediff(sol1->t,sol2->t);
    v[1]=(sol1->pos[1]-sol2->pos[1])/timediff(sol1->t,sol2->t);
    v[2]=(sol1->pos[2]-sol2->pos[2])/timediff(sol1->t,sol2->t);
}
/* initial ins states--------------------------------------------------------*/
static int init_ins(const imud_t *imu,const gsof_t *pos,const insopt_t *iopt,
                    insstate_t *ins)
{
    static gsof_t sols[NP]={0};
    double vr[3]={0};
    int i;

    trace(3,"init_ins: time=%s\n",time_str(imu->time,3));

    /* check position measurement status */
    if (!chksol(pos)) {
        trace(2,"invalid position meas.\n");
        return 0;
    }
    /* save pvt solution buffer */
    for (i=0;i<NP-1;i++) sols[i]=sols[i+1]; sols[i]=*pos;
    for (i=0;i<NP;i++) {
        if (sols[i].solq>=iopt->iisu||sols[i].solq==SOLQ_NONE) return 0;
    }
    /* compute velocity from solutions */
    sol2vel(sols+NP-1,sols+NP-2,vr);

    /* check velocity ok? */
    if (norm(vr,3)<MINVEL||norm(imu->gyro,3)>MAXGYRO) {
        return 0;
    }
    /* initial ins states */
    if (!ant2inins(sols[NP-1].t,sols[NP-1].pos,vr,iopt,NULL,
                   ins,NULL)) {
        trace(2,"initial ins state fail\n");
        return 0;
    }
    ins->time=sols[NP-1].t;
    ins->stat=INSS_INIT;

    /* update ins state in n-frame */
    updinsn(ins);

    trace(3,"initial ins state ok\n");
    return 1;
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
/* thread to send keep alive for monitor port --------------------------------*/
static void *sendkeepalive(void *arg)
{
    trace(3,"sendkeepalive: start\n");

    while (keepalive) {
        strwrite(&moni,(unsigned char *)"\r",1);
        sleepms(INTKEEPALIVE);
    }
    trace(3,"sendkeepalive: stop\n");
    return NULL;
}
/* open monitor port ---------------------------------------------------------*/
static int openmoni(int port)
{
    pthread_t thread;
    char path[64];

    trace(3,"openmomi: port=%d\n",port);

    sprintf(path,":%d",port);

    if (!stropen(&moni,STR_TCPSVR,STR_MODE_RW,path)) return 0;
    strsettimeout(&moni,timeout,reconnect);
    keepalive=1;
    pthread_create(&thread,NULL,sendkeepalive,NULL);
    return 1;
}
/* close monitor port---------------------------------------------------------*/
static int close_moni(stream_t *moni)
{
    trace(3,"closemoni:\n");
    keepalive=0;

    /* send disconnect message */
    strwrite(moni,(unsigned char *)MSG_DISCONN,strlen(MSG_DISCONN));

    /* wait fin from clients */
    sleepms(1000);
    strclose(moni);
}
/* open output solution file--------------------------------------------------*/
static int open_solfile(const char *file)
{
    trace(3,"open_solfile: file=%s\n",file);
    if (!stropen(&frst,STR_FILE,STR_MODE_W,file)) return 0;
#if !FORWARD_IN_MEMO
    if (!(fp_fwd_sol=fopen(solfile,"wb"))) return 0;
#endif
    return 1;
}
/* write solution header to output stream ------------------------------------*/
static void writesolhead(stream_t *stream, const solopt_t *solopt)
{
    unsigned char buff[1024];
    int n;

    n=outsolheads(buff,solopt);
    strwrite(stream,buff,n);
}
/* time of output solutions to monitor----------------------------------------*/
static gtime_t out_time={0};
/* adjust backward solution time for plot better------------------------------*/
static void adj_bcksol_time(rtk_t *rtk)
{
    double dt;
    dt=timediff(out_time,rtk->sol.time);
    rtk->sol.time=timeadd(out_time,dt);
}
/* write solution status output stream ---------------------------------------*/
static int wrt_solution(rtk_t *rtk,const solopt_t *solopt,int type)
{
    unsigned char buff[1024];
    static int c=0,out_head=0;
    int n=0;

    if (rtk->sol.time.time==0) return 0;

#if SOL_OUTPUT_FILE
    /* output solution to file */
    n=outsols(buff,&rtk->sol,rtk->rb,solopt,&rtk->ins,&rtk->opt.insopt,0);
    if (frst.port) {
        strwrite(&frst,buff,n);
    }
#endif
    if (type==0) out_time=rtk->sol.time;
    if (!out_head) {
        writesolhead(&frst,solopt);
        out_head=1;
    }
    /* output solution to monitor */
    if (moni.port) {

        /* adjust backward solution time */
        if (type==1) {
            adj_bcksol_time(rtk);
        }
        n=outsols(buff,&rtk->sol,rtk->rb,solopt,&rtk->ins,&rtk->opt.insopt,1);
        if (c++>OUTSOLFRQ) {
            strwrite(&moni,buff,n);
            c=0;
        }
    }
    return n;
}
/* free ins states struct----------------------------------------------------*/
extern void freeins(insstate_t *ins)
{
    trace(3,"freeins:\n");
    if (ins->x ) free(ins->x ); ins->x =NULL;
    if (ins->P ) free(ins->P ); ins->P =NULL;
    if (ins->Pa) free(ins->Pa); ins->Pa=NULL;
    if (ins->Pb) free(ins->Pb); ins->Pb=NULL;
    if (ins->xb) free(ins->xb); ins->xb=NULL;
    if (ins->xa) free(ins->xa); ins->xa=NULL;
    if (ins->P0) free(ins->P0); ins->P0=NULL;
    if (ins->F ) free(ins->F ); ins->F =NULL;
    if (ins->gmeas.data) free(ins->gmeas.data); ins->gmeas.data=NULL;

    ins->nx=ins->nb=0;
    ins->gmeas.n=ins->gmeas.nmax=0;
}
/* initial ins solution temporary--------------------------------------------*/
static void init_insol(ins_sol_t *sol,int nx)
{
    gtime_t t0={0};
    sol->time=t0;
    sol->Pc=mat(nx,nx); sol->Pp=mat(nx,nx);
    sol->Ps=mat(nx,nx); sol->F =mat(nx,nx); sol->nx=nx;
}
/* free ins solution temporary-----------------------------------------------*/
static void free_insol(ins_sol_t *sol)
{
    if (sol->Pc) free(sol->Pc); sol->Pc=NULL;
    if (sol->Pp) free(sol->Pp); sol->Pp=NULL;
    if (sol->Ps) free(sol->Ps); sol->Ps=NULL;
    if (sol->F ) free(sol->F ); sol->F =NULL;
}
/* save ins solution to RTS data struct--------------------------------------*/
static void torts(ins_sol_t *sol,const insstate_t *ins,const insopt_t *opt,
                  int type)
{
    if (sol->Pc==NULL||sol->Pp==NULL) return;
    if (type==1) { /* updated/predicted ins state */
        matcpy(sol->pCbe,ins->Cbe,3,3);
        matcpy(sol->pre ,ins->re ,3,1);
        matcpy(sol->pve ,ins->ve ,3,1);
        matcpy(sol->pae ,ins->ae ,3,1);
        matcpy(sol->pba ,ins->ba ,3,1); matcpy(sol->pbg,ins->bg,3,1);
        matcpy(sol->pMa ,ins->Ma ,3,3); matcpy(sol->pMg,ins->Mg,3,3);
        matcpy(sol->Pp  ,ins->P  ,ins->nx,ins->nx);
        matcpy(sol->F   ,ins->F  ,ins->nx,ins->nx);

        matcpy(sol->plever,ins->lever,3,1);

        sol->stat=ins->stat;
    }
    else if (type==2) { /* current ins state */
        matcpy(sol->cCbe,ins->Cbe,3,3);
        matcpy(sol->cre ,ins->re ,3,1);
        matcpy(sol->cve ,ins->ve ,3,1);
        matcpy(sol->cae ,ins->ae ,3,1);
        matcpy(sol->cba ,ins->ba ,3,1); matcpy(sol->cbg,ins->bg,3,1);
        matcpy(sol->cMa ,ins->Ma ,3,3); matcpy(sol->cMg,ins->Mg,3,3);
        matcpy(sol->Pc  ,ins->P  ,ins->nx,ins->nx);
        matcpy(sol->plever,ins->lever,3,1);

        sol->time=ins->time;
    }
    sol->ns   =ins->ns;
    sol->gstat=ins->gstat;
}
/* backup ins states information---------------------------------------------
 * args:  insstate_t *ins  IO  ins states
 *        insopt_t *opt    I   ins options
 *        int type         I   information type
 * return: 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int bckupinsinfo(insstate_t *ins, const insopt_t *opt, int type)
{
    trace(3,"bckupinsinfo:\n");
    torts(&insol,ins,opt,type);
}
/* write forward solution to file in binary----------------------------------*/
static int wrt_fwdsol_bin(const ins_sol_t *data)
{
    int i;
    fwrite(data,sizeof(ins_sol_t),1,fp_fwd_sol);

    for (i=0;i<data->nx*data->nx;i++) {
        fwrite(&data->Pc[i],sizeof(double),1,fp_fwd_sol);
        fwrite(&data->Pp[i],sizeof(double),1,fp_fwd_sol);
        fwrite(&data->F [i],sizeof(double),1,fp_fwd_sol);
    }
    return 1;
}
/* add ins solutions to buffer-----------------------------------------------*/
static int add_ins_sol(ins_solbuf_t *buf,const ins_sol_t *data)
{
#if FORWARD_IN_MEMO
    int i; ins_sol_t *pins;
    if (data->time.time==0) return 0;
    if (buf->nmax<=buf->n) {
        buf->nmax+=1024;
        if (!(pins=(ins_sol_t *)realloc(buf->data,sizeof(ins_sol_t)*buf->nmax))) {
            trace(1,"add_ins_sol malloc error: n=%d\n",buf->nmax);
            for (i=0;i<buf->n;i++) free_insol(buf->data+i);
            return 0;
        }
        buf->data=pins;
    }
    pins=&buf->data[buf->n++]; *pins=*data;

    pins->Ps=zeros(data->nx,data->nx);
    pins->Pp=zeros(data->nx,data->nx);
    pins->Pc=zeros(data->nx,data->nx);
    pins->F =eye(data->nx);

    matcpy(pins->Pp,data->Pp,data->nx,data->nx);
    matcpy(pins->Ps,data->Ps,data->nx,data->nx);
    matcpy(pins->Pc,data->Pc,data->nx,data->nx);
    matcpy(pins->F ,data->F ,data->nx,data->nx);
#else
    /* write solution to file in binary */
    wrt_fwdsol_bin(data);
    insbuf.n++;
#endif
    return 1;
}
/* output solution------------------------------------------------------------*/
static void out_stat(const gmea_t *gm,insstate_t *ins)
{
    /* update ins solution status */
    if (ins->stat==INSS_LCUD&&gm) {
        ins->ns=gm->ns; ins->gstat=gm->stat;
    }
    else if (ins->stat==INSS_MECH) {
        ins->ns=0;
        ins->gstat=SOLQ_NONE;
    }
}
/* free forward/backward solution--------------------------------------------*/
static void free_fbsol(ins_solbuf_t *insbuf)
{
#if FORWARD_IN_MEMO
    int i;
    for (i=0;i<insbuf->n;i++) {
        free_insol(&insbuf->data[i]);
    }
#endif
}
/* forward filter of rts-----------------------------------------------------
 * args:    imu_t *imu        I  imu measurement data
 *          gsof_data_t *pos  I  position measurement data
 *          insopt_t *opt     I  ins options
 * return:  number of epochs
 *---------------------------------------------------------------------------*/
static int fwdfilt(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                   const solopt_t *solopt,rtk_t *rtk)
{
    imud_t imus={0},*imuz=NULL;
    gsof_t poss={0};
    gmea_t gmea={0};
    insstate_t *ins=&rtk->ins;
    const insopt_t *iopt=&rtk->opt.insopt;
    int ws,init=0;

    trace(3,"fwdfilt: ni=%d  np=%d\n",imu->n,pos->n);

    /* static detect window size */
    ws=iopt->zvopt.ws<=0?5:iopt->zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* rtk init. */
    rtkinit(rtk,popt);

    /* initial ins solution temporary */
    init_insol(&insol,rtk->ins.nx);

    while (inputimu(imu,&imus,popt,imuz,ws)) {
        if (inputpos(pos,&poss,imus.time,popt)) {

            if (!gsof2gnss(&poss,&gmea)||outagegsof(popt,&poss)) continue;
            if (!init) {
                /* initialization */
                if (!(init=init_ins(&imus,&poss,iopt,ins))) continue;
            }
            else {
                /* coupled. */
                lcigpos(iopt,&imus,ins,&gmea,INSUPD_MEAS);
            }
        }
        else if (init==1) {
            /* ins mech. */
            lcigpos(iopt,&imus,ins,NULL,INSUPD_TIME);
        }
        else continue;

        /* motion constraint update */
        motion(iopt,imuz,ins,&imus,ws);

        /* odometry velocity aid */
        if (iopt->odo) {
            odo(iopt,&imus,&imus.odo,ins);
        }
        /* sol. status */
        out_stat(&gmea,ins);

        /* solution. */
        ins2sol(&rtk->ins,iopt,&rtk->sol);

#if OUT_MONITOR
        /* write solution */
        if (!wrt_solution(rtk,solopt,0)) continue;
#endif
        /* ins solution add. */
        if (!add_ins_sol(&insbuf,&insol)) {
            trace(2,"add ins solution fail\n");
            continue;
        }
    }
#if !FORWARD_IN_MEMO
    fclose(fp_fwd_sol); fp_fwd_sol=NULL;
#endif
    free_insol(&insol);
    free(imuz);
    return insbuf.n>1;
}
/* get error correction of smoothed state------------------------------------*/
static void err_corr(const ins_sol_t *pre, ins_sol_t *cur,
                     const insopt_t *opt,double *x,double *P)
{
    double omg[3],dC[9];
    int i=0,iba,nba,ibg,nbg,isg,nsg,isa,nsa,nla,ila;

    matmul("TN",3,3,3,1.0,cur->pCbe,pre->sCbe,0.0,dC);
    so3_log(dC,omg,NULL);
    for (i=0;i<3;i++) omg[i]=-omg[i];

    iba=xiBa(opt); nba=xnBa(opt);
    ibg=xiBg(opt); nbg=xnBg(opt);
    isg=xiSg(opt); nsg=xnSg(opt);
    isa=xiSa(opt); nsa=xnSa(opt);
    ila=xiLa(opt); nla=xnLa(opt);

    matcpy(x,omg,1,3);

    /* residual scale factors of gyro. and accl. */
    if (nsg) {
        for (i=isg;i<isg+nsg;i++) {
            x[i]=pre->sMg[i-isg+(i-isg)*3]-cur->pMg[i-isg+(i-isg)*3];
        }
    }
    if (nsa) {
        for (i=isa;i<isa+nsa;i++) {
            x[i]=pre->sMa[i-isa+(i-isa)*3]-cur->pMa[i-isa+(i-isa)*3];
        }
    }
    x[3]=-pre->sve[0]+cur->pve[0];
    x[4]=-pre->sve[1]+cur->pve[1];
    x[5]=-pre->sve[2]+cur->pve[2]; /* velocity */

    x[6]=-pre->sre[0]+cur->pre[0];
    x[7]=-pre->sre[1]+cur->pre[1];
    x[8]=-pre->sre[2]+cur->pre[2]; /* position */

    /* accl. bias */
    if (nba) {
        x[iba+0]=pre->sba[0]-cur->pba[0];
        x[iba+1]=pre->sba[1]-cur->pba[1];
        x[iba+2]=pre->sba[2]-cur->pba[2];
    }
    /* gyro. bias */
    if (nbg) {
        x[ibg+0]=pre->sbg[0]-cur->pbg[0];
        x[ibg+1]=pre->sbg[1]-cur->pbg[1];
        x[ibg+2]=pre->sbg[2]-cur->pbg[2];
    }
    /* correction for lever arm */
    if (nla) {
        x[ila+0]=pre->slever[0]-cur->plever[0];
        x[ila+1]=pre->slever[1]-cur->plever[1];
        x[ila+2]=pre->slever[2]-cur->plever[2];
    }
    for (i=0;i<cur->nx*cur->nx;i++) P[i]=pre->Ps[i]-cur->Pp[i];
    matcpy(cur->sba,cur->cba,3,1);
    matcpy(cur->sbg,cur->cbg,3,1);
    matcpy(cur->sMa,cur->sMa,3,3);
    matcpy(cur->sMg,cur->sMg,3,3);

    matcpy(cur->slever,cur->slever,3,1);
}
/* updates ins state---------------------------------------------------------*/
static void upd_ins_state(ins_sol_t *cur,insstate_t *ins)
{
    matcpy(ins->lever,cur->slever,1,3);
    matcpy(ins->P,cur->Ps,cur->nx,cur->nx);

    matcpy(ins->Cbe,cur->sCbe,3,3);
    matcpy(ins->re ,cur->sre ,3,1);
    matcpy(ins->ve ,cur->sve ,3,1);
    matcpy(ins->ae ,cur->sae ,3,1);
    matcpy(ins->ba ,cur->sba ,3,1);
    matcpy(ins->bg ,cur->sbg ,3,1);
    matcpy(ins->Ma ,cur->sMa ,3,3);
    matcpy(ins->Mg ,cur->sMg ,3,3);
}
/* validation of combined solutions ------------------------------------------*/
static int valsmth(const double *Ps,const double *Pc,const double *dx,int nx)
{
    double var[3],factor=10.0;
    int i,flag=0;

    /* 1.position */
    for (i=6;i<9;i++) var[i-6]=Ps[i+i*nx]+Pc[i+i*nx];
    for (i=6;i<9;i++) {
        if (SQR(dx[i])<=SQR(factor)*var[i-6]) continue;  /* ok if in 4-sigma */
        flag=3; /* position error too large */
    }
    /* 2.velocity */
    for (i=3;i<6;i++) var[i-3]=Ps[i+i*nx]+Pc[i+i*nx];
    for (i=3;i<6;i++) {
        if (SQR(dx[i])<=SQR(factor)*var[i-3]) continue;
        flag++;
    }
    /* 3.attitude */
    for (i=0;i<3;i++) var[i]=Ps[i+i*nx]+Pc[i+i*nx];
    for (i=0;i<3;i++) {
        if (SQR(dx[i])<=SQR(factor)*var[i]) continue;
        flag++;
    }
    return flag<3;
}
/* update ins solution through RTS ------------------------------------------*/
static int upd_ins_rts(const ins_sol_t *pre,ins_sol_t *cur,insstate_t *ins,
                       const insopt_t *opt)
{
    double *dx,*Ak,*dP,*Pk_1,*xs,*Ps,factor=0.9999;
    int nx=cur->nx,i;

    dP=zeros(nx,nx); Ak=zeros(nx,nx); Pk_1=zeros(nx,nx);
    dx=zeros(nx,1); xs=zeros(nx,1); Ps=zeros(nx,nx);

    err_corr(pre,cur,opt,dx,dP);

    matcpy(Pk_1,cur->Pp,nx,nx);
    if (!matinv(Pk_1,nx)) {
        matmul33("NTN",cur->Pc,cur->F,Pk_1,nx,nx,nx,nx,Ak);
    }
    else {
        free(dx); free(dP); free(Ak); free(Pk_1);
        free(xs); 
        return 0;
    }
    matmul("NN",nx,1,nx,factor,Ak,dx,0.0,xs);
    matmul33("NNT",Ak,dP,Ak,nx,nx,nx,nx,Pk_1);

    for (i=0;i<nx*nx;i++) Ps[i]=Pk_1[i]+cur->Pc[i];
    matcpy(cur->Ps,Ps,nx,nx);

    if (valsmth(Ps,cur->Pc,xs,nx)) {
        matcpy(cur->Ps,Ps,nx,nx);
    }
    else {
        fprintf(stderr,"%s: validation of backward solutions fail\n",time_str(cur->time,4));

        matcpy(cur->Ps,cur->Pc,nx,nx);
        setzero(xs,1,nx);
    }
    lcclp(xs,cur->cCbe,cur->cre,cur->cve,NULL,NULL,NULL,
          cur->sre,cur->sve,cur->sae,cur->sba,cur->sbg,
          cur->sMa,cur->sMg,cur->slever,cur->sCbe,
          NULL,NULL,opt);
    
    upd_ins_state(cur,ins);

    ins->time =cur->time;
    ins->ns   =cur->ns;
    ins->gstat=cur->gstat;
    ins->stat =INSS_RTS;

    free(dx); free(dP); free(Ak); free(Pk_1);
    free(xs); free(Ps);
    return 1;
}
/* initial first epoch of backward smoother----------------------------------*/
static void init_bcksmh(ins_sol_t *cur)
{
    matcpy(cur->sCbe,cur->cCbe,3,3);
    matcpy(cur->sre ,cur->cre ,3,1);
    matcpy(cur->sve ,cur->cve ,3,1);
    matcpy(cur->sae ,cur->cae ,3,3);
    matcpy(cur->sba ,cur->cba ,3,1);
    matcpy(cur->sbg ,cur->cbg ,3,1);
    matcpy(cur->sMa ,cur->cMa ,3,3);
    matcpy(cur->sMg ,cur->cMg ,3,3);
    matcpy(cur->Ps,cur->Pc,cur->nx,cur->nx);
    matcpy(cur->slever,cur->clever,3,1);
}
/* get ins state from forwad solution binary file----------------------------*/
static int get_ins_state(int n,int bsize,ins_sol_t *ins,const ins_sol_t *pins)
{
    int i,gsize=0;

    fseek(fp_fwd_sol,(long)bsize*(n-1),SEEK_SET);

    gsize+=fread(ins,sizeof(ins_sol_t),1,fp_fwd_sol);

    ins->Ps=pins->Ps;
    ins->Pp=pins->Pp;
    ins->Pc=pins->Pc;
    ins->F =pins->F ;

    for (i=0;i<ins->nx*ins->nx;i++) {
        gsize+=fread(&ins->Pc[i],sizeof(double),1,fp_fwd_sol);
        gsize+=fread(&ins->Pp[i],sizeof(double),1,fp_fwd_sol);
        gsize+=fread(&ins->F [i],sizeof(double),1,fp_fwd_sol);
    }
    return gsize==(1+ins->nx*ins->nx*3);
}
/* backward smoother solution------------------------------------------------
 * args:    rtk_t *rtk       I  rtk data struct
 *          prcopt_t *opt    I  process option
 *          solopt_t *solopt I  solution output option
 * return: numbers of epochs
 * --------------------------------------------------------------------------*/
static int bwdsmh(rtk_t *rtk,const prcopt_t *popt,const solopt_t *solopt)
{
    ins_sol_t *cur,*pre;
    ins_sol_t fins[2],pfins[2]={0};
    int i,j=0,bsize;

#if FORWARD_IN_MEMO
    cur=&insbuf.data[insbuf.n-1];
#else
    bsize=sizeof(ins_sol_t)+SQR(rtk->ins.nx)*sizeof(double)*3;

    if (!(fp_fwd_sol=fopen(solfile,"rb"))) return 0;

    init_insol(&fins[0],rtk->ins.nx);
    init_insol(&fins[1],rtk->ins.nx);

    pfins[0]=fins[0];
    pfins[1]=fins[1];

    if (!get_ins_state(insbuf.n,bsize,&fins[0],&pfins[0])) {
        trace(2,"read ins state fail\n");
        free_insol(&fins[0]);
        free_insol(&fins[1]); return 0;
    }
    cur=&fins[0];
#endif
    /* first smoother epoch */
    init_bcksmh(cur);

    for (i=insbuf.n-2;i>=0;i--) {
#if FORWARD_IN_MEMO
        pre=cur;
        cur=&insbuf.data[i];
#else
        pre=&fins[j%2]; j++;
        if (!get_ins_state(i+1,bsize,&fins[j%2],&pfins[j%2])) continue;
        cur=&fins[j%2]; 
#endif
        if (!upd_ins_rts(pre,cur,&rtk->ins,&popt->insopt)) {
            fprintf(stderr,"rts fail,time=%s\n",time_str(cur->time,4));
        }
        /* solution. */
        ins2sol(&rtk->ins,&rtk->opt.insopt,&rtk->sol);

#if OUT_MONITOR
        /* write solution */
        if (!wrt_solution(rtk,solopt,1)) {
            continue;
        }
#endif
    }
    free_insol(&fins[0]);
    free_insol(&fins[1]);
    return cur==insbuf.data;
}
/* set the temporary path saved by the forward solution file-----------------*/
extern void set_fwd_soltmp_file(const char *file)
{
#if 1
    /* just for debug */
    strcpy(solfile,"/media/sujinglan/Files/fwd_sol.tmp");
#else
    if (file==NULL) {
        strcpy(solfile,"./fwd_sol.tmp");
    }
    else {
        strcpy(solfile,file);
    }
#endif
}
/* rts smoother for ins/gnss loosely coupled---------------------------------
 * args:  imu_t *imup       I  imu measurement data
 *        gsof_data_t *posp I  position measurement data
 *        prcopt_t *popt    I  options
 *        solopt_t *solopt  I  solution options
 *        int port          I  port for solution monitor (0: none)
 *        char *file        I  file path of output solution (NULL: no output)
 * return: number of epochs
 * --------------------------------------------------------------------------*/
extern int lcrts(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                 const solopt_t *solopt,int port,const char *file)
{
    rtk_t rtk={0};
    int n;

    trace(3,"lcrts: port=%d  file=%s\n",port,file);

    /* open result file */
    if (file&&!open_solfile(file)) {
        trace(2,"open file fail\n");
        return 0;
    }
    /* open monitor */
    if (port&&!openmoni(port)) {
        trace(2,"open monitor fail\n");
        return 0;
    }
    /* forward filter */
    if (fwdfilt(imu,pos,popt,solopt,&rtk)) bwdsmh(&rtk,popt,solopt);
    if (!(n=insbuf.n)) {
        trace(2,"rts solution fail\n");
    }
    free_fbsol(&insbuf);
    rtkfree(&rtk);

    /* close monitor/file */
    if (fp_fwd_sol) fclose(fp_fwd_sol);
    if (port) close_moni(&moni);
    if (file) strclose(&frst);
    return n;
}