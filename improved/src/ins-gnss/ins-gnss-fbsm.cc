/*-----------------------------------------------------------------------------
* ins-gnss-fbsm.cc : ins-gnss combine forward and backward filters by
*                    fixed-point smoother
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
* history : 2018/09/25 1.0 new
*----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants-----------------------------------------------------------------*/
#define MAXTIMEDIFF  3.0                 /* max time difference for RTS smoother */
#define NP           5                   /* number of position solutions */
#define MINVEL       5.0                 /* min velocity for initial ins states */
#define MAXGYRO      (30.0*D2R)          /* max rotation speed value for initial */
#define MAXDIFF      10.0                /* max time difference between solution */
#define INTKEEPALIVE 1000                /* keep alive interval (ms) */
#define OUTSOLFRQ    100                 /* frequency of output ins solutions */
#define OUT_MONITOR  1                   /* output solution to monitor */
#define SOL_OUTPUT_FILE 0                /* output solution to file */
#define DEGRADE_FACTOR  0.618            /* degrade factor to correction if combined solution is inconsistent */

typedef struct {                         /* forward ins solution data type */
    gtime_t time;                        /* ins solution time */
    int nx;                              /* number of error states */
    int ns,stat,gstat;                   /* number of valid satellite/ins solution status/gps status */

    double *x,*P;                        /* forward ins solution data */
} ins_fsol_t;

/* constants/global variables -----------------------------------------------*/
static stream_t    moni={0};             /* monitor stream */
static stream_t    frst={0};             /* solution result file */
static prcopt_t prcopt ={0};             /* processing options */
static solopt_t solopt ={0};             /* solution options */
static filopt_t filopt ={""};            /* file options */

static int ipos=0;                       /* current gsof message index */
static int iimu=0;                       /* current imu measurement data */
static int keepalive=0;                  /* keep alive flag */
static int timeout  =10000;              /* timeout time (ms) */
static int reconnect=10000;              /* reconnect interval (ms) */
static int week=0;                       /* GPS week */
static char solfile[1024];               /* solution output file path */
static FILE *fp_fwd_sol=NULL;            /* foeward solution file pointer */

/* initial ins forward solution ---------------------------------------------*/
static void init_fsol(ins_fsol_t *fsol,const insopt_t *opt)
{
    gtime_t t0={0};
    fsol->ns=fsol->gstat=fsol->stat=0;
    fsol->time=t0; 
    fsol->nx=xnX(opt);

    fsol->x=mat(fsol->nx,       1);
    fsol->P=mat(fsol->nx,fsol->nx);
}
/* free ins forward solution-------------------------------------------------*/
static void free_fsol(ins_fsol_t *fsol)
{
    if (fsol->x) free(fsol->x); fsol->x=NULL;
    if (fsol->P) free(fsol->P); fsol->P=NULL;
}
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
                    imud_t *imuz,int ws,int type)
{
    int i,k;
    if (0<=iimu&&iimu<imu->n) {
        trace(3,"input imu measurement data\n");
    }
    else return 0;

    if (type==0) { /* forward */
        for (k=0,i=(iimu-ws<0?0:iimu-ws);k<ws&&i>=0&&i<imu->n;k++,i++) {
            /* adjust time */
            adj_imutime(&imu->data[i],opt);
            imuz[k]=imu->data[i];
        }
        adj_imutime(&imu->data[iimu],opt);
        *imudata=imu->data[iimu++];
    }
    else if (type==1) { /* backward */
        for (k=ws-1,i=(iimu+ws>=imu->n?imu->n-1:iimu+ws);
             k>=0&&i>=0&&i<imu->n;k--,i--) {
            imuz[k]=imu->data[i];
        }
        *imudata=imu->data[iimu--];
    }
    return 1;
}
/* input position measurement data-------------------------------------------*/
static int inputpos(const gsof_data_t *pos,gsof_t *gsofdata,gtime_t timu,
                    const prcopt_t* opt,int type)
{
    int flag=0,i;
    static double ts=0.5/opt->insopt.hz,sowi,sowp;

    trace(3,"inputpos:\n");

    if (0<=ipos&&ipos<pos->n) {
        trace(3,"input position measurement data\n");
    }
    else return 0;

    sowi=time2gpst(timu,NULL);

    if (type==0) { /* forward */
        for (i=ipos>5?ipos-5:0;i<pos->n&&i>=0;i++) {
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
    }
    else if (type==1) {
        for (i=ipos+5>=pos->n?pos->n-1:ipos+5;i>=0&&i<pos->n;i--) {
            if (fabs(sowi-time2gpst(pos->data[i].t,NULL))<ts) {
                flag=1;
                ipos=i;
                *gsofdata=pos->data[i];

                /* GPS week */
                time2gpst(gsofdata->t,&week);
                break;
            }
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
#if 1   /* for debugs */
        fprintf(stderr,"%s",buff);
#endif
    }
    return n;
}
/* update solution status----------------------------------------------------*/
static void update_stat(const gmea_t *gm,insstate_t *ins)
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
/* add ins solutions to buffer-----------------------------------------------*/
static int add_ins_sol(const insstate_t *ins)
{
    double omg[3];
    int i;

    /* write solution to file in binary */
    fwrite(&ins->time,sizeof(gtime_t),1,fp_fwd_sol);
    fwrite(&ins->ns  ,sizeof(int),1,fp_fwd_sol);
    fwrite(&ins->nx  ,sizeof(int),1,fp_fwd_sol);

    fwrite(&ins->stat ,sizeof(int),1,fp_fwd_sol);
    fwrite(&ins->gstat,sizeof(int),1,fp_fwd_sol);

    so3_log(ins->Cbe,omg,NULL);
    fwrite(omg,sizeof(double),3,fp_fwd_sol);

    fwrite(ins->ve,sizeof(double),3,fp_fwd_sol);
    fwrite(ins->re,sizeof(double),3,fp_fwd_sol);
    fwrite(ins->ba,sizeof(double),3,fp_fwd_sol);
    fwrite(ins->bg,sizeof(double),3,fp_fwd_sol);

    for (i=0;i<ins->nx*ins->nx;i++) {
        fwrite(&ins->P[i],sizeof(double),1,fp_fwd_sol);
    }
    return 1;
}
/* get forward ins solution data from file---------------------------------- */
static int get_fwd_sol(ins_fsol_t *sol,int i,int bsize)
{
    fread(&sol->time,sizeof(gtime_t),1,fp_fwd_sol);
    fread(&sol->ns,sizeof(int),1,fp_fwd_sol);
    fread(&sol->nx,sizeof(int),1,fp_fwd_sol);

    fread(&sol->stat ,sizeof(int),1,fp_fwd_sol);
    fread(&sol->gstat,sizeof(int),1,fp_fwd_sol);

    fread(sol->x+ 0,sizeof(double),3,fp_fwd_sol);
    fread(sol->x+ 3,sizeof(double),3,fp_fwd_sol);
    fread(sol->x+ 6,sizeof(double),3,fp_fwd_sol);
    fread(sol->x+ 9,sizeof(double),3,fp_fwd_sol);
    fread(sol->x+12,sizeof(double),3,fp_fwd_sol);
    for (i=0;i<sol->nx*sol->nx;i++) {
        fread(&sol->P[i],sizeof(double),1,fp_fwd_sol);
    }
    fseek(fp_fwd_sol,-bsize,SEEK_CUR);
}
/* get solution time from forward solution file -----------------------------*/
static int get_time_fsol(FILE *fp_sol,gtime_t *time)
{
    gtime_t t0={0};
    if (!fread(time,sizeof(gtime_t),1,fp_sol)) {
        *time=t0;
        return 0;
    }
    fseek(fp_sol,-sizeof(gtime_t),SEEK_CUR);
    return 1;
}
/* search forward ins solution by time---------------------------------------*/
static int search_fsol(int bsize,gtime_t time,ins_fsol_t *sol)
{
    gtime_t ts,t0={0};
    int i=0;

    sol->time=t0;
    while (true) {
        if      (i== 1&& feof (fp_fwd_sol)) break;
        else if (i==-1&&!ftell(fp_fwd_sol)) break;

        fseek(fp_fwd_sol,(i==0?-1:i)*bsize,SEEK_CUR);

        if (!get_time_fsol(fp_fwd_sol,&ts)) continue;
        if (fabs(timediff(ts,time)<=DTTOL)) {
            get_fwd_sol(sol,i,bsize);
            break;
        }
        if (timediff(ts,time)>0.0) i= 1;
        if (timediff(ts,time)<0.0) i=-1;
    }
    return sol->time.time>0;
}
/* backward initial----------------------------------------------------------*/
static int bwdinit(const imu_t *imu,const gsof_data_t *pos,rtk_t *rtk)
{
    insstate_t *ins=&rtk->ins;
    int i;
    for (i=imu->n-1;i>=0;i--) {
        if (fabs(timediff(rtk->ins.time,imu->data[i].time))<DTTOL) break;
    }
    iimu=i<=0?0:i;
    ipos=ipos+10>=pos->n?pos->n-1:ipos+10;

    setzero(ins->x ,ins->nx,1);
    setzero(ins->xa,ins->nx,1);
#if 0
    getP0(&rtk->opt.insopt,ins->P );
    getP0(&rtk->opt.insopt,ins->Pa);
#endif
    rtk->opt.insopt.soltype=1;
    return 1;
}
/* validation of combined solutions ------------------------------------------*/
static int valsmth(const double *Pb,const double *Pf,const double *dx,int nx,
                   const insopt_t *opt)
{
    double var[3],factor=10.0;
    int i,flag=0;

    /* 1.position */
    for (i=6;i<9;i++) var[i-6]=Pb[i+i*nx]+Pf[i+i*nx];
    for (i=6;i<9;i++) {
        if (SQR(dx[i])<=SQR(factor)*var[i-6]) continue;  /* ok if in 4-sigma */
        flag=3; /* position error too large */
    }
    /* 2.velocity */
    for (i=3;i<6;i++) var[i-3]=Pb[i+i*nx]+Pf[i+i*nx];
    for (i=3;i<6;i++) {
        if (SQR(dx[i])<=SQR(factor)*var[i-3]) continue;
        flag++;
    }
    /* 3.attitude */
    for (i=0;i<3;i++) var[i]=Pb[i+i*nx]+Pf[i+i*nx];
    for (i=0;i<3;i++) {
        if (SQR(dx[i])<=SQR(factor)*var[i]) continue;
        flag++;
    }
    /* accl. bias */
    if (xnBa(opt)) {
        for (i=0;i<3;i++) var[i]=Pb[i+xiBa(opt)+(i+xiBa(opt))*nx]+Pf[i+xiBa(opt)+(i+xiBa(opt))*nx];
        for (i=0;i<3;i++) {
            if (SQR(dx[i+xiBa(opt)])<=SQR(factor)*var[i]) continue;
            flag++;
        }
    }
    /* gyro. bias */
    if (xnBg(opt)) {
        for (i=0;i<3;i++) var[i]=Pb[i+xiBg(opt)+(i+xiBg(opt))*nx]+Pf[i+xiBg(opt)+(i+xiBg(opt))*nx];
        for (i=0;i<3;i++) {
            if (SQR(dx[i+xiBg(opt)])<=SQR(factor)*var[i]) continue;
            flag++;
        }
    }
    return flag<3;
}
/* copy ins stats------------------------------------------------------------*/
static void copy_ins_state(const insstate_t *pins,insstate_t *inss)
{
    insstate_t inst=*inss;
    *inss=*pins;

    inss->x =inst.x ; inss->P =inst.P ;
    inss->xa=inst.xa; inss->Pa=inst.Pa;
    inss->xb=inst.xb; inss->Pb=inst.Pb;
    inss->F =inst.F ; inss->P0=inst.P0;

    inss->rtkp=inst.rtkp;
    inss->gmeas.data=inst.gmeas.data;
    inss->gmeas.n   =inst.gmeas.n;
    inss->gmeas.nmax=inst.gmeas.nmax;

    matcpy(inss->F,pins->F,inss->nx,inss->nx);
    matcpy(inss->P,pins->P,inss->nx,inss->nx);

    matcpy(inss->P0,pins->P0,inss->nx,inss->nx);
    matcpy(inss->Pa,pins->Pa,inss->nx,inss->nx);
    matcpy(inss->Pb,pins->Pb,inss->nx,inss->nx);

    matcpy(inss->x ,pins->x ,inss->nx,1);
    matcpy(inss->xa,pins->xa,inss->nx,1);
    matcpy(inss->xb,pins->xb,inss->nx,1);
}
/* combine forward/backward solutions  --------------------------------------*/
static int combres(insstate_t *ins,const insopt_t *opt,insstate_t *inss)
{
    static int bszie=sizeof(gtime_t)+4*sizeof(int)+(15+ins->nx*ins->nx)*sizeof(double);
    ins_fsol_t fsol={0};
    double *dx,omg[3],fCbe[9],dCbe[9],*Ps,*Pf,*Pb,*dxs,factor=0.9999;
    int i,nx=ins->nx;

    copy_ins_state(ins,inss);
#if 0
    /* only output backward filter solution if enable */
    return 1;
#endif
    init_fsol(&fsol,opt);
    if (!search_fsol(bszie,ins->time,&fsol)) {
        fprintf(stderr,"%s: combine forward/backward fail\n",time_str(ins->time,4));
        free_fsol(&fsol);
        return 0;
    }
    dx=mat(nx, 1); Ps=mat(nx,nx);
    Pf=mat(nx,nx); Pb=mat(nx,nx);

    if (xnBa(opt)) {
        for (i=0;i<3;i++) dx[xiBa(opt)+i]=ins->ba[i]-fsol.x[xiBa(opt)+i];
    }
    if (xnBg(opt)) {
        for (i=0;i<3;i++) dx[xiBg(opt)+i]=ins->bg[i]-fsol.x[xiBg(opt)+i];
    }
    so3_exp(fsol.x,fCbe);
    matmul("NT",3,3,3,1.0,ins->Cbe,fCbe,0.0,dCbe);
    so3_log(dCbe,omg,NULL);
    
    for (i=0;i<3;i++) dx[i]=omg[i];
    for (i=3;i<6;i++) dx[i]=ins->ve[i-3]-fsol.x[i];
    for (i=6;i<9;i++) dx[i]=ins->re[i-6]-fsol.x[i];

    if (!valsmth(ins->P,fsol.P,dx,ins->nx,opt)) {
        fprintf(stderr,"%s: combine forward/backward degrade\n",time_str(ins->time,4));
        factor=DEGRADE_FACTOR;
       
        inss->stat=INSS_DEGRADE;
    }
    matcpy(Pf,fsol.P,nx,nx);
    matcpy(Pb,ins->P,nx,nx);
    if (!matinv(Pf,nx)&&!matinv(Pb,nx)) {
        for (i=0;i<nx*nx;i++) Ps[i]=Pf[i]+Pb[i];
    }
    else {
        free(dx); free(Ps); free(Pf); free(Pb);
        free_fsol(&fsol);
        return 0;
    }
    if (matinv(Ps,nx)) {
        free(dx); free(Ps); free(Pf); free(Pb);
        free_fsol(&fsol);
        return 0;
    }
    dxs=mat(nx,1);
    matmul33("NNN",Ps,Pb,dx,nx,nx,nx,1,dxs);
    for (i=0;i<nx;i++) dxs[i]*=-factor;

    tracemat_std(3,dxs,1,nx,10,6);

    if (xnBa(opt)) matcpy(inss->ba   ,&fsol.x[xiBa(opt)],3,1);
    if (xnBg(opt)) matcpy(inss->bg   ,&fsol.x[xiBg(opt)],3,1);
    if (xnLa(opt)) matcpy(inss->lever,&fsol.x[xiLa(opt)],3,1);

    lcclp(dxs,fCbe,fsol.x+6,fsol.x+3,ins->fb,ins->omgb,ins->Gg,
          inss->re,inss->ve,inss->ae,
          inss->ba,inss->bg,inss->Ma,inss->Mg,
          inss->lever,inss->Cbe,inss->fb,inss->omgb,opt);

    matcpy(inss->P,Ps,nx,nx);
    if (ins->stat!=INSS_DEGRADE) inss->stat=INSS_FBCOMB;

    free(dx); free(Ps); free(Pf);
    free(Pb); free(dxs);
    free_fsol(&fsol);
    return 1;
}
/* set the temporary path saved by the forward solution file-----------------*/
extern void set_fwdtmp_file(const char *file)
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
/* open forward solution binary file-----------------------------------------*/
static int open_fwdsol()
{
    fp_fwd_sol=fopen(solfile,"rb"); if (fp_fwd_sol) fseek(fp_fwd_sol,0,SEEK_END);
    return fp_fwd_sol==NULL;
}
/* loosely coupled filter of smoother----------------------------------------
 * args:    imu_t *imu        I  imu measurement data
 *          gsof_data_t *pos  I  position measurement data
 *          insopt_t *opt     I  ins options
 *          int type          I  process direction (forward/backward)
 * return:  number of epochs
 *---------------------------------------------------------------------------*/
static int fbfilt(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                  const solopt_t *solopt,rtk_t *rtk,int type)
{
    imud_t imus={0},*imuz=NULL;
    gsof_t poss={0};
    gmea_t gmea={0};
    const insopt_t *iopt=&rtk->opt.insopt;
    rtk_t rtks={0};
    int ws,init=0,n=0;

    trace(3,"fwdfilt: ni=%d  np=%d\n",imu->n,pos->n);

    /* static detect window size */
    ws=iopt->zvopt.ws<=0?5:iopt->zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* initialization */
    if (type==0) rtkinit(rtk,popt);
    if (type==1) {
        rtkinit(&rtks,popt); bwdinit(imu,pos,rtk); init=1;
        if (open_fwdsol()) {
            goto exit;
        }
    }
    /* start loosely coupled */
    while (inputimu(imu,&imus,popt,imuz,ws,type)) {
        if (inputpos(pos,&poss,imus.time,popt,type)) {

            if (!gsof2gnss(&poss,&gmea)||outagegsof(popt,&poss)) continue;
            if (!init) {
                /* initialization */
                if (!(init=init_ins(&imus,&poss,iopt,&rtk->ins))) continue;
            }
            else {
                /* coupled. */
                lcigpos(iopt,&imus,&rtk->ins,&gmea,INSUPD_MEAS);
            }
        }
        else if (init==1) {
            /* ins mech. */
            lcigpos(iopt,&imus,&rtk->ins,NULL,INSUPD_TIME);
        }
        else continue;

        /* motion constraint update */
        motion(iopt,imuz,&rtk->ins,&imus,ws);

        /* odometry velocity aid */
        if (iopt->odo) {
            odo(iopt,&imus,&imus.odo,&rtk->ins);
        }
        /* sol. status */
        update_stat(&gmea,&rtk->ins);

        /* forward/backward combined solution */
        if (type&&!combres(&rtk->ins,iopt,&rtks.ins)) {
            trace(2,"%s: forward/backward combined solution fail\n",time_str(rtk->ins.time,4));
            continue;
        }
        /* solution. */
        if (type==0) ins2sol(&rtk->ins,iopt,&rtk->sol);
        else {
            ins2sol(&rtks.ins,iopt,&rtks.sol);
        }
#if OUT_MONITOR
        /* write solution */
        if (!wrt_solution(type?&rtks:rtk,solopt,type)) continue;
#endif
        /* ins solution add. in forward */
        if (type==0&&!add_ins_sol(&rtk->ins)) {
            trace(2,"add ins solution fail\n");
            continue;
        }
        n++;
    }
exit:
    fclose(fp_fwd_sol); fp_fwd_sol=NULL;
    rtkfree(&rtks);
    free(imuz);
    return n;
}
/* forward/backward smoother for ins/gnss loosely coupled--------------------
 * args:  imu_t *imup       I  imu measurement data
 *        gsof_data_t *posp I  position measurement data
 *        prcopt_t *popt    I  options
 *        solopt_t *solopt  I  solution options
 *        int port          I  port for solution monitor (0: none)
 *        char *file        I  file path of output solution (NULL: no output)
 * return: number of epochs
 * --------------------------------------------------------------------------*/
extern int lcfbsm(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                  const solopt_t *solopt,int port,const char *file)
{
    rtk_t rtk={0};
    int n=0;

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
    /* forward/backward solution */
    if (fbfilt(imu,pos,popt,solopt,&rtk,0)) n=fbfilt(imu,pos,popt,solopt,&rtk,1);
    if (n==0) {
        trace(2,"forward/backward combined solution fail\n");
    }
    rtkfree(&rtk);

    /* close monitor/file */
    if (fp_fwd_sol) fclose(fp_fwd_sol);
    if (port) close_moni(&moni);
    if (file) strclose(&frst);
    return 0;
}


