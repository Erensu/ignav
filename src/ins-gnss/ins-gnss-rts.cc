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
#include <navlib.h>

/* constants-----------------------------------------------------------------*/
#define NP           3              /* number of position solutions */
#define MINVEL       5.0            /* min velocity for initial ins states */
#define MAXGYRO      (30.0*D2R)     /* max rotation speed value for initial */
#define MAXDIFF      10.0           /* max time difference between solution */
#define INTKEEPALIVE 1000           /* keep alive interval (ms) */
#define OUTSOLFRQ    100            /* frequency of output ins solutions */
#define OUT_MONITOR  1              /* ouput solution to monitor */

/* type struct---------------------------------------------------------------*/
typedef struct {                    /* rts forward filter result type */
    gtime_t time;                   /* time of update ekf states */
    double *x,*P;                   /* updated ekf states and covariance */
    double *x0,*P0,*F;              /* predict ekf states and covariance/transformation matrix */
} fwd_t;

typedef struct {                    /* rts backward solution type */
    gtime_t time;                   /* time of backward solutions */
    double *Ks;                     /* backward smoother gain matrix*/
    double *xs,*Ps;                 /* backward smoother solutions */
} bwd_t;

typedef struct {                    /* forward solutions type */
    fwd_t *data;                    /* forward solutions data */
    int n,nmax;                     /* number and max number of solutions */
} fwdd_t;

typedef struct {                    /* backward solutions type */
    bwd_t *data;                    /* backward solutions data */
    int n,nmax;                     /* number and max numbers of solutions */
} bwdd_t;

typedef struct {                    /* ins solution buffer type */
    int n,nmax;                     /* number and max number of solutions */
    insstate_t  *data;              /* ins solution data */
} ins_solbuf_t;

/* constants/global variables -----------------------------------------------*/
static imu_t       imu ={0};        /* imu measurement data */
static gsof_data_t pos ={0};        /* position measurement data for ins/gnss coupled */
static fwdd_t      fwd ={0};        /* forward solutions */
static bwdd_t      bwd ={0};        /* backward solutions */
static stream_t    moni={0};        /* monitor stream */
static prcopt_t prcopt={0 };        /* processing options */
static solopt_t solopt={0 };        /* solution options */
static filopt_t filopt={""};        /* file options */
static ins_solbuf_t insbuf={0};     /* ins solution buffer */

static int ipos=0;                  /* current gsof message index */
static int iimu=0;                  /* current imu measurement data */
static int keepalive=0;             /* keep alive flag */
static int timeout  =10000;         /* timeout time (ms) */
static int reconnect=10000;         /* reconnect interval (ms) */

/* input imu measurement data------------------------------------------------*/
static int inputimu(imud_t *imudata,const prcopt_t* opt,imud_t *imuz,int ws)
{
    int i,k;

    if (0<=iimu&&iimu<imu.n) {
        trace(3,"input imu measurement data\n");
    }
    else return 0; 

    /* imu data for static detect */
    for (k=0,i=iimu;k<ws&&i>=0&&i<imu.n;
         k++,i++) {
        imuz[k]=imu.data[i];
    }
    /* forward input */
    *imudata=imu.data[iimu++];
    return 1;
}
/* input position measurement data-------------------------------------------*/
static int inputpos(gsof_t *gsofdata,gtime_t timu,const prcopt_t* opt)
{
    int flag=0,i;
    static double ts=0.5/opt->insopt.hz;

    trace(3,"inputpos:\n");

    if (0<=ipos&&ipos<pos.n) {
        trace(3,"input position measurement data\n");
    }
    else return 0;

    /* find the closest gsof message for imu data */
    for (i=ipos>5?ipos-5:0;i<pos.n;i++) {
        if (fabs(timediff(timu,pos.data[i].t))<ts) {
            flag=1;
            ipos=i;
            *gsofdata=pos.data[i];
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
        trace(2,"invalid position measurement\n");
        return 0;
    }
    /* save pvt solution buffer */
    for (i=0;i<NP-1;i++) sols[i]=sols[i+1]; sols[i]=*pos;
    for (i=0;i<NP;i++) {
        if (sols[i].solq>iopt->iisu||sols[i].solq==SOLQ_NONE) return 0;
    }
    /* compute velocity from solutions */
    sol2vel(sols+MAXSOL-1,sols+MAXSOL-2,vr);

    /* check velocity ok? */
    if (norm(vr,3)<MINVEL||norm(imu->gyro,3)>MAXGYRO) {
        return 0;
    }
    /* initial ins states */
    if (!ant2inins(sols[MAXSOL-1].t,sols[MAXSOL-1].pos,vr,
                   iopt,
                   NULL,
                   ins,NULL)) {
        trace(2,"initial ins state fail\n");
        return 0;
    }
    ins->time=sols[MAXSOL-1].t;

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
/* write solution status output stream ---------------------------------------*/
static int wrt_solution(const rtk_t *rtk,const solopt_t *solopt)
{
    unsigned char buff[1024];
    static int c=0;
    int n=0;

    if (moni.port) {
        n=outsols(buff,&rtk->sol,rtk->rb,solopt,
                  rtk->opt.mode>=PMODE_INS_UPDATE?rtk->ins:NULL,
                  &rtk->opt.insopt,1);
        if (rtk->opt.mode>=PMODE_INS_UPDATE&&rtk->opt.mode<=PMODE_INS_TGNSS) {
            if (c++>OUTSOLFRQ) {
                strwrite(&moni,buff,n); c=0;
            }
        }
        else {
            strwrite(&moni,buff,n);
        }
    }
    return n;
}
/* free forward/backward solutions-------------------------------------------*/
static int free_sol(fwdd_t *fwd,bwdd_t *bwd)
{
    int i; if (fwd) {
        for (i=0;i<fwd->n;i++) {
            free(fwd->data[i].x ); free(fwd->data[i].P);
            free(fwd->data[i].x0); free(fwd->data[i].P0);
        }
        free(fwd->data);
        fwd->n=0;
        fwd->nmax=0;
    }
    if (bwd) {
        for (i=0;i<bwd->n;i++) {
            free(bwd->data[i].xs); free(bwd->data[i].Ps);
            free(bwd->data[i].Ks);
        }
        free(bwd->data);
        bwd->n=0;
        bwd->nmax=0;
    }
}
/* add forward solution to buffer--------------------------------------------*/
static int add_fwd_sol(fwdd_t *fwd,const fwd_t *data)
{
    fwd_t *pfwd;

    if (fwd->nmax<=fwd->n) {
        fwd->nmax+=1024;
        if (!(pfwd=(fwd_t *)realloc(fwd->data,sizeof(fwd_t)*fwd->nmax))) {
            trace(1,"add_fwd_sol malloc error: n=%d\n",fwd->nmax);
            free_sol(fwd,NULL);
            return 0;
        }
        fwd->data=pfwd;
    }
    fwd->data[fwd->n++]=*data;
    return 1;
}
/* add backward solution to buffer-------------------------------------------*/
static int add_bwd_sol(bwdd_t *bwd,const bwd_t *data)
{
    bwd_t *pbwd;

    if (bwd->nmax<=bwd->n) {
        bwd->nmax+=1024;
        if (!(pbwd=(bwd_t *)realloc(bwd->data,sizeof(bwd_t)*bwd->nmax))) {
            trace(1,"add_bwd_sol malloc error: n=%d\n",bwd->nmax);
            free_sol(NULL,bwd);
            return 0;
        }
        bwd->data=pbwd;
    }
    bwd->data[bwd->n++]=*data;
    return 1;
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
    if (ins->x0) free(ins->x0); ins->x0=NULL;
    if (ins->P0) free(ins->P0); ins->P0=NULL;
    if (ins->F ) free(ins->F ); ins->F =NULL;
    if (ins->gmeas.data) free(ins->gmeas.data); ins->gmeas.data=NULL;

    ins->nx=ins->nb=0;
    ins->gmeas.n=ins->gmeas.nmax=0;
}
/* add ins solutions to buffer-----------------------------------------------*/
static int add_ins_sol(ins_solbuf_t *buf,const insstate_t *data)
{
    int i; insstate_t *pins;
    if (buf->nmax<buf->n) {
        buf->n+=1024;
        if (!(pins=(insstate_t *)realloc(buf->data,sizeof(insstate_t)*buf->nmax))) {

            trace(1,"add_ins_sol malloc error: n=%d\n",buf->nmax);
            for (i=0;i<buf->n;i++) {
                freeins(buf->data+i); free(buf->data+i);
            }
            return 0;
        }
        buf->data=pins;
    }
    pins=&buf->data[buf->n++];
    *pins=*data;

    pins->x =zeros(1,data->nx);
    pins->xa=zeros(1,data->nx);
    pins->xb=zeros(1,data->nx);

    pins->P =zeros(data->nx,data->nx);
    pins->Pa=zeros(data->nx,data->nx);
    pins->Pb=zeros(data->nx,data->nx);
    pins->x0=NULL;
    pins->P0=NULL;
    pins->F =NULL;

    /* copy ins states */
    matcpy(pins->P ,data->P ,data->nx,data->nx);
    matcpy(pins->Pa,data->Pa,data->nx,data->nx);
    matcpy(pins->Pb,data->Pb,data->nx,data->nx);

    matcpy(pins->x ,data->x ,1,data->nx);
    matcpy(pins->xa,data->xa,1,data->nx);
    matcpy(pins->xb,data->xb,1,data->nx);
    return 1;
}
/* initial forward/backward solution-----------------------------------------*/
static void init_fbsol(const insstate_t *ins,fwd_t *fwds,bwd_t *bwds)
{
    int nx=ins->nx;
    gtime_t t0={0};
    if (fwds) {
        fwds->x =zeros(1,nx); fwds->P =zeros(nx,nx);
        fwds->x0=zeros(1,nx); fwds->P0=zeros(nx,nx);
        fwds->time=t0;
    }
    if (bwds) {
        bwds->Ks=zeros(nx,nx);
        bwds->xs=zeros(nx,1 );
        bwds->Ps=zeros(nx,nx);
        bwds->time=t0;
    }
}
/*  save forward solutions---------------------------------------------------*/
static int sve_solution(const rtk_t *rtk,const prcopt_t *opt)
{
    int nx=rtk->ins.nx;
    fwd_t fwds={0};

    init_fbsol(&rtk->ins,&fwds,NULL);

    matcpy(fwds.x,rtk->ins.x, 1,nx);
    matcpy(fwds.P,rtk->ins.P,nx,nx);

    matcpy(fwds.x0,rtk->ins.x0, 1,nx);
    matcpy(fwds.P0,rtk->ins.P0,nx,nx);

    fwds.time=rtk->ins.time;
    return add_fwd_sol(&fwd,&fwds);
}
/* forward filter of rts-----------------------------------------------------
 * args:    imu_t *imu        I  imu measurement data
 *          gsof_data_t *pos  I  position measurement data
 *          insopt_t *opt     I  ins options
 * return:  number of epochs
 *---------------------------------------------------------------------------*/
static int fwdfilt(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                   const solopt_t *solopt)
{
    imud_t imus={0},*imuz=NULL;
    gsof_t poss={0};
    gmea_t gmea={0};
    rtk_t  rtk={{0}};
    insopt_t *iopt=&popt->insopt;
    int ws,init=0;

    trace(3,"fwdfilt: ni=%d  np=%d\n",imu->n,pos->n);

    /* static detect window size */
    ws=popt->insopt.zvopt.ws<=0?5:popt->insopt.zvopt.ws;
    imuz=(imud_t*)malloc(sizeof(imud_t)*ws);

    /* rtk init. */
    rtkinit(&rtk,popt);

    while (inputimu(&imus,popt,imuz,ws)) {
        if (inputpos(&poss,imus.time,popt)) {

            if (!gsof2gnss(&poss,&gmea)) continue;
            if ( outagegsof(popt,&poss)) continue;

            if (!init) {
                /* initialization */
                init=init_ins(&imus,&poss,&popt->insopt,&rtk.ins);
            }
            else {

                /* coupled. */
                lcigpos(&popt->insopt,&imus,&rtk.ins,&gmea,
                        INSUPD_MEAS);
            }
        }
        else if (init==1) { /* ins mech. */
            lcigpos(&popt->insopt,&imus,&rtk.ins,
                    &gmea,INSUPD_TIME);
        }
        else continue;

        /* motion constraint update */
        motion(iopt,imuz,&rtk.ins,&imus,ws);

        /* odometry velocity aid */
        if (iopt->odo) {
            odo(iopt,&imus,&imus.odo,&rtk.ins);
        }
        /* solution. */
        ins2sol(&rtk.ins,&popt->insopt,&rtk.sol);

#if OUT_MONITOR
        /* write solution */
        wrt_solution(&rtk,solopt);
#endif
        /* save solution */
        if (!sve_solution(&rtk,popt)) {
            trace(2,"add fwd. solution fail\n");
        }
        /* ins solution add. */
        if (!add_ins_sol(&insbuf,&rtk.ins)) {
            trace(2,"add ins solution fail\n");
            continue;
        }
    }
    return fwd.n;
}
/* rts-----------------------------------------------------------------------*/
static int rts(const rtk_t *rtk,const fwd_t *fwd,const double *xsk,
               const double *Psk)
{
    int i,nx=rtk->ins.nx;
    double *dx,*dP,*A,*P0,*dxs,*dPs;
    bwd_t bwds={0};

    dx=zeros(1 ,nx);
    dP=zeros(nx,nx);
    A =zeros(nx,nx);
    P0=zeros(nx,nx);

    for (i=0;i<nx   ;i++) dx[i]=xsk[i]-fwd->x0[i];
    for (i=0;i<nx*nx;i++) dP[i]=Psk[i]-fwd->P0[i];

    matcpy(P0,fwd->P0,nx,nx);
    if (!matinv(P0,nx)) {
        matmul33("NTN",fwd->P,fwd->F,P0,nx,nx,nx,nx,A);

        dxs=zeros(1,nx);
        matmul("NN",nx,1,nx,1.0,A,dx,0.0,dxs);
        dPs=zeros(nx,nx);
        matmul33("NNT",A,dP,A,nx,nx,nx,nx,dPs);

        for (i=0;i<nx   ;i++) dx[i]=fwd->x[i]+dxs[i];
        for (i=0;i<nx*nx;i++) dP[i]=fwd->P[i]+dPs[i];
        free(dxs);
        free(dPs);
    }
    else {
        free(dx); free(dP);
        free(A ); free(P0);
        return 0;
    }
    init_fbsol(&rtk->ins,NULL,&bwds);

    matcpy(bwds.Ks,A ,nx,nx);
    matcpy(bwds.xs,dx, 1,nx);
    matcpy(bwds.Ps,dP,nx,nx);

    bwds.time=fwd->time;

    add_bwd_sol(&bwd,&bwds);

    free(dx); free(dP);
    free(A ); free(P0);
    return 1;
}
/* update ins states in backward smoother------------------------------------*/
static void upd_bwdsmh(insstate_t *ins,const double *xs,const double *Ps,
                       const insopt_t *opt)
{
    /* update states/cov. */
    matcpy(ins->P,Ps,ins->nx,ins->nx);
    matcpy(ins->x,xs,ins->nx,1);

    /* close loop */
    clp(ins,opt,xs);
}
/* find ins solution---------------------------------------------------------*/
static int find_ins_sol(const ins_solbuf_t *buf,gtime_t time,double dt)
{
    int index=timediff(time,buf->data[0].time)/dt;

    index=index<0?0:index;

    while (true) {
        if (index<0||index>=buf->n) return -1;

        if      (timediff(buf->data[index].time,time)<0.0) index++;
        else if (timediff(buf->data[index].time,time)>0.0) index--;
        else if (fabs(timediff(buf->data[index].time,time))<DTTOL) {
            return index;
        }
    }
}
/* backward smoother solution------------------------------------------------
 * args:    rtk_t *rtk       I  rtk data struct
 *          prcopt_t *opt    I  process option
 *          solopt_t *solopt I  solution output option
 * return: numbers of epochs
 * --------------------------------------------------------------------------*/
static int bwdsmh(const rtk_t *rtk,const prcopt_t *popt,const solopt_t *solopt)
{
    double *xs=NULL,*Ps=NULL,dt=1.0/popt->insopt.nhz;
    int i,nx=rtk->ins.nx,j;
    bwd_t bwds={0};
    rtk_t rtk0={0};

    trace(3,"bwdsmh:\n");

    init_fbsol(&rtk->ins,NULL,&bwds);

    matcpy(bwds.xs,fwd.data[fwd.n-1].x,nx,nx);
    matcpy(bwds.Ps,fwd.data[fwd.n-1].P,nx,nx);

    /* first backward solution */
    xs=bwds.xs;
    Ps=bwds.Ps;

    /* initial. */
    rtk0=*rtk;

    for (i=fwd.n-2;i>=0;i++) {
        if (!rts(rtk,fwd.data+i,xs,Ps)) continue;

        /* next backward solution */
        xs=bwd.data[bwd.n-1].xs;
        Ps=bwd.data[bwd.n-1].Ps;

        /* find ins solution */
        if ((j=find_ins_sol(&insbuf,fwd.data[i].time,dt))<0) continue;

        /* updates. */
        upd_bwdsmh(&insbuf.data[j],xs,Ps,&popt->insopt);

        /* solution. */
        ins2sol(&insbuf.data[j],&popt->insopt,&rtk0.sol);

#if OUT_MONITOR
        /* write solution */
        wrt_solution(&rtk0,solopt);
#endif
        /* save solution */
        if (!sve_solution(&rtk0,popt)) {
            trace(2,"add fwd. solution fail\n");
        }
    }
    return bwd.n;
}