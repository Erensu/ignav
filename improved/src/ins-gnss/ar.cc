/*------------------------------------------------------------------------------
* ar.cc : ambiguity resolution functions
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2012/02/03 1.0  new
*-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants/macros ----------------------------------------------------------*/
#define FIX_THRES       0.129       /* fix threshold (cycle): p0=0.9999 */
#define MIN_ARC_GAP     300.0       /* min arc gap (s) */
#define MIN_FIX_CNT     10          /* min fix count to fix ambiguity */

#define SQR(x)          ((x)*(x))
#define ROUND(x)        (int)floor((x)+0.5)

#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)      /* phase bias (s:satno,f:freq) */

/* wave length of LC (m) -----------------------------------------------------*/
static double lam_LC(int i, int j, int k)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    return CLIGHT/(i*f1+j*f2+k*f5);
}
/* carrier-phase LC (m) ------------------------------------------------------*/
static double L_LC(int i, int j, int k, const double *Li, const double *Lj)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    double L1,L2,L5;

    if ((i&&(!Li[0]||!Lj[0]))||(j&&(!Li[1]||!Lj[1]))||(k&&(!Li[2]||!Lj[2]))) {
        return 0.0;
    }
    L1=CLIGHT/f1*(Li[0]-Lj[0]);
    L2=CLIGHT/f2*(Li[1]-Lj[1]);
    L5=CLIGHT/f5*(Li[2]-Lj[2]);
    return (i*f1*L1+j*f2*L2+k*f5*L5)/(i*f1+j*f2+k*f5);
}
/* pseudorange LC (m) --------------------------------------------------------*/
static double P_LC(int i, int j, int k, const double *Pi, const double *Pj)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    double P1,P2,P5;

    if ((i&&(!Pi[0]||!Pj[0]))||(j&&(!Pi[1]||!Pj[1]))||(k&&(!Pi[2]||!Pj[2]))) {
        return 0.0;
    }
    P1=Pi[0]-Pj[0];
    P2=Pi[1]-Pj[1];
    P5=Pi[2]-Pj[2];
    return (i*f1*P1+j*f2*P2+k*f5*P5)/(i*f1+j*f2+k*f5);
}
/* noise variance of LC (m) --------------------------------------------------*/
static double var_LC(int i, int j, int k, double sig)
{
    const double f1=FREQ1,f2=FREQ2,f5=FREQ5;
    return (SQR(i*f1)+SQR(j*f2)+SQR(k*f5))/SQR(i*f1+j*f2+k*f5)*SQR(sig);
}
/* single-difference noise variance ------------------------------------------*/
static double SD_var(double var, double el)
{
    double sinel=sin(el);
    return 2.0*(var+var/sinel/sinel);
}
/* average LC ----------------------------------------------------------------*/
static void average_LC(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *azel)
{
    ambc_t *amb;
    double LC1,LC2,LC3,var1,var2,var3,err=rtk->opt.err[1]*rtk->opt.eratio[0];
    int i,j,k;

    for (i=0;i<ns;i++) {
        if (satsys(sat[i],NULL)!=SYS_GPS) continue;
        j=iu[i]; k=ir[i];

        /* triple-freq carrier and code LC (m) */
        LC1=L_LC(1,-1, 0,obs[j].L,obs[k].L)-P_LC(1,1,0,obs[j].P,obs[k].P);
        LC2=L_LC(0, 1,-1,obs[j].L,obs[k].L)-P_LC(0,1,1,obs[j].P,obs[k].P);
        LC3=L_LC(1,-6, 5,obs[j].L,obs[k].L)-P_LC(1,1,0,obs[j].P,obs[k].P);

        /* measurement noise variance (m) */
        var1=SD_var(var_LC(1,1,0,err),azel[1+2*iu[i]]);
        var2=SD_var(var_LC(0,1,1,err),azel[1+2*iu[i]]);
        var3=SD_var(var_LC(1,1,0,err),azel[1+2*iu[i]]);

        amb=rtk->ambc+sat[i]-1;

        if (rtk->ssat[sat[i]-1].slip[0]||rtk->ssat[sat[i]-1].slip[1]||
            rtk->ssat[sat[i]-1].slip[2]||amb->n[0]==0||
            fabs(timediff(amb->epoch[0],obs[0].time))>MIN_ARC_GAP) {

            amb->fixcnt=0;
            if (LC1) {
                amb->n[0]=1; amb->LC[0]=LC1; amb->LCv[0]=var1;
            }
            if (LC2) {
                amb->n[1]=1; amb->LC[1]=LC2; amb->LCv[1]=var2;
            }
            if (LC3) {
                amb->n[2]=1; amb->LC[2]=LC3; amb->LCv[2]=var3;
            }
        }
        else { /* averaging */
            if (LC1) {
                amb->LC [0]+=(LC1-amb->LC[0])/(++amb->n[0]);
                amb->LCv[0]+=(var1-amb->LCv[0])/amb->n[0];
            }
            if (LC2) {
                amb->LC [1]+=(LC2-amb->LC[1])/(++amb->n[1]);
                amb->LCv[1]+=(var2-amb->LCv[1])/amb->n[1];
            }
            if (LC3) {
                amb->LC [2]+=(LC3-amb->LC[2])/(++amb->n[2]);
                amb->LCv[2]+=(var3-amb->LCv[2])/amb->n[2];
            }
        }
        amb->epoch[0]=obs[0].time;
    }
}
/* fix solution --------------------------------------------------------------*/
static int fixsol(rtk_t *rtk, const double *v, const double *H, int nv)
{
    double *R; int i,j,info;

    if (nv<=0) return 0;

    R=zeros(nv,nv);

    for (i=0;i<nv;i++) R[i+i*nv]=0.0001;

    /* update states with constraints */
    if (!(info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv))) {

        /* set solution */
        for (i=0;i<rtk->na;i++) {
            rtk->xa[i]=rtk->x[i]; for (j=0;j<rtk->na;j++) {
                rtk->Pa[i+j*rtk->na]=rtk->Pa[j+i*rtk->na]=rtk->P[i+j*rtk->nx];
            }
        }
    }
    else {
        trace(1,"filter error (info=%d)\n",info);
        nv=0;
    }
    free(R);
    return nv;
}
/* resolve integer ambiguity by WL-NL ----------------------------------------*/
extern int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *azel)
{
    ambc_t *ambi,*ambj;
    double lam_N,lam_W,C1,C2,BW,vW,BC,vC,B1,v1,*H,*v;
    int i,j,k,l,NW,N1,N2,nv=0,fix;

    if (ns<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.nf<2) return 0;

    trace(3,"resamb_WLNL: time=%s ns=%d\n",time_str(obs[0].time,0),ns);

    lam_N=lam_LC(1, 1,0);
    lam_W=lam_LC(1,-1,0);
    C1= SQR(lam_carr[1])/(SQR(lam_carr[1])-SQR(lam_carr[0]));
    C2=-SQR(lam_carr[0])/(SQR(lam_carr[1])-SQR(lam_carr[0]));

    v=zeros(ns,1); H=zeros(rtk->nx,ns);

    /* average LC */
    average_LC(rtk,obs,sat,iu,ir,ns,nav,azel);

    /* search reference */
    for (i=0,j=1;j<ns;j++) {
        ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
        if (ambj->n[0]>ambi->n[0]) i=j;
    }
    /* resolve double-difference ambiguity */
    for (j=0;j<ns;j++) {
        if (i==j) continue;
        ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
        k=IB(sat[i],0,&rtk->opt);
        l=IB(sat[j],0,&rtk->opt);

        if (!ambi->n[0]||!ambj->n[0]) continue;

        /* wide lane ambiguity (cycle) */
        BW=(ambi->LC[0]-ambj->LC[0])/lam_W;
        vW=(ambi->LCv[0]/ambi->n[0]+ambj->LCv[0]/ambj->n[0])/SQR(lam_W);
        NW=ROUND(BW);

        /* narrow lane ambiguity (cycle) */
        BC=rtk->x[k]-rtk->x[l];
        vC=rtk->P[k+k*rtk->nx]+rtk->P[l+l*rtk->nx]-2.0*rtk->P[k+l*rtk->nx];
        B1=(BC+C2*lam_carr[1]*NW)/lam_N;
        v1=vC/SQR(lam_N);
        N1=ROUND(B1); N2=N1-NW;

        fix=fabs(NW-BW)<=FIX_THRES*1.5&&sqrt(vW)<=FIX_THRES&&
            fabs(N1-B1)<=FIX_THRES*1.5&&sqrt(v1)<=FIX_THRES;

        trace(3,"%s sat=%2d-%2d B=%13.3f %13.3f sig=%7.3f %7.3f %s\n",
              time_str(obs[0].time,0),sat[i],sat[j],BW,B1,sqrt(vW),
              sqrt(v1),fix?"FIX":"---");

        if (!fix||++ambj->fixcnt<MIN_FIX_CNT) continue;

        /* constraint to dd-ambiguity */
        v[nv]=(C1*lam_carr[0]*N1+C2*lam_carr[1]*N2)-(rtk->x[k]-rtk->x[l]);
        H[k+nv*rtk->nx]= 1.0;
        H[l+nv*rtk->nx]=-1.0;
        nv++;
    }
    /* fix solution */
    nv=fixsol(rtk,v,H,nv);

    free(v); free(H);
    return nv>=ns/2; /* fix if a half ambiguities fixed */
}
/* resolve integer ambiguity by TCAR -----------------------------------------*/
extern int resamb_TCAR(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *azel)
{
    ambc_t *ambi,*ambj;
    double lam_E,lam_F,lam_N,C1,C2,BE,BF,vE,vF,BC,vC,B1,v1,*H,*v;
    int i,j,k,l,nv=0,NE,NF,NW,N1,N2,fix;

    if (ns<=0||rtk->opt.ionoopt!=IONOOPT_IFLC||rtk->opt.nf<3) return 0;

    trace(3,"resamb_TCAR: time=%s ns=%d\n",time_str(obs[0].time,0),ns);

    lam_E=lam_LC(0,1,-1);
    lam_F=lam_LC(1,-6,5);
    lam_N=lam_LC(1, 1,0);
    C1= SQR(lam_carr[1])/(SQR(lam_carr[1])-SQR(lam_carr[0]));
    C2=-SQR(lam_carr[0])/(SQR(lam_carr[1])-SQR(lam_carr[0]));

    v=zeros(ns,1); H=zeros(rtk->nx,ns);

    /* average LC */
    average_LC(rtk,obs,sat,iu,ir,ns,nav,azel);

    /* search reference */
    for (i=0,j=1;j<ns;j++) {
        ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
        if (ambj->n[0]>ambi->n[0]) i=j;
    }
    /* resolve double-difference bias */
    for (j=0;j<ns;j++) {
        if (i==j) continue;
        ambi=rtk->ambc+sat[i]-1;
        ambj=rtk->ambc+sat[j]-1;
        k=IB(sat[i],0,&rtk->opt);
        l=IB(sat[j],0,&rtk->opt);

        if (!ambi->n[1]||!ambj->n[1]||!ambi->n[2]||!ambj->n[2]) continue;

        /* extra wide lane ambiguity (cycle) */
        BE=(ambi->LC[1]-ambj->LC[1])/lam_E;
        BF=(ambi->LC[2]-ambj->LC[2])/lam_F;
        vE=(ambi->LCv[1]/ambi->n[1]+ambj->LCv[1]/ambi->n[1])/SQR(lam_E);
        vF=(ambi->LCv[2]/ambi->n[2]+ambj->LCv[2]/ambi->n[2])/SQR(lam_F);
        NE=ROUND(BE); NF=ROUND(BF); NW=5*NE+NF;

        /* narrow lane ambiguity (cycle) */
        BC=rtk->x[k]-rtk->x[l];
        vC=rtk->P[k+k*rtk->nx]+rtk->P[l+l*rtk->nx]-2.0*rtk->P[k+l*rtk->nx];
        B1=(BC+C2*lam_carr[1]*NW)/lam_N;
        v1=vC/SQR(lam_N);
        N1=ROUND(B1); N2=N1-NW;

        fix=fabs(NE-BE)<=FIX_THRES*1.5&&sqrt(vE)<=FIX_THRES&&
            fabs(NF-BF)<=FIX_THRES*1.5&&sqrt(vF)<=FIX_THRES&&
            fabs(N1-B1)<=FIX_THRES*1.5&&sqrt(v1)<=FIX_THRES;

        trace(1,"%s sat=%2d-%2d B=%13.3f %13.3f %13.3f sig=%7.3f %7.3f %7.3f %s\n",
              time_str(obs[0].time,0),sat[i],sat[j],BE,BF,B1,sqrt(vE),
              sqrt(vF),sqrt(v1),fix?"FIX":"---");

        if (!fix||++ambj->fixcnt<MIN_FIX_CNT) continue;

        /* constraint to dd-ambiguity */
        v[nv]=(C1*lam_carr[0]*N1+C2*lam_carr[1]*N2)-(rtk->x[k]-rtk->x[l]);
        H[k+nv*rtk->nx]= 1.0;
        H[l+nv*rtk->nx]=-1.0;
        nv++;
    }
    /* fix solution */
    nv=fixsol(rtk,v,H,nv);

    free(v); free(H);
    return nv>=ns/2; /* fix if a half ambiguities fixed */
}