/*----------------------------------------------------------------------------
* ins-vel-match-align.cc : ins velocity matching alignment functions
*
* reference :
*    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
*        Navigation System, Artech House, 2008
*    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
*        for IMU calibration without external equipments,2014.
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/11/05 1.0 new
*----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants -----------------------------------------------------------------*/
#define MINVELE     5.0          /* min velocity for ins velocity matching alignment (m/s) */
#define MAXROT      (10.0*D2R)   /* max rotation of vehicle when velocity matching alignment */
#define MINNUMP     5            /* min number of rtk-fix position from message */
#define MININTERNAL 3.0          /* min time internal of rtk-fix position */

/* coarse alignment for velocity matching alignment--------------------------*/
static int cvmali(const double *pos,int np,double *dt,double *Cbn,int *ind,
                  double *ve)
{
    int i,j,k,*index=imat(np,1);
    double *vel=mat(np-1,3),Cnb[9],rpy[3]={0},vn[3],C[9],llh[3];

    for (k=0,i=0;i<np-1;i++) {
        for (j=0;j<3;j++) {
            vel[3*k+j]=(pos[3*(i+1)+j]-pos[3*i+j])/dt[i];
        }
        if (norm(vel+3*k,3)<MINVELE||
                dt[i]>MININTERNAL) continue;

        index[k++]=i; /* index of position */
    }
    if (k>3&&chksdri(vel,k)) {

        *ind=index[k/2]; /* middle index of position */

        matcpy(ve,vel+3*(k/2),1,3);
        ecef2pos(pos+3**ind,llh);
        ned2xyz(llh,C);
        matmul("TN",3,1,3,1.0,C,ve,0.0,vn); /* velocity in navigation frame */

        rpy[2]=vel2head(vn); /* yaw */
    }
    rpy2dcm(rpy,Cnb);
    matt(Cnb,3,3,Cbn);

    free(vel); free(index); return norm(rpy,3)>0.0;
}
/* given gsof message time and return synchronization imu measurement index-*/
static int gt2ii(gtime_t time,int iis,int iin,const imu_t *imu)
{
    int i;
    for (i=iis>iin?iis-iin:0;i<imu->n;i++) {
        if (fabs(timediff(imu->data[i].time,time))<DTTOL) return i;
    }
    return 0;
}
/* coarse alignment for velocity matching alignment through gsof message-----
 * args    :  gsof_data_t *gsof  I  gsof message data
 *            imu_t *imu         I  imu measurement data
 *            int iis            I  start align imu measurement index
 *            int igs            I  start align gsof measurement index
 *            insopt_t *opt      I  ins options
 *            int solq           I  position status
 *            insstate_t *ins    IO ins state
 * return : 1: ok,0: fail
 * -------------------------------------------------------------------------*/
extern int cvmalign(const gsof_data_t *gsof,int igs,const imu_t *imu,int iis,
                    const insopt_t *opt,int solq,insstate_t *ins)
{
    int i,k=0,p,ind=-1,stat=0,iie=iis,*index=imat(MINNUMP,1);
    int iin=(int)opt->hz*2;
    
    double pos[MINNUMP*3],dt[MINNUMP],Cbn[9],Cne[9],ve[3];
    gtime_t time[MINNUMP]={0};

    trace(3,"cvmalign:\n");

    /* initial ins attitude from gsof message */
    for (i=igs;i<gsof->n;i++) {

        if (gsof->data[i].solq==solq) {
            
            iie=gt2ii(gsof->data[i].t,iie,iin,imu);
            if (norm(imu->data[iie].gyro,3)>90.0*D2R) continue;

            matcpy(pos+3*k,gsof->data[i].pos,1,3);

            index[k]=i; /* index of message */
            time[k++]=gsof->data[i].t; /* message time */
        }
        if (k==MINNUMP) {
            for (p=0;p<k-1;p++) dt[p]=timediff(time[p+1],time[p]);
            if ((stat=cvmali(pos,k,dt,Cbn,&ind,ve))) break;
            k=0;
        }
    }
    ind=index[ind]; /* index of gsof message */

    if (stat) {

        /* approximate attitude */
        ned2xyz(gsof->data[ind].llh,Cne);
        matcpy(ins->Cbn,Cbn,3,3);
        matmul("NN",3,3,3,1.0,Cne,Cbn,0.0,ins->Cbe);

        /* ins position */
        iie=gt2ii(gsof->data[ind].t,iis,iin,imu);

        gapv2ipv(gsof->data[ind].pos,ve,ins->Cbe,
                 opt->lever,&imu->data[iie],ins->re,ins->ve);

        /* alignment time */
        ins->time=gsof->data[ind].t;
    }
    free(index); return stat;
}
/* easy velocity matching alignment------------------------------------------
 * args    : gsof_data_t *gsof  I  gsof message data
 *           int igs            I  start align imu measurement index
 *           imu_t *imu         I  imu measurement data
 *           int iis            I  end align imu measurement index
 *           insopt_t *opt      I  ins options
 *           int solq           I  gnss positioning status
 *           insstate_t *ins    O  ins instate
 * return  : 1 (ok) or 0 (fail)
 * -------------------------------------------------------------------------*/
extern int easyvmali(const gsof_data_t *gsof,int igs,const imu_t *imu,int iis,
                     const insopt_t *opt,int solq,insstate_t *ins)
{
    int i,j,k,l;
    int in=(int)opt->hz*2;
    double pos[2*3],dt,ve[3],vn[3],llh[3],C[9],rpy[3]={0};
    gtime_t time[2];

    for (k=0,i=igs;i<gsof->n;i++) {

        /* obtain position from gsof message */
        if (gsof->data[i].solq==solq) {

            pos[3*k  ]=gsof->data[i].pos[0];
            pos[3*k+1]=gsof->data[i].pos[1];
            pos[3*k+2]=gsof->data[i].pos[2];
            
            time[k++]=gsof->data[i].t;
        }
        if (k==2) {
            if (fabs(dt=timediff(time[1],time[0]))
                >MININTERNAL) {
                k=0; continue;
            }
            /* determine ins initial state */
            for (j=0;j<3;j++) ve[j]=(pos[j+3]-pos[j])/dt;

            if (fabs(norm(ve,3))<MINVELE) {k=0; continue;}

            ecef2pos(pos,llh);
            ned2xyz(llh,C);
            matmul("TN",3,1,3,1.0,C,ve,0.0,vn); /* velocity in navigation frame */

            rpy[2]=vel2head(vn); /* yaw */

            rpy2dcm(rpy,C);
            matt(C,3,3,ins->Cbn);

            ned2xyz(llh,C);
            matmul("NN",3,3,3,1.0,C,ins->Cbn,0.0,ins->Cbe);

            if (!(l=gt2ii(time[1],iis,in,imu))) {k=0; continue;}

            gapv2ipv(pos,ve,ins->Cbe,
                     opt->lever,&imu->data[l],ins->re,ins->ve);

            ins->time=time[1]; /* alignment time */
            break;
        }
    }
    return k==2&&i<gsof->n;
}
/* input observation data----------------------------------------------------
 * args   :  int revs        I  input observation data direction
 *           obs_t *obss     I  input all observation data
 *           obsd *obs       O  output observatino data
 *           prcopt_t *popt  I  rtk options
 *           int *iobsu      O  index of next rover station observation data
 *           int *iobsr      O  index of next reference station observation data
 *           int *nuo        O  number of rover station observation data
 *           int *nro        O  number of base station observation data
 * return : number of observation data
 * --------------------------------------------------------------------------*/
extern int inputobs(int revs,const obs_t *obss,obsd_t *obs,const prcopt_t *popt,
                    int *iobsu,int *iobsr,int *nuo,int *nro)
{
    int i,nu,nr,n=0;

    if (!revs) { /* input forward data */
        if ((nu=nextobsf(obss,iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsf(obss,iobsr,2))>0;iobsr+=nr)
                if (timediff(obss->data[*iobsr].time,
                             obss->data[*iobsu].time)>-DTTOL) break;
        }
        else {
            for (i=*iobsr;(nr=nextobsf(obss,&i,2))>0;*iobsr=i,i+=nr)
                if (timediff(obss->data[i].time,
                             obss->data[*iobsu].time)>DTTOL) break;
        }
        nr=nextobsf(obss,iobsr,2);
        if (nr<=0) {
            nr=nextobsf(obss,iobsr,2);
        }
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss->data[*iobsu+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss->data[*iobsr+i];
        *iobsu+=nu;
    }
    else{ /* input backward data */
        if ((nu=nextobsb(obss,iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(nr=nextobsb(obss,iobsr,2))>0;iobsr-=nr)
                if (timediff(obss->data[*iobsr].time,
                             obss->data[*iobsu].time)<DTTOL) break;
        }
        else {
            for (i=*iobsr;(nr=nextobsb(obss,&i,2))>0;*iobsr=i,i-=nr)
                if (timediff(obss->data[i].time,
                             obss->data[*iobsu].time)<-DTTOL) break;
        }
        nr=nextobsb(obss,iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss->data[*iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss->data[*iobsr-nr+1+i];
        *iobsu-=nu;
    }
    if (nuo) *nuo=nu;
    if (nro) *nro=nr;
    return n;
}
/* use single positioning to initial ins states------------------------------
 * args   :  obs_t *obs      I  observation data
 *           imu_t *imu      I  imu measurement data (NULL: no input)
 *           nav_t *nav      I  navigation data
 *           prcopt_t *opt   I  options
 *           insstate_t *ins IO ins states
 *           int *iobss      O  alignment time index for observation
 *           int *iimu       O  alignment time index for imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int alisgpos(const obs_t *obs,const nav_t *nav,const imu_t *imu,
                    const prcopt_t *opt,insstate_t *ins,int *iobss,int *iimu)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t *sol,sol0={{0}};
    const insopt_t *insopt=&opt->insopt;
    int i,j,k,iobs,m,np,info=0,*ind,n;
    char msg[128];
    double *dt,rr[3],vr[3];

    trace(3,"alisgpos:\n");

    n=opt->insopt.minp?opt->insopt.minp:MINNUMP;

    ind=imat(n,1); dt=mat(n,1);

    sol=(sol_t*)malloc(sizeof(sol_t)*n);

    if (sol==NULL||ind==NULL||dt==NULL) {
        trace(2,"malloc error\n");
        return 0;
    }
    for (i=0;i<n;i++) sol[i]=sol0;

    /* input observation data */
    for (np=0,iobs=0;(m=nextobsf(obs,&iobs,1))>0;iobs+=m) {

        /* collect observation data */
        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */

        /* standard positioning */
        if (!pntpos(data,j,nav,opt,sol+np,NULL,NULL,NULL,msg)) continue;

        ind[np++]=iobs;

        if (np==n) {
            for (i=0;i<np-1;i++) {
                dt[i]=timediff(sol[i+1].time,sol[i].time);
            }
            /* check time continuity */
            if (norm(dt,np-1)>SQRT(np-1)
                ||norm(dt,np-1)==0.0) {
                np=0;
                continue;
            }
            else {
                k=(np-1)/2;
                matcpy(rr,sol[k+1].rr+0,1,3);
                matcpy(vr,sol[k+1].rr+3,1,3);

                if (norm(vr,3)==0.0) {
                    for (i=0;i<3;i++) {
                        vr[i]=(sol[k+1].rr[i]-sol[k].rr[i])/dt[k];
                    }
                }
                if (norm(vr,3)<MINVELE) {
                    np=0;
                    continue;
                }
                /* initial ins state use single positioning */
                if (!ant2inins(sol[k+1].time,rr,vr,insopt,imu,ins,iimu)) {
                    np=0;
                    continue;
                }
                /* initial receiver clock bias */
                matcpy(ins->dtr,sol[k+1].dtr,1,4);

                /* initial receiver clock drift */
                ins->dtrr=sol[k+1].dtrr;

                /* alignment time */
                ins->time=sol[k+1].time;

                /* index of align observation data */
                if (iobss) *iobss=ind[k+1];
                info=1;
                break; /* ok */
            }
        }
    }
    if (iobs>=obs->n) {
        trace(2,"no enough observation to align ins\n");
        info=0;
    }
    free(sol);
    free(ind); free(dt);
    return info;
}
/* initial ambiguity for ins update states-----------------------------------*/
static int initamb(insstate_t *ins,const rtk_t *rtk)
{
    int i,na1,na2,na,nb,nx;

    trace(3,"initamb:\n");

    na1=ins->nx-ins->nb;
    na2=rtk->nx-rtk->na;
    na=rtk->na;
    nb=ins->nb;
    nx=rtk->nx;

    /* check number of ambiguity */
    if (na1!=na2) {
        trace(2,"ambiguity initial for ins update fail\n ");
        return 0;
    }
    /* initial ambiguity value and variance */
    for (i=0;i<na1;i++) {
        if (rtk->x[rtk->na+i]!=0.0) {
            insinitx(ins,rtk->x[na+i],rtk->P[(na+i)+(na+i)*nx],nb+i);
        }
    }
    return 1;
}
/* use rtk positioning for initialing ins state------------------------------
 * args   :  obs_t *obs      I  observation data
 *           imu_t *imu      I  imu measurement data (NULL: no input)
 *           nav_t *nav      I  navigation data
 *           prcopt_t *opt   I  options
 *           insstate_t *ins IO ins states
 *           int *iobsu      O  alignment time index for rover observation
 *           int *iobsr      O  alignment time index for base observation
 *           int *iimu       O  alignment time index for imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int alirtkpos(const obs_t *obs,const nav_t *nav,const imu_t *imu,
                     const prcopt_t *opt,insstate_t *ins,int *iobsu,int *iobsr,
                     int *iimu)
{
    gtime_t ts={0};
    sol_t *solp;
    rtk_t rtk={{0}};
    obsd_t data[MAXOBS];
    prcopt_t popt=*opt;
    int nobs,i,np,n,k=0,flag,info=0,nu,nr,ir=0,iu=0;
    double rr[3],vr[3],*dt;

    trace(3,"alirtkpos:\n");

    /* position mode */
    if      (opt->insopt.alimethod==INSALIGN_DGPS) popt.mode=PMODE_DGPS  ;
    else if (opt->insopt.alimethod==INSALIGN_RTK ) popt.mode=PMODE_KINEMA;

    /* ionosphere and troposphere option */
    popt.ionoopt=IONOOPT_IFLC;
    popt.tropopt=TROPOPT_SAAS;

    rtkinit(&rtk,&popt);

    np=opt->insopt.minp?opt->insopt.minp:MINNUMP;
    solp=(sol_t*)malloc(sizeof(sol_t)*np);

    dt=mat(1,np);

    if (solp==NULL||dt==NULL) {
        trace(2,"malloc error\n");
        return 0;
    }
    while ((nobs=inputobs(0,obs,data,&popt,&iu,&ir,&nu,&nr))>=0) {

        /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(data[i].sat,NULL)&popt.navsys)&&
                popt.exsats[data[i].sat-1]!=1) data[n++]=data[i];
        }
        if (n<=0) continue;
#if 1
        if (!screent(data[0].time,ts,ts,1.0)) continue;
#endif
        /* carrier-phase bias correction */
        if (nav->nf>0) {
            corr_phase_bias_fcb(data,n,nav);
        }
        else if (!strstr(popt.pppopt,"-DIS_FCB")) {
            corr_phase_bias_ssr(data,n,nav);
        }
        /* rtk positioning */
        if (!rtkpos(&rtk,data,n,nav)) continue;

        solp[k++]=rtk.sol;

        if (k==np) {
            /* check solution status */
            for (flag=1,i=0;i<np;i++) {
                flag&=solp[i].stat<=SOLQ_DGPS&&solp[i].stat!=SOLQ_NONE;
            }
            /* check ok? */
            if (!flag) {k=0; continue;}

            for (i=0;i<np-1;i++) {
                dt[i]=timediff(solp[i+1].time,solp[i].time);
            }
            /* check time continuity */
            if (norm(dt,np-1)>SQRT(np-1)
                ||norm(dt,np-1)==0.0) {
                k=0; continue;
            }
            matcpy(rr,solp[k-1].rr+0,1,3);
            matcpy(vr,solp[k-2].rr+3,1,3);

            if (norm(vr,3)==0.0) {
                for (i=0;i<3;i++) {
                    vr[i]=(solp[k-1].rr[i]-solp[k-2].rr[i])/dt[k-1];
                }
            }
            if (norm(vr,3)<MINVELE) {
                k=0; continue;
            }
            /* initial ins state using rtk positioning */
            if (!ant2inins(solp[k-1].time,rr,vr,&popt.insopt,
                           imu,ins,iimu)) {
                k=0; continue;
            }
#if 0
            /* ambiguity initial */
            if (!initamb(ins,&rtk)) {k=0; continue;}
#endif
            /* alignment time */
            ins->time=solp[k-1].time;

            /* index of align observation data */
            if (iobsu) *iobsu=iu-nu;
            if (iobsr) *iobsr=ir-nr;

            /* ok */
            info=1; break;
        }
    }
    rtkfree(&rtk);

    free(solp); free(dt);
    return info;
}
/* use ppk positioning for initialing ins state------------------------------
 * args   :  obs_t *obs      I  observation data
 *           imu_t *imu      I  imu measurement data (NULL: no input)
 *           nav_t *nav      I  navigation data
 *           prcopt_t *opt   I  options
 *           insstate_t *ins IO ins states
 *           int *iobsu      O  alignment time index for rover observation
 *           int *iimu       O  alignment time index for imu measurement data
 * return : 1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int alippkpos(const obs_t *obs,const nav_t *nav,const imu_t *imu,
                     const prcopt_t *opt,insstate_t *ins,int *iobsu,int *iimu)
{
    trace(3,"alirtkpos:\n");
    return 1;
}