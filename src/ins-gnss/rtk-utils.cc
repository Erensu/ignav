/*------------------------------------------------------------------------------
* rtk-utils.cc : rtk precise positioning utils
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2017/12/25 1.0  new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* add double-difference ambiguity-------------------------------------------*/
extern int add_ddamb(amb_t *amb)
{
    ddamb_t *data;

    if (amb->nmax<=amb->nb) {
        if (amb->nmax<=0) amb->nmax=MAXSAT; else amb->nmax*=2;
        if (!(data=(ddamb_t *)realloc(amb->amb,sizeof(ddamb_t)*amb->nmax))) {
            trace(2,"add double-difference ambiguity: realloc error n=%dx%d\n",
                  sizeof(ddamb_t),amb->nmax);
            free(amb->amb); amb->amb=NULL; amb->nb=amb->nmax=0;
            return -1;
        }
        amb->amb=data;
    }
    return 1;
}
/* get double-difference ambiguity -------------------------------------------*/
extern ddamb_t *get_ddamb(amb_t *bias,int sat1,int sat2,int f)
{
    int i; for (i=0;i<bias->nb;i++)
        if (sat1==bias->amb[i].sat1&&
            sat2==bias->amb[i].sat2&&f==bias->amb[i].f) return bias->amb+i;
    return NULL;
}
/* extract double-difference ambiguity----------------------------------------*/
extern int extract_ddamb(rtk_t *rtk, const ddsat_t *ddsat, int nb,
                         const double *bias)
{
    int i,j,sat1,sat2,f;
    ddamb_t *amb=NULL,amb0={0};
    double p,pr;

    trace(3,"extract_ddamb:\n");

    /* get double-difference ambiguity list */
    for (j=0,i=0;i<nb;i++) {

        f=ddsat[i].f;
        sat1=ddsat[i].sat1;
        sat2=ddsat[i].sat2;

        if ((amb=get_ddamb(&rtk->bias,sat1,sat2,f))==NULL) {
            add_ddamb(&rtk->bias);

            /* new double-difference ambiguity */
            amb=&rtk->bias.amb[rtk->bias.nb++];
            *amb=amb0;
        }
        /* set observation time and numbers of ambiguity */
        amb->pt=amb->time; amb->time=rtk->sol.time;
        p=amb->b;
        pr=amb->ratio;

        /* updates double-difference ambiguity */
        amb->pb=amb->b; amb->b=bias[i]; /* double-difference ambiguity */
        amb->sat1=sat1; amb->sat2=sat2; /* double-difference satellite */
        amb->f=ddsat[i].f;              /* frequency no. */

        /* first epoch ambiguity */
        if (amb->c==0) {
            amb->c++; continue;
        }
        /* check dd-ambiguity */
        if      (fabs(p-amb->b)<1E-3) amb->c++;
        else if (fabs(p-amb->b)>1E-3) {

            if (fabs(p-amb->b)<5.0
                &&pr>rtk->sol.ratio) {

                /* check cycle slip */
                if (!rtk->ssat[sat1-1].slip[f]&&
                    !rtk->ssat[sat2-1].slip[f]) {

                    /* no cycle slip,no change ambiguity */
                    amb->b=p; amb->pb=p; amb->c++;
                }
            }
            else {
                amb->c=0; /* reset fix counts */
            }
        }
        /* ratio */
        if (rtk->sol.stat==SOLQ_INHERIT) {
            amb->ratio=pr; /* inherit precious ratio */
        }
        else {
            amb->ratio=rtk->sol.ratio; /* current ratio */
        }
    }
#ifdef TRACE
    for (i=0;i<rtk->bias.nb;i++) {
        trace(3,"sat=%3d - %3d (frq=%d) : %s  %12.4lf  %3d  %12.4lf\n",
              rtk->bias.amb[i].sat1,
              rtk->bias.amb[i].sat2,rtk->bias.amb[i].f+1,
              time_str(rtk->bias.amb[i].time,2),
              rtk->bias.amb[i].b,rtk->bias.amb[i].c,
              rtk->sol.ratio);
    }
#endif
    return j; /* numbers of double-difference ambiguity */
}

