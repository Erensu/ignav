/*----------------------------------------------------------------------------
* ground-truth.cc : read the ground truth solution for evaluating results
*
* version : $Revision:$ $Date:$
* history : 2017/03/20  1.0  new
*----------------------------------------------------------------------------*/
#include <carvig.h>

/* search ground truth solution from solution buffer-------------------------*/
extern int findgtsols(solbuf_t *solbuf,gtime_t time)
{
    int i;
    if (solbuf->dt==0.0&&solbuf->n>10) {
        for (i=0;i<10;i++) {
            solbuf->dt+=timediff(solbuf->data[i+1].time,solbuf->data[i].time);
        }
        solbuf->dt/=10.0;
    }
    for (i=solbuf->start;i<solbuf->end;i++) {
        if (fabs(timediff(solbuf->data[i].time,time))<(0.5+0.001)*solbuf->dt) {
            solbuf->start=i; return i;
        }
    }
    return -1;
}
/* read karlsruhe dataset ground truth solution------------------------------*/
static int readkarl(FILE *fp,solbuf_t *solbuf)
{
    char buff[MAXBUFF];
    double stamp,pos[3],x0,y0,z0,rpy[3];
    sol_t sol={0};

    trace(3,"readkarl:\n");

    /* read record */
    while (fgets(buff,MAXBUFF,fp)) {
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                   &stamp,&pos[0],&pos[1],&pos[2],&x0,&y0,&z0,&rpy[0],
                   &rpy[1],&rpy[2])<10) {
            continue;
        }
        stamp=stamp*1E-9;
        sol.time.sec =stamp-(time_t)stamp;
        sol.time.time=(time_t)stamp;
        sol.stat=SOLQ_GRTH;

        pos[0]*=D2R; pos[1]*=D2R;
        pos2ecef(pos,sol.rr);
        addsol(solbuf,&sol);
    }
    solbuf->start=0; solbuf->end=solbuf->n;
    return solbuf->n;
}
/* read ground truth solution from file--------------------------------------
 * args  :  char* file       I   input ground truth solution file
 *          solbuf_t *solbuf IO  solution data buffer
 *          int format       I   ground truth solution format
 * return: numbers of ground truth solutions
 * --------------------------------------------------------------------------*/
extern int read_gt_sols(const char *file,solbuf_t *solbuf,int format)
{
    FILE *fp=NULL;

    trace(3,"read_gt_sols: format=%d\n",format);

    if ((fp=fopen(file,"r"))==NULL) {
        trace(2,"open file fail\n");
        return 0;
    }
    initsolbuf(solbuf,0,1024);

    switch (format) {
        case GROUND_TRUTH_KARL : return readkarl(fp,solbuf);
        case GROUND_TRUTH_EUROC: return 0;
        case GROUND_TRUTH_MALA : return 0;
    }
    return 0;
}
