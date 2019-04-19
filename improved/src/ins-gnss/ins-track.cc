/*------------------------------------------------------------------------------
 * ins-track.cc : feature points tracking functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Bruce D. Lucas and Takeo Kanade. An Iterative Image Registration Technique
 *        with an Application to Stereo Vision. International Joint Conference on
 *        Artificial Intelligence, pages 674-679, 1981.
 *    [3] Carlo Tomasi and Takeo Kanade. Detection and Tracking of Point Features.
 *        Carnegie Mellon University Technical Report CMU-CS-91-132, April 1991.
 *    [4] Jianbo Shi and Carlo Tomasi. Good Features to Track. IEEE Conference on
 *        Computer Vision and Pattern Recognition, pages 593-600, 1994.
 *    [5] Stan Birchfield. Derivation of Kanade-Lucas-Tomasi Tracking Equation.
 *        Unpublished, January 1997.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/10/25 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* constants------------------------------------------------------------------*/
#define REFINE                 0         /* tracking refine feature points */
#define OUTPUT_PPM             1         /* output tracking data to ppm file */
#define HASH_FIND_TRACK        1         /* find track index from hash table */
#define REPLACE_LOST_FEATURE   1         /* replace lost track feature in `struct: track_t' */
#define MAX_NUM_IMG            16        /* max number of image buffer */

/* type definitions ----------------------------------------------------------*/
typedef struct hashtable {        /* hash table type */
    int index,last_idx;           /* index of track feature for matching */
    UT_hash_handle hh;            /* makes this structure hashable */
} hashtable_t;

typedef struct {                  /* image data buffer type */
    int index;                    /* index of current image */
    img_t imgbuf[MAX_NUM_IMG];    /* image data buffer */
} imgbuf_t;

/* global variables-----------------------------------------------------------*/
static long int id_seed=1;        /* generate a new feature id */
static hashtable_t *hash=NULL;    /* hash table storing all track feature index in `struct:track_t' */
static imgbuf_t imgbuf={0};       /* image data buffer */

/* add element to hash table---------------------------------------------------*/
static int hash_add(hashtable_t **ht,const int index,int last_idx)
{
    struct hashtable *s=NULL;
    HASH_FIND_INT(*ht,&last_idx,s);  /* id already in the hash? */
    if (s==NULL) {

        /* new track feature match index */
        s=(struct hashtable *)malloc(sizeof(struct hashtable));
        s->last_idx=last_idx;
        s->index=index;
        HASH_ADD_INT(*ht,last_idx,s); return -1;
    }
    else {
        trace(2,"hash element have exist\n");

        /* replace this element if it track lost */
        s->index=index;
        s->last_idx=last_idx; return index;
    }
}
/* delete hash table-----------------------------------------------------------*/
static void hash_delete(hashtable_t **ht)
{
    struct hashtable *current,*tmp;

    HASH_ITER(hh,*ht,current,tmp) {
        HASH_DEL(*ht,current);  /* delete; users advances to next */
        free(current);          /* optional- if you want to free  */
    }
    *ht=NULL; /* delete */
}
/* counts of elements in hash table-------------------------------------------*/
static int hash_counts(const struct hashtable *ht)
{
    return HASH_COUNT(ht);
}
/* find track data from hash table--------------------------------------------*/
static int hash_findtrack(hashtable_t **ht,const int last_idx)
{
    struct hashtable *s=NULL;
    int index;

    HASH_FIND_INT(*ht,&last_idx,s);  /* s: output pointer */
    if (s==NULL) return -1;
    index=s->index;
    HASH_DEL(*ht,s); return index;
}
/* find element in hash table-------------------------------------------------*/
static hashtable *hash_find(hashtable_t **ht,const int last_idx)
{
    struct hashtable *s=NULL;

    HASH_FIND_INT(*ht,&last_idx,s);  /* s: output pointer */
    return s;
}
/* initial track---------------------------------------------------------------
 * args:    trackd_t *data  IO  track set data
 *          voopt_t *opt    I   track options
 * return: status (1: ok,0: fail)
 * ----------------------------------------------------------------------------*/
extern int inittrack(trackd_t *data,const voopt_t *opt)
{
    gtime_t t0={0};

    data->ts=data->te=t0;
    data->n=data->nmax=0; data->uid=id_seed++;
    data->first_frame=0;
    data->last_frame =0;
    data->last_idx   =0;
    data->data=NULL;
    return 1;
}
/* find a track in tracking set data------------------------------------------*/
static int findtrack(const track_t *track,const match_point *mp)
{
#if HASH_FIND_TRACK
    return hash_findtrack(&hash,mp->ip);
#else
    int i;
    for (i=0;i<track->n;i++) {
        if (mp->ip==track->data[i].last_idx) return i;
    }
    return -1;
#endif
}
/* add a feature to a frame---------------------------------------------------*/
static void hash_addfeat(feature** feats,int uid,double u,double v,gtime_t time)
{
    feature *s=NULL;
    HASH_FIND_INT(*feats,&uid,s);  /* id already in the hash? */
    if (s==NULL) {

        s=(feature*)malloc(sizeof(feature));
        s->time=time;
        s->uid=uid;
        s->u=u;
        s->v=v;
        HASH_ADD_INT(*feats,uid,s);  /* id: name of key field */
    }
}
/* copy image data------------------------------------------------------------*/
static void copyimg(img_t *out,const img_t *in)
{
    feature *current,*tmp;

    /* update image raw data */
    if (in==NULL) return;
    if (out->data) {
        free(out->data);
        out->data=(unsigned char*)malloc(sizeof(unsigned char)*in->w*in->h);
    }
    else {
        out->data=(unsigned char*)malloc(sizeof(unsigned char)*in->w*in->h);
    }
    memcpy(out->data,in->data,sizeof(unsigned char)*in->w*in->h);
    out->h=in->h;
    out->w=in->w;
    out->id=in->id;
    out->time=in->time;

    /* update feature points */
    HASH_ITER(hh,out->feat,current,tmp) {

        /* delete feature point */
        HASH_DEL(out->feat,current);
        free(current);
    }
    out->feat=NULL;
    HASH_ITER(hh,in->feat,current,tmp) {
        hash_addfeat(&out->feat,current->uid,current->u,current->v,current->time);
    }
}
/* add new feature and image data to track------------------------------------*/
static int addnewfeatimg(trackd_t *track,const feature *feat)
{
    feature *obs_data;

    if (track->nmax<=track->n) {
        if (track->nmax<=0) track->nmax=8; else track->nmax*=2;
        if (!(obs_data=(feature*)realloc(track->data,sizeof(feature)*track->nmax))) {
            free(track->data); track->data=NULL;

            track->n=track->nmax=0;
            return -1;
        }
        track->data=obs_data;
    }
    track->data[track->n++]=*feat;
    return 1;
}
/* add new track to list------------------------------------------------------*/
static int addnewtrack(track_t *track,const trackd_t *data)
{
    trackd_t *obs_data;

    if (track->nmax<=track->n) {
        if (track->nmax<=0) track->nmax=16*2; else track->nmax*=2;
        if (!(obs_data=(trackd_t *)realloc(track->data,sizeof(trackd_t)*track->nmax))) {

            free(track->data);
            track->data=NULL;
            track->n=track->nmax=0; return -1;
        }
        track->data=obs_data;
    }
    track->data[track->n++]=*data;
    return 1;
}
/* replace index of lost track feature----------------------------------------*/
static int indexofrptrack(track_t *track)
{
    int i;
    for (i=0;i<track->n;i++) {
        if (track->data[i].last_idx==-1) return i;
    }
    return -1;
}
/* create a new track---------------------------------------------------------*/
static int newtrack(const match_point_t *mp,gtime_t tp,gtime_t tc,int curr_frame,
                    const voopt_t *opt,track_t *track)
{
    feature fp,fc; trackd_t ntrack,*ptk;
    int index,i;

#if REPLACE_LOST_FEATURE
    /* find index of lost track feature */
    if ((index=indexofrptrack(track))>=0) {

        /* new track. */
        fp.time=tp; fp.u=mp->up; fp.v=mp->vp; fp.valid=1;
        fc.time=tc; fc.u=mp->uc; fc.v=mp->vc; fc.valid=1;
        fp.uid=id_seed;
        fc.uid=id_seed;

        ptk=&track->data[index]; ptk->n=0;
        addnewfeatimg(ptk,&fp);
        addnewfeatimg(ptk,&fc);

        ptk->first_frame=curr_frame-1; /* update frame id */
        ptk->last_frame =curr_frame;
        ptk->last_idx   =mp->ic;
        ptk->uid=id_seed++; /* new feature id */
        ptk->ts=tp;
        ptk->te=tc;
    
        for (i=0;i<3;i++) ptk->xyz[i]=0.0;
        ptk->flag=TRACK_NEW;
        return index;
    }
#endif
    /* new track. */
    fp.time=tp; fp.u=mp->up; fp.v=mp->vp; fp.valid=1;
    fc.time=tc; fc.u=mp->uc; fc.v=mp->vc; fc.valid=1;

    inittrack(&ntrack,opt);

    fp.uid=ntrack.uid;
    fc.uid=ntrack.uid;

    addnewfeatimg(&ntrack,&fp);
    addnewfeatimg(&ntrack,&fc); /* add track feature. */

    ntrack.first_frame=curr_frame-1; /* update frame id */
    ntrack.last_frame =curr_frame;
    ntrack.last_idx   =mp->ic;
    ntrack.ts=tp;
    ntrack.te=tc;
    for (i=0;i<3;i++) ntrack.xyz[i]=0.0;
    ntrack.flag=TRACK_NEW;

    /* add new track to track table */
    if (addnewtrack(track,&ntrack)<=0) return 0;
    return track->n-1;
}
/* add image data to image buffer---------------------------------------------*/
static void addimg2buf(const img_t *data)
{
    copyimg(&imgbuf.imgbuf[imgbuf.index++%MAX_NUM_IMG],data);
}
/* update tracking feature point data-----------------------------------------*/
static int updatetrack(const img_t *pimg,const img_t *cimg,const voopt_t *opt,int curr_frame, gtime_t tp,
                       gtime_t tc,track_t *track)
{
    trackd_t ntrack={0},*ptk;
    feature *current,*tmp;
    int i,find;

    trace(3,"updatetrack:\n");

    HASH_ITER(hh,cimg->feat,current,tmp) {
        for (find=0,i=0;i<track->n;i++) {
            if (current->uid==track->data[i].uid&&track->data[i].last_idx!=-1) {
                find=1;
                break;
            }
        }
        if (find) {
            /* update track */
            addnewfeatimg(&track->data[i],current);

            track->data[i].last_frame=curr_frame;
            track->data[i].te=tc;

            /* update index */
            track->updtrack[track->nupd++]=i;
            track->data[i].flag=TRACK_UPDATED;
        }
        else {
            /* replace lost track */
            for (ptk=NULL,i=0;i<track->n;i++) {
                if (track->data[i].last_idx==-1) {ptk=&track->data[i]; break;}
            }
            if (ptk) {
                ptk->first_frame=curr_frame; /* update frame id */
                ptk->last_frame =curr_frame;
                ptk->last_idx=0;
                ptk->ts=tc;
                ptk->te=tc;
                ptk->uid=current->uid;
                ptk->flag=TRACK_NEW;
                ptk->n=0;

                setzero(ptk->xyz,1,3);
                addnewfeatimg(ptk,current);
            }
            else {
                /* new track */
                inittrack(&ntrack,opt);
                addnewfeatimg(&ntrack,current);

                ntrack.first_frame=curr_frame; /* update frame id */
                ntrack.last_frame =curr_frame;
                ntrack.ts=tc;
                ntrack.te=tc;
                ntrack.flag=TRACK_NEW;
                ntrack.uid =current->uid;

                /* add to track list */
                if (addnewtrack(track,&ntrack)<=0) continue;
            }
            /* update index */
            track->newtrack[track->nnew]=i;
            track->nnew++;
        }
    }
    /* check track lost features */
    for (i=0;i<track->n;i++) {
        if (track->data[i].flag==TRACK_LOST) track->data[i].last_idx=-1;
    }
    /* add image data to buffer */
    addimg2buf(cimg);
    return track->n>0;
}
/* matched feature points set data convert to track set-----------------------
 * args:  match_set *set   I  matched feature points set data
 *        gtime_t tp,tc    I  precious/current time
 *        int curr_frame   I  uid of current frame
 *        img_t *pimg,cimg I  precious/current image data (NULL: no input)
 *        voopt_t *opt     I  options
 *        track_t *track   O  tracking data
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int match2track(const match_set *mset,gtime_t tp,gtime_t tc,int curr_frame,
                       const img_t *pimg,const img_t *cimg,
                       const voopt_t *opt,track_t *track)
{
    register int i,idx=0;
    feature feat;

    trace(3,"match2track: n=%d\n",track->n);

    /* initial feature numbers */
    track->nnew=0;
    track->nupd=0;

    trace(3,"match2track:\n");

    /* initial flag of all track */
    for (i=0;i<track->n;i++) track->data[i].flag=TRACK_LOST;

    /* input detected feature data */
    if (cimg->feat) {
        return updatetrack(pimg,cimg,opt,curr_frame,tp,tc,track);
    }
    /* initial track */
    if (track->n==0) {
        for (i=0;i<mset->n;i++) {

            /* create a new track */
            newtrack(&mset->data[i],tp,tc,curr_frame,opt,track);

            /* update index */
            track->newtrack[track->nnew]=i;
            track->nnew++;
        }
        /* update hash table */
        for (i=0;i<track->nnew;i++) {

            /* add match index to hash table */
            hash_add(&hash,track->newtrack[i],track->data[track->newtrack[i]].last_idx);
        }
        /* add image data to buffer */
        if (pimg) addimg2buf(pimg);
        if (cimg) addimg2buf(cimg);
        return 1;
    }
    /* update feature point tracks data */
    for (i=0;i<mset->n;i++) {
        if ((idx=findtrack(track,&mset->data[i]))<0) {

            /* create a new track */
            if ((idx=newtrack(&mset->data[i],tp,tc,curr_frame,opt,track))<0) continue;

            /* update index */
            track->newtrack[track->nnew]=idx;
            track->nnew++;
            continue;
        }
#if REFINE
        track->data[idx].data[track->data[idx].n-1].u=mset->data[i].up;
        track->data[idx].data[track->data[idx].n-1].v=mset->data[i].vp;
#endif
        /* update track */
        feat.uid=track->data[idx].uid;

        feat.u=mset->data[i].uc;
        feat.v=mset->data[i].vc;

        feat.time=tc;
        feat.valid=1;
        addnewfeatimg(&track->data[idx],&feat);

        /* track flag. */
        track->data[idx].last_idx  =mset->data[i].ic;
        track->data[idx].last_frame=curr_frame;
        track->data[idx].te=tc; /* timestamp */

        /* update index */
        track->updtrack[track->nupd++]=idx;
        track->data[idx].flag=TRACK_UPDATED;
    }
    /* remove all elements in hash table */
    hash_delete(&hash);

    /* add new match index in hash table */
    for (i=0;i<track->nnew;i++) {
        hash_add(&hash,track->newtrack[i],track->data[track->newtrack[i]].last_idx);
    }
    /* update match index in hash table */
    for (i=0;i<track->nupd;i++) {
        hash_add(&hash,track->updtrack[i],track->data[track->updtrack[i]].last_idx);
    }
    for (i=0;i<track->n;i++) {
        if (track->data[i].flag!=TRACK_LOST) continue;
        track->data[i].last_idx=-1;
        track->nlost++;
    }
    /* add image data to buffer */
    addimg2buf(cimg);
    return track->n>0;
}
/* free tracking--------------------------------------------------------------*/
extern void freetrack(trackd_t *track)
{
    trace(3,"freetrack:\n");
    if (track->data) {
        free(track->data); track->data=NULL;
    }
    track->n=track->nmax=0;
    track->first_frame=0;
    track->last_frame =0; track->last_idx=0; track->uid=0;
}
/* free tracking set data ----------------------------------------------------*/
extern void freetrackset(track_t *track)
{
    int i;
    trace(3,"freetrackset:\n");

    if (track->data) {
        for (i=0;i<track->n;i++) freetrack(&track->data[i]);
    }
    track->n=track->nmax=0;
    if (hash) hash_delete(&hash);
}
/* write feature points to PPM file-------------------------------------------
 * args:    feature *featurelist    I  feature points list
 *          int nfeature            I  number of feature points
 *          unsigned char *greyimg  I  image data
 *          int ncols,nrows         I  size of image
 *          char *filename          I  ppm file path
 * return: none
 * ---------------------------------------------------------------------------*/
extern void feat2ppm(feature *featurelist, int nfeature,unsigned char *greyimg,
                     int ncols,int nrows,char *filename)
{
    trace(3,"feat2ppm:\n");

    int nbytes=ncols*nrows*sizeof(char);
    unsigned char *redimg,*grnimg,*bluimg;
    int offset;
    int x,y,xx,yy;
    int i;

    /* Allocate memory for component images */
    redimg=(unsigned char*)malloc(nbytes);
    grnimg=(unsigned char*)malloc(nbytes);
    bluimg=(unsigned char*)malloc(nbytes);
    if (redimg==NULL||grnimg==NULL||bluimg==NULL) return;

    memcpy(redimg,greyimg,nbytes);
    memcpy(grnimg,greyimg,nbytes);
    memcpy(bluimg,greyimg,nbytes);

    /* Overlay features in red */
    for (i=0;i<nfeature;i++) {
        x=(int)(featurelist[i].u+0.5);
        y=(int)(featurelist[i].v+0.5);
        for (yy=y-1;yy<=y+1;yy++) for (xx=x-1;xx<=x+1;xx++) {
            if (xx>=0&&yy>=0&&xx<ncols&&yy<nrows)  {
                offset=yy*ncols+xx;
                *(redimg+offset)=255;
                *(grnimg+offset)=0;
                *(bluimg+offset)=0;
            }
        }
    }
    /* write to PPM file */
    ppmWriteFileRGB(filename,redimg,grnimg,bluimg,ncols,nrows);

    /* free memory */
    free(redimg);
    free(grnimg);
    free(bluimg);
}
/* trace tracking data--------------------------------------------------------*/
extern void tracetrack(const track_t *track)
{
    char dir[126],path[126];
    trackd_t *pta; img_t *imgp;
    int i,j;

    trace(3,"tracetrack: n=%d\n",track->n);

    for (i=0;i<track->n;i++) {
        pta=&track->data[i]; trace(3,"\ntrack-uid=%4d:  \n");
        for (j=0;j<pta->n;j++) {

            /* [time  frame_id  (u,v)] */
            trace(3,"    %s  %4d  (%8.3lf  %8.3lf)  \n",
                  time_str(pta->data[j].time,4),
                  pta->first_frame+j,pta->data[j].u,pta->data[j].v);
        }
        /* output to PPM file */
#if OUTPUT_PPM&TRACE_TRACK
        sprintf(dir,"/media/sujinglan/Files/carvig-debug/test_track/track_%ld",pta->uid);

        /* create new dir. */
        if (access(dir,F_OK)==-1) {
            if (mkdir(dir,0777)!=0) continue;
        }
        for (j=0;j<pta->n;j++) {

            if (!(imgp=getimgdata(pta->data[j].time))) continue;

            sprintf(path,"%s/%d.ppm",dir,j);
            feat2ppm(&pta->data[j],1,imgp->data,imgp->w,
                     imgp->h,path);
        }
#endif
    }
}
/* get track given track-uid--------------------------------------------------
 * args:    track_t *track  I  tracking data
 *          int uid         I  track uid need to find
 * return: pointer of track (NULL: no found)
 * ---------------------------------------------------------------------------*/
extern trackd_t *gettrack(const track_t *track,int uid)
{
    trace(3,"gettrack:\n");
    int i;
    for (i=0;i<track->n;i++) {
        if (track->data[i].uid==uid) return &track->data[i];
    }
    return NULL;
}
/* check feature point whether is out of view or track lost-------------------
 * args:    track_t *track  I  tracking data
 *          voopt_t *opt    I  visual odomtery options
 *          int id          I  feature point id
 *          gtime_t time    I  tracking timestamp
 * return:  status (1: out of view,0: in view,-1: no found)
 * ---------------------------------------------------------------------------*/
extern int outofview(const track_t *track,const voopt_t *opt,int id,gtime_t time)
{
    trackd_t *pt=NULL;
    int i;

    for (i=0;i<track->n;i++) {
        if (track->data[i].uid==id) {pt=&track->data[i]; break;}
    }
    if (!pt) return -1;
    for (i=0;i<pt->n;i++) {
        if (fabs(timediff(time,pt->data[i].time))<1E-6) return 1;
    }
    return 0;
}
/* initial track image data buffer--------------------------------------------*/
extern void inittrackimgbuf(const voopt_t *opt)
{
    gtime_t t0={0};
    int i;
    trace(3,"inittrackimgbuf:\n");

    for (i=0;i<MAX_NUM_IMG;i++) {
        initimg(imgbuf.imgbuf+i,opt->match.img_w,opt->match.img_h,t0);
    }
}
/* free track image data buffer-----------------------------------------------*/
extern void freetrackimgbuf()
{
    trace(3,"freetrackimgbuf:\n");
    int i;
    for (i=0;i<MAX_NUM_IMG;i++) freeimg(imgbuf.imgbuf+i);
}
/* get image data pointer from image data buffer------------------------------*/
extern img_t* getimgdata(gtime_t time)
{
    int i;
    for (i=0;i<MAX_NUM_IMG;i++) {
        if (fabs(timediff(time,imgbuf.imgbuf[i].time))<1E-5) return &imgbuf.imgbuf[i];
    }
    return NULL;
}




