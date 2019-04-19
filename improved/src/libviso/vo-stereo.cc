
#include "vo-stereo.h"

using namespace std;

VisualOdometryStereo::VisualOdometryStereo(parameters param)
                     :param(param),VisualOdometry(param)
{
    matcher->setIntrinsics(param.calib.f,param.calib.cu,param.calib.cv,param.base);
}

VisualOdometryStereo::VisualOdometryStereo() {}

VisualOdometryStereo::~VisualOdometryStereo() {}

bool VisualOdometryStereo::process(uint8_t *I1,uint8_t *I2,int32_t* dims,bool replace) {
  
    /* push back images */
    matcher->pushBack(I1,I2,dims,replace);
  
    /* bootstrap motion estimate if invalid */
    if (!Tr_valid) {
        matcher->matchFeatures(2);
        matcher->bucketFeatures(param.bucket.max_features,
                                param.bucket.bucket_width,
                                param.bucket.bucket_height);
        p_matched=matcher->getMatches();
        updateMotion();
    }
    /* match features and update motion */
    if (Tr_valid) matcher->matchFeatures(2,&Tr_delta);
    else          matcher->matchFeatures(2);
    matcher->bucketFeatures(param.bucket.max_features,
                            param.bucket.bucket_width,
                            param.bucket.bucket_height);
    p_matched=matcher->getMatches();
    return updateMotion();
}

vector<double> VisualOdometryStereo::estimateMotion(vector<Matcher::p_match> p_matched) {

    bool success=true;

    /* compute minimum distance for RANSAC samples */
    double width=0,height=0;
    for (vector<Matcher::p_match>::iterator
                 it=p_matched.begin();it!=p_matched.end();it++) {
        if (it->u1c>width)  width =it->u1c;
        if (it->v1c>height) height=it->v1c;
    }
    double min_dist=min(width,height)/3.0;
  
    /* get number of matches */
    int32_t N=p_matched.size();
    if (N<6) return vector<double>();

    X         =new double[N];
    Y         =new double[N];
    Z         =new double[N];
    J         =new double[4*N*6];
    p_predict =new double[4*N];
    p_observe =new double[4*N];
    p_residual=new double[4*N];

    /* project matches of previous image into 3d */
    for (int32_t i=0;i<N;i++) {
        double d=max(p_matched[i].u1p-p_matched[i].u2p,0.0001f);
        X[i]=(p_matched[i].u1p-param.calib.cu)*param.base/d;
        Y[i]=(p_matched[i].v1p-param.calib.cv)*param.base/d;
        Z[i]=param.calib.f*param.base/d;
    }
    /* loop variables */
    vector<double> tr_delta;
    vector<double> tr_delta_curr;
    tr_delta_curr.resize(6);
  
    /* clear parameter vector */
    inliers.clear();

    /* initial RANSAC estimate */
    for (int32_t k=0;k<param.ransac_iters;k++) {

        /* draw random sample set */
        vector<int32_t> active=getRandomSample(N,3);

        /* clear parameter vector */
        for (int32_t i=0;i<6;i++) tr_delta_curr[i]=0;

        /* minimize reprojection errors */
        VisualOdometryStereo::result result=UPDATED;
        int32_t iter=0;
        while (result==UPDATED) {
            result=updateParameters(p_matched,active,tr_delta_curr,1,1e-6);
            if (iter++>20||result==CONVERGED) break;
        }
        /* overwrite best parameters if we have more inliers */
        if (result!=FAILED) {
            vector<int32_t> inliers_curr=getInlier(p_matched,tr_delta_curr);
            if (inliers_curr.size()>inliers.size()) {
                inliers=inliers_curr;
                tr_delta=tr_delta_curr;
            }
        }
    }
    /* final optimization (refinement) */
    if (inliers.size()>=6) {
        int32_t iter=0;
        VisualOdometryStereo::result result=UPDATED;
        while (result==UPDATED) {
            result=updateParameters(p_matched,inliers,tr_delta,1,1e-8);
            if (iter++>100||result==CONVERGED) break;
        }
        /* not converged */
        if (result!=CONVERGED) success=false;
    }
    else {
        success=false;
    }
    /* release dynamic memory */
    delete[] X;
    delete[] Y;
    delete[] Z;
    delete[] J;
    delete[] p_predict;
    delete[] p_observe;
    delete[] p_residual;
  
    /* parameter estimate succeeded? */
    if (success) return tr_delta;
    else         return vector<double>();
}

vector<int32_t> VisualOdometryStereo::getInlier(vector<Matcher::p_match> &p_matched,vector<double> &tr) {

    /* mark all observations active */
    vector<int32_t> active;
    for (int32_t i=0;i<(int32_t)p_matched.size();i++) active.push_back(i);

    /* extract observations and compute predictions */
    computeObservations(p_matched,active);
    computeResidualsAndJacobian(tr,active);

    /* compute inliers */
    vector<int32_t> inliers;
    for (int32_t i=0;i<(int32_t)p_matched.size();i++)
        if (pow(p_observe[4*i+0]-p_predict[4*i+0],2)+pow(p_observe[4*i+1]-p_predict[4*i+1],2)+
            pow(p_observe[4*i+2]-p_predict[4*i+2],2)+pow(p_observe[4*i+3]-p_predict[4*i+3],2)
            <param.inlier_threshold*param.inlier_threshold)
        inliers.push_back(i);
    return inliers;
}

VisualOdometryStereo::result VisualOdometryStereo::updateParameters(vector<Matcher::p_match> &p_matched,
                                                                    vector<int32_t> &active,
                                                                    vector<double> &tr,
                                                                    double step_size,double eps) {
    /* we need at least 3 observations */
    if (active.size()<3) return FAILED;
  
    /* extract observations and compute predictions */
    computeObservations(p_matched,active);
    computeResidualsAndJacobian(tr,active);

    Matrix A(6,6);
    Matrix B(6,1);

    /* fill matrices A and B */
    for (int32_t m=0;m<6;m++) {
        for (int32_t n=0;n<6;n++) {
            double a=0;
            for (int32_t i=0;i<4*(int32_t)active.size();i++) {
                a+=J[i*6+m]*J[i*6+n];
            }
            A.val[m][n]=a;
        }
        double b=0;
        for (int32_t i=0;i<4*(int32_t)active.size();i++) {
            b+=J[i*6+m]*(p_residual[i]);
        }
        B.val[m][0]=b;
    }
    /* perform elimination */
    if (B.solve(A)) {
        bool converged=true;
        for (int32_t m=0; m<6; m++) {
            tr[m]+=step_size*B.val[m][0];
            if (fabs(B.val[m][0])>eps) converged=false;
        }
        if (converged) return CONVERGED;
        else           return UPDATED;
    }
    else {
        return FAILED;
    }
}

void VisualOdometryStereo::computeObservations(vector<Matcher::p_match> &p_matched,
                                               vector<int32_t> &active) {

    /* set all observations */
    for (int32_t i=0; i<(int32_t)active.size(); i++) {
        p_observe[4*i+0]=p_matched[active[i]].u1c; /* u1 */
        p_observe[4*i+1]=p_matched[active[i]].v1c; /* v1 */
        p_observe[4*i+2]=p_matched[active[i]].u2c; /* u2 */
        p_observe[4*i+3]=p_matched[active[i]].v2c; /* v2 */
    }
}

void VisualOdometryStereo::computeResidualsAndJacobian(vector<double> &tr,
                                                       vector<int32_t> &active) {

    /* extract motion parameters */
    double rx=tr[0]; double ry=tr[1]; double rz=tr[2];
    double tx=tr[3]; double ty=tr[4]; double tz=tr[5];

    /* precompute sine/cosine */
    double sx=sin(rx); double cx=cos(rx); double sy=sin(ry);
    double cy=cos(ry); double sz=sin(rz); double cz=cos(rz);

    /* compute rotation matrix and derivatives : R=Rx*Ry*Rz */
    double r00=+cy*cz;          double r01=-cy*sz;          double r02= +sy;
    double r10=+sx*sy*cz+cx*sz; double r11=-sx*sy*sz+cx*cz; double r12=-sx*cy;
    double r20=-cx*sy*cz+sx*sz; double r21=+cx*sy*sz+sx*cz; double r22=+cx*cy;

    /* --------------------------------------------------------------------------------------------------------
       [                                      0,        -cos(z)*sin(y),                         -cos(y)*sin(z)]
       [   cos(x)*cos(z)*sin(y) - sin(x)*sin(z),  cos(y)*cos(z)*sin(x),   cos(x)*cos(z) - sin(x)*sin(y)*sin(z)]
       [   cos(x)*sin(z) + cos(z)*sin(x)*sin(y), -cos(x)*cos(y)*cos(z),   cos(z)*sin(x) + cos(x)*sin(y)*sin(z)]
       [                                      0,         sin(y)*sin(z),                         -cos(y)*cos(z)]
       [ - cos(z)*sin(x) - cos(x)*sin(y)*sin(z), -cos(y)*sin(x)*sin(z), - cos(x)*sin(z) - cos(z)*sin(x)*sin(y)]
       [   cos(x)*cos(z) - sin(x)*sin(y)*sin(z),  cos(x)*cos(y)*sin(z),   cos(x)*cos(z)*sin(y) - sin(x)*sin(z)]
       [                                      0,                cos(y),                                      0]
       [                         -cos(x)*cos(y),         sin(x)*sin(y),                                      0]
       [                         -cos(y)*sin(x),        -cos(x)*sin(y),                                      0]
     ---------------------------------------------------------------------------------------------------------*/
    double rdrx10=+cx*sy*cz-sx*sz; double rdrx11=-cx*sy*sz-sx*cz; double rdrx12=-cx*cy;
    double rdrx20=+sx*sy*cz+cx*sz; double rdrx21=-sx*sy*sz+cx*cz; double rdrx22=-sx*cy;
    double rdry00=-sy*cz;          double rdry01=+sy*sz;          double rdry02=+cy;
    double rdry10=+sx*cy*cz;       double rdry11=-sx*cy*sz;       double rdry12=+sx*sy;
    double rdry20=-cx*cy*cz;       double rdry21=+cx*cy*sz;       double rdry22=-cx*sy;
    double rdrz00=-cy*sz;          double rdrz01=-cy*cz;
    double rdrz10=-sx*sy*sz+cx*cz; double rdrz11=-sx*sy*cz-cx*sz;
    double rdrz20=+cx*sy*sz+sx*cz; double rdrz21=+cx*sy*cz-sx*sz;

    /* loop variables */
    double X1p,Y1p,Z1p;
    double X1c,Y1c,Z1c,X2c;
    double X1cd,Y1cd,Z1cd;

    /* for all observations do */
    for (int32_t i=0; i<(int32_t)active.size(); i++) {

        /* get 3d point in previous coordinate system */
        X1p=X[active[i]];
        Y1p=Y[active[i]];
        Z1p=Z[active[i]];

        /* compute 3d point in current left coordinate system */
        X1c=r00*X1p+r01*Y1p+r02*Z1p+tx;
        Y1c=r10*X1p+r11*Y1p+r12*Z1p+ty;
        Z1c=r20*X1p+r21*Y1p+r22*Z1p+tz;
    
        /* weighting */
        double weight=1.0;
        if (param.reweighting)
            weight=1.0/(fabs(p_observe[4*i+0]-param.calib.cu)/fabs(param.calib.cu)+0.05);
    
        /* compute 3d point in current right coordinate system */
        X2c=X1c-param.base;

        /* for all paramters do (rx,ry,rz,tx,ty,tz)*/
        for (int32_t j=0;j<6;j++) {

            /* derivatives of 3d pt. in curr. left coordinates wrt. param j */
            switch (j) {
                case 0:
                    X1cd=0;
                    Y1cd=rdrx10*X1p+rdrx11*Y1p+rdrx12*Z1p;
                    Z1cd=rdrx20*X1p+rdrx21*Y1p+rdrx22*Z1p;
                    break;
                case 1:
                    X1cd=rdry00*X1p+rdry01*Y1p+rdry02*Z1p;
                    Y1cd=rdry10*X1p+rdry11*Y1p+rdry12*Z1p;
                    Z1cd=rdry20*X1p+rdry21*Y1p+rdry22*Z1p;
                    break;
                case 2:
                    X1cd=rdrz00*X1p+rdrz01*Y1p;
                    Y1cd=rdrz10*X1p+rdrz11*Y1p;
                    Z1cd=rdrz20*X1p+rdrz21*Y1p;
                    break;
                case 3: X1cd=1; Y1cd=0; Z1cd=0; break;
                case 4: X1cd=0; Y1cd=1; Z1cd=0; break;
                case 5: X1cd=0; Y1cd=0; Z1cd=1; break;
            }
            /* set jacobian entries (project via K) */
            J[(4*i+0)*6+j]=weight*param.calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c); /* left u' */
            J[(4*i+1)*6+j]=weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c); /* left v' */
            J[(4*i+2)*6+j]=weight*param.calib.f*(X1cd*Z1c-X2c*Z1cd)/(Z1c*Z1c); /* right u' */
            J[(4*i+3)*6+j]=weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c); /* right v' */
        }
        /* set prediction (project via K) */
        p_predict[4*i+0]=param.calib.f*X1c/Z1c+param.calib.cu; /* left u */
        p_predict[4*i+1]=param.calib.f*Y1c/Z1c+param.calib.cv; /* left v */
        p_predict[4*i+2]=param.calib.f*X2c/Z1c+param.calib.cu; /* right u */
        p_predict[4*i+3]=param.calib.f*Y1c/Z1c+param.calib.cv; /* right v */
    
        /* set residuals */
        p_residual[4*i+0]=weight*(p_observe[4*i+0]-p_predict[4*i+0]);
        p_residual[4*i+1]=weight*(p_observe[4*i+1]-p_predict[4*i+1]);
        p_residual[4*i+2]=weight*(p_observe[4*i+2]-p_predict[4*i+2]);
        p_residual[4*i+3]=weight*(p_observe[4*i+3]-p_predict[4*i+3]);
    }
}

