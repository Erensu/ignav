#include "vo.h"
using namespace std;

VisualOdometry::VisualOdometry(parameters param):param(param) {
    J        =0;
    p_observe=0;
    p_predict=0;
    matcher  =new Matcher(param.match);
    Tr_delta =Matrix::eye(4);
    Tr_valid =false;
    srand(0);
}

VisualOdometry::VisualOdometry() {
    J        =0;
    p_observe=0;
    p_predict=0;
    matcher  =new Matcher(param.match);
    Tr_delta =Matrix::eye(4);
    Tr_valid =false;
    srand(0);
}

VisualOdometry::~VisualOdometry() {
    delete matcher;
}

bool VisualOdometry::updateMotion() {
  
    /* estimate motion */
    std::vector<double> tr_delta=estimateMotion(p_matched);
  
    /* on failure */
    if (tr_delta.size()!=6) return false;
  
    /* set transformation matrix (previous to current frame) */
    Tr_delta=transformationVectorToMatrix(tr_delta);
    Tr_valid=true;
    return true;
}

Matrix VisualOdometry::transformationVectorToMatrix(std::vector<double> tr) {

    /* extract parameters */
    double rx=tr[0],ry=tr[1],rz=tr[2],tx=tr[3],ty=tr[4],tz=tr[5];

    /* precompute sine/cosine */
    double sx=sin(rx),cx=cos(rx),sy=sin(ry),cy=cos(ry),sz=sin(rz),cz=cos(rz);

    /* compute transformation */
    Matrix Tr(4,4);
    Tr.val[0][0]=+cy*cz;          Tr.val[0][1]=-cy*sz;          Tr.val[0][2]=+sy;    Tr.val[0][3]=tx;
    Tr.val[1][0]=+sx*sy*cz+cx*sz; Tr.val[1][1]=-sx*sy*sz+cx*cz; Tr.val[1][2]=-sx*cy; Tr.val[1][3]=ty;
    Tr.val[2][0]=-cx*sy*cz+sx*sz; Tr.val[2][1]=+cx*sy*sz+sx*cz; Tr.val[2][2]=+cx*cy; Tr.val[2][3]=tz;
    Tr.val[3][0]=0;               Tr.val[3][1]=0;               Tr.val[3][2]=0;      Tr.val[3][3]=1;
    return Tr;
}

vector<int32_t> VisualOdometry::getRandomSample(int32_t N,int32_t num) {

    /* init sample and totalset */
    vector<int32_t> sample;
    vector<int32_t> totalset;
  
    /* create vector containing all indices */
    for (int32_t i=0;i<N;i++) totalset.push_back(i);

    /* add num indices to current sample */
    sample.clear();
    for (int32_t i=0;i<num;i++) {
        int32_t j=rand()%totalset.size();
        sample.push_back(totalset[j]);
        totalset.erase(totalset.begin()+j);
    }
    return sample;
}
