/*------------------------------------------------------------------------------
 * ins-pnp.cc : PnP pose estimation algorithm
 *
 * reference :
 *     [1] Lepetit V,Moreno-Noguer F,Fua P. EPnP: An AccurateO(n) Solution
 *         to the PnP Problem[J]. International Journal of Computer Vision,2009,
 *         81(2):155-166.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/11/27 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/*-----------------------------------------------------------------------------*/
static int solve_deg2(double a, double b, double c, double &x1, double &x2)
{
    double delta=b*b-4*a*c;
    if (delta<0) return 0;

    double inv_2a=0.5/a;
    if (delta==0) {
        x1=-b*inv_2a; x2=x1; return 1;
    }
    double sqrt_delta=sqrt(delta);
    x1=(-b+sqrt_delta)*inv_2a;
    x2=(-b-sqrt_delta)*inv_2a;
    return 2;
}
/* ----------------------------------------------------------------------------
 * Reference : Eric W. Weisstein. "Cubic Equation." From MathWorld--A Wolfram
 *             Web Resource.
 * http://mathworld.wolfram.com/CubicEquation.html
 * \return Number of real roots found.
 * ---------------------------------------------------------------------------*/
static int solve_deg3(double a, double b, double c, double d,
                      double &x0, double &x1, double &x2)
{
    if (a==0) {
        if (b==0)	{
            if (c==0) return 0;
            x0=-d/c; return 1;
        }
        x2=0;
        return solve_deg2(b,c,d,x0,x1);
    }
    double inv_a=1.0/a;
    double b_a=inv_a*b,b_a2=b_a*b_a;
    double c_a=inv_a*c;
    double d_a=inv_a*d;
    double Q =(3*c_a-b_a2)/9;
    double R =(9*b_a*c_a-27*d_a-2*b_a*b_a2)/54;
    double Q3=Q*Q*Q;
    double D =Q3+R*R;
    double b_a_3=(1.0/3.0)*b_a;

    if (Q==0) {
        if (R==0) {
            x0=x1=x2=-b_a_3; return 3;
        }
        else {
            x0=pow(2*R,1.0/3.0)-b_a_3;
            return 1;
        }
    }
    if (D<=0) {
        double theta =acos(R/sqrt(-Q3));
        double sqrt_Q=sqrt(-Q);
        x0=2*sqrt_Q*cos( theta      /3.0)-b_a_3;
        x1=2*sqrt_Q*cos((theta+2*PI)/3.0)-b_a_3;
        x2=2*sqrt_Q*cos((theta+4*PI)/3.0)-b_a_3;
        return 3;
    }
    double AD=pow(fabs(R)+sqrt(D),1.0/3.0)*(R>0?1:(R<0?-1:0));
    double BD=(AD==0)?0:-Q/AD;
    x0=AD+BD-b_a_3;
    return 1;
}
/* ----------------------------------------------------------------------------
 * Reference : Eric W. Weisstein. "Quartic Equation." From MathWorld--A Wolfram
 *             Web Resource.
 * http://mathworld.wolfram.com/QuarticEquation.html
 * \return Number of real roots found.
 * ---------------------------------------------------------------------------*/
static int solve_deg4(double a, double b, double c, double d, double e,
                      double &x0, double &x1, double &x2, double &x3)
{
    if (a==0) {
        x3=0; return solve_deg3(b,c,d,e,x0,x1,x2);
    }
    double inv_a=1.0/a;
    b*=inv_a; c*=inv_a; d*=inv_a; e*=inv_a;
    double b2=b*b,bc=b*c,b3=b2*b;

    double r0,r1,r2;
    int n=solve_deg3(1,-c,d*b-4*e,4*c*e-d*d-b2*e,r0,r1,r2);
    if (n==0) return 0;

    double R2=0.25*b2-c+r0,R;
    if (R2<0) return 0;

    R=sqrt(R2);
    double inv_R=1.0/R;
    int nb_real_roots=0;

    double D2,E2;
    if (R<10E-12) {
        double temp=r0*r0-4*e; if (temp<0) D2=E2=-1;
        else {
            double sqrt_temp=sqrt(temp);
            D2=0.75*b2-2*c+2*sqrt_temp; E2=D2-4*sqrt_temp;
        }
    }
    else {
        double u=0.75*b2-2*c-R2,
               v=0.25*inv_R*(4*bc-8*d-b3);
        D2=u+v;
        E2=u-v;
    }
    double b_4=0.25*b,R_2=0.5*R;
    if (D2>=0) {
        double D=sqrt(D2);
        nb_real_roots=2;

        double D_2=0.5*D;
        x0=R_2+D_2-b_4; x1=x0-D;
    }
    if (E2>=0) {
        double E  =sqrt(E2);
        double E_2=0.5*E;
        if (nb_real_roots==0) {
            x0=-R_2+E_2-b_4; x1=x0-E; nb_real_roots=2;
        }
        else {
            x2=-R_2+E_2-b_4; x3=x2-E;
            nb_real_roots=4;
        }
    }
    return nb_real_roots;
}
/*----------------------------------------------------------------------------*/
static int jacobi_4x4(double *A, double *D, double *U)
{
    static double Id[16]={1.0,0.0,0.0,0.0,
                          0.0,1.0,0.0,0.0,
                          0.0,0.0,1.0,0.0,
                          0.0,0.0,0.0,1.0
    };
    double B[4],Z[4],g,h,sum;
    int i,j,k;

    memcpy(U,Id,16*sizeof(double));

    B[0]=A[0]; B[1]=A[5]; B[2]=A[10]; B[3]=A[15];
    memcpy(D,B,4*sizeof(double));
    memset(Z,0,4*sizeof(double));

    for (int iter=0;iter<50;iter++) {
        sum=fabs(A[1])+fabs(A[2])+fabs(A[3])+fabs(A[6])+fabs(A[7])+fabs(A[11]);
        if (sum==0.0) return 1;

        double tresh=(iter<3)?0.2*sum/16.0:0.0;
        for (i=0;i<3;i++) {
            double *pAij=A+5*i+1;

            for (j=i+1;j<4;j++) {
                double Aij=*pAij,eps_machine=100.0*fabs(Aij);

                if (iter>3&&fabs(D[i])+eps_machine==fabs(D[i])
                    &&fabs(D[j])+eps_machine==fabs(D[j])) {
                    *pAij=0.0;
                }
                else if (fabs(Aij)>tresh) {

                    double hh=D[j]-D[i],t;
                    if (fabs(hh)+eps_machine==fabs(hh)) t=Aij/hh;
                    else {
                        double theta=0.5*hh/Aij;
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta<0.0) t=-t;
                    }
                    hh=t*Aij;
                    Z[i]-=hh; Z[j]+=hh;
                    D[i]-=hh; D[j]+=hh; *pAij=0.0;

                    double c=1.0/sqrt(1+t*t),s=t*c,tau=s/(1.0+c);
                    for (k=0;k<=i-1;k++) {
                        g=A[k*4+i]; h=A[k*4+j];
                        A[k*4+i]=g-s*(h+g*tau); A[k*4+j]=h+s*(g-h*tau);
                    }
                    for (k=i+1; k<=j-1;k++) {
                        g=A[i*4+k]; h=A[k*4+j];
                        A[i*4+k]=g-s*(h+g*tau); A[k*4+j]=h+s*(g-h*tau);
                    }
                    for (k=j+1;k<4;k++) {
                        g=A[i*4+k]; h=A[j*4+k];
                        A[i*4+k]=g-s*(h+g*tau); A[j*4+k]=h+s*(g-h*tau);
                    }
                    for (k=0;k<4;k++) {
                        g=U[k*4+i]; h=U[k*4+j];
                        U[k*4+i]=g-s*(h+g*tau);
                        U[k*4+j]=h+s*(g-h*tau);
                    }
                }
                pAij++;
            }
        }
        for (i=0;i<4;i++) B[i]+=Z[i];
        memcpy(D,B,4*sizeof(double));
        memset(Z,0,4*sizeof(double));
    }
    return 0;
}
/* pnp align------------------------------------------------------------------*/
static int pnpalign(double M_end[3][3],
                    double X0,double Y0,double Z0,double X1,double Y1,double Z1,
                    double X2,double Y2,double Z2,double R[3][3],double T[3])
{
    double Cs[3],Ce[3],s[3*3],Qs[16],evs[4],U[16];
    double q[4],ev_max,q02,q0_1,q1_2,q2_3;
    double q12,q0_2,q1_3,q22,q0_3,q32;
    int i,j,i_ev;

    for (i=0;i<3;i++) {
        Ce[i]=(M_end[0][i]+M_end[1][i]+M_end[2][i])/3;
    }
    Cs[0]=(X0+X1+X2)/3;
    Cs[1]=(Y0+Y1+Y2)/3;
    Cs[2]=(Z0+Z1+Z2)/3;

    /* Covariance matrix s: */
    for (j=0;j<3;j++) {
        s[0*3+j]=(X0*M_end[0][j]+X1*M_end[1][j]+X2*M_end[2][j])/3-Ce[j]*Cs[0];
        s[1*3+j]=(Y0*M_end[0][j]+Y1*M_end[1][j]+Y2*M_end[2][j])/3-Ce[j]*Cs[1];
        s[2*3+j]=(Z0*M_end[0][j]+Z1*M_end[1][j]+Z2*M_end[2][j])/3-Ce[j]*Cs[2];
    }
    Qs[0*4+0]=s[0*3+0]+s[1*3+1]+s[2*3+2];
    Qs[1*4+1]=s[0*3+0]-s[1*3+1]-s[2*3+2];
    Qs[2*4+2]=s[1*3+1]-s[2*3+2]-s[0*3+0];
    Qs[3*4+3]=s[2*3+2]-s[0*3+0]-s[1*3+1];

    Qs[1*4+0]=Qs[0*4+1]=s[1*3+2]-s[2*3+1];
    Qs[2*4+0]=Qs[0*4+2]=s[2*3+0]-s[0*3+2];
    Qs[3*4+0]=Qs[0*4+3]=s[0*3+1]-s[1*3+0];
    Qs[2*4+1]=Qs[1*4+2]=s[1*3+0]+s[0*3+1];
    Qs[3*4+1]=Qs[1*4+3]=s[2*3+0]+s[0*3+2];
    Qs[3*4+2]=Qs[2*4+3]=s[2*3+1]+s[1*3+2];

    jacobi_4x4(Qs,evs,U);

    /* Looking for the largest eigen value: */
    i_ev=0;
    ev_max=evs[i_ev];
    for (i=1;i<4;i++) if (evs[i]>ev_max) ev_max=evs[i_ev=i];

    /* Quaternion: */
    for (i=0;i<4;i++) q[i]=U[i*4+i_ev];

    q02 =q[0]*q[0],q12 =q[1]*q[1],q22 =q[2]*q[2],q32=q[3]*q[3];
    q0_1=q[0]*q[1],q0_2=q[0]*q[2],q0_3=q[0]*q[3];
    q1_2=q[1]*q[2],q1_3=q[1]*q[3];
    q2_3=q[2]*q[3];

    R[0][0]=q02+q12-q22-q32;
    R[0][1]=2.0*(q1_2-q0_3);
    R[0][2]=2.0*(q1_3+q0_2);

    R[1][0]=2.0*(q1_2+q0_3);
    R[1][1]=q02+q22-q12-q32;
    R[1][2]=2.0*(q2_3-q0_1);

    R[2][0]=2.0*(q1_3-q0_2);
    R[2][1]=2.0*(q2_3+q0_1);
    R[2][2]=q02+q32-q12-q22;

    for (i=0;i<3;i++) {
        T[i]=Ce[i]-(R[i][0]*Cs[0]+R[i][1]*Cs[1]+R[i][2]*Cs[2]);
    }
    return 1;
}
/* ----------------------------------------------------------------------------
 * Given 3D distances between three points and cosines of 3 angles at the apex,
 * calculates. the lentghs of the line segments connecting projection center (P)
 * and the three 3D points (A, B, C).
 * Returned distances are for |PA|, |PB|, |PC| respectively.
 * Only the solution to the main branch.
 * Reference : X.S. Gao, X.-R. Hou, J. Tang, H.-F. Chang; "Complete Solution
 *             Classification for the Perspective-Three-Point Problem"
 * IEEE Trans. on PAMI, vol. 25, No. 8, August 2003
 * \param lengths3D Lengths of line segments up to four solutions.
 * \param dist3D Distance between 3D points in pairs |BC|, |AC|, |AB|.
 * \param cosines Cosine of the angles /_BPC, /_APC, /_APB.
 * \returns Number of solutions.
 * WARNING: NOT ALL THE DEGENERATE CASES ARE IMPLEMENTED
 * ----------------------------------------------------------------------------*/
static int solve_for_lengths(double lengths[4][3],double distances[3],
                              double cosines[3])
{
    double p=cosines[0]*2;
    double q=cosines[1]*2;
    double r=cosines[2]*2;

    double inv_d22=1.0/(distances[2]*distances[2]);
    double a=inv_d22*(distances[0]*distances[0]);
    double b=inv_d22*(distances[1]*distances[1]);

    double a2=a*a,b2=b*b,p2=p*p,q2=q*q,r2=r*r;
    double pr=p*r,pqr=q*pr;

    /* check reality condition (the four points should not be coplanar) */
    if (p2+q2+r2-pqr-1==0) {
        return 0;
    }
    double ab=a*b,a_2=2*a;
    double A=-2*b+b2+a2+1+ab*(2-r2)-a_2;

    /* check reality condition */
    if (A==0) return 0;

    double a_4=4*a;
    double B=q*(-2*(ab+a2+1-b)+r2*ab+a_4)+pr*(b-b2+ab);
    double C=q2+b2*(r2+p2-2)-b*(p2+pqr)-ab*(r2+pqr)+(a2-a_2)*(2+q2)+2;
    double D=pr*(ab-b2+b)+q*((p2-2)*b+2*(ab-a2)+a_4-2);
    double E=1+2*(b-a-ab)+b2-b*p2+a2;

    double temp=(p2*(a-1+b)+r2*(a-1-b)+pqr-a*pqr);
    double b0=b*temp*temp;

    /* check reality condition */
    if (b0==0) return 0;

    double real_roots[4];
    int n=solve_deg4(A,B,C,D,E,real_roots[0],real_roots[1],real_roots[2],real_roots[3]);

    if (n==0) return 0;

    int nb_solutions=0;
    double r3=r2*r,pr2=p*r2,r3q=r3*q;
    double inv_b0=1.0/b0;

    /* for each solution of x */
    for(int i=0;i<n;i++) {
        double x=real_roots[i];

        /* check reality condition */
        if (x<=0) continue;

        double x2=x*x;
        double b1=
                ((1-a-b)*x2+(q*a-q)*x+1-a+b)*(((r3*(a2+ab*(2-r2)-a_2+b2-2*b+1))*x+

                 (r3q*(2*(b-a2)+a_4+ab*(r2-2) -2)+pr2*(1+a2 +2*(ab-a-b) +r2*(b-b2) +b2)))*x2 +
                 (r3*(q2*(1-2*a+a2)+r2*(b2-ab)-a_4+2*(a2-b2)+2)+r*p2*(b2+2*(ab-b-a)+1+a2)+pr2*q*(a_4+2*(b-ab-a2)-2-r2*b))*x+

                 2*r3q*(a_2-b-a2+ab-1)+pr2*(q2-a_4+2*(a2-b2)+r2*b+q2*(a2-a_2)+2)+
                 p2*(p*(2*(ab-a-b)+a2+b2+1)+2*q*r*(b+a_2-a2-ab-1)));

        /* check reality condition */
        if (b1<=0) continue;
        double y=inv_b0*b1;
        double v=x2+y*y-x*y*r;

        if (v<=0) continue;
        double Z=distances[2]/sqrt(v);
        double X=x*Z;
        double Y=y*Z;

        lengths[nb_solutions][0]=X;
        lengths[nb_solutions][1]=Y;
        lengths[nb_solutions][2]=Z;

        nb_solutions++;
    }
    return nb_solutions;
}
/* internal solve for p3p ----------------------------------------------------*/
static int p3pinternal(double R[4][3][3], double t[4][3], const voopt_t *opt,
                       double mu0, double mv0, double X0, double Y0, double Z0,
                       double mu1, double mv1, double X1, double Y1, double Z1,
                       double mu2, double mv2, double X2, double Y2, double Z2)
{
    double mk0,mk1,mk2;
    double norm,cx_fx,cy_fy;
    double inv_fx,inv_fy,M_orig[3][3];;
    double distances[3],cosines[3],lengths[4][3];
    int n,i,nsols=0;

    cx_fx=opt->cam.K[6]/opt->cam.K[0];
    cy_fy=opt->cam.K[7]/opt->cam.K[4];
    inv_fx=1.0/opt->cam.K[0]; inv_fy=1.0/opt->cam.K[4];

    mu0=inv_fx*mu0-cx_fx;
    mv0=inv_fy*mv0-cy_fy;

    norm=sqrt(mu0*mu0+mv0*mv0+1);
    mk0 =1.0/norm; mu0*=mk0; mv0*=mk0;

    mu1=inv_fx*mu1-cx_fx;
    mv1=inv_fy*mv1-cy_fy;

    norm=sqrt(mu1*mu1+mv1*mv1+1);
    mk1 =1.0/norm; mu1*=mk1; mv1*=mk1;

    mu2=inv_fx*mu2-cx_fx;
    mv2=inv_fy*mv2-cy_fy;

    norm=sqrt(mu2*mu2+mv2*mv2+1);
    mk2 =1.0/norm; mu2*=mk2; mv2*=mk2;

    distances[0]=sqrt((X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2));
    distances[1]=sqrt((X0-X2)*(X0-X2)+(Y0-Y2)*(Y0-Y2)+(Z0-Z2)*(Z0-Z2));
    distances[2]=sqrt((X0-X1)*(X0-X1)+(Y0-Y1)*(Y0-Y1)+(Z0-Z1)*(Z0-Z1));

    /* calculate angles */
    cosines[0]=mu1*mu2+mv1*mv2+mk1*mk2;
    cosines[1]=mu0*mu2+mv0*mv2+mk0*mk2;
    cosines[2]=mu0*mu1+mv0*mv1+mk0*mk1;

    n=solve_for_lengths(lengths,distances,cosines);

    for (i=0;i<n;i++) {
        M_orig[0][0]=lengths[i][0]*mu0;
        M_orig[0][1]=lengths[i][0]*mv0;
        M_orig[0][2]=lengths[i][0]*mk0;

        M_orig[1][0]=lengths[i][1]*mu1;
        M_orig[1][1]=lengths[i][1]*mv1;
        M_orig[1][2]=lengths[i][1]*mk1;

        M_orig[2][0]=lengths[i][2]*mu2;
        M_orig[2][1]=lengths[i][2]*mv2;
        M_orig[2][2]=lengths[i][2]*mk2;

        if (!pnpalign(M_orig,X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,R[nsols],t[nsols])) {
            continue;
        }
        nsols++;
    }
    return nsols;
}
/* p3p solver-----------------------------------------------------------------*/
static int p3psolve(double R[3][3], double t[3], const voopt_t *opt,
                    double mu0, double mv0, double X0, double Y0, double Z0,
                    double mu1, double mv1, double X1, double Y1, double Z1,
                    double mu2, double mv2, double X2, double Y2, double Z2,
                    double mu3, double mv3, double X3, double Y3, double Z3)
{
    double Rs[4][3][3],ts[4][3];
    int ns=0,i,j,n;
    double min_reproj=0;

    n=p3pinternal(Rs,ts,opt,mu0,mv0,X0,Y0,Z0,mu1,mv1,X1,Y1,Z1,mu2,mv2,X2,Y2,Z2);

    if (n==0) return 0;
    for (i=0;i<n;i++) {
        double X3p=Rs[i][0][0]*X3+Rs[i][0][1]*Y3+Rs[i][0][2]*Z3+ts[i][0];
        double Y3p=Rs[i][1][0]*X3+Rs[i][1][1]*Y3+Rs[i][1][2]*Z3+ts[i][1];
        double Z3p=Rs[i][2][0]*X3+Rs[i][2][1]*Y3+Rs[i][2][2]*Z3+ts[i][2];

        double mu3p=opt->cam.K[6]+opt->cam.K[0]*X3p/Z3p;
        double mv3p=opt->cam.K[7]+opt->cam.K[4]*Y3p/Z3p;
        double reproj=(mu3p-mu3)*(mu3p-mu3)+(mv3p-mv3)*(mv3p-mv3);
        if (i==0||min_reproj>reproj) {
            ns=i; min_reproj=reproj;
        }
    }
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) R[i][j]=Rs[ns][i][j];
        t[i]=ts[ns][i];
    }
    return 1;
}
/* p3p pose estimation use only three points----------------------------------
 * args:    feature *feats  I  feature points
 *          int nf          I  number of feature points
 *          double *xp      I  points of world frame
 *          int np          I  number of points in world frame
 *          voopt_t *opt    I  options
 *          double *R,*t    O  rotation matrix and translation vector
 * return: status (1: ok,0: fail)
 * note: input more than four points, use three to solve, use others to validate
 * ---------------------------------------------------------------------------*/
extern int p3pthree(const feature *feats,int nf,double *xp,int np,
                    const voopt_t *opt,double *R,double *t)
{
    double Rot[3][3],trans[3];

    trace(3,"p3pthree:\n");

    if (nf<4||np<4) {
        trace(2,"need more than four feature points\n");
        return 0;
    }
    double mp[2*4];
    int i,j;

    for (i=0;i<4;i++) {
        mp[2*i+0]=feats[i].u; mp[2*i+1]=feats[i].v;
    }
    if (!p3psolve(Rot,trans,opt,mp[0],mp[1],xp[0],xp[1],xp[2],
                  mp[2],mp[3],xp[3],xp[ 4],xp[ 5],
                  mp[4],mp[5],xp[6],xp[ 7],xp[ 8],
                  mp[6],mp[7],xp[9],xp[10],xp[11])) {
        trace(2,"p3pthree solve fail\n");
        seteye(R,3); setzero(t,1,3);
        return 0;
    }
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) R[i+3*j]=Rot[i][j];
        t[i]=trans[i];
    }
    return 1;
}
/* p3p pose estimation--------------------------------------------------------
 * args:    feature *feats  I  feature points
 *          int nf          I  number of feature points
 *          double *xp      I  points of world frame
 *          int np          I  number of points in world frame
 *          voopt_t *opt    I  options
 *          double *R,*t    O  rotation matrix and translation vector
 * return: status (1: ok,0: fail)
 * ---------------------------------------------------------------------------*/
extern int p3p(const feature *feats,int nf,const double *xp,int np,
               const voopt_t *opt,double *R,double *t)
{
    trace(3,"p3p:\n");
}