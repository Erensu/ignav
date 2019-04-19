/*------------------------------------------------------------------------------
* quaternion.cc : quaternion common functions
*
* reference :
*    [1] http://www.euclideanspace.com/maths/geometry/rotations
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2017/10/06 1.0 new
*-----------------------------------------------------------------------------*/
#include <carvig.h>

static const double ZERO_TOLERANCE=0.000001;
extern const quat_t identity_quat={{1.0,0.0,0.0,0.0}};

/* copy vector vi to vo -----------------------------------------------------*/
static void vec3_copy(vec3_t *vo, vec3_t *vi)
{
    memcpy(vo,vi,sizeof(vec3_t));
}
/* normalize angle -----------------------------------------------------------*/
static double normalize_euler_0_2pi(double a)
{
    while (a<0) a+=(2*M_PI); return a;
}
/* normalize angle-----------------------------------------------------------*/
extern double NORMANG(double ang)
{
    while (ang<0.0) ang+=360.0; return ang;
}
/* returns square of the magnitude of v -------------------------------------*/
static double vec3_len2(const vec3_t *v)
{
    return v->x*v->x+v->y*v->y+v->z*v->z;
}
/* vo = normalized vi -------------------------------------------------------*/
static vec3_t *vec3_normalize(vec3_t *vo, const vec3_t *vi)
{
    double len=sqrt(vec3_len2(vi));
    vo->x=vi->x/len;
    vo->y=vi->y/len;
    vo->z=vi->z/len;
    return vo;
}
/* return dot product of v1 and v2 ------------------------------------------*/
static double vec3_dot(const vec3_t *v1, const vec3_t *v2)
{
    return v1->vec[0]*v2->vec[0]+v1->vec[1]*v2->vec[1]+v1->vec[2]*v2->vec[2];
}
/* vo = v1 x v2 (cross product) ---------------------------------------------*/
static vec3_t *vec3_cross(vec3_t *vo, const vec3_t *v1, const vec3_t *v2)
{
    vo->vec[0]=v1->vec[1]*v2->vec[2]-v1->vec[2]*v2->vec[1];
    vo->vec[1]=v1->vec[2]*v2->vec[0]-v1->vec[0]*v2->vec[2];
    vo->vec[2]=v1->vec[0]*v2->vec[1]-v1->vec[1]*v2->vec[0];
    return vo;
}
/* init orientation quaternion from measurements ----------------------------*/
extern void quat_init(quat_t *q, const vec3_t *acc, const vec3_t *mag)
{
    double ax=acc->x;
    double ay=acc->y;
    double az=acc->z;
    double mx=mag->x;
    double my=mag->y;
    double mz=mag->z;

    double init_roll=atan2(-ay,-az);
    double init_pitch=atan2(ax,-az);

    double cos_roll=cos(init_roll);
    double sin_roll=sin(init_roll);
    double cos_pitch=cos(init_pitch);
    double sin_pitch=sin(init_pitch);

    double mag_x=mx*cos_pitch+my*sin_roll*sin_pitch+mz*cos_roll*sin_pitch;
    double mag_y=my*cos_roll-mz*sin_roll;

    double init_yaw=atan2(-mag_y,mag_x);

    cos_roll=cos(init_roll*0.5);
    sin_roll=sin(init_roll*0.5);

    cos_pitch=cos(init_pitch*0.5);
    sin_pitch=sin(init_pitch*0.5);

    double cosHeading=cos(init_yaw*0.5);
    double sinHeading=sin(init_yaw*0.5);

    q->q0=cos_roll*cos_pitch*cosHeading+sin_roll*sin_pitch*sinHeading;
    q->q1=sin_roll*cos_pitch*cosHeading-cos_roll*sin_pitch*sinHeading;
    q->q2=cos_roll*sin_pitch*cosHeading+sin_roll*cos_pitch*sinHeading;
    q->q3=cos_roll*cos_pitch*sinHeading-sin_roll*sin_pitch*cosHeading;
}
/* initialize quaternion from axis angle using floats -----------------------*/
extern void quat_init_axis(quat_t *q, double x, double y, double z, double a)
{
    /* see: http://www.euclideanspace.com/maths/geometry/rotations
            /conversions/angleToQuaternion/index.htm */
    double a2=a*0.5;
    double s=sin(a2);
    q->x=x*s;
    q->y=y*s;
    q->z=z*s;
    q->w=cos(a2);
}
/* initialize quaternion from axis angle using a vector ---------------------*/
extern void quat_init_axis_v(quat_t *q, const vec3_t *v, double a)
{
    quat_init_axis(q,v->x,v->y,v->z,a);
}
/* extract axis and angle from a quaternion */
extern void quat_to_axis(const quat_t *q, double *x, double *y, double *z, double *a)
{
    /* see: http://www.euclideanspace.com/maths/geometry/rotations
            /conversions/quaternionToAngle/index.htm
    */
    double angle=2.0*acos(q->w);
    double s=sqrt(1.0-q->w*q->w);
    if (s<ZERO_TOLERANCE) {
        /* if s close to zero then direction of axis not important */
        *a=0;
        *x=1;
        *y=0;
        *z=0;
    } else {
        *a=angle;
        *x=*a*q->x/s; /* normalise axis */
        *y=*a*q->y/s;
        *z=*a*q->z/s;
    }
}
/* extract axis in vector form and angle from a quaternion ------------------*/
extern void quat_to_axis_v(const quat_t *q, vec3_t *v, double *a)
{
    quat_to_axis(q,&v->x,&v->y,&v->z,a);
}
/* rotate vector v_in in-place via unit quaternion quat ---------------------*/
extern void quat_rot_vec_self(vec3_t *v, const quat_t *q)
{
    vec3_t vo;
    quat_rot_vec(&vo,v,q);
    vec3_copy(v,&vo);
}
/* rotate vector vi via unit quaternion q and put result into vector vo -----*/
extern void quat_rot_vec(vec3_t *vo, const vec3_t *vi, const quat_t *q)
{
    /* see: https://github.com/qsnake/ase/blob/master/ase/quaternions.py */
    const double vx=vi->x,vy=vi->y,vz=vi->z;
    const double qw=q->w,qx=q->x,qy=q->y,qz=q->z;
    const double qww=qw*qw,qxx=qx*qx,qyy=qy*qy,qzz=qz*qz;
    const double qwx=qw*qx,qwy=qw*qy,qwz=qw*qz,qxy=qx*qy;
    const double qxz=qx*qz,qyz=qy*qz;
    vo->x=(qww+qxx-qyy-qzz)*vx+2*((qxy-qwz)*vy+(qxz+qwy)*vz);
    vo->y=(qww-qxx+qyy-qzz)*vy+2*((qxy+qwz)*vx+(qyz-qwx)*vz);
    vo->z=(qww-qxx-qyy+qzz)*vz+2*((qxz-qwy)*vx+(qyz+qwx)*vy);
}
/* copy quaternion qi to qo -------------------------------------------------*/
extern void quat_copy(quat_t *qo, const quat_t *qi)
{
    memcpy(qo,qi,sizeof(quat_t));
}
/* returns len of quaternion ------------------------------------------------*/
extern double quat_len(const quat_t *q)
{
    int i;
    double s=0.0;
    for (i=0;i<4;i++) s+=q->vec[i]*q->vec[i];
    return sqrt(s);
}
/* conjugate quaternion -----------------------------------------------------*/
extern void quat_conj(quat_t *q_out, const quat_t *q_in)
{
    q_out->x=-q_in->x;
    q_out->y=-q_in->y;
    q_out->z=-q_in->z;
    q_out->w= q_in->w;
}
/* convert quaternion to euler angles ---------------------------------------*/
extern void quat_to_euler(euler_t *euler, const quat_t *quat)
{
    const double x=quat->x,y=quat->y,z=quat->z,w=quat->w;
    const double ww=w*w,xx=x*x,yy=y*y,zz=z*z;
    euler->yaw=normalize_euler_0_2pi(atan2(2.0*(x*y+z*w),xx-yy-zz+ww));
    euler->pitch=asin(-2.0*(x*z-y*w));
    euler->roll=atan2(2.0*(y*z+x*w),-xx-yy+zz+ww);
}
/* o = q1 * q2 --------------------------------------------------------------*/
extern void quat_mul(quat_t *o, const quat_t *q1, const quat_t *q2)
{
    /* see: http://www.euclideanspace.com/maths/algebra/
            realNormedAlgebra/quaternions/code/index.htm#mul */
    o->w=-q1->x*q2->x-q1->y*q2->y-q1->z*q2->z+q1->w*q2->w;
    o->x= q1->x*q2->w+q1->y*q2->z-q1->z*q2->y+q1->w*q2->x;
    o->y=-q1->x*q2->z+q1->y*q2->w+q1->z*q2->x+q1->w*q2->y;
    o->z= q1->x*q2->y-q1->y*q2->x+q1->z*q2->w+q1->w*q2->z;
}
/* o = q1 + 12 --------------------------------------------------------------*/
extern void quat_add(quat_t *o, const quat_t *q1, const quat_t *q2)
{
    /* see: http://www.euclideanspace.com/maths/algebra/
            realNormedAlgebra/quaternions/code/index.htm#add */
    o->x=q1->x+q2->x;
    o->y=q1->y+q2->y;
    o->z=q1->z+q2->z;
    o->w=q1->w+q2->w;
}
/* qo = qi * f --------------------------------------------------------------*/
extern void quat_scale(quat_t *o, const quat_t *q, double f)
{
    /* see: http://www.euclideanspace.com/maths/algebra/
            realNormedAlgebra/quaternions/code/index.htm#scale*/
    int i;
    for (i=0;i<4;i++) o->vec[i]=q->vec[i]*f;
}
/* normalizes quaternion q and puts result into o ---------------------------*/
extern void quat_normalize(quat_t *o, const quat_t *q)
{
    /* see: http://www.euclideanspace.com/maths/algebra/
            realNormedAlgebra/quaternions/code/index.htm#normalise */
    quat_scale(o,q,1.0/quat_len(q));
}
/* normalize q in-place -----------------------------------------------------*/
extern void quat_normalize_self(quat_t *q)
{
    quat_normalize(q,q);
}
/* /* Convert quaternion to right handed rotation matrix. m is a pointer-----
 * to 16 floats in column major order.
 * m is pointer to array of 16 floats in column major order -----------------*/
extern void quat_to_rh_rot_matrix(const quat_t *q, double *m)
{
    quat_t qn;
    double qw,qx,qy,qz;

    quat_normalize(&qn,q);

    qw=qn.w; qx=qn.x; qy=qn.y; qz=qn.z;

    m[0]=1.0-2.0*qy*qy-2.0*qz*qz;
    m[1]=2.0*qx*qy-2.0*qz*qw;
    m[2]=2.0*qx*qz+2.0*qy*qw;

    m[3]=2.0*qx*qy+2.0*qz*qw;
    m[4]=1.0-2.0*qx*qx-2.0*qz*qz;
    m[5]=2.0*qy*qz-2.0*qx*qw;

    m[6]=2.0*qx*qz-2.0*qy*qw;
    m[7]=2.0*qy*qz+2.0*qx*qw;
    m[8]=1.0-2.0*qx*qx-2.0*qy*qy;
}
/* vec3 rotate self by axis and angle ---------------------------------------*/
extern vec3_t *vec3_rot_axis_self(vec3_t *vo, double x, double y, double z,
                                  double angle)
{
    quat_t rotate;
    quat_init_axis(&rotate,x,y,z,angle);
    quat_rot_vec_self(vo,&rotate);
    return vo;
}
/* see http://gamedev.stackexchange.com/questions/15070/orienting-a-model-to-face-a-target */
/* Calculate the quaternion to rotate from vector u to vector v */
extern void quat_from_u2v(quat_t *q, const vec3_t *u, const vec3_t *v,
                          const vec3_t *up)
{
    vec3_t un,vn,axis,axisn;
    double dot;
    double angle;

    vec3_normalize(&un,u);
    vec3_normalize(&vn,v);
    dot=vec3_dot(&un,&vn);
    if (fabs(dot-(-1.0))<ZERO_TOLERANCE) {
        /* vector a and b point exactly in the opposite direction
         * so it is a 180 degrees turn around the up-axis
         */
        vec3_t default_up={{0,1,0}};
        if (!up)
            up=&default_up;
        quat_init_axis(q,up->x,up->y,up->z,M_PI);
        return;
    }
    if (fabs(dot-1.0)<ZERO_TOLERANCE) {
        /* vector a and b point exactly in the same direction
         * so we return the identity quaternion
         */
        *q=identity_quat;
        return;
    }
    angle=acos(dot);
    vec3_cross(&axis,&un,&vn);
    vec3_normalize(&axisn,&axis);
    quat_init_axis(q,axisn.x,axisn.y,axisn.z,angle);
}
/* quaternion dot product q1 . q2 -------------------------------------------*/
extern double quat_dot(const quat_t *q1, const quat_t *q2)
{
    return q1->vec[0]*q2->vec[0]+q1->vec[1]*q2->vec[1]+
           q1->vec[2]*q2->vec[2]+q1->vec[3]*q2->vec[3];
}
static quat_t *quat_lerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto,
                         double t)
{
    double cosom=quat_dot(qfrom,qto);
    double scale0=1.0-t;
    double scale1=t;

    /* qto=qfrom or qto=-qfrom so no rotation to slerp */
    if (cosom>=1.0) {
        quat_copy(qo,qfrom);
        return qo;
    }
    /* adjust for shortest path */
    quat_t to1;
    if (cosom<0.0) {
        to1.x=-qto->x;
        to1.y=-qto->y;
        to1.z=-qto->z;
        to1.w=-qto->w;
    } else {
        quat_copy(&to1,qto);
    }
    /* calculate final values */
    qo->x=scale0*qfrom->x+scale1*to1.x;
    qo->y=scale0*qfrom->y+scale1*to1.y;
    qo->z=scale0*qfrom->z+scale1*to1.z;
    qo->w=scale0*qfrom->w+scale1*to1.w;
    return qo;
}
/* calculate normalized linear quaternion interpolation ---------------------*/
extern quat_t *quat_nlerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto,
                          double t)
{
    quat_lerp(qo,qfrom,qto,t);
    quat_normalize_self(qo);
    return qo;
}
/* calculate spherical quaternion interpolation -----------------------------*/
extern quat_t *quat_slerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto,
                          double t)
{
    /* calc cosine */
    double cosom=quat_dot(qfrom,qto);

    /* qto=qfrom or qto=-qfrom so no rotation to slerp */
    if (cosom>=1.0) {
        quat_copy(qo,qfrom);
        return qo;
    }
    /* adjust for shortest path */
    quat_t to1;
    if (cosom<0.0) {
        cosom=-cosom;
        to1.x=-qto->x;
        to1.y=-qto->y;
        to1.z=-qto->z;
        to1.w=-qto->w;
    } else {
        quat_copy(&to1,qto);
    }
    /* calculate coefficients */
    double scale0,scale1;
    if (cosom<0.99995) {
        /* standard case (slerp) */
        double omega=acos(cosom);
        double sinom=sin(omega);
        scale0=sin((1.0-t)*omega)/sinom;
        scale1=sin(t*omega)/sinom;
    } else {
        /* "from" and "to" quaternions are very close
         *  ... so we can do a linear interpolation
         */
        scale0=1.0-t;
        scale1=t;
    }
    /* calculate final values */
    qo->x=scale0*qfrom->x+scale1*to1.x;
    qo->y=scale0*qfrom->y+scale1*to1.y;
    qo->z=scale0*qfrom->z+scale1*to1.z;
    qo->w=scale0*qfrom->w+scale1*to1.w;
    return qo;
}
/* quaternion inverse----------------------------------------------------*/
extern void quat_inv(const quat_t *q,quat_t *qi)
{
    qi->x=-q->x; qi->y=-q->y; qi->z=-q->z; qi->w=q->w;
}
/* quaternion inverse----------------------------------------------------*/
extern quat_t invquat(const quat_t *q)
{
    quat_t quat;
    quat.x=-q->x; quat.y=-q->y; quat.z=-q->z; quat.w=q->w;
    return quat;
}
/* euler angles convert to quaternion------------------------------------*/
extern void euler_to_quat(const euler_t *euler,quat_t *q)
{
    quat_t quat;

    double cosRoll=cos(euler->roll*0.5);
    double sinRoll=sin(euler->roll*0.5);

    double cosPitch=cos(euler->pitch*0.5);
    double sinPitch=sin(euler->pitch*0.5);

    double cosHeading=cos(euler->yaw*0.5);
    double sinHeading=sin(euler->yaw*0.5);

    quat.w=cosRoll*cosPitch*cosHeading+sinRoll*sinPitch*sinHeading;
    quat.x=sinRoll*cosPitch*cosHeading-cosRoll*sinPitch*sinHeading;
    quat.y=cosRoll*sinPitch*cosHeading+sinRoll*cosPitch*sinHeading;
    quat.z=cosRoll*cosPitch*sinHeading-sinRoll*sinPitch*cosHeading;

    quat_normalize(q,&quat);
}
/* Apply incremental yaw, pitch and roll relative to the quaternion.
 * For example, if the quaternion represents an orientation of a ship,
 * this will apply yaw/pitch/roll *in the ship's local coord system to the
 * orientation
 * --------------------------------------------------------------------------*/
extern quat_t *quat_apply_relative_yaw_pitch_roll(quat_t *q,
                                                  double yaw, double pitch,
                                                  double roll)
{
    quat_t qyaw,qpitch,qroll,qrot,q1,q2,q3,q4;

    /* calculate amount of yaw to impart this iteration... */
    quat_init_axis(&qyaw,0.0,0.0,1.0,yaw);
    /* Calculate amount of pitch to impart this iteration... */
    quat_init_axis(&qpitch,0.0,1.0,0.0,pitch);
    /* Calculate amount of roll to impart this iteration... */
    quat_init_axis(&qroll,1.0,0.0,0.0,roll);
    /* Combine pitch, roll and yaw */
    quat_mul(&q1,&qyaw,&qpitch);
    quat_mul(&qrot,&q1,&qroll);

    /* Convert rotation to local coordinate system */
    quat_mul(&q1,q,&qrot);
    quat_conj(&q2,q);
    quat_mul(&q3,&q1,&q2);
    /* Apply to local orientation */
    quat_mul(&q4,&q3,q);
    quat_normalize_self(&q4);
    *q=q4;
    return q;
}
/* euler angles convert to quaternion----------------------------------------*/
extern void rpy2quat(const double *rpy,quat_t *q)
{
    quat_t quat;
    euler_t euler;

    euler.roll=rpy[0]; euler.pitch=rpy[1]; euler.yaw=rpy[2];

    double cosRoll=cos(euler.roll*0.5);
    double sinRoll=sin(euler.roll*0.5);

    double cosPitch=cos(euler.pitch*0.5);
    double sinPitch=sin(euler.pitch*0.5);

    double cosHeading=cos(euler.yaw*0.5);
    double sinHeading=sin(euler.yaw*0.5);

    quat.w=cosRoll*cosPitch*cosHeading+sinRoll*sinPitch*sinHeading;
    quat.x=sinRoll*cosPitch*cosHeading-cosRoll*sinPitch*sinHeading;
    quat.y=cosRoll*sinPitch*cosHeading+sinRoll*cosPitch*sinHeading;
    quat.z=cosRoll*cosPitch*sinHeading-sinRoll*sinPitch*cosHeading;

    quat_normalize(q,&quat);
}

