/*------------------------------------------------------------------------------
 * ins-gnss-iekf.cc : ins-gnss loosely coupled common functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *        for IMU calibration without external equipments,2014.
 *    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *        INS 2008.
 *    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [5] Barrau A , Bonnabel S . The Invariant Extended Kalman Filter as a Stable
 *        Observer[J]. IEEE Transactions on Automatic Control, 2016:1-1.
 *    [6] Martin P , Salaun E . Invariant observers for attitude and heading
 *        estimation from low-cost inertial and magnetic sensors[C]// IEEE
 *        Conference on Decision & Control. IEEE, 2010.
 *    [7] Barczyk M , Lynch A F . Invariant Extended Kalman Filter design for a
 *        magnetometer-plus-GPS aided inertial navigation system[C]// IEEE
 *        Conference on Decision & Control & European Control Conference.
 *        IEEE, 2012.
 *    [8] Bonnabel S . Left-invariant extended Kalman filter and attitude
 *        estimation[C]// IEEE Conference on Decision & Control. IEEE, 2007.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2019/01/22 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

