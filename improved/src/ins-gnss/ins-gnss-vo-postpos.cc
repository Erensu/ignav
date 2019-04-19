/*------------------------------------------------------------------------------
 * ins-gnss-vo-postpos.cc : ins-gnss-vo coupled post-processing common functions
 *
 * reference :
 *    [1] P.D.Groves, Principles of GNSS, Intertial, and Multisensor Integrated
 *        Navigation System, Artech House, 2008
 *    [2] Tedaldi D, Pretto A, Menegatti E. A robust and easy to implement method
 *        for IMU calibration without external equipments,2014.
 *    [3] Skog I, Handel P. Effects of time synchronization errors in GNSS-aided
 *        INS 2008.
 *    [4] Shin E H. Accuracy Improvement of Low Cost INS/GPS for Land Applications
 *    [5] Li M, Mourikis A I. High-precision, consistent EKF-based visual–inertial
 *        odometry[J].International Journal of Robotics Research,2013,32(6):690-711.
 *    [6] Monocular Visual Inertial Odometry on a Mobile Device.
 *    [7] Mourikis A / , Roumeliotis S / . A Multi-State Constraint Kalman Filter
 *        for Vision-aided Inertial Navigation[C]// IEEE International Conference
 *        on Robotics & Automation. IEEE, 2007.
 *    [8] Forster C , Carlone L , Dellaert F , et al. On-Manifold Preintegration
 *        for Real-Time Visual-Inertial Odometry[J]. IEEE Transactions on Robotics,
 *        2015, 33(1):1-21.
 *    [9] Observability-constrained vision-aided inertial navigation.
 *    [10] Li M ,Mourikis A I. Improving the accuracy of EKF-based visual-inertial
 *         odometry[C]// IEEE International Conference on Robotics & Automation.
 *         IEEE, 2012.
 *    [11] Li M ,Mourikis A I. High-precision, consistent EKF-based visual-inertial
 *         odometry[M].Sage Publications, Inc. 2013.
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/02 1.0 new
 *-----------------------------------------------------------------------------*/
#include "carvig.h"

/* ins-gnss-vo coupled post-processing-----------------------------------------
 * args:    gtime_t ts       I   processing start time (ts.time==0: no limit)
 *          gtime_t te       I   processing end time   (te.time==0: no limit)
 *          double ti        I   processing interval  (s) (0:all)
 *          double tu        I   processing unit time (s) (0:all)
 *          prcopt_t *popt   I   processing options
 *          solopt_t *sopt   I   solution options
 *          filopt_t *fopt   I   file options
 *          char   **infile  I   input files (see below)
 *          int    n         I   number of input files
 *          char   *outfile  I   output file ("":stdout, see below)
 * return: status (0:ok,0>:error,1:aborted)
 * notes : input files should contain observation data, navigation data, precise
 *          ephemeris/clock (optional), sbas log file (optional), ssr message
 *          log file (optional) and tec grid file (optional). only the first
 *          observation data file in the input files is recognized as the rover
 *          data.
 *
 *          the type of an input file is recognized by the file extention as ]
 *          follows:
 *              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
 *              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
 *              .lex,.LEX            : qzss lex message log files
 *              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
 *              .*i,.*I              : tec grid files (ionex)
 *              .fcb,.FCB            : satellite fcb
 *              .imu,.bin            : imu measurement data
 *              .gsof                : gsof measurements data from trimble
 *              .gps                 : NovAtel solution measurement data
 *              .m39                 : M39 timestamp file of tx2-image and M39-imu data
 *              others               : rinex obs, nav, gnav, hnav, qnav or clock
 * ---------------------------------------------------------------------------*/
extern int igvopostpos(const gtime_t ts, gtime_t te, double ti, double tu,
                       const prcopt_t *popt, const solopt_t *sopt,
                       const filopt_t *fopt, char **infile, int n,int port,
                       const char *outfile)
{
    trace(3,"igvopostpos:\n");
    return 0;
}