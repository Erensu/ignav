/*------------------------------------------------------------------------------
* carvig.h : computer vision, ins and gnss coupled for positioning header file
*
* version : $Revision:$ $Date:$
* history : 2018/11/09 1.0  new (first version)
*-----------------------------------------------------------------------------*/
#ifndef CARVIG_H
#define CARVIG_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include <sys/stat.h>
#include <stddef.h>
#include <byteswap.h>
#include <termios.h>
#include <queue.h>
#include <dirent.h>
#include <unistd.h>
#include <malloc.h>

#include <structure/uthash.h>

/* define fixed-width datatypes for Visual Studio projects */
#if defined(_MSC_VER)&&(_MSC_VER<1600)
typedef __int8            int8_t;
typedef __int16           int16_t;
typedef __int32           int32_t;
typedef __int64           int64_t;
typedef unsigned __int8   uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t;
#else
#include <stdint.h>
#endif

#ifndef uchar
typedef unsigned char uchar;
#endif

#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <pthread.h>
#endif
#ifdef __cplusplus
extern "C"{
#endif

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport)    /* for Windows DLL */
#else
#define EXPORT
#endif
#define ENACMP
#define ENAOPENCV   0                   /* enable opencv library */

/* constants -----------------------------------------------------------------*/
#define VER_CARVIG  "1.0.0"             /* library version */

#define PATCH_LEVEL "b1"                /* patch level */
#define INSPOSSTR   "@INS@"             /* string to express ins+gnss position mode */
#define GNSPOSSTR   "@GNS@"             /* string to express gnss position mode */
#define VOSPOSSTR   "@VOS@"             /* string to express vo position mode */
#define GTHPOSSTR   "@GROUNDTRUTH@"     /* string to express ground truth mode */

#define OBS_200HZ
#define CPUTIME_IN_GPST
#define EXTGSI
#define LAPACK
#define TRACE                           /* trace information for debug */
#define TRACE_INS     1                 /* trace ins updates information */
#define TRACE_STDERR  0                 /* trace information to stderr if set */
#define VIG_TRACE_MAT 1                 /* trace matrix for debugs */
#define MAXBUFF     4096                /* size of line buffer */
#define MAXHIST     256                 /* size of history buffer */
#define DISFLAG     123456              /* flag of disable estimated state */
#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define M2R         (PI/(60.0*180.0))   /* min to rad */
#define R2M         ((60.0*180.0)/PI)   /* rad to min */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */
#define Mg2M        9.80665E-6          /* micro-g to meters per second squared */

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */

#define HION        350000.0            /* ionosphere height (m) */

#define MAXFREQ     7                   /* max NFREQ */

#define FREQ1       1.57542E9           /* L1/E1  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2     frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/LEX frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9          /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ1_CMP   1.561098E9          /* BeiDou B1 frequency (Hz) */
#define FREQ2_CMP   1.20714E9           /* BeiDou B2 frequency (Hz) */
#define FREQ3_CMP   1.26852E9           /* BeiDou B3 frequency (Hz) */

#define FREQOCXO    1E8                 /* crystal frequency (100Mhz) */

#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   1.5                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_QZS   1.0                 /* error factor: QZSS */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_IRN   1.5                 /* error factor: IRNSS */
#define EFACT_SBS   3.0                 /* error factor: SBAS */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_CMP     0x20                /* navigation system: BeiDou */
#define SYS_IRN     0x40                /* navigation system: IRNS */
#define SYS_LEO     0x80                /* navigation system: LEO */
#define SYS_ALL     0xFF                /* navigation system: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
#define TSYS_IRN    6                   /* time system: IRNSS time */

#ifndef NFREQ
#define NFREQ       3                   /* number of carrier frequencies */
#endif
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */

#ifndef NEXOBS
#define NEXOBS      10                  /* number of extended obs codes */
#endif

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#ifdef ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif
#ifdef ENAGAL
#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   30                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1
#else
#define MINPRNGAL   0
#define MAXPRNGAL   0
#define NSATGAL     0
#define NSYSGAL     0
#endif
#ifdef ENAQZS
#define MINPRNQZS   193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS   199                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                 /* min satellite PRN number of QZSS SAIF */
#define MAXPRNQZS_S 189                 /* max satellite PRN number of QZSS SAIF */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif
#ifdef ENACMP
#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   35                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1
#else
#define MINPRNCMP   0
#define MAXPRNCMP   0
#define NSATCMP     0
#define NSYSCMP     0
#endif
#ifdef ENAIRN
#define MINPRNIRN   1                   /* min satellite sat number of IRNSS */
#define MAXPRNIRN   7                   /* max satellite sat number of IRNSS */
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1) /* number of IRNSS satellites */
#define NSYSIRN     1
#else
#define MINPRNIRN   0
#define MAXPRNIRN   0
#define NSATIRN     0
#define NSYSIRN     0
#endif
#ifdef ENALEO
#define MINPRNLEO   1                   /* min satellite sat number of LEO */
#define MAXPRNLEO   10                  /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif
#define NSYS        (NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP+NSYSIRN+NSYSLEO) /* number of systems */
#define NUMSYS      7

#define MINPRNSBS   120                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   142                 /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */

#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
                                        /* max satellite number (1 to MAXSAT) */
#define MAXSTA      255

#ifndef MAXOBS
#define MAXOBS      64                  /* max number of obs in an epoch */
#endif
#define MAXRCV      64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */

#ifdef  OBS_200HZ
#define DTTOL       0.0025              /* tolerance of time difference (s) */
#else
#define DTTOL       0.005               /* tolerance of time difference (s) */
#endif

#define MAXDTOE     7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 10800.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_SBS 360.0               /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S   86400.0             /* max time difference to ephem toe (s) for other */
#define MAXGDOP     300.0               /* max GDOP */

#define INT_SWAP_TRAC 86400.0           /* swap interval of trace file (s) */
#define INT_SWAP_STAT 86400.0           /* swap interval of solution status file (s) */

#define MAXEXFILE   1024                /* max number of expanded files */
#define MAXSBSAGEF  30.0                /* max age of SBAS fast correction (s) */
#define MAXSBSAGEL  1800.0              /* max age of SBAS long term corr (s) */
#define MAXSBSURA   8                   /* max URA of SBAS satellite */
#define MAXBAND     10                  /* max SBAS band of IGP */
#define MAXNIGP     201                 /* max number of IGP in SBAS band */
#define MAXNGEO     4                   /* max number of GEO satellites */
#define MAXCOMMENT  10                  /* max number of RINEX comments */
#define MAXSTRPATH  1024                /* max length of stream path */
#define MAXSTRMSG   1024                /* max length of stream message */
#define MAXSTRRTK   12                  /* max number of stream in RTK server */
#define MAXSBSMSG   32                  /* max number of SBAS msg in RTK server */
#define MAXSOLMSG   8191                /* max length of solution message */
#define MAXRAWLEN   4096                /* max length of receiver raw message */
#define MAXERRMSG   4096                /* max length of error/warning message */
#define MAXANT      64                  /* max length of station name/antenna type */
#define MAXSOLBUF   500                 /* max number of solution buffer */
#define MAXOBSBUF   128                 /* max number of observation data buffer */
#define MAXIMU      500                 /* max number of imu measurement data buffer */
#define MAXSOL      100                 /* max number of solution data buffer */
#define MAXIMUBUF   1000                /* max number of imu measurement data buffer */
#define MAXIMGBUF   500                 /* max number of image data buffer */
#define MAXPOSEBUF  1000                /* max number of pose measurement buffer */
#define MAXFRAMEBUF 100                 /* max number of image frames buffer */
#define MAXPOSE     50                  /* max number of input pose measurement data */
#define MAXIMG      50                  /* max number of input image data */
#define MAXNRPOS    16                  /* max number of reference positions */
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define MAXGISLAYER 32                  /* max number of GIS data layers */
#define MAXRCVCMD   4096                /* max length of receiver commands */

#define RNX2VER     2.10                /* RINEX ver.2 default output version */
#define RNX3VER     3.00                /* RINEX ver.3 default output version */

#define OBSTYPE_PR  0x01                /* observation type: pseudorange */
#define OBSTYPE_CP  0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP 0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR 0x08                /* observation type: SNR */
#define OBSTYPE_ALL 0xFF                /* observation type: all */

#define FREQTYPE_L1 0x01                /* frequency type: L1/E1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/B1 */
#define FREQTYPE_L5 0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_L6 0x08                /* frequency type: E6/LEX/B3 */
#define FREQTYPE_L7 0x10                /* frequency type: E5b/B2 */
#define FREQTYPE_L8 0x20                /* frequency type: E5(a+b) */
#define FREQTYPE_L9 0x40                /* frequency type: S */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */

#define FEAT_CREATE           0         /* feature point status: created */
#define FEAT_INIT_FAIL        1         /* feature point status: feature point position initial fail */

#define TRACK_UPDATED         0         /* feature point track flag: updated */
#define TRACK_NEW             1         /* feature point track flag: new track */
#define TRACK_LOST            2         /* feature point track flag: track lost */
#define TRACK_BAD             3         /* feature point track flag: bad track */
#define TRACK_FILTER          4         /* feature point track flag: have used to filter */
#define TRACK_INIT_POS_OK     5         /* feature point track flag: initial feature point position ok */
#define TRACK_INIT_POS_FAIL   6         /* feature point track flag: initial feature point position fail */
#define TRACK_BUCKET          7         /* feature point track flag: track have used to bucket robust estimate */

#define KLT_INIT              1         /* KLT track status: initial */
#define KLT_TRACKED           0         /* KLT track status: tracked */
#define KLT_NOT_FOUND        -1         /* KLT track status: not found */
#define KLT_SMALL_DET        -2         /* KLT track status: small motion */
#define KLT_MAX_ITERATIONS   -3         /* KLT track status: max iterations */
#define KLT_OOB              -4         /* KLT track status: out of border */
#define KLT_LARGE_RESIDUE    -5         /* KLT track status: large residue */

#define TRACE_TRACK           1         /* trace track data for debugs */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P    (GPS,GLO) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless (GPS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A        (GAL) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P) (GAL,QZS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1SAIF (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1I+Q (GPS,QZS,CMP) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5/E5aI    (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5/E5aQ    (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5/E5aI+Q/L5B+C (GPS,GAL,QZS,IRN,SBS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2I   (GAL,CMP) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2Q   (GAL,CMP) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2I+Q (GAL,CMP) */
#define CODE_L6A    30                  /* obs code: E6A        (GAL) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,CMP) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C    (GAL) */
#define CODE_L6S    35                  /* obs code: LEXS       (QZS) */
#define CODE_L6L    36                  /* obs code: LEXL       (QZS) */
#define CODE_L8I    37                  /* obs code: E5(a+b)I   (GAL) */
#define CODE_L8Q    38                  /* obs code: E5(a+b)Q   (GAL) */
#define CODE_L8X    39                  /* obs code: E5(a+b)I+Q (GAL) */
#define CODE_L2I    40                  /* obs code: B1I        (BDS) */
#define CODE_L2Q    41                  /* obs code: B1Q        (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define MAXCODE     55                  /* max number of obs code */

#define PMODE_SINGLE 0                  /* positioning mode: single */
#define PMODE_DGPS   1                  /* positioning mode: DGPS/DGNSS */
#define PMODE_KINEMA 2                  /* positioning mode: kinematic */
#define PMODE_STATIC 3                  /* positioning mode: static */
#define PMODE_MOVEB  4                  /* positioning mode: moving-base */
#define PMODE_FIXED  5                  /* positioning mode: fixed */
#define PMODE_PPP_KINEMA 6              /* positioning mode: PPP-kinemaric */
#define PMODE_PPP_STATIC 7              /* positioning mode: PPP-static */
#define PMODE_PPP_FIXED  8              /* positioning mode: PPP-fixed */
#define PMODE_INS_UPDATE 9              /* positioning mode: ins mechanization */
#define PMODE_INS_LGNSS  10             /* positioning mode: ins/gnss loosely coupled */
#define PMODE_INS_TGNSS  11             /* positioning mode: ins/gnss tightly coupled */
#define PMODE_INS_LGNSS_VO  12          /* positioning mode: ins/gnss-lc/vo coupled */
#define PMODE_INS_TGNSS_VO  13          /* positioning mode: ins/gnss-tc/vo coupled */

#define CAM_PINHOLE      0              /* camera model: pinhole */
#define DISTORT_RT       0              /* distortion model: radial-tangential */

#define INSTC_SINGLE   1                /* ins and gnss tightly coupled mode: single */
#define INSTC_PPK      2                /* ins and gnss tightly coupled mode: ppp-kinematic */
#define INSTC_DGPS     3                /* ins and gnss tightly coupled mode: DGPS/DGNSS */
#define INSTC_RTK      4                /* ins and gnss tightly coupled mode: RTK */

#define INSLC_SINGLE   1                /* ins and gnss loosely coupled mode: with single point position */
#define INSLC_PPK      2                /* ins and gnss loosely coupled mode: with ppk position */
#define INSLC_DGPS     3                /* ins and gnss loosely coupled mode: with dgps position */
#define INSLC_RTK      4                /* ins and gnss loosely coupled mode: with rtk position */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */
#define SOLF_ENU    2                   /* solution format: e/n/u-baseline */
#define SOLF_NMEA   3                   /* solution format: NMEA-183 */
#define SOLF_STAT   4                   /* solution format: solution status */
#define SOLF_GSIF   5                   /* solution format: GSI F1/F2 */
#define SOLF_INS    6                   /* solution format: output ins states */
#define SOLF_VO     7                   /* solution format: output vo states */
#define SOLF_GRTH   8                   /* solution format: ground truth */

#define SOLQ_NONE    0                  /* solution status: no solution */
#define SOLQ_FIX     1                  /* solution status: fix */
#define SOLQ_FLOAT   2                  /* solution status: float */
#define SOLQ_SBAS    3                  /* solution status: SBAS */
#define SOLQ_DGPS    4                  /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE  5                  /* solution status: single */
#define SOLQ_PPP     6                  /* solution status: PPP */
#define SOLQ_DR      7                  /* solution status: dead reconing */
#define SOLQ_DOP     8                  /* solution status: doppler measurement aid */
#define SOLQ_INHERIT 9                  /* solution status: ambiguity inherit fix status */
#define SOLQ_ROUND   10                 /* solution status: ambiguity round fix */
#define SOLQ_VO      11                 /* solution status: visual odometry status aid */
#define SOLQ_GRTH    12                 /* solution status: ground truth status */
#define SOLQ_WAAS    13                 /* solution status: solution calculated using corrections from an WAAS */
#define SOLQ_PROP    14                 /* solution status: propagated by a kalman filter without new observations */
#define SOLQ_OMIN    15                 /* solution status: OmniSTAR VBS position */
#define SOLQ_MONO    16                 /* solution status: monocular camera pose for pose measurement */
#define SOLQ_STEREO  17                 /* solution status: stereo camera pose for pose measurement */
#define MAXSOLQ      17                 /* max number of solution status */

#define TIMES_GPST  0                   /* time system: gps time */
#define TIMES_UTC   1                   /* time system: utc */
#define TIMES_JST   2                   /* time system: jst */

#define IONOOPT_OFF 0                   /* ionosphere option: correction off */
#define IONOOPT_BRDC 1                  /* ionosphere option: broadcast model */
#define IONOOPT_SBAS 2                  /* ionosphere option: SBAS model */
#define IONOOPT_IFLC 3                  /* ionosphere option: L1/L2 or L1/L5 iono-free LC */
#define IONOOPT_EST 4                   /* ionosphere option: estimation */
#define IONOOPT_TEC 5                   /* ionosphere option: IONEX TEC model */
#define IONOOPT_QZS 6                   /* ionosphere option: QZSS broadcast model */
#define IONOOPT_LEX 7                   /* ionosphere option: QZSS LEX ionospehre */
#define IONOOPT_STEC 8                  /* ionosphere option: SLANT TEC model */

#define TROPOPT_OFF 0                   /* troposphere option: correction off */
#define TROPOPT_SAAS 1                  /* troposphere option: Saastamoinen model */
#define TROPOPT_SBAS 2                  /* troposphere option: SBAS model */
#define TROPOPT_EST 3                   /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG 4                  /* troposphere option: ZTD+grad estimation */
#define TROPOPT_ZTD 5                   /* troposphere option: ZTD correction */

#define EPHOPT_BRDC 0                   /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC 1                   /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS 2                   /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC 3                 /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM 4                 /* ephemeris option: broadcast + SSR_COM */
#define EPHOPT_LEX  5                   /* ephemeris option: QZSS LEX ephemeris */

#define ARMODE_OFF  0                   /* AR mode: off */
#define ARMODE_CONT 1                   /* AR mode: continuous */
#define ARMODE_INST 2                   /* AR mode: instantaneous */
#define ARMODE_FIXHOLD 3                /* AR mode: fix and hold */
#define ARMODE_PPPAR 4                  /* AR mode: PPP-AR */
#define ARMODE_PPPAR_ILS 5              /* AR mode: PPP-AR ILS */
#define ARMODE_WLNL 6                   /* AR mode: wide lane/narrow lane */
#define ARMODE_TCAR 7                   /* AR mode: triple carrier ar */
#define ARMODE_WLNLC 8                  /* AR mode: wide lane/narrow lane */
#define ARMODE_TCARC 9                  /* AR mode: triple carrier ar */

#define SBSOPT_LCORR 1                  /* SBAS option: long term correction */
#define SBSOPT_FCORR 2                  /* SBAS option: fast correction */
#define SBSOPT_ICORR 4                  /* SBAS option: ionosphere correction */
#define SBSOPT_RANGE 8                  /* SBAS option: ranging */

#define POSOPT_POS   0                  /* pos option: LLH/XYZ */
#define POSOPT_SINGLE 1                 /* pos option: average of single pos */
#define POSOPT_FILE  2                  /* pos option: read from pos file */
#define POSOPT_RINEX 3                  /* pos option: rinex header pos */
#define POSOPT_RTCM  4                  /* pos option: rtcm station pos */
#define POSOPT_RAW   5                  /* pos option: raw station pos */

#define STR_NONE     0                  /* stream type: none */
#define STR_SERIAL   1                  /* stream type: serial */
#define STR_FILE     2                  /* stream type: file */
#define STR_TCPSVR   3                  /* stream type: TCP server */
#define STR_TCPCLI   4                  /* stream type: TCP client */
#define STR_NTRIPSVR 6                  /* stream type: NTRIP server */
#define STR_NTRIPCLI 7                  /* stream type: NTRIP client */
#define STR_FTP      8                  /* stream type: ftp */
#define STR_HTTP     9                  /* stream type: http */
#define STR_NTRIPC_S 10                 /* stream type: NTRIP caster server */
#define STR_NTRIPC_C 11                 /* stream type: NTRIP caster client */
#define STR_UDPSVR   12                 /* stream type: UDP server */
#define STR_UDPCLI   13                 /* stream type: UDP server */
#define STR_MEMBUF   14                 /* stream type: memory buffer */

#define STRFMT_RTCM2 0                  /* stream format: RTCM 2 */
#define STRFMT_RTCM3 1                  /* stream format: RTCM 3 */
#define STRFMT_OEM4  2                  /* stream format: NovAtel OEMV/4 */
#define STRFMT_OEM3  3                  /* stream format: NovAtel OEM3 */
#define STRFMT_UBX   4                  /* stream format: u-blox LEA-*T */
#define STRFMT_SS2   5                  /* stream format: NovAtel Superstar II */
#define STRFMT_CRES  6                  /* stream format: Hemisphere */
#define STRFMT_STQ   7                  /* stream format: SkyTraq S1315F */
#define STRFMT_GW10  8                  /* stream format: Furuno GW10 */
#define STRFMT_JAVAD 9                  /* stream format: JAVAD GRIL/GREIS */
#define STRFMT_NVS   10                 /* stream format: NVS NVC08C */
#define STRFMT_BINEX 11                 /* stream format: BINEX */
#define STRFMT_RT17  12                 /* stream format: Trimble RT17 */
#define STRFMT_SEPT  13                 /* stream format: Septentrio */
#define STRFMT_CMR   14                 /* stream format: CMR/CMR+ */
#define STRFMT_TERSUS 15                /* stream format: TERSUS */
#define STRFMT_LEXR  16                 /* stream format: Furuno LPY-10000 */
#define STRFMT_RINEX 17                 /* stream format: RINEX */
#define STRFMT_SP3   18                 /* stream format: SP3 */
#define STRFMT_RNXCLK 19                /* stream format: RINEX CLK */
#define STRFMT_SBAS  20                 /* stream format: SBAS messages */
#define STRFMT_NMEA  21                 /* stream format: NMEA 0183 */
#define STRFMT_GSOF  22                 /* stream format: GSOF */
#define STRFMT_UBXM8 23                 /* stream format: ublox evk-m8u for input imu data */
#define STRFMT_UBXSOL 24                /* stream format: ublox evk-m8u for input PVT solution data */
#define STRFMT_M39   25                 /* stream format: M39 imu raw data */
#define STRFMT_RINEX_RT  26             /* stream format: real-time rinex data input stream  */
#define STRFMT_M39MIX    27             /* stream format: M39 mix raw data (include image data) */
#define STRFMT_EUROCIMU  28             /* stream format: EUROC imu measurement data */
#define STRFMT_EUROCIMG  29             /* stream format: EUROC image measurement data */
#define STRFMT_KARLIMG   30             /* stream format: Karlsruhe dataset image measurement data */
#define STRFMT_MAGGNSS   31             /* stream format: malaga dataset GNSS position measurement data */
#define STRFMT_MAGIMU    32             /* stream format: malaga dataset imu measurement data */
#define STRFMT_MAGIMG    33             /* stream format: malaga dataset image measurement data */
#define STRFMT_OEM6_SOL  34             /* stream format: NovAtel OEM6 solution data */
#define STRFMT_OEM6_POSE 35             /* stream format: NovAtel OEM6 dual ant. pose measurement */
#define STRFMT_OEM6_RAW  36             /* stream format: NovAtel OEM6 raw observation data */
#define STRFMT_IGVSIM_FEAT 37           /* stream format: ins-gnss-vo multisensor simulator feature point measurement data */
#define STRFMT_IGVSIM_IMU  38           /* stream format: ins-gnss-vo multisensor simulator imu measurement data  */
#define STRFMT_IGVSIM_GNSS 39           /* stream format: ins-gnss-vo multisensor simulator gnss measurement data*/
#define STRFMT_IGVSIM_FEATALL 40        /* stream format: read all ins-gnss-vo multisensor simulator feature point measurement data */
#define STRFMT_STIM300        41        /* stream format: STIM300 imu raw data */

#define GROUND_TRUTH_KARL  1            /* Karlsruhe dataset ground truth solution format */
#define GROUND_TRUTH_EUROC 2            /* EuRoC MAV dataset ground truth solution format */
#define GROUND_TRUTH_MALA  3            /* Malaga Urban dataset ground truth solution format */

#ifndef EXTLEX
#define MAXRCVFMT    15                 /* max number of receiver format */
#else
#define MAXRCVFMT    16
#endif

#define STR_MODE_R  0x1                 /* stream mode: read */
#define STR_MODE_W  0x2                 /* stream mode: write */
#define STR_MODE_RW 0x3                 /* stream mode: read/write */

#define GEOID_EMBEDDED    0             /* geoid model: embedded geoid */
#define GEOID_EGM96_M150  1             /* geoid model: EGM96 15x15" */
#define GEOID_EGM2008_M25 2             /* geoid model: EGM2008 2.5x2.5" */
#define GEOID_EGM2008_M10 3             /* geoid model: EGM2008 1.0x1.0" */
#define GEOID_GSI2000_M15 4             /* geoid model: GSI geoid 2000 1.0x1.5" */
#define GEOID_RAF09       5             /* geoid model: IGN RAF09 for France 1.5"x2" */

#define COMMENTH    "%"                 /* comment line indicator for solution */
#define MSG_DISCONN "$_DISCONNECT\r\n"  /* disconnect message */

#define LLI_SLIP    0x01                /* LLI: cycle-slip */
#define LLI_HALFC   0x02                /* LLI: half-cycle not resovled */
#define LLI_BOCTRK  0x04                /* LLI: boc tracking of mboc signal */
#define LLI_HALFA   0x40                /* LLI: half-cycle added */
#define LLI_HALFS   0x80                /* LLI: half-cycle subtracted */

#define IMUFMT_KVH    1                 /* imu data format KVH */
#define IMUFMT_GI310  2                 /* imu data format GI310 */
#define IMUFMT_UBX    3                 /* imu data format ublox: EVK-M8U-0-00 */

#define INSUPD_TIME    0                /* only ins mechanization to updates ins states and updates states and covariance*/
#define INSUPD_MEAS    1                /* use gnss measurements to updates ins states */
#define INSUPD_INSS    2                /* only ins mechanization to updates ins states */

#define INSS_NONE      0                /* ins updates status: no update */
#define INSS_MECH      1                /* ins updates status: ins mechanization */
#define INSS_TIME      2                /* ins updates status: ins mechanization and propagate states and covariance */
#define INSS_LCUD      3                /* ins updates status: ins-gnss loosely-coupled updates  */
#define INSS_ZVU       4                /* ins updates status: ins zero velocity update */
#define INSS_ZARU      5                /* ins updates status: ins zero angular rate update */
#define INSS_NHC       6                /* ins updates status: ins non-holonomic constraint updates */
#define INSS_ODO       7                /* ins updates status: use odometry velocity measurement to aid ins */
#define INSS_TCUD      8                /* ins updates status: ins/gnss tightly coupled */
#define INSS_REBOOT    9                /* ins updates status: reboot */
#define INSS_DOPP      10               /* ins updates status: doppler measurement aid */
#define INSS_LACK      11               /* ins updates status: lack satellite tightly-coupled */
#define INSS_INIT      12               /* ins updates status: initialization ins states */
#define INSS_REINIT    13               /* ins updates status: re-initialization ins states */
#define INSS_VO        14               /* ins updates status: ins and visual odometry loosely coupled */
#define INSS_POSE      15               /* ins updates status: pose measurement fusion */
#define INSS_MAGH      16               /* ins updates status: magnetic heading auxiliary */
#define INSS_RTS       17               /* ins updates status: RTS smoother */
#define INSS_DEGRADE   18               /* ins updates status: degrade solution in forward/backward combined */
#define INSS_FBCOMB    19               /* ins updates status: forward/backward combined solution */

#define UPDINT_IMU     1                /* ins updates states time internal by imu data in ins-gnss coupled */
#define UPDINT_GNSS    2                /* ins updates states time internal by gnss data in ins-gnss coupled */
#define UPD_IN_EULER   0                /* update attitude in euler angles space when coupled measurement */

#define INS_BAEST      1                /* estimate accl bias in ins/gnss coupled */
#define INS_BAFIX      2                /* fix accl bias,and no estimate in ins/gnss coupled */
#define INS_BGEST      1                /* estimate gyro bias in ins/gnss coupled */
#define INS_BGFIX      2                /* fix gyro bias,and no estimate in ins/gnss coupled */
#define INS_DTEST      1                /* estimate time synchronization error between ins and gnss measurements for coupling */
#define INS_SGEST      1                /* estimate residual scale factors of gyroscopes */
#define INS_SGFIX      2                /* fix residual scale factors of gyroscopes */
#define INS_SAEST      1                /* estimate residual scale factors of accl */
#define INS_SAFIX      2                /* fix residual scale factors of accl */
#define INS_RAEST      1                /* estimate non-orthogonal between sensor axes for accl */
#define INS_RAFIX      2                /* fix non-orthogonal between sensor axes for accl */
#define INS_RGEST      1                /* estimate non-orthogonal between sensor axes for gyro. */
#define INS_RGFIX      2                /* fix non-orthogonal between sensor axes for gyro */
#define INS_LAEST      1                /* estimate lever arm for body to ant. */
#define INS_LAFIX      2                /* fix lever arm for body to ant. */

#define INS_RANDOM_WALK    1            /* stochastic process settings: random walk */
#define INS_RANDOM_CONS    2            /* stochastic process settings: random const */
#define INS_GAUSS_MARKOV   3            /* stochastic process settings: gauss-markov */

#define IMUDECFMT_RATE     1            /* imu measurements data format: angular rate/acceleration */
#define IMUDECFMT_INCR     2            /* imu measurements data format: angular/velocity increment */

#define IGCOM_USEGSOF      1            /* ins/gnss loosely couple using gsof message data */
#define IGCOM_USEOBS       2            /* ins/gnss loosely conple using code/phase observation data */
#define IGCOM_USESOL       3            /* ins/gnss loosely couple using pvt solution data */

#define INSALIGN_CORSE     1            /* ins initial corse alignment */
#define INSALIGN_FINE      2            /* ins initial fine alignment included fn and vn measurements */
#define INSALIGN_FINEEX    3            /* ins initial fine alignment extend method */
#define INSALIGN_LARGE     4            /* ins initial fine alignment with large yaw misalignment angle */
#define INSALIGN_DEFAULT   5            /* ins initial alignment defalut method */
#define INSALIGN_VELMATCH  6            /* ins initial velocity match alignment */
#define INSALIGN_SINGLE    7            /* ins initial use single point position */
#define INSALIGN_PPK       8            /* ins initial use ppk position */
#define INSALIGN_DGPS      9            /* ins initial use dgps position */
#define INSALIGN_RTK       10           /* ins initial use rtk position */
#define INSALIGN_OFF       11           /* ins initial no alignment */

#define IMUCOOR_FRD        1            /* imu body coordinate frame: frd */
#define IMUCOOR_RFU        2            /* imu body coordinate frame: rfu */

#define IMUVALFMT_DEG      1            /* imu gyro measurement data value format: degree */
#define IMUVALFMT_RAD      2            /* imu gyro measurement data value format: rad */

#define IMUDETST_GLRT      1            /* runs the generalized likelihood test for detect static imu measurement */
#define IMUDETST_MV        2            /* runs the acceleration moving variance detector */
#define IMUDETST_MAG       3            /* runs the acceleration magnitude detector */
#define IMUDETST_ARE       4            /* runs the angular rate energy detector */
#define IMUDETST_ALL       5            /* runs all static detector */

#define NPOS               3            /* # of precious epochs gps antenna position */
#define BASE_Q_EMPTY       0
#define BASE_Q_NO_EMPTY    1

#define VO_MONO            1            /* mono camera visual odometry mode */
#define VO_STEREO          2            /* stereo camera visual odometry mode */
#define DESCR_MAXLEN       128          /* max length of feature descriptor */

#define POSE_CAM           1            /* fusion external pose measurement from camera visual odometry */
#define POSE_DUAL_ANT      2            /* fusion external pose measurement from dual antennas */

#define COPLED_RB          0            /* ins/gnss coupled with robocentric formulation */

#define P2_5        0.03125               /* 2^-5 */
#define P2_6        0.015625              /* 2^-6 */
#define P2_10       9.765625000000000E-04 /* 2^-10 */
#define P2_11       4.882812500000000E-04 /* 2^-11 */
#define P2_12       2.441406250000000E-04 /* 2^-12 */
#define P2_15       3.051757812500000E-05 /* 2^-15 */
#define P2_17       7.629394531250000E-06 /* 2^-17 */
#define P2_19       1.907348632812500E-06 /* 2^-19 */
#define P2_20       9.536743164062500E-07 /* 2^-20 */
#define P2_21       4.768371582031250E-07 /* 2^-21 */
#define P2_23       1.192092895507810E-07 /* 2^-23 */
#define P2_24       5.960464477539063E-08 /* 2^-24 */
#define P2_27       7.450580596923828E-09 /* 2^-27 */
#define P2_29       1.862645149230957E-09 /* 2^-29 */
#define P2_30       9.313225746154785E-10 /* 2^-30 */
#define P2_31       4.656612873077393E-10 /* 2^-31 */
#define P2_32       2.328306436538696E-10 /* 2^-32 */
#define P2_33       1.164153218269348E-10 /* 2^-33 */
#define P2_35       2.910383045673370E-11 /* 2^-35 */
#define P2_38       3.637978807091710E-12 /* 2^-38 */
#define P2_39       1.818989403545856E-12 /* 2^-39 */
#define P2_40       9.094947017729280E-13 /* 2^-40 */
#define P2_43       1.136868377216160E-13 /* 2^-43 */
#define P2_48       3.552713678800501E-15 /* 2^-48 */
#define P2_50       8.881784197001252E-16 /* 2^-50 */
#define P2_55       2.775557561562891E-17 /* 2^-55 */

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define MAX(x,y)    ((x)>=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#ifdef WIN32
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#define FILEPATHSEP '\\'
#else
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define lock(f)     pthread_mutex_lock(f)
#define unlock(f)   pthread_mutex_unlock(f)
#define FILEPATHSEP '/'
#endif

/* FEATURE_FWD_MATCH <BR> FEATURE_BCK_MATCH <BR> FEATURE_MDL_MATCH ----------*/
enum feature_match_type {
    FEATURE_FWD_MATCH,
    FEATURE_BCK_MATCH,
};
/* type definitions ----------------------------------------------------------*/
typedef struct {                /* time struct */
    time_t time;                /* time (s) expressed by standard time_t */
    double sec;                 /* fraction of second under 1 s */
} gtime_t;

typedef struct {                /* camera model struct type */
    double k1,k2,k3,k4,k5,k6;   /* distortion parameter */
    double p1,p2,s1,s2,s3,s4;
    double K[9];                /* intrinsic parameter matrix */
} cam_t;

typedef struct feature {        /* feature data type */
    gtime_t time;               /* feature point extract time */
    double u,v;                 /* pixel u,v coord */
    int valid,uid;              /* valid flag and uid of this feature point */
    int descrlen;               /* length of descriptor */
    double descr[DESCR_MAXLEN]; /* descriptor */
    UT_hash_handle hh;          /* makes this structure hashable */
} feat_data_t;

typedef struct  {               /* image data type */
    gtime_t time;               /* data time */
    int w,h,id;                 /* width and height of image/frame id of image */
    unsigned char *data;        /* image data buffer */
    feature *feat;              /* feature point measurement data */
    UT_hash_handle hh;          /* makes this structure hashable */
} img_t;

typedef struct trackd {         /* feature points track record */
    gtime_t ts,te;              /* start/end timestamp */
    int n,nmax;                 /* number and max number of feature points */
    int last_idx,first_frame,last_frame;
                                /* last feature point index/first frame index/last frame index of this track */
    unsigned char flag;         /*
                                 * 0: updated track                  --TRACK_UPDATED
                                 * 1: new track,                     --TRACK_NEW
                                 * 2: track lost,                    --TRACK_LOST
                                 * 3: bad track,                     --TRACK_BAD
                                 * 4: have used to filter            --TRACK_FILTER
                                 * 5: initial feature position ok    --TRACK_INIT_POS_OK
                                 * 6: initial feature position fail  --TRACK_INIT_POS_FAIL
                                 * */
    int uid;                    /* a unique identifier of this track */
    struct feature *data;       /* track feature data */
    double xyz[3];              /* feature position in ecef */
} trackd_t;

typedef struct track {          /* store all feature track data type */
    int n,nmax;                 /* number and max number of track data */
    int newtrack[MAXBUFF],updtrack[MAXBUFF],losttrack[MAXBUFF],exceedtrack[MAXBUFF],nnew,nupd,nlost,nexceed;
    int nless,lesstrack[MAXBUFF];
                                /* new/updated/lost/exceed-max-track-length feature track index in `.data' */
    trackd_t *data;             /* track data */
} track_t;

typedef struct {                /* bucketing parameters */
    int nmax;                   /* maximal number of features per bucket */
    double w;                   /* width of bucket */
    double h;                   /* height of bucket */
} bucketopt_t;

typedef struct {                /* camera parameters (all are mandatory/need to be supplied) for matching feature points */
    double f,fu,fv;             /* focal length (in pixels) */
    double cu;                  /* principal point (u-coordinate) */
    double cv;                  /* principal point (v-coordinate) */
} calib_t;

typedef struct {                /* visual odometry matching options */
    int img_w,img_h;            /* image width and height in pixel */
    int nms_n;                  /* non-max-suppression: min. distance between maxima (in pixels) */
    int nms_tau;                /* non-max-suppression: interest point peakiness threshold */
    int match_binsize;          /* matching bin width/height (affects efficiency only) */
    int match_radius;           /* matching radius (du/dv in pixels) */
    int match_disp_tol;         /* dv tolerance for stereo matches (in pixels) */
    int outlier_disp_tol;       /* outlier removal: disparity tolerance (in pixels) */
    int outlier_flow_tol;       /* outlier removal: flow tolerance (in pixels) */
    int multi_stage;            /* 0=disabled,1=multistage matching (denser and faster) */
    int half_res;               /* 0=disabled,1=match at half resolution, refine at full resolution */
    int refine;                 /* refinement (0=none,1=pixel,2=sub-pixel) */
    double f,fu,fv,cu,cv,base;  /* calibration parameters (only for match prediction) */
    int roi[2][2];              /* match image ROI options:
                                 * roi[0]: (left-top),
                                 * roi[1]: (right-bottom)
                                 * */
    bucketopt_t bucket;         /* bucketing parameters */
} matchopt_t;

typedef struct match_point {    /* match feature point type */
    float up,vp;                /* u,v-coordinates in previous image */
    float uc,vc;                /* u,v-coordinates in current image */
    int ip,ic;                  /* feature id/feature index in precious/current image(for tracking) */
    int kltstat;                /* KLT track status (if KLT enable) */
} match_point_t;

typedef struct match_set {      /* store all matched feature points */
    int n,nmax;                 /* max number and number of matched points */
    match_point_t *data;        /* matched points data */
} match_set_t;

typedef struct {                /* match struct type */
    gtime_t time,pt;            /* current/precious match time */
    img_t Ip,Ic;                /* precious/current image data */
    match_set_t mp_dense;       /* dense macthed feature points data */
    match_set_t mp_sparse;      /* sparse macthed feature points data */
    match_set_t mp_bucket;      /* bucket matched feature points data */
    matchopt_t opt;             /* matched options */
} match_t;

typedef struct {                /* observation data record */
    gtime_t time;               /* receiver sampling time (GPST) */
    unsigned char sat,rcv;      /* satellite/receiver number */
    unsigned char SNR [NFREQ+NEXOBS]; /* signal strength (0.25 dBHz) */
    unsigned char LLI [NFREQ+NEXOBS]; /* loss of lock indicator */
    unsigned char code[NFREQ+NEXOBS]; /* code indicator (CODE_???) */
    double L[NFREQ+NEXOBS]; /* observation data carrier-phase (cycle) */
    double P[NFREQ+NEXOBS]; /* observation data pseudorange (m) */
    float  D[NFREQ+NEXOBS]; /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {            /* magnetometer measurement data type */
    gtime_t time;           /* magnetometer measurement data time */
    double val[3];          /* magnetometer measurement data expressed in body frame */
} mag_t;

typedef union {             /* generic 3d vector */
    struct {
        double x,y,z;
    };
    double vec[3];
} vec3_t;

typedef union {             /* quaternion type */
    struct {
        double q0,q1,q2,q3;
    };
    struct {
        double w,x,y,z;
    };
    double vec[4];
} quat_t;

typedef union {             /* euler angle */
    struct {
        double yaw;
        double pitch;
        double roll;
    };
    double vec[3];
} euler_t;

typedef struct {            /* odometry measurement record type */
    gtime_t time;           /* measurement time */
    double dt;              /* time difference related to increment distance */
    double dr;              /* increment of distance (m) */
    double vr[3];           /* wheel velocity in vehicle rear frame */
} odod_t;

typedef struct {            /* odometry measurement data */
    int n,nmax;             /* number/max number of imu data record */
    odod_t *data;           /* odometry measurement data records */
} odo_t;

typedef struct {            /* imu data record type */
    gtime_t time;           /* time */
    double gyro[3];         /* angular rate measurements in b-frame (rad/s) */
    double accl[3];         /* force measurements in b-frame (m/s^2) */
    double temp;            /* temperature (C) */
    int stat;               /* sensor status (0:ok) */

    unsigned int pps;       /* imu pps */
    unsigned int imuc;      /* imu count */

    short int odoc;         /* odometer count */
    odod_t odo;             /* odometry measurement data */
} imud_t;

typedef struct {            /* imu data type */
    imud_t *data;           /* imu data records */
    int n,nmax;             /* number/max number of imu data record */
    int decfmt;             /* imu measurement data decode format */
    int format;             /* imu format (IMUFORMAT_???) */
    int coor;               /* imu body coordinate frame (IMUCOOR_???) */
    int valfmt;             /* imu gyro measurement data format (IMUVALFMT_???) */
} imu_t;

typedef struct {            /* m39/image time tag */
    int sowc;               /* counts of sow increments */
    double sow,week;        /* GPS sow (s)/GPS week */
    gtime_t time,zda;       /* frame time */
    timespec pps,ppsp,fts;  /* pps/frame time */
} m39_mix_t;

typedef struct {            /* gnss positioning result for ins-gnss loosely coupled */
    gtime_t t;              /* measurements time in GPST (s)*/
    unsigned char ns;       /* number of valid satellites */
    unsigned char stat;     /* status of PVT solution */
    double pe[3];           /* position of gnss antenna in ecef-frame (m) */
    double ve[3];           /* velocity of gnss antenna in ecef-frame (m/s) */
    double std[6];          /* position/velocity variance from gnss positioning {m,m/s} {xx,yy,zz} */
    double covp[9],covv[9]; /* position/velocity covariance {m^2,m/s^2} {xx,yy,zz,xy,xz,yz}*/
} gmea_t;

typedef struct {            /* pose measurement data struct */
    gtime_t time;           /* time of pose measurement data */
    int type;               /* pose measurement data type (POSE_FUSION_???) */
    int stat;               /* pose measurement data status */
    double rpy[3],len;      /* pose measurement data expressed in euler angle (Cm=R(x)*R(y)*R(z)) and baseline length (m) */
    double C[9];            /* pose measurement data expressed in dcm format */
    double var[3];          /* pose measurement variance expressed in tangent space or euler angle */
} pose_meas_t;

typedef struct {            /* pose measurement data buffer type */
    int n,nmax;             /* number of solution/max number of buffer */
    int cyclic;             /* cyclic buffer flag */
    int start,end;          /* start/end index */
    double dt;              /* time difference between solutions */
    gtime_t time;           /* current solution time */
    pose_meas_t *data;      /* solution data */
} posebuf_t;

typedef struct {            /* gnss position/velocity measurement */
    int n,nmax;             /* number and max number of position/velocity measurements */
    gmea_t *data;           /* gnss position/velocity measurement data */
} gmeas_t;

typedef struct {               /* camera states type */
    gtime_t time;              /* camera states time */
    int status;                /* monocular visual odometry estimate status */
    double rc[3],Cce[9],rc0[3],Cce0[9];
                               /* current camera position/attitude in e-frame */
    double Cbc[9],lbc[3];      /* dcm and lever arm of camera to imu body frame */
    double dT[16],ratio;       /* camera transformation matrix from precious to current */
    double T[16];              /* camera transformation matrix in current epoch */
    double phi[3],P[9];        /* covariance of camera pose expressed in tangent space */
} vostate_t;

typedef struct {             /* ins states type */
    gtime_t time,ptime;      /* ins states time and precious time of ins states*/
    gtime_t plct,ptct;       /* precious time of ins-gnss loosely/tightly coupled (0: no coupled)*/
    double dt;               /* time difference of precious ins updates and current time */
    double Cbe[9];           /* b-frame to e-frame trans matrix */
    double re[3],ve[3],ae[3],oxyz[3];
                             /* position/velocity/acceleration (e-frame) */
    double dtr[6],dtrr;      /* receiver clock bias to time systems (s) and clock‐drift (m/s) */
    
    double rn[3],vn[3],an[3];/* position/velocity/acceleration (n-frame) */
    double Cbn[9],dvn[3];    /* transform matrix of b-frame to n-frame/velocity increment at the precious epoch */

    double ba[3],bg[3];     /* accelemeter-bias and gyro-bias */
    double Ma[9];           /* non-orthogonal between sensor axes and body frame for accelerometers
                             *    | sax    ra_xy  ra_xz |      | sgx    rg_xy  rg_xz |
                             * Ma=| ra_yx  say    ra_yz |   Mg=| rg_yx  sgy    rg_yz |
                             *    | ra_zx  ra_zy  saz   |      | rg_zx  rg_zy  sgz   |
                             * */
    double Mg[9];           /* non-orthogonal between sensor axes and body frame for gyroscopes */
    double Gg[9];           /* g-dependent bias for a gyro triad */
    double fb0[3],omgb0[3]; /* uncorrected specific-force (b-frame)/angular rate (b-frame) */
    double fb[3],omgb[3];   /* corrected specific-force (b-frame)/angular rate (b-frame) */

    double fbp[3],omgbp[3]; /* precious epoch corrected specific-force (b-frame)/angular rate (b-frame) */
    double lever[3];        /* lever arm for body to ant. (m) */
    double lbc[3],Cbc[9];   /* lever arm/misalignment for body to camera (for additional) */

    double Cbr[9];          /* transform matrix from odometry rear frame to body frame */
    double os;              /* odometry scale factor */
    double rbl[3];          /* lever arm of odometry rear frame to body frame */
    double vr[3];           /* odometry velocity in rear frame */

    double fx,fy,ox,oy,k1,k2,p1,p2;
                            /* camera calibration parameters */
    double Cvb[9],len;      /* misalignment from v-frame to b-frame (defined at dual ant.) and length of dual ant. */
    double dopv[3];         /* doppler velocity (ecef, m/s) */
    double pnull[3],vnull[3],Cnull[9];
                            /* these three variables should have the same physical
                             * interpretation with `orientation`, `position`, and
                             * `velocity` in `ecef'. there three variables are used to modify
                             * the transition matrices to make the observability matrix
                             * have proper null space.
                             * */
    double *x,*P;           /* ekf estimated states/covariance matrix */
    double *xa,*Pa;         /* estimated states and covariance for ins-gnss loosely coupled */
    double *xb,*Pb;         /* fixed states and covariance (except phase bias) */
    double *P0,*F;          /* predict error states correction and its covariance matrix/transmit matrix */

    double pins[9],pCbe[9]; /* ins states (position/velocity/acceleration/attitude) of precious epoch in ecef-frame */
    gmeas_t gmeas;          /* gps position/velocity measurements */

    double age,ratio;       /* age of differential of ins and gnss (s)/ambiguity fix ratio */
    int stat,gstat,pose;    /* ins updates stat,gnss updates status and pose fusion status */
    int ns;                 /* numbers valid satellite for loosely coupled */
    int nx,nb;              /* numbers of estimated states,fixed states (except phase bias) */
    void *rtkp;             /* pointer rtk struct data */
    vostate_t vo;           /* vosidual odometry states */
} insstate_t;

typedef struct {            /* PSD for ins-gnss loosely coupled ekf states */
    double gyro;            /* gyro noise PSD (rad^2/s) */
    double accl;            /* accelerometer noise PSD (m^2 s^-3) */
    double ba;              /* accelerometer bias random walk PSD (m^2 s^-5) */
    double bg;              /* gyro bias random walk PSD (rad^2 s^-3) */
    double clk;             /* receiver clock phase-drift PSD (m^2/s) */
    double clkr;            /* receiver clock drift PSD */
} psd_t;

typedef struct {            /* initial uncertainty for ins-gnss loosely coupled */
    double att;             /* initial attitude uncertainty per axis (rad) (ecef) */
    double vel;             /* initial velocity uncertainty per axis (m/s) (ecef) */
    double pos;             /* initial position uncertainty per axis (m) (ecef) */
    double ba;              /* initial accl bias uncertainty (m/s^2) */
    double bg;              /* initial gyro bias uncertainty (rad/s) */
    double dt;              /* initial time synchronization error uncertainty (m^2) */
    double sg;              /* initial residual scale factors of gyroscopes uncertainty */
    double sa;              /* initial residual scale factors of accl uncertainty */
    double ra;              /* initial non-orthogonal between sensor axes for accl uncertainty */
    double rg;              /* initial non-orthogonal between sensor axes for gyro uncertainty */
    double lever,mla;       /* initial lever arm for body to ant. uncertainty {m/rad} */
    double os;              /* initial odometry scale factor uncertainty */
    double ol;              /* initial odometry lever arm uncertainty */
    double oa;              /* initial odometry misalignment uncertainty */
    double rc;              /* initial receiver clock uncertainty (m) */
    double rr;              /* initial receiver clock drift uncertainty (s/s) */
    double cma,lma;         /* initial misalignment/lever-arm from camera to imu body uncertainty (rad/m) */
    double fo[4],kp[4];     /* initial camera calibration parameters({fx,fy,ox,oy},{k1,k2,p1,p2}) uncertainty ({pixel,1.0}) */
    double vma;             /* initial misalignment from v-frame to imu body uncertainty (rad) */
} unc_t;

typedef struct {            /* imu error model type */
    double bg[3];           /* gyro constant bias (rad/sec) */
    double ba[3];           /* accl constant bias (m/s^2) */
    double Ma[9],Mg[9];     /* non-orthogonal between sensor axes (ppm) */
    double Gg[9];           /* g-dependent bias for a gyro triad */
    double sR0G,TauG;       /* gyro correlated bias (TauG: time for Markov process (s),sROG: deg/h */
    double sROA,TauA;       /* accl correlated bias (TauA: time for Markov process (s),sROA: ug */
    double wgn[3];          /* angular random walk (rad s^-0.5) */
    double wan[3];          /* velocity random walk (m s^-1.5) */
    double wbg[3];          /* gyro bias random walk PSD (rad^2 s^-3) */
    double wba[3];          /* accelerometer bias random walk PSD (m^2 s^-5) */
} imu_err_t;

typedef struct {            /* ins navigation initial fine alignments type */
    double eb[3];           /* gyro constant bias (deg/h) */
    double db[3];           /* accl constant bias (ug) */
    double web[3];          /* angular random walk (deg/sqrt(h)) */
    double wdb[3];          /* velocity random walk (ug/sqrt(Hz)) */
    int ns;                 /* numbers of imu data for ins alignment */
    int chkstatic;          /* check static imu data when static alignment */
    double dt;              /* imu sampling interval */
    double phi0[3];         /* initial misalignment angles estimation */
    double wvn[3];          /* velocity measurement noise (3x1 vector) */
} ins_align_t;

typedef struct {            /* zero velocity detector setting */
    int ws;                 /* static detector window size */
    double mt;              /* min span time of static imu measurement data */
    double sp;              /* zero velocity detect time span from first epoch */
    double gthres;          /* threshold of static detector by gravity  */
    double athres[3];       /* threshold of static detector by accl measurement for checking its variance */
    double gyrothres[3];    /* threshold of static detector by gyro measurement */
    double odt;             /* time difference to convert odometry measurement to velocity */

    double sig_a,sig_g;     /* sig_a: standard deviation of the acceleromter noise (m/s^2). this is used to
                             * control the zero-velocity detectors trust in the accelerometer data.
                             * sig_g: standard deviation of the gyroscope noise (rad/s]). this is used to
                             * control the zero-velocity detectors trust in the gyroscope data */
    double gamma[4];        /* threshold used in the zero-velocity detector. If the test statistics are
                             * below this value the zero-velocity hypothesis is chosen,
                             * [0]: GLRT, [1]: MV, [2]: MAG, [3]: ARE */
} ins_zv_t;

typedef struct {            /* odometry options types */
    int dir;                /* odometry direction */
    int all;                /* 0: only use forward velocity,1: use side and forward velocity */
    double res;             /* resolution of odometry */
    double s;               /* scale factor of odometry */
    double d;               /* diameter of wheel */
    double lever[3];        /* lever arm from odometry to ins,expressed in body-frame */
    double odt;             /* time difference for convert odometry increment distance to velocity */
    double ostd;            /* odometry standard deviation for ins and odometry coupled */
} odopt_t;

typedef struct {            /* visual odometry aid ins options (here for rectified image) */
    matchopt_t match;       /* feature match options */
    calib_t calib;          /* camera parameters for matching feature points */
    cam_t   cam;            /* camera intrinsic parameters */
    double height;          /* camera height above ground (meters) */
    double inlier_thres;    /* fundamental matrix inlier threshold */
    double motion_thres;    /* directly return false on small motions */
    double hz;              /* camera frame update frequency */
    double lbc[3],ebc[3];   /* lever arm/euler angle of camera to imu (m/rad)
                             * Cbc=R(ebc[0])*R(ebc[1])*R(ebc[2])
                             * */
    int    ransac_iters;    /* number of RANSAC iterations */
} voopt_t;

typedef struct {            /* magnetometer options type */
    double sx,sy,sz;        /* magnetometer data scale factor */
    double ox,oy,oz;        /* magnetometer data offset */
} magopt_t;

typedef struct {            /* ins options type */
    int exinserr;           /* extend ins error model */
    int inithead;           /* use the roll, pitch and gyro measurements to initial head */
    int gravityex;          /* precious gravity model */
    int updint;             /* ins states and covariance propagate time internal (1: IMU,2:GNSS) */
    int baopt;              /* accelemeter-bias option */
    int bgopt;              /* gyro-bias option */
    int cnscl;              /* coning and sculling compensation */
    int estsg;              /* estimate residual scale factors of gyroscopes */
    int estsa;              /* estimate residual scale factors of accelemeter */
    int estdt;              /* estimate time synchronization error between ins and gnss */
    int estrg;              /* estimate non-orthogonal between sensor axes for gyro. */
    int estra;              /* estimate non-orthogonal between sensor axes for accl. */
    int estlever;           /* estimate lever arm for body to ant. */
    int estodos;            /* estimate odometry scale factor */
    int estodoa;            /* estimate odometry misalignment of body frame and rear frame */
    int estodol;            /* estimate odometry lever arm of rear frame to body frame */
    int estmisv;            /* estimate misalignment from v-frame (defined at dual ant.) to b-frame */
    int estcama;            /* estimate misalignment from b-frame to c-frame */
    int estcaml;            /* estimate lever arm from b-frame to c-frame */
    int estcam_fo;          /* estimate camera calibration parameters: fx,fy,ox,oy */
    int estcam_kp;          /* estimate camera calibration parameters: k1,k2,p1,p2 */

    int baproopt;           /* accl. bias stochastic process settings (INS_RANDOM_WALK,...) */
    int bgproopt;           /* gryo. bias stochastic process settings (INS_RANDOM_WALK,...) */
    int sgproopt;           /* residual scale factors of gyroscopes stochastic process settings (INS_RANDOM_WALK,...) */
    int saproopt;           /* residual scale factors of accl. stochastic process settings (INS_RANDOM_WALK,...) */
    int dtproopt;           /* ins-gnss time synchronization error stochastic process setting (INS_RANDOM_WALK) */
    int rgproopt;           /* non-orthogonal between sensor axes for gyro stochastic process setting */
    int raproopt;           /* non-orthogonal between sensor axes for accl stochastic process setting */
    int osproopt;           /* odomatry scale factor stochastic process option */
    int olproopt;           /* odometry lever arm stochastic process option */
    int oaproopt;           /* odometry misalignment stochastic process option */
    int cmaopt,vmaopt;      /* camera misalignment and v-frame misalignment to b-frame stochastic process option */
    int claopt;             /* camera lever arm from b-frame to c-frame stochastic process option */
    int cfoopt,ckpopt;      /* camera calibration parameters stochastic process option */

    int odo;                /* use odometry velocity measurement to aid ins */
    int pose_aid;           /* use pose measurement from dual ant. to aid ins navigation */

    int align_vn;           /* ins initial align uses kalman filter with vn as measurement in static-alignment mode */
    int align_fn;           /* ins initial align uses kalman filter with fn as measurement in static-alignment mode */
    int align_dualants;     /* ins initial use dual antennas pose measurement data */
    int align_given;        /* ins initial use given attitude/velocity/position */

    int imuformat;          /* imu measurement data format (IMUFORMAT_???) */
    int imudecfmt;          /* imu measurement decode format (IMUDECFMT_???) */
    int imucoors;           /* imu body coordinate frame (IMUCOOR_???) */
    int imuvalfmt;          /* imu measurement gyro data value format (IMUVALFMT_???) */
    int lcopt;              /* ins/gnss loosely coupled use gnss observation or gsof message */

    int exprn;              /* use extend method to get process noise matrix */
    int exphi;              /* use precise system propagate matrix for ekf */
    int exvm;               /* extend method to velocity matching for ins navigation initial */
    int iisu;               /* initial ins state use rtk options (SOLQ_???)*/

    int nhc;                /* non-holonomic constraint options */
    int zvu;                /* zero velocity update options */
    int zaru;               /* zero augular rate update options */
    int detst;              /* detect static imu measurement data */
    int magh;               /* magnetometer auxiliary */

    int lc,tc;                 /* ins-gnss loosely/tightly coupled mode (INSLC_???/INSTC_???) */
    int dopp;               /* use doppler measurement to aid ins updates states */
    int intpref;            /* time-interpolation for observation when tightly coupled */
    int minp;               /* min number position for ins alignment */
    int soltype;            /* solution type (0:forward,1:backward,2:combined,3:RTS) */
    int transmit_corr;      /* transmit error state correction (dx_t=(I+F*dt)dx_t_1) */

    gtime_t ext[16][2];     /* exclude time for processing ins measurement data,[0]: start time,[1]: end time */

    double lever[3];        /* lever arm terms of body-frame to gnss antenna */
    double mis_euler[3];    /* misalignment euler angle from v-frame (defined in dual ant.) to b-frame */
    double len;             /* length of dual antennas (m) */
    double hz;              /* imu measurement sampling frequency */
    double nhz;             /* non-holonomic constraint update frequency */
    double rn0[3],vn0[3],att0[3],t0;
                            /* ins initial position/velocity/attitude in n-frame {deg, m/s, rad}
                             * and week seconds from user giving
                             * */
    psd_t psd;              /* PSD for ins-gnss loosely coupled ekf states */
    unc_t unc;              /* initial uncertainty for ins-gnss loosely coupled */
    imu_err_t imuerr;       /* imu error model */
    ins_align_t align;      /* ins navigation initial alignment options */
    ins_zv_t zvopt;         /* ins zero velocity detector options */
    odopt_t odopt;          /* odometry options */
    voopt_t voopt;          /* visual odometry aid ins options */
    magopt_t magopt;        /* magnetometer options */
    void *gopt;             /* gnss-rtk options */
} insopt_t;

/* type definition -----------------------------------------------------------*/
typedef struct {                        /* signal index type */
    int n;                              /* number of index */
    int frq[MAXOBSTYPE];                /* signal frequency (1:L1,2:L2,...) */
    int pos[MAXOBSTYPE];                /* signal index in obs data (-1:no) */
    unsigned char pri [MAXOBSTYPE];     /* signal priority (15-0) */
    unsigned char type[MAXOBSTYPE];     /* type (0:C,1:L,2:D,3:S) */
    unsigned char code[MAXOBSTYPE];     /* obs code (CODE_L??) */
    double shift[MAXOBSTYPE];           /* phase shift (cycle) */
} sigind_t;

typedef struct {        /* observation data */
    int n,nmax;         /* number of observation data/allocated */
    sigind_t sind[2][7];/* signal index for rover and base station,0: rover,1: base
                         * sind[rcv][0]: GPS
                         * sind[rcv][1]: GLONASS
                         * sind[rcv][2]: Galileo
                         * sind[rcv][3]: QZSS
                         * sind[rcv][4]: SBAS
                         * sind[rcv][5]: BDS
                         * sind[rcv][6]: IRNS */
    obsd_t *data;       /* observation data records */
} obs_t;

typedef struct {        /* earth rotation parameter data type */
    double mjd;         /* mjd (days) */
    double xp,yp;       /* pole offset (rad) */
    double xpr,ypr;     /* pole offset rate (rad/day) */
    double ut1_utc;     /* ut1-utc (s) */
    double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
    int n,nmax;         /* number and max number of data */
    erpd_t *data;       /* earth rotation parameter data */
} erp_t;

typedef struct {        /* antenna parameter type */
    int sat;            /* satellite number (0:receiver) */
    char type[MAXANT];  /* antenna type */
    char code[MAXANT];  /* serial number or satellite code */
    gtime_t ts,te;      /* valid time start and end */
    double off[NFREQ][ 3]; /* phase center offset e/n/u or x/y/z (m) */
    double var[NFREQ][19]; /* phase center variation (m) */
                        /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} pcv_t;

typedef struct {        /* antenna parameters type */
    int n,nmax;         /* number of data/allocated */
    pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

typedef struct {        /* almanac type */
    int sat;            /* satellite number */
    int svh;            /* sv health (0:ok) */
    int svconf;         /* as and sv config */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    gtime_t toa;        /* Toa */
                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,OMGd;
    double toas;        /* Toa (s) in week */
    double f0,f1;       /* SV clock parameters (af0,af1) */
} alm_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode,iodc;      /* IODE,IODC */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    int code;           /* GPS/QZS: code on L2, GAL/CMP: data sources */
    int flag;           /* GPS/QZS: L2 P data flag, CMP: nav type */
    gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    double crc,crs,cuc,cus,cic,cis;
    double toes;        /* Toe (s) in week */
    double fit;         /* fit interval (h) */
    double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[4];      /* group delay parameters */
                        /* GPS/QZS:tgd[0]=TGD */
                        /* GAL    :tgd[0]=BGD E5a/E1,tgd[1]=BGD E5b/E1 */
                        /* CMP    :tgd[0]=BGD1,tgd[1]=BGD2 */
    double Adot,ndot;   /* Adot,ndot for CNAV */
} eph_t;

typedef struct {        /* GLONASS broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode;           /* IODE (0-6 bit of tb field) */
    int frq;            /* satellite frequency number */
    int svh,sva,age;    /* satellite health, accuracy, age of operation */
    gtime_t toe;        /* epoch of epherides (gpst) */
    gtime_t tof;        /* message frame time (gpst) */
    double pos[3];      /* satellite position (ecef) (m) */
    double vel[3];      /* satellite velocity (ecef) (m/s) */
    double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;   /* SV clock bias (s)/relative freq bias */
    double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

typedef struct {        /* precise ephemeris type */
    gtime_t time;       /* time (GPST) */
    int index;          /* ephemeris index for multiple files */
    double pos[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
    float  std[MAXSAT][4]; /* satellite position/clock std (m|s) */
    double vel[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
    float  vst[MAXSAT][4]; /* satellite velocity/clk-rate std (m/s|s/s) */
    float  cov[MAXSAT][3]; /* satellite position covariance (m^2) */
    float  vco[MAXSAT][3]; /* satellite velocity covariance (m^2) */
} peph_t;

typedef struct {        /* precise clock type */
    gtime_t time;       /* time (GPST) */
    int index;          /* clock index for multiple files */
    double clk[MAXSAT][1]; /* satellite clock (s) */
    float  std[MAXSAT][1]; /* satellite clock std (s) */
} pclk_t;

typedef struct {        /* SBAS ephemeris type */
    int sat;            /* satellite number */
    gtime_t t0;         /* reference epoch time (GPST) */
    gtime_t tof;        /* time of message frame (GPST) */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    double pos[3];      /* satellite position (m) (ecef) */
    double vel[3];      /* satellite velocity (m/s) (ecef) */
    double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
    double af0,af1;     /* satellite clock-offset/drift (s,s/s) */
} seph_t;

typedef struct {        /* norad two line element data type */
    char name [32];     /* common name */
    char alias[32];     /* alias name */
    char satno[16];     /* satellilte catalog number */
    char satclass;      /* classification */
    char desig[16];     /* international designator */
    gtime_t epoch;      /* element set epoch (UTC) */
    double ndot;        /* 1st derivative of mean motion */
    double nddot;       /* 2st derivative of mean motion */
    double bstar;       /* B* drag term */
    int etype;          /* element set type */
    int eleno;          /* element number */
    double inc;         /* orbit inclination (deg) */
    double OMG;         /* right ascension of ascending node (deg) */
    double ecc;         /* eccentricity */
    double omg;         /* argument of perigee (deg) */
    double M;           /* mean anomaly (deg) */
    double n;           /* mean motion (rev/day) */
    int rev;            /* revolution number at epoch */
} tled_t;

typedef struct {        /* norad two line element type */
    int n,nmax;         /* number/max number of two line element data */
    tled_t *data;       /* norad two line element data */
} tle_t;

typedef struct {        /* TEC grid type */
    gtime_t time;       /* epoch time (GPST) */
    int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
    double rb;          /* earth radius (km) */
    double lats[3];     /* latitude start/interval (deg) */
    double lons[3];     /* longitude start/interval (deg) */
    double hgts[3];     /* heights start/interval (km) */
    double *data;       /* TEC grid data (tecu) */
    float *rms;         /* RMS values (tecu) */
} tec_t;

typedef struct {        /* satellite fcb data type */
    gtime_t ts,te;      /* start/end time (GPST) */
    double bias[MAXSAT][3]; /* fcb value   (cyc) */
    double std [MAXSAT][3]; /* fcb std-dev (cyc) */
} fcbd_t;

typedef struct {        /* SBAS message type */
    int week,tow;       /* receiption time */
    int prn;            /* SBAS satellite PRN number */
    unsigned char msg[29]; /* SBAS message (226bit) padded by 0 */
} sbsmsg_t;

typedef struct {        /* SBAS messages type */
    int n,nmax;         /* number of SBAS messages/allocated */
    sbsmsg_t *msgs;     /* SBAS messages */
} sbs_t;

typedef struct {        /* SBAS fast correction type */
    gtime_t t0;         /* time of applicability (TOF) */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range-rate correction (RRC) (m/s) */
    double dt;          /* range-rate correction delta-time (s) */
    int iodf;           /* IODF (issue of date fast corr) */
    short udre;         /* UDRE+1 */
    short ai;           /* degradation factor indicator */
} sbsfcorr_t;

typedef struct {        /* SBAS long term satellite error correction type */
    gtime_t t0;         /* correction time */
    int iode;           /* IODE (issue of date ephemeris) */
    double dpos[3];     /* delta position (m) (ecef) */
    double dvel[3];     /* delta velocity (m/s) (ecef) */
    double daf0,daf1;   /* delta clock-offset/drift (s,s/s) */
} sbslcorr_t;

typedef struct {        /* SBAS satellite correction type */
    int sat;            /* satellite number */
    sbsfcorr_t fcorr;   /* fast correction */
    sbslcorr_t lcorr;   /* long term correction */
} sbssatp_t;

typedef struct {        /* SBAS satellite corrections type */
    int iodp;           /* IODP (issue of date mask) */
    int nsat;           /* number of satellites */
    int tlat;           /* system latency (s) */
    sbssatp_t sat[MAXSAT]; /* satellite correction */
} sbssat_t;

typedef struct {        /* SBAS ionospheric correction type */
    gtime_t t0;         /* correction time */
    short lat,lon;      /* latitude/longitude (deg) */
    short give;         /* GIVI+1 */
    float delay;        /* vertical delay estimate (m) */
} sbsigp_t;

typedef struct {        /* IGP band type */
    short x;            /* longitude/latitude (deg) */
    const short *y;     /* latitudes/longitudes (deg) */
    unsigned char bits; /* IGP mask start bit */
    unsigned char bite; /* IGP mask end bit */
} sbsigpband_t;

typedef struct {        /* SBAS ionospheric corrections type */
    int iodi;           /* IODI (issue of date ionos corr) */
    int nigp;           /* number of igps */
    sbsigp_t igp[MAXNIGP]; /* ionospheric correction */
} sbsion_t;

typedef struct {        /* DGPS/GNSS correction type */
    gtime_t t0;         /* correction time */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range rate correction (RRC) (m/s) */
    int iod;            /* issue of data (IOD) */
    double udre;        /* UDRE */
} dgps_t;

typedef struct {        /* SSR correction type */
    gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    double udi[6];      /* SSR update interval (s) */
    int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;           /* issue of data */
    int iodcrc;         /* issue of data crc for beidou/sbas */
    int ura;            /* URA indicator */
    int refd;           /* sat ref datum (0:ITRF,1:regional) */
    double deph [3];    /* delta orbit {radial,along,cross} (m) */
    double ddeph[3];    /* dot delta orbit {radial,along,cross} (m/s) */
    double dclk [3];    /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    double hrclk;       /* high-rate clock correction (m) */
    float  cbias[MAXCODE]; /* code biases (m) */
    double pbias[MAXCODE]; /* phase biases (m) */
    float  stdpb[MAXCODE]; /* std-dev of phase biases (m) */
    double yaw_ang,yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
    unsigned char update; /* update flag (0:no update,1:update) */
} ssr_t;

typedef struct {        /* QZSS LEX message type */
    int prn;            /* satellite PRN number */
    int type;           /* message type */
    int alert;          /* alert flag */
    unsigned char stat; /* signal tracking status */
    unsigned char snr;  /* signal C/N0 (0.25 dBHz) */
    unsigned int ttt;   /* tracking time (ms) */
    unsigned char msg[212]; /* LEX message data part 1695 bits */
} lexmsg_t;

typedef struct {        /* QZSS LEX messages type */
    int n,nmax;         /* number of LEX messages and allocated */
    lexmsg_t *msgs;     /* LEX messages */
} lex_t;

typedef struct {        /* QZSS LEX ephemeris type */
    gtime_t toe;        /* epoch time (GPST) */
    gtime_t tof;        /* message frame time (GPST) */
    int sat;            /* satellite number */
    unsigned char health; /* signal health (L1,L2,L1C,L5,LEX) */
    unsigned char ura;  /* URA index */
    double pos[3];      /* satellite position (m) */
    double vel[3];      /* satellite velocity (m/s) */
    double acc[3];      /* satellite acceleration (m/s2) */
    double jerk[3];     /* satellite jerk (m/s3) */
    double af0,af1;     /* satellite clock bias and drift (s,s/s) */
    double tgd;         /* TGD */
    double isc[8];      /* ISC */
} lexeph_t;

typedef struct {        /* QZSS LEX ionosphere correction type */
    gtime_t t0;         /* epoch time (GPST) */
    double tspan;       /* valid time span (s) */
    double pos0[2];     /* reference position {lat,lon} (rad) */
    double coef[3][2];  /* coefficients lat x lon (3 x 2) */
} lexion_t;

typedef struct {        /* stec data type */
    gtime_t time;       /* time (GPST) */
    unsigned char sat;  /* satellite number */
    double ion;         /* slant ionos delay (m) */
    float std;          /* std-dev (m) */
    float azel[2];      /* azimuth/elevation (rad) */
    unsigned char flag; /* fix flag */
} stec_t;

typedef struct {        /* trop data type */
    gtime_t time;       /* time (GPST) */
    double trp[3];      /* zenith tropos delay/gradient (m) */
    float std[3];       /* std-dev (m) */
} trop_t;

typedef struct {        /* ppp corrections type */
    int nsta;           /* number of stations */
    char stas[MAXSTA][8]; /* station names */
    double rr[MAXSTA][3]; /* station ecef positions (m) */
    int ns[MAXSTA],nsmax[MAXSTA]; /* number of stec data */
    int nt[MAXSTA],ntmax[MAXSTA]; /* number of trop data */
    stec_t *stec[MAXSTA]; /* stec data */
    trop_t *trop[MAXSTA]; /* trop data */
} pppcorr_t;

typedef struct {        /* navigation data type */
    int n,nmax;         /* number of broadcast ephemeris */
    int ng,ngmax;       /* number of glonass ephemeris */
    int ns,nsmax;       /* number of sbas ephemeris */
    int ne,nemax;       /* number of precise ephemeris */
    int nc,ncmax;       /* number of precise clock */
    int na,namax;       /* number of almanac data */
    int nt,ntmax;       /* number of tec grid data */
    int nf,nfmax;       /* number of satellite fcb data */
    eph_t *eph;         /* GPS/QZS/GAL ephemeris */
    geph_t *geph;       /* GLONASS ephemeris */
    seph_t *seph;       /* SBAS ephemeris */
    peph_t *peph;       /* precise ephemeris */
    pclk_t *pclk;       /* precise clock */
    alm_t *alm;         /* almanac data */
    tec_t *tec;         /* tec grid data */
    fcbd_t *fcb;        /* satellite fcb data */
    erp_t  erp;         /* earth rotation parameters */
    double utc_gps[4];  /* GPS delta-UTC parameters {A0,A1,T,W} */
    double utc_glo[4];  /* GLONASS UTC GPS time parameters */
    double utc_gal[4];  /* Galileo UTC GPS time parameters */
    double utc_qzs[4];  /* QZS UTC GPS time parameters */
    double utc_cmp[4];  /* BeiDou UTC parameters */
    double utc_irn[4];  /* IRNSS UTC parameters */
    double utc_sbs[4];  /* SBAS UTC parameters */
    double ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int leaps;          /* leap seconds (s) */
    double lam[MAXSAT][NFREQ+NEXOBS*2];  /* carrier wave lengths (m) */
    double cbias[MAXSAT][3];    /* satellite dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
    double rbias[MAXRCV][2][3]; /* receiver dcb (0:p1-p2,1:p1-c1,2:p2-c2) (m) */
    double wlbias[MAXSAT];      /* wide-lane bias (cycle) */
    double glo_cpbias[4];       /* glonass code-phase bias {1C,1P,2C,2P} (m) */
    char glo_fcn[MAXPRNGLO+1];  /* glonass frequency channel number + 8 */
    pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
    sbssat_t sbssat;    /* SBAS satellite corrections */
    sbsion_t sbsion[MAXBAND+1]; /* SBAS ionosphere corrections */
    dgps_t dgps[MAXSAT]; /* DGPS corrections */
    ssr_t ssr[MAXSAT];   /* SSR corrections */
    lexeph_t lexeph[MAXSAT]; /* LEX ephemeris */
    lexion_t lexion;     /* LEX ionosphere correction */
    pppcorr_t pppcorr;   /* ppp corrections */
    sigind_t sind[2][7]; /* observation data signal index
                          * sind[rcv][0]: GPS
                          * sind[rcv][1]: GLONASS
                          * sind[rcv][2]: Galileo
                          * sind[rcv][3]: QZSS
                          * sind[rcv][4]: SBAS
                          * sind[rcv][5]: BDS
                          * sind[rcv][6]: IRNS
                          * rcv=0: rover,rcv=1: base*/
} nav_t;

typedef struct {          /* station parameter type */
    char name   [MAXANT]; /* marker name */
    char marker [MAXANT]; /* marker number */
    char antdes [MAXANT]; /* antenna descriptor */
    char antsno [MAXANT]; /* antenna serial number */
    char rectype[MAXANT]; /* receiver type descriptor */
    char recver [MAXANT]; /* receiver firmware version */
    char recsno [MAXANT]; /* receiver serial number */
    int antsetup;         /* antenna setup id */
    int itrf;             /* ITRF realization year */
    int deltype;          /* antenna delta type (0:enu,1:xyz) */
    double pos[3];        /* station position (ecef) (m) */
    double del[3];        /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;           /* antenna height (m) */
} sta_t;

typedef struct {        /* NovAtel OEM6 velocity solution type */
    gtime_t time;       /* solution time */
    double v[3],qv[6];  /* velocity and it covariance in ecef */
    int type;           /* velocity solution type (1: instantaneous doppler,2: difference from successive position) */
} solvel_t;

typedef struct {        /* solution type */
    gtime_t time;       /* time (GPST) */
    double rr[9];       /* position/velocity/acceleration (m|m/s|m/s^2) */
                        /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    float  qr[6],pqr[6];/* position variance/covariance (m^2) */
                        /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                        /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    float  qv[6];       /* velocity variance/covariance (m^2/s^2) */
    float  qa[6];       /* vehicle attitude variance/covariance (rad^2)
　　　　　　　　　　　　 * {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} */
    float  qc[6];       /* receiver clock bias to time systems variance */
    float  qcr;         /* receiver clock drift about GPST variance */
    double dtr[6];      /* receiver clock bias to time systems (s) */
    double dtrr;        /* receiver clock drift (m/s) about GPST */
    unsigned char type; /* type (0:xyz-ecef,1:enu-baseline) */
    unsigned char stat,pstat; /* solution status (SOLQ_???) */
    unsigned char ns;   /* number of valid satellites */
    unsigned char ista; /* ins update status (INSS_???) */
    float age;          /* age of differential (s) */
    float ratio;        /* AR ratio factor for valiation */
    float thres;        /* AR ratio threshold for valiation */
    float wlratio;      /* AR ratio for WL ambiguity */

    double att[3];      /* vehicle attitude (roll,pitch,yaw/deg) */
    double ab[3],vb[3]; /* vehicle acceleration and velocity expressed in frd-frame */
    double ba[3],bg[3]; /* acceleration bias and gyro bias */
    double Ma[9],Mg[9]; /* Ma and Mg */
    double os,vr[3];    /* odometry scale factor and odometry velocity in rear frame */
    double dv[3];       /* doppler velocity expressed in ecef */
    double oxyz[3];     

    solvel_t sol_vel;   /* NovAtel OEM6 velocity solution */
    imud_t imu;         /* correctioned imu measurement data */
    vostate_t vo;       /* vo states */
} sol_t;

typedef struct {        /* solution buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    int cyclic;         /* cyclic buffer flag */
    int start,end;      /* start/end index */
    gtime_t time;       /* current solution time */
    sol_t *data;        /* solution data */
    double rb[3];       /* reference position {x,y,z} (ecef) (m) */
    double dt;          /* time difference between solutions */
    unsigned char buff[MAXSOLMSG+1]; /* message buffer */
    int nb;             /* number of byte in message buffer */
} solbuf_t;

typedef struct {        /* solution status type */
    gtime_t time;       /* time (GPST) */
    unsigned char sat;  /* satellite number */
    unsigned char frq;  /* frequency (1:L1,2:L2,...) */
    float az,el;        /* azimuth/elevation angle (rad) */
    float resp;         /* pseudorange residual (m) */
    float resc;         /* carrier-phase residual (m) */
    unsigned char flag; /* flags: (vsat<<5)+(slip<<3)+fix */
    unsigned char snr;  /* signal strength (0.25 dBHz) */
    unsigned short lock;  /* lock counter */
    unsigned short outc;  /* outage counter */
    unsigned short slipc; /* slip counter */
    unsigned short rejc;  /* reject counter */
} solstat_t;

typedef struct {        /* solution status buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    solstat_t *data;    /* solution status data */
} solstatbuf_t;

typedef struct {        /* RTCM control struct type */
    int staid;          /* station id */
    int stah;           /* station health */
    int seqno;          /* sequence number for rtcm 2 or iods msm */
    int outtype;        /* output message type */
    gtime_t time;       /* message time */
    gtime_t time_s;     /* message start time */
    obs_t obs;          /* observation data (uncorrected) */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    dgps_t *dgps;       /* output of dgps corrections */
    ssr_t ssr[MAXSAT];  /* output of ssr corrections */
    char msg[128];      /* special message */
    char msgtype[256];  /* last message type */
    char msmtype[6][128]; /* msm signal types */
    int obsflag;        /* obs data complete flag (1:ok,0:not complete) */
    int ephsat;         /* update satellite of ephemeris */
    double cp[MAXSAT][NFREQ+NEXOBS]; /* carrier-phase measurement */
    unsigned short lock[MAXSAT][NFREQ+NEXOBS]; /* lock time */
    unsigned short loss[MAXSAT][NFREQ+NEXOBS]; /* loss of lock count */
    gtime_t lltime[MAXSAT][NFREQ+NEXOBS]; /* last lock time */
    int nbyte;          /* number of bytes in message buffer */ 
    int nbit;           /* number of bits in word buffer */ 
    int len;            /* message length (bytes) */
    unsigned char buff[1200]; /* message buffer */
    unsigned int word;  /* word buffer for rtcm 2 */
    unsigned int nmsg2[100]; /* message count of RTCM 2 (1-99:1-99,0:other) */
    unsigned int nmsg3[400]; /* message count of RTCM 3 (1-299:1001-1299,300-399:2000-2099,0:ohter) */
    char opt[256];      /* RTCM dependent options */
} rtcm_t;

typedef struct {        /* rinex control struct type */
    gtime_t time;       /* message time */
    double ver;         /* rinex version */
    char   type;        /* rinex file type ('O','N',...) */
    int    sys;         /* navigation system */
    int    tsys;        /* time system */
    char   tobs[7][MAXOBSTYPE][4]; /* rinex obs types */
    obs_t  obs;         /* observation data */
    nav_t  nav;         /* navigation data */
    sta_t  sta;         /* station info */
    int    ephsat;      /* ephemeris satellite number */
    char   opt[256];    /* rinex dependent options */
} rnxctr_t;

typedef struct {        /* option type */
    const char *name;   /* option name */
    int format;         /* option format (0:int,1:double,2:string,3:enum) */
    void *var;          /* pointer to option variable */
    const char *comment; /* option comment/enum labels/unit */
} opt_t;

typedef struct {        /* extended receiver error model */
    int ena[4];         /* model enabled */
    double cerr[4][NFREQ*2]; /* code errors (m) */
    double perr[4][NFREQ*2]; /* carrier-phase errors (m) */
    double gpsglob[NFREQ];   /* gps-glonass h/w bias (m) */
    double gloicb [NFREQ];   /* glonass interchannel bias (m/fn) */
} exterr_t;

typedef struct {        /* SNR mask type */
    int ena[2];         /* enable flag {rover,base} */
    double mask[NFREQ][9]; /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;

typedef struct {        /* processing options type */
    int mode;           /* positioning mode (PMODE_???) */
    int soltype;        /* solution type (0:forward,1:backward,2:combined) */
    int nf;             /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
    int navsys;         /* navigation system */
    double elmin;       /* elevation mask angle (rad) */
    snrmask_t snrmask;  /* SNR mask */
    int sateph;         /* satellite ephemeris/clock (EPHOPT_???) */
    int modear;         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
    int glomodear;      /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
    int bdsmodear;      /* BeiDou AR mode (0:off,1:on) */
    int maxout;         /* obs outage count to reset bias */
    int minlock;        /* min lock count to fix ambiguity */
    int minfix;         /* min fix count to hold ambiguity */
    int armaxiter;      /* max iteration to resolve ambiguity */
    int ionoopt;        /* ionosphere option (IONOOPT_???) */
    int tropopt;        /* troposphere option (TROPOPT_???) */
    int dynamics;       /* dynamics model (0:none,1:velociy,2:accel) */
    int tidecorr;       /* earth tide correction (0:off,1:solid,2:solid+otl+pole) */
    int niter;          /* number of filter iteration */
    int codesmooth;     /* code smoothing window size (0:none) */
    int intpref;        /* interpolate reference obs (for post mission) */
    int sbascorr;       /* SBAS correction options */
    int sbassatsel;     /* SBAS satellite selection (0:all) */
    int rovpos;         /* rover position for fixed mode */
    int refpos;         /* base position for relative mode */
                        /* (0:pos in prcopt,  1:average of single pos, */
                        /*  2:read from file, 3:rinex header, 4:rtcm pos) */
    double eratio[NFREQ]; /* code/phase error ratio */
    double err[5];      /* measurement error factor */
                        /* [0]:reserved */
                        /* [1-3]:error factor a/b/c of phase (m) */
                        /* [4]:doppler frequency (hz) */
    double std[3];      /* initial-state std [0]bias,[1]iono [2]trop */
    double prn[6];      /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv [5] pos */
    double sclkstab;    /* satellite clock stability (sec/sec) */
    double thresar[8];  /* AR validation threshold */
    double elmaskar;    /* elevation mask of AR for rising satellite (deg) */
    double elmaskhold;  /* elevation mask to hold ambiguity (deg) */
    double thresslip;   /* slip threshold of geometry-free phase (m) */
    double maxtdiff;    /* max difference of time (sec) */
    double maxinno;     /* reject threshold of innovation (m) */
    double maxgdop;     /* reject threshold of gdop */
    double baseline[2]; /* baseline length constraint {const,sigma} (m) */
    double ru[3];       /* rover position for fixed mode {x,y,z} (ecef) (m) */
    double rb[3];       /* base position for relative mode {x,y,z} (ecef) (m) */
    char anttype[2][MAXANT]; /* antenna types {rover,base} */
    double antdel[2][3]; /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
    double eps;          /* time internal between epochs */
    pcv_t pcvr[2];       /* receiver antenna parameters {rov,base} */
    unsigned char exsats[MAXSAT]; /* excluded satellites (1:excluded,2:included) */
    int  maxaveep;       /* max averaging epoches */
    int  initrst;        /* initialize by restart */
    int  outsingle;      /* output single by dgps/float/fix/ppp outage */
    char rnxopt[2][256]; /* rinex options {rover,base} */
    int  posopt[6];      /* positioning options */
    int  syncsol;        /* solution sync mode (0:off,1:on) */
    double odisp[2][6*11]; /* ocean tide loading parameters {rov,base} */
    exterr_t exterr;     /* extended receiver error model */
    int freqopt;         /* disable L2-AR */
    char pppopt[256];    /* ppp option */
    char monodir[MAXSTRPATH]; /* mono camera image file directory */
    int aropt;           /* AR opt for WLNL-fix ambiguity (0: ROUND fix NL-amb, 1: ILS fix NL-amb) */
    int adjobs;          /* adjust observation data */
    int gtfmt;           /* ground truth solution format */

    int inherit;          /* ambiguity inherit */
    int inherit_age;      /* inherit double-difference ambiguity age (times of epoch) */
    int inherit_fixc;     /* counts of double-difference ambiguity inherit fix */
    double inherit_thres; /* threshold of inherixing double-difference ambiguity */
    double inherit_IF;    /* inherit double-difference ambiguity impact factor */

    insopt_t insopt;      /* ins option */
    gtime_t ext[16][2];   /* exclude gnss measurement data (included gsof+observation data) for processing */
    sigind_t sind[2][7];  /* observation signal information,0: rover,1: base */
} prcopt_t;

typedef struct {        /* solution options type */
    int posf;           /* solution format (SOLF_???) */
    int times;          /* time system (TIMES_???) */
    int timef;          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;          /* time digits under decimal point */
    int degf;           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outhead;        /* output header (0:no,1:yes) */
    int outopt;         /* output processing options (0:no,1:yes) */
    int outvel;         /* output velocity options (0:no,1:yes) */
    int outatt;         /* output attitude options (0:no,1:yes) */
    int outacc;         /* output acceleration in body frame options (0:no,1:yes) */
    int outvb;          /* output velocity in body frame options (0:no,1:yes) */
    int outba,outbg;    /* output gyro and acceleration bias options (0:no,1:yes) */
    int datum;          /* datum (0:WGS84,1:Tokyo) */
    int height;         /* height (0:ellipsoidal,1:geodetic) */
    int geoid;          /* geoid model (0:EGM96,1:JGD2000) */
    int solstatic;      /* solution of static mode (0:all,1:single) */
    int sstat;          /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;          /* debug trace level (0:off,1-5:debug) */
    double nmeaintv[2]; /* nmea output interval (s) (<0:no,0:all) */
                        /* nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
    char sep[64];       /* field separator */
    char prog[64];      /* program name */
    double maxsolstd;   /* max std-dev for solution output (m) (0:all) */
    int ins_posf;       /* ins position output format: xyz/llh */
    int odo;            /* odometry output options */
    int outclk;         /* output receiver clock bias */
    int dopp;           /* doppler output options */
    int wlratio;        /* WL ambiguity fix ratio */
    int outimuraw;      /* output imu raw data option (0:no,1:yes) */
    int outMa,outMg;    /* output Ma and Mg options (0:no,1:yes) */
    int outenu;         /* output ENU of ins position */
    int outvo;          /* output vo states */
} solopt_t;

typedef struct {              /* file options type */
    char satantp[MAXSTRPATH]; /* satellite antenna parameters file */
    char rcvantp[MAXSTRPATH]; /* receiver antenna parameters file */
    char navfile[MAXSTRPATH]; /* gps navigation data file */
    char bdsfile[MAXSTRPATH]; /* bds navigation data file */
    char glofile[MAXSTRPATH]; /* glo navigation data file */
    char mixfile[MAXSTRPATH]; /* mix navigation data file */
    char stapos [MAXSTRPATH]; /* station positions file */
    char geoid  [MAXSTRPATH]; /* external geoid data file */
    char iono   [MAXSTRPATH]; /* ionosphere data file */
    char dcb    [MAXSTRPATH]; /* dcb data file */
    char eop    [MAXSTRPATH]; /* eop data file */
    char blq    [MAXSTRPATH]; /* ocean tide loading blq file */
    char tempdir[MAXSTRPATH]; /* ftp/http temporaly directory */
    char geexe  [MAXSTRPATH]; /* google earth exec file */
    char solstat[MAXSTRPATH]; /* solution statistics file */
    char trace  [MAXSTRPATH]; /* debug trace file */
    char gtfile [MAXSTRPATH]; /* ground truth file */
    char magfile[MAXSTRPATH]; /* geomagnetic field model coefficients */
    
} filopt_t;

typedef struct {        /* RINEX options type */
    gtime_t ts,te;      /* time start/end */
    double tint;        /* time interval (s) */
    double ttol;        /* time tolerance (s) */
    double tunit;       /* time unit for multiple-session (s) */
    double rnxver;      /* RINEX version */
    int navsys;         /* navigation system */
    int obstype;        /* observation type */
    int freqtype;       /* frequency type */
    char mask[7][64];   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    char staid [32];    /* station id for rinex file name */
    char prog  [32];    /* program */
    char runby [32];    /* run-by */
    char marker[64];    /* marker name */
    char markerno[32];  /* marker number */
    char markertype[32]; /* marker type (ver.3) */
    char name[2][32];   /* observer/agency */
    char rec [3][32];   /* receiver #/type/vers */
    char ant [3][32];   /* antenna #/type */
    double apppos[3];   /* approx position x/y/z */
    double antdel[3];   /* antenna delta h/e/n */
    char comment[MAXCOMMENT][64]; /* comments */
    char rcvopt[256];   /* receiver dependent options */
    unsigned char exsats[MAXSAT]; /* excluded satellites */
    int scanobs;        /* scan obs types */
    int outiono;        /* output iono correction */
    int outtime;        /* output time system correction */
    int outleaps;       /* output leap seconds */
    int autopos;        /* auto approx position */
    int halfcyc;        /* half cycle correction */
    int sep_nav;        /* separated nav files */
    gtime_t tstart;     /* first obs time */
    gtime_t tend;       /* last obs time */
    gtime_t trtcm;      /* approx log start time for rtcm */
    char tobs[7][MAXOBSTYPE][4]; /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    int nobs[7];        /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
} rnxopt_t;

typedef struct {        /* satellite status type */
    unsigned char sys;  /* navigation system */
    unsigned char vs;   /* valid satellite flag single */
    double azel[2];     /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ]; /* residuals of pseudorange (m) */
    double resc[NFREQ]; /* residuals of carrier-phase (m) */
    unsigned char vsat[NFREQ]; /* valid satellite flag */
    unsigned char vsatc[NFREQ];/* valid satellite flag for code */
    unsigned char snr [NFREQ]; /* signal strength (0.25 dBHz) */
    unsigned char fix [NFREQ]; /* ambiguity fix flag (1:fix,2:float,3:hold) */
    unsigned char slip[NFREQ]; /* cycle-slip flag */
    unsigned char half[NFREQ]; /* half-cycle valid flag */
    int lock [NFREQ];   /* lock counter of phase */
    unsigned int index[NFREQ]; /* index of satellite in phase double difference residuals vector */
    unsigned int cind [NFREQ]; /* index of satellite in code double difference residuals vector */
    unsigned int outc [NFREQ]; /* obs outage counter of phase */
    unsigned int slipc[NFREQ]; /* cycle-slip counter */
    unsigned int rejc [NFREQ]; /* reject counter */
    unsigned int news [2*NFREQ];
                        /* new phase double difference satellite flag */
    double  gf;         /* geometry-free phase L1-L2 (m) */
    double  gf2;        /* geometry-free phase L1-L5 (m) */
    double  mw;         /* MW-LC (m) */
    double  phw;        /* phase windup (cycle) */
    gtime_t pt[2][NFREQ]; /* previous carrier-phase time */
    double  ph[2][NFREQ]; /* previous carrier-phase observable (cycle) */
    double  sdi[NFREQ];   /* single-differenced pseudorange observable by INS */
    double  sdg[NFREQ];   /* single-differenced pseudorange observable by GNSS */
} ssat_t;

typedef struct {        /* ambiguity control type */
    gtime_t epoch[4];   /* last epoch */
    int n[4];           /* number of epochs */
    double LC [4];      /* linear combination average */
    double LCv[4];      /* linear combination variance */
    int fixcnt;         /* fix count */
    char flags[MAXSAT]; /* fix flags */
} ambc_t;

typedef struct {        /* double-difference ambiguity type */
    gtime_t time;       /* observation time */
    int sat1,sat2;      /* double-difference ambiguity satellite */
    int f;              /* frequency no. */
    int c;              /* count of fix */
    double ratio;       /* LAMBDA ratio */
    double bias;        /* ambiguity value */
} ddamb_t;

typedef struct {        /* double-difference satellite */
    int sat1,sat2;      /* double difference satellite no. */
    int f;              /* frequency no. */
    int flag;
} ddsat_t;

typedef struct {
    int nb;             /* numbers of double-difference ambiguity */
    int nmax;           /* max numbers of double-difference ambiguity */
    ddamb_t *amb;       /* double difference ambiguity list */
} amb_t;

typedef struct {        /* RTK control/result type */
    sol_t  sol;         /* RTK solution */
    double rb[6];       /* base position/velocity (ecef) (m|m/s) */
    double tt;          /* time difference between current and previous (s) */
    double *x, *P;      /* float states and their covariance */
    double *xa,*Pa;     /* fixed states and their covariance */
    int nx,na;          /* number of float states/fixed states */
    int nfix;           /* number of continuous fixes of ambiguity */
    int neb;            /* bytes in error message buffer */
    int ns;             /* number of double difference satellites */
    int refsat[NUMSYS][2*NFREQ]; /* reference satellite of double-difference residuals* (0:gps/qzs/sbs,1:glo,2:gal,3:bds) */
    int failc;                   /* counts of solution valid fail */
    char errbuf[MAXERRMSG];      /* error message buffer */
    prcopt_t opt;                /* processing options */
    ambc_t ambc[MAXSAT];         /* ambiguity control */
    ssat_t ssat[MAXSAT];         /* satellite status */
    insstate_t ins;              /* ins states */
    amb_t bias;                  /* double-difference ambiguity list */
    amb_t wlbias;                /* WL double-difference ambiguity list */
    ddsat_t sat[MAXSAT];         /* double difference satellite list */
} rtk_t;

typedef struct half_cyc_tag {  /* half-cycle correction list type */
    unsigned char sat;         /* satellite number */
    unsigned char freq;        /* frequency number (0:L1,1:L2,2:L5) */
    unsigned char valid;       /* half-cycle valid flag */
    char corr;                 /* half-cycle corrected (x 0.5 cyc) */
    gtime_t ts,te;             /* time start, time end */
    struct half_cyc_tag *next; /* pointer to next correction */
} half_cyc_t;

typedef struct {        /* GSOF: the base station position and its quality */
    gtime_t t;          /* GPS time */
    double pos[3];      /* latitude,longitude and height */
    int quality;        /* The quality of the base station position:
                           0: Fix not available or invalid
                           1: Autonomous GPS fix
                           2: Differential SBAS or OmniSTAR VBS
                           4: RTK Fixed,xFill
                           5: OmniSTAR XP,OmniSTAR HP,CenterPoint RTX,Float RTK,or Location RTK */
} gsof_base_t;

typedef struct {         /* GSOF message: the clock information. it contains the following data:clock offset,frequency offset */
    int flags;           /* Provides information related to the clock fix process. Defined values are:
                            bit 0 set: clock offset is valid
                            bit 1 set: frequency offset is valid
                            bit 2 set: receiver is in anywhere fix mode */
    double off;          /* the current clock offset in milliseconds */
    double foff;         /* the offset of the local oscillator from the nominal GPS L1 frequency in parts per million */
} gsof_clk_t;

typedef struct {         /* GSOF message data type */
    gtime_t t;           /* position time */
    int ns;              /* number of satellites used to determine the position */
    int status;          /* receiver status code */
    int no;              /* unique number assigned to a group of record packet pages.
                          * prevents page mismatches when multiple sets of record packets exist in output stream */
    int solq;            /* position status */
    int velf;            /* flag of velocity */
    double pos[3];       /* ecef position (m) */
    double llh[3];       /* latitude,longitude and height of rover station (rad,m) */
    double vel[3];       /* enu velocity (m/s) */
    double delta[3];     /* earth-centered, earth-fixed x,y,z deltas between the rover and base position (m)*/
    float cov[8];        /* earth-centered, earth-fixed (ecef) position covariance
                          * [0]: position rms,
                          * [1]: cov-xx(m),
                          * [2]: cov-xy(m),
                          * [3]: cov-xz(m),
                          * [4]: cov-yy(m),
                          * [5]: cov-yz(m),
                          * [6]: cov-zz(m),
                          * [7]: unit variance(m^2)
                          * */
    float sig[6];        /* the position sigma information
                          * [0]: position rms (m)
                          * [1]: east (m)
                          * [2]: north (m)
                          * [3]: east-north (m)
                          * [4]: up
                          * [5]: unit variance,valid only for over-determined solutions.
                          *      unit variance should approach 1.0 value.
                          *      A value of less than 1.0 indicates that apriori variances are too pessimistic
                          * */
    double dop[4];       /* PDOP,HDOP,VDOP,TDOP */
    gsof_base_t base;    /* base station information */
    gsof_clk_t  clk;     /* clock information of rover station */
} gsof_t;

typedef struct {         /* gsof observation data type */
    int n,nmax;          /* number of obervation data/allocated */
    gsof_t *data;        /* gsof observation data */
} gsof_data_t;

typedef struct {         /* temp variable for real time decode rinex obs data */
    gtime_t time;        /* observation data time */
    int start;           /* flag of start decode obs data */
    int i,n,nsat;        /* index/number of satellites */
    int signal;          /* signal index set */
    int endhead;         /* flag of rinex obs header */
    char tobs[NUMSYS][MAXOBSTYPE][4];
                         /* type of observation data */
    obsd_t obs[MAXOBS];  /* observation data */
    sigind_t index[7];   /* signal index */
} rinex_t;

typedef struct {         /* vehicle attitude solution data */
    gtime_t time;        /* attitude solution data time */
    double val[3],std[3];/* attitude solution data and its variance (roll,pitch,yaw/rad) */
} att_t;

typedef struct {        /* receiver raw data control type */
    gtime_t time;       /* message time */
    gtime_t tobs[MAXSAT][NFREQ+NEXOBS]; /* observation data time */
    obs_t obs;          /* observation data */
    sol_t sol;          /* solution data from ublox message */
    obs_t obuf;         /* observation data buffer */
    img_t img;          /* image raw data */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    gsof_t gsof;        /* position information from a trimble proprietary format
                         * and can be used to send information such as position and status to a third-party device */
    gsof_data_t gsofb;  /* gsof data for backward solution */
    att_t att;          /* vehicle attitude solution */
    imud_t imu;         /* imu measurement data (use in decode M39 data) */
    imu_t imut;         /* imu measurement data for an epoch of GNSS measurement */
    imu_t imub;         /* imu measurement data for backward solution */
    odo_t odo;          /* odometry measurement data */
    m39_mix_t m39;      /* m39 mix raw data (inclued image/imu) */
    pose_meas_t pose;   /* pose measurement data from camera or dual ant. */
    int ephsat;         /* sat number of update ephemeris (0:no satellite) */
    sbsmsg_t sbsmsg;    /* SBAS message */
    char msgtype[256];  /* last message type */
    char monodir[PATH_MAX]; /* mono camera image data directory */
    unsigned char subfrm[MAXSAT][380];  /* subframe buffer */
    lexmsg_t lexmsg;    /* LEX message */
    double lockt[MAXSAT][NFREQ+NEXOBS]; /* lock time (s) */
    double icpp[MAXSAT],off[MAXSAT],icpc; /* carrier params for ss2 */
    double prCA[MAXSAT],dpCA[MAXSAT]; /* L1/CA pseudrange/doppler for javad */
    unsigned char halfc[MAXSAT][NFREQ+NEXOBS]; /* half-cycle add flag */
    unsigned char tflg; /* time valid flag (only use decode ublox message) */
    unsigned char dire; /* direction of input raw data (0:forward,1:backward,2:combined) */
    char freqn[MAXOBS]; /* frequency number for javad */
    int nbyte;          /* number of bytes in message buffer */
    int len;            /* message length (bytes) */
    int iod;            /* issue of data */
    int tod;            /* time of day (ms) */
    int tbase;          /* time base (0:gpst,1:utc(usno),2:glonass,3:utc(su) */
    int flag;           /* general purpose flag */
    int outtype;        /* output message type */
    int curb;           /* index of raw data for backward solution */
    unsigned char buff[MAXRAWLEN]; /* message buffer */
    char opt[256];      /* receiver dependent options */
    half_cyc_t *half_cyc; /* half-cycle correction list */
    rinex_t rinex;        /* temp variable for real time decode rinex obs data */
    int format;         /* receiver stream format */
    void *rcv_data;     /* receiver dependent data */
    void *optp;         /* process options */
    void *strp;         /* stream pointer */
    int imufmt;         /* imu data type */
} raw_t;

/* type definitions ----------------------------------------------------------*/
typedef struct vt_tag { /* virtual console type */
    int state;          /* state(0:close,1:open) */
    int type;           /* type (0:dev,1:telnet) */
    int in,out;         /* input/output file descriptor */
    int n,nesc;         /* number of line buffer/escape */
    int cur;            /* cursor position */
    int cur_h;          /* current history */
    int brk;            /* break status */
    int blind;          /* blind inpu mode */
    int debug;          /* debug mode */
    struct termios tio; /* original terminal attribute */
    char buff[MAXBUFF]; /* line buffer */
    char esc[8];        /* escape buffer */
    char *hist[MAXHIST];/* history buffer */
    FILE *logfp;        /* log file pointer */
} vt_t;

typedef struct {        /* stream type */
    int type;           /* type (STR_???) */
    int mode;           /* mode (STR_MODE_?) */
    int state;          /* state (-1:error,0:close,1:open) */
    unsigned int inb,inr;   /* input bytes/rate */
    unsigned int outb,outr; /* output bytes/rate */
    unsigned int tick_i; /* input tick tick */
    unsigned int tick_o; /* output tick */
    unsigned int tact;   /* active tick */
    unsigned int inbt,outbt; /* input/output bytes at tick */
    lock_t lock;        /* lock flag */
    void *port;         /* type dependent port control struct */
    char path[MAXSTRPATH]; /* stream path */
    char msg [MAXSTRMSG];  /* stream message */
} stream_t;

typedef struct {        /* stream converter type */
    int itype,otype;    /* input and output stream type */
    int nmsg;           /* number of output messages */
    int msgs[64];       /* output message types */
    double tint[64];    /* output message intervals (s) */
    unsigned int tick[64]; /* cycle tick of output message */
    int ephsat[64];     /* satellites of output ephemeris */
    int stasel;         /* station info selection (0:remote,1:local) */
    rtcm_t rtcm;        /* rtcm input data buffer */
    raw_t raw;          /* raw  input data buffer */
    rtcm_t out;         /* rtcm output data buffer */
} strconv_t;

typedef struct {         /* stream server type */
    int state;           /* server state (0:stop,1:running) */
    int cycle;           /* server cycle (ms) */
    int buffsize;        /* input/monitor buffer size (bytes) */
    int nmeacycle;       /* NMEA request cycle (ms) (0:no) */
    int relayback;       /* relay back of output streams (0:no) */
    int nstr;            /* number of streams (1 input + (nstr-1) outputs */
    int npb;             /* data length in peek buffer (bytes) */
    char cmds_periodic[16][MAXRCVCMD]; /* periodic commands */
    double nmeapos[3];   /* NMEA request position (ecef) (m) */
    unsigned char *buff; /* input buffers */
    unsigned char *pbuf; /* peek buffer */
    unsigned int tick;   /* start tick */
    stream_t stream[16]; /* input/output streams */
    strconv_t *conv[16]; /* stream converter */
    thread_t thread;     /* server thread */
    lock_t lock;         /* lock flag */
} strsvr_t;

typedef struct {         /* time synchronization index in buffer */
    int ni,nip,imu;      /* current/precious number and sync index of imu measurement data  */
    int nr,nrp,rover;    /* current/precious number and sync index of observation data */
    int ns,nsp,pvt;      /* current/precious number and sync index of pvt solution data */
    int nb,nbp,base;     /* current/precious number and sync index of base observation data */
    int nm,nmp,img;      /* current/precious number and sync index of image raw data */
    int np,npp,pose;     /* current/precious number and sync index of pose measurement */
    int of[6];           /* overflow flag (rover,base,imu,pvt,image,pose) */
    unsigned int tali[6];/* time alignment for data (rover-base,pvt-imu,rover-base-imu,imu-image-pvt,pose-imu) */
    unsigned int ws;     /* search window size */
    double dt[6];        /* time difference between input streams */
    gtime_t time[6];     /* current time of rover,base,imu,pvt,image and pose measurement */
} syn_t;

typedef struct {        /* RTK server type */
    int pause;          /* pause program (0:off,1:on ) */
    int reinit;         /* re-initial ins states (0:off,1:on) */
    int state;          /* server state (0:stop,1:running) */
    int cycle;          /* processing cycle (ms) */
    int nmeacycle;      /* NMEA request cycle (ms) (0:no req) */
    int nmeareq;        /* NMEA request (0:no,1:nmeapos,2:single sol) */
    double nmeapos[3];  /* NMEA request position (ecef) (m) */
    int buffsize;       /* input buffer size (bytes) */
    int format[7];      /* input format {rov,base,corr,sol,imu,image} */
    int navsel;         /* ephemeris select (0:all,1:rover,2:base,3:corr) */
    int nsbs;           /* number of sbas message */
    int nsol;           /* number of solution buffer */
    int nb [7];         /* bytes in input buffers {rov,base,corr,sol,imu,image} */
    int nsb[2];         /* bytes in solution buffers */
    int npb[7];         /* bytes in input peek buffers {rov,base,corr,sol,imu,image} */
    int cputime;        /* CPU time (ms) for a processing cycle */
    int prcout;         /* missing observation data count */
    int iprcout;        /* imu missing observation data count */
    int gprcout;        /* solution missing data count */
    int nave;           /* number of averaging base pos */
    int week;           /* GPS week */
    double rb_ave[3];   /* averaging base pos */
    char cmds_periodic[7][MAXRCVCMD]; /* periodic commands */
    char cmd_reset[MAXRCVCMD];/* reset command */
    double bl_reset;          /* baseline length to reset (km) */
    unsigned char *buff[7];   /* input buffers {rov,base,corr,sol,imu,image,pose} */
    unsigned char *sbuf[2];   /* output buffers {sol1,sol2} */
    unsigned char *pbuf[7];   /* peek buffers {rov,base,corr,sol,imu,image,pose} */
    unsigned int nmsg[7][15]; /* input message counts */
    unsigned int tick;        /* start tick */
    char files[3][MAXSTRPATH];/* download paths {rov,base,corr} */
    rtk_t rtk;                /* RTK control/result struct */
    raw_t  raw [7];           /* receiver raw control {rov,base,corr,sol,imu,image,pose} */
    rtcm_t rtcm[7];           /* RTCM control {rov,base,corr,reserve,reserve} */
    gtime_t ftime[3];         /* download time {rov,base,corr} */
    sol_t  solbuf[MAXSOLBUF]; /* solution buffer */
    obs_t  obs[3][MAXOBSBUF]; /* observation data {rov,base,corr} */
    sol_t  pvt[MAXSOLBUF];    /* PVT solutions for ins/gnss couple */
    imud_t imu[MAXIMUBUF];    /* imu measurement data of rover station */
    pose_meas_t pose[MAXPOSEBUF]; /* pose measurement from camera or dual ant. */
    img_t  img[MAXIMGBUF];        /* image raw data from camera measurement */
    syn_t syn;                    /* time synchronization index of buffer */
    nav_t nav;                    /* navigation data */
    sbsmsg_t sbsmsg[MAXSBSMSG];   /* SBAS message buffer */
    stream_t stream[14];/* streams {rov,base,corr,sol,imu,image,pose,sol1,sol2,logr,logb,logc,logi,logs} */
    stream_t *moni;     /* monitor stream */
    stream_t *groundtruth; /* ground-truth path stream */
    vt_t *vt;           /* virtual console */
    solopt_t solopt[2]; /* output solution options {sol1,sol2} */
    solbuf_t gtsols;    /* ground-truth solution data */
    thread_t thread;    /* server thread */
    lock_t lock;        /* lock flag */
} rtksvr_t;

typedef struct {        /* gis data point type */
    double pos[3];      /* point data {lat,lon,height} (rad,m) */
} gis_pnt_t;

typedef struct {        /* gis data polyline type */
    int npnt;           /* number of points */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
    double *pos;        /* position data (3 x npnt) */
} gis_poly_t;

typedef struct {        /* gis data polygon type */
    int npnt;           /* number of points */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
    double *pos;        /* position data (3 x npnt) */
} gis_polygon_t;

typedef struct gisd_tag { /* gis data list type */
    int type;           /* data type (1:point,2:polyline,3:polygon) */
    void *data;         /* data body */
    struct gisd_tag *next; /* pointer to next */
} gisd_t;

typedef struct {        /* gis type */
    char name[MAXGISLAYER][256]; /* name */
    int flag[MAXGISLAYER];       /* flag */
    gisd_t *data[MAXGISLAYER];   /* gis data list */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
} gis_t;

typedef void fatalfunc_t(const char *); /* fatal callback function type */
typedef int8_t s8;                      /* Signed 8-bit integer */
typedef int16_t s16;                    /* Signed 16-bit integer */
typedef int32_t s32;                    /* Signed 32-bit integer */
typedef int64_t s64;                    /* Signed 64-bit integer */
typedef uint8_t u8;                     /* Unsigned 8-bit integer */
typedef uint16_t u16;                   /* Unsigned 16-bit integer */
typedef uint32_t u32;                   /* Unsigned 32-bit integer */
typedef uint64_t u64;                   /* Unsigned 64-bit integer */

/* global variables ----------------------------------------------------------*/
extern const double chisqr[];             /* chi-sqr(n) table (alpha=0.001) */
extern const double lam_carr[];           /* carrier wave length (m) {L1,L2,...} */
extern const prcopt_t prcopt_default;     /* default positioning options */
extern const solopt_t solopt_default;     /* default solution output options */
extern const solopt_t solopt_ins_default; /* default ins solution output options */
extern const solopt_t solopt_vo_default;  /* default vo solution output options */
extern const filopt_t fileopt_default;    /* default file options */
extern const sbsigpband_t igpband1[9][8]; /* SBAS IGP band 0-8 */
extern const sbsigpband_t igpband2[2][5]; /* SBAS IGP band 9-10 */
extern const char *formatstrs[];          /* stream format strings */
extern opt_t sysopts[];                   /* system options table */
extern opt_t insopts[];                   /* ins options table */
extern const double Omge[9];              /* earth rotation matrix in i/e-frame */
extern const quat_t identity_quat;        /* identity quaternion */
extern const double Cen[9];               /* transform matrix of enu-frame convert to ned-frame */
extern const double Crf[9];               /* transform matrix of rfu-frame convert to frd-frame */
extern const char *solqstrs[];            /* solution status strings */

/* satellites, systems, codes functions --------------------------------------*/
EXPORT int  satno   (int sys, int prn);
EXPORT int  satsys  (int sat, int *prn);
EXPORT int  satid2no(const char *id);
EXPORT void satno2id(int sat, char *id);
EXPORT unsigned char obs2code(const char *obs, int *freq);
EXPORT char *code2obs(unsigned char code, int *freq);
EXPORT int  satexclude(int sat, int svh, const prcopt_t *opt);
EXPORT int  testsnr(int base, int freq, double el, double snr,
                    const snrmask_t *mask);
EXPORT void setcodepri(int sys, int freq, const char *pri);
EXPORT int  getcodepri(int sys, unsigned char code, const char *opt);

/* matrix and vector functions -----------------------------------------------*/
EXPORT double *mat  (int n, int m);
EXPORT float  *fmat (int n, int m);
EXPORT int    *imat (int n, int m);
EXPORT double *zeros(int n, int m);
EXPORT float *fzeros(int n, int m);
EXPORT double *eye  (int n);
EXPORT float *feye(int n);
EXPORT double dot (const double *a, const double *b, int n);
EXPORT float dotf(const float *a, const float *b, int n);
EXPORT double norm(const double *a, int n);
EXPORT float normf(const float *a, int n);
EXPORT void cross3(const double *a, const double *b, double *c);
EXPORT int  normv3(const double *a, double *b);
EXPORT void matcpy(double *A, const double *B, int n, int m);
EXPORT void fmatcpy(float *A, const float *B, int n, int m);
EXPORT void matcpy_d2f(float *A, const double *B, int n, int m);
EXPORT void matcpy_f2d(double *A, const float *B, int n, int m);
EXPORT void imatcpy(int *A, const int *B, int n, int m);
EXPORT void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
EXPORT void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D);
EXPORT int  matinv(double *A, int n);
EXPORT int  solve (const char *tr, const double *A, const double *Y, int n,
                   int m, double *X);
EXPORT int  lsq   (const double *A, const double *y, int n, int m, double *x,
                   double *Q);
EXPORT int  filter(double *x, double *P, const double *H, const double *v,
                   const double *R, int n, int m);
EXPORT int  smoother(const double *xf, const double *Qf, const double *xb,
                     const double *Qb, int n, double *xs, double *Qs);
EXPORT void lsmooth3(double *in, double *out, int N);
EXPORT void lsmooth5(double *in, double *out, int N);
EXPORT void lsmooth7(double *in, double *out, int N);
EXPORT void quadsmooth5(double *in, double *out, int N);
EXPORT void matprint (const double *A, int n, int m, int p, int q);
EXPORT void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);

EXPORT void matt(const double *A,int n,int m,double *At);
EXPORT int svd(const double *A,int m,int n,double *U,double *W,double *V);
EXPORT int null(const double *A,int m,int n,double *N,int *p,int *q);
EXPORT void matpow(const double *A,int m,int p,double *B);
EXPORT void dialog(const double *v,int n,double *D);
EXPORT void mat2dmat(const double *a,double **b,int m,int n);
EXPORT void dmat2mat(double **a,double *b,int m,int n);
EXPORT int qr(const double *A,int m,int n,double *Q,double *R,int flag);

EXPORT void asi_blk_mat(double *A,int m,int n,const double *B,int p ,int q,
                        int isr,int isc);
EXPORT int mateigenvalue(const double* A,int n,double *u,double *v);
EXPORT int matdet(const double*A,int n,double*det);
EXPORT double det(const double *A,int n);
EXPORT double stds(const double *val,int n);
EXPORT double avg(const double *val,int n);
EXPORT double norm_distri(const double u);
EXPORT double re_norm(double p);
EXPORT void add_fatal(fatalfunc_t *func);
EXPORT double** dmat(int m,int n);
EXPORT void matrix_udu(u32 n, double *M, double *U, double *D);
/* time and string functions -------------------------------------------------*/
EXPORT double  str2num(const char *s, int i, int n);
EXPORT int     str2time(const char *s, int i, int n, gtime_t *t);
EXPORT void    time2str(gtime_t t, char *str, int n);
EXPORT gtime_t epoch2time(const double *ep);
EXPORT void    time2epoch(gtime_t t, double *ep);
EXPORT gtime_t gpst2time(int week, double sec);
EXPORT double  time2gpst(gtime_t t, int *week);
EXPORT gtime_t gst2time(int week, double sec);
EXPORT double  time2gst(gtime_t t, int *week);
EXPORT gtime_t bdt2time(int week, double sec);
EXPORT double  time2bdt(gtime_t t, int *week);
EXPORT char    *time_str(gtime_t t, int n);

EXPORT gtime_t timeadd  (gtime_t t, double sec);
EXPORT double  timediff (gtime_t t1, gtime_t t2);
EXPORT gtime_t gpst2utc (gtime_t t);
EXPORT gtime_t utc2gpst (gtime_t t);
EXPORT gtime_t gpst2bdt (gtime_t t);
EXPORT gtime_t bdt2gpst (gtime_t t);
EXPORT gtime_t timeget  (void);
EXPORT void    timeset  (gtime_t t);
EXPORT double  time2doy (gtime_t t);
EXPORT double  utc2gmst (gtime_t t, double ut1_utc);
EXPORT int read_leaps(const char *file);
EXPORT double  time2secs(gtime_t t);

EXPORT int adjgpsweek(int week);
EXPORT int adjsind(const prcopt_t *opt,const obsd_t *obs,int *i,int *j,int *k);
EXPORT unsigned int tickget(void);
EXPORT void sleepms(int ms);

EXPORT int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
                   const char *base);
EXPORT int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts,
                    gtime_t te, const char *rov, const char *base);

/* geiod models --------------------------------------------------------------*/
EXPORT int opengeoid(int model, const char *file);
EXPORT void closegeoid(void);
EXPORT double geoidh(const double *pos);

/* coordinates transformation ------------------------------------------------*/
EXPORT void ecef2pos(const double *r, double *pos);
EXPORT void pos2ecef(const double *pos, double *r);
EXPORT void ecef2enu(const double *pos, const double *r, double *e);
EXPORT void enu2ecef(const double *pos, const double *e, double *r);
EXPORT void covenu  (const double *pos, const double *P, double *Q);
EXPORT void covecef (const double *pos, const double *Q, double *P);
EXPORT void xyz2enu (const double *pos, double *E);
EXPORT void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
EXPORT void deg2dms (double deg, double *dms, int ndec);
EXPORT double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
EXPORT void readpos(const char *file, const char *rcv, double *pos);
EXPORT int  sortobs(obs_t *obs);
EXPORT int  sortgsof(gsof_data_t *gsof);
EXPORT void uniqnav(nav_t *nav);
EXPORT int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
EXPORT int  readnav(const char *file, nav_t *nav);
EXPORT int  savenav(const char *file, const nav_t *nav);
EXPORT void freeobs(obs_t *obs);
EXPORT void freenav(nav_t *nav, int opt);
EXPORT int  readblq(const char *file, const char *sta, double *odisp);
EXPORT int  readerp(const char *file, erp_t *erp);
EXPORT int  geterp (const erp_t *erp, gtime_t time, double *val);
EXPORT void freeimudata(imu_t *imu);
EXPORT void freegsofdata(gsof_data_t *data);
EXPORT void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav);
EXPORT void corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav);
/* debug trace functions -----------------------------------------------------*/
EXPORT void traceopen(const char *file);
EXPORT void traceclose(void);
EXPORT void tracelevel(int level);
EXPORT void trace    (int level, const char *format, ...);
EXPORT void tracet   (int level, const char *format, ...);
EXPORT void tracemat (int level, const double *A, int n, int m, int p, int q);
EXPORT void tracemat_std(int level, const double *A, int n, int m, int p, int q);
EXPORT void traceimat(const int *A,int n,int m,int p,int q);
EXPORT void traceobs (int level, const obsd_t *obs, int n);
EXPORT void tracenav (int level, const nav_t *nav);
EXPORT void tracegnav(int level, const nav_t *nav);
EXPORT void tracehnav(int level, const nav_t *nav);
EXPORT void tracepeph(int level, const nav_t *nav);
EXPORT void tracepclk(int level, const nav_t *nav);
EXPORT void traceb   (int level, const unsigned char *p, int n);
EXPORT void tracemr(double **A,int m,int n,int p,int q);
EXPORT void traceobsbuff(rtksvr_t *svr);
EXPORT void tracesync(rtksvr_t *svr);
/* platform dependent functions ----------------------------------------------*/
EXPORT int execcmd(const char *cmd);
EXPORT int expath (const char *path, char *paths[], int nmax);
EXPORT void createdir(const char *path);

/* positioning models --------------------------------------------------------*/
EXPORT double satwavelen(int sat, int frq, const nav_t *nav);
EXPORT double satazel(const double *pos, const double *e, double *azel);
EXPORT double geodist(const double *rs, const double *rr, double *e);
EXPORT void dops(int ns, const double *azel, double elmin, double *dop);
EXPORT void csmooth(obs_t *obs, int ns);
EXPORT double uravalue(int sys, int sva);
EXPORT int uraindex(double value, int sys);

/* atmosphere models ---------------------------------------------------------*/
EXPORT double ionmodel(gtime_t t, const double *ion, const double *pos,
                       const double *azel);
EXPORT double ionmapf(const double *pos, const double *azel);
EXPORT double ionppp(const double *pos, const double *azel, double re,
                     double hion, double *pppos);
EXPORT double tropmodel(gtime_t time, const double *pos, const double *azel,
                        double humi);
EXPORT double tropmapf(gtime_t time, const double *pos, const double *azel,
                       double *mapfw);
EXPORT int iontec(gtime_t time, const nav_t *nav, const double *pos,
                  const double *azel, int opt, double *delay, double *var);
EXPORT void readtec(const char *file, nav_t *nav, int opt);
EXPORT int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var);
EXPORT int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var);
EXPORT double gdelaycorr(const int sys, const double *rr,const double *rs);
/* antenna models ------------------------------------------------------------*/
EXPORT int  readpcv(const char *file, pcvs_t *pcvs);
EXPORT pcv_t *searchpcv(int sat, const char *type, gtime_t time,
                        const pcvs_t *pcvs);
EXPORT void antmodel(const pcv_t *pcv, const double *del, const double *azel,
                     int opt, double *dant);
EXPORT void antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
EXPORT void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
                       double *rmoon, double *gmst);
EXPORT void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr);

/* datum transformation ------------------------------------------------------*/
EXPORT int loaddatump(const char *file);
EXPORT int tokyo2jgd(double *pos);
EXPORT int jgd2tokyo(double *pos);

/* rinex functions -----------------------------------------------------------*/
EXPORT int readrnx (const char *file, int rcv, const char *opt, obs_t *obs,
                    nav_t *nav, sta_t *sta);
EXPORT int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
                    double tint, const char *opt, obs_t *obs, nav_t *nav,
                    sta_t *sta);
EXPORT int readrnxc(const char *file, nav_t *nav);
EXPORT int readgsoff(const char *file,gsof_data_t *gsof);
EXPORT int readimub(const char *file,imu_t* imu,int decfmt,int imufmt,int coor,
                    int valfmt);
EXPORT int rtk_uncompress(const char *file, char *uncfile);
EXPORT int convrnx(int format, rnxopt_t *opt, const char *file, char **ofile);
EXPORT int  init_rnxctr (rnxctr_t *rnx);
EXPORT void free_rnxctr (rnxctr_t *rnx);
EXPORT int  open_rnxctr (rnxctr_t *rnx, FILE *fp);
EXPORT int  input_rnxctr(rnxctr_t *rnx, FILE *fp);
EXPORT int addobsdata(obs_t *obs, const obsd_t *data);
EXPORT int addimudata(imu_t *imu, const imud_t *data);
/* ephemeris and clock functions ---------------------------------------------*/
EXPORT double eph2clk (gtime_t time, const eph_t  *eph);
EXPORT double geph2clk(gtime_t time, const geph_t *geph);
EXPORT double seph2clk(gtime_t time, const seph_t *seph);
EXPORT void eph2pos (gtime_t time, const eph_t  *eph,  double *rs, double *dts,
                     double *var);
EXPORT void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var);
EXPORT void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var);
EXPORT int  peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                     double *rs, double *dts, double *var);
EXPORT void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant);
EXPORT int  satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                   const nav_t *nav, double *rs, double *dts, double *var,
                   int *svh);
EXPORT void satposs(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                    int sateph, double *rs, double *dts, double *var, int *svh);
EXPORT void readsp3(const char *file, nav_t *nav, int opt);
EXPORT int  readsap(const char *file, gtime_t time, nav_t *nav);
EXPORT int  readdcb(const char *file, nav_t *nav, const sta_t *sta);
EXPORT int  readfcb(const char *file, nav_t *nav);
EXPORT void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts);

EXPORT int tle_read(const char *file, tle_t *tle);
EXPORT int tle_name_read(const char *file, tle_t *tle);
EXPORT int tle_pos(gtime_t time, const char *name, const char *satno,
                   const char *desig, const tle_t *tle, const erp_t *erp,
                   double *rs);

/* receiver raw data functions -----------------------------------------------*/
EXPORT unsigned int getbitu(const unsigned char *buff, int pos, int len);
EXPORT int          getbits(const unsigned char *buff, int pos, int len);
EXPORT void setbitu(unsigned char *buff, int pos, int len, unsigned int data);
EXPORT void setbits(unsigned char *buff, int pos, int len, int data);
EXPORT long bin2dec(const unsigned char *b,int n);
EXPORT unsigned int rtk_crc32  (const unsigned char *buff, int len);
EXPORT unsigned int rtk_crc24q (const unsigned char *buff, int len);
EXPORT unsigned short rtk_crc16(const unsigned char *buff, int len);
EXPORT int decode_word (unsigned int word, unsigned char *data);
EXPORT int decode_frame(const unsigned char *buff, eph_t *eph, alm_t *alm,
                        double *ion, double *utc, int *leaps);
EXPORT int test_glostr(const unsigned char *buff);
EXPORT int decode_glostr(const unsigned char *buff, geph_t *geph);
EXPORT int decode_bds_d1(const unsigned char *buff, eph_t *eph);
EXPORT int decode_bds_d2(const unsigned char *buff, eph_t *eph);
EXPORT int decode_gal_inav(const unsigned char *buff, eph_t *eph);

EXPORT int init_raw   (raw_t *raw, int format);
EXPORT void free_raw  (raw_t *raw);
EXPORT int input_raw  (raw_t *raw, int format, unsigned char data);
EXPORT int input_rawf (raw_t *raw, int format, FILE *fp);

EXPORT int init_rt17  (raw_t *raw);
EXPORT int init_cmr   (raw_t *raw);
EXPORT void free_rt17 (raw_t *raw);
EXPORT void free_cmr  (raw_t *raw);
EXPORT void freeigvfp();
EXPORT int update_cmr (raw_t *raw, rtksvr_t *svr, obs_t *obs);

EXPORT int input_oem6_sol (raw_t *raw, unsigned char data);
EXPORT int input_oem6_pose(raw_t *raw, unsigned char data);
EXPORT int input_oem6_raw (raw_t *raw, unsigned char data);
EXPORT int input_igvsim_imu(raw_t *raw,unsigned char data);
EXPORT int input_igvsim_gnss(raw_t *raw, unsigned char data);
EXPORT int input_igvsim_feat(raw_t *raw, unsigned char data);

EXPORT int input_oem4  (raw_t *raw, unsigned char data);
EXPORT int input_oem3  (raw_t *raw, unsigned char data);
EXPORT int input_ubx   (raw_t *raw, unsigned char data);
EXPORT int input_ss2   (raw_t *raw, unsigned char data);
EXPORT int input_cres  (raw_t *raw, unsigned char data);
EXPORT int input_stq   (raw_t *raw, unsigned char data);
EXPORT int input_gw10  (raw_t *raw, unsigned char data);
EXPORT int input_javad (raw_t *raw, unsigned char data);
EXPORT int input_nvs   (raw_t *raw, unsigned char data);
EXPORT int input_bnx   (raw_t *raw, unsigned char data);
EXPORT int input_rt17  (raw_t *raw, unsigned char data);
EXPORT int input_sbf   (raw_t *raw, unsigned char data);
EXPORT int input_cmr   (raw_t *raw, unsigned char data);
EXPORT int input_tersus(raw_t *raw, unsigned char data);
EXPORT int input_lexr  (raw_t *raw, unsigned char data);
EXPORT int input_gsof  (raw_t *raw, unsigned char data);
EXPORT int input_m39   (raw_t *raw, unsigned char data);
EXPORT int input_ubxm8 (raw_t *raw, unsigned char data);
EXPORT int input_ubxsol(raw_t *raw, unsigned char data);
EXPORT int input_rinex (raw_t *raw, unsigned char data);
EXPORT int input_m39_mix(raw_t *raw, unsigned char data);
EXPORT int input_sbp(raw_t *raw, uint8_t data);
EXPORT int input_cnav(raw_t *raw, unsigned char data);
EXPORT int input_igvsim_featall(raw_t *raw, unsigned char data);

EXPORT int input_oem6f_sol (raw_t *raw, FILE *fp);
EXPORT int input_oem6f_pose(raw_t *raw, FILE *fp);
EXPORT int input_oem6f_raw (raw_t *raw, FILE *fp);
EXPORT int input_oem4f (raw_t *raw, FILE *fp);
EXPORT int input_oem3f (raw_t *raw, FILE *fp);
EXPORT int input_ubxf  (raw_t *raw, FILE *fp);
EXPORT int input_ss2f  (raw_t *raw, FILE *fp);
EXPORT int input_cresf (raw_t *raw, FILE *fp);
EXPORT int input_stqf  (raw_t *raw, FILE *fp);
EXPORT int input_gw10f (raw_t *raw, FILE *fp);
EXPORT int input_javadf(raw_t *raw, FILE *fp);
EXPORT int input_nvsf  (raw_t *raw, FILE *fp);
EXPORT int input_bnxf  (raw_t *raw, FILE *fp);
EXPORT int input_rt17f (raw_t *raw, FILE *fp);
EXPORT int input_sbff  (raw_t *raw, FILE *fp);
EXPORT int input_cmrf  (raw_t *raw, FILE *fp);
EXPORT int input_tersusf(raw_t *raw, FILE *fp);
EXPORT int input_lexrf  (raw_t *raw, FILE *fp);
EXPORT int input_gsoff  (raw_t *raw, FILE *fp);
EXPORT int input_m39f   (raw_t *raw, FILE *fp);
EXPORT int input_ubxm8f (raw_t *raw, FILE *fp);
EXPORT int input_ubxsolf (raw_t *raw,FILE *fp);
EXPORT int input_m39_mixf(raw_t *raw,FILE *fp);
EXPORT int input_sbpjsonf(raw_t *raw, FILE *fp);
EXPORT int input_sbpf(raw_t *raw, FILE *fp);
EXPORT int input_cnavf(raw_t *raw, FILE *fp);
EXPORT int input_igvsim_featf(raw_t  *raw,FILE *fp);

EXPORT int gen_ubx (const char *msg, unsigned char *buff);
EXPORT int gen_stq (const char *msg, unsigned char *buff);
EXPORT int gen_nvs (const char *msg, unsigned char *buff);
EXPORT int gen_lexr(const char *msg, unsigned char *buff);

/* rtcm functions ------------------------------------------------------------*/
EXPORT int init_rtcm   (rtcm_t *rtcm);
EXPORT void free_rtcm  (rtcm_t *rtcm);
EXPORT int input_rtcm2 (rtcm_t *rtcm, unsigned char data);
EXPORT int input_rtcm3 (rtcm_t *rtcm, unsigned char data);
EXPORT int input_rtcm2f(rtcm_t *rtcm, FILE *fp);
EXPORT int input_rtcm3f(rtcm_t *rtcm, FILE *fp);
EXPORT int gen_rtcm2   (rtcm_t *rtcm, int type, int sync);
EXPORT int gen_rtcm3   (rtcm_t *rtcm, int type, int sync);

/* solution functions --------------------------------------------------------*/
EXPORT void initsolbuf(solbuf_t *solbuf, int cyclic, int nmax);
EXPORT void freesolbuf(solbuf_t *solbuf);
EXPORT void freesolstatbuf(solstatbuf_t *solstatbuf);
EXPORT int addpose(posebuf_t *posebuf, const pose_meas_t *pose);
EXPORT void initposebuf(posebuf_t *posebuf,int cyclic, int nmax);
EXPORT void freeposebuf(posebuf_t *posebuf);
EXPORT sol_t *getsol(solbuf_t *solbuf, int index);
EXPORT int addsol(solbuf_t *solbuf, const sol_t *sol);
EXPORT int readsol (char *files[], int nfile, solbuf_t *sol);
EXPORT int readsolt(char *files[], int nfile, gtime_t ts, gtime_t te,
                    double tint, int qflag, solbuf_t *sol);
EXPORT int readsolx(char *files[], int nfile, gtime_t ts, gtime_t te,
                    double tint, int qflag, const solopt_t *opt,
                    solbuf_t *solbuf);
EXPORT int readsolstat(char *files[], int nfile, solstatbuf_t *statbuf);
EXPORT int readsolstatt(char *files[], int nfile, gtime_t ts, gtime_t te,
                        double tint, solstatbuf_t *statbuf);
EXPORT int read_gt_sols(const char *file,solbuf_t *solbuf,int format);
EXPORT int findgtsols(solbuf_t *solbuf,gtime_t time);
EXPORT int inputsol(unsigned char data, gtime_t ts, gtime_t te, double tint,
                    int qflag, const solopt_t *opt, solbuf_t *solbuf);
EXPORT int inputsolx(unsigned char data, gtime_t ts, gtime_t te, double tint,
                     int qflag, solbuf_t *solbuf);
EXPORT int outgroundtruth(unsigned char *buff,const sol_t *sol,
                          const solopt_t *opt);
EXPORT int outprcopts(unsigned char *buff, const prcopt_t *opt);
EXPORT int outsolheads(unsigned char *buff, const solopt_t *opt);
EXPORT int outsols  (unsigned char *buff, const sol_t *sol, const double *rb,
                     const solopt_t *opt,const insstate_t *ins,const insopt_t *insopt,
                     int type);
EXPORT int outsolexs(unsigned char *buff, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt,int type);
EXPORT void outprcopt(FILE *fp, const prcopt_t *opt);
EXPORT void outsolhead(FILE *fp, const solopt_t *opt,const insopt_t *insopt);
EXPORT void outsol  (FILE *fp, const sol_t *sol, const double *rb,
                     const solopt_t *opt,const insstate_t *ins,const insopt_t *insopt);
EXPORT void outsolex(FILE *fp, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt);
EXPORT int outnmea_rmc(unsigned char *buff, const sol_t *sol);
EXPORT int outnmea_gga(unsigned char *buff, const sol_t *sol);
EXPORT int outnmea_gsa(unsigned char *buff, const sol_t *sol,
                       const ssat_t *ssat);
EXPORT int outnmea_gsv(unsigned char *buff, const sol_t *sol,
                       const ssat_t *ssat);

EXPORT int outrnxobsh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxobsb(FILE *fp, const rnxopt_t *opt, const obsd_t *obs, int n,
                      int epflag);
EXPORT int outrnxnavh (FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxgnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxhnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxlnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxqnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxcnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxinavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxnavb (FILE *fp, const rnxopt_t *opt, const eph_t *eph);
EXPORT int outrnxgnavb(FILE *fp, const rnxopt_t *opt, const geph_t *geph);
EXPORT int outrnxhnavb(FILE *fp, const rnxopt_t *opt, const seph_t *seph);

/* google earth kml converter ------------------------------------------------*/
EXPORT int convkml(const char *infile, const char *outfile, gtime_t ts,
                   gtime_t te, double tint, int qflg, double *offset,
                   int tcolor, int pcolor, int outalt, int outtime);
EXPORT int savekml(const char *file, const solbuf_t *solbuf, int tcolor,int pcolor,
                   int outalt, int outtime);

/* gpx converter -------------------------------------------------------------*/
EXPORT int convgpx(const char *infile, const char *outfile, gtime_t ts,
                   gtime_t te, double tint, int qflg, double *offset,
                   int outtrk, int outpnt, int outalt, int outtime);

/* sbas functions ------------------------------------------------------------*/
EXPORT int  sbsreadmsg (const char *file, int sel, sbs_t *sbs);
EXPORT int  sbsreadmsgt(const char *file, int sel, gtime_t ts, gtime_t te,
                        sbs_t *sbs);
EXPORT void sbsoutmsg(FILE *fp, sbsmsg_t *sbsmsg);
EXPORT int  sbsdecodemsg(gtime_t time, int prn, const unsigned int *words,
                         sbsmsg_t *sbsmsg);
EXPORT int sbsupdatecorr(const sbsmsg_t *msg, nav_t *nav);
EXPORT int sbssatcorr(gtime_t time, int sat, const nav_t *nav, double *rs,
                      double *dts, double *var);
EXPORT int sbsioncorr(gtime_t time, const nav_t *nav, const double *pos,
                      const double *azel, double *delay, double *var);
EXPORT double sbstropcorr(gtime_t time, const double *pos, const double *azel,
                          double *var);

/* options functions ---------------------------------------------------------*/
EXPORT opt_t *searchopt(const char *name, const opt_t *opts);
EXPORT int str2opt(opt_t *opt, const char *str);
EXPORT int opt2str(const opt_t *opt, char *str);
EXPORT int opt2buf(const opt_t *opt, char *buff);
EXPORT int loadopts(const char *file, opt_t *opts);
EXPORT int saveopts(const char *file, const char *mode, const char *comment,
                    const opt_t *opts);
EXPORT void resetsysopts(void);
EXPORT void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);
EXPORT void setsysopts(const prcopt_t *popt, const solopt_t *sopt,
                       const filopt_t *fopt);
EXPORT void rtk_debug_opt(prcopt_t *prcopt_,solopt_t *solopt_,
                          filopt_t *filopt_);

/* stream data input and output functions ------------------------------------*/
EXPORT void strinitcom(void);
EXPORT void strinit  (stream_t *stream);
EXPORT void strlock  (stream_t *stream);
EXPORT void strunlock(stream_t *stream);
EXPORT int  stropen  (stream_t *stream, int type, int mode, const char *path);
EXPORT void strclose (stream_t *stream);
EXPORT int  strread  (stream_t *stream, unsigned char *buff, int n);
EXPORT int  strwrite (stream_t *stream, unsigned char *buff, int n);
EXPORT void strsync  (stream_t *stream1, stream_t *stream2);
EXPORT int  strstat  (stream_t *stream, char *msg);
EXPORT int  strstatx (stream_t *stream, char *msg);
EXPORT void strsum   (stream_t *stream, int *inb, int *inr, int *outb, int *outr);
EXPORT int  strgetsel(stream_t *stream, char *sel);
EXPORT int  strsetsel(stream_t *stream, const char *sel);
EXPORT int  strsetsrctbl(stream_t *stream, const char *file);
EXPORT void strsetopt(const int *opt);
EXPORT gtime_t strgettime(stream_t *stream);
EXPORT void strsendnmea(stream_t *stream, const sol_t *sol);
EXPORT void strsendcmd(stream_t *stream, const char *cmd);
EXPORT void strsettimeout(stream_t *stream, int toinact, int tirecon);
EXPORT void strsetdir(const char *dir);
EXPORT void strsetproxy(const char *addr);

/* integer ambiguity resolution ----------------------------------------------*/
EXPORT int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s);
EXPORT int lambda_reduction(int n, const double *Q, double *Z);
EXPORT int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s);
EXPORT int bootstrap(int n,const double *a, const double *Q, double *F,double *Ps);
EXPORT int plambda(const double *a,const double *Qa,int n,int m,double *F,
                   double *s,double p0);
/* standard positioning ------------------------------------------------------*/
EXPORT int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, insstate_t *ins,double *azel,
                  ssat_t *ssat, char *msg);

/* precise positioning -------------------------------------------------------*/
EXPORT void rtkinit(rtk_t *rtk, const prcopt_t *opt);
EXPORT void rtkfree(rtk_t *rtk);
EXPORT int  rtkpos (rtk_t *rtk, const obsd_t *obs, int nobs, const nav_t *nav);
EXPORT int  rtkopenstat(const char *file, int level);
EXPORT void rtkclosestat(void);
EXPORT int  rtkoutstat(rtk_t *rtk, char *buff);
EXPORT int  insoutstat(rtk_t *rtk, char *buff);
EXPORT void initP(int is,int ni,int nx,double unc,double unc0,double *P0);
EXPORT void initx(rtk_t *rtk, double xi, double var, int i);
EXPORT void insinitx(insstate_t *ins,double xi,double var,int i);
EXPORT int insinirtobs(rtksvr_t *svr,const obsd_t *obs,int n,const imud_t *imu);
/* precise point positioning -------------------------------------------------*/
EXPORT void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);
EXPORT int pppnx(const prcopt_t *opt);
EXPORT int pppoutstat(rtk_t *rtk, char *buff);

EXPORT int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                  const nav_t *nav, const double *azel, double *x, double *P);

EXPORT int pppcorr_read(pppcorr_t *corr, const char *file);
EXPORT void pppcorr_free(pppcorr_t *corr);
EXPORT int pppcorr_trop(const pppcorr_t *corr, gtime_t time, const double *pos,
                        double *ztd, double *std);
EXPORT int pppcorr_stec(const pppcorr_t *corr, gtime_t time, const double *pos,
                        double *ion, double *std);

/* stream server functions ---------------------------------------------------*/
EXPORT void strsvrinit (strsvr_t *svr, int nout);
EXPORT int  strsvrstart(strsvr_t *svr, int *opts, int *strs, char **paths,
                        strconv_t **conv, char **cmds, char **cmds_priodic,
                        const double *nmeapos);
EXPORT void strsvrstop (strsvr_t *svr, char **cmds);
EXPORT void strsvrstat (strsvr_t *svr, int *stat, int *byte, int *bps, char *msg);
EXPORT strconv_t *strconvnew(int itype, int otype, const char *msgs, int staid,
                             int stasel, const char *opt);
EXPORT void strconvfree(strconv_t *conv);
EXPORT void strsvrsetsrctbl(strsvr_t *svr, const char *file);

/* gis data functions --------------------------------------------------------*/
EXPORT int gis_read(const char *file, gis_t *gis, int layer);
EXPORT void gis_free(gis_t *gis);

/* qzss lex functions --------------------------------------------------------*/
EXPORT int lexupdatecorr(const lexmsg_t *msg, nav_t *nav, gtime_t *tof);
EXPORT int lexreadmsg(const char *file, int sel, lex_t *lex);
EXPORT void lexoutmsg(FILE *fp, const lexmsg_t *msg);
EXPORT int lexconvbin(int type, int format, const char *infile,
                      const char *outfile);
EXPORT int lexeph2pos(gtime_t time, int sat, const nav_t *nav, double *rs,
                      double *dts, double *var);
EXPORT int lexioncorr(gtime_t time, const nav_t *nav, const double *pos,
                      const double *azel, double *delay, double *var);
/* ins functions--------------------------------------------------------------*/
EXPORT void matmul3(const char *tr, const double *A, const double *B, double *C);
EXPORT void matmul3v(const char *tr, const double *A, const double *b, double *c);
EXPORT void seteye(double* A,int n);
EXPORT void setzero(double *A,int n,int m);
EXPORT void skewsym3(const double *ang, double *C);
EXPORT void skewsym3x(double x,double y,double z,double *C);
EXPORT int resize(double **A,int m,int n,int p,int q);
EXPORT double gravity0(const double *pos);
EXPORT void gravity(const double *re, double *ge);
EXPORT void initins(insstate_t *ins, const double *re, double angh,
                    const imud_t *data, int n,const insopt_t *opt);
EXPORT int updateins(const insopt_t *insopt,insstate_t *ins, const imud_t *data);
EXPORT int updateinsbe(const insopt_t *insopt,insstate_t *ins,const imud_t *data);
EXPORT int updateinsbn(const insopt_t *insopt,insstate_t *ins,const imud_t *data);
EXPORT int updateinsn(const insopt_t *insopt,insstate_t *ins,const imud_t *data);
EXPORT int updateinsb(const insopt_t *insopt,insstate_t *ins,const imud_t *data);
EXPORT void correctatt(const double *dphi,const double *C,double *Cc);

EXPORT int kinematicsecef(const double dt,const double *Cbe1,const double *Cbe0,
                          const double *ve1,const double *ve0,
                          const double *re1,
                          double *fb,double *omgb);
EXPORT int simimumeas(double *fb,double *omgb,const imu_err_t *err,double dt);
EXPORT int generatepath(const char *file,const cam_t *cam,const imu_err_t *err,
                        const insopt_t *opt,
                        const char *imufile,
                        const char *gpsfile,const char *vofile);
EXPORT int sim_imu_static(const double *rpy,const double *pos,const double ts,
                          const double T,const imu_err_t *err,imu_t *data);

EXPORT void traceins(int level, const insstate_t *ins);
EXPORT void traceinss(int level, const double *Cbe,const double *re,
                      const double *ve,gtime_t time);
EXPORT void ned2xyz(const double *pos,double *Cne);
EXPORT void rpy2dcm(const double *rpy,double *Cnb);
EXPORT void dcm2rpy(const double *Cnb,double *rpy);
EXPORT void rov2dcm(const double *rv,double *C);
EXPORT void rov2qua(const double *rv,quat_t *q);
EXPORT void quat2rpy(const quat_t *quat,double *rpy);

EXPORT void quat2rot(const double *q,double *v);
EXPORT void quatrot(const double *qab,double *va,int dir,double *vb);
EXPORT void rv2quat(const double *rv,double *q);
EXPORT void rot2dcm(const double *rot,double *dcm);
EXPORT void quat2dcmx(const double *qba, double *Cba);
EXPORT void dcm2quatx(const double *dcm, double* quat);
EXPORT void quatmulx(const double *qab,const double *qca,double *qcb);

EXPORT void normdcm(double *C);
EXPORT void quatupd(const double *vi,quat_t *qo);
EXPORT void qmulv(const double *vi,const quat_t *quat,double *vo);
EXPORT void rfu2frd(const double *rfu,double *frd,double *C);
EXPORT void enu2ned(const double *enu,double *ned,double *C);
EXPORT void dcm2quat(const double *C,quat_t *q);
EXPORT void quat2dcm(const quat_t *q,double *C);
EXPORT void ins_errmodel(const double *accl,const double *gyro,double *cor_accl,
                         double *cor_gyro,insstate_t* ins);
EXPORT void rotscull_corr(insstate_t *ins,const insopt_t *opt,double dt,
                          double *dv,double *da);
EXPORT void rp2head(const double roll,const double picth,const double *gyro,double *head);
EXPORT void getatt(const insstate_t *ins,double *rpy);
EXPORT double NORMANG(double ang);
EXPORT void pregrav(const double *pos, double *g);
EXPORT int  insupdate_ned(const insopt_t *insopt,insstate_t *ins, const imud_t *data);
EXPORT void radii(const double *rn,double *R_N,double *R_E);
EXPORT void gravity_ned(const double *pos,double *gn);
EXPORT void estatt(const imud_t *data,int n,double *Cbn);
EXPORT double georadi(const double *pos);
EXPORT void rmlever(const double *pos, const double *re, const double *ve,
                    const double *lever, const double *Cbe, const double *omgb,
                    double *rec, double *vec);
EXPORT void gapv2ipv(const double *pos,const double *vel,const double *Cbe,
                     const double *lever,const imud_t *imu,double *posi,double *veli);
EXPORT void insp2antp(const insstate_t *ins,double *rr);
EXPORT int ant2inins(gtime_t time,const double *rr,const double *vr,
                     const insopt_t *opt,const imu_t *imu,insstate_t *ins,int *iimu);
EXPORT void propinss(insstate_t *ins,const insopt_t *opt,double dt,
                     double *x,double *P);
EXPORT void getPhi(const insopt_t *opt,const insstate_t *ins,double *Phi);
EXPORT void getQn(const insopt_t *opt,const insstate_t *ins,double *Qn);
EXPORT void getP0(const insopt_t *opt,double *P0);
EXPORT int lcigpos(const insopt_t *opt, const imud_t *data, insstate_t *ins,
                   gmea_t *gnss, int upd);
EXPORT int tcigpos(const prcopt_t *opt,const obsd_t *obs,int n,const nav_t *nav,
                   const imud_t *imu,rtk_t *rtk,insstate_t *ins,int upd);

EXPORT int lcigposrb(const insopt_t *opt, const imud_t *data, insstate_t *ins,
                     gmea_t *gnss, int upd);

EXPORT int doppler(const obsd_t *obs,int n,const nav_t *nav,const prcopt_t *opt,insstate_t *ins);
EXPORT void initlc(insopt_t *insopt, insstate_t *ins);
EXPORT void freelc(insstate_t *ins);
EXPORT void freeins(insstate_t *ins);
EXPORT void getaccl(const double *fib,const double *Cbe,const double *re,
                    const double *ve,double *ae);
EXPORT void cnscl(const double *gyro,const double *accl,double *phim,double *dvbm,
                  const double *wm0,const double *vm0,double *dphim,double *rotm,
                  double *scullm);
EXPORT void ins_errmodel2(const double *accl,const double *gyro,const double *Ma,
                          const double *Mg,const double *ba,const double *bg,
                          const double *Gg,double *cor_accl,double *cor_gyro);
EXPORT int coarse_align(insstate_t *ins,const imud_t *data,int n,const insopt_t *opt);
EXPORT int fine_align(insstate_t *ins,const imud_t *data,int n,const insopt_t *opt);
EXPORT int fine_alignex(insstate_t *ins,const imud_t *data,int n,const insopt_t *opt);
EXPORT int alignfnex(const double *imu,const double *q0,const double *pos,
                     const double *phi0,const ins_align_t *pasi,double dt,
                     int n,double *att0,double *qo);
EXPORT int  alignvnex(const double *imu,int n,double dt,const double *q0,
                      const double *pos,const double *phi0,const double *wvn,
                      const ins_align_t *pas,double *att0,double *qo);
EXPORT int fine_align_lym(insstate_t *ins,const imud_t *data,int n,const insopt_t *opt);
EXPORT int readimu(const char *file, imu_t *imu,int decfmt,int format,int coor,int valfmt);
EXPORT int readstim300(const char *file,imu_t *imu);
EXPORT int sortimudata(imu_t *imu);
EXPORT void adjimudata(const prcopt_t *opt,imu_t *imu);
EXPORT void adjustimu(const prcopt_t *opt,imud_t *imu);
EXPORT int addgmea(gmeas_t *gmeas, const gmea_t *data);
EXPORT int savegmeas(insstate_t *ins,const sol_t *sol,const gmea_t *gmea);
EXPORT int chksdri(const double *vel,int n);
EXPORT int rechkatt(insstate_t *ins,const imud_t *imu);

EXPORT int get_m39_img(const char *dir,const uint64_t fts,uint64_t ftns,
                       char *imgfile);

EXPORT int detstaticw(const insstate_t *ins,const imud_t *imu,int n,const insopt_t *opt);
EXPORT int detstatic_GLRT(const imud_t *imu,int n,const insopt_t *opt,const double *pos);
EXPORT int detstatic_MV(const imud_t *imu,int n,const insopt_t *opt);
EXPORT int detstatic_MAG(const imud_t *imu,int n,const insopt_t *opt,const double *pos);
EXPORT int detstatic_ARE(const imud_t *imu,int n,const insopt_t *opt);
EXPORT int detstatic_ODO(const insopt_t *opt,const odod_t *odo);
EXPORT int detstc(const imud_t *imu,int n,const insopt_t *opt,const double *pos);
EXPORT double vel2head(const double *vel);
EXPORT void corratt(const double *dx,double *C);
EXPORT int insinitrt(rtksvr_t *svr,const sol_t *sol,const imud_t *imu);
EXPORT void insstatocov(const insopt_t *opt,const insstate_t *ins,double *Pa,
                        double *Pv,double *Pp);
EXPORT void ins2sol(insstate_t *ins,const insopt_t *opt,sol_t *sol);
EXPORT void rvec2quat(const double *rvec,double *q);
EXPORT void dcm2rot(const double *C,double *rv);
EXPORT void getvn(const insstate_t *ins,double *vn);
EXPORT void updinsn(insstate_t *ins);
EXPORT void updinse(insstate_t *ins);
EXPORT int insinitdualant(rtksvr_t *svr,const pose_meas_t *pose,const sol_t *sol,
                          const imud_t *imu);
EXPORT int insinitgiven(rtksvr_t *svr,const imud_t *imu);
EXPORT int calibdualant(const insopt_t *opt,const solbuf_t *solbuf,
                        posebuf_t *posebuf,gtime_t ts,gtime_t te,
                        double *Cvb,double *ave,double *std);
EXPORT int calibdualantx(const insopt_t *opt,const solbuf_t *solbuf,
                         posebuf_t *posebuf,gtime_t ts,gtime_t te,
                         double *Cvb,double *std);

/* close-loop for states ----------------------------------------------------*/
EXPORT void lcclp(double *x,double *Cbe,double *re,double *ve,double *fib,
                  double *omgb,double *Gg, double *rec,double *vec,double *aec,
                  double *bac, double *bgc,double *Mac,double *Mgc,
                  double *leverc,double *Cbec,double *fibc,double *omgbc,
                  const insopt_t *opt);

/* motion constraint---------------------------------------------------------*/
EXPORT int nhc(insstate_t *ins,const insopt_t *opt,const imud_t *imu);
EXPORT int zvu(insstate_t *ins,const insopt_t *opt,const imud_t *imu,int flag);
EXPORT int zaru(insstate_t *ins,const insopt_t *opt,const imud_t *imu,int flag);
EXPORT void clp(insstate_t *ins,const insopt_t *opt,const double *x);

/* odometry options-----------------------------------------------------------*/
EXPORT void odores(double rese);
EXPORT void odod(double de);
EXPORT void vr2vb(const double *vr,const double *lever,const double *Cbe,
                  const imud_t *imu,const double s,double *vb);
EXPORT int ve2vr(const double *ve,const double *Cbe,const double *Cbr,
                 const double *lever,const imud_t *imu,const double s,
                 double *vr);
EXPORT int odo(const insopt_t *opt,const imud_t *imu,const odod_t *odo,
               insstate_t *ins);
EXPORT void initodo(const odopt_t *opt,insstate_t *ins);

/* magnetic heading-----------------------------------------------------------*/
EXPORT int magmodel(const char *file);
EXPORT double maghead(const magopt_t *opt,const mag_t *data,
                      const double *Cbn,
                      const double *pos);
EXPORT int magnetometer(insstate_t *ins,const insopt_t *opt,
                        const mag_t *data);

/* quaternion common function-------------------------------------------------*/
EXPORT void quat_init(quat_t *quat, const vec3_t *acc, const vec3_t *mag);
EXPORT void quat_init_axis(quat_t *q, double x, double y, double z, double a);
EXPORT void quat_init_axis_v(quat_t *q, const vec3_t *v, double a);
EXPORT void quat_to_axis(const quat_t *q, double *x, double *y, double *z, double *a);
EXPORT void quat_to_axis_v(const quat_t *q, vec3_t *v, double *a);
EXPORT void quat_rot_vec(vec3_t *vo, const vec3_t *vi, const quat_t *q);
EXPORT void quat_rot_vec_self(vec3_t *v, const quat_t *q);
EXPORT double quat_len(const quat_t *q);
EXPORT void quat_copy(quat_t *qo, const quat_t *qi);
EXPORT void quat_scale(quat_t *qo, const quat_t *qi, double f);
EXPORT void quat_conj(quat_t *qo, const quat_t *qi);
EXPORT void quat_add(quat_t *qo, const quat_t *q1, const quat_t *q2);
EXPORT void quat_mul(quat_t *o, const quat_t *q1, const quat_t *q2);
EXPORT void quat_normalize(quat_t *qo, const quat_t *qi);
EXPORT void quat_normalize_self(quat_t *q);
EXPORT void quat_to_euler(euler_t *e, const quat_t *q);
EXPORT quat_t *quat_apply_relative_yaw_pitch_roll(quat_t *q,double yaw, double pitch, double roll);
EXPORT void quat_to_rh_rot_matrix(const quat_t *q, double *m);
EXPORT vec3_t *vec3_rot_axis_self(vec3_t *vo, double x, double y, double z, double angle);
EXPORT void quat_from_u2v(quat_t *q, const vec3_t *u, const vec3_t *v, const vec3_t *up);
EXPORT double quat_dot(const quat_t *q1, const quat_t *q2);
EXPORT quat_t *quat_nlerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto, double t);
EXPORT quat_t *quat_slerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto, double t);
EXPORT void quat_inv(const quat_t *q,quat_t *qi);
EXPORT void euler_to_quat(const euler_t *euler,quat_t *q);
EXPORT void rpy2quat(const double *rpy,quat_t *q);
EXPORT quat_t invquat(const quat_t *q);

/* ins common functions in enu-rfu-frame--------------------------------------*/
EXPORT void q2mat(const double *qnb,double *C);
EXPORT void qdelphi(const double *phi,double *qnb);
EXPORT void rv2q(const double *rv,double *q);
EXPORT void q2att(const double *q,double *att);
EXPORT void qupdt(const double *rv,double *q);
EXPORT void qmul(const double *q1,const double *q2,double *q);
EXPORT void qmulve(const double *q,const double *vi,double *vo);
EXPORT void qconj(const double *qi,double *qo);
EXPORT void rv2m(const double *rv,double *m);
EXPORT void a2qua(const double *att,double *q);

/* so3 lie group functions----------------------------------------------------*/
EXPORT void so3_exp(const double *omega,double *C);
EXPORT void so3_log(const double *C,double *omega,double *theta);
EXPORT void liebracket(const double *omega1,const double *omega2,double *omega);
EXPORT void so3_hat(const double *omega,double *Omg);
EXPORT void so3_vee(const double *Omg,double *omega);
EXPORT void expmat(const double *A,int n,double *E);
EXPORT void rpy2c(const double *rpy,double *C);
EXPORT void c2rpy(const double *C,double *rpy);
EXPORT void addlie(const double *a1,const double *a2,double *a);
EXPORT void so3_jac(const double *phi,double *Jr,double *Jl);
EXPORT void so3jac(const double *phi,double *Jri);
extern void so3exp(const double *phi,double *C);
extern void so3log(const double *C,double *phi);

/* se3 lie group functions----------------------------------------------------*/
EXPORT void se3_log(const double *R,const double *t,double *omega);
EXPORT void se3_exp(const double *omg,double *R,double *t);
EXPORT void se3_hat(const double *omg,double *Omg);
EXPORT void se3_vee(const double *Omg,double *omg);

/* pose fusion for ins navigation---------------------------------------------*/
EXPORT int posefusion(const insopt_t *opt,const pose_meas_t *data,
                      insstate_t *ins,int flag);

/* RTS functions--------------------------------------------------------------*/
EXPORT int lcrts(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                 const solopt_t *solopt,int port,const char *file);
EXPORT int lcfbsm(const imu_t *imu,const gsof_data_t *pos,const prcopt_t *popt,
                  const solopt_t *solopt,int port,const char *file);
EXPORT void set_fwd_soltmp_file(const char *file);
EXPORT void set_fwdtmp_file(const char *file);
EXPORT int bckupinsinfo(insstate_t *ins, const insopt_t *opt, int type);
EXPORT int tcpostpos(prcopt_t *popt, const solopt_t *solopt, int port,
                     const char *outfile, char** infiles);

/* estimated states interface-------------------------------------------------*/
EXPORT int xnP (const insopt_t *opt);
EXPORT int xnV (const insopt_t *opt);
EXPORT int xnA (const insopt_t *opt);
EXPORT int xnBa(const insopt_t *opt);
EXPORT int xnBg(const insopt_t *opt);
EXPORT int xnDt(const insopt_t *opt);
EXPORT int xnSg(const insopt_t *opt);
EXPORT int xnSa(const insopt_t *opt);
EXPORT int xnRg(const insopt_t *opt);
EXPORT int xnRa(const insopt_t *opt);
EXPORT int xnLa(const insopt_t *opt);
EXPORT int xnOs(const insopt_t *opt);
EXPORT int xnOl(const insopt_t *opt);
EXPORT int xnOa(const insopt_t *opt);
EXPORT int xnCl(const insopt_t *opt);
EXPORT int xnIF(const insopt_t *opt);
EXPORT int xnRc(const insopt_t *opt);
EXPORT int xnRr(const insopt_t *opt);
EXPORT int xnI (const insopt_t *opt);
EXPORT int xnT (const insopt_t *opt);
EXPORT int xnL (const insopt_t *opt);
EXPORT int xnD (const insopt_t *opt);
EXPORT int xnB (const insopt_t *opt);
EXPORT int xnRx(const insopt_t *opt);
EXPORT int xnCm(const insopt_t *opt);
EXPORT int xnCla(const insopt_t *opt);
EXPORT int xnCfo(const insopt_t *opt);
EXPORT int xnCkp(const insopt_t *opt);
EXPORT int xnVm(const insopt_t *opt);
EXPORT int xiVm(const insopt_t *opt);
EXPORT int xnX (const insopt_t *opt);
EXPORT int xiA (const insopt_t *opt);
EXPORT int xiV (const insopt_t *opt);
EXPORT int xiP (const insopt_t *opt);
EXPORT int xiBa(const insopt_t *opt);
EXPORT int xiBg(const insopt_t *opt);
EXPORT int xiDt(const insopt_t *opt);
EXPORT int xiSg(const insopt_t *opt);
EXPORT int xiSa(const insopt_t *opt);
EXPORT int xiRg(const insopt_t *opt);
EXPORT int xiRa(const insopt_t *opt);
EXPORT int xiLa(const insopt_t *opt);
EXPORT int xiOs(const insopt_t *opt);
EXPORT int xiOl(const insopt_t *opt);
EXPORT int xiOa(const insopt_t *opt);
EXPORT int xiCm(const insopt_t *opt);
EXPORT int xiCl(const insopt_t *opt);
EXPORT int xiCfo(const insopt_t *opt);
EXPORT int xiCkp(const insopt_t *opt);
EXPORT int xiRc(const insopt_t *opt);
EXPORT int xiRr(const insopt_t *opt);
EXPORT int xiIo(const insopt_t *opt,int s);
EXPORT int xiTr(const insopt_t *opt,int r);
EXPORT int xiLl(const insopt_t *opt,int f);
EXPORT int xiDl(const insopt_t *opt);
EXPORT int xiBs(const insopt_t *opt,int s,int f);
EXPORT int unusex(const insopt_t *opt,int index,insstate_t *ins,double *x);

/* application defined functions ---------------------------------------------*/
EXPORT int showmsg(char *format,...);
EXPORT void settspan(gtime_t ts, gtime_t te);
EXPORT void settime(gtime_t time);
EXPORT void none(void);

/* image filter functions-----------------------------------------------------*/
EXPORT void blob5x5(const uint8_t* in, int16_t* out, int w, int h);
EXPORT void checkerboard5x5(const uint8_t* in, int16_t* out, int w, int h);
EXPORT void sobel5x5(const uint8_t* in, uint8_t* out_v, uint8_t* out_h,
                     int w, int h);
EXPORT int initimg(img_t *data,int w,int h,gtime_t time);
EXPORT void freeimg(img_t *data);

/* match feature points functions---------------------------------------------*/
EXPORT int matchfeats(match_t *pmatch,const img_t *img);
EXPORT void init_match(match_t *match,const matchopt_t *opt);
EXPORT void free_match(match_t *match);
EXPORT void free_match_set(match_set_t *mset);

/* visual odometry estimator--------------------------------------------------*/
EXPORT int estmonort(const voopt_t *opt,const match_set_t *feat,double *Tr,double *ratio);

/* pgm file read/write--------------------------------------------------------*/
EXPORT void ppmWriteFileRGB(char *fname, unsigned char *redimg,
                            unsigned char *greenimg,unsigned char *blueimg,
                            int ncols, int nrows);
EXPORT void ppmWrite(FILE *fp,unsigned char *redimg,unsigned char *greenimg,
                     unsigned char *blueimg,
                     int ncols, int nrows);
EXPORT void pgmWriteFile(char *fname, unsigned char *img, int ncols,
                         int nrows);
EXPORT void pgmWrite(FILE *fp,unsigned char *img, int ncols, int nrows);
EXPORT unsigned char* pgmReadFile(char *fname,unsigned char *img,
                                  int *ncols, int *nrows);

EXPORT unsigned char* pgmRead(FILE *fp,unsigned char *img,
                              int *ncols, int *nrows);

/* jpeg file read/write------------------------------------------------------*/
EXPORT void savejpg(const char* filename, uchar* body, int h, int w, int ch,
                    int quality);
EXPORT int loadjpg(const char* filename, uchar* &body, int &h, int &w,
                   int &ch);
EXPORT int readjpeg(const char *imgfile,gtime_t time,img_t *img,int flag);

/* klt track-----------------------------------------------------------------*/
EXPORT void initklt();
EXPORT void freeklt();
EXPORT int kltstatus(match_point_t *matchp,const img_t *pimg,const img_t *cimg,
                     const voopt_t *opt);

/* tracking functions--------------------------------------------------------*/
EXPORT void tracetrack(const track_t *track);
EXPORT void feat2ppm(feature *featurelist, int nfeature,unsigned char *greyimg,
                     int ncols,int nrows,char *filename);
EXPORT void freetrackset(track_t *track);
EXPORT void freetrack(trackd_t *track);
EXPORT int match2track(const match_set *mset,gtime_t tp,gtime_t tc,int curr_frame,
                       const img_t *pimg,const img_t *cimg,
                       const voopt_t *opt,track_t *track);
EXPORT int inittrack(trackd_t *data,const voopt_t *opt);
EXPORT void drawtrack(const track_t *track,const voopt_t *opt);
EXPORT void drawtrackd(const trackd_t *track,const voopt_t *opt);
EXPORT trackd_t *gettrack(const track_t *track,int uid);
EXPORT int outofview(const track_t *track,const voopt_t *opt,
                     int id,gtime_t time);
EXPORT void inittrackimgbuf(const voopt_t *opt);
EXPORT void freetrackimgbuf();
EXPORT img_t* getimgdata(gtime_t time);
EXPORT int getgmsmatches(const voopt_t *opt,match_set *match,const img_t *Ip,const img_t *Ic);

/* pnp pose estimate function------------------------------------------------*/
EXPORT int p3pthree(const feature *feats,int nf,double *xp,int np,
                    const voopt_t *opt,double *R,double *t);
EXPORT int p3p(const feature *feats,int nf,const double *xp,int np,
               const voopt_t *opt,double *R,double *t);

/* visual odometry utils-----------------------------------------------------*/
EXPORT int iscolinear(const double *p1,const double *p2,const double *p3,
                      int dim, int flag);
EXPORT int triangulate3D(const double *C21,const double *t12_1,
                         const double *obs1,const double *obs2,
                         const double *K,double *X);

/* camera distort/undistort function-----------------------------------------*/
EXPORT void distortradtan(const cam_t *cam,const double *in,double *out,
                          double *J);
EXPORT int undistortradtan(const cam_t *cam,const double *in,double *out,
                           double *J);

/* vo aid to ins/gnss coupled functions---------------------------------------*/
EXPORT int voigposmsckf(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                        const img_t *img,int flag);

EXPORT void initvoaid(insopt_t *opt);
EXPORT void freevoaid();
EXPORT void freevoaidlc();
EXPORT void initvoaidlc(insopt_t *opt);
EXPORT int kalibrrosbag(const char *datfile,const char *imufile,
                        const char *imgdir,
                        const char *output);
EXPORT int voigposlc(const insopt_t *opt,insstate_t *ins,const imud_t *imu,
                     const img_t *img,int flag);
EXPORT void updcamposevar(const insopt_t *opt,const double *var,vostate_t *vo);
EXPORT void initcamposevar(const insopt_t *opt,insstate_t *ins,vostate_t *vo);
EXPORT int updatecamera(insstate_t *ins,const insopt_t *opt,const double *dT,gtime_t time);
EXPORT void ins2camera(const insstate_t *ins,double *rc,double *Cce);
EXPORT void camera2ins(const double *Cce,const double *rc,const double *Cbc,
                       const double *lbc,double *Cbe,double *re);
EXPORT int initcampose(insstate_t *ins,const insopt_t *opt);
EXPORT int zvucam(const insopt_t *opt,const double *dT,int status,const imud_t *imu,
                  gtime_t tc,insstate_t *ins);

EXPORT void rt2tf(const double *R,const double *t,double *T);
EXPORT void tf2rt(const double *T,double *R,double *t);

EXPORT int predictfeat(const double *R,const double *t,const double *K,
                       const double *uv,double *uvp);
EXPORT int inroi(const float u,const float v,const matchopt_t *opt);
EXPORT int fastfeats(const unsigned char *img,int img_w,int img_h,short barrier,
                     const matchopt_t *opt,int *uv,int *num);
EXPORT int solveRt(const voopt_t *opt,const insstate_t *ins,const match_set_t *mf,
                   double *dT,double *ratio);
EXPORT int solveRt5p(const voopt_t *opt,const insstate_t *ins,const match_set_t *mf,
                     double *dT,double *ratio);
/* ins-gnss-vo coupled post-processing----------------------------------------*/
EXPORT int igvopostpos(const gtime_t ts, gtime_t te, double ti, double tu,
                       const prcopt_t *popt, const solopt_t *sopt,
                       const filopt_t *fopt, char **infile, int n,int port,
                       const char *outfile);

/* car-vig position server functions------------------------------------------*/
EXPORT int carvigsvrmark(rtksvr_t *svr, const char *name, const char *comment);
EXPORT void rtksvrsstat(rtksvr_t *svr, int *sstat, char *msg);
EXPORT int carvigsvrostat(rtksvr_t *svr, int rcv, gtime_t *time, int *sat,
                          double *az, double *el, int **snr, int *vsat);

EXPORT void carvigsvrclosestr(rtksvr_t *svr, int index);
EXPORT int carvigsvropenstr(rtksvr_t *svr, int index, int str, const char *path,
                            const solopt_t *solopt);

EXPORT int carvigsvrstart(rtksvr_t *svr, int cycle, int buffsize, int *strs,
                          char **paths, int *formats, int navsel, char **cmds,
                          char **cmds_periodic, char **rcvopts, int nmeacycle,
                          int nmeareq, const double *nmeapos, prcopt_t *prcopt,
                          solopt_t *solopt, stream_t *moni, char *errmsg);

EXPORT void carvigsvrstop(rtksvr_t *svr, char **cmds);
EXPORT void rtksvrlock  (rtksvr_t *svr);
EXPORT void rtksvrunlock(rtksvr_t *svr);
EXPORT int carvigsvrinit(rtksvr_t *svr);
EXPORT void carvigsvrfree(rtksvr_t *svr);
EXPORT void carvigsvrsetimgpath(const char *dirpath);

/* virtual console functions--------------------------------------------------*/
EXPORT vt_t *vt_open(int sock, const char *dev);
EXPORT void vt_close(vt_t *vt);
EXPORT int vt_getc(vt_t *vt, char *c);
EXPORT int vt_gets(vt_t *vt, char *buff, int n);
EXPORT int vt_putc(vt_t *vt, char c);
EXPORT int vt_puts(vt_t *vt, const char *buff);
EXPORT int vt_printf(vt_t *vt, const char *format, ...);
EXPORT int vt_chkbrk(vt_t *vt);
EXPORT int vt_openlog(vt_t *vt, const char *file);
EXPORT void vt_closelog(vt_t *vt);

/* interface for opencv library-----------------------------------------------*/
EXPORT void dipsplyimg(const img_t *img);
EXPORT void drawalltrack(const track_t *track);

EXPORT int estvo(const voopt_t *opt,const img_t *img,double *dTr,double *ratio);
EXPORT void resetmonoa();
EXPORT void freemonoa();

#ifdef __cplusplus
}
#endif
#endif 
