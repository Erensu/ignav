/*------------------------------------------------------------------------------
* odo.c : read odometry measurement data
*
* version : $Revision:$ $Date:$
* history : 2017/11/18  1.0  new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* read odometry measurement data--------------------------------------------
 * args   :  char *file  I   input odometry measurement data file
 *           odo_t *odo  IO  odometry measurement data records
 * return :  1 (ok) or 0 (fail)
 * --------------------------------------------------------------------------*/
extern int readodo(const char *file,odo_t *odo)
{
    return 1;
}