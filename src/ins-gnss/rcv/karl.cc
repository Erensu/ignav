/*----------------------------------------------------------------------------
* karl.cc : read the Karlsruhe Dataset
*
* version : $Revision:$ $Date:$
* history : 2017/03/30  1.0  new
*----------------------------------------------------------------------------*/
#include <navlib.h>

/* sync a newline------------------------------------------------------------*/
static int syncnewline(unsigned char *buff, int nb)
{
    if (buff[nb-1]=='\n'||(buff[nb-2]=='\r'&&buff[nb-1]=='\n')) return 1;
    return 0;
}
/* input karlsruhe dataset format image--------------------------------------
 * args  : raw_t *raw         IO  raw struct
 *         unsigned char data I   raw data
 * return: >1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int input_img_karl(raw_t *raw, unsigned char data)
{
    return 11;
}