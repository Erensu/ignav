/*------------------------------------------------------------------------------
* io_png.h : decode png image raw data header file
*
* version : $Revision:$ $Date:$
* history : 2018/03/21  1.0  new
*-----------------------------------------------------------------------------*/
#ifndef _IO_PNG_H
#define _IO_PNG_H

#ifdef __cplusplus
extern "C" {
#endif

#define IO_PNG_VERSION "0.20110825"
#include <stddef.h>

/* io_png.c------------------------------------------------------------------*/
extern char *io_png_info(void);
extern unsigned char *io_png_read_u8(const char *fname, size_t *nxp, size_t *nyp,
                                     size_t *ncp);
extern unsigned char *io_png_read_u8_rgb(const char *fname, size_t *nxp,size_t *nyp);
extern unsigned char *io_png_read_u8_gray(const char *fname,size_t *nxp,size_t *nyp);
extern float *io_png_read_f32(const char *fname,size_t *nxp,size_t *nyp,size_t *ncp);
extern float *io_png_read_f32_rgb(const char *fname, size_t *nxp,size_t *nyp);
extern float *io_png_read_f32_gray(const char *fname,size_t *nxp,size_t *nyp);
extern int io_png_write_u8(const char *fname,const unsigned char *data,size_t nx,
                           size_t ny,size_t nc);
extern int io_png_write_f32(const char *fname,const float *data,size_t nx,size_t ny,
                            size_t nc);

#ifdef __cplusplus
}
#endif

#endif /* !_IO_PNG_H */
