/*---------------------------------------------------------------------------
 * serialisation.h
 *
 *  Created on: 02.02.2013
 *      Author: Thomas Hartmann
 * -------------------------------------------------------------------------*/

#ifndef SERIALISATION_INLINE_H_
#define SERIALISATION_INLINE_H_
/* --------------------------------------------------------------------------
 * \brief This header file implements a number of serialisation functions.
 *
 *  instream is treated as little endian. So for little endian streams use
 *  plain getter/setter function, for bit endian streams use reverse setter/ getter
 *  functions.
 *  If BIG_ENDIAN is defined, the platform is supposed to be big endian and
 *  the appropriate functions will be called to reverse or not reverse byte order.
 *
 *  If UNALIGNED_ACCESS is defined, casts are used for non reversing operations,
 *  otherise memcpy or for loops are used.
 * --------------------------------------------------------------------------*/
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* one byte getter/setter functions */
#define U1(__p)				(*((unsigned char  *)(__p)))
#define I1(__p)				(*((char   *)(__p)))
#define setU1(__p,__v)		((*((unsigned char*)(__p)))=__v)
#define setI1(__p,__v)		(((*((char*)(__p)))=__v)

/* multibyte setter/getter functions */
#ifdef BIG_ENDIAN
	/* instream is treated as little endian, so revert byte order */
	#define U2(__p)				__U2r(__p)
	#define U4(__p)				__U4r(__p)
	#define I2(__p)				__I2r(__p)
	#define I4(__p)				__I4r(__p)
	#define R4(__p)				__R4r(__p)
	#define R8(__p)				__R8r(__p)

	/* plattform is also big endian, so don't revert byte order */
	#define U2r(__p)			__U2(__p)
	#define U4r(__p)			__U4(__p)
	#define I2r(__p)			__I2(__p)
	#define I4r(__p)			__I4(__p)
	#define R4r(__p)			__R4(__p)
	#define R8r(__p)			__R8(__p)

	#define setR8(__p,__v)		__setR8r(__p,__v)
	#define setR4(__p,__v)		__setR4r(__p,__v)
	#define setI4(__p,__v)		__setI4r(__p,__v)
	#define setI2(__p,__v)		__setI2r(__p,__v)
	#define setU4(__p,__v)		__setU4r(__p,__v)
	#define setU2(__p,__v)		__setU2r(__p,__v)

	#define setR8r(__p,__v)  	__setR8(__p,__v)
	#define setR4r(__p,__v)	    __setR4(__p,__v)
	#define setI4r(__p,__v) 	__setI4(__p,__v)
	#define setI2r(__p,__v) 	__setI2(__p,__v)
	#define setU4r(__p,__v) 	__setU4(__p,__v)
	#define setU2r(__p,__v) 	__setU2(__p,__v)
#else
	/* LITTLE ENDIAN  getter fucntions */
	#define U2(__p)				__U2(__p)
	#define U4(__p)				__U4(__p)
	#define I2(__p)				__I2(__p)
	#define I4(__p)				__I4(__p)
	#define R4(__p)				__R4(__p)
	#define R8(__p)				__R8(__p)
	/* reverse get */
	#define U2r(__p)			__U2r(__p)
	#define U4r(__p)			__U4r(__p)
	#define I2r(__p)			__I2r(__p)
	#define I4r(__p)			__I4r(__p)
	#define R4r(__p)			__R4r(__p)
	#define R8r(__p)			__R8r(__p)

	/* LITTLE ENDIAN  setter functions */
	#define setR8(__p,__v)		__setR8(__p, __v)
	#define setR4(__p,__v)		__setR4(__p, __v)
	#define setI4(__p,__v)		__setI4(__p, __v)
	#define setI2(__p,__v)		__setI2(__p, __v)
	#define setU4(__p,__v)		__setU4(__p, __v)
	#define setU2(__p,__v)		__setU2(__p, __v)

    #define setR8r(__p,__v) 	__setR8r(__p, __v)
	#define setR4r(__p,__v) 	__setR4r(__p, __v)
	#define setI4r(__p,__v) 	__setI4r(__p, __v)
	#define setI2r(__p,__v) 	__setI2r(__p, __v)
	#define setU4r(__p,__v) 	__setU4r(__p, __v)
	#define setU2r(__p,__v) 	__setU2r(__p, __v)
#endif

#ifdef UNALIGNED_ACCESS
	#define __U2(p)       (*((unsigned short *)(p)))
	#define __U4(p)       (*((unsigned int   *)(p)))
	#define __I2(p)       (*((short  *)(p)))
	#define __I4(p)       (*((int    *)(p)))
	#define __R4(p)       (*((float  *)(p)))

	static double __R8(const unsigned char * const pa_pnInData)
	{
		double nValue;
		memcpy(&nValue,pa_pnInData,sizeof(nValue));
		return nValue;
	}

	static void __setR8(unsigned char * p_pnDest, double pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}

	#define __setR4(__p,__v)		((*((float*)(__p)))=__v)
	#define __setI4(__p,__v)		((*((int32_t*)(__p)))=__v)
	#define __setI2(__p,__v)		((*((int16_t*)(__p)))=__v)
	#define __setU4(__p,__v)		((*((uint32_t*)(__p)))=__v)
	#define __setU2(__p,__v)		((*((uint16_t*)(__p)))=__v)
#else
	static uint16_t __U2(const unsigned char * const pa_pnSrc)
	{
		uint32_t nValue;
		memcpy(&nValue,pa_pnSrc,sizeof(nValue));
		return nValue;
	}
	static uint32_t __U4(const unsigned char * const pa_pnSrc)
	{
		uint32_t nValue;
		memcpy(&nValue,pa_pnSrc,sizeof(nValue));
		return nValue;
	}
	static int16_t __I2(const unsigned char * const pa_pnSrc)
	{
		int16_t nValue;
		memcpy(&nValue,pa_pnSrc,sizeof(nValue));
		return nValue;
	}
	static int32_t __I4(const unsigned char * const pa_pnSrc)
	{
		int32_t nValue;
		memcpy(&nValue,pa_pnSrc,sizeof(nValue));
		return nValue;
	}
	static float __R4(const unsigned char * const pa_pnSrc)
	{
		float nValue;
		memcpy(&nValue,pa_pnSrc,sizeof(nValue));
		return nValue;
	}
	static double __R8(const unsigned char * const pa_pnSrc)
	{
		double nValue;
		memcpy(&nValue,pa_pnSrc,sizeof(nValue));
		return nValue;
	}

	/* setter functions */
	static void __set_plain(unsigned char * p_pnDest, uint32_t pa_nValue, int pa_nSize)
	{
		memcpy(p_pnDest,&pa_nValue,pa_nSize);
	}
	static void __set_plainR8(unsigned char * p_pnDest, double pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}

	static void __setR8(unsigned char * p_pnDest, double pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void __setR4(unsigned char * p_pnDest, float pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void __setI4(unsigned char * p_pnDest, int32_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void __setI2(unsigned char * p_pnDest, int16_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void __setU4(unsigned char * p_pnDest, uint32_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void __setU2(unsigned char * p_pnDest, uint16_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
#endif

/* reverse byte order functions, no casts available, so code it... */
/* getters */
static uint16_t __U2r(const unsigned char * pa_pnSrc)
{
	uint16_t nValue;
	unsigned char * pnValue=(unsigned char *)&nValue+sizeof(nValue);
	int i;
	for (i=0;i<sizeof(nValue);i++) {
		*pnValue--=*pa_pnSrc++;
	}
	return nValue;
}
static uint32_t __U4r(const unsigned char * pa_pnSrc)
{
	uint32_t nValue;
	unsigned char * pnValue=(unsigned char *)&nValue+sizeof(nValue);
	int i;
	for (i=0;i<sizeof(nValue);i++) {
		*pnValue--=*pa_pnSrc++;
	}
	return nValue;
}
static int16_t __I2r(const unsigned char * pa_pnSrc)
{
	int16_t nValue;
	unsigned char * pnValue=(unsigned char *)&nValue+sizeof(nValue);
	int i;
	for (i=0;i<sizeof(nValue);i++) {
		*pnValue--=*pa_pnSrc++;
	}
	return nValue;
}
static int32_t __I4r(const unsigned char * pa_pnSrc)
{
	int32_t nValue;
	unsigned char * pnValue=(unsigned char *)&nValue+sizeof(nValue);
	int i;
	for (i=0;i<sizeof(nValue);i++) {
		*pnValue--=*pa_pnSrc++;
	}
	return nValue;
}
static float __R4r(const unsigned char * pa_pnSrc)
{
	float nValue;
	unsigned char * pnValue=(unsigned char *)&nValue+sizeof(nValue);
	int i;
	for (i=0;i<sizeof(nValue);i++) {
		*pnValue--=*pa_pnSrc++;
	}
	return nValue;
}
static double __R8r(const unsigned char * pa_pnSrc)
{
	double nValue;
	unsigned char *pnValue=(unsigned char *)&nValue+sizeof(nValue);
	int i;
	for (i=0;i<sizeof(nValue);i++) {
		*pnValue--=*pa_pnSrc++;
	}
	return nValue;
}

/* setters */
static void __set_reverse(unsigned char * pa_pnDest, uint32_t pa_nValue, int pa_nSize)
{
	unsigned char * pnValue=((unsigned char *)&pa_nValue)+pa_nSize;
	int i;
	for (i=0;i<pa_nSize;i++) {
		*pa_pnDest++=*pnValue--;
	}
}

static void __setU2r(unsigned char * pa_pnDest, uint16_t pa_nValue)
{
	unsigned char * pnValue=(unsigned char *)&pa_nValue+sizeof(pa_nValue);
	int i;
	for (i=0;i<sizeof(pa_nValue);i++) {
		*pa_pnDest++=*pnValue--;
	}
}
static void __setU4r(unsigned char * pa_pnDest, uint32_t pa_nValue)
{
	unsigned char * pnValue=(unsigned char *)&pa_nValue+sizeof(pa_nValue);
	int i;
	for (i=0;i<sizeof(pa_nValue);i++) {
		*pa_pnDest++=*pnValue--;
	}
}
static void __setI2r(unsigned char * pa_pnDest, int16_t pa_nValue)
{
	unsigned char * pnValue=(unsigned char *)&pa_nValue+sizeof(pa_nValue);
	int i;
	for (i=0;i<sizeof(pa_nValue);i++) {
		*pa_pnDest++=*pnValue--;
	}
}
static void __setI4r(unsigned char * pa_pnDest, int32_t pa_nValue)
{
	unsigned char * pnValue=(unsigned char *)&pa_nValue+sizeof(pa_nValue);
	int i;
	for (i=0;i<sizeof(pa_nValue);i++) {
		*pa_pnDest++=*pnValue--;
	}
}
static void __setR4r(unsigned char * pa_pnDest, float pa_nValue)
{
	unsigned char * pnValue=(unsigned char *)&pa_nValue+sizeof(pa_nValue);
	int i;
	for (i=0;i<sizeof(pa_nValue);i++) {
		*pa_pnDest++=*pnValue--;
	}
}
static void __setR8r(unsigned char * pa_pnDest, double pa_nValue)
{
	unsigned char * pnValue=(unsigned char *)&pa_nValue+sizeof(pa_nValue);
	int i;
	for (i=0;i<sizeof(pa_nValue);i++) {
		*pa_pnDest++=*pnValue--;
	}
}
#endif /* SERIALISATION_INLINE_H_ */
