/* --------------------------------------------------------------------------
 * serialisation.cc
 *
 *  Created on: 02.02.2013
 *      Author: thomas
 *
 *      Brief : Helper functions for serialisation
 *-------------------------------------------------------------------------*/
#include <carvig.h>

/* Getter functions -------------------------------------------------------*/
#ifdef ALIGN_STRICT
	#define U1(p)       (*((unsigned char  *)(p)))
	#define I1(p)       (*((char   *)(p)))

	static uint16_t U2(const unsigned char * const pa_pnInData)
	{
		uint32_t nValue;
		memcpy(&nValue,pa_pnInData,sizeof(nValue));
		return nValue;
	}
	static uint32_t U4(const unsigned char * const pa_pnInData)
	{
		uint32_t nValue;
		memcpy(&nValue,pa_pnInData,sizeof(nValue));
		return nValue;
	}
	static int16_t I2(const unsigned char * const pa_pnInData)
	{
		int16_t nValue;
		memcpy(&nValue,pa_pnInData,sizeof(nValue));
		return nValue;
	}
	static int32_t I4(const unsigned char * const pa_pnInData)
	{
		int32_t nValue;
		memcpy(&nValue,pa_pnInData,sizeof(nValue));
		return nValue;
	}
	static float R4(const unsigned char * const pa_pnInData)
	{
		float nValue;
		memcpy(&nValue,pa_pnInData,sizeof(nValue));
		return nValue;
	}
	static double R8(const unsigned char * p)
	{
		double value;
		unsigned char *q=(unsigned char *)&value;
		int i;
		for (i=0;i<8;i++) *q++=*p++;
		return value;
	}

	static void setR8(unsigned char * p_pnDest, double pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void setR4(unsigned char * p_pnDest, float pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void setI4(unsigned char * p_pnDest, int32_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void setI2(unsigned char * p_pnDest, int16_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void setU4(unsigned char * p_pnDest, uint32_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	static void setU2(unsigned char * p_pnDest, uint16_t pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}

#else
	#define U1(p)       (*((unsigned char  *)(p)))
	#define U2(p)       (*((unsigned short *)(p)))
	#define U4(p)       (*((unsigned int   *)(p)))
	#define I1(p)       (*((char   *)(p)))
	#define I2(p)       (*((short  *)(p)))
	#define I4(p)       (*((int    *)(p)))
	#define R4(p)       (*((float  *)(p)))

	static double R8(const unsigned char *p)
	{
		double value;
		unsigned char *q=(unsigned char *)&value;
		int i;
		for (i=0;i<8;i++) *q++=*p++;
		return value;
	}
	static void setR8(unsigned char * p_pnDest, double pa_nValue)
	{
		memcpy(p_pnDest,&pa_nValue,sizeof(pa_nValue));
	}
	#define setR4(__p,__v)		R4(__p)
	#define setI4(__p,__v)		I4(__p)
	#define setI2(__p,__v)		I2(__p)
	#define setU4(__p,__v)		U4(__p)
	#define setU2(__p,__v)		U2(__p)
#endif


