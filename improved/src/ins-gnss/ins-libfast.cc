/*-----------------------------------------------------------------------------
 * ins-libfast.cc : FAST corner detector
 *
 * reference :
 *    [1] Rosten E , Drummond T . Machine Learning for High-Speed Corner
 *        Detection[J]. 2006..
 *
 * version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
 * history : 2018/12/28 1.0 new
 *----------------------------------------------------------------------------*/
#include <carvig.h>
#include <vector>
using namespace std;

#ifndef __SSE2__
#  error "This file requires SSE2 support. Check your compiler flags."
#else
#  include <emmintrin.h>
#endif

#if __SSE2__
#define CHECK_BARRIER(lo, hi, other, flags)       \
  {                                               \
  __m128i diff = _mm_subs_epu8(lo, other);        \
  __m128i diff2 = _mm_subs_epu8(other, hi);       \
  __m128i z = _mm_setzero_si128();                \
  diff = _mm_cmpeq_epi8(diff, z);                 \
  diff2 = _mm_cmpeq_epi8(diff2, z);               \
  flags = ~(_mm_movemask_epi8(diff) | (_mm_movemask_epi8(diff2) << 16)); \
  }

template <bool Aligned> inline __m128i load_si128(const void* addr)
{
    return _mm_loadu_si128((const __m128i*)addr);
}
template <> inline __m128i load_si128<true>(const void* addr)
{
    return _mm_load_si128((const __m128i*)addr);
}
#endif

/// Check if the pointer is aligned to the specified byte granularity
template<int bytes> bool is_aligned(const void* ptr);
template<> inline bool is_aligned<8>(const void* ptr)
{
    return ((reinterpret_cast<std::size_t>(ptr)) & 0x7) == 0;
}
template<> inline bool is_aligned<16>(const void* ptr)
{
    return ((reinterpret_cast<std::size_t>(ptr)) & 0xF) == 0;
}
struct Less
{
    template <class T1, class T2> static bool eval(const T1 a, const T2 b)
    {
        return a < b;
    }
    static short prep_t(short pixel_val, short barrier)
    {
        return pixel_val - barrier;
    }
};
struct Greater
{
    template <class T1, class T2> static bool eval(const T1 a, const T2 b)
    {
        return a > b;
    }
    static short prep_t(short pixel_val, short barrier)
    {
        return pixel_val + barrier;
    }
};
struct fast_xy {
    short x,y;
    short val;
    fast_xy(short x_,short y_,short val_) :x(x_),y(y_),val(val_) {}
};
typedef unsigned char fast_byte;

static void fast_corner_detect_10(const fast_byte* img, int img_width, int img_height,
                                  int img_stride,short barrier,
                                  vector<fast_xy>& corners)
{
    int y, cb, c_b;
    const fast_byte  *line_max, *line_min;
    const fast_byte* cache_0;

    int pixel[16] = {
            0 + img_stride * 3,
            1 + img_stride * 3,
            2 + img_stride * 2,
            3 + img_stride * 1,
            3 + img_stride * 0,
            3 + img_stride * -1,
            2 + img_stride * -2,
            1 + img_stride * -3,
            0 + img_stride * -3,
            -1 + img_stride * -3,
            -2 + img_stride * -2,
            -3 + img_stride * -1,
            -3 + img_stride * 0,
            -3 + img_stride * 1,
            -2 + img_stride * 2,
            -1 + img_stride * 3,
    };
    for (y = 3 ; y < img_height - 3; y++)
    {
        cache_0 = img + y*img_stride + 3;
        line_min = cache_0 - 3;
        line_max = img + y*img_stride + img_width - 3;

        for(; cache_0 < line_max;cache_0++)
        {
            cb = *cache_0 + barrier;
            c_b= *cache_0 - barrier;

            if(*(cache_0 + pixel[0]) > cb)
                if(*(cache_0 + pixel[8]) > cb)
                    if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_0 + pixel[5]) > cb)
                            if(*(cache_0 + pixel[2]) > cb)
                                if(*(cache_0 + pixel[6]) > cb)
                                    if(*(cache_0 + 3) > cb)
                                        if(*(cache_0 + pixel[7]) > cb)
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[9]) > cb)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[9]) > cb)
                                                    if(*(cache_0 + pixel[10]) > cb)
                                                        if(*(cache_0 + pixel[11]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + pixel[10]) > cb)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[13]) > cb)
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[13]) > cb)
                                            if(*(cache_0 + pixel[14]) > cb)
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + pixel[1]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + 3) < c_b)
                                        if(*(cache_0 + pixel[10]) > cb)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + -3) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[1]) > cb)
                                                                if(*(cache_0 + pixel[15]) > cb)
                                                                    goto success;
                                                                else
                                                                if(*(cache_0 + pixel[7]) > cb)
                                                                    if(*(cache_0 + pixel[9]) > cb)
                                                                        goto success;
                                                                    else
                                                                        continue;
                                                                else
                                                                    continue;
                                                            else
                                                            if(*(cache_0 + pixel[7]) > cb)
                                                                if(*(cache_0 + pixel[9]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + -3) > cb)
                                        if(*(cache_0 + pixel[14]) > cb)
                                            if(*(cache_0 + pixel[10]) > cb)
                                                if(*(cache_0 + pixel[11]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[1]) > cb)
                                                            if(*(cache_0 + pixel[7]) > cb)
                                                                if(*(cache_0 + pixel[9]) > cb)
                                                                    goto success;
                                                                else
                                                                if(*(cache_0 + pixel[15]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else if(*(cache_0 + pixel[1]) < c_b)
                                                            if(*(cache_0 + pixel[7]) > cb)
                                                                if(*(cache_0 + pixel[9]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                        if(*(cache_0 + pixel[9]) > cb)
                                                            if(*(cache_0 + pixel[7]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[6]) < c_b)
                                    if(*(cache_0 + -3) > cb)
                                        if(*(cache_0 + pixel[13]) > cb)
                                            if(*(cache_0 + pixel[14]) > cb)
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + pixel[1]) > cb)
                                                        if(*(cache_0 + 3) > cb)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[10]) > cb)
                                                            if(*(cache_0 + pixel[11]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[7]) > cb)
                                                        if(*(cache_0 + pixel[9]) > cb)
                                                            if(*(cache_0 + pixel[10]) > cb)
                                                                if(*(cache_0 + pixel[11]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + -3) > cb)
                                    if(*(cache_0 + pixel[14]) > cb)
                                        if(*(cache_0 + pixel[15]) > cb)
                                            if(*(cache_0 + pixel[13]) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[10]) > cb)
                                                        if(*(cache_0 + pixel[11]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else if(*(cache_0 + pixel[1]) < c_b)
                                                    if(*(cache_0 + pixel[7]) > cb)
                                                        if(*(cache_0 + pixel[9]) > cb)
                                                            if(*(cache_0 + pixel[10]) > cb)
                                                                if(*(cache_0 + pixel[11]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[7]) > cb)
                                                    if(*(cache_0 + pixel[10]) > cb)
                                                        if(*(cache_0 + pixel[11]) > cb)
                                                            if(*(cache_0 + pixel[9]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + pixel[2]) < c_b)
                                if(*(cache_0 + -3) > cb)
                                    if(*(cache_0 + pixel[9]) > cb)
                                        if(*(cache_0 + pixel[10]) > cb)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + pixel[7]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        if(*(cache_0 + 3) > cb)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            if(*(cache_0 + pixel[14]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                            if(*(cache_0 + pixel[11]) > cb)
                                if(*(cache_0 + pixel[10]) > cb)
                                    if(*(cache_0 + -3) > cb)
                                        if(*(cache_0 + pixel[9]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[6]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        goto success;
                                                    else if(*(cache_0 + 3) < c_b)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            if(*(cache_0 + pixel[14]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else if(*(cache_0 + pixel[6]) < c_b)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[14]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[14]) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[5]) < c_b)
                            if(*(cache_0 + pixel[13]) > cb)
                                if(*(cache_0 + pixel[11]) > cb)
                                    if(*(cache_0 + -3) > cb)
                                        if(*(cache_0 + pixel[14]) > cb)
                                            if(*(cache_0 + pixel[15]) > cb)
                                                if(*(cache_0 + pixel[10]) > cb)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        if(*(cache_0 + pixel[1]) > cb)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[1]) > cb)
                                                        if(*(cache_0 + pixel[2]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[2]) > cb)
                                                        if(*(cache_0 + 3) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + -3) > cb)
                            if(*(cache_0 + pixel[14]) > cb)
                                if(*(cache_0 + pixel[11]) > cb)
                                    if(*(cache_0 + pixel[15]) > cb)
                                        if(*(cache_0 + pixel[10]) > cb)
                                            if(*(cache_0 + pixel[13]) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[2]) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[7]) > cb)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[10]) < c_b)
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[2]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + 3) > cb)
                                            if(*(cache_0 + pixel[2]) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + pixel[3]) < c_b)
                        if(*(cache_0 + -3) > cb)
                            if(*(cache_0 + pixel[10]) > cb)
                                if(*(cache_0 + pixel[13]) > cb)
                                    if(*(cache_0 + pixel[9]) > cb)
                                        if(*(cache_0 + pixel[11]) > cb)
                                            if(*(cache_0 + pixel[14]) > cb)
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + pixel[7]) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[1]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + -3) > cb)
                        if(*(cache_0 + pixel[10]) > cb)
                            if(*(cache_0 + pixel[14]) > cb)
                                if(*(cache_0 + pixel[11]) > cb)
                                    if(*(cache_0 + pixel[13]) > cb)
                                        if(*(cache_0 + pixel[9]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + pixel[14]) < c_b)
                                if(*(cache_0 + 3) > cb)
                                    if(*(cache_0 + pixel[5]) > cb)
                                        if(*(cache_0 + pixel[6]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[9]) > cb)
                                                    if(*(cache_0 + pixel[11]) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                            if(*(cache_0 + 3) > cb)
                                if(*(cache_0 + pixel[13]) > cb)
                                    if(*(cache_0 + pixel[6]) > cb)
                                        if(*(cache_0 + pixel[11]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else if(*(cache_0 + pixel[8]) < c_b)
                    if(*(cache_0 + pixel[11]) > cb)
                        if(*(cache_0 + pixel[2]) > cb)
                            if(*(cache_0 + pixel[15]) > cb)
                                if(*(cache_0 + pixel[1]) > cb)
                                    if(*(cache_0 + pixel[14]) > cb)
                                        if(*(cache_0 + pixel[13]) > cb)
                                            if(*(cache_0 + pixel[3]) > cb)
                                                if(*(cache_0 + -3) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[10]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + 3) > cb)
                                                    if(*(cache_0 + pixel[5]) > cb)
                                                        if(*(cache_0 + pixel[6]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[10]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[3]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[2]) < c_b)
                            if(*(cache_0 + pixel[1]) < c_b)
                                if(*(cache_0 + pixel[3]) < c_b)
                                    if(*(cache_0 + 3) < c_b)
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[6]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        if(*(cache_0 + pixel[10]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + pixel[11]) < c_b)
                        if(*(cache_0 + pixel[6]) > cb)
                            if(*(cache_0 + pixel[14]) > cb)
                                if(*(cache_0 + pixel[3]) > cb)
                                    if(*(cache_0 + pixel[1]) > cb)
                                        if(*(cache_0 + pixel[2]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[15]) > cb)
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[6]) < c_b)
                            if(*(cache_0 + pixel[10]) > cb)
                                if(*(cache_0 + pixel[1]) > cb)
                                    if(*(cache_0 + pixel[2]) > cb)
                                        if(*(cache_0 + pixel[3]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            if(*(cache_0 + pixel[14]) > cb)
                                                                if(*(cache_0 + pixel[15]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + pixel[10]) < c_b)
                                if(*(cache_0 + pixel[5]) > cb)
                                    if(*(cache_0 + pixel[7]) > cb)
                                        if(*(cache_0 + pixel[1]) > cb)
                                            if(*(cache_0 + pixel[2]) > cb)
                                                if(*(cache_0 + pixel[3]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        if(*(cache_0 + -3) > cb)
                                                            if(*(cache_0 + pixel[13]) > cb)
                                                                if(*(cache_0 + pixel[14]) > cb)
                                                                    if(*(cache_0 + pixel[15]) > cb)
                                                                        goto success;
                                                                    else
                                                                        continue;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[7]) < c_b)
                                        if(*(cache_0 + pixel[14]) > cb)
                                            if(*(cache_0 + -3) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[2]) > cb)
                                                        if(*(cache_0 + pixel[3]) > cb)
                                                            if(*(cache_0 + 3) > cb)
                                                                if(*(cache_0 + pixel[13]) > cb)
                                                                    if(*(cache_0 + pixel[15]) > cb)
                                                                        goto success;
                                                                    else
                                                                        continue;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[14]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + -3) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + -3) > cb)
                                        if(*(cache_0 + pixel[1]) > cb)
                                            if(*(cache_0 + pixel[2]) > cb)
                                                if(*(cache_0 + pixel[3]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            if(*(cache_0 + pixel[14]) > cb)
                                                                if(*(cache_0 + pixel[15]) > cb)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[5]) < c_b)
                                    if(*(cache_0 + -3) > cb)
                                        if(*(cache_0 + pixel[2]) < c_b)
                                            if(*(cache_0 + pixel[3]) < c_b)
                                                if(*(cache_0 + 3) < c_b)
                                                    if(*(cache_0 + pixel[7]) < c_b)
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + -3) < c_b)
                                        if(*(cache_0 + pixel[9]) < c_b)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        if(*(cache_0 + pixel[14]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[3]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[14]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[7]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[2]) < c_b)
                                        if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[3]) < c_b)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    if(*(cache_0 + 3) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[15]) < c_b)
                                    if(*(cache_0 + pixel[14]) < c_b)
                                        if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + -3) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                            if(*(cache_0 + -3) > cb)
                                if(*(cache_0 + pixel[1]) > cb)
                                    if(*(cache_0 + pixel[2]) > cb)
                                        if(*(cache_0 + pixel[3]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + -3) > cb)
                            if(*(cache_0 + pixel[3]) > cb)
                                if(*(cache_0 + pixel[1]) > cb)
                                    if(*(cache_0 + pixel[2]) > cb)
                                        if(*(cache_0 + 3) > cb)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[13]) > cb)
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + pixel[3]) > cb)
                        if(*(cache_0 + pixel[5]) > cb)
                            if(*(cache_0 + pixel[14]) > cb)
                                if(*(cache_0 + pixel[15]) > cb)
                                    if(*(cache_0 + pixel[13]) > cb)
                                        if(*(cache_0 + pixel[1]) > cb)
                                            if(*(cache_0 + pixel[2]) > cb)
                                                if(*(cache_0 + 3) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + -3) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[6]) > cb)
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[2]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[7]) > cb)
                                        if(*(cache_0 + pixel[1]) > cb)
                                            if(*(cache_0 + pixel[2]) > cb)
                                                if(*(cache_0 + 3) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + pixel[3]) < c_b)
                        if(*(cache_0 + pixel[1]) < c_b)
                            if(*(cache_0 + pixel[10]) < c_b)
                                if(*(cache_0 + pixel[2]) < c_b)
                                    if(*(cache_0 + 3) < c_b)
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[6]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                if(*(cache_0 + pixel[3]) > cb)
                    if(*(cache_0 + pixel[14]) > cb)
                        if(*(cache_0 + -3) > cb)
                            if(*(cache_0 + pixel[2]) > cb)
                                if(*(cache_0 + 3) > cb)
                                    if(*(cache_0 + pixel[15]) > cb)
                                        if(*(cache_0 + pixel[1]) > cb)
                                            if(*(cache_0 + pixel[13]) > cb)
                                                if(*(cache_0 + pixel[11]) > cb)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[13]) < c_b)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        if(*(cache_0 + pixel[7]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[6]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + 3) < c_b)
                                    if(*(cache_0 + pixel[1]) > cb)
                                        if(*(cache_0 + pixel[10]) > cb)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + pixel[13]) > cb)
                                                    if(*(cache_0 + pixel[15]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[10]) > cb)
                                    if(*(cache_0 + pixel[13]) > cb)
                                        if(*(cache_0 + pixel[11]) > cb)
                                            if(*(cache_0 + pixel[15]) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + -3) < c_b)
                            if(*(cache_0 + pixel[6]) > cb)
                                if(*(cache_0 + pixel[1]) > cb)
                                    if(*(cache_0 + pixel[2]) > cb)
                                        if(*(cache_0 + 3) > cb)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + pixel[7]) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + pixel[6]) > cb)
                            if(*(cache_0 + pixel[2]) > cb)
                                if(*(cache_0 + pixel[5]) > cb)
                                    if(*(cache_0 + pixel[13]) > cb)
                                        if(*(cache_0 + pixel[15]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[1]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[7]) > cb)
                                                    if(*(cache_0 + pixel[15]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[7]) > cb)
                                        if(*(cache_0 + pixel[15]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else if(*(cache_0 + pixel[3]) < c_b)
                    if(*(cache_0 + pixel[2]) > cb)
                        if(*(cache_0 + pixel[9]) > cb)
                            if(*(cache_0 + pixel[1]) > cb)
                                if(*(cache_0 + pixel[10]) > cb)
                                    if(*(cache_0 + pixel[11]) > cb)
                                        if(*(cache_0 + -3) > cb)
                                            if(*(cache_0 + pixel[13]) > cb)
                                                if(*(cache_0 + pixel[14]) > cb)
                                                    if(*(cache_0 + pixel[15]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                if(*(cache_0 + pixel[9]) > cb)
                    if(*(cache_0 + pixel[2]) > cb)
                        if(*(cache_0 + -3) > cb)
                            if(*(cache_0 + pixel[14]) > cb)
                                if(*(cache_0 + pixel[11]) > cb)
                                    if(*(cache_0 + pixel[13]) > cb)
                                        if(*(cache_0 + pixel[15]) > cb)
                                            if(*(cache_0 + pixel[10]) > cb)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                    continue;
            else if(*(cache_0 + pixel[0]) < c_b)
                if(*(cache_0 + pixel[8]) > cb)
                    if(*(cache_0 + pixel[2]) > cb)
                        if(*(cache_0 + pixel[10]) > cb)
                            if(*(cache_0 + pixel[6]) > cb)
                                if(*(cache_0 + pixel[7]) > cb)
                                    if(*(cache_0 + pixel[9]) > cb)
                                        if(*(cache_0 + pixel[5]) > cb)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + 3) > cb)
                                                    if(*(cache_0 + pixel[3]) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + -3) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + -3) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[3]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + -3) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[13]) > cb)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + -3) > cb)
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + pixel[2]) < c_b)
                        if(*(cache_0 + pixel[13]) > cb)
                            if(*(cache_0 + pixel[6]) > cb)
                                if(*(cache_0 + pixel[11]) > cb)
                                    if(*(cache_0 + pixel[9]) > cb)
                                        if(*(cache_0 + pixel[7]) > cb)
                                            if(*(cache_0 + pixel[10]) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        if(*(cache_0 + 3) > cb)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + pixel[6]) < c_b)
                                if(*(cache_0 + pixel[7]) < c_b)
                                    if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[13]) < c_b)
                            if(*(cache_0 + pixel[3]) > cb)
                                if(*(cache_0 + pixel[10]) > cb)
                                    if(*(cache_0 + pixel[7]) > cb)
                                        if(*(cache_0 + 3) > cb)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[6]) > cb)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        if(*(cache_0 + pixel[11]) > cb)
                                                            if(*(cache_0 + -3) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[10]) < c_b)
                                    if(*(cache_0 + pixel[9]) < c_b)
                                        if(*(cache_0 + pixel[1]) < c_b)
                                            if(*(cache_0 + pixel[11]) < c_b)
                                                if(*(cache_0 + -3) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + pixel[3]) < c_b)
                                if(*(cache_0 + pixel[15]) < c_b)
                                    if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[5]) > cb)
                                            if(*(cache_0 + pixel[10]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[11]) < c_b)
                                                        if(*(cache_0 + -3) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[11]) < c_b)
                                                    if(*(cache_0 + -3) < c_b)
                                                        if(*(cache_0 + pixel[14]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[6]) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + -3) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[10]) < c_b)
                                                if(*(cache_0 + pixel[11]) < c_b)
                                                    if(*(cache_0 + -3) < c_b)
                                                        if(*(cache_0 + pixel[14]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[11]) < c_b)
                                            if(*(cache_0 + pixel[10]) > cb)
                                                if(*(cache_0 + 3) < c_b)
                                                    if(*(cache_0 + -3) < c_b)
                                                        if(*(cache_0 + pixel[14]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[10]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + -3) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + -3) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                            if(*(cache_0 + pixel[9]) < c_b)
                                if(*(cache_0 + pixel[11]) < c_b)
                                    if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[10]) < c_b)
                                            if(*(cache_0 + -3) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + pixel[7]) > cb)
                            if(*(cache_0 + pixel[3]) > cb)
                                if(*(cache_0 + pixel[10]) > cb)
                                    if(*(cache_0 + 3) > cb)
                                        if(*(cache_0 + pixel[5]) > cb)
                                            if(*(cache_0 + pixel[6]) > cb)
                                                if(*(cache_0 + pixel[9]) > cb)
                                                    if(*(cache_0 + pixel[11]) > cb)
                                                        if(*(cache_0 + -3) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[7]) < c_b)
                            if(*(cache_0 + pixel[1]) < c_b)
                                if(*(cache_0 + pixel[3]) < c_b)
                                    if(*(cache_0 + 3) < c_b)
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[6]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + -3) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                            if(*(cache_0 + pixel[11]) > cb)
                                if(*(cache_0 + pixel[9]) > cb)
                                    if(*(cache_0 + pixel[10]) > cb)
                                        if(*(cache_0 + pixel[13]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + 3) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[3]) > cb)
                                            if(*(cache_0 + 3) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    if(*(cache_0 + pixel[7]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else if(*(cache_0 + pixel[8]) < c_b)
                    if(*(cache_0 + 3) > cb)
                        if(*(cache_0 + -3) < c_b)
                            if(*(cache_0 + pixel[10]) < c_b)
                                if(*(cache_0 + pixel[14]) < c_b)
                                    if(*(cache_0 + pixel[15]) < c_b)
                                        if(*(cache_0 + pixel[13]) < c_b)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[11]) < c_b)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        if(*(cache_0 + pixel[2]) < c_b)
                                                            if(*(cache_0 + pixel[3]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else if(*(cache_0 + pixel[9]) < c_b)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[3]) < c_b)
                                                        if(*(cache_0 + pixel[2]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    if(*(cache_0 + pixel[11]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[5]) < c_b)
                                        if(*(cache_0 + pixel[6]) < c_b)
                                            if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    if(*(cache_0 + pixel[11]) < c_b)
                                                        if(*(cache_0 + pixel[13]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + 3) < c_b)
                        if(*(cache_0 + pixel[2]) > cb)
                            if(*(cache_0 + pixel[10]) < c_b)
                                if(*(cache_0 + -3) < c_b)
                                    if(*(cache_0 + pixel[11]) < c_b)
                                        if(*(cache_0 + pixel[9]) < c_b)
                                            if(*(cache_0 + pixel[13]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[7]) < c_b)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            if(*(cache_0 + pixel[5]) < c_b)
                                                                if(*(cache_0 + pixel[6]) < c_b)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[6]) < c_b)
                                                            if(*(cache_0 + pixel[5]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[1]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[6]) < c_b)
                                                        if(*(cache_0 + pixel[7]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[3]) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[6]) < c_b)
                                                        if(*(cache_0 + pixel[7]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[2]) < c_b)
                            if(*(cache_0 + pixel[6]) > cb)
                                if(*(cache_0 + pixel[13]) < c_b)
                                    if(*(cache_0 + pixel[14]) < c_b)
                                        if(*(cache_0 + pixel[15]) < c_b)
                                            if(*(cache_0 + -3) < c_b)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    if(*(cache_0 + pixel[3]) < c_b)
                                                        if(*(cache_0 + pixel[11]) < c_b)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[5]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        if(*(cache_0 + pixel[10]) < c_b)
                                                            if(*(cache_0 + pixel[11]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        if(*(cache_0 + pixel[10]) < c_b)
                                                            if(*(cache_0 + pixel[11]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + pixel[6]) < c_b)
                                if(*(cache_0 + pixel[3]) > cb)
                                    if(*(cache_0 + pixel[9]) < c_b)
                                        if(*(cache_0 + pixel[10]) < c_b)
                                            if(*(cache_0 + pixel[11]) < c_b)
                                                if(*(cache_0 + -3) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        if(*(cache_0 + pixel[7]) < c_b)
                                                            if(*(cache_0 + pixel[5]) < c_b)
                                                                goto success;
                                                            else
                                                            if(*(cache_0 + pixel[14]) < c_b)
                                                                if(*(cache_0 + pixel[15]) < c_b)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                        if(*(cache_0 + pixel[1]) < c_b)
                                                            if(*(cache_0 + pixel[14]) < c_b)
                                                                if(*(cache_0 + pixel[15]) < c_b)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[3]) < c_b)
                                    if(*(cache_0 + pixel[5]) > cb)
                                        if(*(cache_0 + pixel[11]) < c_b)
                                            if(*(cache_0 + -3) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            if(*(cache_0 + pixel[1]) < c_b)
                                                                goto success;
                                                            else
                                                            if(*(cache_0 + pixel[7]) < c_b)
                                                                if(*(cache_0 + pixel[9]) < c_b)
                                                                    if(*(cache_0 + pixel[10]) < c_b)
                                                                        goto success;
                                                                    else
                                                                        continue;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[5]) < c_b)
                                        if(*(cache_0 + pixel[7]) > cb)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[1]) > cb)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    if(*(cache_0 + pixel[10]) < c_b)
                                                        if(*(cache_0 + pixel[11]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[15]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[11]) < c_b)
                                                if(*(cache_0 + pixel[10]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[13]) < c_b)
                                            if(*(cache_0 + pixel[15]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[1]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + -3) < c_b)
                                        if(*(cache_0 + pixel[14]) < c_b)
                                            if(*(cache_0 + pixel[11]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        if(*(cache_0 + pixel[1]) > cb)
                                                            if(*(cache_0 + pixel[7]) < c_b)
                                                                if(*(cache_0 + pixel[9]) < c_b)
                                                                    if(*(cache_0 + pixel[10]) < c_b)
                                                                        goto success;
                                                                    else
                                                                        continue;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else if(*(cache_0 + pixel[1]) < c_b)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            if(*(cache_0 + pixel[7]) < c_b)
                                                                if(*(cache_0 + pixel[10]) < c_b)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[11]) < c_b)
                                    if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[10]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + -3) < c_b)
                                                    if(*(cache_0 + pixel[7]) > cb)
                                                        if(*(cache_0 + pixel[1]) < c_b)
                                                            if(*(cache_0 + pixel[14]) < c_b)
                                                                if(*(cache_0 + pixel[15]) < c_b)
                                                                    goto success;
                                                                else
                                                                    continue;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else if(*(cache_0 + pixel[7]) < c_b)
                                                        if(*(cache_0 + pixel[5]) < c_b)
                                                            goto success;
                                                        else
                                                        if(*(cache_0 + pixel[14]) < c_b)
                                                            if(*(cache_0 + pixel[15]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        if(*(cache_0 + pixel[1]) < c_b)
                                                            if(*(cache_0 + pixel[14]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                            if(*(cache_0 + -3) < c_b)
                                if(*(cache_0 + pixel[14]) < c_b)
                                    if(*(cache_0 + pixel[15]) < c_b)
                                        if(*(cache_0 + pixel[13]) < c_b)
                                            if(*(cache_0 + pixel[11]) > cb)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    if(*(cache_0 + pixel[3]) < c_b)
                                                        if(*(cache_0 + pixel[5]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[11]) < c_b)
                                                if(*(cache_0 + pixel[1]) > cb)
                                                    if(*(cache_0 + pixel[7]) < c_b)
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            if(*(cache_0 + pixel[10]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else if(*(cache_0 + pixel[1]) < c_b)
                                                    if(*(cache_0 + pixel[3]) > cb)
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            if(*(cache_0 + pixel[10]) < c_b)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else if(*(cache_0 + pixel[3]) < c_b)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[10]) < c_b)
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[10]) < c_b)
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    if(*(cache_0 + pixel[1]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + pixel[11]) < c_b)
                            if(*(cache_0 + pixel[10]) < c_b)
                                if(*(cache_0 + -3) < c_b)
                                    if(*(cache_0 + pixel[9]) < c_b)
                                        if(*(cache_0 + pixel[13]) > cb)
                                            if(*(cache_0 + pixel[3]) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[6]) < c_b)
                                                        if(*(cache_0 + pixel[7]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else if(*(cache_0 + pixel[13]) < c_b)
                                            if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[6]) < c_b)
                                                    if(*(cache_0 + pixel[5]) < c_b)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + pixel[6]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[5]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + -3) < c_b)
                        if(*(cache_0 + pixel[10]) < c_b)
                            if(*(cache_0 + pixel[14]) < c_b)
                                if(*(cache_0 + pixel[11]) < c_b)
                                    if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[15]) < c_b)
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    if(*(cache_0 + pixel[2]) < c_b)
                                                        if(*(cache_0 + pixel[3]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[3]) < c_b)
                                                if(*(cache_0 + pixel[2]) < c_b)
                                                    if(*(cache_0 + pixel[1]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[6]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                if(*(cache_0 + pixel[2]) < c_b)
                    if(*(cache_0 + -3) > cb)
                        if(*(cache_0 + pixel[6]) < c_b)
                            if(*(cache_0 + pixel[14]) < c_b)
                                if(*(cache_0 + pixel[7]) > cb)
                                    if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[7]) < c_b)
                                    if(*(cache_0 + 3) < c_b)
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[13]) < c_b)
                                    if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + 3) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + -3) < c_b)
                        if(*(cache_0 + pixel[3]) > cb)
                            if(*(cache_0 + pixel[9]) < c_b)
                                if(*(cache_0 + pixel[11]) < c_b)
                                    if(*(cache_0 + pixel[14]) < c_b)
                                        if(*(cache_0 + pixel[13]) < c_b)
                                            if(*(cache_0 + pixel[15]) < c_b)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    if(*(cache_0 + pixel[10]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[3]) < c_b)
                            if(*(cache_0 + pixel[14]) < c_b)
                                if(*(cache_0 + 3) > cb)
                                    if(*(cache_0 + pixel[10]) < c_b)
                                        if(*(cache_0 + pixel[15]) < c_b)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[11]) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + 3) < c_b)
                                    if(*(cache_0 + pixel[15]) < c_b)
                                        if(*(cache_0 + pixel[1]) < c_b)
                                            if(*(cache_0 + pixel[13]) > cb)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[6]) < c_b)
                                                        if(*(cache_0 + pixel[7]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[13]) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[11]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[6]) < c_b)
                                                    if(*(cache_0 + pixel[5]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[10]) < c_b)
                                    if(*(cache_0 + pixel[11]) < c_b)
                                        if(*(cache_0 + pixel[15]) < c_b)
                                            if(*(cache_0 + pixel[13]) < c_b)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + pixel[9]) < c_b)
                            if(*(cache_0 + pixel[10]) < c_b)
                                if(*(cache_0 + pixel[14]) < c_b)
                                    if(*(cache_0 + pixel[11]) < c_b)
                                        if(*(cache_0 + pixel[15]) < c_b)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + pixel[6]) < c_b)
                        if(*(cache_0 + pixel[14]) < c_b)
                            if(*(cache_0 + 3) < c_b)
                                if(*(cache_0 + pixel[13]) > cb)
                                    if(*(cache_0 + pixel[7]) < c_b)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[13]) < c_b)
                                    if(*(cache_0 + pixel[5]) < c_b)
                                        if(*(cache_0 + pixel[15]) < c_b)
                                            if(*(cache_0 + pixel[1]) < c_b)
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[7]) < c_b)
                                    if(*(cache_0 + pixel[15]) < c_b)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + pixel[1]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                    continue;
            else
            if(*(cache_0 + pixel[8]) > cb)
                if(*(cache_0 + pixel[10]) > cb)
                    if(*(cache_0 + 3) > cb)
                        if(*(cache_0 + pixel[2]) > cb)
                            if(*(cache_0 + pixel[6]) > cb)
                                if(*(cache_0 + pixel[7]) > cb)
                                    if(*(cache_0 + pixel[11]) > cb)
                                        if(*(cache_0 + pixel[9]) > cb)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[3]) > cb)
                                                    goto success;
                                                else if(*(cache_0 + pixel[3]) < c_b)
                                                    if(*(cache_0 + -3) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[13]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + -3) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        if(*(cache_0 + pixel[14]) > cb)
                                                            if(*(cache_0 + pixel[15]) > cb)
                                                                goto success;
                                                            else
                                                                continue;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[15]) > cb)
                                                if(*(cache_0 + pixel[14]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        if(*(cache_0 + pixel[13]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[1]) > cb)
                                        if(*(cache_0 + pixel[3]) > cb)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[9]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else if(*(cache_0 + pixel[2]) < c_b)
                            if(*(cache_0 + pixel[11]) > cb)
                                if(*(cache_0 + -3) > cb)
                                    if(*(cache_0 + pixel[9]) > cb)
                                        if(*(cache_0 + pixel[6]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[13]) > cb)
                                                    if(*(cache_0 + pixel[5]) > cb)
                                                        goto success;
                                                    else
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[3]) > cb)
                                                    if(*(cache_0 + pixel[5]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                        if(*(cache_0 + -3) > cb)
                            if(*(cache_0 + pixel[6]) > cb)
                                if(*(cache_0 + pixel[11]) > cb)
                                    if(*(cache_0 + pixel[13]) > cb)
                                        if(*(cache_0 + pixel[7]) > cb)
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    goto success;
                                                else if(*(cache_0 + pixel[5]) < c_b)
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        if(*(cache_0 + pixel[15]) > cb)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                if(*(cache_0 + pixel[15]) > cb)
                                                    if(*(cache_0 + pixel[14]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[3]) > cb)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[7]) > cb)
                                                    if(*(cache_0 + pixel[9]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[3]) > cb)
                                        if(*(cache_0 + pixel[7]) > cb)
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[5]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + 3) < c_b)
                        if(*(cache_0 + pixel[6]) > cb)
                            if(*(cache_0 + pixel[14]) > cb)
                                if(*(cache_0 + pixel[13]) > cb)
                                    if(*(cache_0 + pixel[7]) > cb)
                                        if(*(cache_0 + pixel[15]) > cb)
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[11]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                        if(*(cache_0 + pixel[5]) > cb)
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[11]) > cb)
                                                    if(*(cache_0 + -3) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + pixel[14]) > cb)
                        if(*(cache_0 + pixel[6]) > cb)
                            if(*(cache_0 + -3) > cb)
                                if(*(cache_0 + pixel[5]) > cb)
                                    if(*(cache_0 + pixel[11]) > cb)
                                        if(*(cache_0 + pixel[9]) > cb)
                                            if(*(cache_0 + pixel[7]) > cb)
                                                if(*(cache_0 + pixel[13]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[5]) < c_b)
                                    if(*(cache_0 + pixel[15]) > cb)
                                        if(*(cache_0 + pixel[7]) > cb)
                                            if(*(cache_0 + pixel[9]) > cb)
                                                if(*(cache_0 + pixel[11]) > cb)
                                                    if(*(cache_0 + pixel[13]) > cb)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[15]) > cb)
                                    if(*(cache_0 + pixel[11]) > cb)
                                        if(*(cache_0 + pixel[9]) > cb)
                                            if(*(cache_0 + pixel[13]) > cb)
                                                if(*(cache_0 + pixel[7]) > cb)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                    continue;
            else if(*(cache_0 + pixel[8]) < c_b)
                if(*(cache_0 + pixel[10]) < c_b)
                    if(*(cache_0 + 3) > cb)
                        if(*(cache_0 + pixel[14]) < c_b)
                            if(*(cache_0 + pixel[6]) < c_b)
                                if(*(cache_0 + -3) < c_b)
                                    if(*(cache_0 + pixel[9]) < c_b)
                                        if(*(cache_0 + pixel[11]) < c_b)
                                            if(*(cache_0 + pixel[15]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[7]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else if(*(cache_0 + 3) < c_b)
                        if(*(cache_0 + pixel[6]) < c_b)
                            if(*(cache_0 + -3) > cb)
                                if(*(cache_0 + pixel[2]) < c_b)
                                    if(*(cache_0 + pixel[1]) > cb)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        if(*(cache_0 + pixel[11]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    if(*(cache_0 + pixel[7]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[11]) < c_b)
                                        if(*(cache_0 + pixel[3]) < c_b)
                                            if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + pixel[7]) < c_b)
                                                    if(*(cache_0 + pixel[9]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else if(*(cache_0 + -3) < c_b)
                                if(*(cache_0 + pixel[7]) < c_b)
                                    if(*(cache_0 + pixel[11]) > cb)
                                        if(*(cache_0 + pixel[1]) < c_b)
                                            if(*(cache_0 + pixel[2]) < c_b)
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    if(*(cache_0 + pixel[5]) < c_b)
                                                        if(*(cache_0 + pixel[9]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else if(*(cache_0 + pixel[11]) < c_b)
                                        if(*(cache_0 + pixel[9]) < c_b)
                                            if(*(cache_0 + pixel[5]) > cb)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[14]) < c_b)
                                                        if(*(cache_0 + pixel[15]) < c_b)
                                                            goto success;
                                                        else
                                                            continue;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else if(*(cache_0 + pixel[5]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    goto success;
                                                else
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                            if(*(cache_0 + pixel[15]) < c_b)
                                                if(*(cache_0 + pixel[14]) < c_b)
                                                    if(*(cache_0 + pixel[13]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                    if(*(cache_0 + pixel[1]) < c_b)
                                        if(*(cache_0 + pixel[2]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + pixel[3]) < c_b)
                                                    if(*(cache_0 + pixel[5]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                            if(*(cache_0 + pixel[2]) < c_b)
                                if(*(cache_0 + pixel[1]) < c_b)
                                    if(*(cache_0 + pixel[3]) < c_b)
                                        if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + pixel[5]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[11]) < c_b)
                                    if(*(cache_0 + pixel[3]) < c_b)
                                        if(*(cache_0 + pixel[5]) < c_b)
                                            if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                    if(*(cache_0 + pixel[14]) < c_b)
                        if(*(cache_0 + pixel[6]) < c_b)
                            if(*(cache_0 + -3) < c_b)
                                if(*(cache_0 + pixel[5]) > cb)
                                    if(*(cache_0 + pixel[9]) < c_b)
                                        if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[11]) < c_b)
                                                if(*(cache_0 + pixel[13]) < c_b)
                                                    if(*(cache_0 + pixel[15]) < c_b)
                                                        goto success;
                                                    else
                                                        continue;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else if(*(cache_0 + pixel[5]) < c_b)
                                    if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[11]) < c_b)
                                            if(*(cache_0 + pixel[7]) < c_b)
                                                if(*(cache_0 + pixel[9]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                if(*(cache_0 + pixel[15]) < c_b)
                                    if(*(cache_0 + pixel[13]) < c_b)
                                        if(*(cache_0 + pixel[7]) < c_b)
                                            if(*(cache_0 + pixel[9]) < c_b)
                                                if(*(cache_0 + pixel[11]) < c_b)
                                                    goto success;
                                                else
                                                    continue;
                                            else
                                                continue;
                                        else
                                            continue;
                                    else
                                        continue;
                                else
                                    continue;
                            else
                                continue;
                        else
                            continue;
                    else
                        continue;
                else
                    continue;
            else
                continue;
success:
            corners.push_back(fast_xy(static_cast<short>(cache_0-line_min), /* x. */
                                      static_cast<short>(y),                /* y. */
                                      img[y*img_width+cache_0-line_min]));  /* val. */
        }
    }
}
static inline bool test_gt_set(int a, int b, int& min_diff)
{
    if(a > b)
    {
        if(a-b < min_diff)
            min_diff = a-b;

        return 1;
    }
    return 0;
}
static inline int fast_corner_score_10(const fast_byte* cache_0,
                                       const int offset[], int b)
{
    b++;
    //This function computes the score for a pixel which is known to be
    //a corner at barrier b. So we start looking at b+1 and above to
    //establish where it stops becoming a corner.
    for(;;)
    {
        int cb = *cache_0 + b;
        int c_b= *cache_0 - b;
        int min_diff = INT_MAX;
        if(test_gt_set(*(cache_0 + offset[0]), cb, min_diff))
            if(test_gt_set(*(cache_0 + offset[8]), cb, min_diff))
                if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                                    b += min_diff;
                                                                else
                                                                    break;
                                                            else
                                                                break;
                                                        else
                                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    b += min_diff;
                                                else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else if(test_gt_set(c_b, *(cache_0 + offset[8]), min_diff))
                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                                    b += min_diff;
                                                                else
                                                                    break;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                                    b += min_diff;
                                                                else
                                                                    break;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                        if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
                break;
        else if(test_gt_set(c_b, *(cache_0 + offset[0]), min_diff))
            if(test_gt_set(*(cache_0 + offset[8]), cb, min_diff))
                if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else if(test_gt_set(c_b, *(cache_0 + offset[8]), min_diff))
                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                    if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                            b += min_diff;
                                                        else
                                                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                            b += min_diff;
                                                        else
                                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                                    b += min_diff;
                                                                else
                                                                    break;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                                    b += min_diff;
                                                                else
                                                                    break;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                                b += min_diff;
                                                            else
                                                                break;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                        b += min_diff;
                                                    else
                                                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                        if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                            if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
            if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
                break;
        else
        if(test_gt_set(*(cache_0 + offset[8]), cb, min_diff))
            if(test_gt_set(*(cache_0 + offset[10]), cb, min_diff))
                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[2]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                b += min_diff;
                                            else if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                            b += min_diff;
                                                        else
                                                            break;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                    if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                b += min_diff;
                                            else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(*(cache_0 + offset[3]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                    if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(*(cache_0 + offset[14]), cb, min_diff))
                    if(test_gt_set(*(cache_0 + offset[6]), cb, min_diff))
                        if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                                if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(*(cache_0 + offset[15]), cb, min_diff))
                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                    if(test_gt_set(*(cache_0 + offset[9]), cb, min_diff))
                                        if(test_gt_set(*(cache_0 + offset[13]), cb, min_diff))
                                            if(test_gt_set(*(cache_0 + offset[7]), cb, min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
                break;
        else if(test_gt_set(c_b, *(cache_0 + offset[8]), min_diff))
            if(test_gt_set(c_b, *(cache_0 + offset[10]), min_diff))
                if(test_gt_set(*(cache_0 + offset[4]), cb, min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else if(test_gt_set(c_b, *(cache_0 + offset[4]), min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                        if(test_gt_set(*(cache_0 + offset[12]), cb, min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[1]), cb, min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                if(test_gt_set(*(cache_0 + offset[11]), cb, min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                        if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                    if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                        b += min_diff;
                                                    else
                                                        break;
                                                else
                                                    break;
                                            else
                                                break;
                                        else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                b += min_diff;
                                            else
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                        if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                        if(test_gt_set(c_b, *(cache_0 + offset[2]), min_diff))
                            if(test_gt_set(c_b, *(cache_0 + offset[1]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[3]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                if(test_gt_set(c_b, *(cache_0 + offset[14]), min_diff))
                    if(test_gt_set(c_b, *(cache_0 + offset[6]), min_diff))
                        if(test_gt_set(c_b, *(cache_0 + offset[12]), min_diff))
                            if(test_gt_set(*(cache_0 + offset[5]), cb, min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                                if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                                    b += min_diff;
                                                else
                                                    break;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else if(test_gt_set(c_b, *(cache_0 + offset[5]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                            if(test_gt_set(c_b, *(cache_0 + offset[15]), min_diff))
                                if(test_gt_set(c_b, *(cache_0 + offset[13]), min_diff))
                                    if(test_gt_set(c_b, *(cache_0 + offset[7]), min_diff))
                                        if(test_gt_set(c_b, *(cache_0 + offset[9]), min_diff))
                                            if(test_gt_set(c_b, *(cache_0 + offset[11]), min_diff))
                                                b += min_diff;
                                            else
                                                break;
                                        else
                                            break;
                                    else
                                        break;
                                else
                                    break;
                            else
                                break;
                        else
                            break;
                    else
                        break;
                else
                    break;
            else
                break;
        else
            break;
    }
    return b-1;
}
static void fast_corner_score_10(const fast_byte* img,const int img_stride,
                                 const vector<fast_xy>& corners,
                                 const int threshold,
                                 vector<int>& scores)
{
    scores.resize(corners.size());
    int pixel[16] = {
            0 + img_stride * 3,
            1 + img_stride * 3,
            2 + img_stride * 2,
            3 + img_stride * 1,
            3 + img_stride * 0,
            3 + img_stride * -1,
            2 + img_stride * -2,
            1 + img_stride * -3,
            0 + img_stride * -3,
            -1 + img_stride * -3,
            -2 + img_stride * -2,
            -3 + img_stride * -1,
            -3 + img_stride * 0,
            -3 + img_stride * 1,
            -2 + img_stride * 2,
            -1 + img_stride * 3,
    };
    for (unsigned int n=0; n < corners.size(); n++)
        scores[n] = fast_corner_score_10(img + corners[n].y*img_stride + corners[n].x, pixel, threshold);
}
static void fast_nonmax_3x3(const vector<fast_xy>& corners,const vector<int>& scores,
                            vector<int>& nonmax_corners)
{
    nonmax_corners.clear();
    nonmax_corners.reserve(corners.size());
    if (corners.size() < 1) return;

    // Find where each row begins
    // (the corners are output in raster scan order). A beginning of -1 signifies
    // that there are no corners on that row.
    int last_row = corners.back().y;
    vector<int> row_start(last_row + 1, -1);

    int prev_row = -1;
    for (unsigned int i=0; i< corners.size(); i++)
        if (corners[i].y != prev_row)
        {
            row_start[corners[i].y] = i;
            prev_row = corners[i].y;
        }

    //Point above points (roughly) to the pixel above the one of interest, if there
    //is a feature there.
    int point_above = 0;
    int point_below = 0;

    const int sz = (int)corners.size();

    for(int i=0; i < sz; i++)
    {
        int score = scores[i];
        fast_xy pos = corners[i];

        //Check left
        if (i > 0)
            //if(corners[i-1] == pos-ImageRef(1,0) && (scores[i-1] >= score))
            if (corners[i-1].x == pos.x-1 && corners[i-1].y == pos.y && scores[i-1] >= score)
                continue;

        //Check right
        if (i < (sz - 1))
            //if(corners[i+1] == pos+ImageRef(1,0) &&  (scores[i+1] >= score))
            if (corners[i+1].x == pos.x+1 && corners[i+1].y == pos.y && scores[i+1] >= score)
                continue;

        //Check above (if there is a valid row above)
        if(pos.y != 0 && row_start[pos.y - 1] != -1)
        {
            //Make sure that current point_above is one
            //row above.
            if(corners[point_above].y < pos.y - 1)
                point_above = row_start[pos.y-1];

            //Make point_above point to the first of the pixels above the current point,
            //if it exists.
            for(; corners[point_above].y < pos.y && corners[point_above].x < pos.x - 1; point_above++)
            {}


            for(int i=point_above; corners[i].y < pos.y && corners[i].x <= pos.x + 1; i++)
            {
                int x = corners[i].x;
                if( (x == pos.x - 1 || x ==pos.x || x == pos.x+1) && (scores[i] >= score))
                    goto cont;
            }

        }
        //Check below (if there is anything below)
        if(pos.y != last_row && row_start[pos.y + 1] != -1 && point_below < sz) //Nothing below
        {
            if (corners[point_below].y < pos.y + 1)
                point_below = row_start[pos.y+1];

            // Make point below point to one of the pixels belowthe current point, if it
            // exists.
            for(; point_below < sz && corners[point_below].y == pos.y+1 && corners[point_below].x < pos.x - 1; point_below++)
            {}
            for(int i=point_below; i < sz && corners[i].y == pos.y+1 && corners[i].x <= pos.x + 1; i++)
            {
                int x = corners[i].x;
                if( (x == pos.x - 1 || x ==pos.x || x == pos.x+1) && (scores[i] >= score))
                    goto cont;
            }
        }
        nonmax_corners.push_back(i);
cont:
        ;
    }
}
template <class C> inline bool is_corner_10(const unsigned char* p, const int w,
                                            const int barrier) {
    const int w3 = 3*w;
    const int t = C::prep_t(*p, barrier);
    if (!C::eval(p[-1-w3],t)) { // ???????????????-
        if (!C::eval(p[3+w],t)) { // ?????-?????????-
            return false;
        } // ?????@?????????-
        if (!C::eval(p[2+2*w],t)) { // ?????@-????????-
            return false;
        } // ?????@@????????-
        if (!C::eval(p[-1+w3],t)) { // ?????@@??-?????-
            return false;
        } // ?????@@??@?????-
        if (!C::eval(p[1+w3],t)) { // ?????@@-?@?????-
            return false;
        } // ?????@@@?@?????-
        if (!C::eval(p[w3],t)) { // ?????@@@-@?????-
            return false;
        } // ?????@@@@@?????-
        if (!C::eval(p[-2+2*w],t)) { // ?????@@@@@-????-
            if (!C::eval(p[-w3],t)) { // -????@@@@@-????-
                return false;
            } // @????@@@@@-????-
            if (!C::eval(p[3],t)) { // @???-@@@@@-????-
                return false;
            } // @???@@@@@@-????-
            if (!C::eval(p[1-w3],t)) { // @-??@@@@@@-????-
                return false;
            } // @@??@@@@@@-????-
            if (!C::eval(p[2-2*w],t)) { // @@-?@@@@@@-????-
                return false;
            } // @@@?@@@@@@-????-
            if (!C::eval(p[3-w],t)) { // @@@-@@@@@@-????-
                return false;
            } // @@@@@@@@@@-????-
            return true;
        } // ?????@@@@@@????-
        if (!C::eval(p[-3+w],t)) { // ?????@@@@@@-???-
            if (!C::eval(p[3],t)) { // ????-@@@@@@-???-
                return false;
            } // ????@@@@@@@-???-
            if (!C::eval(p[1-w3],t)) { // ?-??@@@@@@@-???-
                return false;
            } // ?@??@@@@@@@-???-
            if (!C::eval(p[2-2*w],t)) { // ?@-?@@@@@@@-???-
                return false;
            } // ?@@?@@@@@@@-???-
            if (!C::eval(p[3-w],t)) { // ?@@-@@@@@@@-???-
                return false;
            } // ?@@@@@@@@@@-???-
            return true;
        } // ?????@@@@@@@???-
        if (!C::eval(p[3],t)) { // ????-@@@@@@@???-
            if (!C::eval(p[-3],t)) { // ????-@@@@@@@-??-
                return false;
            } // ????-@@@@@@@@??-
            if (!C::eval(p[-3-w],t)) { // ????-@@@@@@@@-?-
                return false;
            } // ????-@@@@@@@@@?-
            if (!C::eval(p[-2-2*w],t)) { // ????-@@@@@@@@@--
                return false;
            } // ????-@@@@@@@@@@-
            return true;
        } // ????@@@@@@@@???-
        if (!C::eval(p[3-w],t)) { // ???-@@@@@@@@???-
            if (!C::eval(p[-3],t)) { // ???-@@@@@@@@-??-
                return false;
            } // ???-@@@@@@@@@??-
            if (!C::eval(p[-3-w],t)) { // ???-@@@@@@@@@-?-
                return false;
            } // ???-@@@@@@@@@@?-
            return true;
        } // ???@@@@@@@@@???-
        if (!C::eval(p[-3],t)) { // ???@@@@@@@@@-??-
            if (!C::eval(p[2-2*w],t)) { // ??-@@@@@@@@@-??-
                return false;
            } // ??@@@@@@@@@@-??-
            return true;
        } // ???@@@@@@@@@@??-
        return true;
    } // ???????????????@
    if (!C::eval(p[-2-2*w],t)) { // ??????????????-@
        if (!C::eval(p[3],t)) { // ????-?????????-@
            return false;
        } // ????@?????????-@
        if (!C::eval(p[3+w],t)) { // ????@-????????-@
            return false;
        } // ????@@????????-@
        if (!C::eval(p[w3],t)) { // ????@@??-?????-@
            return false;
        } // ????@@??@?????-@
        if (!C::eval(p[1+w3],t)) { // ????@@?-@?????-@
            return false;
        } // ????@@?@@?????-@
        if (!C::eval(p[2+2*w],t)) { // ????@@-@@?????-@
            return false;
        } // ????@@@@@?????-@
        if (!C::eval(p[3-w],t)) { // ???-@@@@@?????-@
            if (!C::eval(p[-1+w3],t)) { // ???-@@@@@-????-@
                return false;
            } // ???-@@@@@@????-@
            if (!C::eval(p[-3-w],t)) { // ???-@@@@@@???--@
                return false;
            } // ???-@@@@@@???@-@
            if (!C::eval(p[-2+2*w],t)) { // ???-@@@@@@-??@-@
                return false;
            } // ???-@@@@@@@??@-@
            if (!C::eval(p[-3+w],t)) { // ???-@@@@@@@-?@-@
                return false;
            } // ???-@@@@@@@@?@-@
            if (!C::eval(p[-3],t)) { // ???-@@@@@@@@-@-@
                return false;
            } // ???-@@@@@@@@@@-@
            return true;
        } // ???@@@@@@?????-@
        if (!C::eval(p[2-2*w],t)) { // ??-@@@@@@?????-@
            if (!C::eval(p[-3],t)) { // ??-@@@@@@???-?-@
                return false;
            } // ??-@@@@@@???@?-@
            if (!C::eval(p[-1+w3],t)) { // ??-@@@@@@-??@?-@
                return false;
            } // ??-@@@@@@@??@?-@
            if (!C::eval(p[-2+2*w],t)) { // ??-@@@@@@@-?@?-@
                return false;
            } // ??-@@@@@@@@?@?-@
            if (!C::eval(p[-3+w],t)) { // ??-@@@@@@@@-@?-@
                return false;
            } // ??-@@@@@@@@@@?-@
            return true;
        } // ??@@@@@@@?????-@
        if (!C::eval(p[1-w3],t)) { // ?-@@@@@@@?????-@
            if (!C::eval(p[-1+w3],t)) { // ?-@@@@@@@-????-@
                return false;
            } // ?-@@@@@@@@????-@
            if (!C::eval(p[-2+2*w],t)) { // ?-@@@@@@@@-???-@
                return false;
            } // ?-@@@@@@@@@???-@
            if (!C::eval(p[-3+w],t)) { // ?-@@@@@@@@@-??-@
                return false;
            } // ?-@@@@@@@@@@??-@
            return true;
        } // ?@@@@@@@@?????-@
        if (!C::eval(p[-w3],t)) { // -@@@@@@@@?????-@
            if (!C::eval(p[-1+w3],t)) { // -@@@@@@@@-????-@
                return false;
            } // -@@@@@@@@@????-@
            if (!C::eval(p[-2+2*w],t)) { // -@@@@@@@@@-???-@
                return false;
            } // -@@@@@@@@@@???-@
            return true;
        } // @@@@@@@@@?????-@
        return true;
    } // ??????????????@@
    if (!C::eval(p[-3-w],t)) { // ?????????????-@@
        if (!C::eval(p[1+w3],t)) { // ???????-?????-@@
            return false;
        } // ???????@?????-@@
        if (!C::eval(p[3-w],t)) { // ???-???@?????-@@
            return false;
        } // ???@???@?????-@@
        if (!C::eval(p[3],t)) { // ???@-??@?????-@@
            return false;
        } // ???@@??@?????-@@
        if (!C::eval(p[3+w],t)) { // ???@@-?@?????-@@
            return false;
        } // ???@@@?@?????-@@
        if (!C::eval(p[2+2*w],t)) { // ???@@@-@?????-@@
            return false;
        } // ???@@@@@?????-@@
        if (!C::eval(p[2-2*w],t)) { // ??-@@@@@?????-@@
            if (!C::eval(p[w3],t)) { // ??-@@@@@-????-@@
                return false;
            } // ??-@@@@@@????-@@
            if (!C::eval(p[-3],t)) { // ??-@@@@@@???--@@
                return false;
            } // ??-@@@@@@???@-@@
            if (!C::eval(p[-1+w3],t)) { // ??-@@@@@@-??@-@@
                return false;
            } // ??-@@@@@@@??@-@@
            if (!C::eval(p[-2+2*w],t)) { // ??-@@@@@@@-?@-@@
                return false;
            } // ??-@@@@@@@@?@-@@
            if (!C::eval(p[-3+w],t)) { // ??-@@@@@@@@-@-@@
                return false;
            } // ??-@@@@@@@@@@-@@
            return true;
        } // ??@@@@@@?????-@@
        if (!C::eval(p[1-w3],t)) { // ?-@@@@@@?????-@@
            if (!C::eval(p[-3+w],t)) { // ?-@@@@@@???-?-@@
                return false;
            } // ?-@@@@@@???@?-@@
            if (!C::eval(p[w3],t)) { // ?-@@@@@@-??@?-@@
                return false;
            } // ?-@@@@@@@??@?-@@
            if (!C::eval(p[-1+w3],t)) { // ?-@@@@@@@-?@?-@@
                return false;
            } // ?-@@@@@@@@?@?-@@
            if (!C::eval(p[-2+2*w],t)) { // ?-@@@@@@@@-@?-@@
                return false;
            } // ?-@@@@@@@@@@?-@@
            return true;
        } // ?@@@@@@@?????-@@
        if (!C::eval(p[-w3],t)) { // -@@@@@@@?????-@@
            if (!C::eval(p[w3],t)) { // -@@@@@@@-????-@@
                return false;
            } // -@@@@@@@@????-@@
            if (!C::eval(p[-1+w3],t)) { // -@@@@@@@@-???-@@
                return false;
            } // -@@@@@@@@@???-@@
            if (!C::eval(p[-2+2*w],t)) { // -@@@@@@@@@-??-@@
                return false;
            } // -@@@@@@@@@@??-@@
            return true;
        } // @@@@@@@@?????-@@
        return true;
    } // ?????????????@@@
    if (!C::eval(p[-w3],t)) { // -????????????@@@
        if (!C::eval(p[2+2*w],t)) { // -?????-??????@@@
            return false;
        } // -?????@??????@@@
        if (!C::eval(p[1+w3],t)) { // -?????@-?????@@@
            return false;
        } // -?????@@?????@@@
        if (!C::eval(p[-2+2*w],t)) { // -?????@@??-??@@@
            return false;
        } // -?????@@??@??@@@
        if (!C::eval(p[w3],t)) { // -?????@@-?@??@@@
            return false;
        } // -?????@@@?@??@@@
        if (!C::eval(p[-1+w3],t)) { // -?????@@@-@??@@@
            return false;
        } // -?????@@@@@??@@@
        if (!C::eval(p[-3+w],t)) { // -?????@@@@@-?@@@
            if (!C::eval(p[1-w3],t)) { // --????@@@@@-?@@@
                return false;
            } // -@????@@@@@-?@@@
            if (!C::eval(p[3+w],t)) { // -@???-@@@@@-?@@@
                return false;
            } // -@???@@@@@@-?@@@
            if (!C::eval(p[2-2*w],t)) { // -@-??@@@@@@-?@@@
                return false;
            } // -@@??@@@@@@-?@@@
            if (!C::eval(p[3-w],t)) { // -@@-?@@@@@@-?@@@
                return false;
            } // -@@@?@@@@@@-?@@@
            if (!C::eval(p[3],t)) { // -@@@-@@@@@@-?@@@
                return false;
            } // -@@@@@@@@@@-?@@@
            return true;
        } // -?????@@@@@@?@@@
        if (!C::eval(p[-3],t)) { // -?????@@@@@@-@@@
            if (!C::eval(p[3+w],t)) { // -????-@@@@@@-@@@
                return false;
            } // -????@@@@@@@-@@@
            if (!C::eval(p[2-2*w],t)) { // -?-??@@@@@@@-@@@
                return false;
            } // -?@??@@@@@@@-@@@
            if (!C::eval(p[3-w],t)) { // -?@-?@@@@@@@-@@@
                return false;
            } // -?@@?@@@@@@@-@@@
            if (!C::eval(p[3],t)) { // -?@@-@@@@@@@-@@@
                return false;
            } // -?@@@@@@@@@@-@@@
            return true;
        } // -?????@@@@@@@@@@
        return true;
    } // @????????????@@@
    if (!C::eval(p[-3],t)) { // @???????????-@@@
        if (!C::eval(p[2+2*w],t)) { // @?????-?????-@@@
            return false;
        } // @?????@?????-@@@
        if (!C::eval(p[2-2*w],t)) { // @?-???@?????-@@@
            return false;
        } // @?@???@?????-@@@
        if (!C::eval(p[3-w],t)) { // @?@-??@?????-@@@
            return false;
        } // @?@@??@?????-@@@
        if (!C::eval(p[3+w],t)) { // @?@@?-@?????-@@@
            return false;
        } // @?@@?@@?????-@@@
        if (!C::eval(p[3],t)) { // @?@@-@@?????-@@@
            return false;
        } // @?@@@@@?????-@@@
        if (!C::eval(p[1-w3],t)) { // @-@@@@@?????-@@@
            if (!C::eval(p[1+w3],t)) { // @-@@@@@-????-@@@
                return false;
            } // @-@@@@@@????-@@@
            if (!C::eval(p[-3+w],t)) { // @-@@@@@@???--@@@
                return false;
            } // @-@@@@@@???@-@@@
            if (!C::eval(p[w3],t)) { // @-@@@@@@-??@-@@@
                return false;
            } // @-@@@@@@@??@-@@@
            if (!C::eval(p[-1+w3],t)) { // @-@@@@@@@-?@-@@@
                return false;
            } // @-@@@@@@@@?@-@@@
            if (!C::eval(p[-2+2*w],t)) { // @-@@@@@@@@-@-@@@
                return false;
            } // @-@@@@@@@@@@-@@@
            return true;
        } // @@@@@@@?????-@@@
        return true;
    } // @???????????@@@@
    if (!C::eval(p[1-w3],t)) { // @-??????????@@@@
        if (!C::eval(p[1+w3],t)) { // @-?????-????@@@@
            return false;
        } // @-?????@????@@@@
        if (!C::eval(p[-3+w],t)) { // @-?????@???-@@@@
            return false;
        } // @-?????@???@@@@@
        if (!C::eval(p[w3],t)) { // @-?????@-??@@@@@
            return false;
        } // @-?????@@??@@@@@
        if (!C::eval(p[-1+w3],t)) { // @-?????@@-?@@@@@
            return false;
        } // @-?????@@@?@@@@@
        if (!C::eval(p[-2+2*w],t)) { // @-?????@@@-@@@@@
            return false;
        } // @-?????@@@@@@@@@
        return true;
    } // @@??????????@@@@
    if (!C::eval(p[2-2*w],t)) { // @@-?????????@@@@
        if (!C::eval(p[-3+w],t)) { // @@-????????-@@@@
            return false;
        } // @@-????????@@@@@
        if (!C::eval(p[w3],t)) { // @@-?????-??@@@@@
            return false;
        } // @@-?????@??@@@@@
        if (!C::eval(p[-1+w3],t)) { // @@-?????@-?@@@@@
            return false;
        } // @@-?????@@?@@@@@
        if (!C::eval(p[-2+2*w],t)) { // @@-?????@@-@@@@@
            return false;
        } // @@-?????@@@@@@@@
        return true;
    } // @@@?????????@@@@
    if (!C::eval(p[-3+w],t)) { // @@@????????-@@@@
        if (!C::eval(p[3-w],t)) { // @@@-???????-@@@@
            return false;
        } // @@@@???????-@@@@
        if (!C::eval(p[3],t)) { // @@@@-??????-@@@@
            return false;
        } // @@@@@??????-@@@@
        if (!C::eval(p[3+w],t)) { // @@@@@-?????-@@@@
            return false;
        } // @@@@@@?????-@@@@
        return true;
    } // @@@????????@@@@@
    if (!C::eval(p[-2+2*w],t)) { // @@@???????-@@@@@
        if (!C::eval(p[3-w],t)) { // @@@-??????-@@@@@
            return false;
        } // @@@@??????-@@@@@
        if (!C::eval(p[3],t)) { // @@@@-?????-@@@@@
            return false;
        } // @@@@@?????-@@@@@
        return true;
    } // @@@???????@@@@@@
    if (!C::eval(p[3-w],t)) { // @@@-??????@@@@@@
        if (!C::eval(p[-1+w3],t)) { // @@@-?????-@@@@@@
            return false;
        } // @@@-?????@@@@@@@
        return true;
    } // @@@@??????@@@@@@
    return true;
}
template <bool Aligned> void faster_corner_detect_10(const fast_byte* img,
                                                     int img_width, int img_height,
                                                     int img_stride,
                                                     short barrier,
                                                     vector<fast_xy>& corners)
{
    const int w = img_stride; //cfo: before: I.size().x;
    const int stride = 3*w;

    // The compiler refuses to reserve a register for this
    register const __m128i barriers = _mm_set1_epi8((fast_byte)barrier);

    int xend = img_width - 3;
    xend -= (img_width-3) % 16;

    for(int y=3; y < img_height - 3; y++)
    {
        for(int x=3; x < 16; x++)
            if(is_corner_10<Less>( (fast_byte*)img + y*img_stride + x, img_stride, barrier) ||
               is_corner_10<Greater>( (fast_byte*)img + y*img_stride + x, img_stride, barrier))
                corners.push_back(fast_xy(x, y,img[x+y*img_width]));

        for(int x=16; x < xend; x+=16)
        {
            const fast_byte* p = (fast_byte*)img + y*img_stride + x;
            __m128i lo, hi;
            {
                const __m128i here = load_si128<Aligned>((const __m128i*)(p));
                lo = _mm_subs_epu8(here, barriers);
                hi = _mm_adds_epu8(barriers, here);
            }
            unsigned int ans_b, ans_e;
            {
                __m128i top = load_si128<Aligned>((const __m128i*)(p-stride));
                __m128i bottom = load_si128<Aligned>((const __m128i*)(p+stride));

                CHECK_BARRIER(lo, hi, top, ans_b);
                CHECK_BARRIER(lo, hi, bottom, ans_e);
                if (!(ans_b | ans_e))
                    continue;
            }

            unsigned int ans_m, ans_p, possible;
            {
                __m128i ul = _mm_loadu_si128((const __m128i*)(p-2-2*w));
                __m128i lr = _mm_loadu_si128((const __m128i*)(p+2+2*w));
                CHECK_BARRIER(lo, hi, ul, ans_m);
                CHECK_BARRIER(lo, hi, lr, ans_p);
                possible = (ans_m & ans_b) | (ans_e & ans_p);
                if (!possible)
                    continue;
            }

            unsigned int ans_o, ans_n;
            {
                __m128i ll = _mm_loadu_si128((const __m128i*)(p-2+2*w));
                __m128i ur = _mm_loadu_si128((const __m128i*)(p+2-2*w));
                CHECK_BARRIER(lo, hi, ll, ans_o);
                CHECK_BARRIER(lo, hi, ur, ans_n);
                possible &= ans_o | (ans_b & ans_n);
                possible &= ans_n | (ans_e & ans_o);
                if (!possible)
                    continue;
            }

            unsigned int ans_h, ans_k;
            {
                __m128i left = _mm_loadu_si128((const __m128i*)(p-3));
                __m128i right = _mm_loadu_si128((const __m128i*)(p+3));
                CHECK_BARRIER(lo, hi, left, ans_h);
                CHECK_BARRIER(lo, hi, right, ans_k);
                possible &= ans_h | (ans_n & ans_k & ans_p);
                possible &= ans_k | (ans_m & ans_h & ans_o);
                if (!possible)
                    continue;
            }

            unsigned int ans_a, ans_c;
            {
                __m128i a = _mm_loadu_si128((const __m128i*)(p-1-stride));
                __m128i c = _mm_insert_epi16(_mm_srli_si128(a,2), *(const unsigned short*)(p+15-stride), 7);
                //__m128i c = _mm_loadu_si128((const __m128i*)(p+1-stride));
                CHECK_BARRIER(lo, hi, a, ans_a);
                CHECK_BARRIER(lo, hi, c, ans_c);
                possible &= ans_a | (ans_e & ans_p);
                possible &= ans_c | (ans_o & ans_e);
                if (!possible)
                    continue;
            }

            unsigned int ans_d, ans_f;
            {
                __m128i d = _mm_loadu_si128((const __m128i*)(p-1+stride));
                __m128i f = _mm_insert_epi16(_mm_srli_si128(d,2), *(const unsigned short*)(p+15+stride), 7);
                //__m128i f = _mm_loadu_si128((const __m128i*)(p+1+stride));
                CHECK_BARRIER(lo, hi, d, ans_d);
                CHECK_BARRIER(lo, hi, f, ans_f);
                const unsigned int ans_abc = ans_a & ans_b & ans_c;
                possible &= ans_d | (ans_abc & ans_n);
                possible &= ans_f | (ans_m & ans_abc);
                if (!possible)
                    continue;
            }

            unsigned int ans_g, ans_i;
            {
                __m128i g = _mm_loadu_si128((const __m128i*)(p-3-w));
                __m128i ii = _mm_loadu_si128((const __m128i*)(p-3+w));
                CHECK_BARRIER(lo, hi, g, ans_g);
                CHECK_BARRIER(lo, hi, ii, ans_i);
                possible &= ans_g | (ans_f & ans_p & ans_k);
                possible &= ans_i | (ans_c & ans_n & ans_k);
                if (!possible)
                    continue;
            }

            unsigned int ans_j, ans_l;
            {
                __m128i jj = _mm_loadu_si128((const __m128i*)(p+3-w));
                __m128i l = _mm_loadu_si128((const __m128i*)(p+3+w));
                CHECK_BARRIER(lo, hi, jj, ans_j);
                CHECK_BARRIER(lo, hi, l, ans_l);
                const unsigned int ans_ghi = ans_g & ans_h & ans_i;
                possible &= ans_j | (ans_d & ans_o & ans_ghi);
                possible &= ans_l | (ans_m & ans_a & ans_ghi);
                if (!possible)
                    continue;
            }

            possible |= (possible >> 16);
            //if(possible & 0x0f) //Does this make it faster?
            {
                if(possible & (1<< 0))
                    corners.push_back(fast_xy(x + 0, y,img[x+0+y*img_width]));
                if(possible & (1<< 1))
                    corners.push_back(fast_xy(x + 1, y,img[x+1+y*img_width]));
                if(possible & (1<< 2))
                    corners.push_back(fast_xy(x + 2, y,img[x+2+y*img_width]));
                if(possible & (1<< 3))
                    corners.push_back(fast_xy(x + 3, y,img[x+3+y*img_width]));
                if(possible & (1<< 4))
                    corners.push_back(fast_xy(x + 4, y,img[x+4+y*img_width]));
                if(possible & (1<< 5))
                    corners.push_back(fast_xy(x + 5, y,img[x+5+y*img_width]));
                if(possible & (1<< 6))
                    corners.push_back(fast_xy(x + 6, y,img[x+6+y*img_width]));
                if(possible & (1<< 7))
                    corners.push_back(fast_xy(x + 7, y,img[x+7+y*img_width]));
            }
            //if(possible & 0xf0) //Does this mak( ,  fast)r?
            {
                if(possible & (1<< 8))
                    corners.push_back(fast_xy(x + 8, y,img[x+8+y*img_width]));
                if(possible & (1<< 9))
                    corners.push_back(fast_xy(x + 9, y,img[x+9+y*img_width]));
                if(possible & (1<<10))
                    corners.push_back(fast_xy(x +10, y,img[x+10+y*img_width]));
                if(possible & (1<<11))
                    corners.push_back(fast_xy(x +11, y,img[x+11+y*img_width]));
                if(possible & (1<<12))
                    corners.push_back(fast_xy(x +12, y,img[x+12+y*img_width]));
                if(possible & (1<<13))
                    corners.push_back(fast_xy(x +13, y,img[x+13+y*img_width]));
                if(possible & (1<<14))
                    corners.push_back(fast_xy(x +14, y,img[x+14+y*img_width]));
                if(possible & (1<<15))
                    corners.push_back(fast_xy(x +15, y,img[x+15+y*img_width]));
            }
        }

        for(int x=xend; x < img_width - 3; x++)
            if(is_corner_10<Less>(img+y*img_stride+x, img_stride, barrier) ||
               is_corner_10<Greater>(img+y*img_stride+x, img_stride, barrier))
                corners.push_back(fast_xy(x, y,img[x+y*img_width]));
    }
}
static void fast_corner_detect_10_sse2(const fast_byte* img, int img_width,
                                       int img_height, int img_stride,
                                       short barrier,
                                       vector<fast_xy>& corners)
{
    if (img_width<22) {
        fast_corner_detect_10(img,img_width,img_height,img_stride,
                              barrier,corners);
        return;
    }
    else if (img_width<22||img_height<7) {
        return;
    }
    if (is_aligned<16>(img)&&is_aligned<16>(img+img_stride)) {
        faster_corner_detect_10<true>(img,img_width,img_height,img_stride,
                                      barrier,
                                      corners);
    }
    else {
        faster_corner_detect_10<false>(img,img_width,img_height,img_stride,
                                       barrier,
                                       corners);
    }
}
static float shitomasiscore(const unsigned char *img, const int img_w, const int img_h,
                            int u, int v)
{
    float dXX=0.0;
    float dYY=0.0;
    float dXY=0.0;
    int halfbox_size=4;
    int box_size=2*halfbox_size;
    int box_area=box_size*box_size;
    int x_min=u-halfbox_size;
    int x_max=u+halfbox_size;
    int y_min=v-halfbox_size;
    int y_max=v+halfbox_size;

    if (x_min<1||x_max>=img_w-1||y_min<1||y_max>=img_h-1)
        return 0.0; /* patch is too close to the boundary */

    int stride=img_w;
    for (int y=y_min;y<y_max;++y) {

        const uint8_t* ptr_left  =img+stride*y+x_min-1;
        const uint8_t* ptr_right =img+stride*y+x_min+1;
        const uint8_t* ptr_top   =img+stride*(y-1)+x_min;
        const uint8_t* ptr_bottom=img+stride*(y+1)+x_min;
        for (int x=0;x<box_size;++x,++ptr_left,++ptr_right,++ptr_top,++ptr_bottom) {
            float dx=*ptr_right -*ptr_left;
            float dy=*ptr_bottom-*ptr_top;
            dXX+=dx*dx;
            dYY+=dy*dy;
            dXY+=dx*dy;
        }
    }
    /* find and return smaller eigenvalue: */
    dXX=dXX/(2.0f*box_area);
    dYY=dYY/(2.0f*box_area);
    dXY=dXY/(2.0f*box_area);
    return 0.5f*(dXX+dYY-sqrtf((dXX+dYY)*(dXX+dYY)-4.0f*(dXX*dYY-dXY*dXY)));
}
/* FAST corner detect--------------------------------------------------------
 * args:    unsigned char* img  I  image data
 *          int img_w,img_h     I  size of image
 *          int wthstep         I  width step of reading image data
 *          short barrier       I  barrier options
 *          matchopt_t *opt     I  match options
 *          int *uv             O  FAST corners
 *          int *num            O  number of FAST corners
 * return: status (1: ok, 0: fail)
 * --------------------------------------------------------------------------*/
extern int fastfeats(const unsigned char *img,int img_w,int img_h,short barrier,
                     const matchopt_t *opt,int *uv,int *num)
{
    static const float threshold=20.0;
    vector<int> scores,nm_corners;
    vector<fast_xy> corners;
    
    trace(3,"fastfeats:\n");

#ifdef __SSE2__
    /* FAST corner detect */
    fast_corner_detect_10_sse2(img,img_w,img_h,img_w,barrier,corners);
#else
    /* FAST corner detect */
    fast_corner_detect_10(img,img_w,img_h,img_w,20,corners);
#endif
    /* non-maximum suppression */
    fast_corner_score_10(img,img_w,corners,10,scores);
    fast_nonmax_3x3(corners,scores,nm_corners);

    *num=0; /* number of FAST corners */

    /* extract high score features */
    for (auto it:nm_corners) {
        fast_xy& xy=corners.at(it);

        float score=shitomasiscore(img,img_w,img_h,xy.x,xy.y);
        if (score>threshold) {
            uv[3**num+0]=xy.x;
            uv[3**num+1]=xy.y;
            uv[3**num+2]=xy.val; (*num)++;
        }
    }
    return *num;
}
