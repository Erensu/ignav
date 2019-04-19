/* ---------------------------------------------------------------------------
** The MIT license:
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is furnished
** to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all
** copies or substantial portions of the Software.

** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
** CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**----------------------------------------------------------------------------*/
#include "sse2neonbinding.h"
#ifdef WIN32
#include <xmmintrin.h>
#include <emmintrin.h>
#include <malloc.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#include <stdio.h>
#endif
#ifdef WIN32
void* platformAlignedAlloc(size_t size)
    {
        return _aligned_malloc(size,16);
    }
    void platformAlignedFree(void* ptr)
    {
        _aligned_free(ptr);
    }
#else
void* platformAlignedAlloc(size_t size)
{
    //return ::memalign(16,size);
    void *address;
    int ret=posix_memalign(&address,16,size);
    if (ret!=0) {
        fprintf(stderr,"Error at File %s line number %d\n",__FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    return address;
}
void platformAlignedFree(void* ptr)
{
    //::free(ptr);
    free(ptr);
}
#endif
