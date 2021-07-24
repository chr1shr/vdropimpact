#ifndef VDI_COMMON_HH
#define VDI_COMMON_HH

#include <cstdio>

#ifdef _OPENMP
#include "omp.h"
#else
#include <ctime>
#endif

void fatal_error(const char *p,int code);
FILE* safe_fopen(const char* filename,const char* mode);
void safe_fread(void *ptr,size_t size,size_t count,FILE *fp,const char* p);

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
inline double wtime() {return omp_get_wtime();}
#else
inline double wtime() {return (1./CLOCKS_PER_SEC)*static_cast<double>(clock());}
#endif

#endif
