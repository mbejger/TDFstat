#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <complex.h>
#include <omp.h>

#include "settings.h"

#if defined(SLEEF)
// uncomment to use sleef headers with inline versions of functions
/*
#define SLEEF_ALWAYS_INLINE __attribute__((always_inline))
#define SLEEF_INLINE static inline
#define SLEEF_CONST const
#include <x86intrin.h>
#include <float.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
*/

#include <sleef.h>
#define SLEEF_ENABLE_OMP_SIMD
#ifdef USE_AVX2
//#include <sleefinline_avx2.h>
#define VLEN 4
typedef __m256d vdouble;
typedef Sleef___m256d_2 vdouble2;
#define v_sincos Sleef_sincosd4_u35avx2
#endif
#ifdef USE_AVX512
//#include <sleefinline_avx512f.h>
#define VLEN 8
typedef __m512d vdouble;
typedef Sleef___m512d_2 vdouble2;
//typedef vdouble2_avx512f_sleef vdouble2;
#define v_sincos Sleef_sincosd8_u35avx512f
#endif
#ifdef USE_AVX2_SINCOS
//#include <sleefinline_avx2.h>
#define VLEN 4
typedef __m256d vdouble;
typedef Sleef___m256d_2 vdouble2;
static vdouble vload_vd_p(const double *ptr) { return _mm256_load_pd(ptr); }
static void vstore_v_p_vd(double *ptr, vdouble v) { _mm256_store_pd(ptr, v); }
#define v_sincos Sleef_sincosd4_u35
#endif

#endif


void spindown_modulation(int nifo, int N,
			 double het1, double *sgnlt, double _tmp1[][N],
			 fftw_complex *fxa, fftw_complex *fxb) {
  int i, j, n;
  
  for(n=0; n<nifo; ++n) {
    
#if defined(SLEEF)
    // use simd sincos from the SLEEF library;
    int Nvec = (N/VLEN)*VLEN; // fully vectorizable range
#if defined(USE_AVX2)
#pragma omp parallel for schedule(static) default(shared) private(j)
    for (i=0; i<Nvec; i+=VLEN) {
      vdouble2 v;
      vdouble a;
      a = _mm256_setr_pd( het1*(i)   + sgnlt[1]*_tmp1[n][i],
			  het1*(i+1) + sgnlt[1]*_tmp1[n][i+1],
			  het1*(i+2) + sgnlt[1]*_tmp1[n][i+2],
			  het1*(i+3) + sgnlt[1]*_tmp1[n][i+3]
			  );
      v = v_sincos(a);
      // exph = _c[j] - I*_p[j];
      // fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph;
      __m256d mare = _mm256_setr_pd( creal(ifo[n].sig.xDatma[i]),
				     creal(ifo[n].sig.xDatma[i+1]),
				     creal(ifo[n].sig.xDatma[i+2]),
				     creal(ifo[n].sig.xDatma[i+3]) );
      __m256d maim = _mm256_setr_pd( cimag(ifo[n].sig.xDatma[i]),
				     cimag(ifo[n].sig.xDatma[i+1]),
				     cimag(ifo[n].sig.xDatma[i+2]),
				     cimag(ifo[n].sig.xDatma[i+3]) );
      __m256d vec1, vre, vim;
      vec1 = _mm256_mul_pd(mare, v.y);
      vre = _mm256_fmadd_pd(maim, v.x, vec1);
      vec1 = _mm256_mul_pd(mare, v.x);
      vim = _mm256_fmsub_pd(maim, v.y, vec1);
      if (n==0)
	for(j=0; j<VLEN; ++j)
	  fxa[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
      else
	for(j=0; j<VLEN; ++j)
	  fxa[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];
      
      mare = _mm256_setr_pd( creal(ifo[n].sig.xDatmb[i]),
			     creal(ifo[n].sig.xDatmb[i+1]),
			     creal(ifo[n].sig.xDatmb[i+2]),
			     creal(ifo[n].sig.xDatmb[i+3]) );
      maim = _mm256_setr_pd( cimag(ifo[n].sig.xDatmb[i]),
			     cimag(ifo[n].sig.xDatmb[i+1]),
			     cimag(ifo[n].sig.xDatmb[i+2]),
			     cimag(ifo[n].sig.xDatmb[i+3]) );
      vec1 = _mm256_mul_pd(mare, v.y);
      vre = _mm256_fmadd_pd(maim, v.x, vec1);
      vec1 = _mm256_mul_pd(mare, v.x);
      vim = _mm256_fmsub_pd(maim, v.y, vec1);
      if (n==0)
	for(j=0; j<VLEN; ++j)
	  fxb[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
      else
	for(j=0; j<VLEN; ++j)
	  fxb[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];
    } // i
#elif defined(USE_AVX512)
#pragma omp parallel for schedule(static) default(shared) private(j)
    for (i=0; i<Nvec; i+=VLEN) {
      vdouble2 v;
      vdouble a;
      a = _mm512_setr_pd( het1*(i)   + sgnlt[1]*_tmp1[n][i],
			  het1*(i+1) + sgnlt[1]*_tmp1[n][i+1],
			  het1*(i+2) + sgnlt[1]*_tmp1[n][i+2],
			  het1*(i+3) + sgnlt[1]*_tmp1[n][i+3],
			  het1*(i+4) + sgnlt[1]*_tmp1[n][i+4],
			  het1*(i+5) + sgnlt[1]*_tmp1[n][i+5],
			  het1*(i+6) + sgnlt[1]*_tmp1[n][i+6],
			  het1*(i+7) + sgnlt[1]*_tmp1[n][i+7]
			  );
      v = v_sincos(a);
      // exph = _c[j] - I*_p[j];
      // fxa[i+j] = ifo[0].sig.xDatma[i+j]*exph;
      __m512d mare = _mm512_setr_pd( creal(ifo[n].sig.xDatma[i]),
				     creal(ifo[n].sig.xDatma[i+1]),
				     creal(ifo[n].sig.xDatma[i+2]),
				     creal(ifo[n].sig.xDatma[i+3]),
				     creal(ifo[n].sig.xDatma[i+4]),
				     creal(ifo[n].sig.xDatma[i+5]),
				     creal(ifo[n].sig.xDatma[i+6]),
				     creal(ifo[n].sig.xDatma[i+7])
				     );
      __m512d maim = _mm512_setr_pd( cimag(ifo[n].sig.xDatma[i]),
				     cimag(ifo[n].sig.xDatma[i+1]),
				     cimag(ifo[n].sig.xDatma[i+2]),
				     cimag(ifo[n].sig.xDatma[i+3]),
				     cimag(ifo[n].sig.xDatma[i+4]),
				     cimag(ifo[n].sig.xDatma[i+5]),
				     cimag(ifo[n].sig.xDatma[i+6]),
				     cimag(ifo[n].sig.xDatma[i+7])
				     );
      __m512d vec1, vre, vim;
      vec1 = _mm512_mul_pd(mare, v.y);
      vre = _mm512_fmadd_pd(maim, v.x, vec1);
      vec1 = _mm512_mul_pd(mare, v.x);
      vim = _mm512_fmsub_pd(maim, v.y, vec1);
      if (n==0)
	for(j=0; j<VLEN; ++j)
	  fxa[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
      else
	for(j=0; j<VLEN; ++j)
	  fxa[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];
      
      mare = _mm512_setr_pd( creal(ifo[n].sig.xDatmb[i]),
			     creal(ifo[n].sig.xDatmb[i+1]),
			     creal(ifo[n].sig.xDatmb[i+2]),
			     creal(ifo[n].sig.xDatmb[i+3]),
			     creal(ifo[n].sig.xDatmb[i+4]),
			     creal(ifo[n].sig.xDatmb[i+5]),
			     creal(ifo[n].sig.xDatmb[i+6]),
			     creal(ifo[n].sig.xDatmb[i+7])
			     );
      maim = _mm512_setr_pd( cimag(ifo[n].sig.xDatmb[i]),
			     cimag(ifo[n].sig.xDatmb[i+1]),
			     cimag(ifo[n].sig.xDatmb[i+2]),
			     cimag(ifo[n].sig.xDatmb[i+3]),
			     cimag(ifo[n].sig.xDatmb[i+4]),
			     cimag(ifo[n].sig.xDatmb[i+5]),
			     cimag(ifo[n].sig.xDatmb[i+6]),
			     cimag(ifo[n].sig.xDatmb[i+7])
			     );
      vec1 = _mm512_mul_pd(mare, v.y);
      vre = _mm512_fmadd_pd(maim, v.x, vec1);
      vec1 = _mm512_mul_pd(mare, v.x);
      vim = _mm512_fmsub_pd(maim, v.y, vec1);
      if (n==0)
	for(j=0; j<VLEN; ++j)
	  fxb[i+j] = ((double *)&vre)[j] + I*((double *)&vim)[j];
      else
	for(j=0; j<VLEN; ++j)
	  fxb[i+j] += ((double *)&vre)[j] + I*((double *)&vim)[j];
    } // i

#elif defined(USE_AVX2_SINCOS)
    // use only SIMD sincos from sleef
#pragma omp parallel for schedule(static) default(shared) private(j)    
    for (i=0; i<Nvec; i+=VLEN) {
      double _p[VLEN] __attribute__((aligned(32))),
	     _c[VLEN] __attribute__((aligned(32)));
      complex double exph;
      vdouble2 v;
      vdouble a;
      for(j=0; j<VLEN; ++j)
        _p[j] =  het1*(i+j) + sgnlt[1]*_tmp1[n][i+j];
      a = vload_vd_p(_p);
      v = v_sincos(a);
      vstore_v_p_vd(_p, v.x); // reuse _p for sin
      vstore_v_p_vd(_c, v.y);
      
      if (n==0)
	for(j=0; j<VLEN; ++j){
	  exph = _c[j] - I*_p[j];
	  fxa[i+j] = ifo[n].sig.xDatma[i+j]*exph;
	  fxb[i+j] = ifo[n].sig.xDatmb[i+j]*exph;
	}
      else
	for(j=0; j<VLEN; ++j){
	  exph = _c[j] - I*_p[j];
	  fxa[i+j] += ifo[n].sig.xDatma[i+j]*exph;
	  fxb[i+j] += ifo[n].sig.xDatmb[i+j]*exph;
	}
    } // i
#endif
    
    // for all vectorized versions:
    // calculate remaining elements if N is not a multiple of VLEN
    for(i=Nvec; i<N; ++i) {
      double phase, sp, cp;
      complex double exph;
      phase = het1*i + sgnlt[1]*_tmp1[n][i];
      sincos(phase, &sp, &cp);
      exph = cp - I*sp;
      if (n==0){
	fxa[i] = ifo[n].sig.xDatma[i]*exph;
	fxb[i] = ifo[n].sig.xDatmb[i]*exph;
      } else {
	fxa[i] += ifo[n].sig.xDatma[i]*exph;
	fxb[i] += ifo[n].sig.xDatmb[i]*exph;
      }	    
    } // i
    
#elif defined(GNUSINCOS)
    double phase, cp, sp;
    complex double exph;
#pragma omp parallel for schedule(static, 1) private(phase, exph, cp, sp) default(shared)
    for(i=N-1; i>-1; --i) {
      phase = het1*i + sgnlt[1]*_tmp1[n][i];
      sincos(phase, &sp, &cp);
      exph = cp - I*sp;
      if (n==0){
	fxa[i] = ifo[n].sig.xDatma[i]*exph;
	fxb[i] = ifo[n].sig.xDatmb[i]*exph;
      } else {
	fxa[i] += ifo[n].sig.xDatma[i]*exph;
	fxb[i] += ifo[n].sig.xDatmb[i]*exph;
      }	    
      
    }
#else
    // generic version for systems without sincos function
    double phase, cp, sp;
    complex double exph;
#pragma omp parallel for schedule(static, 1) private(phase, exph, cp, sp) default(shared)
    for(i=N-1; i!=-1; --i) {
      phase = het1*i + sgnlt[1]*_tmp1[n][i];
      cp = cos(phase);
      sp = sin(phase);
      exph = cp - I*sp;
      if (n==0){
	fxa[i] = ifo[n].sig.xDatma[i]*exph;
	fxb[i] = ifo[n].sig.xDatmb[i]*exph;
      } else {
	fxa[i] += ifo[n].sig.xDatma[i]*exph;
	fxb[i] += ifo[n].sig.xDatmb[i]*exph;
      }	    
    }
#endif

  } // n
  
} // phase_mod
