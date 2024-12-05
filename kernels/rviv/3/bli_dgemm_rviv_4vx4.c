/*

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2023, The University of Texas at Austin

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


*/
#include "bli_rviv_utils.h"

////// HERO_1 includes /////
/*#ifdef __HERO_1
#include "encoding.h"
#include "inttypes.h"
#include "matvec_dev.h"
extern volatile uint32_t dma_wait_cycles;
////// HOST includes /////
#else
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <unistd.h>

#include <libhero/hero_api.h>
#include <omp.h>
#define snrt_printf printf

static inline void fence() { asm volatile("fence" ::: "memory"); }
*/
void bli_dgemm_rviv_asm_4vx4
    (
             intptr_t   k,
       const void*      alpha,
       const void*      a,
       const void*      b,
       const void*      beta,
             void*      c, intptr_t rs_c, intptr_t cs_c
    );

void bli_dgemm_rviv_4vx4
     (
             dim_t      m,
             dim_t      n,
             dim_t      k,
       const void*      alpha,
       const void*      a,
       const void*      b,
       const void*      beta,
             void*      c, inc_t rs_c, inc_t cs_c,
       const auxinfo_t* data,
       const cntx_t*    cntx
     )
{

	// The assembly kernels always take native machine-sized integer arguments.
	// dim_t and inc_t are normally defined as being machine-sized. If larger, assert.
	bli_static_assert( sizeof(dim_t) <= sizeof(intptr_t) &&
	                   sizeof(inc_t) <= sizeof(intptr_t) );

	// Extract vector-length dependent mr, nr that are fixed at configure time.
	const inc_t mr = bli_cntx_get_blksz_def_dt( BLIS_DOUBLE, BLIS_MR, cntx );
	const inc_t nr = 4;

	GEMM_UKR_SETUP_CT( d, mr, nr, false );

	// The kernel assumes rs_c == 1, and the context should not deviate from it.
	assert( rs_c == 1 );

/*#pragma omp target device(1) map(to : alpha, a, b, beta, c, k, cs_c)
    {
        volatile *void a_ = a;
        volatile *void b_ = b;
        volatile *void alpha_ = alpha;
        volatile *void beta_ = beta;
        volatile *void c_ = c;
	volatile dim_t k_ = k;

//	bli_dgemm_rviv_asm_4vx4( k, alpha, a, b, beta, c,
//	    			get_vlenb(), cs_c * sizeof(double) );
*/
	const int M = 4; // Assuming 4x4 blocks as per the assembly
        const int N = 4; // Assuming 4x4 blocks as per the assembly
        //double XC[4][4] = {{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0}}; // Accumulator for the result
   	
       	double * XA = (double *)a;
       	double * XB = (double *)b;
       	double * alpha_ = (double *)alpha;
       	double * beta_ = (double *)beta;
	double * XC = (double *)c;
	int lda = M;
	int ldb = k;
	int ldc;
	#ifdef DGEMM_MR
    	ldc = ( ( m - 1 ) / DGEMM_MR + 1 ) * DGEMM_MR;
	#else
    	ldc     = m;
	#endif

	for ( int j = 0; j < N; j ++ ) {
                for ( int i = 0; i < M; i ++ ) {
                        for ( int p = 0; p < k; p ++ ) {
                                XC[ j * ldc + i ] *= *beta_; 
				XA[ p * lda + i ] *= *alpha_;
			}
		}
	}
	for ( int j = 0; j < N; j ++ ) {
        	for ( int i = 0; i < M; i ++ ) {
            		for ( int p = 0; p < k; p ++ ) {
                		XC[ j * ldc + i ] += XA[ p * lda + i ] * XB[ j * ldb + p ];
			}
		}
	}
/*		
    // Matrix multiplication
    for (intptr_t l = 0; l < k; l++) {
        // Compute each block row of A and corresponding B row
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                AB[i][j] += a_[i + l * M] * b_[l * N + j];
            }
        }
    }

    // Multiply accumulators by alpha
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            AB[i][j] *= *alpha_;
        }
    }

    // Apply beta and store in C
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            c_[i * rs_c + j * cs_c] = AB[i][j] + (*beta_ * c_[i * rs_c + j * cs_c]);
        }
    }
}
*/	
//}
	GEMM_UKR_FLUSH_CT( d );
}
