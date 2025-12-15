/* vim: set fdm=marker:
 *
 *
 * HilbertTransform optimized using r2c and c2r DFTI
 * Warning: Slight Faster than HilbertTransformCRCDFTIOpt and uses only 3/4 of the memory
 * but the program randomly gets hung up and waste time
 *
 * Uses size_t to index large arrays
 *
 * Note
 * INPUT: -(1 / PI) * Im{u(t)}
 * OUTPUT: Re{u(t)}
 *
 */

#ifndef _HILBERTTRANSFORMR2CDFTI_H
#define _HILBERTTRANSFORMR2CDFTI_H

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "mkl.h"
#include "mkl_dfti.h"
#include "mkl_service.h"

#include "C_extension.h"


 
#define MKL_Complex16 double complex 


void Hilbert_sign_fct(void *in, size_t *Nhalfcmplx) 
{
     double complex *product = (double complex *) in;
     int k;
     for (k = 0; k < Nhalfcmplx[2]; k++) {
        product[k] *= -I * M_PI * intSgn(k) * intSign((int) Nhalfcmplx[2] - k - 2);
    }
}

struct HilbertParams
{
    size_t *N,
           *padHT,
           *NpadHT;
};


int Hilbert_transform(double *in, double *out, void *params)
{
    // ON -> 1, OFF -> 0 defined C_extension.h
    int PRINT = OFF,
        CHECKSUM = OFF;
    if (PRINT) {
        printf("PRINT ON\n");
    }
    if (CHECKSUM) {
        printf("CHECKSUM ON\n");
    }


    struct HilbertParams *p = (struct HilbertParams *) params;
    size_t *N      = p->N,
           *padHT  = p->padHT,
           *NpadHT = p->NpadHT;

    size_t i, j, k, l, n, m;
    
    size_t Nhalfcmplx[3] = {N[0], N[1], NpadHT[2] / 2 + 1},
           Nhalfreal[3]  = {N[0], N[1], 2 * (NpadHT[2] / 2 + 1)};
    size_t N3          = N[0] * N[1] * N[2],
           N3halfcmplx = N[0] * N[1] * Nhalfcmplx[2],
           N3padHT     = N[0] * N[1] * NpadHT[2],
           N3halfreal  = N[0] * N[1] * Nhalfreal[2];


    double *bigin  = calloc(N3halfreal, sizeof(double));


    // pad array double 3d
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                l = k + padHT[2] + Nhalfreal[2] * m; 
                bigin[l] = in[n];
            }
        }
    }

    MKL_LONG status = 0,
             dim = 1,
             size = NpadHT[2],
             size_real = NpadHT[2] + 2, // size of real domain
             size_cmplx = NpadHT[2] / 2 + 1; // size of complex domain
    DFTI_DESCRIPTOR_HANDLE hand_forward,
                           hand_backward;

    double scale_forwards = 1,
           scale_backwards = 1 / (double) size;
    status = DftiCreateDescriptor(&hand_forward, DFTI_DOUBLE, DFTI_REAL, dim, size);
    status = DftiSetValue(hand_forward, DFTI_PLACEMENT, DFTI_INPLACE);
    status = DftiSetValue(hand_forward, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(hand_forward, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
    status = DftiSetValue(hand_forward, DFTI_FORWARD_SCALE, scale_forwards);
    status = DftiCommitDescriptor(hand_forward);

    status = DftiCreateDescriptor(&hand_backward, DFTI_DOUBLE, DFTI_REAL, dim, size);
    status = DftiSetValue(hand_backward, DFTI_PLACEMENT, DFTI_INPLACE);
    status = DftiSetValue(hand_backward, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(hand_backward, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
    status = DftiSetValue(hand_backward, DFTI_BACKWARD_SCALE, scale_backwards);
    status = DftiCommitDescriptor(hand_backward);


    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = Nhalfreal[2] * (j + N[1] * i); // stride length between domains
            status = DftiComputeForward(hand_forward, &bigin[m]);
            Hilbert_sign_fct(&bigin[m], Nhalfcmplx);
            status = DftiComputeBackward(hand_backward, &bigin[m]);
        } 
    }
    status = DftiFreeDescriptor(&hand_forward);
    status = DftiFreeDescriptor(&hand_backward);


    // unpad double array 3d and normalize
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                l = k + padHT[2] + Nhalfreal[2] * m; 
                out[n] = bigin[l];
            }
        }
    }


    free(bigin);

    return status;
}


#endif /* _HILBERTTRANSFORMR2CDFTIOPTII_H */

