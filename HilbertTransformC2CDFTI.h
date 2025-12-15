/* vim: set fdm=marker:
 *
 *
 * HilbertTransform optimized using c2c DFTI
 * Uses size_t to index large arrays
 *
 *
 *
 * Note
 * INPUT: -(1 / PI) * Im{u(t)}
 * OUTPUT: Re{u(t)}
 *
 *
 *
 */

#ifndef _HILBERTTRANSFORMC2CDFTI_H
#define _HILBERTTRANSFORMC2CDFTI_H

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "mkl.h"
#include "mkl_dfti.h"
#include "mkl_service.h"

#include "C_extension.h"


 
#define MKL_Complex16 double complex 


MKL_Complex16 hilbert_sign_fct(int i, size_t *Npad) {
    return -I * M_PI * intSgn(i) * intSign( (int) NpadHT[2] / 2 - i - 1);
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
    size_t *N    = p->N,
           *pad  = p->pad,
           *Npad = p->Npad;

    size_t i, j, k;
    size_t N3    = N[0] * N[1] * N[2],
           N3pad = N[0] * N[1] * Npad[2];

    MKL_Complex16 *hilbertsign = mkl_calloc(Npad[2], sizeof(MKL_Complex16), 64),
                  *big =         mkl_calloc(N3pad,   sizeof(MKL_Complex16), 64);


//    for (i = 0; i < Npad[2]; i++) {
//        hilbertsign[i] = -I * M_PI * intSign(i) * intSgn((int) N[2] / 2 + (int) pad[2] - i - 1);
//    }


    // pad double to double complex array 3d
    size_t n, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                m = k + pad[2] + Npad[2] * (j + N[1] * i); 
                big[m] = in[n];
            }
        }
    }

    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    MKL_LONG size = Npad[2];
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 1, size);
    status = DftiCommitDescriptor(hand);

    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = Npad[2] * (j + N[1] * i);
            status = DftiComputeForward(hand, &big[n]);
            for (k = 0; k < Npad[2]; k++) {
                big[k + n] *= hilbert_sign_fct(k, Npad);
            }
            status = DftiComputeBackward(hand, &big[n]);
        } 
    }
    status = DftiFreeDescriptor(&hand);
    // unpad double complex to double array 3d and normalize
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                m = k + pad[2] + Npad[2] * (j + N[1] * i); 
                out[n] = big[m] / Npad[2];
            }
        }
    }

    mkl_free(big);
    mkl_free(hilbertsign);

    return status;
}


#endif /* _HILBERTTRANSFORMC2CDFTI_H */

