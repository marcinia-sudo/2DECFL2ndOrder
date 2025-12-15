/* vim: set fdm=marker:
 *
 *
 *
 * ConvolutionchiBubWW.h
 *
 *
 *
 */


#ifndef _CONVOLUTIONCHIBUBWW_H
#define _CONVOLUTIONCHIBUBWW_H

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <mkl.h>
#include <mkl_dfti.h>
#include <mkl_service.h>

#include "C_extension.h"
#include "Functionschi.h"


#define MKL_Complex16 double complex

// convolution products afermi * conj(bfermi) - conj(afermi) * bfermi
void product_chiBubWW(void *in1,
                      void *in2,
                      void *in3,
                      void *in4,
                      void *out,
                      size_t *N)
{
    // for a full complex array of complex type the dimensions of the array are {N[0], N[1], N[2]}
    // for a half complex array of complex type the dimensions of the array are {N[0], N[1], N[2] / 2 + 1}
    double complex *a      = (double complex *) in1,
                   *bfermi = (double complex *) in2,
                   *afermi = (double complex *) in3,
                   *b      = (double complex *) in4,
                   *product = (double complex *) out;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) { // (Npad[2] / 2 + 1)
                n = k + N[2] * (j + N[1] * i);
                // conj[ a(p) ] * bf(p+q) - conj[ af(p) ] * b(p+q)
                product[n] = (conj(a[n]) * bfermi[n] - conj(afermi[n]) * b[n]) * pow(-1, k);
            }
        }
    }
}


struct ChiBubWWParams
{
    size_t *N,
           *padFT,
           *NpadFT;
    double DELTAomega;
    double *en,
           *Fermi_array,
           *rho_bigG,
           *rho_chiBubWW;
};


MKL_INT chiBubWW_function(void *params)
{

    // ON = 1, OFF = 0 defined in C_extention.h
    int PRINT = OFF,
        CHECKSUM = OFF;
    if (PRINT == ON) {
        printf("PRINT ON\n");
    }
    if (CHECKSUM == ON) {
        printf("CHECKSUM ON\n");
    }

    struct ChiBubWWParams *p = (struct ChiBubWWParams *) params;
    size_t *N = p->N,
           *padFT = p->padFT,
           *NpadFT = p->NpadFT;
    double DELTAomega = p->DELTAomega;
    double *en = p->en,
           *fermi = p->Fermi_array,
           *rho_bigG = p->rho_bigG,
           *rho_chiBubWW = p->rho_chiBubWW;

    size_t i, j, k, l, m, n;


    size_t Nhalfcmplx[3]  = {N[0], N[1], (NpadFT[2] / 2 + 1)}, // complex array
           Nhalfreal[3] = {N[0], N[1], 2 * (NpadFT[2] / 2 + 1)}; // real to complex array
    size_t Ns      = N[0] * N[1],
           N3      = N[0] * N[1] * N[2],
           N3halfcmplx  = N[0] * N[1] * Nhalfcmplx[2],
           N3padFT   = N[0] * N[1] * NpadFT[2],
           N3halfreal = N[0] * N[1] * Nhalfreal[2];

    if (PRINT) {
        printf("%zu %zu %zu %zu\n", N3, N3halfcmplx, N3padFT, N3halfreal);
    }


    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand_forward = 0,
                           hand_backward = 0;
    char version[DFTI_VERSION_LENGTH];
    DftiGetValue(0, DFTI_VERSION, version);
    if (PRINT) {
        printf(version);
        printf("\n");
    }

    MKL_LONG halfcmplx_strides[4]  = {0, Nhalfcmplx[2]  * Nhalfcmplx[1],  Nhalfcmplx[2],  1}, // strides for half-complex padded-array
             full_strides[4]  = {0, NpadFT[2]   * NpadFT[1],   NpadFT[2],   1}, // strides for full padded-array
             halfreal_strides[4] = {0, Nhalfreal[2] * Nhalfreal[1], Nhalfreal[2], 1}; // strides for half-real padded-array


    double scale_forwards = 1,
           scale_backwards = 2 * (DELTAomega) / (Ns * N3padFT);

    status = DftiCreateDescriptor(&hand_forward, DFTI_DOUBLE, DFTI_REAL, 3, NpadFT);
    if (status != DFTI_NO_ERROR) goto failed_forward;
    status = DftiSetValue(hand_forward, DFTI_PLACEMENT, DFTI_INPLACE); // Default 
    if (status != DFTI_NO_ERROR) goto failed_forward;
    status = DftiSetValue(hand_forward, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed_forward;
    status = DftiSetValue(hand_forward, DFTI_FORWARD_SCALE, scale_forwards);
    if (status != DFTI_NO_ERROR) goto failed_forward;
    status = DftiSetValue(hand_forward, DFTI_INPUT_STRIDES, halfreal_strides);
    if (status != DFTI_NO_ERROR) goto failed_forward;
    status = DftiSetValue(hand_forward, DFTI_OUTPUT_STRIDES, halfcmplx_strides);
    if (status != DFTI_NO_ERROR) goto failed_forward;
    status = DftiCommitDescriptor(hand_forward);
    if (status != DFTI_NO_ERROR) goto failed_forward;


    status = DftiCreateDescriptor(&hand_backward, DFTI_DOUBLE, DFTI_REAL, 3, NpadFT);
    if (status != DFTI_NO_ERROR) goto failed_backward;
    status = DftiSetValue(hand_backward, DFTI_PLACEMENT, DFTI_INPLACE); // Default 
    if (status != DFTI_NO_ERROR) goto failed_backward;
    status = DftiSetValue(hand_backward, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != DFTI_NO_ERROR) goto failed_backward;
    status = DftiSetValue(hand_backward, DFTI_BACKWARD_SCALE, scale_backwards);
    if (status != DFTI_NO_ERROR) goto failed_backward;
    status = DftiSetValue(hand_backward, DFTI_INPUT_STRIDES, halfcmplx_strides);
    if (status != DFTI_NO_ERROR) goto failed_backward;
    status = DftiSetValue(hand_backward, DFTI_OUTPUT_STRIDES, halfreal_strides);
    if (status != DFTI_NO_ERROR) goto failed_backward;
    status = DftiCommitDescriptor(hand_backward);
    if (status != DFTI_NO_ERROR) goto failed_backward;


    double *out1 = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out2 = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out3 = mkl_calloc(N3halfreal, sizeof(double), 64);

    double *rho_bigG_forward           = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_rho_bigG_forward        = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en2_rho_bigG_forward       = mkl_calloc(N3halfreal, sizeof(double), 64),
           *rho_bigG_Fermi_forward     = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_rho_bigG_Fermi_forward  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en2_rho_bigG_Fermi_forward = mkl_calloc(N3halfreal, sizeof(double), 64);


    double rho_bigG_temp;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                rho_bigG_temp = rho_bigG[n];
                rho_bigG_forward[l] = rho_bigG_temp;
                en_rho_bigG_forward[l] = en[m] * rho_bigG_temp;
                en2_rho_bigG_forward[l] = en[m] * en[m] * rho_bigG_temp;
                rho_bigG_Fermi_forward[l] = rho_bigG_temp * fermi[k];
                en_rho_bigG_Fermi_forward[l] = en[m] * rho_bigG_temp * fermi[k];
                en2_rho_bigG_Fermi_forward[l] = en[m] * en[m] * rho_bigG_temp * fermi[k];
            }
        }
    }



    status = DftiComputeForward(hand_forward, rho_bigG_forward);
    if (status != DFTI_NO_ERROR) goto failed_1;
    status = DftiComputeForward(hand_forward, en_rho_bigG_forward);
    if (status != DFTI_NO_ERROR) goto failed_1;
    status = DftiComputeForward(hand_forward, en2_rho_bigG_forward);
    if (status != DFTI_NO_ERROR) goto failed_1;

   
    
    status = DftiComputeForward(hand_forward, rho_bigG_Fermi_forward);
    if (status != DFTI_NO_ERROR) goto failed_1;
    status = DftiComputeForward(hand_forward, en_rho_bigG_Fermi_forward);
    if (status != DFTI_NO_ERROR) goto failed_1;
    status = DftiComputeForward(hand_forward, en2_rho_bigG_Fermi_forward);
    if (status != DFTI_NO_ERROR) goto failed_1;


    // conj[a(p)] * bf(p+q) - conj[af(p)] * b(p+q)
    product_chiBubWW(en2_rho_bigG_forward, rho_bigG_Fermi_forward, en2_rho_bigG_Fermi_forward, rho_bigG_forward, out1, Nhalfcmplx);
    product_chiBubWW(rho_bigG_forward, en2_rho_bigG_Fermi_forward, rho_bigG_Fermi_forward, en2_rho_bigG_forward, out2, Nhalfcmplx);
    product_chiBubWW(en_rho_bigG_forward, en_rho_bigG_Fermi_forward, en_rho_bigG_Fermi_forward, en_rho_bigG_forward, out3, Nhalfcmplx);

    status = DftiComputeBackward(hand_backward, out1);
    if (status != DFTI_NO_ERROR) goto failed_2;
    status = DftiComputeBackward(hand_backward, out2);
    if (status != DFTI_NO_ERROR) goto failed_2;
    status = DftiComputeBackward(hand_backward, out3);
    if (status != DFTI_NO_ERROR) goto failed_2;
 



// padFT = 3 Nw / 2
// NpadFT =  2 * padFT + Nw
//    double scale1 = 2 * (DELTAomega) / (Ns * N3padFT);
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                rho_chiBubWW[n] = (out1[l] + out2[l] - 2 * out3[l]);
            }
        }
    }
    
    cleanup:
 
    mkl_free(out1);
    mkl_free(out2);
    mkl_free(out3);
    mkl_free(rho_bigG_Fermi_forward);
    mkl_free(en_rho_bigG_Fermi_forward);
    mkl_free(en2_rho_bigG_Fermi_forward);
    mkl_free(rho_bigG_forward);
    mkl_free(en_rho_bigG_forward);
    mkl_free(en2_rho_bigG_forward);


    cleanup_backward:
    DftiFreeDescriptor(&hand_backward);

    cleanup_forward:
    DftiFreeDescriptor(&hand_forward);


    return status;

    failed_forward:
    printf("chiBubWW\n");
    printf("ERROR: %s\n", DftiErrorMessage(status));
    printf("Forward ERROR, status = %li\n", status); status = 1; goto cleanup_forward;

    failed_backward:
    printf("chiBubWW\n");
    printf("ERROR: %s\n", DftiErrorMessage(status));
    printf("Backward ERROR, status = %li\n", status); status = 1; goto cleanup_backward;
 
    failed_1:
    printf("chiBubWW\n");
    printf("ERROR: %s\n", DftiErrorMessage(status));
    printf("failed_1 ERROR, status = %li\n", status); status = 1; goto cleanup;
 
    failed_2:
    printf("chiBubWW\n");
    printf("ERROR: %s\n", DftiErrorMessage(status));
    printf("failed_2 ERROR, status = %li\n", status); status = 1; goto cleanup;
 

}


#endif
