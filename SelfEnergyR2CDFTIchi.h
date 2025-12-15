/* vim: set fdm=marker:
 *
 *
 * SelfEnergyR2CDFTIchi.h
 *
 * This module computes the self energies of the 2d second order ECFL theory. The self energies are found by computing
 * a set of double convolution integrals. This is achieved using convolution thereom which invovles taking an FFT of
 * the real signal, computing the product in frequency space and applying inverse FFT to arrive at the solution.
 *
 * con(p, q, p+q-k) = -ggg
 *
 * Here we compute chi[1] and psi[1] of the lambda expansion
 */


#ifndef _SELFENERGYR2CDFTICHI_H
#define _SELFENERGYR2CDFTICHI_H

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <mkl.h>
#include <mkl_dfti.h>
#include <mkl_service.h>

#include "C_extension.h"


#define MKL_Complex16 double complex

// convolution products afermi * bfermi conj(cfermibar) + afermibar * bfermibar * conj(cfermi)
void product1(void *in1,
              void *in2,
              void *in3,
              void *in4,
              void *in5,
              void *in6,
              void *out,
              size_t *N)
{
    // for a full complex array of complex type the dimensions of the array are {N[0], N[1], NpadFT[2]}
    // for a half complex array of complex type the dimensions of the array are {N[0], N[1], NpadFT[2] / 2 + 1}
    double complex *afermi    = (double complex *) in1,
                   *bfermi    = (double complex *) in2,
                   *cfermi    = (double complex *) in3,
                   *afermibar = (double complex *) in4,
                   *bfermibar = (double complex *) in5,
                   *cfermibar = (double complex *) in6,
                   *product   = (double complex *) out;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) { // (NpadFT[2] / 2 + 1)
                n = k + N[2] * (j + N[1] * i);
                product[n] = afermi[n] * bfermi[n] * conj(cfermibar[n])
                           + afermibar[n] * bfermibar[n] * conj(cfermi[n]);
            }
        }
    }
}


// convolution products afermi * bfermi conj(cfermibar) + afermibar * bfermibar * conj(cfermi)
void product2(void *in1,
              void *in2,
              void *in3,
              void *in4,
              void *in5,
              void *in6,
              void *in7,
              void *in8,
              void *in9,
              void *in10,
              void *in11,
              void *in12,
              void *out1,
              void *out2,
              void *out3,
              void *out4,
              size_t *N)
{
    double complex *a1fermi    = (double complex *) in1,
                   *a2fermi    = (double complex *) in2,
                   *a3fermi    = (double complex *) in3,
                   *a4fermi    = (double complex *) in4,
                   *bfermi     = (double complex *) in5,
                   *cfermi     = (double complex *) in6,
                   *a1fermibar = (double complex *) in7,
                   *a2fermibar = (double complex *) in8,
                   *a3fermibar = (double complex *) in9,
                   *a4fermibar = (double complex *) in10,
                   *bfermibar  = (double complex *) in11,
                   *cfermibar  = (double complex *) in12,
                   *prod1      = (double complex *) out1,
                   *prod2      = (double complex *) out2,
                   *prod3      = (double complex *) out3,
                   *prod4      = (double complex *) out4;

    double complex tmp1, tmp2;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp1 = bfermi[n]    * conj(cfermibar[n]);
                tmp2 = bfermibar[n] * conj(cfermi[n]);
                prod1[n] = a1fermi[n] * tmp1 + a1fermibar[n] * tmp2;
                prod2[n] = a2fermi[n] * tmp1 + a2fermibar[n] * tmp2;
                prod3[n] = a3fermi[n] * tmp1 + a3fermibar[n] * tmp2;
                prod4[n] = a4fermi[n] * tmp1 + a4fermibar[n] * tmp2;
            }
        }
    }

}


// convolution products afermi * bfermi * conj(cfermibar) + afermibar * bfermibar * conj(cfermi)
void product3(void *in1,
              void *in2,
              void *in3,
              void *in4,
              void *in5,
              void *in6,
              void *in7,
              void *in8,
              void *in9,
              void *in10,
              void *out1,
              void *out2,
              void *out3,
              void *out4,
              void *out5,
              void *out6,
              void *out7,
              void *out8,
              void *out9,
              void *out10,
              size_t *N)
{

              double complex *a1fermi    = (double complex *) in1,
                             *a2fermi    = (double complex *) in2,
                             *a3fermi    = (double complex *) in3,
                             *a4fermi    = (double complex *) in4,
                             *cfermi     = (double complex *) in5,
                             *a1fermibar = (double complex *) in6,
                             *a2fermibar = (double complex *) in7,
                             *a3fermibar = (double complex *) in8,
                             *a4fermibar = (double complex *) in9,
                             *cfermibar  = (double complex *) in10,
                             *prod1      = (double complex *) out1,
                             *prod2      = (double complex *) out2,
                             *prod3      = (double complex *) out3,
                             *prod4      = (double complex *) out4,
                             *prod5      = (double complex *) out5,
                             *prod6      = (double complex *) out6,
                             *prod7      = (double complex *) out7,
                             *prod8      = (double complex *) out8,
                             *prod9      = (double complex *) out9,
                             *prod10     = (double complex *) out10;


    size_t i, j, k, n;
    double complex tmp1, tmp2;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp1 = conj(cfermibar[n]);
                tmp2 = conj(cfermi[n]);
                prod1[n]  = a1fermi[n] * a1fermi[n] * tmp1 + a1fermibar[n] * a1fermibar[n] * tmp2;
                prod2[n]  = a1fermi[n] * a2fermi[n] * tmp1 + a1fermibar[n] * a2fermibar[n] * tmp2;
                prod3[n]  = a1fermi[n] * a3fermi[n] * tmp1 + a1fermibar[n] * a3fermibar[n] * tmp2;
                prod4[n]  = a1fermi[n] * a4fermi[n] * tmp1 + a1fermibar[n] * a4fermibar[n] * tmp2;
                prod5[n]  = a2fermi[n] * a2fermi[n] * tmp1 + a2fermibar[n] * a2fermibar[n] * tmp2;
                prod6[n]  = a2fermi[n] * a3fermi[n] * tmp1 + a2fermibar[n] * a3fermibar[n] * tmp2;
                prod7[n]  = a2fermi[n] * a4fermi[n] * tmp1 + a2fermibar[n] * a4fermibar[n] * tmp2;
                prod8[n]  = a3fermi[n] * a3fermi[n] * tmp1 + a3fermibar[n] * a3fermibar[n] * tmp2;
                prod9[n]  = a3fermi[n] * a4fermi[n] * tmp1 + a3fermibar[n] * a4fermibar[n] * tmp2;
                prod10[n] = a4fermi[n] * a4fermi[n] * tmp1 + a4fermibar[n] * a4fermibar[n] * tmp2;
            }
        }
    }


}

// convolution products afermi * bfermi conj(cfermibar) + afermibar * bfermibar * conj(cfermi)
void product4(void *in1,
              void *in2,
              void *in3,
              void *in4,
              void *in5,
              void *in6,
              void *in7,
              void *in8,
              void *in9,
              void *in10,
              void *in11,
              void *in12,
              void *in13,
              void *in14,
              void *in15,
              void *in16,
              void *in17,
              void *in18,
              void *in19,
              void *in20,
              void *in21,
              void *in22,
              void *in23,
              void *in24,
              void *out1,
              void *out2,
              void *out3,
              void *out4,
              void *out5,
              void *out6,
              void *out7,
              void *out8,
              void *out9,
              void *out10,
              size_t *N)
{
              double complex *a1fermi     = (double complex *) in1,
                             *a2fermi     = (double complex *) in2,
                             *a3fermi     = (double complex *) in3,
                             *a4fermi     = (double complex *) in4,
                             *a5fermi     = (double complex *) in5,
                             *a6fermi     = (double complex *) in6,
                             *a7fermi     = (double complex *) in7,
                             *a8fermi     = (double complex *) in8,
                             *a9fermi     = (double complex *) in9,
                             *a10fermi    = (double complex *) in10,
                             *bfermi      = (double complex *) in11,
                             *cfermi      = (double complex *) in12,
                             *a1fermibar  = (double complex *) in13,
                             *a2fermibar  = (double complex *) in14,
                             *a3fermibar  = (double complex *) in15,
                             *a4fermibar  = (double complex *) in16,
                             *a5fermibar  = (double complex *) in17,
                             *a6fermibar  = (double complex *) in18,
                             *a7fermibar  = (double complex *) in19,
                             *a8fermibar  = (double complex *) in20,
                             *a9fermibar  = (double complex *) in21,
                             *a10fermibar = (double complex *) in22,
                             *bfermibar   = (double complex *) in23,
                             *cfermibar   = (double complex *) in24,
                             *prod1       = (double complex *) out1,
                             *prod2       = (double complex *) out2,
                             *prod3       = (double complex *) out3,
                             *prod4       = (double complex *) out4,
                             *prod5       = (double complex *) out5,
                             *prod6       = (double complex *) out6,
                             *prod7       = (double complex *) out7,
                             *prod8       = (double complex *) out8,
                             *prod9       = (double complex *) out9,
                             *prod10      = (double complex *) out10;


    size_t i, j, k, n;
    double complex tmp1, tmp2;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp1 = bfermi[n]    * conj(cfermibar[n]);
                tmp2 = bfermibar[n] * conj(cfermi[n]);
                prod1[n]  = a1fermi[n]  * tmp1 + a1fermibar[n]  * tmp2;
                prod2[n]  = a2fermi[n]  * tmp1 + a2fermibar[n]  * tmp2;
                prod3[n]  = a3fermi[n]  * tmp1 + a3fermibar[n]  * tmp2;
                prod4[n]  = a4fermi[n]  * tmp1 + a4fermibar[n]  * tmp2;
                prod5[n]  = a5fermi[n]  * tmp1 + a5fermibar[n]  * tmp2;
                prod6[n]  = a6fermi[n]  * tmp1 + a6fermibar[n]  * tmp2;
                prod7[n]  = a7fermi[n]  * tmp1 + a7fermibar[n]  * tmp2;
                prod8[n]  = a8fermi[n]  * tmp1 + a8fermibar[n]  * tmp2;
                prod9[n]  = a9fermi[n]  * tmp1 + a9fermibar[n]  * tmp2;
                prod10[n] = a10fermi[n] * tmp1 + a10fermibar[n] * tmp2;
            }
        }
    }

}


struct SelfEnergyParams
{
    size_t *N,
           *padFT,
           *NpadFT;
    double lambda,
           J,
           DELTAomega;
    double *Fermi_array,
           *Fermibar_array,
           *en,
           *c,
           *s,
           *rho_g_temp,
           *rho_psi0,
           *rho_psi1,
           *rho_chi0,
           *rho_chi1,
           *rho_chi2;
};


MKL_INT self_energy(void *params)
{

    // ON = 1, OFF = 0 defined in C_extention.h
    int PRINT = OFF,
        CHECKSUM = OFF;
    if (PRINT == ON) {
        printf("PRINT DATA ON\n");
    }
    if (CHECKSUM == ON) {
        printf("CHECKSUM ON\n");
    }

    struct SelfEnergyParams *p = (struct SelfEnergyParams *) params;
    size_t *N = p->N,
           *padFT = p->padFT,
           *NpadFT = p->NpadFT;
    double lambda = p->lambda,
           J = p->J,
           DELTAomega = p->DELTAomega;
    double *fermi = p->Fermi_array,
           *fermibar = p->Fermibar_array,
           *en = p->en,
           *c = p->c,
           *s = p->s,
           *rho_g_temp = p->rho_g_temp,
           *rho_psi0 = p->rho_psi0,
           *rho_psi1 = p->rho_psi1,
           *rho_chi0 = p->rho_chi0,
           *rho_chi1 = p->rho_chi1,
           *rho_chi2 = p->rho_chi2;

    size_t i, j, k, l, m, n;


    size_t Nhalfcmplx[3] = {N[0], N[1], (NpadFT[2] / 2 + 1)},
           Nhalfreal[3]  = {N[0], N[1], 2 * (NpadFT[2] / 2 + 1)};
    size_t Ns           = N[0] * N[1],
           N3           = N[0] * N[1] * N[2],
           N3halfcmplx  = N[0] * N[1] * Nhalfcmplx[2],
           N3padFT      = N[0] * N[1] * NpadFT[2],
           N3halfreal   = N[0] * N[1] * Nhalfreal[2];

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

    MKL_LONG halfcmplx_strides[4]  = {0, Nhalfcmplx[2] * Nhalfcmplx[1], Nhalfcmplx[2],  1}, // strides for half-complex padded-array
             full_strides[4]  = {0, NpadFT[2] * NpadFT[1], NpadFT[2],   1}, // strides for full padded-array
             halfreal_strides[4] = {0, Nhalfreal[2] * Nhalfreal[1], Nhalfreal[2], 1}; // strides for half-real padded-array


    status = DftiCreateDescriptor(&hand_forward, DFTI_DOUBLE, DFTI_REAL, 3, NpadFT);
    status = DftiSetValue(hand_forward, DFTI_PLACEMENT, DFTI_INPLACE); // Defalut 
    status = DftiSetValue(hand_forward, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(hand_forward, DFTI_INPUT_STRIDES, halfreal_strides);
    status = DftiSetValue(hand_forward, DFTI_OUTPUT_STRIDES, halfcmplx_strides);
    status = DftiCommitDescriptor(hand_forward);
 
    status = DftiCreateDescriptor(&hand_backward, DFTI_DOUBLE, DFTI_REAL, 3, NpadFT);
    status = DftiSetValue(hand_backward, DFTI_PLACEMENT, DFTI_INPLACE); // Defalut 
    status = DftiSetValue(hand_backward, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    status = DftiSetValue(hand_backward, DFTI_INPUT_STRIDES, halfcmplx_strides);
    status = DftiSetValue(hand_backward, DFTI_OUTPUT_STRIDES, halfreal_strides);
    status = DftiCommitDescriptor(hand_backward);



    double *out1  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out2  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out3  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out4  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out5  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out6  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out7  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out8  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out9  = mkl_calloc(N3halfreal, sizeof(double), 64),
           *out10 = mkl_calloc(N3halfreal, sizeof(double), 64);


    // g(p)g(q)g(p+q-k) -> ggg
    double *temp1  = mkl_calloc(N3, sizeof(double), 64), // -e_p ggg
           *temp2  = mkl_calloc(N3, sizeof(double), 64), // -J_{p-k}/2 ggg
           *temp3  = mkl_calloc(N3, sizeof(double), 64), // -ggg
           *temp4  = mkl_calloc(N3, sizeof(double), 64), // -e_p e_{p+q-k} ggg
           *temp5  = mkl_calloc(N3, sizeof(double), 64), // -e_p J_{k-q}/2 ggg
           *temp6  = mkl_calloc(N3, sizeof(double), 64), // -e_p J_{k-p}/2 ggg
           *temp7  = mkl_calloc(N3, sizeof(double), 64), // -e_{p+q-k}J_{k-p}ggg
           *temp8  = mkl_calloc(N3, sizeof(double), 64), // -J_{k-p}J_{p-k}/4 ggg
           *temp9  = mkl_calloc(N3, sizeof(double), 64), // -J_{k-p}J_{k-q}/4 ggg
           *temp10 = mkl_calloc(N3, sizeof(double), 64); // -e_{p+q-k} ggg 


    double *rho_g_fermi       = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_rho_g_fermi    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cx_rho_g_fermi    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx_rho_g_fermi    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cy_rho_g_fermi    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sy_rho_g_fermi    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_cx_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_sx_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_cy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_sy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64);


    double *rho_g_fermibar       = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_rho_g_fermibar    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cx_rho_g_fermibar    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx_rho_g_fermibar    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cy_rho_g_fermibar    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sy_rho_g_fermibar    = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_cx_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_sx_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_cy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *en_sy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64);

        double rho_g_fermi_tmp,
               rho_g_fermibar_tmp;
        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m  = (j + N[1] * i);
                for (k = 0; k < N[2]; k++) {
                    l = k + padFT[2] + Nhalfreal[2] * m;
                    n = k + N[2] * m;

                    rho_g_fermi_tmp      = rho_g_temp[n] * fermi[k];
                    rho_g_fermi[l]       = rho_g_fermi_tmp;
                    en_rho_g_fermi[l]    = en[m] * rho_g_fermi_tmp;
                    cx_rho_g_fermi[l]    = c[i]  * rho_g_fermi_tmp;
                    sx_rho_g_fermi[l]    = s[i]  * rho_g_fermi_tmp;
                    cy_rho_g_fermi[l]    = c[j]  * rho_g_fermi_tmp;
                    sy_rho_g_fermi[l]    = s[j]  * rho_g_fermi_tmp;
                    en_cx_rho_g_fermi[l] = en[m] * c[i] * rho_g_fermi_tmp;
                    en_sx_rho_g_fermi[l] = en[m] * s[i] * rho_g_fermi_tmp;
                    en_cy_rho_g_fermi[l] = en[m] * c[j] * rho_g_fermi_tmp;
                    en_sy_rho_g_fermi[l] = en[m] * s[j] * rho_g_fermi_tmp;

                    rho_g_fermibar_tmp      = rho_g_temp[n] * fermibar[k];
                    rho_g_fermibar[l]       = rho_g_fermibar_tmp;
                    en_rho_g_fermibar[l]    = en[m] * rho_g_fermibar_tmp;
                    cx_rho_g_fermibar[l]    = c[i]  * rho_g_fermibar_tmp;
                    sx_rho_g_fermibar[l]    = s[i]  * rho_g_fermibar_tmp;
                    cy_rho_g_fermibar[l]    = c[j]  * rho_g_fermibar_tmp;
                    sy_rho_g_fermibar[l]    = s[j]  * rho_g_fermibar_tmp;
                    en_cx_rho_g_fermibar[l] = en[m] * c[i] * rho_g_fermibar_tmp;
                    en_sx_rho_g_fermibar[l] = en[m] * s[i] * rho_g_fermibar_tmp;
                    en_cy_rho_g_fermibar[l] = en[m] * c[j] * rho_g_fermibar_tmp;
                    en_sy_rho_g_fermibar[l] = en[m] * s[j] * rho_g_fermibar_tmp;
            }
        }
    }

    if (CHECKSUM) {
        check_sum_half_real_3d("rho_g_fermi", rho_g_fermi, NpadFT);
        check_sum_half_real_3d("en_rho_g_fermi", en_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("cx_rho_g_fermi", cx_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("sx_rho_g_fermi", sx_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("cy_rho_g_fermi", cy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("sy_rho_g_fermi", sy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("en_cx_rho_g_fermi", en_cx_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("en_sx_rho_g_fermi", en_sx_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("en_cy_rho_g_fermi", en_cy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("en_sx_rho_g_fermi", en_sx_rho_g_fermi, NpadFT);
    }

    if (CHECKSUM) {
        check_sum_half_real_3d("rho_g_fermibar", rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("en_rho_g_fermibar", en_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("cx_rho_g_fermibar", cx_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("sx_rho_g_fermibar", sx_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("cy_rho_g_fermibar", cy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("sy_rho_g_fermibar", sy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("en_cx_rho_g_fermibar", en_cx_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("en_sx_rho_g_fermibar", en_sx_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("en_cy_rho_g_fermibar", en_cy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("en_sx_rho_g_fermibar", en_sx_rho_g_fermibar, NpadFT);
    }

    status = DftiComputeForward(hand_forward, rho_g_fermi);
    status = DftiComputeForward(hand_forward, en_rho_g_fermi);
    status = DftiComputeForward(hand_forward, cx_rho_g_fermi);
    status = DftiComputeForward(hand_forward, sx_rho_g_fermi);
    status = DftiComputeForward(hand_forward, cy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, sy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, en_cx_rho_g_fermi);
    status = DftiComputeForward(hand_forward, en_sx_rho_g_fermi);
    status = DftiComputeForward(hand_forward, en_cy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, en_sy_rho_g_fermi);


    if (CHECKSUM) {
        check_sum_half_complex_3d("rho_g_fermi", rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("en_rho_g_fermi", en_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("cx_rho_g_fermi", cx_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("sx_rho_g_fermi", sx_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("cy_rho_g_fermi", cy_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("sy_rho_g_fermi", sy_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("en_cx_rho_g_fermi", en_cx_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("en_sx_rho_g_fermi", en_sx_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("en_cy_rho_g_fermi", en_cy_rho_g_fermi, NpadFT);
        check_sum_half_complex_3d("en_sx_rho_g_fermi", en_sx_rho_g_fermi, NpadFT);
    }


    status = DftiComputeForward(hand_forward, rho_g_fermibar);
    status = DftiComputeForward(hand_forward, en_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, cx_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, sx_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, cy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, sy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, en_cx_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, en_sx_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, en_cy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, en_sy_rho_g_fermibar);

    if (CHECKSUM) {
        check_sum_half_complex_3d("rho_g_fermibar", rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("en_rho_g_fermibar", en_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("cx_rho_g_fermibar", cx_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("sx_rho_g_fermibar", sx_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("cy_rho_g_fermibar", cy_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("sy_rho_g_fermibar", sy_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("en_cx_rho_g_fermibar", en_cx_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("en_sx_rho_g_fermibar", en_sx_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("en_cy_rho_g_fermibar", en_cy_rho_g_fermibar, NpadFT);
        check_sum_half_complex_3d("en_sx_rho_g_fermibar", en_sx_rho_g_fermibar, NpadFT);
    }


// psi: {{{

    double scale1 = (pow(DELTAomega, 2)) / (pow(Ns, 2) * N3padFT);

    // psi_66b(k) = - 2 lambda (epsilon_p - u0/2 + 1/2 J_{k-p})
    // temp1 = -epsilon_p * ggg
    // out1 = afermi * bfermi * conj(cfermibar) + afermibar * bfermibar * conj(cfermi)
    product1(en_rho_g_fermi,    rho_g_fermi,    rho_g_fermi,
             en_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out1, Nhalfcmplx);
    status = DftiComputeBackward(hand_backward, out1);
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp1[n] = scale1 * out1[l];
            }
        }
    }


    /*
     * temp2   = -(1/2) J_{k-p} ggg
     * tmp1    = rho_g_fermi[n]    * conj(rho_g_fermibar[n]);
     * tmp2    = rho_g_fermibar[n] * conj(rho_g_fermi[n]);
     * out1[n] = cx_rho_g_fermi[n] * tmp1 + cx_rho_g_fermibar[n] * tmp2;
     * out2[n] = sx_rho_g_fermi[n] * tmp1 + sx_rho_g_fermibar[n] * tmp2;
     * out3[n] = cy_rho_g_fermi[n] * tmp1 + cy_rho_g_fermibar[n] * tmp2;
     * out4[n] = sy_rho_g_fermi[n] * tmp1 + sy_rho_g_fermibar[n] * tmp2;
     */
    product2(cx_rho_g_fermi, sx_rho_g_fermi, cy_rho_g_fermi,
             sy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cx_rho_g_fermibar, sx_rho_g_fermibar, cy_rho_g_fermibar,
             sy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar,
             out1, out2, out3, out4, Nhalfcmplx);

    status = DftiComputeBackward(hand_backward, out1);
    status = DftiComputeBackward(hand_backward, out2);
    status = DftiComputeBackward(hand_backward, out3);
    status = DftiComputeBackward(hand_backward, out4);

    double scale2 = (2 * J) * scale1 / 2;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp2[n] = scale2 * (c[i] * out1[l] + s[i] * out2[l] + c[j] * out3[l] + s[j] * out4[l]);
            }
        }
    }


    // temp3 = -ggg
    product1(rho_g_fermi,    rho_g_fermi,    rho_g_fermi,
             rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out1, Nhalfcmplx);
    status = DftiComputeBackward(hand_backward, out1);
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp3[n] = scale1 * out1[l];
            }
        }
    }

    if (CHECKSUM) {
        check_sum_double_array_1d("temp1", temp1, N3);
        check_sum_double_array_1d("temp2", temp2, N3);
        check_sum_double_array_1d("temp3", temp3, N3);
    }

    double rho_psi0_65a,
           rho_psi1_65a;
    /* 
     *
     *
     * rho_psi_65a = - 2(epsilon_p - u0 / 2 + 1/2 J_{k-p}) ggg
     *
     * rho_psi power series expansion in (-u0/2)
     * rho_psi = rho_psi0 + (-u0 / 2) * rho_psi1 + ...
     *
     */
    for (i = 0; i < N3; i++) {
        rho_psi0_65a = 2 * temp1[i] + 2 * temp2[i];
        rho_psi1_65a = 2 * temp3[i];
        rho_psi0[i] = lambda * lambda * rho_psi0_65a;
        rho_psi1[i] = lambda * lambda * rho_psi1_65a;
    }


    if (CHECKSUM) {
        check_sum_double_array_1d("rho_psi0", rho_psi0, N3);
        check_sum_double_array_1d("rho_psi1", rho_psi1, N3);
    }


// }}}


// chi: {{{


    // chi_66b(k) = -2 lambda (epsilon_p - u0/2 + 1/2 J_{k-p}) * (epsilon_{p + q - k} - u0/2 + 1/2 J_{k-q}) ggg
    // chi_66a(k) = -2 lambda (epsilon_{p + q - k} - u0/2 + 1/2 J_{k-p}) (1/2 J_{p-k}) ggg
    // chi(k) = -(2 epsilon_p - u0 + J_{k-p)) * (epsilon_{p + q - k} - u0/2 + 1/2(J_{k-q} + J_{k-p}) ggg


    // temp4 = -epsilon_p epsilon_{p+q-k} ggg
    product1(en_rho_g_fermi, rho_g_fermi, en_rho_g_fermi,
             en_rho_g_fermibar, rho_g_fermibar, en_rho_g_fermibar, out1, Nhalfcmplx);
    status = DftiComputeBackward(hand_backward, out1);
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp4[n] = scale1 * out1[l];
            }
        }
    }


    /*
     * temp5 = -epsilon_p J_{k-q} / 2 ggg
     * tmp1 = en_rho_g_fermi[n]    * conj(rho_g_fermibar[n]);
     * tmp2 = en_rho_g_fermibar[n] * conj(rho_g_fermi[n]);
     * out1[n] = cx_rho_g_fermi[n] * tmp1 + cx_rho_g_fermibar[n] * tmp2;
     * out2[n] = sx_rho_g_fermi[n] * tmp1 + sx_rho_g_fermibar[n] * tmp2;
     * out3[n] = cy_rho_g_fermi[n] * tmp1 + cy_rho_g_fermibar[n] * tmp2;
     * out4[n] = sy_rho_g_fermi[n] * tmp1 + sy_rho_g_fermibar[n] * tmp2;
     * note that product3 can do the job of product2
     */
    product2(cx_rho_g_fermi, sx_rho_g_fermi, cy_rho_g_fermi,
             sy_rho_g_fermi, en_rho_g_fermi, rho_g_fermi,
             cx_rho_g_fermibar, sx_rho_g_fermibar, cy_rho_g_fermibar,
             sy_rho_g_fermibar, en_rho_g_fermibar, rho_g_fermibar,
             out1, out2, out3, out4, Nhalfcmplx);

    status = DftiComputeBackward(hand_backward, out1);
    status = DftiComputeBackward(hand_backward, out2);
    status = DftiComputeBackward(hand_backward, out3);
    status = DftiComputeBackward(hand_backward, out4);

    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp5[n] = scale2 * (c[i] * out1[l] + s[i] * out2[l] + c[j] * out3[l] + s[j] * out4[l]);
            }
        }
    }

    /*
     * temp6 = -epsilon_p J_{k-p} / 2 ggg
     * tmp1 = rho_g_fermi[n]    * conj(rho_g_fermibar[n]);
     * tmp2 = rho_g_fermibar[n] * conj(rho_g_fermi[n]);
     * out1[n] = en_cx_rho_g_fermi[n] * tmp1 + en_cx_rho_g_fermibar[n] * tmp2;
     * out2[n] = en_sx_rho_g_fermi[n] * tmp1 + en_sx_rho_g_fermibar[n] * tmp2;
     * out3[n] = en_cy_rho_g_fermi[n] * tmp1 + en_cy_rho_g_fermibar[n] * tmp2;
     * out4[n] = en_sy_rho_g_fermi[n] * tmp1 + en_sy_rho_g_fermibar[n] * tmp2;
     */
    product2(en_cx_rho_g_fermi, en_sx_rho_g_fermi, en_cy_rho_g_fermi,
             en_sy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             en_cx_rho_g_fermibar, en_sx_rho_g_fermibar, en_cy_rho_g_fermibar,
             en_sy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar,
             out1, out2, out3, out4, Nhalfcmplx);
    status = DftiComputeBackward(hand_backward, out1);
    status = DftiComputeBackward(hand_backward, out2);
    status = DftiComputeBackward(hand_backward, out3);
    status = DftiComputeBackward(hand_backward, out4);
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp6[n] = scale2 * (c[i] * out1[l] + s[i] * out2[l] + c[j] * out3[l] + s[j] * out4[l]);
            }
        }
    }

    /*
     * temp7 = -1/2 J_{k-p} epsilon_{p + q - k} ggg
     * tmp1 = en_rho_g_fermi[n]    * conj(rho_g_fermibar[n]);
     * tmp2 = en_rho_g_fermibar[n] * conj(rho_g_fermi[n]);
     * out1[n] = cx_rho_g_fermi[n] * tmp1 + cx_rho_g_fermibar[n] * tmp2;
     * out2[n] = sx_rho_g_fermi[n] * tmp1 + sx_rho_g_fermibar[n] * tmp2;
     * out3[n] = cy_rho_g_fermi[n] * tmp1 + cy_rho_g_fermibar[n] * tmp2;
     * out4[n] = sy_rho_g_fermi[n] * tmp1 + sy_rho_g_fermibar[n] * tmp2;
     */
    product2(cx_rho_g_fermi, sx_rho_g_fermi, cy_rho_g_fermi,
             sy_rho_g_fermi, rho_g_fermi, en_rho_g_fermi,
             cx_rho_g_fermibar, sx_rho_g_fermibar, cy_rho_g_fermibar,
             sy_rho_g_fermibar, rho_g_fermibar, en_rho_g_fermibar,
             out1, out2, out3, out4, Nhalfcmplx);
    status = DftiComputeBackward(hand_backward, out1);
    status = DftiComputeBackward(hand_backward, out2);
    status = DftiComputeBackward(hand_backward, out3);
    status = DftiComputeBackward(hand_backward, out4);
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp7[n] = scale2 * (c[i] * out1[l] + s[i] * out2[l] + c[j] * out3[l] + s[j] * out4[l]);
            }
        }
    }

    /*
     * temp8 = -J_{k-p} J_{q-k} / 4 ggg
     * a1 = cx_rho_g, a2 = sx_rho_g, a3 = cy_rho_g, a4 = sy_rho_g, c = rho_g
     * tmp1 = conj(cfermibar[n]);
     * tmp2 = conj(cfermi[n]);
     * out1[n]  = a1fermi[n] * a1fermi[n] * tmp1 + a1fermibar[n] * a1fermibar[n] * tmp2;
     * out2[n]  = a1fermi[n] * a2fermi[n] * tmp1 + a1fermibar[n] * a2fermibar[n] * tmp2;
     * out3[n]  = a1fermi[n] * a3fermi[n] * tmp1 + a1fermibar[n] * a3fermibar[n] * tmp2;
     * out4[n]  = a1fermi[n] * a4fermi[n] * tmp1 + a1fermibar[n] * a4fermibar[n] * tmp2;
     * out5[n]  = a2fermi[n] * a2fermi[n] * tmp1 + a2fermibar[n] * a2fermibar[n] * tmp2;
     * out6[n]  = a2fermi[n] * a3fermi[n] * tmp1 + a2fermibar[n] * a3fermibar[n] * tmp2;
     * out7[n]  = a2fermi[n] * a4fermi[n] * tmp1 + a2fermibar[n] * a4fermibar[n] * tmp2;
     * out8[n]  = a3fermi[n] * a3fermi[n] * tmp1 + a3fermibar[n] * a3fermibar[n] * tmp2;
     * out9[n]  = a3fermi[n] * a4fermi[n] * tmp1 + a3fermibar[n] * a4fermibar[n] * tmp2;
     * out10[n] = a4fermi[n] * a4fermi[n] * tmp1 + a4fermibar[n] * a4fermibar[n] * tmp2;
     */
    product3(cx_rho_g_fermi,    sx_rho_g_fermi,    cy_rho_g_fermi,    sy_rho_g_fermi,    rho_g_fermi,
             cx_rho_g_fermibar, sx_rho_g_fermibar, cy_rho_g_fermibar, sy_rho_g_fermibar, rho_g_fermibar,
             out1, out2, out3, out4, out5, out6, out7, out8, out9, out10, Nhalfcmplx);

    status = DftiComputeBackward(hand_backward, out1);
    status = DftiComputeBackward(hand_backward, out2);
    status = DftiComputeBackward(hand_backward, out3);
    status = DftiComputeBackward(hand_backward, out4);
    status = DftiComputeBackward(hand_backward, out5);
    status = DftiComputeBackward(hand_backward, out6);
    status = DftiComputeBackward(hand_backward, out7);
    status = DftiComputeBackward(hand_backward, out8);
    status = DftiComputeBackward(hand_backward, out9);
    status = DftiComputeBackward(hand_backward, out10);

    double scale3 = pow(2*J, 2) * scale1 / 4;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp8[n] = scale3
                 * (c[i]*c[i] * out1[l]  + 2*c[i]*s[i] * out2[l] + 2*c[i]*c[j] * out3[l]  + 2*c[i]*s[j] * out4[l]
                                         +   s[i]*s[i] * out5[l] + 2*s[i]*c[j] * out6[l]  + 2*s[i]*s[j] * out7[l]
                                                                 +   c[j]*c[j] * out8[l]  + 2*c[j]*s[j] * out9[l]
                                                                                          +   s[j]*s[j] * out10[l]);
            }
        }
    }


    // temp10 = epsilon_{p + q -k} ggg

    product1(rho_g_fermi,    rho_g_fermi,    en_rho_g_fermi,
             rho_g_fermibar, rho_g_fermibar, en_rho_g_fermibar, out1, Nhalfcmplx);
    status = DftiComputeBackward(hand_backward, out1);

    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp10[n] = scale1 * out1[l];
            }
        }
    }


    mkl_free(en_rho_g_fermi);
    mkl_free(cx_rho_g_fermi);
    mkl_free(sx_rho_g_fermi);
    mkl_free(cy_rho_g_fermi);
    mkl_free(sy_rho_g_fermi);
    mkl_free(en_cx_rho_g_fermi);
    mkl_free(en_sx_rho_g_fermi);
    mkl_free(en_cy_rho_g_fermi);
    mkl_free(en_sy_rho_g_fermi);

    mkl_free(en_rho_g_fermibar);
    mkl_free(cx_rho_g_fermibar);
    mkl_free(sx_rho_g_fermibar);
    mkl_free(cy_rho_g_fermibar);
    mkl_free(sy_rho_g_fermibar);
    mkl_free(en_cx_rho_g_fermibar);
    mkl_free(en_sx_rho_g_fermibar);
    mkl_free(en_cy_rho_g_fermibar);
    mkl_free(en_sy_rho_g_fermibar);


    double *cx2_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx2_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cy2_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sy2_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64);

    double *cx2_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx2_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cy2_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sy2_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64);


    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m  = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;

                rho_g_fermi_tmp      = rho_g_temp[n] * fermi[k];
                cx2_rho_g_fermi[l]   = c[i] * c[i] * rho_g_fermi_tmp;
                sx2_rho_g_fermi[l]   = s[i] * s[i] * rho_g_fermi_tmp;
                cy2_rho_g_fermi[l]   = c[j] * c[j] * rho_g_fermi_tmp;
                sy2_rho_g_fermi[l]   = s[j] * s[j] * rho_g_fermi_tmp;

                rho_g_fermibar_tmp      = rho_g_temp[n] * fermibar[k];
                cx2_rho_g_fermibar[l]   = c[i] * c[i] * rho_g_fermibar_tmp;
                sx2_rho_g_fermibar[l]   = s[i] * s[i] * rho_g_fermibar_tmp; // error found s[i] * c[i] --> s[i] * s[i]
                cy2_rho_g_fermibar[l]   = c[j] * c[j] * rho_g_fermibar_tmp;
                sy2_rho_g_fermibar[l]   = s[j] * s[j] * rho_g_fermibar_tmp;
            }
        }
    }

    if (CHECKSUM) {
        check_sum_half_real_3d("cx2_rho_g_fermi", cx2_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("sx2_rho_g_fermi", sx2_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("cy2_rho_g_fermi", cy2_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("sy2_rho_g_fermi", sy2_rho_g_fermi, NpadFT);

        check_sum_half_real_3d("cx2_rho_g_fermibar", cx2_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("sx2_rho_g_fermibar", sx2_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("cy2_rho_g_fermibar", cy2_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("sy2_rho_g_fermibar", sy2_rho_g_fermibar, NpadFT);
    }


    status = DftiComputeForward(hand_forward, cx2_rho_g_fermi);
    status = DftiComputeForward(hand_forward, sx2_rho_g_fermi);
    status = DftiComputeForward(hand_forward, cy2_rho_g_fermi);
    status = DftiComputeForward(hand_forward, sy2_rho_g_fermi);

    status = DftiComputeForward(hand_forward, cx2_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, sx2_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, cy2_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, sy2_rho_g_fermibar);

    /*
     * temp9 =  J_{k-p} J_{p-k} / 4
     * a1 = cx2_rho_g,   a2 = sx2_rho_g,    a3 = cy2_rho_g,   a4 = sy2_rho_g,
     * a5 = cx_sx_rho_g, a6 = cx_cy_rho_g,  a7 = cx_sy_rho_g, a8 = sx_cy_rho_g,
     * a9 = sx_sy_rho_g, a10 = cy_sy_rho_g, b = rho_g, c = rho_g
     * tmp1 = bfermi[n]    * conj(cfermibar[n]);
     * tmp2 = bfermibar[n] * conj(cfermi[n]);
     * out1[n]  = a1fermi[n]  * tmp1 + a1fermibar[n]  * tmp2;
     * out2[n]  = a2fermi[n]  * tmp1 + a2fermibar[n]  * tmp2;
     * out3[n]  = a3fermi[n]  * tmp1 + a3fermibar[n]  * tmp2;
     * out4[n]  = a4fermi[n]  * tmp1 + a4fermibar[n]  * tmp2;
     * out5[n]  = a5fermi[n]  * tmp1 + a5fermibar[n]  * tmp2;
     * out6[n]  = a6fermi[n]  * tmp1 + a6fermibar[n]  * tmp2;
     * out7[n]  = a7fermi[n]  * tmp1 + a7fermibar[n]  * tmp2;
     * out8[n]  = a8fermi[n]  * tmp1 + a8fermibar[n]  * tmp2;
     * out9[n]  = a9fermi[n]  * tmp1 + a9fermibar[n]  * tmp2;
     * out10[n] = a10fermi[n] * tmp1 + a10fermibar[n] * tmp2;
     */

    product1(cx2_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cx2_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out1, Nhalfcmplx);
    product1(sx2_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             sx2_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out2, Nhalfcmplx);
    product1(cy2_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cy2_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out3, Nhalfcmplx);
    product1(sy2_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             sy2_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out4, Nhalfcmplx);


    mkl_free(cx2_rho_g_fermi);
    mkl_free(sx2_rho_g_fermi);
    mkl_free(cy2_rho_g_fermi);
    mkl_free(sy2_rho_g_fermi);

    mkl_free(cx2_rho_g_fermibar);
    mkl_free(sx2_rho_g_fermibar);
    mkl_free(cy2_rho_g_fermibar);
    mkl_free(sy2_rho_g_fermibar);


    double *cx_sx_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cx_cy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cx_sy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx_cy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx_sy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cy_sy_rho_g_fermi = mkl_calloc(N3halfreal, sizeof(double), 64);

    double *cx_sx_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cx_cy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cx_sy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx_cy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *sx_sy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64),
           *cy_sy_rho_g_fermibar = mkl_calloc(N3halfreal, sizeof(double), 64);


    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m  = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;

                rho_g_fermi_tmp      = rho_g_temp[n] * fermi[k];
                cx_sx_rho_g_fermi[l] = c[i] * s[i] * rho_g_fermi_tmp;
                cx_cy_rho_g_fermi[l] = c[i] * c[j] * rho_g_fermi_tmp;
                cx_sy_rho_g_fermi[l] = c[i] * s[j] * rho_g_fermi_tmp;
                sx_cy_rho_g_fermi[l] = s[i] * c[j] * rho_g_fermi_tmp;
                sx_sy_rho_g_fermi[l] = s[i] * s[j] * rho_g_fermi_tmp;
                cy_sy_rho_g_fermi[l] = c[j] * s[j] * rho_g_fermi_tmp;

                rho_g_fermibar_tmp      = rho_g_temp[n] * fermibar[k];
                cx_sx_rho_g_fermibar[l] = c[i] * s[i] * rho_g_fermibar_tmp;
                cx_cy_rho_g_fermibar[l] = c[i] * c[j] * rho_g_fermibar_tmp;
                cx_sy_rho_g_fermibar[l] = c[i] * s[j] * rho_g_fermibar_tmp;
                sx_cy_rho_g_fermibar[l] = s[i] * c[j] * rho_g_fermibar_tmp;
                sx_sy_rho_g_fermibar[l] = s[i] * s[j] * rho_g_fermibar_tmp;
                cy_sy_rho_g_fermibar[l] = c[j] * s[j] * rho_g_fermibar_tmp;
            }
        }
    }

    if (CHECKSUM) {
        check_sum_half_real_3d("cx_sx_rho_g_fermi", cx_sx_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("cx_cy_rho_g_fermi", cx_cy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("cx_sy_rho_g_fermi", cx_sy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("sx_cy_rho_g_fermi", sx_cy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("sx_sy_rho_g_fermi", sx_sy_rho_g_fermi, NpadFT);
        check_sum_half_real_3d("cy_sy_rho_g_fermi", cy_sy_rho_g_fermi, NpadFT);

        check_sum_half_real_3d("cx_sx_rho_g_fermibar", cx_sx_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("cx_cy_rho_g_fermibar", cx_cy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("cx_sy_rho_g_fermibar", cx_sy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("sx_cy_rho_g_fermibar", sx_cy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("sx_sy_rho_g_fermibar", sx_sy_rho_g_fermibar, NpadFT);
        check_sum_half_real_3d("cy_sy_rho_g_fermibar", cy_sy_rho_g_fermibar, NpadFT);
    }


    status = DftiComputeForward(hand_forward, cx_sx_rho_g_fermi);
    status = DftiComputeForward(hand_forward, cx_cy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, cx_sy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, sx_cy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, sx_sy_rho_g_fermi);
    status = DftiComputeForward(hand_forward, cy_sy_rho_g_fermi);

    status = DftiComputeForward(hand_forward, cx_sx_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, cx_cy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, cx_sy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, sx_cy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, sx_sy_rho_g_fermibar);
    status = DftiComputeForward(hand_forward, cy_sy_rho_g_fermibar);


    // temp9 continued
    product1(cx_sx_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cx_sx_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out5, Nhalfcmplx);
    product1(cx_cy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cx_cy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out6, Nhalfcmplx);
    product1(cx_sy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cx_sy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out7, Nhalfcmplx);
    product1(sx_cy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             sx_cy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out8, Nhalfcmplx);
    product1(sx_sy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             sx_sy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out9, Nhalfcmplx);
    product1(cy_sy_rho_g_fermi, rho_g_fermi, rho_g_fermi,
             cy_sy_rho_g_fermibar, rho_g_fermibar, rho_g_fermibar, out10, Nhalfcmplx);


    mkl_free(rho_g_fermi);
    mkl_free(cx_sx_rho_g_fermi);
    mkl_free(cx_cy_rho_g_fermi);
    mkl_free(cx_sy_rho_g_fermi);
    mkl_free(sx_cy_rho_g_fermi);
    mkl_free(sx_sy_rho_g_fermi);
    mkl_free(cy_sy_rho_g_fermi);

    mkl_free(rho_g_fermibar);
    mkl_free(cx_sx_rho_g_fermibar);
    mkl_free(cx_cy_rho_g_fermibar);
    mkl_free(cx_sy_rho_g_fermibar);
    mkl_free(sx_cy_rho_g_fermibar);
    mkl_free(sx_sy_rho_g_fermibar);
    mkl_free(cy_sy_rho_g_fermibar);


    status = DftiComputeBackward(hand_backward, out1);
    status = DftiComputeBackward(hand_backward, out2);
    status = DftiComputeBackward(hand_backward, out3);
    status = DftiComputeBackward(hand_backward, out4);
    status = DftiComputeBackward(hand_backward, out5);
    status = DftiComputeBackward(hand_backward, out6);
    status = DftiComputeBackward(hand_backward, out7);
    status = DftiComputeBackward(hand_backward, out8);
    status = DftiComputeBackward(hand_backward, out9);
    status = DftiComputeBackward(hand_backward, out10);


    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                l = k + padFT[2] + Nhalfreal[2] * m;
                n = k + N[2] * m;
                temp9[n] = scale3
                         * (c[i]*c[i] * out1[l] + s[i]*s[i] * out2[l] + c[j]*c[j] * out3[l] + s[j]*s[j] * out4[l]
                         + 2 * c[i]*s[i] * out5[l]  + 2 * c[i]*c[j] * out6[l] + 2 * c[i]*s[j] * out7[l]
                         + 2 * s[i]*c[j] * out8[l]  + 2 * s[i]*s[j] * out9[l] + 2 * c[j]*s[j] * out10[l]);
            }
        }
    }




    if (CHECKSUM) {
        check_sum_double_array_1d("temp4", temp4, N3);
        check_sum_double_array_1d("temp5", temp5, N3);
        check_sum_double_array_1d("temp6", temp6, N3);
        check_sum_double_array_1d("temp7", temp7, N3);
        check_sum_double_array_1d("temp8", temp8, N3);
        check_sum_double_array_1d("temp9", temp9, N3);
        check_sum_double_array_1d("temp10", temp10, N3);
    }


    // chi_66b(k) = -2(epsilon_p - u0/2 + 1/2 J_{k-p}) * (epsilon_{p + q - k} - u0/2 + 1/2 J_{k-q}) ggg
    // chi_66a(k) = -2(epsilon_{p + q - k} - u0/2 + 1/2 J_{k-p}) (1/2 J_{p-k}) ggg


    /*
     * rho_chi power series expansion in (-u0/2) 
     *
     * rho_chi = rho_chi0 + (-u0/2) rho_chi1 + (-u0/2)^2
     *
     */
    double rho_chi0_66b,
           rho_chi1_66b,
           rho_chi2_66b,
           rho_chi0_67a,
           rho_chi1_67a;

    for (i = 0; i < N3; i++) {
        rho_chi0_66b = 2 * (temp4[i] + temp5[i] + temp7[i] + temp8[i]);
        rho_chi1_66b = 2 * (temp1[i] + temp10[i] + 2 * temp2[i]);
        rho_chi2_66b = 2 * temp3[i];
        rho_chi0_67a = 2 * (temp6[i] + temp9[i]);
        rho_chi1_67a = 2 * temp2[i];
        rho_chi0[i] = pow(lambda, 2) * (rho_chi0_66b + rho_chi0_67a);
        rho_chi1[i] = pow(lambda, 2) * (rho_chi1_66b + rho_chi1_67a);
        rho_chi2[i] = pow(lambda, 2) * rho_chi2_66b;
    }


    if (CHECKSUM) {
        check_sum_double_array_1d("rho_chi0", rho_chi0, N3);
        check_sum_double_array_1d("rho_chi1", rho_chi1, N3);
        check_sum_double_array_1d("rho_chi2", rho_chi2, N3);
    }

// }}}


// phi: {{{
/*
{
    double *rho_phi0 = mkl_calloc(N3, sizeof(double), 64),
           *rho_phi1 = mkl_calloc(N3, sizeof(double), 64),
           *rho_phi2 = mkl_calloc(N3, sizeof(double), 64);

    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                rho_phi0[n] = rho_chi0[n] + en[m] * rho_psi0[n];
                rho_phi1[n] = (-1.0/2.0) * rho_chi1[n] + (-1.0/2.0) * rho_psi0[n] + en[m] * rho_psi1[n] *(-1.0/2.0); 
                rho_phi2[n] = (1/4.0)*rho_chi2[n] + (1.0/4.0) * rho_psi1[n];
            }
        }
    }


    if (CHECKSUM) {
        check_sum_double_array_3d("rho_phi0", rho_phi0, N);
        check_sum_double_array_3d("rho_phi1", rho_phi1, N);
        check_sum_double_array_3d("rho_phi2", rho_phi2, N);
    }

       mkl_free(rho_phi0);
       mkl_free(rho_phi1);
       mkl_free(rho_phi2);

 }
*/

// }}}


   goto cleanup;
    cleanup:
    DftiFreeDescriptor(&hand_forward);
    DftiFreeDescriptor(&hand_backward);

    mkl_free(out1);
    mkl_free(out2);
    mkl_free(out3);
    mkl_free(out4);
    mkl_free(out5);
    mkl_free(out6);
    mkl_free(out7);
    mkl_free(out8);
    mkl_free(out9);
    mkl_free(out10);

    mkl_free(temp1);
    mkl_free(temp2);
    mkl_free(temp3);
    mkl_free(temp4);
    mkl_free(temp5);
    mkl_free(temp6);
    mkl_free(temp7);
    mkl_free(temp8);
    mkl_free(temp9);
    mkl_free(temp10);

    return status;
}


#endif /* _SELFENERGYR2CDFTI_H */
