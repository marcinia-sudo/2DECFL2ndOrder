/* vim: set fdm=marker:
 *
 *
 *  MKL_EXTENSION.h
 *
 *
 */

#ifndef _MKL_EXTENSION_H
#define _MKL_EXTENSION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mkl.h> 



// pad_MKL_Complex16_array: {{{


void pad_MKL_Complex16_array_1d(MKL_Complex16 *in, MKL_Complex16 *out, size_t N, size_t pad)
{
    size_t Npad = N + 2 * pad;
    int n;
    for (n = 0; n < Npad; n++) {
        out[n].real = 0.0;
        out[n].imag = 0.0;
    }
    for (n = 0; n < N; n++) {
        out[n + pad].real = in[n].real;
        out[n + pad].imag = in[n].imag;
    }
}

void unpad_MKL_Complex16_array_1d(MKL_Complex16 *in, MKL_Complex16 *out, size_t N, size_t pad)
{
    int n;
    for (n = 0; n < N; n++) {
        out[n].real = in[n + pad].real;
        out[n].imag = in[n + pad].imag;
    }
}


void pad_MKL_Complex16_array_3d(MKL_Complex16 *in, MKL_Complex16 *out, size_t *N, size_t *pad)
{
    size_t rank = 3,
           N3pad = 1;
    size_t Npad[3];
    int i, j, k, n, m;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
        N3pad *= Npad[i];
    }
    for (i = 0; i < N3pad; i++) {
        out[i].real = 0;
        out[i].imag = 0;
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                m = k + N[2] * (j + N[1] * i);
                out[n].real = in[m].real;
                out[n].imag = in[n].imag;
            }
        }
    }
}

void unpad_MKL_Complex16_array_3d(MKL_Complex16 *in, MKL_Complex16 *out, size_t *N, size_t *pad)
{
    size_t rank = 3;
    size_t Npad[3];
    int i, j, k, n, m;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                m = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                out[n].real = in[m].real;
                out[n].imag = in[m].imag;
            }
        }
    }
}

// }}}


// pad_double_to_MKL_Complex16_array: {{{


void pad_double_to_MKL_Complex16_array_1d(double *in, MKL_Complex16 *out, size_t N, size_t pad)
{
    size_t Npad = N + 2 * pad;
    int n;
    for (n = 0; n < Npad; n++) {
        out[n].real = 0.0;
        out[n].imag = 0.0;
    }
    for (n = 0; n < N; n++) {
        out[n + pad].real = in[n];
        out[n + pad].imag = 0.0;
    }
}


void unpad_MKL_Complex16_to_double_array_1d(MKL_Complex16 *in, double *out, size_t N, size_t pad)
{
    int n;
    for (n = 0; n < N; n++) {
        out[n] = in[n + pad].real;
    }
}


void pad_double_to_MKL_Complex16_array_3d(double *in, MKL_Complex16 *out, size_t *N, size_t *pad)
{
    size_t rank = 3,
           N3pad = 1;
    size_t Npad[3];
    int i, j, k, n, m;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
        N3pad *= Npad[i];
    }
    for (i = 0; i < N3pad; i++) {
        out[i].real = 0.0;
        out[i].imag = 0.0;
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                m = k + N[2] * (j + N[1] * i);
                out[n].real = in[m];
                out[n].imag = 0.0;
            }
        }
    }
}


void unpad_MKL_Complex16_to_double_array_3d(MKL_Complex16 *in, double *out, size_t *N, size_t *pad)
{
    size_t rank = 3;
    size_t Npad[3];
    int i, j, k, n, m;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                m = k + N[2] * (j + N[1] * i); 
                out[m] = in[n].real;
            }
        }
    }
}


// }}}


// copy_MKL_Complex16_array: {{{


void copy_MKL_Complex16_array_1d(double complex *in, double complex *out, size_t N)
{
    int n;
    for (n = 0; n < N; n++) {
        in[n] = out[n];
    }
}


void copy_MKL_Compelx16_array_2d(double complex *in, double complex *out, size_t *N)
{
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            out[n] = in[n];
        }
    }
}


void copy_MKL_Complex16_array_3d(double complex *in, double complex *out, size_t *N)
{
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                out[n] = in[n];
            }
        }
    }
}


// }}} copy_double_complex_array


// print_MKL_Complex16_array: {{{


void print_MKL_Complex16_array_1d(MKL_Complex16 *in, size_t N)
{
    double re, im;
    int n;
    for (n = 0; n < N; n++) {
        re = in[n].real;
        im = in[n].imag;
        printf("%02d %+.8f%+.8fi \n", n, re, im);
    }
}


void print_MKL_Complex16_array_2d(MKL_Complex16 *in, size_t *N)
{
    double re, im;
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            re = in[n].real;
            im = in[n].imag;
            printf("%+.8f%+.8fi \n", re, im);
        }
        printf("\n");
    }
}


void print_MKL_Complex16_array_3d(MKL_Complex16 *in, size_t *N)
{
    double re, im;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                re = in[n].real;
                im = in[n].imag;
                printf("%02d %02d %02d %+.8f%+.8fi\n", i, j, k, re, im);
            }
            printf("\n\n");
        }
        printf("\n\n\n");
    }
}


// }}} print_double_complex_array


// write_to_file_MKL_Complex16_array: {{{


void write_to_file_MKL_Complex16_1d(char *name, MKL_Complex16 *in, MKL_LONG N)
{
    FILE *fdata = fopen(name, "w");
    double re, im;
    int n;
    for (n = 0; n < N; n++) {
        re = in[n].real;
        im = in[n].imag;
        fprintf(fdata, "%5.8lf%+5.8lf\n", re, im);
    }
    fclose(fdata);
}


void write_to_file_MKL_Complex16_array_2d(char *name, MKL_Complex16 *in, MKL_LONG *N)
{
    FILE *fdata = fopen(name, "w");
    double re, im;
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            re = in[n].real;
            im = in[n].imag;
            fprintf(fdata, "%5.8lf%+5.8lf\n", re, im);
        }
    }
    fclose(fdata);
}


void write_to_file_MKL_Complex16_array_3d(char *name, MKL_Complex16 *in, MKL_LONG *N)
{
    FILE *fdata = fopen(name, "w");
    double re, im;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                re = in[n].real;
                im = in[n].imag;
                fprintf(fdata, "%5.8lf%+5.8lf\n", re, im);
            }
        }
    }
    fclose(fdata);
}


// }}} write_to_file_double_array


#endif /* _MKL_EXTENSION_H */

