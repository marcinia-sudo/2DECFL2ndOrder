/* vim: set fdm=marker:
 *
 *
 *  C_EXTENSION.h
 *
 *
 */

#ifndef _C_EXTENSION_H
#define _C_EXTENSION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mkl.h> 


// we should create a library with the following header files:
//#include "Conversions.h" // Contains useful conversions
// etc
// and search for it in the makefile

#define Re(z) (__real__ z)
#define Im(z) (__imag__ z)
#define TRUE 1
#define FALSE 0
#define ON 1
#define OFF 0

// clock {{{

// Reports CPU_TIME and WALL_TIME and exit program if WALL_TIME exceeds TIME_LIMIT
void timer (char *name, clock_t tic, struct timeval tv1, int TIME_LIMIT) {

        clock_t toc = clock();
        double time_elapsed = (double) -(tic - toc) / CLOCKS_PER_SEC;
        printf("%s: CPU_TIME elapsed = %+5.15lf\n", name, time_elapsed);

        struct timeval tv2;
        gettimeofday(&tv2, NULL);
        double time_elapsed_usec = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000;
        double time_elapsed_sec = (double) (tv2.tv_sec - tv1.tv_sec);
        printf("%s: WALL_TIME = %f seconds, %f millsecs\n", name, time_elapsed_sec, time_elapsed_usec);

        if (time_elapsed_sec > TIME_LIMIT) {
            printf("%s: Time Limit\n", name);
            exit(0);
        }

}

// }}}


// power_fcts: {{{


int int_pow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
           result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

///    }}}


// sign_fcts: {{{

int intSign(int x)
{
    int res;
    if (x >= 0) {
        res = 1;
    } else if (x < 0) {
        res = -1;
    } else {
        printf("nan");
        // something ridiculous
        return 123456789;
    }
    return res;
}

// signum
int intSgn(int x)
{
    int res = 0;
    if (x > 0) {
        res = 1;
    } else if (x < 0) {
        res = -1;
    } else if (x == 0) {
        res = 0;
    } else {
        printf("nan\n");
        // something ridiculous
        return 123456789;
    }
    return res;
}


double Sign(double x)
{
    double res;
    if (x >= 0.0) {
        res = 1.0;
    } else if (x < 0.0 ) {
        res = -1.0;
    } else {
        printf("auxillary error sign fct nan\n");
        // something ridiculous
        return 123456789.0;
    }
    return res;
}


// signum
int Sgn(double x)
{
    int res = 0,
        eps = pow(10, -15);
    if (x > eps) {
        res = 1.0;
    } else if (x < -eps ) {
        res = -1.0;
    } else if (x >= -eps && x <= eps) {
        res = 0.0;
    } else {
        printf("auxillary error sgn fct nan\n");
        // something ridiculous
        return 123456789.0;
    }
    return res;
}

double SgnExact(double x)
{
    double res = 0;
//           eps = pow(10, -15);
    if (x > 0.0) {
        res = 1.0;
    } else if (x < 0.0 ) {
        res = -1.0;
    } else if (x == 0.0) {
        res = 0.0;
    } else {
        printf("nan\n");
        // something ridiculous
        return 123456789.0;
    }
    return res;
}


double THETASign(double x)
{
     return ((1.0 + Sign(x)) / 2.0);
}

double THETASgn(double x)
{
     return ((1.0 + Sgn(x)) / 2.0);
}
// }}}


// pad_array_double: {{{

void pad_double_array_1d(double *in, double *out, size_t N, size_t pad)
{
    size_t Npad = N + 2 * pad;
    int n, m;
    for (n = 0; n < N; n++) {
        m = n + (int) pad;
        out[m] = in[n];
    }
}

void unpad_double_array_1d(double *in, double *out, size_t N, size_t pad)
{
    int n;
    for (n = 0; n < N; n++) {
        out[n] = in[n + pad];
    }
}


void pad_double_array_3d(double *in, double *out, size_t *N, size_t *pad) {
// pad array does not initialize the out array
// malloc + memset can do this but calloc is more efficient
    size_t rank = 3,
           N3pad = 1;
    size_t Npad[3];
    int i, j, k, m, n;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                m = k + N[2] * (j + N[1] * i); 
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                out[n] = in[m];
            }
        }
    }

}


void unpad_double_array_3d(double *in, double *out, size_t *N, size_t *pad) {

    size_t rank = 3;
    size_t Npad[3];
    int i, j, k, m, n;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                m = k + N[2] * (j + N[1] * i); 
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                out[m] = in[n];
            }
        }
    }

}
// }}}


// pad_double_complex_array: {{{


void pad_double_complex_array_1d(double complex *in, double complex *out, size_t N, size_t pad)
{
    size_t Npad = N + 2 * pad;
    int n;
    for (n = 0; n < Npad; n++) {
        out[n] = 0.0;
    }
    for (n = 0; n < N; n++) {
        out[n + pad] = in[n];
    }
}

void unpad_double_complex_array_1d(double complex *in, double complex *out, size_t N, size_t pad)
{
    int n;
    for (n = 0; n < N; n++) {
        out[n] = in[n + pad];
    }
}


void pad_double_complex_array_3d(double complex *in, double complex *out, size_t *N, size_t *pad)
{
    // pad array does not initialize the out array
    // malloc + memset can do this but calloc is more efficient
    size_t rank = 3,
           N3pad = 1;
    size_t Npad[rank];
    int i, j, k, m, n;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
        N3pad *= Npad[i];
    }
    for (i = 0; i < N3pad; i++) {
        out[i] = 0;
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                m = k + N[2] * (j + N[1] * i); 
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                out[n] = in[m];
            }
        }
    }
}

void unpad_double_complex_array_3d(double complex *in, double complex *out, size_t *N, size_t *pad)
{

    size_t rank = 3;
    size_t Npad[3];
    int i, j, k, m, n;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                m = k + N[2] * (j + N[1] * i); 
                n = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0]));
                out[m] = in[n];
            }
        }
    }
}

// }}}


// pad_double_to_double_complex_array: {{{


void pad_double_to_double_complex_array_1d(double *in, double complex *out, size_t N, size_t pad)
{
    size_t Npad = N + 2 * pad;
    int n;
    for (n = 0; n < N; n++) {
        out[n + pad] = in[n];
    }
}


void unpad_double_complex_to_double_array_1d(double complex *in, double *out, size_t N, size_t pad)
{
    int n;
    for (n = 0; n < N; n++) {
        out[n] = in[n + pad];
    }
}


void pad_double_to_double_complex_array_3d(double *in, double complex *out, size_t *N, size_t *pad)
{
    size_t rank = 3,
           N3pad = 1;
    size_t Npad[rank];
    size_t i, j, k, m, n;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
        N3pad *= Npad[i];
    }
    for (i = 0; i < N3pad; i++) {
        out[i] = 0;
    } 
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                m = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0])); 
                out[m] = in[n];
            }
        }
    }
}

void unpad_double_complex_to_double_array_3d(double complex *in, double *out, size_t *N, size_t *pad)
{
    size_t rank = 3;
    size_t Npad[rank];
    size_t i, j, k, m, n;
    for (i = 0; i < rank; i++) {
        Npad[i] = N[i] + 2 * pad[i];
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                m = k + pad[2] + Npad[2] * (j + pad[1] + Npad[1] * (i + pad[0])); 
                out[n] = in[m];
            }
        }
    }
}


// }}}


// print_double_array: {{{


void print_double_array_1d(double *in, size_t N)
{
    double tmp;
    int i;
    for (i = 0; i < N; i++) {
        tmp = in[i];
        printf("%+5.15lf     ", tmp);
    }
    printf("\n");
}


void print_double_array_2d(double *in, size_t *N)
{
    double tmp;
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            tmp = in[n];
            printf("%+5.15lf     ", tmp);
        }
        printf("\n");
    }
}


void print_double_array_3d(double *in, size_t *N)
{
    double tmp;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp = in[n];
                printf("%+5.15lf     ", tmp);
            }
            printf("\n");
        }
        printf("\n");
    }
}


// }}} print_double_array


// print_double_complex_array: {{{


void print_double_complex_array_1d(double complex *in, size_t N)
{
    double re, im;
    int n;
    for (n = 0; n < N; n++) {
        re = creal(in[n]);
        im = cimag(in[n]);
        printf("%+.8f%+.8fi \n", re, im);
    }
}


void print_double_complex_array_2d(double complex *in, size_t *N)
{
    double re, im;
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            re = creal(in[n]);
            im = cimag(in[n]);
            printf("%+.8f%+.8fi \n", re, im);
        }
        printf("\n");
    }
}


void print_double_complex_array_3d(double complex *in, size_t *N)
{
    double re, im;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                re = creal(in[n]);
                im = cimag(in[n]);
                printf("%+.8f%+.8fi ", re, im);
            }
            printf("\n");
        }
        printf("\n");
    }
}


// }}} print_double_complex_array


// copy_double_array: {{{


void copy_double_array_1d(double *in, double *out, size_t N)
{
    int n;
    for (n = 0; n < N; n++) {
        in[n] = out[n];
    }
}


void copy_double_array_2d(double *in, double *out, size_t *N)
{
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            out[n] = in[n];
        }
    }
}


void copy_double_array_3d(double *in, double *out, size_t *N) 
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


// }}} copy_double_array


// copy_double_complex_array: {{{


void copy_double_complex_array_1d(double complex *in, double complex *out, size_t N)
{
    int n;
    for (n = 0; n < N; n++) {
        in[n] = out[n];
    }
}


void copy_double_complex_array_2d(double complex *in, double complex *out, size_t *N)
{
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            out[n] = in[n];
        }
    }
}


void copy_double_complex_array_3d(double complex *in, double complex *out, size_t *N)
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


// write_to_file_double_array: {{{

void write_to_file_double(char *name, double in)
{
    FILE *fdata = fopen(name, "w");
        fprintf(fdata, "%+5.6lf\n", in);
    fclose(fdata);
}



void write_to_file_double_data_array_1d(char *name, double *in, size_t N)
{
    FILE *fdata = fopen(name, "w");
    double tmp;
    int n;
    for (n = 0; n < N; n++) {
        tmp = in[n];
        fprintf(fdata, "%+5.6lf\n", tmp);
    }
    fclose(fdata);
}



void write_to_file_double_array_1d(char *name, double *in, size_t N)
{
    FILE *fdata = fopen(name, "w");
    double tmp;
    int n;
    for (n = 0; n < N; n++) {
        tmp = in[n];
        fprintf(fdata, "%+5.15lf\n", tmp);
    }
    fclose(fdata);
}


void write_to_file_double_array_2d(char *name, double *in, size_t *N)
{
    FILE *fdata = fopen(name, "w");
    double tmp;
    int i, j, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            n = j + N[1] * i;
            tmp = in[n];
            fprintf(fdata, "%+5.15lf\n", tmp);
        }
    }
    fclose(fdata);
}


void write_to_file_double_array_3d(char *name, double *in, size_t *N)
{
    FILE *fdata = fopen(name, "w");
    double tmp;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp = in[n];
                fprintf(fdata, "%+5.15lf\n", tmp);
            }
        }
    }
    fclose(fdata);
}



void write_to_file_format_double_array_3d(char *name, char *format, double *in, size_t *N)
{
    FILE *fdata = fopen(name, "w");
    double tmp;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp = in[n];
                fprintf(fdata, format, tmp);
            }
        }
    }
    fclose(fdata);
}




void write_to_file_double_data_array_3d(char *name, double *in, size_t *N)
{
    FILE *fdata = fopen(name, "w");
    double tmp;
    int i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                tmp = in[n];
                fprintf(fdata, "%+5.6lf\n", tmp);
            }
        }
    }
    fclose(fdata);
}


// }}} write_double_array_to_file


// check_sum: {{{

void check_sum_double_array_1d(char *name, double *array, size_t N)
{
    double sum = 0;
    size_t i;
    for (i = 0; i < N; i++) {
        sum += array[i];
    }
    printf("Check Sum: %s = %5.15lf\n", name, sum);
}


void check_sum_double_array_2d(char *name, double *array, size_t *N)
{
    double sum = 0;
    size_t i, j, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            m = j + N[1] * i;
            sum += array[m];
        }
    }
    printf("Check Sum: %s = %5.15lf\n", name, sum);
}


void check_sum_half_complex_3d(char *name, double *array, size_t *N)
{
    // we can fourier transform a 3d data array of length N[0] * N[1] * N[2] that is entirely real
    // using inplace real-to-complex 3d-FFT where the input array 
    // is prepared as real numbers 2 * (N[2] / 2 + 1) * N[1] * N[0] long in memory
    // where the last two numbers in the final index are padding and are ignored.
    // The output array contains complex numbers 
    // stored sequentially as pairs of numbers 2 * (N[2]/2 + 1) * N[1] * N[0] long in memory
    // check_sum_half_complex_3d will sum over the 1st and 2nd number separately
    // and print out the total of each as a complex number
    double sumre = 0,
           sumim = 0;
    size_t i, j, k, n, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < (N[2] / 2 + 1); k++) {
                n = 0 + 2 * (k + (N[2] / 2 + 1) * (j + N[1] * i));
                m = 1 + 2 * (k + (N[2] / 2 + 1) * (j + N[1] * i));
                sumre += array[n];
                sumim += array[m];
            }
        }
    }
    printf("Check Sum: %s = %5.15lf%+.15lf\n", name, sumre, sumim);
}


void check_sum_half_real_3d(char *name, double *array, size_t *N)
{
    // we can fourier transform a 3d data array of length N[0] * N[1] * N[2] of complex numbers 
    // whose postive frequencies are the complex conjugate of the negative frequencies
    // using inplace complex-to-real FFT 3d where the input array is the zero through positive complex numbers
    // stored sequentially as pairs of numbers 2 * (N[2] / 2 + 1) * N[1] * N[0] long in memory.
    // The output array is a an array of real numbers 2 * (N[2] / 2 + 1) * N[1] * N[0] long in memory
    // where the last two numbers in the final index of the output array are padding and should be ignored
    // check_sum_half_real_3d will sum over the entire array skipping the padding 
    // and print out the sum of the array
    double sum = 0;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                n = k + 2 * (N[2] / 2 + 1) * (j + N[1] * i);
                sum += array[n];
            }
        }
    }
    printf("Check Sum: %s = %5.15lf\n", name, sum);
}


void check_sum_double_array_3d(char *name, double *array, size_t *N)
{
    // for data arrays that store pairs of number
    // check sum will sum over the entire data array
    double sum = 0;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                sum += array[n];
            }
        }
    }
    printf("Check Sum: %s = %5.15lf\n", name, sum);
}



void check_sum_double_strides_3d(char *name, double *array, size_t *stride, size_t *N)
{
    // for data arrays that store pairs of numbers
    // stride[i] jumps over blocks of memory
    // stride[3] = 0 (1) selects the 1st number (2nd number) in the pair
    double sum = 0;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                n = stride[3] + stride[2] * k + stride[1] * j + stride[0] * i;
                sum += array[n];
            }
        }
    }
    printf("Check Sum: %s = %5.15lf\n", name, sum);
}

void check_sum_double_pair_strides_3d(char *name, double *array, size_t *stride, size_t *N)
{
    // for data arrays that store pairs of numbers
    // stride[i] jumps over blocks of memory
    // 0 (1) selects the 1st number (2nd number) in the pair
    double sum = 0,
           sum1 = 0;
    size_t i, j, k, n, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                n = stride[2] * k + stride[1] * j + stride[0] * i;
                m = 1 + stride[2] * k + stride[1] * j + stride[0] * i;
                sum  += array[n];
                sum1 += array[m];
            }
        }
    }
    printf("Check Sum: %s = %5.15lf%+.15lf\n", name, sum, sum1);
}


void check_sum_double_4d(char *name, double *array, size_t *N)
{
    // for data arrays that store pairs of number
    // check sum will sum over the entire data array
    double sum = 0;
    size_t i, j, k, l, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                for (l = 0; l < N[3]; l++) {
                    n = l + N[3] * (k + N[2] * (j + N[0] * i));
                    sum += array[n];
                }
            }
        }
    }
    printf("Check Sum: %s = %5.15lf\n", name, sum);
}


void check_sum_double_pair_3d(char *name, double *array, size_t *N)
{
    // data that stores a pair of numbers
    // stride = 0 select the 1st number in the pair
    // stride = 1 select the 2nd number in the pair
    double sum = 0,
           sum1 = 0;
    size_t i, j, k, n, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                n = N[3] * (k + N[2] * (j + N[1] * i));
                m = 1 + N[3] * (k + N[2] * (j + N[1] * i));
                sum  += array[n];
                sum1 += array[m];
            }
        }
    }
    printf("Check Sum: %s = %5.15lf%+.15lf\n", name, sum, sum1);
}

void check_sum_double_complex_1d(char *name, double complex *array, size_t N)
{
    double complex sum = 0;
    size_t i;
    for (i = 0; i < N; i++) {
        sum += array[i];
    }
    printf("Check Sum: %s = %.15lf%+.15lfi\n", name, Re(sum), Im(sum));
}


void check_sum_double_complex_3d(char *name, double complex *array, size_t *N)
{
    double complex sum = 0;
    size_t i, j, k, n, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) { 
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                sum += array[n];
            }
        }
    }
    printf("check sum: %s = %.15lf%+.15lfi\n", name, Re(sum), Im(sum));
}


///    }}}


// error_handling_functions: {{{

void check_pointer(void *ptr, char *file, const char *fct, int line ) 
{
    if (ptr != NULL) {
        printf("Pointer is not NULL\n");
        printf("File: %s fct: %s(), line: %d\n", file, fct, line);
        exit(0);
    }
}


void check_status(int status)
{
    if (status != 0) {
        printf("Status Fail\n");
        printf("File: %s fct: %s(), line: %d\n", __FILE__, __FUNCTION__, __LINE__);
    }
}


// }}}


#endif /* _C_EXTENSION_H */

