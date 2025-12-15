/* vim: set fdm=marker:
 *
 *
 * GSL_EXTENSION.h
 *
 *
 */

#ifndef _GSL_EXTENSION_H
#define _GSL_EXTENSION_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>


// #include "Constants.h"


#define GSL_ARRAY(z, i) ((z)[(i)])
#define GSL_VECTOR(z, i) ((z)->data[(i) * (z)->stride])


#define GSL_ARRAY_REAL(z,i) ((z)[2 * (i)])
#define GSL_ARRAY_IMAG(z,i) ((z)[2 * (i) + 1])

#define GSL_I (gsl_complex_rect(0.0, 1.0))
#define GSL_COMPLEX_I (gsl_complex_rect(0.0, 1.0))
#define GSL_IMAGINARY (gsl_complex_rect(0.0, 1.0))
#define GSL_NI (gsl_complex_rect(0.0, -1.0))
#define GSL_COMPLEX_NI (gsl_complex_rect(0.0, -1.0))
#define GSL_NIMAGINARY (gsl_complex_rect(0.0, -1.0))


// gsl_vector_complex: {{{


gsl_vector_complex * copy_vector_complex(gsl_vector_complex *cvec, size_t size)
{
    gsl_vector_complex *cvecresized = gsl_vector_complex_calloc(size);
    int i,
        imax = size;
    if(size <= 0) {
        printf("INVALID RESIZE");
    } else {
        for (i = 0; i < 2 * imax; i++) {
            GSL_VECTOR(cvecresized, i) = GSL_VECTOR(cvec, i);
        }  
    }
    return cvecresized;
}


/*
 *  copy vector to complex vector and vice versa
 */
gsl_vector_complex * copy_vector_to_vector_complex_real(gsl_vector *vec)
{
    size_t size = vec->size;
    gsl_vector_complex *cvec = gsl_vector_complex_calloc(size);
    int i,
        imax = size;
    if(size <= 0) {
        printf("INVALID RESIZE");
    } else {
        for (i = 0; i < imax; i++) {
            GSL_VECTOR_REAL(cvec, i) = GSL_VECTOR(vec, i);
        }  
    }
    return cvec;
}


gsl_vector_complex * copy_vector_to_vector_complex_imag(gsl_vector *vec)
{
    size_t size = vec->size;
    gsl_vector_complex *cvec = gsl_vector_complex_calloc(size);
    int i,
        imax = size;
    if(size <= 0) {
        printf("INVALID RESIZE");
    } else {
        for (i = 0; i < imax; i++) {
            GSL_VECTOR_IMAG(cvec, i) = GSL_VECTOR(vec, i);
        }  
    }
    return cvec;
}


gsl_vector * copy_vector_complex_real_to_vector(gsl_vector_complex *cvec)
{
    size_t size = cvec->size;
    gsl_vector *vec = gsl_vector_calloc(size);
    int i,
        imax = size;
    if(size <= 0) {
        printf("INVALID RESIZE");
    } else {
        for (i = 0; i < imax; i++) {
            GSL_VECTOR(vec, i) = GSL_VECTOR_REAL(cvec, i);
        }
    }
    return vec;
}


gsl_vector * copy_vector_complex_imag_to_vector(gsl_vector_complex *cvec)
{
    size_t size = cvec->size;
    gsl_vector *vec = gsl_vector_calloc(size);
    int i,
        imax = size;
    if(size <= 0) {
        printf("INVALID RESIZE");
    } else {
        for (i = 0; i < imax; i++) {
            GSL_VECTOR(vec, i) = GSL_VECTOR_IMAG(cvec, i);
        }
    }
    return vec;
}


gsl_vector_complex * view_vector_complex(gsl_vector_complex *cvec, int size, int offset)
{
    gsl_vector_complex *cvecresized = gsl_vector_complex_calloc(size);
    int i,
        imax = size - offset;
    if(offset + size <= 0) {
        printf("INVALID RESIZE\n");
        return NULL;
    } else {
        for (i = 0; i < 2 * imax; i++) {
            GSL_VECTOR(cvecresized, i) = GSL_VECTOR(cvec, i);
        }  
    }
    return cvecresized;
}


gsl_vector * copy_subvector(gsl_vector *vec, size_t offset, size_t size) 
{
    gsl_vector_view view = gsl_vector_subvector(vec, offset, size);
    gsl_vector *newvec = gsl_vector_alloc(size);
    gsl_vector_memcpy(newvec, &view.vector);
    return newvec;
}


gsl_vector_complex * copy_subvector_complex(gsl_vector_complex *cvec, size_t offset, size_t size) 
{
    gsl_vector_complex *newcvec = gsl_vector_complex_alloc(size);
    gsl_complex ctmp;
    int i;
    for (i = 0; i < size; i++) {
        ctmp = gsl_vector_complex_get(cvec, i + offset);
        gsl_vector_complex_set(newcvec, i, ctmp);
    }
    return newcvec;
}


// }}} gsl_vector_complex


// read_write_vectors: {{{


void write_double_to_file(FILE *fdata, char *name, double tmp) 
{
    fdata = fopen(name, "w+");
    fprintf(fdata, "%-5.15lf\n", tmp);
}


void write_vector_to_file(FILE *fdata, char *name, gsl_vector *vec) 
{
    fdata = fopen(name, "w+");
    size_t vecsize = vec->size;
    double tmpx;
    int i;
    for (i = 0; i < vecsize; i++) {
        tmpx = GSL_VECTOR(vec, i);
        fprintf(fdata, "%-5.15lf\n", tmpx);
    }
}


void write_vector_complex_to_file(FILE *fdata, char *name, gsl_vector_complex *cvec) 
{
    fdata = fopen(name, "w+");
    size_t vecsize = cvec->size;
    double tmpx,
           tmpy;
    int i;
    for (i = 0; i < vecsize; i++) {
        tmpx = GSL_VECTOR_REAL(cvec, i);
        tmpy = GSL_VECTOR_IMAG(cvec, i);
        fprintf(fdata, "%-5.15lf %-5.15lf \n", tmpx, tmpy);
    }
}


// }}} read_write_vectors


// print_data_types: {{{ 

void print_vector(gsl_vector *vec) 
{
    size_t vecsize = vec->size;
    double tmpx;
    int i;
    for (i = 0; i < vecsize; i++) {
        tmpx = GSL_VECTOR(vec, i);
        printf("%-5.15lf\n", tmpx);
    }
}


void print_vector_address(gsl_vector *vec) 
{
    size_t vecsize = vec->size;
    int i;
    for (i = 0; i < vecsize; i++) {
        printf("%p\n", &GSL_VECTOR(vec, i));
    }
}


void print_vector_complex(gsl_vector_complex *cvec) 
{
    size_t vecsize = cvec->size;
    double tmpx,
           tmpy;
    int i;
    for (i = 0; i < vecsize; i++) {
        tmpx = GSL_VECTOR_REAL(cvec, i);
        tmpy = GSL_VECTOR_IMAG(cvec, i);
        printf("%-5.15lf %-5.15lf \n", tmpx, tmpy);
    }
}


void print_vector_complex_address(gsl_vector_complex *cvec) 
{
    size_t vecsize = cvec->size;
    int i;
    for (i = 0; i < vecsize; i++) {
        printf("%p %p \n", &GSL_VECTOR_REAL(cvec, i), &GSL_VECTOR_IMAG(cvec, i));
    }
}


// }}} print_data_types


// plot_functions: {{{ 


gsl_vector * plot_function(gsl_function F, double xmin, double xmax, size_t Nmax) 
{

            double tmp;
            gsl_vector *vec = gsl_vector_alloc(Nmax);
            int i;
            for (i = 0; i < Nmax; i++) {
                tmp = xmin + (xmax - xmin) * i / (Nmax - 1);
                GSL_VECTOR(vec, i) = GSL_FN_EVAL(&F, tmp);
            }
            return vec;
}


// }}} plot_ functions


// gsl_vector: {{{ 


double gsl_complex_real_get(gsl_complex ctmp)
{
    return GSL_REAL(ctmp);
}


double gsl_complex_imag_get(gsl_complex ctmp)
{
    return GSL_IMAG(ctmp);
}


void gsl_complex_real_set(gsl_complex ctmp, double dtmp)
{
    GSL_SET_REAL(&ctmp, dtmp);
}


void gsl_complex_imag_set(gsl_complex ctmp, double dtmp)
{
    GSL_SET_IMAG(&ctmp, dtmp);
}


// FIXME these function are BASE depedent that needs to be remedied
// but maybe that's fine we just need to make more functions
/****/
double gsl_vector_complex_real_get(gsl_vector_complex *cvec, size_t i)
{
    gsl_complex ctmp;
    ctmp = gsl_vector_complex_get(cvec, i);
    return GSL_REAL(ctmp);
}


double gsl_vector_complex_imag_get(gsl_vector_complex *cvec, size_t i)
{
    gsl_complex ctmp;
    ctmp = gsl_vector_complex_get(cvec, i);
    return GSL_IMAG(ctmp);
}


void gsl_vector_complex_real_set(gsl_vector_complex *cvec, size_t i, double dtmp)
{
    gsl_complex ctmp;
    ctmp = gsl_vector_complex_get(cvec, i);
    GSL_SET_REAL(&ctmp, dtmp);
    gsl_vector_complex_set(cvec, i, ctmp);
}


void gsl_vector_complex_imag_set(gsl_vector_complex *cvec, size_t i, double dtmp)
{
    gsl_complex ctmp;
    ctmp = gsl_vector_complex_get(cvec, i);
    GSL_SET_IMAG(&ctmp, dtmp);
    gsl_vector_complex_set(cvec, i, ctmp);
}


double gsl_vector_complex_REAL_get(gsl_vector_complex *cvec, size_t i)
{    
    gsl_vector_view view = gsl_vector_complex_real(cvec);
    return gsl_vector_get(&view.vector, i);
}


double gsl_vector_complex_IMAG_get(gsl_vector_complex *cvec, size_t i)
{
    gsl_vector_view view = gsl_vector_complex_imag(cvec);
    return gsl_vector_get(&view.vector, i);
}


void gsl_vector_complex_REAL_set(gsl_vector_complex *cvec, size_t i, double dtmp)
{
    gsl_vector_view view = gsl_vector_complex_real(cvec);
    gsl_vector_set(&view.vector, i, dtmp);
}


void gsl_vector_complex_IMAG_set(gsl_vector_complex *cvec, size_t i, double dtmp)
{
    gsl_vector_view view = gsl_vector_complex_imag(cvec);
    gsl_vector_set(&view.vector, i, dtmp);
}
/*****/


double * gsl_vector_complex_REAL_ptr(gsl_vector_complex *cvec, size_t i)
{
    gsl_vector_view view = gsl_vector_complex_real(cvec);
    return gsl_vector_ptr(&view.vector, i);

}


double * gsl_vector_complex_IMAG_ptr(gsl_vector_complex *cvec, size_t i)
{
    gsl_vector_view view = gsl_vector_complex_imag(cvec);
    return gsl_vector_ptr(&view.vector, i);
}


double * gsl_vector_complex_real_ptr(gsl_vector_complex *cvec, size_t i)
{    
    return &cvec->data[i * cvec->stride];
}


double * gsl_vector_complex_imag_ptr(gsl_vector_complex *cvec, size_t i)
{
    return (double *) (cvec->data + i * cvec->stride + 1);
}


// pad the ends of a complex vector
gsl_vector_complex * gsl_vector_complex_pad(gsl_vector_complex *cvec, size_t padsize)
{
    size_t vecsize = cvec->size,
           padvecsize = 2 * padsize  + cvec->size;
    gsl_vector_complex *cvecres = gsl_vector_complex_calloc(padvecsize);
    gsl_complex ctmp;
    int i;
    for (i = 0; i < vecsize; i++) {
        ctmp = gsl_vector_complex_get(cvec, i);
        gsl_vector_complex_set(cvecres, i + padsize, ctmp);
    }
    return cvecres;
}


// unpad the ends of a complex vector 
gsl_vector_complex * gsl_vector_complex_unpad(gsl_vector_complex *cvec, size_t padsize)
{
    size_t unpadvecsize = cvec->size - 2 * padsize;
    gsl_vector_complex *cvectmp = gsl_vector_complex_calloc(unpadvecsize);
    gsl_complex ctmp;
    int i;
    for (i = 0; i < unpadvecsize; i++) {
        ctmp = gsl_vector_complex_get(cvec, i + padsize);
        gsl_vector_complex_set(cvectmp, i, ctmp);
    }
    return cvectmp;
}


// pad gsl 
gsl_vector * gsl_vector_pad(gsl_vector *vec, size_t padsize)
{
    size_t vecsize = vec->size,
           padvecsize = 2 * padsize + vec->size;
    gsl_vector *vectmp = gsl_vector_calloc(padvecsize);
    double tmp;
    int i,
        imax = vecsize;
    for (i = 0; i < imax; i++) {
        tmp = gsl_vector_get(vec, i);
        gsl_vector_set(vectmp, i + padsize, tmp);
    }
    return vectmp;
}


gsl_vector * gsl_vector_unpad(gsl_vector *vec, size_t padsize)
{
    size_t unpadvecsize = vec->size - 2 * padsize;
    gsl_vector *vectmp = gsl_vector_calloc(unpadvecsize);
    double tmp;
    int i,
        imax = unpadvecsize;
    for (i = 0; i < imax; i++) {
        tmp = gsl_vector_get(vec, i + padsize);
        gsl_vector_set(vectmp, i, tmp);
    }
    return vectmp;
}


// }}} gsl_vector


#endif /* _GSL_EXTENSION_H */

