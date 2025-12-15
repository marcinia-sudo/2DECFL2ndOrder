/* vim: set fdm=marker:
 *
 * MultirootFinder.h
 * Computes the roots of a set of equations
 *
 */

#ifndef _MULTIROOT_FINDER_H
#define _MULTIROOT_FINDER_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

//#include "Constants.h"
#include "C_extension.h"
#include "gsl_extension.h"


void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
    printf ("iter = %3zu x = % .3f x2 = % .3f\t"
            "f1(x1, x2) = % .3e f2(x1, x2) = % .3e\n",
            iter,
            gsl_vector_get(s->x, 0),
            gsl_vector_get(s->x, 1),
            gsl_vector_get(s->f, 0),
            gsl_vector_get(s->f, 1));
}

void print_state_n(size_t iter, gsl_multiroot_fsolver * s, size_t n)
{
    int i,
        imax = n;
    printf("iter = %3zu x = ", iter);
    for (i = 0; i < imax; i++) {
        printf("% .3f ", gsl_vector_get(s->x, i));
    }
    printf("\t");
    printf("f(x) = ");
    for (i = 0; i < imax; i++) {
        printf("% .3e ", gsl_vector_get(s->f, i));
    }
    printf("\n");
}


void multiroot_finder(gsl_vector *x, gsl_multiroot_function *F, gsl_vector *y)
{
    // ON = 1, OFF = 0 defined in C_extension.h
    int PRINT = ON;
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    int status;
    size_t iter = 0;
    const size_t n = F->n;

    gsl_multiroot_function f = *F;
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, n);
    gsl_multiroot_fsolver_set(s, &f, x);
    if (PRINT) {
        print_state_n(iter, s, n);
    }

    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        if (PRINT)
            print_state_n (iter, s, n);
        if (status)   /* check if solver is stuck */
            break;
        status = gsl_multiroot_test_residual(s->f, 1e-9);
    } while (status == GSL_CONTINUE && iter < 1000);
    if (PRINT) {
        printf ("status = %s\n", gsl_strerror(status));
    }
    GSL_VECTOR(y, 0) = gsl_vector_get(s->x, 0);
    GSL_VECTOR(y, 1) = gsl_vector_get(s->x, 1);
    gsl_multiroot_fsolver_free (s);
}


#endif /* _MULTIROOT_FINDER_H */
