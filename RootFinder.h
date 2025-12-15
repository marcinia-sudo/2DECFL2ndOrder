/* vim: set fdm=marker:
 *
 * RootFinder.h
 *
 *
 */

#ifndef _ROOT_FINDER_H
#define _ROOT_FINDER_H
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include "C_extension.h"


struct FindRootParams 
{ 
    double x_lo, 
           x_hi;
    int SWITCH;
};


double find_root(gsl_function *F, void *params)
{
    struct FindRootParams *p = (struct FindRootParams *) params;

    double x_lo = p->x_lo,
           x_hi = p->x_hi;
 // ON -> 1 and OFF -> 0 is defined C_extensions.h
    int PRINT = p->SWITCH;


    int status;
    int iter = 0,
        max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    T = gsl_root_fsolver_brent;
//    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc(T);

    gsl_set_error_handler_off();
    gsl_root_fsolver_set(s, F, x_lo, x_hi);

    if (PRINT) {
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.00001);
        printf("%d\n", status);
        printf ("using %s method\n", gsl_root_fsolver_name (s));
        printf ("%5s [%9s, %9s] %9s %9s\n", "iter", "lower", "upper", "root", "err(est)");
    }

    double root = 0;
    do {
        status = gsl_root_fsolver_iterate(s);
        root = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.00001);

        if (PRINT) {
            if (status == GSL_SUCCESS) {
                printf ("Converged: \n");
            }
            printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, root, x_hi - x_lo);
        }
        iter++;
     } while (status == GSL_CONTINUE && iter < max_iter);
     gsl_root_fsolver_free(s);
     return root;
}


#endif /* _ROOT_FINDER_H */
