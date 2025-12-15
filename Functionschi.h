/* vim: set fdm=marker:
 *
 *
 * Functionschi.h
 *
 * In this version we optimize all the functions.
 *
 */

#ifndef _FUNCTIONSCHI_H
#define _FUNCTIONSCHI_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>


#include "C_extension.h"


double THETA(double x) {
    return (1 + tanh(x / 0.005)) / 2;
}


double epsilon_k(double t, double tp, double tpp, double kx, double ky)
{
    return - 2 * t * (cos(kx) + cos(ky))
           - 4 * tp * cos(kx) * cos(ky) 
           - 2 * tpp * (cos(2 * kx) + cos(2 * ky));
}

double derivative_epsilon_kx(double t, double tp, double tpp, double kx, double ky)
{
    return   2 * t * sin(kx) 
           + 4 * tp * sin(kx) * cos(ky)
           + 4 * tpp * sin(2 * kx);
}

double derivative_epsilon_ky(double t, double tp, double tpp, double kx, double ky)
{
    return   2 * t * sin(ky) 
           + 4 * tp * cos(kx) * sin(ky)
           + 4 * tpp * sin(2 * ky);
}

double derivative_epsilon_kxx(double t, double tp, double tpp, double kx, double ky)
{
    return   2 * t * cos(kx) 
           + 4 * tp * cos(kx) * cos(ky)
           + 8 * tpp * cos(2 * kx);
}

double derivative_epsilon_kyy(double t, double tp, double tpp, double kx, double ky)
{
    return   2 * t * cos(ky) 
           + 4 * tp * cos(kx) * cos(ky)
           + 8 * tpp * cos(2 * ky);
}

double derivative_epsilon_kxy(double t, double tp, double tpp, double kx, double ky)
{
    return -4 * tp * sin(kx) * sin(ky);
}

double J_k(double J, double kx, double ky)
{
    return 2.0 * J * (cos(kx) + cos(ky));
}

// Fermi_Energy Function: {{{

double find_density(double Fermi_energy, double *en, size_t *N) 
{
    size_t i, j, m,
           Ns = N[0] * N[1];
    double density = 0;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = j + N[1] * i;
            density += THETA(Fermi_energy - en[m]);
        }
    }
    density /= Ns;
    return density; 
}


// find_Fermi_energy_root_function_params
struct FindFermiEnergyRootFunctionParams
{
    size_t *N;
    double nd;
    double *en;
};


double find_Fermi_energy_root_function(double Fermi_energy, void *params) 
{
    struct FindFermiEnergyRootFunctionParams *p = (struct FindFermiEnergyRootFunctionParams *) params;

    size_t *N = p->N;
    double nd = p->nd;
    double *en = p->en;

    return nd / 2 - find_density(Fermi_energy, en, N); 
}

// }}}


// k_Fermi function: {{{
struct k_FermiRootFunctionParams
{
    double t,
           tp,
           tpp,
           Fermi_energy;
};

double k_Fermi_nodal_root_function(double x, void *params)
{
    struct k_FermiRootFunctionParams *p = (struct k_FermiRootFunctionParams *) params;

    double t = p->t,
           tp = p->tp,
           tpp = p->tpp,
           Fermi_energy = p->Fermi_energy;

    return Fermi_energy - epsilon_k(t, tp, tpp, x, x);
}


double k_Fermi_antinodal_root_function(double x, void *params)
{
    struct k_FermiRootFunctionParams *p = (struct k_FermiRootFunctionParams *) params;

    double t = p->t,
           tp = p->tp,
           tpp = p->tpp,
           Fermi_energy = p->Fermi_energy;

    return Fermi_energy - epsilon_k(t, tp, tpp, x, 0);
}

// }}}


void localize_function(double *rho_g, double *rho_g_loc, size_t *N) {
    size_t i, j, k, n,
           Ns = N[0] * N[1];
    for (k = 0; k < N[2]; k++) {
        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                 n = k + N[2] * (j + N[1] * i);
                 rho_g_loc[k] += rho_g[n];
            }
        }
        rho_g_loc[k] /= Ns;
     }
}


double sum_function(double *rho_g_temp, double DELTAomega, size_t *N)
{
    size_t i, j, k, n,
           Ns = N[0] * N[1];
    double  sum = 0;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                sum += rho_g_temp[n];
            }
        }
    }
    sum *= (DELTAomega / Ns);
    return sum;
}


double Fermi_function(double x, double tau)
{   // tau is the natural temperature
    if (x / tau < -20.0) {
        return 1.0;
    }
    else if (x / tau > 20.0) {
        return 0.0;
    }
    else if (x / tau > -20.0 &&  x / tau < 20.0) {
        return 1.0 / (exp(x / tau) + 1);
    }
    else {
        return 0.0;
    }
}


double derivative_Fermi_function_over_beta(double x, double tau)
{
    if (x / tau < -20.0) {
        return 0.0;
    }
    else if (x / tau > 20.0) {
        return 0.0;
    }
    else if (x / tau > -20.0 && x / tau < 20.0) {
        return -1.0 / pow(2.0 * cosh(x / (2.0 * tau)), 2);
    }
    else {
        return 0.0;
    }
}


double Lorentz_function(double x, double x0, double gamma)
{
    return (gamma / (2.0 * M_PI)) / (pow(x - x0, 2) + pow(gamma / 2.0, 2));
}



double pseudo_Fermi_surface(double *rho_bigG, double *k_array, double *omega_array, double DELTAomega, double tau, size_t *N) {

    size_t i, j, k, m, n;
    size_t Ns = N[0] * N[1];
    double gamma = 0,
           norm = 0,
           density = 0,
           omega = 0,
           temp = 0;
    double *sum_array_2d = calloc(Ns, sizeof(double)),
           *norm_array_2d = calloc(Ns, sizeof(double));
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                omega = omega_array[k];
                temp = rho_bigG[n] / cosh(omega / (2 * tau));
                sum_array_2d[m] += temp * omega;
                norm_array_2d[m] += temp;
            }
            sum_array_2d[m] *= DELTAomega;
            norm_array_2d[m] *= DELTAomega;
            gamma = - sum_array_2d[m] / norm_array_2d[m];
            density += THETA(gamma);
        }
    }
    density /= Ns; 
    free(sum_array_2d);
    free(norm_array_2d);
    return density;
}


// single_root_initial_params
struct RootFunctionParams
{
    double nd,
           DELTAomega;
    int neta;
    double *omega_array,
           *Fermi_array,
           *en,
           *rho_g_initial;
    size_t *N;
};


double root_function(double mup_initial, void *params)
{
    struct RootFunctionParams *p = (struct RootFunctionParams *) params;

    double nd = p->nd,
           DELTAomega = p->DELTAomega;
    int neta = p->neta;
    double *omega_array = p->omega_array,
           *Fermi_array = p->Fermi_array,
           *en = p->en,
           *rho_g = p->rho_g_initial;
    size_t *N = p->N;

    size_t Ns = N[0] * N[1],
           N3 = N[0] * N[1] * N[2];
    size_t i, j, k, n, m;


    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                rho_g[n] = Lorentz_function(omega_array[k], en[m] - mup_initial, pow(2, neta) * DELTAomega);
            }
        }
    }


    double *density_of_states = calloc(N3, sizeof(double));

    for (k = 0; k < N[2]; k++) {
        for (j = 0; j < N[1]; j++) {
            for (i = 0; i < N[0]; i++) {
                n = k + N[2] * (j + N[1] * i);
                density_of_states[k] += rho_g[n];
             }
        }
        density_of_states[k] /= Ns;
    }


    double sum = 0;
    for (k = 0; k < N[2]; k++) {
        sum += density_of_states[k] * Fermi_array[k];
    }
    sum *= DELTAomega;
    free(density_of_states);
    return nd / 2 - sum;
}


void rho_g_and_rho_bigG_formulas(double lambda,
                                 double nd,
                                 double gamma_SE,
                                 double mup,
                                 double u0,
                                 double *omega_array,
                                 double *en,
                                 double *xi,
                                 double *rho_g,
                                 double *rho_bigG,
                                 double *rho_psi0,
                                 double *rho_psi1,
                                 double *rho_chi0,
                                 double *rho_chi1,
                                 double *rho_chi2,
                                 double *re_psi0,
                                 double *re_psi1,
                                 double *re_chi0,
                                 double *re_chi1,
                                 double *re_chi2,
                                 size_t *N)
{
    size_t i, j, k, n, m;
    size_t N3 = N[0] * N[1] * N[2];

    double rho_psi,
           rho_chi,
           rho_phi,
           re_psi,
           re_chi,
           re_phi,
           X = (-u0/2),
           rr,
           denom,
           re_g;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                rho_psi     = rho_psi0[n] + rho_psi1[n] * X;
                re_psi      = re_psi0[n]  + re_psi1[n]  * X;
                rho_chi     = rho_chi0[n] + rho_chi1[n] * X + rho_chi2[n] * pow(X, 2);
                re_chi      = re_chi0[n]  + re_chi1[n]  * X + re_chi2[n]  * pow(X, 2);
                rho_phi     = rho_chi     + (en[m] + X) * rho_psi; 
                re_phi      = re_chi      + (en[m] + X) * re_psi;
                rr          = omega_array[k] + mup - xi[m] - re_phi;
                denom       = pow(rr, 2) + pow(M_PI * rho_phi, 2);
                re_g        = rr / denom;
                rho_g[n]    = rho_phi / denom;
                rho_bigG[n] = rho_g[n] * (1 - lambda * gamma_SE + re_psi) + re_g * rho_psi;
            }
        }
    }

}

double density_function(double *rho_g_temp, double *Fermi_array, double DELTAomega, size_t *N)
{
    size_t i, j, k, n;
    size_t Ns = N[0] * N[1];
    double  density = 0;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                density += rho_g_temp[n] * Fermi_array[k];
            }
        }
    }
    density *= (DELTAomega / Ns);
    return density;
}


double u0_sum_rule(double lambda,
                   double nd,
                   double J,
                   double gamma_SE,
                   double mup,
                   double u0,
                   double DELTAomega,
                   double gfe,
                   double ng,
                   double *omega_array, 
                   double *Fermi_array, 
                   double *en,
                   double *rho_bigG,
                   int DENSITYAPPROX,
                   size_t *N) {

// notes on the definition of mup and xip_k
// e_k = epsilon_k; ep_k = (e_k - u0/2); ebar_k = e_k * aG;
// bare band: e_k - mu
// ECFL band: xi[0]_k - mup[0]: xi[0]_k = (e_k - u0/2) * aG, mup[0] = mu - (u0/2) + (J0/2) * gamma_SE
// ECFL w/ SE: xibar_k - mu, xibar_k = xi[0]_k + chi[0]_k, mubar = mup[0] - chi[0]_c
// ECFL w/ SE & u0-indep energy: xip_k - mup, xip_k = xibar_k + u0/2 * aG, mup = mubar + u0/2 * aG

// xip_k = e_k * aG + chi[0]_k
// mup = mup[0] - lambda * chi[0]_{c} + (u0 / 2) * aG
// mup[0] = mu - (u0/2) + (J0/2) * lambda * gamma_SE, aG = 1 - lambda * gamma_SE
// gamma_SE = nG / 2 
// chi[0] = chi[0]_{c} + chi[0]_k
// chi[0]_{c} = - gfe + (u0 / 2) * (ng / 2), where gfe = \sum_{p} g_{p} e_p and (ng / 2) = \sum_{p} g_{p} 
// chi[0]_{k} = -(1 / 2) \sum_{p} J_{k-p} 

// mu = mup + (u0/2) - lambda * (J0 / 2) * gamma_SE -(u0/2)*aG + lambda * chi[0]_{c}

    double density = 0;
    if (DENSITYAPPROX) {
        density = nd;
    }
    else {
        density = ng;
    }


    size_t Ns = N[0] * N[1];
    double J0 = J_k(J, 0, 0),
           aG = 1 - lambda * gamma_SE,
           chi0c = -gfe + (u0 / 2) * (density / 2),
           mu = mup + (u0/2) - lambda * gamma_SE * (J0/2) - (u0/2) *  aG + lambda * chi0c; 
//    printf("J0 = %lf, aG = %lf, gamma_SE = %lf, chi0c = %lf, mu = %lf\n", J0, aG, gamma_SE, chi0c, mu);


    size_t i, j, k, n, m;
    double nG = 0,
           sum = 0,
           firstmoment = 0,
           averagekineticenergy = 0,
           expectationvaluemu = 0;
           
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                nG += rho_bigG[n] * Fermi_array[k];
                firstmoment += rho_bigG[n] * Fermi_array[k] * omega_array[k];
                averagekineticenergy += rho_bigG[n] * Fermi_array[k] * en[m];
                expectationvaluemu += rho_bigG[n] * Fermi_array[k] * mu;
                sum += rho_bigG[n] * Fermi_array[k] * (omega_array[k] + mu - en[m]); 
            }
        }
    }
    sum *= (DELTAomega / Ns);
    nG *= 2 * (DELTAomega / Ns);
    firstmoment *= (DELTAomega / Ns);
    averagekineticenergy *= (DELTAomega / Ns);
    expectationvaluemu *= (DELTAomega / Ns);

    printf("ng = %5.6lf\n", ng);
    printf("nG = %5.6lf\n", nG);
    printf("mu true = %5.6lf\n", mu);
    printf("firstmoment = %5.6lf\n", firstmoment);
    printf("averagekineticenergy = %5.6lf\n", averagekineticenergy);
    printf("expectationvaluemu = %5.6lf\n", expectationvaluemu);
    printf("total = %5.6lf\n", firstmoment-averagekineticenergy + expectationvaluemu);
    printf("sum = %5.6lf\n", sum);


    return sum;
}



// singleroot_increment_params
struct SinglerootFunctionParams
{
    double lambda,
           nd,
           J,
           gamma_SE,
           u0,
           DELTAomega,
           gfe;
    double *omega_array,
           *Fermi_array,
           *en,
           *xi,
           *rho_psi0,
           *rho_psi1,
           *rho_chi0,
           *rho_chi1,
           *rho_chi2,
           *re_psi0,
           *re_psi1,
           *re_chi0,
           *re_chi1,
           *re_chi2;
    size_t *N;
};


double singleroot_function(double mup, void *params)
{
    struct SinglerootFunctionParams *p = (struct SinglerootFunctionParams *) params;

    double lambda       = p->lambda,
           nd           = p->nd,
           gamma_SE     = p->gamma_SE,
           u0           = p->u0,
           DELTAomega   = p->DELTAomega,
           gfe          = p->gfe;
    double *omega_array = p->omega_array,
           *Fermi_array = p->Fermi_array,
           *en          = p->en,
           *xi          = p->xi,
           *rho_psi0    = p->rho_psi0,
           *rho_psi1    = p->rho_psi1,
           *rho_chi0    = p->rho_chi0,
           *rho_chi1    = p->rho_chi1,
           *rho_chi2    = p->rho_chi2,
           *re_psi0     = p->re_psi0,
           *re_psi1     = p->re_psi1,
           *re_chi0     = p->re_chi0,
           *re_chi1     = p->re_chi1,
           *re_chi2     = p->re_chi2;
    size_t *N           = p->N;
    size_t Ns = N[0] * N[1],
           N3 = N[0] * N[1] * N[2];
    // x0 = mup and x1 = u0

    double *rho_g    = calloc(N3, sizeof(double)),
           *rho_bigG = calloc(N3, sizeof(double));

    rho_g_and_rho_bigG_formulas(lambda,
                                nd,
                                gamma_SE,
                                mup,
                                u0,
                                omega_array,
                                en,
                                xi,
                                rho_g,
                                rho_bigG,
                                rho_psi0,
                                rho_psi1,
                                rho_chi0,
                                rho_chi1,
                                rho_chi2,
                                re_psi0,
                                re_psi1,
                                re_chi0,
                                re_chi1,
                                re_chi2,
                                N);

 

//    double ng = density_function(rho_g,    Fermi_array, DELTAomega, N) - nd / 2;
    double  nbigG = 2 * density_function(rho_bigG, Fermi_array, DELTAomega, N) - nd / 2;
    

//    double y0 = ng - nd;
    double y0 = nbigG  - nd ;

    free(rho_g);
    free(rho_bigG);

    return y0;
}


// multiroot_increment_params
struct MultirootFunctionParams
{
    size_t *N;
    double lambda,
           nd,
           J,
           gamma_SE,
           DELTAomega,
           gfe;
    double *omega_array,
           *Fermi_array,
           *en,
           *xi,
           *rho_psi0,
           *rho_psi1,
           *rho_chi0,
           *rho_chi1,
           *rho_chi2,
           *re_psi0,
           *re_psi1,
           *re_chi0,
           *re_chi1,
           *re_chi2;
    int U0SUMRULE,
        DENSITYAPPROX;
};



int multiroot_function(const gsl_vector *x, void *params, gsl_vector *f)
{
    struct MultirootFunctionParams *p = (struct MultirootFunctionParams *) params;

    size_t *N           = p->N;
    double lambda       = p->lambda,
           nd           = p->nd,
           J            = p->J,
           gamma_SE     = p->gamma_SE,
           DELTAomega   = p->DELTAomega,
           gfe          = p->gfe;
    double *omega_array = p->omega_array,
           *Fermi_array = p->Fermi_array,
           *en          = p->en,
           *xi          = p->xi,
           *rho_psi0    = p->rho_psi0,
           *rho_psi1    = p->rho_psi1,
           *rho_chi0    = p->rho_chi0,
           *rho_chi1    = p->rho_chi1,
           *rho_chi2    = p->rho_chi2,
           *re_psi0     = p->re_psi0,
           *re_psi1     = p->re_psi1,
           *re_chi0     = p->re_chi0,
           *re_chi1     = p->re_chi1,
           *re_chi2     = p->re_chi2;
    int U0SUMRULE    = p->U0SUMRULE,
        DENSITYAPPROX = p->DENSITYAPPROX;

    size_t Ns = N[0] * N[1],
           N3 = N[0] * N[1] * N[2];
    // x0 = mup and x1 = u0
    const double mup = gsl_vector_get(x, 0),
                 u0 = gsl_vector_get(x, 1);

    double *rho_g    = calloc(N3, sizeof(double)),
           *rho_bigG = calloc(N3, sizeof(double));

    rho_g_and_rho_bigG_formulas(lambda,
                                nd,
                                gamma_SE,
                                mup,
                                u0,
                                omega_array,
                                en,
                                xi,
                                rho_g,
                                rho_bigG,
                                rho_psi0,
                                rho_psi1,
                                rho_chi0,
                                rho_chi1,
                                rho_chi2,
                                re_psi0,
                                re_psi1,
                                re_chi0,
                                re_chi1,
                                re_chi2,
                                N);

    double ng = 2 * density_function(rho_g, Fermi_array, DELTAomega, N),
           nbigG = 2 * density_function(rho_bigG, Fermi_array, DELTAomega, N);
    double z0,
           z1,
           u0SumRule;

    if (U0SUMRULE) {

        u0SumRule = u0_sum_rule(lambda,
                                nd,
                                J,
                                gamma_SE,
                                mup,
                                u0,
                                DELTAomega,
                                gfe,
                                ng,
                                omega_array, 
                                Fermi_array, 
                                en,
                                rho_bigG,
                                
                                DENSITYAPPROX,
                                
                                N);

        z0 = nbigG - nd,
        z1 = u0SumRule + nd * nd * J / 2;
    }
    else {
        z0 = ng / 2 - nd / 2;
        z1 = nbigG / 2 - nd / 2;
    }
 

     const double y0 = z0,
                  y1 = z1;

    gsl_vector_set(f, 0, y0);
    gsl_vector_set(f, 1, y1);

    free(rho_g);
    free(rho_bigG);

    return GSL_SUCCESS;
}


double our_window(double omega, double omega0, double omega1)
{
	double uu = 0;
    if (fabs(omega) <= omega0) {
        uu = 1.0;
    } else if (fabs(omega) > omega0 && fabs(omega) < omega1) {
        uu = (1.0 / 2.0) + (1.0 / 2.0) * sin((M_PI / 2.0) * (omega0 + omega1 - 2 * fabs(omega)) / (omega1 - omega0));
    } else if(fabs(omega) >= omega1) {
        uu = 0.0;
    } else {
        uu = 0.0;
    }
    return uu;
}


void tukey_array_1d_function(double *omega_array,
                             double *tukey_array_1d,
                             double OMEGAcminus,
                             double OMEGAcplus,
                             size_t Nomega)
{
    size_t k;
    for (k = 0; k < Nomega; k++) {
        tukey_array_1d[k] = our_window(omega_array[k], OMEGAcminus, OMEGAcplus);
    }
}

void tukey_array_3d_function(double *in,
                             double *out,
                             double *tukey_array_1d,
                             double DELTAomega,
                             size_t *N)
{
    size_t i, j, k, n;
    size_t N3 = N[0] * N[1] * N[2];


    double renorm = 0;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                renorm += in[n] * tukey_array_1d[k];
            }
        }
    }
    renorm *= (DELTAomega / pow(N[0], 2));
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                out[n] = in[n] * tukey_array_1d[k] / renorm;
            }
        }
    }
}

void tukey_max_array_3d_function(double *in,
                                 double *out,
                                 double *omega_array,
                                 double DELTAomega,
                                 double OMEGAcminus,
                                 double OMEGAcplus,
                                 int TUKEY,
                                 size_t *N)
{
    size_t i, j, k, l, m, n;
    size_t Ns = N[0] * N[1],
           N3 = N[0] * N[1] * N[2];

    size_t max_index = 0,
           stride = 1;
    double *omega_max = calloc(Ns, sizeof(double)),
           *tukey_max_array = calloc(N3, sizeof(double)),
           *renorm = calloc(Ns, sizeof(double));

    if (TUKEY) {
        for (i = 0; i < N[0]; i++) {
            l = N[1] * i;
            for (j = 0; j < N[1]; j++) {
                m = j + l;
                max_index = gsl_stats_max_index(in + N[2] * m, stride, N[2]);
                omega_max[m] = omega_array[max_index];
            }
        }

        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = j + N[1] * i;
                for (k = 0; k < N[2]; k++) {
                    n = k + N[2] * m;
                    tukey_max_array[n] = our_window(omega_array[k] - omega_max[m], OMEGAcminus, OMEGAcplus);
                }
            }
        }

        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = j + N[1] * i;
                renorm[m] = 0;
                for (k = 0; k < N[2]; k++) {
                    n = k + N[2] * m;
                    renorm[m] += in[n] * tukey_max_array[n];
                }
                renorm[m] *= DELTAomega;
//       printf("omega_max[%zu]=%2.6lf, renorm[%zu]=%2.6lf\n", m, omega_max[m], m, renorm[m]);
            }
        }


        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = (j + N[1] * i);
                for (k = 0; k < N[2]; k++) {
                    n = k + N[2] * m;
                    out[n] = in[n] * tukey_max_array[n] / renorm[m];
                }
            }
        }
    }
    else {
        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = (j + N[1] * i);
                for (k = 0; k < N[2]; k++) {
                    n = k + N[2] * m;
                    out[n] = in[n];
                }
            }
        }
    }

    free(renorm);
    free(omega_max);
    free(tukey_max_array);
}

// returns weighted average of the 3d array
void weighted_average_array_3d(double *rho_g_temp, double *rho_g_final, double weight, size_t *N)
{


    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                rho_g_temp[n] = (1 - weight) * rho_g_temp[n] + weight * rho_g_final[n];
             }
        }
    }


}




// average two array returning  result is returned in the first parameter
void average_array_3d(double *rho_g_temp, double *rho_g_final, size_t *N)
{
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                rho_g_temp[n] = (rho_g_temp[n] + rho_g_final[n]) / 2;
             }
        }
    }


}


double calculate_gerror(double *rho_g_final, double *rho_g_temp, size_t *N)
{
    double gerror = 0;
    size_t i, j, k, n;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                gerror += fabs(rho_g_final[n] - rho_g_temp[n]);
             }
        }
    }
    gerror /= N[0] * N[1] * N[2];
    return gerror;
}


// in is signal array
// peak_value is an array containing the peak value for each point in k-space
// peak_index is an array containing the index location for each point in k-space
// gamma is an array containing the half-width of the peak for each point in k-space
void find_GAMMA(double *in, double *peak_value, size_t *peak_index, double *GAMMA, double DELTAomega, size_t *N) {

    size_t halfpeak_index_left, 
           halfpeak_index_right,
           max_index,
           stride = 1;

    double halfpeak_value;
    size_t i, j, k, n, m;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = j + N[1] * i;
            max_index = gsl_stats_max_index(in + N[2] * m, stride, N[2]);
            peak_index[m] = max_index;
            peak_value[m] = in[max_index];
            halfpeak_value = peak_value[m] / 2;
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                if (in[n] > halfpeak_value) {
                    halfpeak_index_left = k;
                    break;
                }
            }
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                if (in[n] < halfpeak_value) {
                    halfpeak_index_right = k;
                    break;
                }
            }
            GAMMA[m] = (halfpeak_index_right - halfpeak_index_left) * DELTAomega / 2.0;
        }
    }

}

#endif /* _FUNCTIONS_H */
