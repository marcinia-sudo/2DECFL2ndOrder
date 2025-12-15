/* vim: set fdm=marker:
 *
 *
 *  2D ECFL - 2nd order
 *
 * This program uses Functionchi.h, SelfEnergyR2CDFTIchi.h and HilbertTransformR2CDFTI.h
 * optimizations. We can modify our calculations to minimize either WALL-Time, CPU-Time, or RAM usage.
 *
 *
 *
 */
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_statistics_double.h>

#include "C_extension.h"
#include "SelfEnergyR2CDFTIchi.h"
#include "HilbertTransformR2CDFTI.h"
#include "Functionschi.h"
#include "RootFinder.h"
#include "MultirootFinder.h"
#include "ConvolutionchiBubrhorho.h"
#include "ConvolutionchiBubWW.h"



int main(int argc, char *argv[]) {

    clock_t tic = clock();

    struct timeval tv1;
    gettimeofday(&tv1, NULL);

    int TIME_LIMIT = 60 * 60 * 24; // Time limit in seconds

    // ON = 1 and OFF = 0 defined in C_extension.h
    int WRITE = ON, //  write to file calcuables at EOF
        CHECKSUM = ON, // toggles checksum on and off
        SETUP = OFF, //  toggles on and off print outs from setup section
        FERMI_ENERGY = ON, // show root finder print outs for Fermi energy
        K_FERMI = ON, // show root finder print outs for Fermi wavevector
        INITIAL = OFF, // shows root finder print outs for computation of g0, i.e., initiaization section 
        IMPORT = ON, // imports rho_g_initial from file
        EXPORT = ON, // exports rho_g_final to file
        ITER = OFF, // exports rho_g_loc and rho_bigG_loc to file for each iteration 
        SELFENERGY = OFF, // prints self energy checksum to screen
        CHEM = ON, // prints chem potential to screen
        TUKEY = OFF, // toggles Tukey window on and off
        U0SUMRULE = ON, // use u0 sum rules rather than particle sum rules
        DENSITYAPPROX = OFF; // uses nd in place of ng in mu_true formula

    if (CHECKSUM == ON) { printf("CHECKSUM ON\n"); } else { printf("CHECKSUM OFF\n"); }
    if (WRITE == ON) { printf("WRITE ON\n"); } else { printf("WRITE OFF\n"); }
    if (SETUP == ON) { printf("SETUP ON\n"); } else { printf("SETUP OFF\n"); }
    if (FERMI_ENERGY == ON) { printf("FERMI_ENERGY ON\n");} else { printf("FERMI_ENERGY OFF\n"); }
    if (K_FERMI == ON) { printf("K_FERMI ON\n"); } else { printf("K_FERMI OFF\n"); }
    if (INITIAL == ON) { printf("INITIAL ON\n"); } else { printf("INITIAL OFF\n"); }
    if (IMPORT == ON) { printf("IMPORT ON\n");} else { printf("IMPORT OFF\n");}
    if (EXPORT == ON) { printf("EXPORT ON\n");} else { printf("EXPORT OFF\n");}
    if (ITER == ON) { printf("ITERATION ON\n"); } else  { printf("ITERATION OFF\n");} 
    if (SELFENERGY == ON) { printf("SELFENERGY ON\n");} else { printf("SELFENERGY OFF\n");} 
    if (CHEM == ON) { printf("CHEM ON\n"); } else { printf("CHEM OFF\n"); } 
    if (TUKEY == ON) { printf("TUKEY ON\n");} else { printf("TUKEY OFF\n");}
    if (U0SUMRULE == ON) { printf("U0SUMRULE ON\n"); } else { printf("U0SUMRULE OFF\n"); }
    if (DENSITYAPPROX == ON) { printf("DENSITYAPPROX ON\n"); } else { printf("DENSITYAPPROX OFF\n"); }




    printf("Executable name: %s\n", argv[0]);
    printf("This is %s from %s(), line %d\n", __FILE__, __func__, __LINE__);

//    srand(time(NULL));
    int status = 0;
    char name[250],
         file[250],
         directory[250];

    sprintf(directory, "Data/");


// set_up: {{{


    size_t i, j, k, n, m;

    int d = atoi(argv[1]),
        Nk = atoi(argv[2]), // L = Nk * a0 is the size of the system
        Nomega = pow(2, d);

     size_t N[3] = {Nk, Nk, Nomega}, // is the length of each dimension
            padFT[3] = {0, 0, 3 * Nomega / 2}, // is the padding for each dimension for the convolution function
            padHT[3] = {0, 0, 4 * Nomega}, // is the padding for each dimension for the Hilbert transform
            NpadFT[3],
            NpadHT[3];

    for (n = 0; n < 3; n++) {
        NpadFT[n] = N[n] + 2 * padFT[n];
        NpadHT[n] = N[n] + 2 * padHT[n];
    }

    size_t Ns = N[0] * N[1], // Ns is the total sites in the system
           N3 = N[0] * N[1] * N[2], // N3 is the total length the 3-dimensional array for omega-k-space data
           N3padFT = NpadFT[0] * NpadFT[1] * NpadFT[2], /* N3padFT is the length of the padded
                                                           array for the Fast Fourier Transform */
           N3padHT = NpadHT[0] * NpadHT[1] * NpadHT[2]; /* N3padHT is the length of the padded
                                                           array for the Hilbert transform */


    double lambda = 1, // interpolation parameter
           nd = atof(argv[3]), // density
           tau = atof(argv[4]), // natural Temperature tau = kb * T in units of t
           t = 1.0, // 1st nearest neighbor hopping parameter
           tp = atof(argv[5]) * t, // 2nd next nearest neighbor hopping parameter
           tpp = atof(argv[6]) * t, // 3rd nearest neighbors hopping param
           J = atof(argv[7]) * t, // exchange parameter
           gamma_SE = nd / 2; // gamma self-energy


    double mup_initial = atof(argv[8]), // initial guess for the chemical potential 
           u0_initial = atof(argv[9]); // initial guess for the second chemical potential

    int neta = atoi(argv[10]); // width of the Lorentzian: GAMMA = pow(2, neta) * DELTAomega

    double omega_domain = atof(argv[11]) * t,
           omegai = -omega_domain,
           omegaf = omega_domain,
           DELTAomega = (omegaf - omegai) / Nomega,
           ki = 0,
           kf = 2 * M_PI,
           DELTAk = (kf - ki) / Nk;


    // tukey window bounds
    double OMEGAcminus = 6 * t, // inner bound
           OMEGAcplus = 8 * t; // outer bound
    printf("u0SumRule - Summary:\n");
    printf("Tukey Window: \n");
    printf("OMEGAcminus = %lf\n", OMEGAcminus);
    printf("OMEGAcplus = %lf\n", OMEGAcplus);
    printf("omega_domain =%lf\n", omega_domain);
    printf("omega_max = %lf, omega_min = %lf\n", omegaf, omegai);
    printf("Nomega = %d\n", Nomega);
    printf("Nk (L) = %d\n", Nk);
    printf("lambda = %lf\n", lambda);
    printf("Density = %lf\n", nd);
    printf("tau = %lf\n", tau);
    printf("t = %lf\n", t);
    printf("tp = %lf\n", tp);
    printf("tpp = %lf\n", tpp);
    printf("J = %lf\n", J);
    printf("gamma_SE = %lf\n", gamma_SE);
    printf("mup_initial = %lf\n", mup_initial);
    printf("u0_initial = %lf\n", u0_initial);


    // computing the omega and k arrays
    double *k_array     = calloc(Nk, sizeof(double)),
           *omega_array = calloc(Nomega, sizeof(double));
    for (i = 0; i < Nk; i++) {
         k_array[i] = ki + DELTAk * i;
    }
    for (i = 0; i < Nomega; i++) {
         omega_array[i] = omegai + DELTAomega * i;
    }


    if (CHECKSUM) {
        check_sum_double_array_1d("karray", k_array, Nk);
        check_sum_double_array_1d("omega_array", omega_array, Nomega);
    }


    double *Fermi_array    = calloc(Nomega, sizeof(double)),
           *Fermibar_array = calloc(Nomega, sizeof(double));

    for (k = 0; k < Nomega; k++) {
        Fermi_array[k] = Fermi_function(omega_array[k], tau);
        Fermibar_array[k] = 1 - Fermi_array[k];
    }
    if (CHECKSUM && SETUP) {
        check_sum_double_array_1d("Fermi_array", Fermi_array, Nomega);
        check_sum_double_array_1d("Fermibar_array", Fermibar_array, Nomega);
    }


    double *en = calloc(Ns, sizeof(double)), // band energy
           *xi = calloc(Ns, sizeof(double)), // effective band energy
           *c  = calloc(Nk, sizeof(double)), // trig arrays
           *s  = calloc(Nk, sizeof(double));


    for (i = 0; i < Nk; i++) {
        c[i] = cos(k_array[i]);
        s[i] = sin(k_array[i]);
    }
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            en[m] = epsilon_k(t, tp, tpp, k_array[i], k_array[j]);
        }
    }


    if (CHECKSUM && SETUP) {
        check_sum_double_array_2d("en", en, N);
        check_sum_double_array_2d("xi", xi, N);
        check_sum_double_array_1d("cx", c, Nk);
        check_sum_double_array_1d("sx", s, Nk);
    }


// }}} set_up


// Fermi energy: {{{


    // Find Root
    struct FindFermiEnergyRootFunctionParams find_Fermi_energy_root_function_params = {N,
                                                                                       nd,
                                                                                       en};

    gsl_function F;
    F.function = &find_Fermi_energy_root_function;
    F.params = &find_Fermi_energy_root_function_params;


    // {lower bound, upper bound}
    struct FindRootParams find_root_params = {-5, 5, FERMI_ENERGY};
    double Fermi_energy = find_root(&F, &find_root_params);
    if (FERMI_ENERGY) {
        printf("Fermi_energy = %-5.15lf\n", Fermi_energy);
    }
    // check density
    double density = find_density(Fermi_energy, en, N); 
    if (FERMI_ENERGY) {
        printf("check density = %lf\n", density);
    }


// }}}


// k_fermi: {{{


    // Find Root of Fermi_energy - espilon_k(k,k) = 0
    struct k_FermiRootFunctionParams k_Fermi_root_function_params = {t,
                                                                     tp,
                                                                     tpp,
                                                                     Fermi_energy};

    F.function = &k_Fermi_nodal_root_function;
    F.params = &k_Fermi_root_function_params;

    // {lower bound, upper bound}
    struct FindRootParams find_root_k_Fermi_params = {0.0, M_PI, K_FERMI};
    double k_Fermi_nodal = find_root(&F, &find_root_k_Fermi_params);
    printf("k_Fermi_nodal = %-5.15lf\n", k_Fermi_nodal);


    F.function = &k_Fermi_antinodal_root_function;
    F.params = &k_Fermi_root_function_params;

    double k_Fermi_antinodal = find_root(&F, &find_root_k_Fermi_params);
    printf("k_Fermi_antinodal = %-5.15lf\n", k_Fermi_antinodal);


// }}}


// initialize: {{{

    double *rho_g_initial = calloc(N3, sizeof(double));
    double mup_root = mup_initial;


    // Find Root
    struct RootFunctionParams root_function_params = {nd,
                                                      DELTAomega,
                                                      neta,
                                                      omega_array,
                                                      Fermi_array,
                                                      en,
                                                      rho_g_initial,
                                                      N};

//    gsl_function F;
    F.function = &root_function;
    F.params = &root_function_params;

    // {lower bound, upper bound}
    struct FindRootParams find_root_initial_params = {-5.0, 5.0, INITIAL};
    mup_root = find_root(&F, &find_root_initial_params);
    printf("root_finder mup_initial = %-5.15lf\n", mup_root);


    for (i = 0; i < Nk; i++) {
        for (j = 0; j < Nk; j++) {
            m = (j + Nk * i);
            for (k = 0; k < Nomega; k++) {
                n = k + Nomega * m;
                rho_g_initial[n] = Lorentz_function(omega_array[k], en[m] - mup_initial, pow(2, neta) * DELTAomega);
            }
        }
    }


    double sum = 0;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                sum += rho_g_initial[n] * Fermi_array[k];
            }
        }
    }
    sum *= DELTAomega / Ns;
    printf("sum gf = %lf\n", sum);

    if (CHECKSUM ) {
        check_sum_double_array_1d("rho_g_initial", rho_g_initial, N3);
    }

    if (INITIAL) {

/*
{
        strcpy(name, directory);
        sprintf(file, "Spectrum/rho_g_nodal_initial"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf}.dat"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial);
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (i + N[1] * i);
                fprintf(fdatatmp, "%+5.6lf\t %+5.6lf\t %+5.6lf\n", k_array[i], omega_array[n], rho_g_initial[n]);
            }
        }
        fclose(fdatatmp);
}
*/

        }


// }}} initialize


// Import: {{{
// we could add import arguments here
    if (IMPORT) {
{       
        double importdata;
        FILE *importfile;
        strcpy(name, directory);
        sprintf(file, "Import/rho_g_export{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, dom%3.3lf}.dat"
        , d, Nk, nd, 0.057, -0.2, tpp, J, omega_domain);
        strcat(name, file);
        printf("%s\n", name);
        importfile = fopen(name, "r");
 //       pfile = fopen("/Users/Michael/Dropbox/Exports/im_deltapsiQP.dat", "r");
        if (importfile == NULL) {printf("importfile == NULL\n"); exit(0);} else {printf("importfile != NULL\n");}
        k = 0;
        while (EOF != fscanf(importfile, "%lf", &importdata))
        {
            rho_g_initial[k] = importdata;
//            printf("rho_psiQP[%zu]=%lf\n", k, rho_psiQP[k]);
//            printf("k=%zu\n", k);
            k++;
        }
        fclose(importfile);
}
    }


// }}}


// increment: {{{


    int iter = 1,
        itermax = 99,
        CONVERGENCE = 0;
    double mup_temp = mup_initial,
           u0_temp = u0_initial,
           mup,
           u0;
    double error = pow(10, -5);

    double *rho_g_temp  = calloc(N3, sizeof(double)),
           *rho_g       = calloc(N3, sizeof(double)),
           *rho_psi0    = calloc(N3, sizeof(double)),
           *rho_psi1    = calloc(N3, sizeof(double)),
           *rho_chi0    = calloc(N3, sizeof(double)),
           *rho_chi1    = calloc(N3, sizeof(double)),
           *rho_chi2    = calloc(N3, sizeof(double)),
           *re_psi0     = calloc(N3, sizeof(double)),
           *re_psi1     = calloc(N3, sizeof(double)),
           *re_chi0     = calloc(N3, sizeof(double)),
           *re_chi1     = calloc(N3, sizeof(double)),
           *re_chi2     = calloc(N3, sizeof(double));

    copy_double_array_3d(rho_g_initial, rho_g_temp, N);

    if (CHECKSUM && INITIAL) {
        check_sum_double_array_1d("rho_g_initial", rho_g_initial, N3);
        check_sum_double_array_1d("rho_g_temp", rho_g_temp, N3);
    }

    free(rho_g_initial); rho_g_initial = NULL; check_pointer(rho_g_initial, __FILE__, __func__, __LINE__);


    double gerror = 1;


    do {
        printf("do iter = %d\n", iter);


// self_eneriges: {{{


        struct SelfEnergyParams self_energy_params = {N,
                                                      padFT,
                                                      NpadFT,
                                                      lambda,
                                                      J,
                                                      DELTAomega,
                                                      Fermi_array,
                                                      Fermibar_array,
                                                      en,
                                                      c,
                                                      s,
                                                      rho_g_temp,
                                                      rho_psi0,
                                                      rho_psi1,
                                                      rho_chi0,
                                                      rho_chi1,
                                                      rho_chi2};


        status = self_energy(&self_energy_params);


        if (CHECKSUM && SELFENERGY) {
            check_sum_double_array_1d("rho_psi0", rho_psi0, N3);
            check_sum_double_array_1d("rho_psi1", rho_psi1, N3);
            check_sum_double_array_1d("rho_chi0", rho_chi0, N3);
            check_sum_double_array_1d("rho_chi1", rho_chi1, N3);
            check_sum_double_array_1d("rho_chi2", rho_chi2, N3);
        }

        timer("Self energy", tic, tv1, TIME_LIMIT);

        struct HilbertParams hilbert_params = {N,
                                               padHT,
                                               NpadHT};
        status = Hilbert_transform(rho_psi0, re_psi0, &hilbert_params);
        status = Hilbert_transform(rho_psi1, re_psi1, &hilbert_params);
        status = Hilbert_transform(rho_chi0, re_chi0, &hilbert_params);
        status = Hilbert_transform(rho_chi1, re_chi1, &hilbert_params);
        status = Hilbert_transform(rho_chi2, re_chi2, &hilbert_params);


        timer("Hilbert", tic, tv1, TIME_LIMIT);


        if (CHECKSUM && SELFENERGY) {
            check_sum_double_array_1d("re_psi0", re_psi0, N3);
            check_sum_double_array_1d("re_psi1", re_psi1, N3);
            check_sum_double_array_1d("re_chi0", re_chi0, N3);
            check_sum_double_array_1d("re_chi1", re_chi1, N3);
            check_sum_double_array_1d("re_chi2", re_chi2, N3);
        }


// }}} selfenergies


        double gfe = 0,
               js  = 0,
               jc  = 0;
        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = (j + N[1] * i);
                for (k = 0; k < N[2]; k++) {
                    n = k + N[2] * m;
                    gfe += rho_g_temp[n] * Fermi_array[k] * en[m];
                    js  += rho_g_temp[n] * Fermi_array[k] * s[i];
                    jc  += rho_g_temp[n] * Fermi_array[k] * c[i];
                }
            }
        }
        gfe *= DELTAomega / Ns;
        js  *= DELTAomega / Ns;
        jc  *= DELTAomega / Ns;

        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = (j + N[1] * i);
                xi[m] = en[m] * (1 - lambda * gamma_SE) 
                      - lambda * J * ((c[i] + c[j]) * jc + (s[i] + s[j]) * js);
            }
        }

        if (ITER) {
            printf("gfe = %+5.15lf, js = %+5.15lf, jc = %+5.15lf\n", gfe, js, jc);
        }

// multiroot_finder: {{{


        const size_t num_of_roots = 2;
        struct MultirootFunctionParams multirootFunctionParams = {N,
                                                                  lambda,
                                                                  nd,
                                                                  J,
                                                                  gamma_SE,
                                                                  DELTAomega,
                                                                  gfe,
                                                                  omega_array,
                                                                  Fermi_array,
                                                                  en,
                                                                  xi,
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
                                                                  U0SUMRULE,
                                                                  DENSITYAPPROX};

        gsl_multiroot_function f = {&multiroot_function, num_of_roots, &multirootFunctionParams};

        gsl_vector *x = gsl_vector_alloc(num_of_roots),
                   *y = gsl_vector_alloc(num_of_roots);
        gsl_vector_set(x, 0, mup_temp);
        gsl_vector_set(x, 1, u0_temp);
        if (CHEM) {
            printf("mup_temp = %+5.15lf, u0_temp = %+5.15lf\n", x->data[0], x->data[1]);
        }
        multiroot_finder(x, &f, y);
        mup = gsl_vector_get(y, 0);
        u0  = gsl_vector_get(y, 1);
        if (CHEM) {
            printf("mup = %+5.15lf, u0 = %+5.15lf\n", y->data[0], y->data[1]);
        }
        gsl_vector_free(x);
        gsl_vector_free(y);

        mup_temp = mup;
        u0_temp  = u0;



// }}}  multiroot_finder


// rho_g_tukey: {{{


        double *rho_g_tukey = calloc(N3, sizeof(double));
        double rho_psi,
               rho_chi,
               rho_phi,
               re_psi,
               re_chi,
               re_phi,
               rr,
               X = (-u0/2),
               denom,
               re_g;
        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = (j + N[1] * i);
                for (k = 0; k < N[2]; k++) {
                    n = k + N[2] * m;
                    rho_psi  = rho_psi0[n] + rho_psi1[n] * X;
                    re_psi   = re_psi0[n]  + re_psi1[n]  * X;
                    rho_chi  = rho_chi0[n] + rho_chi1[n] * X + rho_chi2[n] * X * X;
                    re_chi   = re_chi0[n]  + re_chi1[n]  * X + re_chi2[n]  * X * X;
                    rho_phi  = rho_chi + (en[m] + X) * rho_psi; 
                    re_phi   = re_chi  + (en[m] + X) * re_psi;
                    rr       = omega_array[k] + mup - xi[m] - re_phi;
                    denom    = pow(rr, 2) + pow(M_PI * rho_phi, 2);
                    rho_g[n] = rho_phi / denom;
                }
            }
        }

        if (CHECKSUM && OFF) {
            check_sum_double_array_1d("rho_g", rho_g, N3);
        }

        tukey_max_array_3d_function(rho_g, rho_g_tukey, omega_array, DELTAomega, OMEGAcminus, OMEGAcplus, TUKEY, N);

        if (CHECKSUM) {
            check_sum_double_array_1d("rho_g", rho_g, N3);
            check_sum_double_array_1d("rho_g_tukey", rho_g_tukey, N3);
        }


        if (WRITE && ITER) {

{
    double *local = calloc(N3, sizeof(double));
    localize_function(rho_g, local, N);
    strcpy(name, directory);
    sprintf(file, "Spectrum/rho_g_loc"
    "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, iter%02d}.dat"
    , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, iter);
    strcat(name, file);
    printf("%s\n", name);
    write_to_file_double_array_1d(name, local, N[2]);
    free(local);
}

    }

        printf("gerror = %-5.15lf\n", gerror);
        double gerrorcalc = calculate_gerror(rho_g_temp, rho_g_tukey, N);
        printf("gerrorcalc = %-5.15lf\n", gerrorcalc);
        gerror = gerrorcalc;
        
        double weight = 1 / pow(3, 0.5);
        weighted_average_array_3d(rho_g_temp, rho_g_tukey, weight, N);

        free(rho_g_tukey); rho_g_tukey = NULL; check_pointer(rho_g_tukey, __FILE__, __func__, __LINE__);



// }}}





 //       if ((gerrorcalc > gerror) && (iter > 0)) {
//            printf("goto\n");
//            goto label1;
//        }


         if (gerror < error) {
            printf("break\n");
            CONVERGENCE = 1;
            break;
        }



        timer("Iter", tic, tv1, TIME_LIMIT);

//       label1:
        iter++;
    } while (iter <= itermax);


// }}} increment


// compute_data: {{{


// rho_bigG: {{{


    double *re_g         = calloc(N3, sizeof(double)), 
           *rho_bigG     = calloc(N3, sizeof(double)),
           *re_bigG      = calloc(N3, sizeof(double));
    
    double *rho_bigG_loc = calloc(N[2], sizeof(double)),
           *re_bigG_loc  = calloc(N[2], sizeof(double)),
           *rho_g_loc    = calloc(N[2], sizeof(double)),
           *re_g_loc     = calloc(N[2], sizeof(double));


    double rho_psi,
           rho_chi,
           rho_phi,
           re_psi,
           re_chi,
           re_phi,
           rr,
           X = -u0/2,
           denom;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                rho_psi     = rho_psi0[n] + rho_psi1[n] * X;
                re_psi      = re_psi0[n]  + re_psi1[n]  * X;
                rho_chi     = rho_chi0[n] + rho_chi1[n] * X + rho_chi2[n] * X * X;
                re_chi      = re_chi0[n]  + re_chi1[n]  * X + re_chi2[n]  * X * X;
                rho_phi     = rho_chi     + (en[m] + X) * rho_psi; 
                re_phi      = re_chi      + (en[m] + X) * re_psi;
                rr          = omega_array[k] + mup - xi[m] - re_phi;
                denom       = pow(rr, 2) + pow(M_PI * rho_phi, 2);
                re_g[n]     = rr / denom;
                rho_g[n]    = rho_phi / denom;
                re_bigG[n]  = re_g[n] * (1 - lambda * gamma_SE + re_psi) - pow(M_PI, 2) * rho_g[n] * rho_psi;
                rho_bigG[n] = rho_g[n] * (1 - lambda * gamma_SE + re_psi) + re_g[n] * rho_psi;
            }
        }
    }


    free(re_chi0); re_chi0 = NULL; check_pointer(re_chi0, __FILE__, __func__, __LINE__);
    free(re_chi1); re_chi1 = NULL; check_pointer(re_chi1, __FILE__, __func__, __LINE__);
    free(re_chi2); re_chi2 = NULL; check_pointer(re_chi2, __FILE__, __func__, __LINE__);

    free(re_psi0); re_psi0 = NULL; check_pointer(re_psi0, __FILE__, __func__, __LINE__);
    free(re_psi1); re_psi1 = NULL; check_pointer(re_psi1, __FILE__, __func__, __LINE__);

    free(rho_chi0); rho_chi0 = NULL; check_pointer(re_chi0, __FILE__, __func__, __LINE__);
    free(rho_chi1); rho_chi1 = NULL; check_pointer(re_chi1, __FILE__, __func__, __LINE__);
    free(rho_chi2); rho_chi2 = NULL; check_pointer(re_chi2, __FILE__, __func__, __LINE__);


    free(rho_psi0); rho_psi0 = NULL; check_pointer(re_psi0, __FILE__, __func__, __LINE__);
    free(rho_psi1); rho_psi1 = NULL; check_pointer(re_psi1, __FILE__, __func__, __LINE__);

    timer("rho_bigG", tic, tv1, TIME_LIMIT);


    localize_function(rho_g, rho_g_loc, N);
    localize_function(rho_bigG, rho_bigG_loc, N);
    localize_function(re_g, re_g_loc, N);
    localize_function(re_bigG, re_bigG_loc, N);

    timer("localize", tic, tv1, TIME_LIMIT);


    double ng = 2 * density_function(rho_g, Fermi_array, DELTAomega, N),
           nbigG = 2 * density_function(rho_bigG, Fermi_array, DELTAomega, N),
           pseudo_nbigG = 2 * pseudo_Fermi_surface(rho_bigG, k_array, omega_array, DELTAomega, tau, N);

    if (CHECKSUM) {
        check_sum_double_array_1d("rho_g", rho_g, N3);
        check_sum_double_array_1d("rho_bigG", rho_bigG, N3);
    }


    if(EXPORT) {
{
        strcpy(name, directory);
        sprintf(file, "Import/rho_g_export"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, dom%3.3lf}.dat"
            , d, Nk, nd, tau, tp, tpp, J, omega_domain);
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < N3; i++) {
            fprintf(fdatatmp, "%+3.6lf\n", rho_g[i]);
        }
        fclose(fdatatmp);
}
    }


    if (WRITE) {

{
    strcpy(name, directory);
    sprintf(file, "Spectrum/g_loc_cmplx"
    "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, dom%2.0lf}%s"
    , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
    strcat(name, file);
    printf("%s\n", name);
    FILE *fdatatmp = fopen(name, "w");
    for (k = 0; k < Nomega; k++) {
        fprintf(fdatatmp, "%+5.6lf\t%5.6lf\t%+5.6lf\n", omega_array[k], re_g_loc[k], -M_PI * rho_g_loc[k]);
    }
    fclose(fdatatmp);
}


{
    strcpy(name, directory);
    sprintf(file, "Spectrum/bigG_loc_cmplx"
    "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%02d, dom%2.0lf}%s"
    , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
    strcat(name, file);
    printf("%s\n", name);
    FILE *fdatatmp = fopen(name, "w");
    for (k = 0; k < Nomega; k++) {
        fprintf(fdatatmp, "%+5.6lf\t%5.6lf\t%+5.6lf\n", omega_array[k], re_bigG_loc[k] , -M_PI * rho_bigG_loc[k]);
    }
    fclose(fdatatmp);
}


{
        strcpy(name, directory);
        sprintf(file, "Spectrum/g_nodal_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (i + N[1] * i);
                fprintf(fdatatmp, "%+5.6lf\t%5.6lf\t%+5.6lf\t%+5.6lf\n", 
                                   k_array[i], omega_array[k], re_g[n] , -M_PI * rho_g[n]);
            }
        }
        fclose(fdatatmp);
}


{
        strcpy(name, directory);
        sprintf(file, "Spectrum/bigG_nodal_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (i + N[1] * i);
                fprintf(fdatatmp, "%+5.6lf\t%5.6lf\t%+5.6lf\t%+5.6lf\n", 
                                   k_array[i], omega_array[k], re_bigG[n] , -M_PI * rho_bigG[n]);
            }
        }
        fclose(fdatatmp);
}


    }

        timer("write rho_bigG", tic, tv1, TIME_LIMIT);



// }}}


// kinetic energy: {{{


double Tx = 0, 
       Ty = 0, 
       Txy = 0,
       cxave = 0, 
       cyave = 0, 
       cxyave = 0, 
       qx, 
       qy;
    for (i = 0; i < N[0]; i++) {
        qx = k_array[i];
        for (j = 0; j < N[1]; j++) {
            qy = k_array[j];
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                Tx += cos(qx) * rho_bigG[n] * Fermi_array[k];
                Ty += cos(qy) * rho_bigG[n] * Fermi_array[k];
                Txy += cos(qx) * cos(qy) * rho_bigG[n] * Fermi_array[k];
                cxave += cos(qx) * rho_bigG[n] * Fermi_array[k];
                cyave += cos(qy) * rho_bigG[n] * Fermi_array[k];
                cxyave += cos(qx) * cos(qy) * rho_bigG[n] * Fermi_array[k];
            }
        }
    }
    Tx  *= -2 * t * DELTAomega / Ns;
    Ty  *= -2 * t * DELTAomega / Ns;
    Txy *= -4 * tp * DELTAomega / Ns;
    cxave *= DELTAomega / Ns;
    cyave *= DELTAomega / Ns;
    cxyave *= DELTAomega / Ns;

    timer("kinetic energy", tic, tv1, TIME_LIMIT);



// }}}


// rho_SIGMA_Dyson: {{{


    double gfe = 0,
           bigGfe = 0; // total kinetic energy
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                gfe += rho_g[n] * Fermi_array[k] * en[m];
                bigGfe += rho_bigG[n] * Fermi_array[k] * en[m];
            }
        }
    }
    gfe *= DELTAomega / Ns;
    bigGfe *= DELTAomega / Ns;


    double J0 = J_k(J, 0, 0),
           aG = 1 - lambda * gamma_SE,
           chi0c = -gfe + (u0/2) * (ng/2),
           mu = mup + u0 / 2 - lambda * (J0/2) * gamma_SE - u0/2 * aG + lambda * chi0c; 
    
//    mu1 = mu;
    printf("mu = %lf \n", mu);



    double *rho_SIGMA = calloc(N3, sizeof(double)),
           *re_SIGMA  = calloc(N3, sizeof(double));

    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * (j + N[1] * i);
                rr = omega_array[k] + mu - en[m];
                denom = pow(re_bigG[n], 2) + pow(M_PI * rho_bigG[n], 2);
                rho_SIGMA[n] = rho_bigG[n] / denom;
                re_SIGMA[n] = rr - re_bigG[n] / denom;
            }
        }
    }

    timer("Dyson Self energy", tic, tv1, TIME_LIMIT);



// }}}


// Z: {{{



// computes quasi-particle peak height, location and width
    double *peak_value = calloc(Ns, sizeof(double)),
           *GAMMA      = calloc(Ns, sizeof(double));
    size_t *peak_index = calloc(Ns, sizeof(double));
    find_GAMMA(rho_bigG, peak_value, peak_index, GAMMA, DELTAomega, N);


    free(GAMMA);
    free(peak_value);
    free(peak_index);


// compute quasi-particle weight
    size_t peak_size = 11;
    double eta = DELTAomega; // small parameter for peak
    double *theta  = calloc(peak_size, sizeof(double)),
           *Z = calloc(Ns, sizeof(double));
    size_t l;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            for (k = 0; k < peak_size; k++) {
                n = k - 4 + N[2] / 2 + N[2] * m;
                l = k - 6 + N[2] / 2 + N[2] * m;
                denom = 1.0 - (re_SIGMA[n] - re_SIGMA[l]) / (2 * eta);
                theta[k] = 1.0 / denom;
                Z[m] += theta[k];
           }
           Z[m] /= peak_size;
        }
    }



    if (WRITE) {
{
        strcpy(name, directory);
        sprintf(file, "Spectrum/Z"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        double kx, ky;
        for (i = 0; i < N[0]; i++) {
            for (j = 0; j < N[1]; j++) {
                m = j + N[1] * i;
                kx = k_array[i];
                ky = k_array[j];
                fprintf(fdatatmp, "%+5.6lf %+5.6lf %+5.6lf %+5.6lf\n", kx, ky, k_Fermi_nodal, Z[m]);
            }
        }
        fclose(fdatatmp);
}
    }


    free(Z); 
    free(theta);

    timer("Z", tic, tv1, TIME_LIMIT);


// }}}


// Polairization: {{{




    double complex *chiBubWW_cmplx = calloc(N3, sizeof(double complex)),
                   *chiBubrhorho_cmplx = calloc(N3, sizeof(double complex));



// chiBubWW: {{{


    double *rho_chiBubWW = calloc(N3, sizeof(double)),
           *re_chiBubWW = calloc(N3, sizeof(double));


    struct ChiBubWWParams chiBubWW_params = {N,
                                             padFT,
                                             NpadFT,
                                             DELTAomega,
                                             en,
                                             Fermi_array,
                                             rho_bigG,
                                             rho_chiBubWW};

    chiBubWW_function(&chiBubWW_params); 


    struct HilbertParams hilbert_params_chiBubWW = {N,
                                                    padHT,
                                                    NpadHT};


    status = Hilbert_transform(rho_chiBubWW, re_chiBubWW, &hilbert_params_chiBubWW);



    if (WRITE) {
{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubWW_GAMMAtoX_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (i + N[1] * i);
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf %5.6lf\n", k_array[i], omega_array[k], 
                        re_chiBubWW[n], -M_PI * rho_chiBubWW[n]);
            }
        }
        fclose(fdatatmp);
}
    }


    if (WRITE) {

{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubWW_GAMMAtoM_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (N[1] * i);
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf %5.6lf\n", k_array[i], omega_array[k], 
                        re_chiBubWW[n], -M_PI * rho_chiBubWW[n]);
            }
        }
        fclose(fdatatmp);
}
    }


    if (WRITE) {

{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubWW_MtoX_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (j = 0; j < Nk / 2 + 1; j++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (j + N[1] * (Nk / 2));
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf %5.6lf\n", k_array[j], omega_array[k], 
                        re_chiBubWW[n], -M_PI * rho_chiBubWW[n]);
            }
        }
        fclose(fdatatmp);
}
    }




    if (WRITE) {

{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubWWq0"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (j = 0; j < Nk / 2 + 1; j++) {
                m = (j + N[1] * i);
                l = N[2]/ 2 + N[2] * m;
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf\n", k_array[i], k_array[j], re_chiBubWW[l]);
            }
        }
        fclose(fdatatmp);
}
    }



    for (i = 0; i < Nk; i++) {
        for (j = 0; j < Nk; j++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (j + N[1] * i);
                chiBubWW_cmplx[n] = re_chiBubWW[n] - M_PI * rho_chiBubWW[n] * I;
            }
        }
    }




    free(rho_chiBubWW);
    free(re_chiBubWW);

    timer("chiBubWW", tic, tv1, TIME_LIMIT);



//    }}}


// chiBubrhorho: {{{


    double *rho_chiBubrhorho = calloc(N3, sizeof(double)),
           *re_chiBubrhorho = calloc(N3, sizeof(double));


    struct ChiBubrhorhoParams chiBubrhorho_params = {N,
                                                     padFT,
                                                     NpadFT,
                                                     DELTAomega,
                                                     Fermi_array,
                                                     rho_bigG,
                                                     rho_chiBubrhorho};

    chiBubrhorho_function(&chiBubrhorho_params); 

    struct HilbertParams hilbert_params_chiBubrhorho = {N,
                                                        padHT,
                                                        NpadHT};

    status = Hilbert_transform(rho_chiBubrhorho, re_chiBubrhorho, &hilbert_params_chiBubrhorho);

    for (i = 0; i < Nk; i++) {
        for (j = 0; j < Nk; j++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (j + N[1] * i);
                chiBubrhorho_cmplx[n] = re_chiBubrhorho[n] - M_PI * rho_chiBubrhorho[n] * I;
            }
        }
    }


    if (WRITE) {


{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubrhorho_GAMMAtoX_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (i + N[1] * i);
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf %5.6lf\n", k_array[i], omega_array[k], 
                        re_chiBubrhorho[n], -M_PI * rho_chiBubrhorho[n]);
            }
        }
        fclose(fdatatmp);
}


{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubrhorho_GAMMAtoM_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (N[1] * i);
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf %5.6lf\n", k_array[i], omega_array[k], 
                        re_chiBubrhorho[n], -M_PI * rho_chiBubrhorho[n]);
            }
        }
        fclose(fdatatmp);
}


{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubrhorho_MtoX_cmplx"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (j = 0; j < Nk / 2 + 1; j++) {
            for (k = 0; k < Nomega; k++) {
                n = k + N[2] * (j + N[1] * (Nk / 2));
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf %5.6lf\n", k_array[j], omega_array[k], 
                        re_chiBubrhorho[n], -M_PI * rho_chiBubrhorho[n]);
            }
        }
        fclose(fdatatmp);
}






{
        strcpy(name, directory);
        sprintf(file, "Polarization/chiBubrhorhoq0"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, " 
            "J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%03d, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial,
            neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        for (i = 0; i < Nk / 2 + 1; i++) {
            for (j = 0; j < Nk / 2 + 1; j++) {
                m = (j + N[1] * i);
                l = N[2]/ 2 + N[2] * m;
                fprintf(fdatatmp, "%5.6lf %5.6lf %+5.6lf\n", k_array[i], k_array[j], re_chiBubrhorho[l]);
            }
        }
        fclose(fdatatmp);
}


    }

    double complex chiBubrhorho00 = chiBubrhorho_cmplx[Nomega/2]; 


    free(rho_chiBubrhorho);
    free(re_chiBubrhorho);


    timer("chiBubrhorho", tic, tv1, TIME_LIMIT);



// }}}


free(chiBubWW_cmplx);
free(chiBubrhorho_cmplx);



    timer("chiBubWW", tic, tv1, TIME_LIMIT);



// }}}


// }}}


// sigma: {{{


    double *phik      = calloc(Ns, sizeof(double)),
           *phik2     = calloc(Ns, sizeof(double));
    double sigmaxx = 0,
           sigmayy = 0,
           sigmaxy = 0,
           sigmac  = 0,
           sigmacn = 0,
           kx,
           ky,
           vx,
           vy,
           vxx,
           vyy,
           vxy,
           ux = 1.0,
           uy = 0.0,
           eta1,
           eta2,
           eta3,
           df,
           tmp;
    double a1g = 0,
           b1g = 0,
           b2g = 0;
    for (i = 0; i < N[0]; i++) {
        for (j = 0; j < N[1]; j++) {
            m = (j + N[1] * i);
            kx = k_array[i];
            ky = k_array[j];
            vx  = derivative_epsilon_kx(t, tp, tpp, kx, ky);
            vy  = derivative_epsilon_ky(t, tp, tpp, kx, ky);
            vxx = derivative_epsilon_kxx(t, tp, tpp, kx, ky);
            vyy = derivative_epsilon_kyy(t, tp, tpp, kx, ky);
            vxy = derivative_epsilon_kxy(t, tp, tpp, kx, ky);
            eta1 = pow(vx * ux + vy * uy, 2);
            eta2 = pow(vx, 2) * vyy - vx * vy * vxy;
            eta3 = pow(vy, 2);
            for (k = 0; k < N[2]; k++) {
                n = k + N[2] * m;
                tmp = rho_bigG[n];
                df = derivative_Fermi_function_over_beta(omega_array[k], tau);
                sigmaxx  += -df * pow(tmp, 2) * eta1;
                sigmayy  += -df * pow(tmp, 2) * eta3;
                sigmaxy  += -df * pow(tmp, 3) * eta2;
                sigmac   += -df * pow(tmp, 2);
                a1g      += -df * pow(tmp, 2) * pow(vxx + vyy, 2) / 4;
                b1g      += -df * pow(tmp, 2) * pow(vxx - vyy, 2) / 4; 
                b2g      += -df * pow(tmp, 2) * pow(vxy, 2);
                phik[m]  += -df * pow(tmp, 2); 
                phik2[m] += -df * pow(tmp, 3);
            }
        }
        phik[m]  *= (4 * tp * pow(2 * M_PI, 2) * DELTAomega / tau);
        phik2[m] *= (DELTAomega / tau);
    }

    for (k = 0; k < N[2]; k++) {
        sigmacn += -derivative_Fermi_function_over_beta(omega_array[k], tau) * pow(rho_bigG_loc[k], 2);
    }


    sigmaxx *= pow(2 * M_PI, 2) * DELTAomega / (tau * N[0] * N[1]);
    sigmayy *= pow(2 * M_PI, 2) * DELTAomega / (tau * N[0] * N[1]);
    sigmaxy *= pow(2 * M_PI, 2) * DELTAomega / (3 * tau * N[0] * N[1]);
    sigmac  *= 4 * M_PI * DELTAomega / (tau * N[0] * N[1]); 
    sigmacn *= 4 * M_PI * DELTAomega / tau;
    a1g     *= pow(2 * M_PI, 2) * DELTAomega / (tau * N[0] * N[1]);
    b1g     *= pow(2 * M_PI, 2) * DELTAomega / (tau * N[0] * N[1]);
    b2g     *= pow(2 * M_PI, 2) * DELTAomega / (tau * N[0] * N[1]);



    if (WRITE) {
{
        strcpy(name, directory);
        sprintf(file, "Spectrum/sigma"
            "{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, dom%2.0lf}%s"
            , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdatatmp = fopen(name, "w");
        fprintf(fdatatmp, "%5.6lf %+5.6lf %+5.6lf %+5.6lf %+5.6lf\n", nd, tau, sigmaxx, sigmayy, sigmaxy);
        fclose(fdatatmp);
 }
    }


free(phik);
free(phik2);

timer("sigma", tic, tv1, TIME_LIMIT);


// }}}


// }}}


    if (WRITE) {
        strcpy(name, directory);
        sprintf(file, "Summary/"
        "Data{Nomega2E%02d, Nk%03d, nd%3.3lf, tau%3.3lf, tp%3.3lf, tpp%3.3lf, J%3.3lf, mup%3.3lf, u0_%3.3lf, neta%02d, dom%2.0lf}%s"
        , d, Nk, nd, tau, tp, tpp, J, mup_initial, u0_initial, neta, omega_domain, CONVERGENCE == 1 ? "SUCCESS.dat" : "FAILED.dat");
        strcat(name, file);
        printf("%s\n", name);
        FILE *fdata = fopen(name, "w");

        if (CHECKSUM == ON) { fprintf(fdata,"CHECKSUM ON\n"); } else { fprintf(fdata,"CHECKSUM OFF\n"); }
        if (WRITE == ON) { fprintf(fdata,"WRITE ON\n"); } else { fprintf(fdata,"WRITE OFF\n"); }
        if (SETUP == ON) { fprintf(fdata,"SETUP ON\n"); } else { fprintf(fdata,"SETUP OFF\n"); }
        if (FERMI_ENERGY == ON) { fprintf(fdata,"FERMI_ENERGY ON\n");} else { fprintf(fdata,"FERMI_ENERGY OFF\n"); }
        if (K_FERMI == ON) { fprintf(fdata,"K_FERMI ON\n"); } else { fprintf(fdata,"K_FERMI OFF\n"); }
        if (INITIAL == ON) { fprintf(fdata,"INITIAL ON\n"); } else { fprintf(fdata,"INITIAL OFF\n"); }
        if (IMPORT == ON) { printf("IMPORT ON\n");} else { printf("IMPORT OFF\n");}
        if (EXPORT == ON) { printf("EXPORT ON\n");} else { printf("EXPORT OFF\n");}
        if (ITER == ON) { fprintf(fdata,"ITERATION ON\n"); } else  { fprintf(fdata,"ITERATION OFF\n");} 
        if (SELFENERGY == ON) { fprintf(fdata,"SELFENERGY ON\n");} else { fprintf(fdata,"SELFENERGY OFF\n");} 
        if (CHEM == ON) { fprintf(fdata,"CHEM ON\n"); } else { fprintf(fdata,"CHEM OFF\n"); } 
        if (TUKEY == ON) { fprintf(fdata,"TUKEY ON\n");} else { fprintf(fdata,"TUKEY OFF\n");}
        if (U0SUMRULE == ON) { fprintf(fdata,"U0SUMRULE ON\n"); } else { fprintf(fdata,"U0SUMRULE OFF\n"); }
        if (DENSITYAPPROX == ON) { fprintf(fdata,"DENSITYAPPROX ON\n"); } else { fprintf(fdata,"DENSITYAPPROX OFF\n"); }


        fprintf(fdata, "u0SumRule Summary: \n");
        fprintf(fdata, "omega_domain= %lf\n", omega_domain);
        fprintf(fdata, "omega_max = %lf, omega_min = %lf\n", omegaf, omegai);
        fprintf(fdata, "Tukey Window: \n");
        fprintf(fdata, "OMEGAcminus = %lf\n", OMEGAcminus);
        fprintf(fdata, "OMEGAcplus = %lf\n", OMEGAcplus);
        fprintf(fdata, "Nomega = %d\n", Nomega);
        fprintf(fdata, "Nk (L) = %d\n", Nk);
        fprintf(fdata, "lambda = %lf\n", lambda);
        fprintf(fdata, "Density = %lf\n", nd);
        fprintf(fdata, "tau = %lf\n", tau);
        fprintf(fdata, "t = %lf\n", t);
        fprintf(fdata, "tp = %lf\n", tp);
        fprintf(fdata, "tpp = %lf\n", tpp);
        fprintf(fdata, "J = %lf\n", J);
        fprintf(fdata, "gamma_SE = %lf\n", gamma_SE);
        fprintf(fdata, "mup_initial = %lf\n", mup_initial);
        fprintf(fdata, "u0_initial = %lf\n", u0_initial);

        fprintf(fdata, "%s\n",  CONVERGENCE == 1 ? "SUCCESS" : "FAILED");
        fprintf(fdata, "iter = %d\n", iter);
        fprintf(fdata, "mup_final = %5.7lg\n", mup);
        fprintf(fdata, "u0_final = %5.7lg\n", u0);
        fprintf(fdata, "ng = %lg\n", ng);
        fprintf(fdata, "nG = %lg\n", nbigG);
        fprintf(fdata, "pseudo_nbigG = %5.7lg\n", pseudo_nbigG);
        fprintf(fdata, "mu = %5.7g\n", mu);
        fprintf(fdata, "kinetic energy = %lg\n", Tx + Ty + Txy);

        fprintf(fdata, "Tx = %lg\n", Tx);
        fprintf(fdata, "Ty = %lg\n", Ty);
        fprintf(fdata, "Txy = %lg\n", Txy);

        fprintf(fdata, "cxave = %lg\n", cxave);
        fprintf(fdata, "cyave = %lg\n", cyave);
        fprintf(fdata, "cxyave = %lg\n", cxyave);


        fprintf(fdata, "Fermi energy = %lg\n", Fermi_energy);
        fprintf(fdata, "mu0 = %lg\n", mup_root);
        fprintf(fdata, "mu = %5.7lg\n", mu);
     
        fprintf(fdata, "k_Fermi_nodal = %lg\n", k_Fermi_nodal);
        fprintf(fdata, "k_Fermi_antinodal = %lg\n", k_Fermi_antinodal);
        fprintf(fdata, "sigmaxx = %lg\n", sigmaxx);
        fprintf(fdata, "sigmayy = %lg\n", sigmayy);
        fprintf(fdata, "sigmaxy = %lg\n", sigmaxy);
        fprintf(fdata, "sigmac = %lg\n", sigmac);
        fprintf(fdata, "sigmacn = %lg\n", sigmacn);
        fprintf(fdata, "a1g = %lg\n", a1g);
        fprintf(fdata, "b1g = %lg\n", b1g);
        fprintf(fdata, "b2g = %lg\n", b2g);
        fclose(fdata);
    }


    if (CHECKSUM == ON) { printf("CHECKSUM ON\n"); } else { printf("CHECKSUM OFF\n"); }
    if (WRITE == ON) { printf("WRITE ON\n"); } else { printf("WRITE OFF\n"); }
    if (SETUP == ON) { printf("SETUP ON\n"); } else { printf("SETUP OFF\n"); }
    if (FERMI_ENERGY == ON) { printf("FERMI_ENERGY ON\n");} else { printf("FERMI_ENERGY OFF\n"); }
    if (K_FERMI == ON) { printf("K_FERMI ON\n"); } else { printf("K_FERMI OFF\n"); }
    if (INITIAL == ON) { printf("INITIAL ON\n"); } else { printf("INITIAL OFF\n"); }
    if (ITER == ON) { printf("ITERATION ON\n"); } else  { printf("ITERATION OFF\n");} 
    if (SELFENERGY == ON) { printf("SELFENERGY ON\n");} else { printf("SELFENERGY OFF\n");} 
    if (CHEM == ON) { printf("CHEM ON\n"); } else { printf("CHEM OFF\n"); } 
    if (TUKEY == ON) { printf("TUKEY ON\n");} else { printf("TUKEY OFF\n");}
    if (U0SUMRULE == ON) { printf("U0SUMRULE ON\n"); } else { printf("U0SUMRULE OFF\n"); }
    if (DENSITYAPPROX == ON) { printf("DENSITYAPPROX ON\n"); } else { printf("DENSITYAPPROX OFF\n"); }





    printf("Tukey Window: \n");
    printf("OMEGAcminus = %lf\n", OMEGAcminus);
    printf("OMEGAcplus = %lf\n", OMEGAcplus);
    printf("omega_domain = %lf\n", omega_domain);
    printf("omega_max = %lf, omega_min = %lf\n", omegaf, omegai);
    printf("Nomega = %d\n", Nomega);
    printf("Nk (L) = %d\n", Nk);
    printf("lambda = %lf\n", lambda);
    printf("Density = %lf\n", nd);
    printf("tau = %lf\n", tau);
    printf("t = %lf\n", t);
    printf("tp = %lf\n", tp);
    printf("tpp = %lf\n", tpp);
    printf("J = %lf\n", J);
    printf("gamma_SE = %lf\n", gamma_SE);
    printf("mup_initial = %lf\n", mup_initial);
    printf("u0_initial = %lf\n", u0_initial);

    printf("%s\n",  CONVERGENCE == 1 ? "SUCCESS" : "FAILED");
    printf("iter = %d\n", iter);
    printf("mup_final = %5.7lg\n", mup);
    printf("u0_final = %5.7lg\n", u0);
    printf("ng = %lg\n", ng);
    printf("nG = %lg\n", nbigG);
    printf("pseudo_nbigG = %5.7lg\n", pseudo_nbigG);
    printf("kinetic energy = %lg\n", Tx + Ty + Txy);
    printf("Tx = %lg\n", Tx);
    printf("Ty = %lg\n", Ty);
    printf("Txy = %lg\n", Txy);
    printf("cxave = %lg\n", cxave);
    printf("cyave = %lg\n", cyave);
    printf("cxyave = %lg\n", cxyave);
    printf("Fermi energy = %lg\n", Fermi_energy);
    printf("mu0 = %lg\n", mup_root);
    printf("mu = %5.7lg\n", mu);
    printf("k_Fermi_nodal = %lg\n", k_Fermi_nodal);
    printf("k_Fermi_antinodal = %lg\n", k_Fermi_antinodal);
    printf("chiBubrhorho00 = %5.6lf%+5.6lfi\n", creal(chiBubrhorho00), cimag(chiBubrhorho00));
    // printf("dn/dmu = %lg\n", dnoverdmu);
    printf("sigmaxx = %lg\n", sigmaxx);
    printf("sigmayy = %lg\n", sigmayy);
    printf("sigmaxy = %lg\n", sigmaxy);
    printf("sigmac = %lg\n", sigmac);
    printf("sigmacn = %lg\n", sigmacn);
    printf("a1g = %lg\n", a1g);
    printf("b1g = %lg\n", b1g);
    printf("b2g = %lg\n", b2g);


    free(rho_bigG_loc);
    free(re_bigG_loc);
    free(rho_g_loc);
    free(re_g_loc);

// }}} compute_data


// }}}


    free(re_g);
    free(re_bigG);

    free(rho_g);
    free(rho_bigG);


    free(rho_g_temp);


    free(Fermi_array);
    free(Fermibar_array);


    free(c);
    free(s);;
    free(xi);
    free(en);

    free(omega_array);
    free(k_array);

    timer("EOF", tic, tv1, TIME_LIMIT);

    return status;
}
