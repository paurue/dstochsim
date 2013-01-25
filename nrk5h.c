/*
 * Copyright (c) 2013 Pau Rué <pau.rue@gmail.com>
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the "Software"), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the Software is 
 * furnished to do so, subject to the following conditions:
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "methods.h"
#include "time.h"
#include<unistd.h>

#define NRK5_A21 3.512099547699939e-02
#define NRK5_A32 9.333225239662185e-02
#define NRK5_A43 1.940500426863389e-01
#define NRK5_A54 3.552631979151575e-01

#ifdef PRINT_RUNTIME
clock_t start, end;
#endif

void sim_nrk5h(Model_t * m, double tt, double tau){
    int i, j, step;
    long seed;
    int nreactions, nspecies;
    int **rs, **stoich, **as;
    double *L, *d, *y, *f;
    double *state;
    propensityFunc * prop;
    double **params;
    int nsteps;
    double *rates;
    const gsl_rng_type * type = gsl_rng_default;
    gsl_rng * r;

    /* Get the pointers */
    nreactions = m->Nreactions;
    nspecies = m->nspecies;
    stoich = imatrix(nspecies, nreactions);
    rs = m->rstoichiometry;
    for(i=0; i< nspecies;i++){
        for(j=0; j< nreactions; j++){
            stoich[i][j] = (m->pstoichiometry[j][i] - rs[j][i]);
        }
    }
    prop = m->prop;
    params = m->params;
    state = dzeros(nspecies);
    rates = dzeros(nreactions);
    for(i=0; i<nspecies; i++) state[i] = (double) m->istate[i];

    as = m->acting_species;


    gsl_rng_env_setup();
    r = gsl_rng_alloc (type);
    seed = time(NULL) * getpid();
    gsl_rng_set (r, seed);                  // set seed

    L = dvector(nreactions);
    d = dvector(nspecies);
    f = dvector(nspecies);
    y = dvector(nspecies);

    #ifdef OUTPUT_SPECIES
    /* Header: column names */
    printf("#time ");
    for(i=0; i<nspecies; i++) printf("%s ", m->species[i]);
    printf("\n");
    printf("0 ");
    for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
    printf("\n");
    #endif
    if(tau == 0){
        report_error("Tau-leap requires a strictly positive time step\n");
        exit(1);
    } else {
        #ifdef PRINT_RUNTIME
        start = clock();
        #endif
        nsteps = (int) ceil(tt / tau);
        for(step=0; step < nsteps; step++) {
            /* Step 0: Compute propensities and L(tau,x) = Pois(tau*x) -tau*x */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](state, nspecies, rs[j], params[j], as[j]);
                L[j] = gsl_ran_poisson (r, tau * rates[j]) - tau * rates[j];
            }
            /* Step 1: compute d = stoichiometry * L
             *                 f(y) = stoich. * propensities and Y2
             */
            for(i=0; i<nspecies; i++){
                d[i] = 0;
                f[i] = 0;
                for(j=0; j<nreactions; j++) {
                    d[i] += (stoich[i][j]) * L[j];
                    f[i] += (stoich[i][j]) * rates[j];
                }
                /* Y2 = y + A21 * (tau * f(y)  + d) */
                y[i] = state[i]  + NRK5_A21 * (tau * f[i] + d[i]);
            }

            /* Step 2:  Compute propensities for Y2 */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](y, nspecies, rs[j], params[j], as[j]);
            }
            /* Step 3: compute f(Y2) = stoich. * propensities
             *                 and Y3
             */
            for(i=0; i<nspecies; i++){
                f[i] = 0;
                for(j=0; j<nreactions; j++) {
                    f[i] += (stoich[i][j]) * rates[j];
                }
                /* Y3 = y + A32 * (tau * f(Y2)  + d) */
                y[i] = state[i]  + NRK5_A32 * (tau * f[i] + d[i]);
            }
            /* Step 4:  Compute propensities for Y3 */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](y, nspecies, rs[j], params[j], as[j]);
            }
            /* Step 5: compute f(Y3) = stoich. * propensities
             *                 and Y4
             */
            for(i=0; i<nspecies; i++){
                f[i] = 0;
                for(j=0; j<nreactions; j++) {
                    f[i] += (stoich[i][j]) * rates[j];
                }
                /* Y4 = y + A43 * (tau * f(Y3)  + d) */
                y[i] = state[i]  + NRK5_A43 * (tau * f[i] + d[i]);
            }
            /* Step 6:  Compute propensities for Y4 */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](y, nspecies, rs[j], params[j], as[j]);
            }
            /* Step 7: compute f(Y4) = stoich. * propensities
             *                 and Y5
             */
            for(i=0; i<nspecies; i++){
                f[i] = 0;
                for(j=0; j<nreactions; j++) {
                    f[i] += (stoich[i][j]) * rates[j];
                }
                /* Y5 = y + A54 * (tau * f(Y4)  + d) */
                y[i] = state[i]  + NRK5_A54 * (tau * f[i] + d[i]);
            }
            /* Step 6:  Compute propensities for Y5 */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](y, nspecies, rs[j], params[j], as[j]);
            }
            /* Step 7: compute f(Y5) = stoich. * propensities and states[n+1]
             */
            for(i=0; i<nspecies; i++){
                f[i] = 0;
                for(j=0; j<nreactions; j++) {
                    f[i] += (stoich[i][j]) * rates[j];
                }
                /* y = ROUND(y + tau * f(Y3)  + d) */
                state[i] += (tau * f[i] + d[i]);
                if(state[i]<0) state[i] = 0;
            }
            #ifdef OUTPUT_SPECIES
            printf("%g ", tau * (step+1));
            for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
            printf("\n");
            #endif
        }
        #ifdef PRINT_RUNTIME
        end = clock();
        printf("%g ", (double) (end - start)/CLOCKS_PER_SEC);
        #endif
        #ifdef OUTPUT_SPECIES
        printf("%g ", tau * (step+1));
        for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
        printf("\n");
        #endif
    }
    gsl_rng_free(r);
    return;
}
