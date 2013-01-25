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

#ifdef PRINT_RUNTIME
clock_t start, end;
#endif

void sim_tleap(Model_t * m, double tt, double tau){
    int i, j, step;
    long seed;
    int nreactions, nspecies;
    int **rs, **stoich, **as, *K;
    double *state;
    propensityFunc * prop;
    double **params;
    int nsteps;
    double rate;
    const gsl_rng_type * type = gsl_rng_default;
    gsl_rng * r;

    /* Get the pointers */
    nreactions = m->Nreactions;
    nspecies = m->nspecies;
    /* Convention: rows correspond to species while columns to reactions, thus
     * stoich[i][j] refers to the stoichiometry of the species i due to reaction j */
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
    for(i=0; i<nspecies; i++) state[i] = (double) m->istate[i];
    as = m->acting_species;

    gsl_rng_env_setup();
    r = gsl_rng_alloc (type);
    seed = time(NULL) * getpid();
    gsl_rng_set (r, seed);                  // set seed

    K = ivector(nreactions);

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
            for(j=0; j< nreactions; j++){
                rate = prop[j](state, nspecies, rs[j], params[j], as[j]);
                K[j] = gsl_ran_poisson (r, tau * rate);
            }
            /* Species update */
            for(i=0; i<nspecies; i++){
                for(j=0; j<nreactions; j++) {
                    state[i] += K[j] * stoich[i][j];
                }
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
    free_ivector(K);
    gsl_rng_free(r);
    return;
}
