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

clock_t start, end;

void sim_direct_method(Model_t * m, double tt, double hurdle){
    int i, j, step;
    long seed;
    int nreactions, nspecies;
    int **rs, **ps, **as;
    double *state;
    propensityFunc * prop;
    double **params;
    int nsteps;
    double t, tau, a0, r1, r2, runningSum, thr, nextHurdle;
    double * rates;
    const gsl_rng_type * type = gsl_rng_default;
    gsl_rng * r;

    /* Get the pointers */
    nreactions = m->Nreactions;
    nspecies = m->nspecies;
    prop = m->prop;
    params = m->params;
    state = dzeros(nspecies);
    rates = dvector(m->Nreactions);
    for(i=0; i<nspecies; i++) state[i] = (double) m->istate[i];
    rs = m->rstoichiometry;
    ps = m->pstoichiometry;
    as = m->acting_species;

    gsl_rng_env_setup();
    r = gsl_rng_alloc (type);
    seed = time(NULL) * getpid();
    gsl_rng_set (r, seed);                  // set seed

    t = 0;

    /* Header: column names */
    printf("#time ");
    for(i=0; i<nspecies; i++) printf("%s ", m->species[i]);
    printf("\n");
    printf("%g ", t);
    for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
    printf("\n");

    if(hurdle == 0){
        /* */
        report_warning("Not yet implemented\n");
    } else {
        nsteps = (int) ceil(tt / hurdle);
        step = 0;
        nextHurdle = hurdle;
        //while(t < time){
        start = clock();
        while(step < nsteps){
            for(i=0; i< nreactions; i++){
                rates[i] = prop[i](state, nspecies, rs[i], params[i], as[i]);
            }
            a0 = dsum(rates, nreactions);

            /* Sample tau and update time*/
            r1 = gsl_rng_uniform_pos (r);
            if(a0 > 0){
                tau = (-1/a0) * log(r1);
            } else {
                /* No more reactions are likely to occur*/
                tau = tt;
                return;
            }

            t = t + tau;

            /* Sample reaction j */
            r2 = gsl_rng_uniform_pos(r);
            thr = a0 * r2;
            runningSum = 0;
            for(j=0; j<nreactions; j++){
                runningSum += rates[j];
                if(runningSum > thr) break;
            }
            /* Hurdles
             * if t reaches nextHurdle, the system state at t=nextHurdle is
             * the system state before updating. That applies for all the
             * following hurdles t reaches.
             */
            i = t > nextHurdle;
            i = t < tt;
            while(t > nextHurdle){
                step += 1;
                printf("%g ",nextHurdle);
                for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
                printf("\n");
                nextHurdle += hurdle;
            }
            /* Species update */
            for(i=0; i<nspecies; i++){
                state[i] += ps[j][i] -rs[j][i];
            }
        }
        end = clock();
		printf("%g ", nextHurdle);
        for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
        printf("\n");
    }
    gsl_rng_free(r);
    return;
}
