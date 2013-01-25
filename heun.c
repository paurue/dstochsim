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

void sim_heun(Model_t * m, double tt, double tau){
    int i, j, step;
    int nreactions, nspecies;
    int **rs, **ps, **as;
    double *y2, *f, *f2;
    double *state;
    propensityFunc * prop;
    double **params;
    int nsteps;
    double *rates;

    /* Get the pointers */
    nreactions = m->Nreactions;
    nspecies = m->nspecies;
    prop = m->prop;
    params = m->params;
    state = dzeros(nspecies);
    rates = dzeros(nreactions);
    for(i=0; i<nspecies; i++) state[i] = (double) m->istate[i];
    rs = m->rstoichiometry;
    ps = m->pstoichiometry;
    as = m->acting_species;


    f  = dvector(nspecies);
    f2 = dvector(nspecies);
    y2  = dvector(nspecies);

    /* Header: column names */
    printf("#time ");
    for(i=0; i<nspecies; i++) printf("%s ", m->species[i]);
    printf("\n");
    printf("0 ");
    for(i=0; i<nspecies; i++) printf("%ld ", (long) state[i]);
    printf("\n");

    if(tau == 0){
        report_error("Heun method requires a strictly positive time step\n");
        exit(1);
    } else {
        nsteps = (int) ceil(tt / tau);
        for(step=0; step < nsteps; step++) {
            /* Step 0: Compute propensities and L(tau,x) = Pois(tau*x) -tau*x */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](state, nspecies, rs[j], params[j], as[j]);
            }
            /* Step 1: compute d = stoichiometry * L
             *                 f(y) = stoich. * propensities
             *                 and Y2
             */
            for(i=0; i<nspecies; i++){
                f[i] = 0;
                for(j=0; j<nreactions; j++) {
                    f[i] += (ps[j][i] - rs[j][i]) * rates[j];
                }
                /* Y2 = y + A21 * (tau * f(y)  + d) */
                y2[i] = state[i]  + tau * f[i];
            }
            /* Step 2:  Compute propensities for Y2 */
            for(j=0; j< nreactions; j++){
                rates[j] = prop[j](y2, nspecies, rs[j], params[j], as[j]);
            }
            /* Step 3: compute f(Y2) = stoich. * propensities
             *                 and state[n+1]
             */
            for(i=0; i<nspecies; i++){
                f2[i] = 0;
                for(j=0; j<nreactions; j++) {
                    f2[i] += (ps[j][i] - rs[j][i]) * rates[j];
                }
                /* y_{n+1} = y_{n} + h/2 * (f(t, Y1) + f(t + h, Y2)) */
                state[i] = state[i]  + 0.5 * tau * (f[i] + f2[i]);
            }
            printf("%g ", tau * (step+1));
            for(i=0; i<nspecies; i++) printf("%g ", state[i]);
            printf("\n");
        }
    }
    return;
}
