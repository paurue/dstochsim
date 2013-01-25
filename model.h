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

#ifndef DSSUTILS_H_
#define DSSUTILS_H_

#include "utils.h"

typedef double (*propensityFunc)(double * state, int nreactants, int * rstoichiometry, double *params, int* acting_species);

typedef struct _Model_t {
	int nspecies;
	int Nreactions;
    int * nparams;
	long * istate;
	double * dstate;
	char ** species;
	long * ics;
	double time;
    propensityFunc * prop;
    double ** params;
    int ** acting_species; /* Some reaction types need these. Such as the propensity depending on another variable */
    int **rstoichiometry, **pstoichiometry;
} Model_t;

Model_t * model_new();

void free_model(Model_t *model);

void model_set_allocate(Model_t * model, int nspecies, int nreactions);

void model_print(Model_t * m);

double prop_MA(double *x , int nx, int *c, double *params, int * acting_species);
double prop_HA(double *x , int nx, int *c, double *params, int * acting_species);
double prop_HI(double *x , int nx, int *c, double *params, int * acting_species);
double prop_MAHI(double *x , int nx, int *c, double *params, int * acting_species);
double prop_HIHA(double *x , int nx, int *c, double *params, int * acting_species);
double prop_CI(double *x , int nx, int *c, double *params, int * acting_species);
double prop_HAHAC(double *x , int nx, int *c, double *params, int * acting_species);
double prop_HAHAHIC(double *x , int nx, int *c, double *params, int * acting_species);
double prop_TEST(double *x , int nx, int *c, double *params, int * acting_species);
#endif /* DSSUTILS_H_ */
