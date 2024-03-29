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

#ifndef METHODS_H_
#define METHODS_H_


#define PRINT_RUNTIME

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include "model.h"


void sim_direct_method(Model_t * m, double tt, double hurdle);
void sim_tleap(Model_t * m, double tt, double tau);
void sim_nrk3l(Model_t * m, double tt, double tau);
void sim_nrk3m(Model_t * m, double tt, double tau);
void sim_nrk3h(Model_t * m, double tt, double tau);
void sim_nrk5l(Model_t * m, double tt, double tau);
void sim_nrk5m(Model_t * m, double tt, double tau);
void sim_nrk5h(Model_t * m, double tt, double tau);

void sim_heun(Model_t * m, double tt, double hurdle);


#endif /* METHODS_H_ */
