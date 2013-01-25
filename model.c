/*
 * Copyright (c) 2013 Pau RuÈ <pau.rue@gmail.com>
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

#include "model.h"
#define MAX_NAME_SIZE 128


Model_t * model_new(){
	Model_t * model;

	model = (Model_t *) malloc(sizeof(Model_t));
	if (!model) {
		report_error("allocation failure in model()");
		exit(1);
	}
	return model;
}

void free_model(Model_t * model){
	/* We should add reactions and species here also*/
	//int i;
/*TODO:Do it right
 	for(i=0; i< model->Nreactions; i++)
        if(model->rstoichiometry[i] != NULL)
            free((char *) model->rstoichiometry[i]);
        if(model->pstoichiometry[i] != NULL)
            free((char *) model->pstoichiometry[i]);
	    if(model->params[i] != NULL)
	        free((char *) model->params[i]);
	printf("=\n");
    if(model->rstoichiometry != NULL)
        free((char *) model->rstoichiometry);
    if(model->pstoichiometry != NULL)
        free((char *) model->pstoichiometry);
    if(model->params != NULL)
        free((char *) model->params);
*/
	free((char *) model);
}

void model_set_allocate(Model_t * model, int nspecies, int nreactions) {
	/* Allocate memory just for Nreactants, Nprodcuts, Nparams and the "rows" of
	 * reactants, prodcuts, params
	 * */
    int i=0;
	int flag_error;
    model->nspecies = nspecies;
    model->Nreactions = nreactions;
	flag_error = 0;
    model->prop = (propensityFunc *) malloc(nreactions * sizeof(propensityFunc *));
	if (model->prop == NULL) flag_error = 1;
	model->nparams = izeros(nreactions);
	if (model->nparams == NULL) flag_error = 1;
    model->rstoichiometry = (int **) calloc(nreactions, sizeof(int *));
    if (model->rstoichiometry == NULL) flag_error = 1;
    model->pstoichiometry = (int **) calloc(nreactions, sizeof(int *));
    if (model->pstoichiometry == NULL) flag_error = 1;
	model->params = (double **) malloc(nreactions * sizeof(double *));
	if (model->params == NULL) flag_error = 1;
    model->acting_species = (int **) calloc(nreactions, sizeof(int *));
    if (model->acting_species == NULL) flag_error = 1;
	model->species = (char **) malloc(nreactions * sizeof(char *));
    if (model->species == NULL) flag_error = 1;
	if(flag_error) {
	    report_error("allocation failure in model_set_allocate()");
	    exit(1);
	 }
	for(i=0; i<nreactions; i++) {
	    model->rstoichiometry[i] = (int *) calloc(nspecies, sizeof(int));
        model->pstoichiometry[i] = (int *) calloc(nspecies, sizeof(int));
        model->species[i] = (char *) malloc(MAX_NAME_SIZE * sizeof(char));
	}
	model->istate =lzeros(nspecies);
	model->ics = lzeros(nspecies);
	return;
}

void model_print(Model_t * m){
    int i, j;
    printf("# [Species] # %d species\n", m->nspecies);
    for(i=0; i<m->nspecies; i++){
        printf("%s = %ld\n", m->species[i], m->ics[i]);
    }
    printf("# [Reactions] # %d reactions\n", m->Nreactions);
    for(i=0; i<m->Nreactions;i++) {
        printf("# ");
        for(j=0; j<m->nspecies-1; j++) {
            if(m->rstoichiometry[i][j] == 1) printf("%s + ", m->species[j]);
            else if(m->rstoichiometry[i][j] > 1) printf("%d*%s + ", m->rstoichiometry[i][j], m->species[j]);
        }
        if(m->rstoichiometry[i][m->nspecies-1] == 1) printf("%s", m->species[m->nspecies-1]);
        else if(m->rstoichiometry[i][m->nspecies-1] > 1) printf("%d*%s", m->rstoichiometry[i][m->nspecies-1], m->species[j]);
        printf(" -> ");
        for(j=0; j<m->nspecies-1; j++) {
            if(m->pstoichiometry[i][j] == 1) printf("%s + ", m->species[j]);
            else if(m->pstoichiometry[i][j] > 1) printf("%d*%s + ", m->pstoichiometry[i][j], m->species[j]);
        }
        if(m->pstoichiometry[i][m->nspecies-1] == 1) printf("%s", m->species[m->nspecies-1]);
        else if(m->pstoichiometry[i][m->nspecies-1] > 1) printf("%d*%s", m->pstoichiometry[i][m->nspecies-1], m->species[j]);
        printf(" | \n");
    }
}
void model_print_state(Model_t * m){
    int i;
    printf("%g ",m->time);
    for(i=0; i<m->nspecies-1;i++){
        printf("%ld ",m->istate[i]);
    }
    printf("%ld\n",m->istate[m->nspecies-1]);

}

/* Propensity reactions*/
double prop_MA(double *x , int nx, int *c, double *params, int * acting_species){
	/* Mass Action Law propensity
	 * For species Xi with stoichiometric coefficients Ci, the propensity is given by:
	 * 	rate * binomial(X1, C1) * ··· * binomial(XN, CN)
	 *
	 * */
	int i;
	double prop;
	prop = params[0];
	for(i=0; i < nx; i++){
		if(c[i]>0)
		    //printf("x(%d)=%d, %d\n",i,x[i],c[i]);
			prop *= dchoose(x[i], c[i]);
	}
//	printf("prop=%g\n",prop);
	return prop;
}

double prop_HA(double *x , int nx, int *c, double *params, int * acting_species){
    /* Hill Activation
     * Reaction rate depends on an extra
     *
     * */
    int y;
    double rate, km, hcoop, prop;
    rate = params[0];
    km = params[1];
    hcoop = params[2];
    y = x[acting_species[0]];
    //printf("rate=%g, km=%g, hcoop=%g  -> %d\n", rate, km, hcoop, y);
    prop = rate / (1 + pow(km/y, hcoop));
    return prop;
}

double prop_HI(double *x , int nx, int *c, double *params, int * acting_species){
    /* Hill Inhibition
     * Reaction rate depends on an extra
     *
     * */
    int y;
    double rate, km, hcoop, prop;
    rate = params[0];
    km = params[1];
    hcoop = params[2];

    y = x[acting_species[0]];
    prop = rate / (1 + pow(y/km, hcoop));
    return prop;
}


double prop_MAHI(double *x , int nx, int *c, double *params, int * acting_species){
    /* Hill Inhibition
     * Reaction rate depends on an extra
     *
     * */
    int ym, yi;
    double rate, kmi, hcoopi, prop;
    rate = params[0];
    kmi = params[1];
    hcoopi = params[2];

    ym = x[acting_species[0]];
    yi = x[acting_species[1]];
    prop = rate * ym / (1 + pow(yi/kmi, hcoopi));
    return prop;
}

double prop_CI(double *x , int nx, int *c, double *params, int * acting_species){
    /* Competitive Inhibition
     *
     * */
    int yi, ya;
    double rate, gamma, hcoopi, kma, hcoopa, prop;
    rate = params[0];
    kma = params[1];
    hcoopa = params[2];
    gamma = params[3];
    hcoopi = params[4];
    //printf("rate=%g, kma=%g, hcoopa=%g, gamma=%g, hcoopi=%g ", rate, kma, hcoopa, gamma, hcoopi);
    ya = x[acting_species[0]];
    yi = x[acting_species[1]];

    if(ya==0) {prop = 0; }
    else { prop = rate * pow(ya, hcoopa) / (pow(kma, hcoopa) + pow(ya, hcoopa) + pow(gamma*yi, hcoopi)); }
    return prop;
}


double prop_HIHA(double *x , int nx, int *c, double *params, int * acting_species){
    /* Hill Inhibition
     * Reaction rate depends on an extra
     *
     * */
    int yi, ya;
    double rate, kmi, hcoopi, kma, hcoopa, prop;
    rate = params[0];
    kmi = params[1];
    hcoopi = params[2];
    kma= params[3];
    hcoopa = params[4];

    yi = x[acting_species[0]];
    ya = x[acting_species[1]];

    if(ya==0) {prop = 0; }
    else { prop = rate / (1 + pow(yi/kmi, hcoopi)) / (1 + pow(kma/ya, hcoopa)); }
    return prop;
}


double prop_HAHAC(double *x , int nx, int *c, double *params, int * acting_species){
    /* 2 competitive activators
     * Reaction rate depends on an extra
     *
     * */
    int ya1, ya2;
    double rate1, rate2, kma1, hcoop1, kma2, hcoop2, prop;
    rate1 = params[0];
    kma1 = params[1];
    hcoop1 = params[2];
    
	rate2 = params[3];
    kma2= params[4];
    hcoop2 = params[5];

    ya1 = x[acting_species[0]];
    ya2 = x[acting_species[1]];

    prop = (rate1*pow(ya1/kma1, hcoop1) + rate2*pow(ya2/kma2, hcoop2) )/ (1 + pow(ya1/kma1, hcoop1) + pow(ya2/kma2, hcoop2));
    return prop;
}

double prop_HAHAHIC(double *x , int nx, int *c, double *params, int * acting_species){
    /* 3 competitive species, 2 activators and 1 inhibitor
     * Reaction rate depends on an extra
     *
     * */
    int ya1, ya2, yi;
    double rate1, rate2, kma1, hcoop1, kma2, hcoop2, kmi, hcoopi, prop;
    rate1 = params[0];
    kma1 = params[1];
    hcoop1 = params[2];

    rate2 = params[3];
    kma2= params[4];
    hcoop2 = params[5];

	kmi = params[6];
	hcoopi = params[7];

    ya1 = x[acting_species[0]];
    ya2 = x[acting_species[1]];
    yi  = x[acting_species[2]];

    prop = (rate1*pow(ya1/kma1, hcoop1) + rate2*pow(ya2/kma2, hcoop2) )/ (1 + pow(ya1/kma1, hcoop1) + pow(ya2/kma2, hcoop2) + pow(yi/kmi, hcoopi));

    return prop;
}



/* Propensity reactions*/
double prop_TEST(double *x , int nx, int *c, double *params, int * acting_species){
	return 42;
}
