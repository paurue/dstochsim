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

#include "parser.h"

#define REACTION_TYPES 4

typedef enum _reaction_type2{
	NOT_VALID,
	MASS_ACTION,
	HILL_ACTIVATION,
	HILL_INHIBITION
} reaction_type;

const char reaction_types[][REACTION_TYPES] = {"\0", "MA", "HA", "HI"};

typedef enum _stoich_type{
	REACTANTS,
	PRODUCTS
} stoich_type;

typedef enum _parsing_section{
	PARSING_NONE = 0,
	PARSING_SPECIES = 1 ,
	PARSING_REACTIONS = 2,
} parsing_section;

void parse_line_species(Model_t * model, List_t * lines) {
	/* Species lines have the format:
	 * SpeciesName1 = 2312 # Some comment
	 * SpeciesName2 = 0
	 * We first create a list of species to ensure all of them are unique!
	 * */
    int i, idx;
    List_t * species_list;
    iList_t * ics_list;
    char * species_str, * ic_str, * saveptr;
    species_list = list_new();
    ics_list = ilist_new();

    for(i=0; i< lines->size; i++) {
        species_str = strtok_r(lines->items[i], "=", &saveptr);
        trim(species_str);
        /*TODO: Allow for more than one species per line
         *     for(p=strsep(&values, ","); p != NULL; p=strsep(&values, ",")) {
         *      value = atoi(p);
         */
        ic_str = strtok_r (NULL, "=", &saveptr);
        trim(ic_str);
        idx = string_find(species_str, model->species, model->nspecies);
        if(idx == -1){ /* Species not in list */
            list_append(species_list, species_str);
            ilist_append(ics_list, atoi(ic_str));
        } else {
            report_warning("Species '%s' is being redefined", species_str);
            ics_list->items[idx] = atoi(ic_str);
        }
    }
    /* Fill the model */
    model->nspecies = species_list->size; 
	/* in case we have less species than lines (repeated species)*/
    for(i=0; i< model->nspecies; i++){
        strcpy(model->species[i], species_list->items[i]);
        model->ics[i] = ics_list->items[i];
        model->istate[i] = ics_list->items[i];
    }
	return;
}

void parse_stoichiometry(char * str, Model_t * model, int ireaction, stoich_type type) {
	/* Format is 2*A + 3*B + C */
	int coeff;
    int nspecies;
	int species[2000], stoichiometry[2000];
	char * aux_str1, * aux_str2, *aux_str3, *saveptr1, *saveptr2;
	int k;
	aux_str1 = strtok_r (str, "+", &saveptr1);
	nspecies = 0;
	while(aux_str1 != NULL) {
		trim(aux_str1);
		//coeff=0;
        //printf(".--.%s",aux_str1);
		if(strlen(aux_str1)==0 || strcmp(aux_str1, "0")!=0){
			aux_str2 = strtok_r(aux_str1, "*", &saveptr2);
			aux_str3 = strtok_r(NULL, "*", &saveptr2);
			if(aux_str3 == NULL) { /* Case "A" instead of "1*A"*/
				coeff = 1;
			} else { /* Case "3*A" */
				coeff = atoi(aux_str2);
				aux_str2 = aux_str3;
			}
			k = string_find(aux_str2, model->species, model->nspecies);
			if(k == -1){
				report_error("Species %s has not been defined", aux_str1);
				exit(1);
			} else if(type==REACTANTS){
			    model->rstoichiometry[ireaction][k] = coeff;
			} else if(type==PRODUCTS) {
	            model->pstoichiometry[ireaction][k] = coeff;
			}
			species[nspecies] = k;
			stoichiometry[nspecies] = coeff;
			nspecies++;
		}
		aux_str1 = strtok_r (NULL, "+", &saveptr1);
	}
    return;
}
void parse_params(char * params_str, Model_t * model, int ireaction, char *rtype) {
    char * saveptr;
    char * aux_str;
    int idx;
    if (strcmp(rtype, "MA") == 0) {
        model->params[ireaction] = (double *) malloc(sizeof(double));
        model->params[ireaction][0] = atof(params_str);
    } else if ((strcmp(rtype, "HA") == 0) || (strcmp(rtype, "HI") == 0) ) {
        model->params[ireaction] = (double *) malloc(3 * sizeof(double)); /* rate, Ks, coop.*/
        model->acting_species[ireaction] = (int *) malloc(sizeof(int));
        /* Find Which species is acting */
        trim(params_str);
        aux_str = strtok_r(params_str, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][0] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate */
        model->params[ireaction][0] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant */
        model->params[ireaction][1] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity */
        model->params[ireaction][2] = atof(aux_str);
    } else if ((strcmp(rtype, "HIHA") == 0) || (strcmp(rtype, "HIHI") == 0) || (strcmp(rtype, "HAHA") == 0) || (strcmp(rtype, "CI") == 0) ) {
        //printf("HIHA\n");
        model->params[ireaction] = (double *) malloc(5 * sizeof(double)); /* rate, Ks1, coop1, Ks2, coop2.*/
        /* for CI: rate, Ks, coop1, gamma, coop2. -> rate * (y1)^coop1 / ((y1)^coop1 + (Ks)^coop2 + (gamma*y2)^coop2) */
		model->acting_species[ireaction] = (int *) malloc(2 * sizeof(int));
        trim(params_str);
        aux_str = strtok_r(params_str, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][0] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][1] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate */
        model->params[ireaction][0] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 1*/
        model->params[ireaction][1] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity 1*/
        model->params[ireaction][2] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 2*/
        model->params[ireaction][3] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity 2 */
        model->params[ireaction][4] = atof(aux_str);
    } else if ((strcmp(rtype, "MAHA") == 0) || (strcmp(rtype, "MAHI") == 0)) {
        model->params[ireaction] = (double *) malloc(3 * sizeof(double)); /* rate, Ks1, coop1.*/
        model->acting_species[ireaction] = (int *) malloc(2 * sizeof(int));
        trim(params_str);
        aux_str = strtok_r(params_str, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][0] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][1] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate */
        model->params[ireaction][0] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 1*/
        model->params[ireaction][1] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity*/
        model->params[ireaction][2] = atof(aux_str);
    } else if ((strcmp(rtype, "HAHAC") == 0) ) {
        model->params[ireaction] = (double *) malloc(6 * sizeof(double)); /* rate1, Ks1, coop1, rate2, Ks2, coop2.*/
        /* for CI: rate, Ks, coop1, gamma, coop2. -> rate * (y1)^coop1 / ((y1)^coop1 + (Ks)^coop2 + (gamma*y2)^coop2) */
		model->acting_species[ireaction] = (int *) malloc(2 * sizeof(int));
        trim(params_str);
        aux_str = strtok_r(params_str, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][0] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][1] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate 1*/
        model->params[ireaction][0] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 1*/
        model->params[ireaction][1] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity 1*/
        model->params[ireaction][2] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate 2*/
        model->params[ireaction][3] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 2*/
        model->params[ireaction][4] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity 2 */
        model->params[ireaction][5] = atof(aux_str);
    } else if ((strcmp(rtype, "HAHAHIC") == 0) ) {
        model->params[ireaction] = (double *) malloc(8 * sizeof(double)); /* rate1, Ks1, coop1, rate2, Ks2, coop2, Ksi, coopi.*/
        /* for CI: rate, Ks, coop1, gamma, coop2. -> rate * (y1)^coop1 / ((y1)^coop1 + (Ks)^coop2 + (gamma*y2)^coop2) */
		model->acting_species[ireaction] = (int *) malloc(3 * sizeof(int));
        trim(params_str);
        aux_str = strtok_r(params_str, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][0] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][1] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* Species name */
        idx = string_find(aux_str, model->species, model->nspecies);
        if(idx == -1) {
            report_error("Species '%s' not found", aux_str);
            exit(1);
        } else {
            model->acting_species[ireaction][2] = idx;
        }
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate 1*/
        model->params[ireaction][0] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 1*/
        model->params[ireaction][1] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity 1*/
        model->params[ireaction][2] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* rate 2*/
        model->params[ireaction][3] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant 2*/
        model->params[ireaction][4] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity 2 */
        model->params[ireaction][5] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* binding constant inhib*/
        model->params[ireaction][6] = atof(aux_str);
        aux_str = strtok_r(NULL, " ", &saveptr); /* cooperativity inhib */
        model->params[ireaction][7] = atof(aux_str);
    } else {
        report_error("Reaction type '%s' not recognised", rtype);
        exit(1);
    }
}

void parse_line_reactions(Model_t * model, List_t * reactions_lines) {
    int i;
    char rtype[100];
    char * aux_str;
	char * saveptr1;
	for(i=0; i<reactions_lines->size; i++) {
        // Reactants
	    aux_str = strtok_r(reactions_lines->items[i], ">", &saveptr1);
	    trim(aux_str);
        parse_stoichiometry(aux_str, model, i, REACTANTS);
        // Products
        aux_str = strtok_r(NULL, "|", &saveptr1);
        trim(aux_str);
        parse_stoichiometry(aux_str, model, i, PRODUCTS);
        // Reaction types
	    aux_str = strtok_r(NULL, "|",&saveptr1);
	    strcpy(rtype, aux_str);
	    trim(rtype);
		// printf("type: %s\n", rtype);
	    if(strcmp(rtype, "MA") == 0) {
	        model->prop[i] = prop_MA;
	    } else if(strcmp(rtype, "HA") == 0) {
            model->prop[i] = prop_HA;
        } else if(strcmp(rtype, "HI") == 0) {
            model->prop[i] = prop_HI;
        } else if(strcmp(rtype, "HIHA") == 0) {
            model->prop[i] = prop_HIHA;
        } else if(strcmp(rtype, "MAHI") == 0) {
            model->prop[i] = prop_MAHI;
        } else if(strcmp(rtype, "CI") == 0) {
            model->prop[i] = prop_CI;
		} else if(strcmp(rtype, "HAHAHIC") == 0) {
            model->prop[i] = prop_HAHAHIC;
		} else if(strcmp(rtype, "HAHAC") == 0) {
            model->prop[i] = prop_HAHAC;
	    }else{
			printf("Warning! Propensity type not detected\n");
	        model->prop[i] = prop_TEST;
	    }
	    //Parameters
	    aux_str = strtok_r (NULL, "|", &saveptr1);
        parse_params(aux_str, model, i, rtype);
	}
    return;
}

Model_t * load_model_from_file(char * fname) {
	FILE * in;
	Model_t * model;
	char line[MAX_LINE_SIZE];
	char section_name[MAX_LINE_SIZE];
	int lastchr; /* Last character read */
	int linenum, len; /* line number and length*/
    parsing_section section;
    List_t * species_lines, * reactions_lines;

	/* Lists initialisation */
    species_lines = list_new();
    reactions_lines = list_new();

    /* Here we get the strings for the species (and its initial conditions),
     * and the reactions.
     * */
	in = fopen(fname, "r");
	if (!in){
		report_error("cannot open file");
		exit(1);
	}
	linenum = 0;
	lastchr = 0;
	section = PARSING_NONE;
	while (fgets(line + lastchr, MAX_LINE_SIZE-lastchr, in) != NULL) {
		linenum++;
		len = (int) strlen(line)-1;
		/* Avoid buffer overflows */
		if(line[len]!='\n'){
			report_error("line too long to be parsed in file");
		fclose(in);
	    exit(1);
		}
		/* Get rid of \n and spaces at end of line */
		while ((len>=0) && ((line[len]=='\n') || (isspace(line[len])))) {
			line[len]='\0' ;
			len-- ;
		}
		len = remove_comments(line, '#');
		if(len < 1) { /* Skip empty line*/
	    } else if (line[0]=='[' && line[len]==']') {
	        /* Section name line
	         * TODO: Allow for different case  and check the output from sscanf
	         * */
	        sscanf(line, "[%[^]]", section_name);
			if(strcmp(section_name, "Species") == 0){
				section = PARSING_SPECIES;
			} else if(strcmp(section_name, "Reactions") == 0) {
				if(section != PARSING_SPECIES){
		        	report_error("In file: %s, line %d: Reactions must be specified after species", fname, linenum);
		        	exit(1);
				}
				section = PARSING_REACTIONS;
	        } else {
	        	report_error("In file: %s, section %s not valid", fname, section_name);
	        	exit(1);
	        }
		} else if(section == PARSING_SPECIES){
            list_append(species_lines, line);
		} else if(section == PARSING_REACTIONS){
			list_append(reactions_lines, line);
		} else {
			report_error("Line not correctly formatted!");
			exit(1);
		}
	}
	fclose(in);

	model = model_new();
    model_set_allocate(model, species_lines->size, reactions_lines->size);
    parse_line_species(model, species_lines);
    parse_line_reactions(model, reactions_lines);

	return model;
}
