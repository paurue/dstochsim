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

#include<unistd.h>
#include "model.h"
#include "parser.h"
#include "methods.h"


int main(int argc, char ** argv){
	Model_t * m;
    int opterr, c;
	char fname[1000], algorithm[100];

	double time = 0, timestep = 1;
    opterr = 0;
    while ((c = getopt (argc, argv, "a:m:n:t:d:")) != -1)
      switch (c)
        {
        case 't':
          time = atof(optarg);
          break;
        case 'd':
          timestep = atof(optarg);
          break;
        case 'm':
          strcpy(fname, optarg);
          break;
        case 'a':
          strcpy(algorithm, optarg);
          break;
        case '?':
          if (optopt == 'c')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          else if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        default:
          abort ();
        }

	m = load_model_from_file(fname);
    if(strcmp(algorithm,"tleap") == 0) {
        sim_tleap(m, time, timestep);
    } else if(strcmp(algorithm,"nrk3l") == 0) {
        sim_nrk3l(m, time, timestep);
    } else if(strcmp(algorithm,"nrk3m") == 0) {
        sim_nrk3m(m, time, timestep);
    } else if(strcmp(algorithm,"nrk3h") == 0) {
        sim_nrk3h(m, time, timestep);
    } else if(strcmp(algorithm,"nrk5l") == 0) {
        sim_nrk5l(m, time, timestep);
    } else if(strcmp(algorithm,"nrk5m") == 0) {
        sim_nrk5m(m, time, timestep);
    } else if(strcmp(algorithm,"nrk5h") == 0) {
        sim_nrk5h(m, time, timestep);
    } else if(strcmp(algorithm,"heun") == 0) {
        sim_heun(m, time, timestep);
    } else {
        sim_direct_method(m, time, timestep);
    }
	free_model(m);
	return 0;
}
