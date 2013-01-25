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


#ifndef UTILS_H_
#define UTILS_H_
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"

#define LIST_MAX_SIZE 100
#define LIST_MAX_ITEM_SIZE 100


typedef struct _List_t
/* List structure */
{
        char items[LIST_MAX_SIZE][LIST_MAX_ITEM_SIZE];
        int size;
} List_t;
typedef struct _iList_t
/* iList structure */
{
        int items[LIST_MAX_SIZE];
        int size;
} iList_t;
/*
 * TODO: Implement a dynamic list
 * */


double *dvector(long n);
/* Allocate a double vector of size n */

double *dzeros(long n);
/* Allocate a double vector of size n */

int *ivector(long n);
/* Allocate a int vector of size n */
int **imatrix(int nr, int nc);
/* Allocate a int vector of size nr x nc */

int *izeros(long n);

long *lvector(long n);
/* Allocate a long int vector of size n */

long *lzeros(long n);

void free_dvector( double *v);
/* free a double vector allocated with dvector() */

void free_ivector( int *v);
/* free an integer vector allocated with ivector() */

void free_lvector( long *v);
/* free a long integer vector allocated with lvector() */

void print_ivector(int * v, int n);

void print_lvector(long * v, int n);


int choose(int n, int k);
double dchoose(int n, int k);
int isum(int *m, int N);
double dsum(double *m, int N);

iList_t * ilist_new();

void ilist_append(iList_t * list, int item);

void print_ilist(iList_t * list);

List_t * list_new();

void list_append(List_t * list, char * item);

int list_has_item(List_t * list, char * item);

void print_list(List_t * list);

void free_list( List_t * list);

int trim(char *s);
int remove_comments(char *s, char cmtsymbol);
int string_find(char * string, char **string_list, int list_size);

void report_error(const char *fmt, ...);
void report_warning(const char *fmt, ...);
/* Output error messages to stderr */
//void error(const char *fmt, ...)

#endif /* UTILS_H_ */
