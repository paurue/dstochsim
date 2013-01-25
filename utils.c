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

#include <stdarg.h>
#include "utils.h"

double *dvector(long n)
/* Allocate a double vector of size n */
{
	double *v;

	v = (double *) malloc( n * sizeof(double));
	if (!v) report_error("allocation failure in dvector()");
	return v;
}

double *dzeros(long n)
/* Allocate a int vector of size n and set it to zero*/
{
    int i;
    double *v;

    v = dvector(n);
    for(i=0; i<n; i++)
        v[i] = 0;
    return v;
}

int *ivector(long n)
/* Allocate a int vector of size n */
{
	int *v;

	v = (int *) malloc( n * sizeof(int));
	if (!v) report_error("allocation failure in ivector()");
	return v;
}

int **imatrix(int nr, int nc)
/* Allocate a int vector of size n */
{
    int i;
    int **m;

    m = (int **) calloc(nr, sizeof(int *));
    if (!m) {
        report_error("allocation failure in imatrix()");
        exit(1);
    }
    for(i=0; i<nr; i++) {
        m[i] = (int *) calloc(nc, sizeof(int));
        if (!m[i]) {
            report_error("allocation failure in imatrix()");
            exit(1);
        }
    }
    return m;
}


int *izeros(long n)
/* Allocate a int vector of size n and set it to zero*/
{
    int i,*v;

    v = ivector(n);
    for(i=0; i<n; i++)
        v[i] = 0;
    return v;
}

void print_ivector(int * v, int n){
    int i;
    printf("(");
    for(i=0;i<n-1;i++){
        printf("%d, ",v[i]);
    }
    printf("%d)\n", v[n]);
    return;
}

long *lvector(long n)
/* Allocate a int vector of size n */
{
	long *v;

	v = (long *) malloc( n * sizeof(long));
	if (!v) report_error("allocation failure in lvector()");
	return v;
}

long *lzeros(long n)
/* Allocate a int vector of size n and set it to zero*/
{
    int i;
    long *v;

    v = lvector(n);
    for(i=0; i<n; i++)
        v[i] = 0;
    return v;
}

void print_lvector(long * v, int n){
    int i;
    printf("(");
    for(i=0;i<n-1;i++){
        printf("%ld, ",v[i]);
    }
    printf("%ld)\n", v[n]);
    return;
}

int choose(int n, int k) {
	/* Binomial coefficient computation - adapted from Wikipedia */
	int i = 0;
    if (k > n)
        return 0;
    if (k > n/2)
        k = n-k; // Take advantage of symmetry
    double accum = 1;
    for (i = 1; i <= k; i++)
         accum = accum * (n-k+i) / i;

    return accum + 0.5; // avoid rounding error
}

int isum(int *m, int N){
    int i;
    int sum;
    sum = 0;
    for(i = 0; i < N; i++){
        sum += m[i];
    }
    return sum;
}

double dsum(double *m, int N){
    int i;
    double sum;
    sum = 0;
    for(i = 0; i < N; i++){
        sum += m[i];
    }
    return sum;
}


double dchoose(int n, int k) {
	/* Binomial coefficient computation - adapted from Wikipedia */
	int i = 0;
    if (k > n)
        return 0;
    if (k > n/2)
        k = n-k; // Take advantage of symmetry
    double accum = 1;
    for (i = 1; i <= k; i++)
         accum = accum * (n-k+i) / i;

    return accum; // avoid rounding error
}
void free_dvector( double *v)
/* free a double vector allocated with dvector() */
{
	free((char *) (v));
}

void free_ivector( int *v)
/* free an integer vector allocated with dvector() */
{
	free((char *) (v));
}

iList_t * ilist_new() {
	iList_t * l;

	l = (iList_t *) malloc(sizeof(iList_t));
	if (!l) report_error("allocation failure in ilist_new()");
	return l;
}

void ilist_append(iList_t * list, int item) {
	if(list->size < LIST_MAX_SIZE) {
		list->items[list->size] = item;
		list->size++;
	}
}

void print_ilist(iList_t * list) {
	int i;

	for(i=0; i<list->size; i++)
		printf("%d\n", list->items[i]);
}
List_t * list_new(){
	List_t * l;

	l = (List_t *) malloc(sizeof(List_t));
	if (!l) report_error("allocation failure in list()");
	return l;
}


void list_append(List_t * list, char * item) {
	if(list->size < LIST_MAX_SIZE) {
		strcpy(list->items[list->size], item);
		list->size++;
	}
}

int list_has_item(List_t * list, char * item) {
	int i, found;
	found = 0;
	for(i=0; i<list->size; i++){
		found = (strcmp(list->items[i], item) == 0);
		if(found) break;
	}

	if(found){
		return i;
	}
	else{
		return -1;
	}
}
void print_list(List_t * list)
{
	int i;

	for(i=0; i<list->size; i++)
		printf("%s\n", list->items[i]);
}

void free_list( List_t * list)
{
	free((char *) (list));
}

int trim(char *s) {
    int start, end, len;
    int i;
    len = strlen(s);
    start = 0;
    end = len - 1;
    while ((start < len) && (s[start] <= ' ')) {
        start++;
    }
    while ((start < end) && (s[end] <= ' ')) {
        end--;
    }
    if (start > end) {
        memset(s, '\0', len);
        return 0;
    }
    for (i = 0; (i + start) <= end; i++) {
        s[i] = s[start + i];
    }
    memset((s + i), '\0', len - i);
    return len;
}

int remove_comments(char *s, char cmtsymbol){
	int end, len;
	int i;
	len = strlen(s);
	end = len-1;
	for(i=end; i>=0; i--) {
		if(s[i] == cmtsymbol){
			s[i] = '\0';
			end = i;
			}
	}
	//memset((s + end + 1), '\0', len-end-1);
	return end;
}
int string_find(char * string, char **string_list, int list_size){
    int i, flag;
    flag = 0;
    for(i=0; i<list_size; i++){
        if(strcmp(string, string_list[i])==0){
            flag=1;
            break;
        }
    }
    if(flag) return i;
    return -1;
}

void report_error(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

void report_warning(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}
