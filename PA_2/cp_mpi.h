#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>


/*
 * Function Declarations
 */


int get_row_proc(int);


void get_matrix_size(char *);


void get_size_of_line(int *);


int get_size_of_buffer(int);


int get_limit_of_row(int);


double * user_matrix_read(char *, int);


double* user_matrix_load();


void print_matrix(double *, int, int);


void print(double *);


double* aloc_matrix(int);


void calculate_RREF(double *, int);


void calculate_RREF_p(double *, int);


void calculate_REF(double *, int, int);


double get_local_maximum(double *, int);


double get_global_maximum(double);


void divide_by_max(double *, int, double);


void write_file(double *, int);


void print_best_acceptance_threshold(double *, int);
