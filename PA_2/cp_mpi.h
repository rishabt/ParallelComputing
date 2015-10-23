#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

void get_matrix_size(char *);
void get_size_of_line(int *);
int get_size_of_buffer(int);
int get_limit_of_row(int);
double * user_matrix_read(char *, int);
double * user_matrix_load();
void print_matrix(double *, int, int);
void print(double *);
void calculate_RREF(double *, int);
void calculate_RREF_p(double *, int);
void calculate_REF(double *, int, int);
double find_local_max(double *, int);
double find_global_max(double);
void divide_by_max(double *, int, double);
void write_clicking_probabilities(double *, int);
void print_best_acceptance_threshold(double *, int);
