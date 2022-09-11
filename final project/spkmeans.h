#include <stdio.h>
#include <stdlib.h>

double** make_mat(int n, int m);
double** k_mean(int k, double **data_points, double **centroids, int vec_num,
                int vec_size);                       
double*** main_by_goal(int k, int goal, double **matrix, int vec_num, int vec_size);
void print_mat_rows(double **mat, int n, int m);
void print_vec_row(double *vec, int n);
