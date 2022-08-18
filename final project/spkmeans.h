#include <stdio.h>
#include <stdlib.h>

double** main_by_goal(int k, int goal, double **matrix, int vec_num, int vec_size);
static double** k_mean(int k, double **data_points, double **centroids,
                        int vec_num,int vec_size);
double** make_mat(int n, int m);
