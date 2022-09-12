#ifndef MAX_INT
#define MAX_INT 2147483647
#include "spkmeans.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* declerations */
double** make_mat(int n, int m);
double euclidean_dist(double *v1, double *v2, int n);
double ** weighted_adjecency_matrix(int n, double **data_points, int vec_size);
double* diagonal_degree_matrix(double **weighted_mat, int n);
double** matrix_multiplication(double **A, double **B, int size_A, int size_mid, int size_B);
double** normalized_graph_laplacian(double *diagonal_mat, double **weighted_mat, int n);
void matrix_sub(double **A, double **B, int n);
int is_diagonal(double **A, double **B, int n, double epsilon);
double** eye(int n);
double** make_P(int i, int j, double c, double s, int n);
int* pivot(double **sym_mat, int n);
double* obtain_c_t(double **sym_mat, int i, int j);
double square_off(double **mat, int n);
double*** jacobi_algorithm(double **sym, int n);
void* invalid_input();
/* int* read_data_points(FILE* file, double** data_points);
int read_sym_mat(FILE* file, double** sym_mat); */
double*** read_file(FILE *file);
void print_mat_rows(double **mat, int n, int m);
void print_mat_cols(double **mat, int n, int m);
void print_vec_row(double *vec, int n);
double** k_mean(int k, double **data_points, double **centroids, int vec_num,
                int vec_size);                       
void avg(double **cluster, int n, double *centroid, int vec_size);
void copy_mat(double** copy_to, double **copy_from, int n);
int eigengap_heuristic(double *eigenvalues, int n);
double*** main_by_goal(int k, int goal, double **matrix, int vec_num, int vec_size);
int comparator(const void *num1, const void *num2);
/*  */

/* works fine! */
double** make_mat(int n, int m)
{
    int i = 0;
    double **mat = (double **)malloc(n * sizeof(double*));
    assert(mat && "An Error Has Accured");
    for (; i < n; i++)
    {
        mat[i] = (double *)calloc(m , sizeof(double));
        assert(mat[i] && "An Error Has Accured");
    }
    return mat;
}


int comparator(const void *num1, const void *num2)
{
    double *n1, *n2;
    n1 = (double *)num1;
    n2 = (double *)num2;
    if(*n1 > *n2)
    {
        return -1;
    }
    if(*n1 < *n2)
    {
        return 1;
    }
    return 0;
}


double euclidean_dist(double *v1, double *v2, int n)
{
    double res = 0;
    int i = 0;
    for (; i < n; i++)
    {
        res += pow((v1[i] - v2[i]), 2);
    }
    res = pow(res , 0.5);
    return res;
}



double ** weighted_adjecency_matrix(int n, double **data_points, int vec_size)
{
    int i = 0, j;
    double **weighted_mat = make_mat(n, n);
    for (i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if( i == j)
            {
                weighted_mat[i][j] = 0;
            }
            else
            {
                weighted_mat[i][j] = exp(-(euclidean_dist(data_points[i],
                                             data_points[j], vec_size))/2);
            }
        }
    }
    return weighted_mat;
}


/* the matrix will be represented by a vector where the i'th element is the ii'th element 
and the rest are considered to be zeros */
double* diagonal_degree_matrix(double **weighted_mat, int n)
{
    int i = 0, z;
    double *ddm = (double *)calloc(n , sizeof(double));
    assert(ddm && "An Error Has Accured");
    for (; i < n; i++)
    {
        for (z = 0; z < n; z++)
        {
            ddm[i] += weighted_mat[i][z];
        }
    }
    return ddm;
}

/* works fine! */
double** matrix_multiplication(double **A, double **B, int size_A, int size_mid, int size_B)
{
    int row_a = 0, col_b, i;
    double **product = make_mat(size_A, size_B);
    for (; row_a < size_A; row_a++)
    {
        for (col_b = 0; col_b < size_B; col_b++)
        {
            for (i = 0; i < size_mid; i++)
            {
                product[row_a][col_b] += A[row_a][i] * B[i][col_b];
            }
        }
    }
    return product;
}


double** normalized_graph_laplacian(double *diagonal_mat, double **weighted_mat, int n)
{
    int i = 0;
    double **l_norm = eye(n);
    double **dwd, **dw;
    double **d = make_mat(n,n);
    for(; i < n; i++)
    {
        d[i][i] = pow(diagonal_mat[i], -0.5);
    }
/*     for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if (i == j)
            {
                d[i][j] = pow(diagonal_mat[i], -0.5);
            }
            else
            {
                d[i][j] = 0;
            }
        }   
    } */
/*     printf("d:\n");
    print_mat_rows(d, n, n); */
    dw = matrix_multiplication(d, weighted_mat, n, n, n);
/*     printf("\ndw:\n");
    print_mat_rows(dw, n, n); */
    dwd = matrix_multiplication(dw, d, n, n, n);
/*     printf("\ndwd:\n");
    print_mat_rows(dwd, n, n); */
    matrix_sub(l_norm, dwd, n);
    for(i = 0; i < n; i++)
    {
        free(d[i]);
        free(dw[i]);
        free(dwd[i]);
    }
    free(d);
    free(dw);
    free(dwd);  
    return l_norm;
}


/* performs A - B inplace */
void matrix_sub(double **A, double **B, int n)
{
    int i = 0, j;
    for(; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] -= B[i][j];
        }
    }
}

/* works fine! */
int is_diagonal(double **A, double **B, int n, double epsilon)
{
    return fabs(square_off(A, n) - square_off(B, n)) <= epsilon;
}

/* works fine! */
double** eye(int n)
{
    int i = 0;
    double **mat = make_mat(n, n);
    for(; i < n; i++)
    {
        mat[i][i] = 1; 
    }
    return mat;
}

/* same as Aviv */
double** make_P(int i, int j, double c, double s, int n)
{
    double **P = eye(n);
    P[i][i] = c;
    P[j][j] = c;
    P[i][j] = s;
    P[j][i] = -s;
    return P;
}


/* res = [i,j] */
int* pivot(double **sym_mat, int n)
{
    int i = 0, j;
    double max_val = -INFINITY;
    int *res = (int *)calloc(2 , sizeof(int));
    assert(res && "An Error Has Accured");
    res[0] = -1;
    res[1] = -1;
    for(; i < n; i++)
    {
        for(j = i + 1; j < n; j++)
        {
            if (i != j && fabs(sym_mat[i][j]) > max_val)
            {
                res[0] = i;
                res[1] = j;
                max_val = fabs(sym_mat[i][j]);
            }
        }
    }
    return res;
}


/* res = [c,t] */
double* obtain_c_t(double **sym_mat, int i, int j)
{
    double *res = (double *)calloc(2 , sizeof(double));
    double theta, t, c;
    int sign;
    theta = (sym_mat[j][j] - sym_mat[i][i]) / (2 * sym_mat[i][j]);
    /* sign */
    if (theta < 0) sign = -1;
    else sign = 1;

    t = sign / (fabs(theta) + pow(pow(theta, 2) + 1, 0.5));
    c = 1 / pow(pow(t, 2) + 1, 0.5);
    res[0] = c;
    res[1] = t;
    return res;
}

/* probably works fine */
double square_off(double **mat, int n)
{
    int i = 0, j;
    double sum_of_squares = 0;
    for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if (i != j)
            {
                sum_of_squares += pow(mat[i][j], 2);
            }
        }
    }
    return sum_of_squares;
}





void copy_mat(double** copy_to, double **copy_from, int n)
{   
    int i = 0, j;
    for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            copy_to[i][j] = copy_from[i][j];
        }
    }
}


/* with relation to the assignment A' will be referred to as B 
res = [eigenvalues, eigenvectors]
*/
double*** jacobi_algorithm(double **sym, int n)
{
    int i, j, iter = 0, r, idx, convergance;
    int *piv;
    double c, t, s;                          
    double *c_t;
    double **B = NULL, **A = NULL, **product, **temp_product, **P, **eigenvalues,
            **eigenvectors;
    double ***res;
    double epsilon = 1.0 *pow(10, -5);
    res = (double ***)malloc(2 * sizeof(double**));
    assert(res && "An Error Has Accured");
    eigenvalues = make_mat(1, n);
    product = eye(n);
    A = sym;
    B = make_mat(n,n);
    convergance = is_diagonal(A, B, n, epsilon) == 1;
    copy_mat(B, A, n);
    while(iter < 100 && convergance == 0)
    {
        piv = pivot(A, n);
        i = piv[0];
        j = piv[1];
        if(A[i][j] != 0)
        {
            c_t = obtain_c_t(A, i, j);
            c = c_t[0];
            t = c_t[1];
            s = c * t;
            free(c_t);
        }
        else
        {
            c = 1;
            s = 0;
        }
        /* calculating A' = P_TAP using formula from the assignment PDF */
        for(r = 0; r < n; r++)
        {
            if (r != i || r != j)
            {
                B[r][i] = c*A[r][i] - s*A[r][j];
                B[i][r] = c*A[r][i] - s*A[r][j];
                B[r][j] = c*A[r][j] + s*A[r][i];
                B[j][r] = c*A[r][j] + s*A[r][i];
            }
        }
        B[i][i] = pow(c, 2)*A[i][i] + pow(s, 2)*A[j][j] - 2*s*c*A[i][j];
        B[j][j] = pow(s, 2)*A[i][i] + pow(c, 2)*A[j][j] + 2*s*c*A[i][j];
        B[i][j] = 0;
        B[j][i] = 0;

        P = make_P(i, j, c, s, n);
        temp_product = matrix_multiplication(product, P, n, n, n);
        for(idx = 0; idx < n; idx++)
        {
            free(P[idx]);
        }
        free(P);
        free(product);
        product = temp_product;
        free(piv);
        iter++;
        convergance = is_diagonal(A, B, n, epsilon) == 1;
        copy_mat(A, B, n);

    } 
    for(idx = 0; idx < n; idx++)
    {
        free(B[idx]);
    }
    free(B);
    eigenvectors = (double **)calloc(n , sizeof(double *));
    assert(eigenvectors && "An Error Has Accured");
    for (idx = 0; idx < n; idx++)
    {
        eigenvalues[0][idx] = A[idx][idx];
        eigenvectors[idx] = product[idx];
    }
    free(A);
    free(product);
    res[0] = eigenvalues;
    res[1] = eigenvectors;
    return res;

}

/* the assignment says to go up to n/2 nםt sure i understood correctly */
int eigengap_heuristic(double *eigenvalues, int n)
{
    int i, k;
    double max_arg, deltha;
    max_arg = -INFINITY;
    qsort(eigenvalues, n, sizeof(double), comparator);
    for(i = 0; i <= (int)n/2; i++)
    {
        deltha = fabs(eigenvalues[i] - eigenvalues[i + 1]); 
        if(deltha > max_arg)
        {
            max_arg = deltha;
            k = i;
        }
    }
    return k + 1;
}

void* invalid_input()
{
    printf("Invalid Input!");
    exit(1);
}


/* res = [matrix, vec_size, vec_num] */

double*** read_file(FILE* file)
{
    char a;
    double num;
    int i = 0, vec_size_int = 0, vec_num_int = 0;
    double *vec;
    double **vec_size, **vec_num, **matrix;
    double ***res = (double***)malloc(3 * sizeof(double**));
    matrix = (double**)malloc(1000 * sizeof(double*));
    assert(matrix && "An Error Has Accured");
    vec = (double*)malloc(10 * sizeof(double));
    assert(vec && "An Error Has Accured");
    vec_size = make_mat(1,1);
    vec_num = make_mat(1,1);
    while(fscanf(file, "%lf", &num) == 1)
    {
        a = fgetc(file);
        if(vec_num_int > 0)
        {
            if(i == 0)
            {
                matrix[vec_num_int] = (double*)malloc(vec_size_int * sizeof(double));
            }
            matrix[vec_num_int][i++] = num;
            if(a == '\n')
            { 
                i = 0;
                vec_num_int++;
            }
        }
        /*this part is only for the first vector*/
        else
        {                
            vec[vec_size_int] = num;
            vec_size_int++;
            if(a != ',')
            {
                vec = (double*)realloc(vec, vec_size_int * sizeof(double));
                matrix[0] = vec;
                vec_num_int = 1;
            }
        }
        
    }
    fclose(file);
    matrix = (double**)realloc(matrix, vec_num_int * sizeof(double*));
    vec_size[0][0] = vec_size_int * 1.0;
    vec_num[0][0] = vec_num_int * 1.0;
    res[0] = matrix;
    res[1] = vec_size;
    res[2] = vec_num;
    return res;
}



void print_vec_row(double *vec, int n)
{
    int i = 0;
    for(; i < n - 1; i++)
    {
        printf("%.4lf,", vec[i]);
    }
    printf("%.4lf\n", vec[n - 1]);
}


void print_mat_rows(double **mat, int n, int m)
{
    int i = 0;
    for(; i < n; i++)
    {
        print_vec_row(mat[i], m);
    }
}


void print_mat_cols(double **mat, int n, int m)
{
    int i = 0, j;
    for(; i < n; i++)
    {
        for(j = 0; j < m - 1; j++)
        {
            printf("%.4lf,", mat[j][i]);
        }
        printf("%.4lf\n",mat[j][i]);
    }
}


double** k_mean(int k, double **data_points, double **centroids, int vec_num,
                int vec_size)
{
    int i, iteration=0, cluster_i = 0, idx, g, more_than_epsilon = 1;
    double min_euclidist, euclidist, norm, dist;
    double *change_vector, *zero;
    double **temp_centroids, **new_centroids;
    double ***clusters;
    int *clusters_sizes;
    new_centroids = make_mat(k , vec_size);
    assert(new_centroids && "An Error Has Accured");
    clusters = (double ***)malloc(k * sizeof(double**));
    assert(clusters && "An Error Has Accured");
    clusters_sizes = (int *)calloc(k, sizeof(int));
    assert(clusters_sizes && "An Error Has Accured");
    zero = (double*)calloc(vec_size, sizeof(double));
    assert(zero && "An Error Has Accured");
    for (i = 0 ; i < k ; i++)
    {
        clusters[i] = (double **)malloc((vec_num - k + 1) * sizeof(double *)); 
        /*largest cluster size can be num of data points - (k-1)*/
        assert(clusters[i] && "An Error Has Accured");
    }
    change_vector = (double *)malloc(k * sizeof(double));
    assert(change_vector && "An Error Has Accured");
    while (more_than_epsilon && iteration < 300)
    {
        iteration++;
        /* set cluster sizes to 0*/
        for (g = 0; g < k ; g++)
        {
            clusters_sizes[g] = 0;
        }
        /* make clusters */
        for (idx = 0 ; idx < vec_num ; idx++)
        {
            min_euclidist = (double)(MAX_INT);
            for (i = 0 ; i < k ; i++)
            {
                euclidist = euclidean_dist(data_points[idx], centroids[i], vec_size);
                if (euclidist < min_euclidist)
                {
                    min_euclidist = euclidist;
                    cluster_i = i;
                }
            }
            clusters[cluster_i][clusters_sizes[cluster_i]] = data_points[idx];
            clusters_sizes[cluster_i] += 1;
        }
        /* make centroids*/
        for (i = 0 ; i < k ; i++)
        {
            avg(clusters[i], clusters_sizes[i], new_centroids[i], vec_size);
        }
        /* make change vector:
        makes a sub vector of each two centroids and
        the norm of this sub is the cordinate in change vector*/
        for (i = 0 ; i < k ; i++)
        {
            norm = euclidean_dist(centroids[i], new_centroids[i], vec_size);
            change_vector[i] = norm;
        }
        dist = euclidean_dist(change_vector, zero, vec_size);
        if (dist < 0.0)
        {
            more_than_epsilon = 0;
        }

        temp_centroids = centroids;
        centroids = new_centroids;
        new_centroids = temp_centroids;
    }
    
    free(change_vector);
    for (i = 0 ; i < k ; i++)
    {
        free(clusters[i]);
    }
    free(centroids);
    free(clusters);
    free(clusters_sizes);
    printf("centroids in kmeans c is:\n");
    print_mat_rows(new_centroids, k, vec_size); 
    return new_centroids;
}


void avg(double **cluster, int n, double *centroid, int vec_size)
{
    int i, j;
    double sum;
    for (i = 0 ; i < vec_size ; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++) 
        {
            sum += cluster[j][i];
        }
        centroid[i] = sum/n;
    }
}

/* goals:
        1 = spk -> res = [T, NULL]
        2 = wam -> res= [weighted_mat, NULL]
        3 = ddg -> res = [ddm, NULL]
        4 = lnorm -> res = [laplace, NULL]
        5 = jacobi -> res = [eigenvalues, eigenvectors]
*/        
double*** main_by_goal(int k, int goal, double **matrix, int vec_num, int vec_size)
{
    int i, j, max_idx;
    double max_val, pow_sum;
    double *eigenvalues, *ddm;
    double **eigenvectors, **weighted_mat, **k_eigenvectors, 
            **full_ddm, **laplace, **U, **T;
    double ***res;
    res = (double ***)malloc(2 * sizeof(double **));
    assert(res && "An Error Has Accured\n");
    if(goal == 5)
    {
        res = jacobi_algorithm(matrix, vec_num);
    }
    else
    {
        res[0] = make_mat(1,1);
        weighted_mat = weighted_adjecency_matrix(vec_num, matrix, vec_size);
        if(goal == 2)
        /* goal is wam */
        {
            res[1] = weighted_mat;
        }
        else
        /* goal is ddg or lnorm */
        {
            ddm = diagonal_degree_matrix(weighted_mat, vec_num);
            if(goal == 3)
            /* goal is ddg */
            {
                /* this part makes a vector into a diagonal matrix */ 
                full_ddm = make_mat(vec_num, vec_num);
                for(i = 0; i < vec_num; i++)
                {
                    full_ddm[i][i] = ddm[i];
                }
                res[1] = full_ddm;
            }
            else 
            /* goal is lnorm or spk */
            {
                laplace = normalized_graph_laplacian(ddm, weighted_mat, vec_num);
                if(goal == 4)
                /* goal is lnorm */
                {
                    res[1] = laplace;
                }
                else
                /* goal is spk*/
                {
                    /* laplace is of size vec_num x vec_num */
                    res = jacobi_algorithm(laplace, vec_num);
                    eigenvalues = res[0][0];
                    eigenvectors = res[1];

                    /*works untill here depanding on the fact that jacobi,lnorm,ddg,wam is currect*/

                    /* determining k */
                    if(k == 0)
                    {
                        k = eigengap_heuristic(eigenvalues, vec_num);
                        res[0][0][0] = (double)k;
                        printf("k:%i\n",k);
                        exit(1);
                    }
                    k_eigenvectors = (double**)calloc(k, sizeof(double*));
                    assert(k_eigenvectors && "An Error Has Accured");
                    U = make_mat(vec_num, k);
                    /* finding the largest k eigenvalues */
                    for(i = 0; i < k; i++)
                    {
                        max_val = -INFINITY;
                        max_idx = -1;
                        for(j = 0; j < vec_num; j++)
                        {
                            if(eigenvalues[j] > max_val)
                            {
                                max_val = eigenvalues[j];
                                max_idx = j;
                            }
                        }
                        k_eigenvectors[i] = eigenvectors[max_idx];
                        eigenvalues[max_idx] = -INFINITY;
                    }
                    /* making U in transpose */
                    for(i = 0; i < k; i++)
                    {
                        for(j = 0; j < vec_num; j++)
                        {
                            U[j][i] = k_eigenvectors[i][j];
                        }
                    }
                    free(k_eigenvectors);
                    T = make_mat(vec_num, k);
                    for(i = 0; i < vec_num; i++)
                    {
                        pow_sum = 0;
                        for(j = 0; j < k; j++)
                        {
                            pow_sum += pow(U[i][j], 2);
                        }
                        pow_sum = sqrt(pow_sum);
                        for(j = 0; j < k; j++)
                        {
                            /* to avoid dividing by zero */
                            if(pow_sum > 0)
                            {
                                T[i][j] = U[i][j] / pow_sum;
                            } 
                            else T[i][j] = 0.0;           
                        }
                    }
                    res[1] = T;
                    /* fails to free */
/*                     for(i = 0; i < vec_num; i++)
                    {
                        free(U[i]);
                        free(eigenvectors[i]);
                        free(laplace[i]); 
                        free(weighted_mat[i]);
                    } */
                    free(U);
                    free(eigenvalues);
                    free(eigenvectors);
                    /* free(laplace); */ /* fails to free */
                    free(weighted_mat);
                }
            }
            /* ddm is a vector */
            free(ddm);
        }
    }   
    /* fails to free */
/*
    for(i = 0; i < vec_num; i++)
    {
        free(matrix[i]);
    }
    free(matrix); */
    return res;
}
/* reading input from CMD */
int main(int argc, char *argv[])
{
    int vec_size = 0, vec_num, i;
    char *goal;
    FILE *file;
    /* matrix is what we read of the file : data points / sym mat */
    double **matrix = NULL, **weighted_mat, **ddm, **laplace, **eigenvectors, **eigenvalues;
    double ***res = NULL;
    double ***data;
    if(argc != 3)
    {
        invalid_input();
    }
    file = fopen(argv[2], "r");
    if(file == NULL)
    {
        invalid_input();
    }
    goal = argv[1];
    if(strcmp(goal, "wam") == 0)
    {
        data = read_file(file);
        matrix = data[0];
        vec_size = (int)data[1][0][0];
        vec_num = (int)data[2][0][0];
        free(data[1][0]);
        free(data[1]);
        free(data[2][0]);
        free(data[2]);
        res = main_by_goal(-1, 2, matrix, vec_num, vec_size);
        weighted_mat = res[1];
        print_mat_rows(weighted_mat, vec_num, vec_num);
        for(i = 0; i < vec_num; i++)
        {
            free(weighted_mat[i]);
            free(matrix[i]);
        }
        free(matrix);
        free(weighted_mat);
        free(res);
        free(data);
    }
    else if(strcmp(goal, "ddg") == 0)
    {
        data = read_file(file);
        matrix = data[0];
        vec_size = (int)data[1][0][0];
        vec_num = (int)data[2][0][0];
        free(data[1][0]);
        free(data[1]);
        free(data[2][0]);
        free(data[2]);
        res = main_by_goal(-1, 3, matrix, vec_num, vec_size);
        ddm = res[1];
        print_mat_rows(ddm, vec_num, vec_num);
        for(i = 0; i < vec_num; i++)
        {
            free(ddm[i]);
            free(matrix[i]);
        }
        free(matrix);
        free(ddm);
        free(res);
        free(data);
    }
    else if(strcmp(goal, "lnorm") == 0)
    {
        data = read_file(file);
        matrix = data[0];
        vec_size = (int)data[1][0][0];
        vec_num = (int)data[2][0][0];
        free(data[1][0]);
        free(data[1]);
        free(data[2][0]);
        free(data[2]);
        res = main_by_goal(-1, 4, matrix, vec_num, vec_size);
        laplace = res[1];
        print_mat_rows(laplace, vec_num, vec_num);
        for(i = 0; i < vec_num; i++)
        {
            free(laplace[i]);
            free(matrix[i]);
        }
        free(matrix);
        free(laplace);
        free(res);
        free(data);
    } 

    else if(strcmp(goal, "jacobi")==0)
    {
        data = read_file(file);
        matrix = data[0];
        vec_size = (int)data[1][0][0];
        vec_num = (int)data[2][0][0];
        free(data[1][0]);
        free(data[1]);
        free(data[2][0]);
        free(data[2]);
        res = main_by_goal(-1, 5, matrix, vec_num, vec_size);
        eigenvalues = res[0];
        eigenvectors = res[1];
        print_vec_row(eigenvalues[0],vec_num);
        print_mat_rows(eigenvectors, vec_num, vec_num);
        /* problem here: */
/*         for(i = 0; i < vec_num; i++)
        {
            free(eigenvectors[i]);
            free(matrix[i]);
        }
        free(matrix); */
        free(eigenvalues[0]);
        free(eigenvalues);
        free(eigenvectors);
        free(res);
        free(data);
        
    }
    else 
    {
        printf("error - skip");
        invalid_input();
    }
    
}
#endif