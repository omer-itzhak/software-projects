#include <math.h>
#include <stdio.h>

/* declerations */
double** make_mat(int n, int m);
double euclidean_dist(double *v1, double *v2, int n);
double ** weighted_adjecency_matrix(int n, double **data_points);
double* diagonal_degree_matrix(double **weighted_mat, int n);
double** matrix_multiplication(double **product, double **A, double **B, int size_A, int size_mid, int size_B);
double** normalized_graph_laplacian(double *diagonal_mat, double **weighted_mat, int n);
void* matrix_sub(double **A, double **B, int n);
int is_diagonal(double **A, double **B, int n, double epsilon);
double** eye(int n);
double** make_P(double **P, double **sym_mat, int n,int i, int j, int c, int s);
int* pivot(double **sym_mat, int n);
double* obtain_c_t(double **sym_mat, int i, int j);
double square_off(double **mat, int n);
double*** jacobi_algorithm(double **sym, int n, double *eigenvalues);
double** eigengap_heuristic(double **l_norm );
/*  */



double** make_mat(int n, int m)
{
    int i = 0;
    double **mat = (double **)calloc(n * sizeof(double**));
    assert(mat);
    for (; i < m; i++)
    {
        mat[i] = (double *)calloc(m * sizeof(double));
        assert(mat[i]);
    }
    return mat;
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


/* the weighted matrix will only be full at the bottom triangle (where i > j)*/
double ** weighted_adjecency_matrix(int n, double **data_points)
{
    int i = 0, j;
    /* creating a lower triangle */
    double **weighted_mat = (double**)calloc(n * sizeof(double*));
    assert(weighted_mat);
    for (; i < n; i++)
    {
        weighted_mat[i] = (double*)calloc((i + 1) * sizeof(double));
        assert(weighted_mat[i]);
    }
    for (i = 0; i < n; i++)
    {
        for(j = 0; j < i; j++)
        {
            weighted_mat[i][j] = exp(-(euclidean_dist(data_points[i], data_points[j], n))/2);
        }
    }
    return weighted_mat;
}


/* the matrix will be represented by a vector where the i'th element is the ii'th element 
and the rest are considered to be zeros */
double* diagonal_degree_matrix(double **weighted_mat, int n)
{
    int i = 0, z;
    double *ddm = (double *)calloc(n * sizeof(double));
    assert(ddm);
    for (; i < n; i++)
    {
        for (z = 0; z < n; z++)
        {
            if (i > z)
            {
                ddm[i] += weighted_mat[i][z];
            }
            else
            {
                ddm[i] += weighted_mat[z][i];
            }
        }
    }
    return ddm;
}


double** matrix_multiplication(double **product, double **A, double **B, int size_A, int size_mid, int size_B)
{
    int row_a = 0, col_b, i;
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
    int i = 0, j;
    double **l_norm = eye(n);
    double **dwd, **dw;
    double **d = make_mat(n,n);
    dw = make_mat(n, n);
    dwd = make_mat(n, n);
    for(; i < n; i++)
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
    }
    dw = matrix_multiplication(dw, d, weighted_mat, n, n, n);
    dwd = matrix_multiplication(dwd, dw, d, n, n, n);
    free(d);
    free(dw);
    matrix_sub(l_norm, dwd, n);
    free(dwd);
    return l_norm;
}


/* performs A - B inplace */
void* matrix_sub(double **A, double **B, int n)
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


int is_diagonal(double **A, double **B, int n, double epsilon)
{
    return square_off(A, n) - square_off(B, n) <= epsilon;
}


double** eye(int n)
{
    int i = 0, j;
    double **mat = make_mat(n, n);
        for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            if (i == j)
            {
                mat[i][j] = 1;
            }
            else
            {
                mat[i][j] = 0;
            }
        }   
    }
    return eye;
}


double** make_P(double **P, double **sym_mat, int n,int i, int j, int c, int s)
{
    P[i][i] = c;
    P[j][j] = c;
    P[i][j] = s;
    P[j][i] = -s;
    return P;
}


/* res = [i,j] */
int* pivot(double **sym_mat, int n)
{
    double max_val = -INFINITY;
    double *res = (double *)calloc(2 * sizeof(double));
    assert(res);
    res[0] = -1;
    res[1] = -1;
    int i = 0, j;
    for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
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
    double *res = (double *)calloc(2 * sizeof(double));
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


/* double** transpose(double **mat_T, double **mat, int n)
{
    int i = 0, j;
    for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            mat_T[j][i] = mat[i][j];
        }
    }
    return mat_T;
} */


/* double** copy_mat(double** mat, int n)
{
    double **copy = make_mat(n, n);
    int i = 0, j;
    for(; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            copy[i][j] = mat[i][j];
        }
    }
    return copy;
} */


/* with relation to assignment A' will be referred to as B 
eigenvectors is etited inplace, res = matrix of eigenvectors*/
double*** jacobi_algorithm(double **sym, int n, double *eigenvalues)
{
    int i, j, iter = 0, r;
    int *piv;
    double c, t, s;
    double *c_t;
    double **B, **A, **product, **P;
    double epsilon = 1.0 *pow(10, -5);
    P = eye(n);
    //P_T = make_mat(n, n);
    product = eye(n);
    B = make_mat(n,n);
//    A = copy_mat(sym, n);
    A = sym;
    while(iter <= 100 || !is_diagonal(A, B, n, epsilon))
    {
        /* obtainning i, j, c, s */
        piv = pivot(A, n);
        i, j = piv[0], piv[1];
        free(piv);
        c_t = obtain_c_t(A, i, j);
        c, t = c_t[0], c_t[1];
        free(c_t);
        s = c * t;
        /* calculating A' = P_TAP using formula */
        for(r = 0; r < n; r++)
        {
            if (r != i && r != j)
            {
                B[r][i] = c*A[r][i] - s*A[r][j];
                B[i][r] = c*A[i][r] - s*A[j][r];
                B[r][j] = c*A[r][j] + s*A[r][i];
                B[j][r] = c*A[j][r] + s*A[i][r];
            }
        }
        B[i][i] = pow(c, 2)*A[i][i] + pow(s, 2)*A[j][j] - 2*s*c*A[i][j];
        B[j][j] = pow(s, 2)*A[i][i] + pow(c, 2)*A[j][j] + 2*s*c*A[i][j];
        B[i][j] = 0;
        B[j][i] = 0;
        /*  */
/*         P_T = transpose(P_T, P, n);
        temp_B = matrix_multiplication(temp_B, P_T, A, n, n, n);
        B = matrix_multiplication(B, temp_B, P, n, n, n);
        A = B; */
        P = make_P(P, A, n, i, j, c, s);
        product = matrix_multiplication(product, product, P, n, n, n);
        iter++;
    }
    free(P);
    free(B);
    //free(P_T);
    for (i = 0; i < n; i++)
    {
        eigenvalues[i] = A[i][i];
    }
    return product;
}

/* didnt understand :( */
double** eigengap_heuristic(double **l_norm ){}


