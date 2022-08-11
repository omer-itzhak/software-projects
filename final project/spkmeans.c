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
double ** weighted_adjecency_matrix(int n, double **data_points, int vec_size)
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
            weighted_mat[i][j] = exp(-(euclidean_dist(data_points[i], data_points[j], vec_size))/2);
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
            /* if i == z then weighted_mat[i][z] doesn't exist
            and is considered to be 0 */
            if (i > z)
            {
                ddm[i] += weighted_mat[i][z];
            }
            else if (i != z)
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
    int i = 0;
    double **mat = make_mat(n, n);
    for(; i < n; i++)
    {
        mat[i][i] = 1; 
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


/* with relation to the assignment A' will be referred to as B 
eigenvectors is edited inplace, res = matrix of eigenvectors*/
double*** jacobi_algorithm(double **sym, int n, double *eigenvalues)
{
    int i, j, iter = 0, r, cnt;
    int *piv;
    double c, t, s;
    double *c_t, *temp_eigenvalues;
    double **B, **A, **product, **P, **res;
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
    temp_eigenvalues = (double *)calloc(n * sizeof(double));
    res = (double **)calloc(n * sizeof(double *));
    for (i = 0; i < n; i++)
    {
        temp_eigenvalues[i] = A[i][i];
    }
    free(A);
    /* sort eigenvalues and corresponding eigen vectors */
    for(i = 0; i < n; i++)
    {
        cnt = 0;
        for(j = 0; j < n; j++)
        {
            if (temp_eigenvalues[i] <= temp_eigenvalues[j]) cnt++;
        }
        eigenvalues[cnt] = temp_eigenvalues[i];
        res[cnt] = product[i];
    }
    free(product);
    free(temp_eigenvalues);
    return res;
}

/* the assignment says to go up to n/2 nut sure i understood correctly */
double** eigengap_heuristic(double *eigenvalues, int n)
{
    int i, k;
    double max_arg = -INFINITY, deltha;
    for(i = 0; i <= n/2; i++)
    {
        deltha = eigenvalues[i] - eigenvalues[i + 1]; 
        if(deltha > max_arg)
        {
            max_arg = deltha;
            k = i;
        }
    }
    return k;
}


void* invalid_input()
{
    printf("Invalid Input!");
    exit(1);
}

/* reading input from CMD */
int main(int argc, char *argv[])
{
    int i, n, vec_size, vec_num, j, g = -1;
    char *goal;
    FILE *file;
    double **sym_mat, **data_points, **eigenvectors, **weighted_mat, 
            **full_weighted_mat, *ddm, **full_ddm, **laplace;
    double *eigenvalues;
    int *size_num;
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
    if(strcmp(goal, "jacobi"))
    {
        n = read_sym_mat(file, sym_mat);
        eigenvalues = (double*)malloc(n, sizeof(double));
        assert(eigenvalues && "An Error Has Accured");
        eigenvectors = jacobi_algorithm(sym_mat, n, eigenvalues);
        print_vec_row(eigenvalues, n);
        print_mat_cols(eigenvectors, n, n);
        for(i = 0; i < n; i++)
        {
            free(sym_mat[i]);
        }
        free(sym_mat);
        free(eigenvalues);
        free(eigenvectors);
    }
    else
    {
        if(strcmp(goal, "wam") == 0) g = 0;
        else if(strcmp(goal, "ddg") == 0) g = 1;
        else if(strcmp(goal, "lnorm") == 0) g = 2;
        else invalid_input();
        size_num = read_data_points(file, data_points);
        vec_size, vec_num = size_num[0], size_num[1];
        weighted_mat = weighted_adjecency_matrix(vec_num, data_points, vec_size);
        if(g == 0)
        /* goal is wam */
        {
            full_weighted_mat = make_mat(vec_num, vec_num);
                for (i = 0; i < n; i++)
                {
                    full_weighted_mat[i][i] = 0;
                    for (j = 0; j < n; j++)
                    {
                        /* if i == j then weighted_mat[i][j] doesn't exist
                        and is considered to be 0 */
                        if (i > j)
                        {
                            full_weighted_mat[i][j] = weighted_mat[i][j];
                        }
                        else if (i != j)
                        {
                            full_weighted_mat[i][j] = weighted_mat[j][i];
                        }
                    }
                }
            print_mat_rows(full_weighted_mat, n, n);
            free(full_weighted_mat);
        }
        else
        /* goal is ddg or lnorm */
        {
            ddm = diagonal_degree_matrix(weighted_mat, vec_num);
            if(g == 1)
            /* goal is ddg */
            {
                full_ddm = (double**)calloc(vec_num, sizeof(double*));
                for(i = 0; i < n; i++) 
                {
                    full_ddm[i][i] = ddm[i];
                }
                print_mat_rows(full_ddm, vec_num, vec_num);
                free(full_ddm);
            }
            else
            /* goal is lnorm */
            {
                laplace = normalized_graph_laplacian(ddm, weighted_mat, vec_num);
                print_mat_rows(laplace, vec_num, vec_num);
                free(laplace);
            }
            free(ddm);
        }
        free(weighted_mat);
        for(i = 0; i < vec_num; i++)
        {
            free(data_points[i]);
        }
        free(data_points);
    }
    return 0;
}

/* res = [vec_size, vec_num] */
int* read_data_points(FILE* file, double** data_points)
{
    char a;
    double num;
    int i = 0, vec_size = 0, vec_num = 0, size_of_data_points = 1024;
    double *vec;
    int *res = (int*)malloc(2, sizeof(int));
    data_points = (double**)malloc(size_of_data_points, sizeof(double*));
    assert(data_points && "An Error Has Accured");
    vec = (double*)malloc(size_of_data_points, sizeof(double));
    assert(vec && "An Error Has Accured");
    while(fscanf(file, "%lf", num) && fgetc(file, &a))
    {
        if(vec_num > 0)
        {
            if(vec_num > size_of_data_points)
            {
                data_points = (double**)realloc(data_points, 2 * size_of_data_points * sizeof(double*));
                size_of_data_points *= 2;
            }
            if(i == 0)
            {
                data_points[vec_num] = (double*)malloc(vec_size * sizeof(double));
            }
            data_points[vec_num][i++] = num;
            if(!a)
            {
                i = 0;
                vec_num++;
            }
        }
        /* this part is only for the first vector */
        else
        {
            vec[vec_size] = num;
            if(a)
            {
                vec_size++;
            }
            else
            {
                vec = (double*)realloc(vec, vec_size * sizeof(double));
                data_points[0] = vec;
                vec_num = 1;
            }
        }
    }
    fclose(file);
    data_points = (double**)realloc(data_points, vec_num * sizeof(double*));
    res[0] = vec_size;
    res[1] = vec_num;
    return res;
}


/* res = n */
double** read_sym_mat(FILE* file, double** sym_mat)
{
    int n = 0, row = 0, col, size_of_mat = 1024;
    double num;
    double *vec;
    char a;
    sym_mat = (double**)malloc(size_of_mat, sizeof(double*));
    assert(sym_mat && "An Error Has Accured");
    vec = (double*)malloc(size_of_mat, sizeof(double));
    assert(vec && "An Error Has Accured");
    while(fscanf(file, "%lf", num) && fgetc(file, &a) && !row)
    {
        if(n > size_of_mat)
        {
            vec = (double*)realloc(vec, size_of_mat * 2 * sizeof(double));
            assert(vec && "An Error Has Accured");
            size_of_mat *= 2;
        }
        vec[n] = num;
        if(a) n++;
        else row = 1;

    }
    sym_mat[0] = (double*)realloc(vec, n * sizeof(double));
    assert(sym_mat[0] && "An Error Has Accured");
    sym_mat = (double**)realloc(sym_mat, n * sizeof(double*));
    assert(sym_mat && "An Error Has Accured");
    for(; row < n; row++)
    {
        sym_mat[row] = (double*)malloc(n , sizeof(double));
        assert(sym_mat && "An Error Has Accured");
        for(col = 0; col < n; col++)
        {
            sym_mat[row][col] = fscanf(file, "%lf");
            fgetc(file);
        }
    }    
    fclose(file);
    return n;
}


void* print_vec_row(double *vec, int n)
{
    int i = 0;
    for(; i < n - 1; i++)
    {
        printf("%lf,", vec[i]);
    }
    printf("%lf\n", vec[n - 1]);
}


print_mat_rows(double **mat, int n, int m)
{
    int i = 0;
    for(; i < n; i++)
    {
        print_vec_row(mat[i], m);
    }
}


print_mat_cols(double **mat, int n, int m)
{
    int i = 0, j;
    for(; i < n; i++)
    {
        for(j = 0; j < m; j++)
        {
            printf("%lf ", mat[j][i]);
        }
        printf("\n");
    }
}


