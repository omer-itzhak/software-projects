#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_INT 2147483647

void avg(double **cluster, int n, double *centroid, int size_vec);
void sub_vector(double *old, double *new, double *sub, int size_vec);
double euclidean_norm(double* vector, int size_vec);
int is_numeric(char *str);
static double** k_mean(int k, int max_iter, double epsilon, double **data_points_c,double **initial_centroids_c,int total_vec_number,int size_vec);
double** create_mat(int vec_num, int vec_size);


static double** k_mean(int k, int max_iter,double epsilon, double ** data_points, double **centroids,int total_vec_number,int size_vec)
{
    int i, iteration=0, cluster_i = 0, idx, g, more_than_epsilon = 1;
    double min_euclidean_dist, euclidean_dist, norm,dist;
    double *change_vector, *sub;
    double **temp_centroids, **new_centroids;
    double ***clusters;
    int *clusters_sizes;
    new_centroids = create_mat(k , size_vec);
    if (!new_centroids)
    {
        return NULL;
    }
    sub = (double *)calloc(size_vec, sizeof(double));
    if (!sub)
    {
        return NULL;
    }
    clusters = (double ***)calloc(k,sizeof(double**));
    if (!clusters)
    {
        return NULL;
    }
    clusters_sizes = (int *)malloc(k*sizeof(int));
    if(!clusters_sizes)
    {
        return NULL;
    }
    for (i = 0 ; i < k ; i++)
    {
        clusters[i] = (double **)malloc((total_vec_number - k + 1)*sizeof(double *)); 
        /*largest cluster size can be num of data points - (k-1)*/
        if (!clusters[i])
        {
            return  NULL;
        }
    }
    change_vector = (double *)malloc(k*sizeof(double));
    if(!change_vector)
    {
        return NULL;
    }


    while (more_than_epsilon && iteration < max_iter)
    {
        iteration += 1;

        /* set cluster sizes to 0*/

        for (g = 0; g < k ; g++)
        {
            clusters_sizes[g] = 0;
        }

        /* make clusters */

        for (idx = 0 ; idx < total_vec_number ; idx++)
        {
            min_euclidean_dist = (double)(MAX_INT);
            for (i = 0 ; i < k ; i++)
            {
                sub_vector(data_points[idx], centroids[i], sub, size_vec);
                euclidean_dist = euclidean_norm(sub, size_vec);
                if (euclidean_dist < min_euclidean_dist)
                {
                    min_euclidean_dist = euclidean_dist;
                    cluster_i = i;
                }
            }
            clusters[cluster_i][clusters_sizes[cluster_i]] = data_points[idx];
            clusters_sizes[cluster_i] += 1;
        }

        /* make centroids*/
        for (i = 0 ; i < k ; i++)
        {
            avg(clusters[i], clusters_sizes[i], new_centroids[i], size_vec);
        }

        /* make change vector:
        makes a sub vector of each two centroids and
        the norm of this sub is the cordinate in change vector*/
        for (i = 0 ; i < k ; i++)
        {
            sub_vector(centroids[i], new_centroids[i], sub, size_vec);
            norm = euclidean_norm(sub, size_vec);
            change_vector[i] = norm;
        }
        dist = euclidean_norm(change_vector, size_vec);
        if (dist < epsilon)
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
    free(sub);
    free(clusters_sizes);
    return new_centroids;
}
void avg(double **cluster, int n, double *centroid, int size_vec)
{
    int i, j;
    double sum;
    for (i = 0 ; i < size_vec ; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++) 
        {
            sum += cluster[j][i];
        }
        centroid[i] = sum/n;
    }
}

void sub_vector(double *old, double *new, double *sub, int size_vec)
{
    int i;
    for(i = 0 ; i < size_vec ; i++)
    {
        sub[i] = old[i] - new[i];
    }
}

double euclidean_norm(double* vector, int size_vec)
{
    int i;
    double res = 0;
    for (i = 0 ; i < size_vec ; i++){
        res += pow(vector[i],2);
    }
    res = pow(res, 0.5);
    return res;
}
int is_numeric(char *str)
{
    int i = 0;
    while (str[i] != 0){
        if(str[i] < 48 || str[i] > 57){
            return 0;
        }
        i++;
    }
    return 1;
}


double** parse_py_table_to_C(PyObject *lst, int vec_num, int vec_size)
{
    int row, col;
    double **data_points_c = create_mat(vec_num, vec_size);
    if (!data_points_c)
    {
        return NULL;
    }
    for (row = 0; row < vec_num ; row++)
    {
        for (col = 0 ; col < vec_size ; col++)
        {
             data_points_c[row][col] = PyFloat_AsDouble(PyList_GetItem(lst, row * vec_size + col));
        }
    }
    return data_points_c;
    
}
double** create_mat(int vec_num, int vec_size){
    int i=0;
    double **mat;
    mat = calloc(vec_num, sizeof(double*));
    if (!mat)
    {
        return NULL;
    }
    for (; i < vec_num; i++)
    {
        mat[i] = calloc(vec_size, sizeof(double));
        if (!mat[i])
        {
            return NULL;
        }
    }
    return mat;
}
// double** create_mat(int vec_num, int vec_size){
//     int i;
//     double *vec;
//     double **mat;
//     vec = calloc(vec_num*vec_size, sizeof(double));
//     if (!vec){
//         return NULL;
//     }
//     mat = calloc(vec_num, sizeof(double*));
//     if (!mat){
//         return NULL;
//     }
//     for(i=0; i<vec_num; i++){
//         mat[i] = vec + i * vec_size;
//     }
//     return mat;
// }

static PyObject* fit(PyObject *self, PyObject *args)
{
    int k, max_iter, total_vec_number, size_vec, i, j;
    double epsilon;
    double **data_points_c, **initial_centroids_c, **centroids_c;
    PyObject *data_points_py, *initial_centroids_py, *centroids_py, *centroid_py;
    if (!PyArg_ParseTuple(args, "iidOOii", &k, &max_iter, &epsilon, &data_points_py, &initial_centroids_py, &total_vec_number, &size_vec))
    {
        return NULL;
    }
    data_points_c = parse_py_table_to_C(data_points_py, total_vec_number, size_vec);
    initial_centroids_c = parse_py_table_to_C(initial_centroids_py, k, size_vec);
    centroids_c = k_mean(k, max_iter, epsilon, data_points_c, initial_centroids_c, total_vec_number, size_vec); // TODO: change
    centroids_py = PyList_New(0);
    /* ask if necessary*/
    if (!centroids_py)
    {
        return NULL;
    }
    for (i = 0; i < k ; i++)
    {
        centroid_py = PyList_New(0);
        if (centroid_py == NULL)
        {
            return NULL;
        }
        for (j = 0 ; j < size_vec ; j++)
        {
            PyList_Append(centroid_py, PyFloat_FromDouble(centroids_c[i][j]));
        }
        PyList_Append(centroids_py, centroid_py);
    }

    /*free vecs: */
   /*  free(data_points_c[0]); */
/* problem with this memory free loop? */
    for (i = 0; i < total_vec_number ; i++)
    {        
        free(data_points_c[i]);
    }
    free(data_points_c); 
    free(centroids_c);
    return centroids_py;
}


/* methods from class*/

static PyMethodDef capiMethods[] = {
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("returns calculated centroids for given matrix of data points")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}