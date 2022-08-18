#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* declerations */
static PyObject* kmeans_C(PyObject *self, PyObject *args);
static PyObject* C_part(PyObject *self, PyObject *args);
double** parse_py_table_to_c(PyObject *table_py, int n, int m);
PyObject* parse_c_table_to_py(double **table_c, int n, int m)
/*  */

static PyObject* C_part(PyObject *self, PyObject *args)
{
    int k, goal, vec_num, vec_size, i;
    PyObject *matrix_py, *T_py;
    double **matrix_c, **T_c;
    if(!PyArg_ParseTuple(args, "iiOii", &k, &goal, &matrix_py, &vec_num, &vec_size))
    {
        return NULL;
    }
    matrix_c = parse_py_table_to_c(matrix_py, vec_num, vec_size);
    T_c = main_by_goal(k, goal, matrix_c, vec_num, vec_size);
    for(i = 0; i < vec_num; i++)
    {
        free(matrix_c[i]);
    }
    free(matrix_c);
    if(T_c != NULL)
    /* goal is spk - go to kmeans++ */
    {
        T_py = parse_c_table_to_py(T_c, vec_num, k);
        free(T_c);
        return T_py;
    }
    return NULL;
}


static PyObject* kmeans_C(PyObject *self, PyObject *args)
{
    int k, vec_num, vec_size, i;
    PyObject *initial_centroids_py, *data_points_py, *centroids_py;
    double **data_points_c, **initial_centroids_c, **centroids_c;
    if(!PyArg_ParseTuple(args, "iOOii", &k, &data_points_py, &initial_centroids_py, &vec_num, &vec_size))
    {
        return NULL;
    }
    data_points_c = parse_py_table_to_c(data_points_py, vec_num, vec_size);
    initial_centroids_c = parse_py_table_to_c(initial_centroids_py, vec_num, vec_size);
    centroids_c = k_mean(k, data_points_c, initial_centroids_c, vec_num, vec_size);
    centroids_py = parse_c_table_to_py(centroids_c, vec_num, vec_size);
    for(i = 0; i < vec_num; i++)
    {
        free(data_points_c[i]);
    }
    free(data_points_c);
    for(i = 0; i < k; i++)
    {
        free(initial_centroids_c[i]);
        free(centroids_c[i]);
    }
    free(initial_centroids_c);
    free(centroids_c);
    return centroids_py;
}


/* static PyObject* heuristic_c(PyObject *self, PyObject *args)
{
    int k, vec_num, vec_size, i;
    PyObject *data_points_py;
    double **data_points_c, **weighted_mat, **diag_mat, **laplace, **eigenvectors;
    double *eigenvalues;
    if(!PyArg_ParseTuple(args, "Oii", &data_points_py, &vec_num, &vec_size))
    {
        return NULL;
    }
    eigenvalues = (double*)malloc(vec_num * sizeof(double));
    assert(eigenvalues && "An Error Has Accured");
    data_points_c = parse_py_table_to_c(data_points_py, vec_num, vec_size);
    weighted_mat = weighted_adjecency_matrix(data_points_c, vec_num, vec_size);
    diag_mat = diagonal_degree_matrix(weighted_mat, vec_num);
    laplace = normalized_graph_laplacian(diag_mat, weighted_mat, vec_num);
    eigenvectors = jacobi_algorithm(laplace, vec_num, eigenvalues);
    k = eigengap_heuristic(eigenvalues, vec_num);
    for(i = 0; i < vec_num; i++)
    {
        free(data_points_c[i]);
        free(weighted_mat);
        free(laplace[i]);
        free(eigenvectors[i]);
    }
    free(data_points_c);
    free(eigenvalues);
    free(weighted_mat);
    free(diag_mat);
    free(laplace);
    free(eigenvectors);
    return k;
}
 */
double** parse_py_table_to_c(PyObject *table_py, int n, int m)
{
    int i, j;
    double** table_c;
    table_c = make_mat(n, m);
    for (i = 0; i < n ; i++)
    {
        for (j = 0 ; j < m ; j++)
        {
            table_c[i][j] = PyFloat_AsDouble(PyList_GetItem(table_py, i * n + j));
        }
    }
    return table_c;
}


PyObject* parse_c_table_to_py(double **table_c, int n, int m)
{
    int i, j;
    PyObject *table_py;
    table_py = PyList_New(0);
    assert(table_py && "An Error Has Accured");
    for (i = 0; i < n ; i++)
    {
        for (j = 0 ; j < m ; j++)
        {
            PyList_Append(table_py, PyFloat_FromDouble(table_c[i][j]));
        }
    }
    return table_py;
}



static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, C_part, "all parts of C except kmeans"),
        FUNC(METH_VARARGS, kmeans_C, "kmeans part of C code"),
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _methods
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

#endif