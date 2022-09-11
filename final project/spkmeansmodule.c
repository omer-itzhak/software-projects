#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* declerations */
static PyObject* fit(PyObject *self, PyObject *args);
static PyObject* C_part(PyObject *self, PyObject *args);
double** parse_py_table_to_c(PyObject *mat_py, int n, int m);
PyObject* parse_c_table_to_py(double **table_c, int first_m, int n, int m);
void print_mat_py(PyObject *mat, int n, int m);
/*  */

static PyObject* C_part(PyObject *self, PyObject *args)
{
    int k, goal, vec_num, vec_size, i;
    PyObject *matrix_py, *res_py = Py_None;
    double *vec;
    double **matrix_c, **mat, **combined_res; 
    double ***res_c;
    if(!PyArg_ParseTuple(args, "iiOii", &k, &goal, &matrix_py, &vec_num, &vec_size))
    {
        return Py_None;
    }
    matrix_c = parse_py_table_to_c(matrix_py, vec_num, vec_size); 
    res_c = main_by_goal(k, goal, matrix_c, vec_num, vec_size);
    
    /* if goal = 5: */
    /* res_c = [eigenvalues, eigenvectors] */
    /* else: */
    /* res_c = [k, matrix] */

    vec = res_c[0][0];
    mat = res_c[1];
    combined_res = (double**)malloc((vec_num + 1) * sizeof(double*));
    assert(res && "An Error Has Accured");
    combined_res[0] = vec;
    for(i = 1; i < vec_num + 1; i++)
    {
        combined_res[i] = mat[i - 1];
    }
    if(goal == 5)
    {
        /* vec = eigenvalues -> size = 1 x vec_num
            mat = eigenvectors -> size = vec_num x vec_num */
        res_py = parse_c_table_to_py(combined_res, vec_num, vec_num, vec_num);
    }
    else if(goal == 1)
    {
        /* vec = k -> size = 1 x 1
            mat = T -> size = vec_num x k */
        res_py = parse_c_table_to_py(combined_res, 1, vec_num, k);
    }
    else
    {
        /* vec = k (= 0) -> size = 1 x 1
            mat = result of goal -> size = vec_num x vec_num */
        res_py = parse_c_table_to_py(combined_res, 1, vec_num, vec_num);
    }
    free(res_c[0]);
    /* fails to free */
/*     for(i = 0; i < vec_num + 1; i++)
    {
        free(res[i]);
    } */
    free(mat);
    free(combined_res);
    return Py_BuildValue("O", res_py); 
}


static PyObject* fit(PyObject *self, PyObject *args)
{
    int k, vec_num, vec_size;
    PyObject *initial_centroids_py, *data_points_py, *centroids_py;
    double **data_points_c, **initial_centroids_c, **centroids_c;
    if(!PyArg_ParseTuple(args, "iOOii", &k, &data_points_py, &initial_centroids_py, &vec_num, &vec_size))
    {
        return NULL;
    }
    printf("line 85 in fit");
    data_points_c = parse_py_table_to_c(data_points_py, vec_num, vec_size);
    printf("line 87 in fit");
    initial_centroids_c = parse_py_table_to_c(initial_centroids_py, vec_num, vec_size);
    printf("line 89 in fit");
    centroids_c = k_mean(k, data_points_c, initial_centroids_c, vec_num, vec_size);
    printf("line 91 in fit");
    centroids_py = parse_c_table_to_py(centroids_c, vec_size, k - 1, vec_size);
    printf("line 93 in fit\n");
    /* fails to free */
/*     for(i = 0; i < vec_num; i++)
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
    free(centroids_c); */
    printf("line 107 in fit\n");
    return Py_BuildValue("O", centroids_py);
}



PyObject* parse_c_table_to_py(double **table_c, int first_m, int n, int m)
{
    int i, j;
    PyObject *table_py = PyList_New(0), *vec_py;
    assert(table_py && "An Error Has Accured");
    vec_py = PyList_New(0);
    assert(vec_py && "An Error Has Accured");
    /* inserting first vec into py table */
    for(i = 0; i < first_m; i++)
    {
        PyList_Append(vec_py, PyFloat_FromDouble((double)table_c[0][i]));
    }
    PyList_Append(table_py, vec_py);
    /* inserting rest of vectors into py table */
    for (i = 1; i < n + 1; i++)
    {
        vec_py = PyList_New(0);
        assert(vec_py && "An Error Has Accured");
        for (j = 0 ; j < m ; j++)
        {
            PyList_Append(vec_py, PyFloat_FromDouble((double)table_c[i][j]));
        }
        PyList_Append(table_py, vec_py);
    }
    return table_py;
}

#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
        FUNC(METH_VARARGS, C_part, "all parts of C except kmeans"),
        FUNC(METH_VARARGS, fit, "kmeans part of C code"),
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


void print_mat_py(PyObject *mat, int n, int m)
{
    int i, j;
    for(i=0; i < n; i++)
    {
        for( j = 0; j < m ; j++)
        {
            printf("%.4lf,", PyFloat_AsDouble(PyList_GetItem(mat, (i * m) + j)));
        }
        printf("\n");
    }
}

double** parse_py_table_to_c(PyObject *mat_py, int n, int m)
{
    int i, j;
    double **mat_c = make_mat(n , m);
    double f;
    for(i=0; i < n; i++)
    {
        for( j = 0; j < m ; j++)
        {
            f =  PyFloat_AsDouble(PyList_GetItem(mat_py, (i * m) + j));
            mat_c[i][j] = f;
        }
    }
    return mat_c;
}



