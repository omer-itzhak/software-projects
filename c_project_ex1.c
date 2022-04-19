#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

double* avg(double **cluster, int cluster_i);
double** read_file(char *input_filename);
double* sub_vector(double *old, double *new);
double** read_file(char *input_filename);
double* sub_vector(double *old, double *new);
double* sub_vector(double *old, double *new);
double** read_file(char *input_filename);
double* sub_vector(double *old, double *new);
double euclidean_norm(double* vector);
int is_numeric(char *str);
int k_mean(int k, int max_iter, char *input_filename,char *output_filename);
int invalid_input();
int write_to_file(char *output_filename, double **centroids, int k);
int an_error();
int size_vec,total_vec_number;
int *clusters_sizes;

void print_vector(double *vector)
{
    int g;
    for (g=0; g<size_vec; g++)
    {
        printf("%Lf,", vector[g]);
    }
    printf("\n"); 
}

int main(int argc,char *argv[])
{
    FILE* opf;
    int k, max_iter, idx_of_file;
    char *input_filename, *output_filename;
    if (argc < 4 || argc > 5)
    {
        invalid_input();
    }
    if (is_numeric(argv[1]) == 1)
    {
        k = atoi(argv[1]);
    }
    else
    {
        invalid_input();
    }
    if (argc == 5)
    {
        if (is_numeric(argv[2]) == 1)
        {
            max_iter = atoi(argv[2]);
        }
        else
        {
            invalid_input();
        }
        idx_of_file = 3;
    } 
    else /* argc = 5*/
    {
        max_iter = 200;
        idx_of_file = 2;
    }
    
    opf = fopen(argv[idx_of_file + 1], "w");
    if (opf == NULL)
    {
        invalid_input();
    }
    fclose(opf);
    input_filename = argv[idx_of_file];
    output_filename = argv[idx_of_file + 1];
    if (k <= 0 || max_iter <= 0)
    {
        invalid_input();
    }
    return k_mean(k,max_iter,input_filename,output_filename);
}
int k_mean(int k, int max_iter, char *input_filename,char *output_filename)
{
    int i,iteration=0,cluster_i = 0, idx, g;
    char more_than_epsilon;
    double epsilon,min_euclidean_dist,euclidean_dist,norm,dist;
    double *data_point,*change_vector,*sub;
    double **data_points,**centroids,**new_centroids, **temp_centroids;
    double ***clusters;
    more_than_epsilon = 1;
    epsilon = 0.001;
    data_points = read_file(input_filename);
    centroids = (double **)calloc(k,sizeof(double*));
    new_centroids = (double **)calloc(k,sizeof(double*));
    assert(centroids && "An Error Has Occurred");
    assert(new_centroids && "An Error Has Occurred");
    for (i = 0 ; i < k ; i++)
    {
        centroids[i] = data_points[i];
    }
    clusters = (double ***)calloc(k,sizeof(double**));
    assert(clusters && "An Error Has Occurred");
    clusters_sizes = (int *)malloc(k*sizeof(int));
    assert(clusters_sizes && "An Error Has Occurred");
    for (i = 0 ; i < k ; i++)
    {
        double **cluster = (double **)malloc((total_vec_number - k + 1)*sizeof(double *)); 
        /*largest cluster size can be num of data points - (k-1)*/
        assert(cluster != NULL && "An Error Has Occurred");
        clusters[i] = cluster;
    }
    change_vector = (double *)malloc(k*sizeof(double));
    assert(change_vector && "An Error Has Occurred");


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
            data_point = data_points[idx];

            min_euclidean_dist = (double)(INT_MAX);
            for (i = 0 ; i < k ; i++)
            {
                euclidean_dist = euclidean_norm(sub_vector(data_point, centroids[i]));
                if (euclidean_dist < min_euclidean_dist)
                {
                    min_euclidean_dist = euclidean_dist;
                    cluster_i = i;
                }
            }
            clusters[cluster_i][clusters_sizes[cluster_i]] = data_point;
            clusters_sizes[cluster_i] += 1;
        }

        /* make centroids*/
        for (i = 0 ; i < k ; i++)
        {
            new_centroids[i] = avg(clusters[i], i);
        }

        /* make change vector*/
        /* makes a sub vector of each two centroids and
        the norm of this sub is the cordinate in change vector*/
        for (i = 0 ; i < k ; i++)
        {
            sub = sub_vector(centroids[i], new_centroids[i]);
            norm = euclidean_norm(sub);
            free(sub);
            change_vector[i] = norm;
        }
        dist = euclidean_norm(change_vector);
        if (dist < epsilon)
        {
            more_than_epsilon = 0;
        }

        temp_centroids = centroids;
        centroids = new_centroids;
        new_centroids = temp_centroids;
        for (g = 0 ; g < k ; g++)
        {
            free(new_centroids[g]);
        }

    }
    for (i=0;i<k;i++)
    {
        printf("clusters_sizes[%d] = %d\n",i,clusters_sizes[i]);
        printf("cluster[%d]=\n",i);
        for (g=0;g<clusters_sizes[i];g++)
        {
            print_vector(clusters[i][g]);
        }
    }
    assert(write_to_file(output_filename, centroids, k) == 0 && "An Error Has Occurred");
    free(change_vector);
    for (i = 0; i < total_vec_number ; i++)
    {
        free(data_points[i]);
    }
    free(data_points);
    for (i = 0 ; i < k ; i++)
    {
        free(centroids[i]);
        free(clusters[i]);
    }
    free(centroids);
    free(new_centroids);
    free(temp_centroids);
    free(clusters);
    return 0;
}
double* avg(double **cluster, int cluster_i)
{
    int i, j, n;
    double sum;
    double *centroid;
    n = clusters_sizes[cluster_i];/* maybe*/
    centroid = (double *)malloc(size_vec*sizeof(double));
    assert(centroid && "An Error Has Occurred");
    for (i = 0 ; i < size_vec ; i++)
    {
        sum = 0;
        for (j = 0; j < n; j++) 
        {
            sum += cluster[j][i];
        }
        centroid[i] = sum/n;
    }
    return centroid;
}
double** read_file(char *input_filename)
{  
    int index=0, line=0;
    double a;
    char b;
    double **data_points;
    FILE *ipf;

    /* count cordinate number in vector */

    ipf = fopen(input_filename, "r");
    assert (ipf != NULL && "An Error Has Occurred");
    size_vec = 0;
    while ((fscanf(ipf,"%lf",&a) == 1)){
        b = fgetc(ipf);
        size_vec++;
        if( b == '\n'){
            break;
        }
    }
    fclose(ipf);

    /* count total number of vectors */

    ipf = fopen(input_filename, "r");
    assert (ipf != NULL && "An Error Has Occurred");
    total_vec_number = 0;
    for(b = fgetc(ipf) ; b != EOF ; b = fgetc(ipf)){
        if(b == '\n'){
            total_vec_number += 1;
        }
    }
    fclose(ipf);

    /* main read section */

    data_points = (double **)calloc(total_vec_number,sizeof(double*));
    assert(data_points && "An Error Has Occurred");
    ipf = fopen(input_filename, "r");
    assert (ipf != NULL && "An Error Has Occurred");
    for(line = 0 ; line < total_vec_number ; line++)
    { 
        double *vector = calloc(size_vec,sizeof(double));
        assert(vector && "An Error Has Occurred");
        for(index = 0 ; index < size_vec ; index++)
        {    
            fscanf(ipf,"%lf",&a);
            vector[index] = a;
            fgetc(ipf);
        }
        data_points[line] = vector;
        printf("data_points[%d] = ",line);
        print_vector(data_points[line]);
    }
    fclose(ipf);  
    return data_points;
}
double* sub_vector(double *old, double *new)
{
    int i;
    double *sub = (double *)calloc(size_vec,sizeof(double));
    assert(sub && "An Error Has Occurred");
    for(i = 0 ; i < size_vec ; i++)
    {
        sub[i] = old[i] - new[i];
    }
    return sub;
}
double euclidean_norm(double* vector)
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
int invalid_input()
{
    printf("Invalid Input!");
    exit(1);
}
int an_error()
{
    printf("An Error Has Occurred");
    exit(1);
}
int write_to_file(char *output_filename, double **centroids, int k)
{
    FILE* opf;
    int j,i;
    opf = fopen(output_filename, "w");
    assert(opf != NULL && "An Error Has Occurred");
    for (j = 0 ; j < k ; j++)
    {
        for (i = 0 ; i < size_vec; i++)
        {
            fprintf(opf,"%.4f" , centroids[j][i]);
            if (i != size_vec - 1)
            {
                fprintf(opf,"%s",",");
            }
        }
        fprintf(opf,"\n");
    }
    fclose(opf);
    return 0;
}
