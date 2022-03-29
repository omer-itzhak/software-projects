#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


double* avg(double **cluster)
{
    int res = 0;
    int n = sizeof(cluster[0]);
    double *centroid = (double *)malloc(n*sizeof(double));
    for (int i=0; i<n;i++)
    {
        double temp = 0;
        for (int j=0; j<sizeof(cluster); j++)
        {
            temp += cluster[j][i];
        }
        centroid[i] = temp/sizeof(cluster);
    }
    return centroid;
}
double** read_file(char *input_filename)
{  
    double **temp_data_points = (double **)malloc(1024*sizeof(double));
    int cnt_coordinates=0, cnt_lines=0;
    double* temp_vector = (double *)malloc(128*sizeof(double));
    if(temp_vector != NULL && temp_data_points != NULL)
    {
        FILE *ipf = NULL;
        ipf = fopen(input_filename, "r");
        if (ipf != NULL)
        {
            double a;
            char b;
            while((fscanf(ipf,"%f%c",&a,&b))!=EOF)
            {    
                temp_vector[cnt_coordinates] = a;
                if (b == ',')
                {
                    cnt_coordinates++;
                }
                else{
                    if (b == '\n')
                    {
                        temp_data_points[cnt_lines] = temp_vector;
                        cnt_coordinates = 0;
                        cnt_lines++;
                    }   
                }             
                
            }
        }
        
    }
    double** data_points = (double**)malloc(cnt_lines*sizeof(double*));
    for (int i=0; i<cnt_lines; i++)
    {
        data_points[i] = temp_data_points[i];
    }
    for (int i=0; i<1024; i++){
        free(temp_data_points[i]);
    }
    free(temp_data_points);
    free(temp_vector);
    return data_points;
}



double* sub_vector(double *old, double *new)
{
    double *sub = (double *)malloc(sizeof(old)*sizeof(double));
    for(int i=0; i<sizeof(old); i++)
    {
        sub[i] = old[i] - new[i];
    }
    return sub;
}

double euclidean_norm(double* vector)
{
    double res = 0;
    for (int i=0; i<sizeof(vector); i++){
        res += pow(vector[i],2);
    }
    res = pow(res, 0.5);
    return res;
}

void* k_mean(int k, int max_iter, char *input_filename,char *output_filename)
{
    if (k <= 0 || max_iter <= 0)
    {
        printf("Invalid Input!");
        exit;
    }
    
    if (fopen(output_filename, "w") == NULL)
    {
        printf("Invalid Input!");
        exit;
    }   
    fclose(output_filename);
    char more_than_epsilon = 1;
    double epsilon = 0.001;
    int i = 0;
    double **data_points = read_file(input_filename);
    double **centroids = (double **)calloc(k,sizeof(double*));
    for (int i=0; i < k; i++)
    {
        centroids[i] = data_points[i];// how to copy address
    }
    int iteration = 0;
    while (more_than_epsilon && iteration < max_iter)
    {
        iteration += 1;
        double ***clusters = (double ***)calloc(k,sizeof(double**));
        int *clusters_indices = (int *)calloc(k,sizeof(int));
        for (double* data_point = data_points[0]; data_point != NULL; data_point+=sizeof(double*)) // foreach?
        {
            int cluster_i = 0;
            double min_euclidean_dist = (double)(INT_MAX);
            for (int i=0; i<k; i++)
            {
                double euclidean_dist = euclidean_norm(sub_vector(data_point, centroids[i]));
                if (euclidean_dist < min_euclidean_dist)
                {
                    min_euclidean_dist = euclidean_dist;
                    cluster_i = i;
                }
            }
            clusters[cluster_i][clusters_indices[cluster_i]] = data_point;
            clusters_indices[cluster_i] += 1;
            double **new_centroids = (double **)calloc(k, sizeof(double*));
            for (int i=0; i<k; i++)
            {
                new_centroids[i] = avg(clusters[i]);
            }
            double *change_vector = (double *)malloc(k*sizeof(double));
            for (int i=0; i<k ; i++)
            {
                double *sub = sub_vector(centroids[i], new_centroids[i]);
                double norm = euclidean_norm(sub);
                change_vector[i] = norm;
            }
            double dist = euclidean_norm(change_vector);
            if (dist < epsilon)
            {
                more_than_epsilon = 0;
            }
            centroids = new_centroids;
        }
    }
    FILE *opf = fopen(output_filename, "w");
    if (opf != NULL)
    {
        for (int j=0; j<k; j++)
    {
        for (int i=0; i<sizeof(centroids[j]); i++)
        {
            fprintf(opf,"%.4f" , centroids[j][i]);        
            if (i != sizeof(centroids[j]-1))
            {
                fprintf(opf,",");
            }
        }
        fprintf(opf,"\n");
    }
    fclose(opf);
    }
    
}

void* main()
{
    k_mean(3,100,"input1.text", "output_file.txt");
}

