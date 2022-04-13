#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

double* avg(double **cluster);
double** read_file(char *input_filename);
double* sub_vector(double *old, double *new);
double* avg(double **cluster);
double** read_file(char *input_filename);
double* sub_vector(double *old, double *new);
double* sub_vector(double *old, double *new);
double** read_file(char *input_filename);
double* sub_vector(double *old, double *new);
double euclidean_norm(double* vector);
int is_numeric(char *ch);
void* k_mean(int k, int max_iter, char *input_filename,char *output_filename);
int invalid_input();
int size_vec;


int main(int argc,char *argv[])
{
    // FILE* opf;
    // int k, max_iter, idx_of_file;
    // char *input_filename, *output_filename;
    // // validity of input check:
    // if (argc < 5 || argc > 6)
    // {
    //     invalid_input();
    // }
    // if (is_numeric(argv[1]) == 1)
    // {
    //     k = (int)(*argv[1]);
    // }
    // else
    // {
    //     invalid_input();
    // }
    // if (argc == 6)
    // {
    //     if (is_numeric(argv[2]) == 1)
    //     {
    //         max_iter = (int)(*argv[2]);
    //     }
    //     else
    //     {
    //         invalid_input();
    //     }
    //     idx_of_file = 3;
    // } 
    // else // argc = 5
    // {
    //     max_iter = 200;
    //     idx_of_file = 2;
    // }
    
    // opf = fopen(*argv[idx_of_file + 1], "w");
    // if (opf == NULL)
    // {
    //     invalid_input();
    // }
    // fclose(*argv[idx_of_file]);
    // input_filename = argv[idx_of_file];
    // output_filename = argv[idx_of_file + 1];
    // if (k <= 0 || max_iter <= 0)
    // {
    //     invalid_input();
    // }
    // return k_mean(k,max_iter,input_filename,output_filename);
    k_mean(3,100,"input_2.txt", "output_file.txt");
}
void* k_mean(int k, int max_iter, char *input_filename,char *output_filename)
{
    int i,j,iteration=0,cluster_i = 0;
    int *clusters_indices;
    FILE *opf;
    char more_than_epsilon;
    double epsilon,min_euclidean_dist,euclidean_dist,norm,dist;
    double *data_point,*change_vector,*sub;
    double **data_points,**centroids,**new_centroids;
    double ***clusters;
    printf("%s","entered k_min");
    // if (k <= 0 || max_iter <= 0)
    // {
    //     printf("Invalid Input!");
    //     exit(0);
    // }
    // check = fopen(output_filename, "w");
    // if (check == NULL)
    // {   
    //     printf("Invalid Input!");
    //     exit(0);
    // }   
    // fclose(check);
   // printf("%s","file_read");
    more_than_epsilon = 1;
    epsilon = 0.001;
    data_points = read_file(input_filename);
    centroids = (double **)calloc(k,sizeof(double*));
    for (i=0; i < k; i++)
    {
        centroids[i] = data_points[i];
    }
    while (more_than_epsilon && iteration < max_iter)
    {
        iteration += 1;
        clusters = (double ***)calloc(k,sizeof(double**));
        clusters_indices = (int *)calloc(k,sizeof(int));
        for (data_point = data_points[0]; data_point != NULL; data_point+=sizeof(double*)) /* foreach? */
        {
            
            min_euclidean_dist = (double)(INT_MAX);
            for (i=0; i<k; i++)
            {
                euclidean_dist = euclidean_norm(sub_vector(data_point, centroids[i]));
                if (euclidean_dist < min_euclidean_dist)
                {
                    min_euclidean_dist = euclidean_dist;
                    cluster_i = i;
                }
            }
            clusters[cluster_i][clusters_indices[cluster_i]] = data_point;
            clusters_indices[cluster_i] += 1;
            new_centroids = (double **)calloc(k, sizeof(double*));
            assert(new_centroids && "An Error Has Occurred");
            for (i=0; i<k; i++)
            {
                new_centroids[i] = avg(clusters[i]);
            }
            change_vector = (double *)malloc(k*sizeof(double));
            assert(change_vector && "An Error Has Occurred");
            for (i=0; i<k ; i++)
            {
                sub = sub_vector(centroids[i], new_centroids[i]);
                norm = euclidean_norm(sub);
                change_vector[i] = norm;
            }
            dist = euclidean_norm(change_vector);
            if (dist < epsilon)
            {
                more_than_epsilon = 0;
            }
            centroids = new_centroids;
        }
    }
    // write to file:
    opf = fopen(output_filename, "w");
    assert(opf && "An Error Has Occurred");
    //if (opf != NULL)
    //{
        for (j=0; j<k; j++)
        {
            for (i=0; i<size_vec; i++)
                {
                    fprintf(opf,"%.4f" , centroids[j][i]);        
                    if (i != size_vec-1)
                        {
                            fprintf(opf,"%s",",");
                        }
                }
            fprintf(opf,"\n");
        }
    fclose(opf);
    //}
    return centroids;
}
double* avg(double **cluster)
{
    int i, j, n;
    double temp;
    double *centroid;
    n = sizeof(cluster[0]);
    centroid = (double *)malloc(n*sizeof(double));
    for (i=0; i<n;i++)
    {
        temp = 0;
        for (j=0; j<(int)sizeof(cluster); j++) 
        {
            temp += cluster[j][i];
        }
        centroid[i] = temp/sizeof(cluster);
    }
    return centroid;
}
double** read_file(char *input_filename)
{  
    int index=0, line=0,total_vec_number = 0;
    double a;
    char b;
    //int i,q;
    double **data_points;

    /* count cordinate number in vector */
    FILE *ipf = fopen(input_filename, "r");
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
    for(b = fgetc(ipf) ; b != EOF ; b = fgetc(ipf)){
        if(b == '\n'){
            total_vec_number += 1;
        }
    }
    fclose(ipf);

    /* main read section */
    data_points = (double **)malloc(total_vec_number*sizeof(double*));
    if(data_points != NULL)
    {
        ipf = fopen(input_filename, "r");
        assert (ipf && "An Error Has Occurred");
        // if (ipf !=NULL)
        for(line=0;line<total_vec_number;line++)
        { 
            double *vector = malloc(size_vec*sizeof(double));
            assert(vector && "An Error Has Occurred");
            for(index=0;index<size_vec;index++)
            {    
                fscanf(ipf,"%lf",&a); 
                vector[index] = a;
                fgetc(ipf);
            }
            data_points[line] = vector;  
        }
    }            
    
    
    // for(i = 0;i <10;i++ )
    // {

    //     for(q=0;q<size_vec;q++)
    //     {
    //         printf("datapoints[%d][%d] is:%f\n" ,i,q,data_points[i][q]);
    //     }
    // }
    fclose(ipf);    
    return data_points;
}
double* sub_vector(double *old, double *new)
{
    int i;
    double *sub = (double *)malloc(sizeof(old)*sizeof(double));
    for(i=0; i<(int)sizeof(old); i++)
    {
        sub[i] = old[i] - new[i];
    }
    return sub;
}
double euclidean_norm(double* vector)
{
    int i;
    double res = 0;
    for (i=0; i<(int)sizeof(vector); i++){
        res += pow(vector[i],2);
    }
    res = pow(res, 0.5);
    return res;
}
int is_numeric(char *ch)
{
    int i = 0;
    while (ch[i] != 0){
        if(ch[i] < 48 || ch[i] > 57){
            return 0;
        }
        i++;
    }
    return 1;
}
int invalid_input()
{
    printf("Invalid Input!");
    exit(0);
}

