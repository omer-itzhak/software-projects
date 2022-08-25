import sys
from enum import Enum
import pandas as pd
import numpy as np
import mykmeanssp as spk
np.random.seed(0)

class Goal(Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5

def args_parsing():
    argv = sys.argv
    if len(argv) != 4:
        invalid_input()
    k = argv[1]
    if not k.isnumeric():
        invalid_input()
    k = int(k)
    if k < 0:
        invalid_input()
    goal = argv[2]
    try:
        goal = Goal[goal]
    except:
        invalid_input()
    input_file_name = argv[3]
    return k, goal, input_file_name

    
def kmeans_pp(k, data_points):
    dists = pd.DataFrame(index=range(data_points.shape[0]),columns=range(k))
    num_of_centroids = 1
    idx_of_first_centroid = np.random.choice(data_points.shape[0])
    centroids_indices = [idx_of_first_centroid]
    centroids = [data_points[idx_of_first_centroid]]
    while num_of_centroids != k:
        weights = []
        centroid = centroids[-1]
        for vector_idx in range(data_points.shape[0]):
            data_point = data_points[vector_idx]
            dists[num_of_centroids - 1][vector_idx] = euclidean_dist(data_point,centroid)
            weights.append(dists.iloc[vector_idx].min())
        sum_of_dists = sum(weights)
        for i in range(data_points.shape[0]):
            weights[i] = weights[i]/sum_of_dists
        idx_of_centroid = np.random.choice(data_points.shape[0], p = weights) 
        centroids_indices.append(idx_of_centroid)
        centroids.append(data_points[idx_of_centroid])
        num_of_centroids += 1
    return (centroids, centroids_indices)
    

def euclidean_dist(v1 , v2):
    return sum([((v1[i] - v2[i])**2) for i in range (len(v1))])


def invalid_input():
    print("Invalid Input!")
    exit(1)


def main():
    k, goal, filename = args_parsing()
    matrix = pd.read_csv(filename)
    vec_num = matrix.shape[0]
    vec_size = matrix.shape[1]
    if k > vec_num:
        invalid_input()
    match goal:
        case 'spk':
            t = spk.C_part(k, 1, matrix, vec_num, vec_size)
            initial_centroids, centroids_indices = kmeans_pp(k, t)
            print_vec(centroids_indices)
            centroids = spk.kmeans_c(k, t, initial_centroids, vec_num, vec_size)
            for i in range(len(centroids)):
                print_vec(centroids[i])
            return 0
        case 'wam':    
            weighted_mat = spk.C_part(k, 2, matrix, vec_num, vec_size)
            print_mat_rows(weighted_mat)
            return 0
        case 'ddg':
            ddm = spk.C_part(k, 3, matrix, vec_num, vec_size)
            print_mat_rows(ddm)
            return 0
        case 'lnorm':
            laplace = spk.C_part(k, 4, matrix, vec_num, vec_size)
            print_mat_rows(laplace)
            return 0
        case 'jacobi':
            jacobi = spk.C_part(k, 5, matrix, vec_num, vec_size)
            eigenvalues, eigenvectors = jacobi[0], jacobi[1]
            print_vec(eigenvalues)
            print_mat_cols(eigenvectors)
            return 0
            

def print_vec(vec):
    for i in range(len(vec) - 1):
        x = float("{0:.4f}".format(vec[i]))
        print(f"{x},")
    print(f"{vec[-1]}\n")
    

def print_mat_rows(mat):
    for i in range(len(mat)):
        print_vec(mat[i])

def print_mat_cols(mat):
    for i in range(len(mat)):
        for j in range(len(mat)):
            x = float("{0:.4f}".format(mat[j][i]))
            print(f"{x}")
        print("\n")


if __name__ == '__main__':
    main()