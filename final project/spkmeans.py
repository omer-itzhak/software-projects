import sys
import pandas as pd
import numpy as np
import mykmeanssp as spk
np.random.seed(0)

def args_parsing():
    goals = {'spk','wam','ddg','lnorm','jacobi'}
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
    if goal not in goals:
        invalid_input()
    input_file_name = argv[3]
    return k, goal, input_file_name

    
def kmeans_pp(k, data_points):
    vec_num = len(data_points)
    dists = pd.DataFrame(index=range(vec_num),columns=range(k))
    num_of_centroids = 1
    idx_lst = [i for i in range(vec_num)]
    idx_of_first_centroid = np.random.randint(vec_num)
    centroids_indices = [idx_of_first_centroid]
    centroids = [data_points[idx_of_first_centroid]]
    while num_of_centroids != k:
        weights = []
        centroid = centroids[-1]
        for vector_idx in range(vec_num):
            data_point = data_points[vector_idx]
            dists[num_of_centroids - 1][vector_idx] = euclidean_dist(data_point,centroid)
            weights.append(dists.iloc[vector_idx].min())
        sum_of_dists = sum(weights)
        for i in range(vec_num):
            weights[i] = weights[i]/sum_of_dists
        idx_of_centroid = (np.random.choice(idx_lst, 1, p = weights))[0]
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
    # matrix is datapoints or sym mat
    matrix = np.loadtxt(filename, delimiter=',')
    vec_num = np.shape(matrix)[0]
    vec_size = np.shape(matrix)[1]
    mat_for_c = matrix.flatten().tolist()
    if goal == 'spk':
        if k > vec_num:
            invalid_input()
        T = spk.C_part(k, 1, mat_for_c, vec_num, vec_size)
        vec_num = len(T)
        k = len(T[0])
        vec_size = k;
        initial_centroids, centroids_indices = kmeans_pp(k, T)
        centroids_for_c = prepare_mat_for_c(initial_centroids)
        t_for_c = prepare_mat_for_c(T)
        sorted_indices = sorted(centroids_indices, reverse = 1)
        print_int_vec(sorted_indices)
        spk.fit(k, t_for_c, centroids_for_c, vec_num, vec_size)
        # centroids = spk.fit(k, t_for_c, centroids_for_c, vec_num, vec_size)
        # print_mat_rows(centroids)
    elif goal == 'wam':
        weighted_mat = spk.C_part(k, 2, mat_for_c, vec_num, vec_size)[1]
        print_mat_rows(weighted_mat)
    elif goal == 'ddg':
        ddm = spk.C_part(k, 3, mat_for_c, vec_num, vec_size)[1]
        print_mat_rows(ddm)
    elif goal == 'lnorm':
        laplace = spk.C_part(k, 4, mat_for_c, vec_num, vec_size)[1]
        print_mat_rows(laplace)
    else: # goal = 'jacobi'
        jacobi = spk.C_part(k, 5, mat_for_c, vec_num, vec_size)
        eigenvalues = jacobi[0][0]
        eigenvectors = jacobi[1]
        print_vec(eigenvalues)
        print_mat_rows(eigenvectors)
        return 0
            

def prepare_mat_for_c(mat):
    res = []
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            res.append(mat[i][j])
    return res

def print_vec(vec):
    str = ""
    for i in range(len(vec)):
        x = "{0:.4f}".format(vec[i])
        str += x
        if i != len(vec) - 1:
            str += ","
    print(str)


def print_int_vec(vec):
    str = ""
    for i in range(len(vec)):
        str += "{}".format(vec[i])
        if i != len(vec) - 1:
            str += ","
    print(str)
    
    

def print_mat_rows(mat):
    for i in range(len(mat)):
        print_vec(mat[i])


def print_mat_cols(mat):
    for i in range(len(mat)):
        str = ""
        for j in range(len(mat)):
            x = "{0:.4f}".format(mat[j][i])
            str += x
            if j != len(mat) - 1:
                str += ","
        print(str)


if __name__ == '__main__':
    main()