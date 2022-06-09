import sys
import pandas as pd
import numpy as np
import mykmeanssp

np.random.seed(0)


def euclidean_dist(v1 , v2):
    return sum([((v1[i] - v2[i])**2) for i in range (len(v1))])

def part_one(k):
    dists = pd.DataFrame(index=range(data_points.shape[0]),columns=range(k))
    num_of_centroids = 1
    idx_of_first_centroid = np.random.choice(data_points.shape[0])
    centroids_indices = [indices[idx_of_first_centroid]]
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
        centroids_indices.append(indices[idx_of_centroid])
        centroids.append(data_points[idx_of_centroid])
        num_of_centroids += 1
    return (centroids, centroids_indices)


def make_data_points_from_files(filename_1, filename_2):
    table_1 = pd.read_csv(filename_1, header=None, index_col = 0)
    table_2 = pd.read_csv(filename_2, header=None, index_col = 0)
    merged_data = table_1.merge(table_2, how = 'inner', left_index = True, right_index = True)
    merged_data.index = merged_data.index.astype('int')
    merged_data = merged_data.sort_index()
    return merged_data


def main():
    k, max_iter, epsilon,  data_points_py, initial_centroids, total_vec_number, size_vec, initial_centroids_indices = args_parsing()
    if k > total_vec_number:
        invalid_input()
    initial_centroids_c = []
    for centroid in initial_centroids:
        for c in centroid:
            initial_centroids_c.append(c)
    data_points_c = []
    for data_point in data_points_py:
        for d in data_point:
            data_points_c.append(d)
    centroids = mykmeanssp.fit(k, max_iter, epsilon, data_points_c,initial_centroids_c, total_vec_number, size_vec)
    if centroids == None:
        print("An Error Has Accured in C!")
    #printing the initial choosen centroids 
    print("".join(f"{idx}," for idx in initial_centroids_indices)[:-1:])
    #printing the final centroids
    for centroid in centroids:
        print("".join(f"%.4f," %coordinate for coordinate in centroid)[:-1:])


def args_parsing():
    args = sys.argv
    if len(args) < 5 or len(args) > 6:
        invalid_input()
    k = args[1]
    if not k.isnumeric():
        invalid_input()
    k = int(k)
    if k <= 0:
        invalid_input()
    if len(args) == 5:
        max_iter = 300
        idx_of_epsilon = 2
    else:
        idx_of_epsilon = 3
        max_iter = args[2]
        if not max_iter.isnumeric():
            invalid_input()
        max_iter = int(max_iter)
        if max_iter <= 0:
            invalid_input()
    epsilon = args[idx_of_epsilon]
    """ if not epsilon.isnumeric():
        print("6\n");
        invalid_input() """
    epsilon = float(epsilon)
    if epsilon < 0:
        invalid_input()
    filename_1 = args[idx_of_epsilon + 1]
    filename_2 = args[idx_of_epsilon + 2] 
    merged_data = make_data_points_from_files(filename_1, filename_2)
    # pd.DataFrame(merged_data).to_numpy();
    total_vec_number = merged_data.shape[0]
    vec_size = merged_data.shape[1]
    global data_points
    data_points = merged_data.values
    global indices
    indices = merged_data.index
    if k > total_vec_number:
        invalid_input()
    initial_centroids, initial_centroids_indices = part_one(k)
    return (k, max_iter, epsilon, data_points, initial_centroids, total_vec_number, vec_size, initial_centroids_indices)

def invalid_input():
    print("Invalid Input!")
    exit(1)

if __name__ == '__main__':
    main()

