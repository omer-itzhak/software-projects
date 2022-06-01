import sys
import pandas as pd
import numpy as np
import mykmeanssp

np.random.seed(0)
merged_data = np.array([[]])

def sub_vector(old, new):
    return [old[i] - new[i] for i in range(len(old))]


def euclidean_norm(vector):
    res = 0
    for x in vector:
        res += x ** 2
    res = res ** 0.5
    return res


def part_one(k):
    i = 1
    idx_of_first_centroid = np.random.choice(len(merged_data), 1, replace = False)[0]
    first_centroid = merged_data[idx_of_first_centroid]
    centroids = [first_centroid]
    weights = []
    indices = [idx_of_first_centroid]
    while i != k:
        sum = 0
        min_dists = []
        for data_point in merged_data:
            min_dist = float("inf")
            for centroid in centroids:
                dist = euclidean_norm(sub_vector(data_point, centroid))
                if dist < min_dist:
                    min_dist = dist
            min_dists.append(min_dist)
            sum += min_dist
        for i in range(merged_data.shape[1]):
            weights.append(min_dists[i] / sum)
        idx_of_centroid = np.random.choice(len(merged_data), 1, replace = False, p = weights)[0]
        centroids.append(merged_data[idx_of_centroid])
        indices.append(int(centroid.index))
        i += 1
    return (centroids, indices)
        

def make_data_points_from_files(filename_1, filename_2):
    table_1 = pd.read_csv(filename_1, header=None, index_col = 0)
    table_2 = pd.read_csv(filename_2, header=None, index_col = 0)
    global merged_data
    merged_data = table_1.merge(table_2, how = 'inner', left_index = True, right_index = True)
    merged_data.index = merged_data.index.astype('int')
    return merged_data


def main():
    k, max_iter, epsilon,  data_points, initial_centroids, total_vec_number, size_vec, initial_centroids_indices = args_parsing()
    initial_centroids_c = []
    for centroid in initial_centroids:
        for c in centroid:
            initial_centroids_c.append(c)
    data_points_c = []
    for data_point in data_points:
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
    if epsilon <= 0:
        invalid_input()
    filename_1 = args[idx_of_epsilon + 1]
    filename_2 = args[idx_of_epsilon + 2] 
    data_points = make_data_points_from_files(filename_1, filename_2)
    pd.DataFrame(data_points).to_numpy();
    total_vec_number = data_points.shape[0]
    vec_size = data_points.shape[1]
    if k > total_vec_number:
        invalid_input()
    initial_centroids, initial_centroids_indices = part_one(k)
    return (k, max_iter, epsilon, data_points, initial_centroids, total_vec_number, vec_size, initial_centroids_indices)

def invalid_input():
    print("Invalid Input!")
    exit(1)

if __name__ == '__main__':
    main()

