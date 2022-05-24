import sys

import numpy as np
import pandas as pd

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
    rand = np.random.seed(0)
    centroids = [rand.choice(merged_data)]
    odds = []
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
            odds.append(min_dists[i] / sum)
        centroids.append(np.random.choice(merged_data, None, False, odds))


def merge_files(filename_1, filename_2):
    table_1 = pd.read_csv(filename_1, header=None)
    table_2 = pd.read_csv(filename_2, header=None)
    table_1.rename({'0': 'indices'},axis=1,inplace=True)
    #table_2.rename(columns={"0": "indices"})

    global merged_data
    merged_data = pd.merge(table_1, table_2, on='0')
    new_file = "../" + filename_1  # TODO : make sure its working
    f = open(new_file, 'w')
    for line in range(merged_data.shape[1]):
        for i, x in enumerate(merged_data.iloc[line]):  # TODO : check if correct
            if i == 0:
                continue
            f.write('%.4f' % x)
            if i != merged_data.shape[0] - 1:  # -2 beacuse we dont need the first collum
                f.write(',')
        f.write('\n')
    f.close()
    return new_file


def main():
    # return args_parsing() # TODO : maybe different return value, depends on C extension
    merge_files(r"C:\Users\omeri\Software project\project 2\ex2\test_data\test_data\input_1_db_1.txt",
                r"C:\Users\omeri\Software project\project 2\ex2\test_data\test_data\input_1_db_2.txt")


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
    if not epsilon.isnumeric():
        invalid_input()
    epsilon = int(epsilon)
    if epsilon <= 0:
        invalid_input()
    input_filename = merge_files(args[idx_of_epsilon + 1], args[idx_of_epsilon + 2])
    output_filename = open("../" + input_filename, 'w')  # TODO : check if its working
    # use C extension: TODO merge with C
    # return k_means(k, max_iter, epsilon, input_filename, output_filename)
    return 0


def invalid_input():
    print("Invalid Input!")
    exit(1)


if __name__ == '__main__':
    main()
